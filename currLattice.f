program Sheglov
  use params
  use potential2d
  implicit none
      INTERFACE
      FUNCTION DSECND()
      DOUBLE PRECISION   DSECND
      END
      END INTERFACE
  double precision time(20)
  character*(20) timename(20)

  complex(WP), allocatable :: sigmaRT(:,:), sigmaLT(:,:), gammaRT(:,:)
  complex(WP), allocatable :: GRnn(:,:,:), GLnn(:,:,:), Gnn(:,:,:)
  real(WP), allocatable :: CuryN(:,:), CurxN(:,:)
  real(WP), allocatable :: CuryEq(:,:), CurxEq(:,:)
  complex(WP), allocatable :: GRNxn(:,:,:)
  complex(WP), allocatable :: GRnNx(:,:,:)
  complex(WP), allocatable :: GLn1(:,:,:), GL1n(:,:,:)
  complex(WP), allocatable :: Gn1(:,:,:), G1n(:,:,:)
  complex(WP), allocatable :: GNx_n(:,:,:), Gn_Nx(:,:,:)
  complex(WP), allocatable :: p(:), pc(:)
  complex(WP), allocatable :: Gmid1(:,:), Gmid2(:,:)
  complex(WP), allocatable :: Gmid(:,:)!, Gnmid(:,:)
  complex(WP), allocatable :: gammaL(:,:),gammaR(:,:)
  complex(WP), allocatable :: sigmaL(:,:),sigmaR(:,:)

  Integer nprocs, myrank, ierr, istart, iend, irank
  type(p2d_t) :: vxy
  type(p2d_t) :: vxy1
  real(WP) :: vmin, btesla, b, be, ef, hmag, dtmp, lds
  integer :: Nx,m,i,j,k,n,ief,iefmin,iefmax,ibt,ibtmax,isim,njmax,mjmax,Mc
  real(WP) :: sumd, sumUp, sumDown, y, x, Lx
  real(WP) :: sum,sled, gtot,Jmax,jx,jy, period, g1x,g2x,g3x,g1y,g2y,g3y
  complex(WP) :: emt, emkt
  real(WP) :: tv,dx, y0, x0, Umin, Umax, Umean, Utemp
  complex(WP), parameter::zzero=(0.0,0.0)
  complex(WP), parameter::zone = (1.0,0.0)
  integer isum, ixshift, iyshift
!
     vxy%xmax = 4*800.
     vxy%ymax = 4*800.
     vxy%xmin = 0
     vxy%ymin = 0
     vxy%nx   = 641!161!81
     vxy%hx   = 5
     vxy%ny   = 641!161!81
     vxy%hy   = 5
     x0 = 0.5*vxy%xmax
     y0 = 0.5*vxy%ymax
     dx = 200!0.25*x0
     call p2d_alloc(vxy,vxy%nx,vxy%ny)
     write(*,*) vxy%nx, vxy%ny,vxy%hx,vxy%hy
     Umin=100000
     Umax=-100000
     do i = 1, vxy%nx
        x = vxy%hx*(i-1)
        do j = 1, vxy%ny
           y = vxy%hy*(j-1)
           vxy%u(i,j) = 5.*exp(-0.5*((x-x0)/dx)**2)+25.*(y/y0-1.)**8
           Utemp = vxy%u(i,j)
           if(Umin.gt.utemp) Umin = Utemp
           if(Umax.lt.utemp) Umax = Utemp
        end do
     end do
     write(*,*) Umin, Umax
     goto 11
     call p2d_read_n(vxy1,'../DOSxy/uidSQad_Ns15-700246.dat')!ideal case
     iyshift = 0
     ixshift = 20
     vxy%xmax = vxy1%xmax-ixshift*vxy1%hx
     vxy%ymax = vxy1%ymax-iyshift*vxy1%hy
     vxy%xmin = vxy1%xmin+ixshift*vxy1%hx
     vxy%ymin = vxy1%ymin+iyshift*vxy1%hy
     vxy%nx   = vxy1%nx - 2*ixshift
     vxy%hx   = vxy1%hx
     vxy%ny   = vxy1%ny - 2*iyshift
     vxy%hy   = vxy1%hy
     call p2d_alloc(vxy,vxy%nx,vxy%ny)
     write(*,*) vxy%nx, vxy%ny,vxy%hx,vxy%hy
     Umin=100000
     Umax=-100000
     do i = 1, vxy%nx
        x = vxy%hx*(i-1)
        do j = 1, vxy%ny
           vxy%u(i,j) = 1000*vxy1%u(i+ixshift,j+iyshift)
           Utemp = vxy%u(i,j)
           if(Umin.gt.utemp) Umin = Utemp
           if(Umax.lt.utemp) Umax = Utemp
        end do
     end do
   call p2d_transpose(vxy)
   write(*,*) vxy%nx, vxy%ny, vxy%hx, vxy%hy
11   call p2d_write_n(vxy,'Uxybar.dat')
!  stop
  open(DTOUT1,file=DATFILE,status='replace',iostat=ierr)!GatE_B
  open(DTOUT,file=DATAFILE,status='replace',iostat=ierr)!LDSxy
  if(ierr.ne.0) then
     write(*,*) 'Error writing file ',DATAFILE
     Stop
  endif
!
  open(DTOUT2,file=DATFILE1,status='replace',iostat=ierr)
  if(ierr.ne.0) then
     write(*,*) 'Error writing file ',DATFILE1
     Stop
  endif

  ! Now vxy%u( 1:vxy%nx, 1:vxy%ny ) contains the rectangular area
  ! with potential. Let m=vxy%ny. Take the 'left column' of points,
  ! vxy%u(1,:), and compute 2*m eigenvalues w(k) and 2*m eigenvectors
  ! vr(:,k) in the uniform channel of shape vxy%u(1,:). For this we
  ! need
  ! dimensionless cell size: be,
  ! magnetic field: hmag,
  ! and Fermi energy: ef.

  b = vxy%hy                    ! or sqrt(vxy%hy * vxy%hx) ?
  be = b**2 / E00
  m = vxy%ny
  Nx = vxy%nx
  allocate(sigmaRT(1:m,1:m))
  allocate(sigmaLT(1:m,1:m))
  allocate(gammaRT(1:m,1:m))
  allocate(gnn(1:m,1:m,1:Nx))
  allocate(gn1(1:m,1:m,1:Nx))
  allocate(g1n(1:m,1:m,1:Nx))
  allocate(gNx_n(1:m,1:m,1:Nx))
  allocate(gn_Nx(1:m,1:m,1:Nx))
  allocate(gLnn(1:m,1:m,1:Nx))
  allocate(gLn1(1:m,1:m,1:Nx))
  allocate(gL1n(1:m,1:m,1:Nx))
  allocate(gRnn(1:m,1:m,1:Nx))
  allocate(gRNxn(1:m,1:m,1:Nx))
  allocate(gRnNx(1:m,1:m,1:Nx))
  allocate(gammaL(1:m,1:m))
  allocate(gammaR(1:m,1:m))
  allocate(sigmaL(1:m,1:m))
  allocate(sigmaR(1:m,1:m))
  write(*,*) "number of points along y=", m
  allocate(p(1:m))
  allocate(pc(1:m))
  tv=1./be
  allocate(Gmid(1:m,1:m))
  allocate(Gmid1(1:m,1:m))
  allocate(Gmid2(1:m,1:m))
  iEfmax = (Emax-Emin)/deltaE+1
  ibtmax = (B_max-B_min)/deltaB+1
  iEfmin=1
   Btesla=B_min
   Ef = Emin
!      if(isim.eq.2) then
!         Btesla=-Btesla
!         Curxold = CurxN
!         Curyold = CuryN
!      end if
15      hmag = btesla * b**2 * 0.006/25.
      emt = exp(cmplx(0.0,2.0*PI*hmag,WP))
      P(1) = emt
      emkt = emt
      do k = 2, M-1
         emkt = emkt * emt
         P(k) = emkt
      end do
      P(M)=emkt*emt
      Pc=conjg(P)
      write(*,*) 'before'
  call calc_gamma(m,vxy%u(1,1:m),tv,hmag,ef,Pc,sigmaL,gammaL)
  deallocate(pc)
  call calc_gamma(m,vxy%u(Nx,1:m),tv,hmag,ef,P,sigmaR,gammaR)
      write(*,*) 'gamma'
  call calc_GreenL(M,Nx,vxy%u(1:Nx,1:M),tv,hmag,ef,P,sigmaL,sigmaR,GLnn,GLn1,GL1n)
      write(*,*) 'GreenL'
  call calc_GreenR(M,Nx,vxy%u(1:Nx,1:M),tv,hmag,ef,P,sigmaL,sigmaR,GRnn,GRNxn,GRnNx)
      write(*,*) 'GreenR'
  call calc_GreenT(M,Nx,P,GLnn,GLn1,GL1n, GRnn,GRNxn,GRnNx, Gnn,Gn1,G1n,Gn_Nx,GNx_n)
      write(*,*) 'GreenT'
!--------conductance-----------------------------------
       call zgemm('N','N',m,m,m,zone,&
       &gammaR(1:m,1:m),size(gammaR,dim=1),&
       &gLn1(1:m,1:m,Nx),size(Gmid1,dim=1),&  !???
       &zzero,Gmid1(1:m,1:m),size(Gmid1,dim=1))
       Gmid2 = transpose(GLn1(1:m,1:m,Nx))
       Gmid = conjg(Gmid2)

       call zgemm('N','N',m,m,m,zone,&
       &gammaL(1:m,1:m),size(gammaL,dim=1),&
       &Gmid(1:m,1:m),size(Gmid,dim=1),&
       &zzero,Gmid2(1:m,1:m),size(Gmid2,dim=1))

       call zgemm('N','N',m,m,m,zone,&
       &Gmid1(1:m,1:m),size(Gmid1,dim=1),&
       &Gmid2(1:m,1:m),size(Gmid2,dim=1),&
       &zzero,Gmid(1:m,1:m),size(Gmid,dim=1))
    gtot=0.
    do j = 1, M
          gtot=gtot+(Gmid(j,j))
    end do
    write(DTOUT1,98) Ef, gtot, Btesla!, grtot+gtot
!    write(DTOUT,98) Ef, gtot, Btesla!, grtot+gtot
!  if(gtot.lt.0) then
!     Ef=Ef+0.01*deltaE
!     goto 15
!  end if
! --------------------------------------------
! LDS(x,y)
    write(*,*) Btesla, gtot, Ef
    lds = 0
    do n = 2, Nx-1
       do i = 1, M
           lds = lds - aimag(Gnn(i,i,n))/PI
          !write(DTOUT1,98) n*b, b*i, -aimag(Gnn(i,i,n))/PI
          write(DTOUT1,'(g13.6,1x,$)') -aimag(Gnn(i,i,n))/PI
       end do
       write (DTOUT1,*)
    end do
  write(DTOUT,98) Ef, Btesla, gtot
  write(DTOUT,98) Btesla, Ef, lds/(Nx-2)/M
  close(DTOUT)
  close(DTOUT1)

! COMPUTE CURRENT

  allocate(curxN(1:m,1:Nx))
  allocate(curyN(1:m,1:Nx))
  allocate(curxEq(1:m,1:Nx))
  allocate(curyEq(1:m,1:Nx))
  curxN=0
  curyN=0
  curxEq=0
  curyEq=0
14    do n = 2, Nx-2
       do j = 1, M
       do i = 1, M
          sigmaRT(i,j)=GRnn(i,j,n+1)*conjg(P(i))*P(j)            !5.16a
          sigmaLT(i,j)=GLnn(i,j,n-1)*conjg(P(j))*P(i)           !5.13
       end do
       end do
       sigmaR=transpose(sigmaRT(1:m,1:m))
       sigmaL=transpose(sigmaLT(1:m,1:m))
       sigmaR=conjg(sigmaR)
       sigmaR=sigmaRT(1:m,1:m)-sigmaR
       gammaRT(1:m,1:m) = cmplx(-aimag(sigmaR(1:m,1:m)),real(sigmaR(1:m,1:m))) !5.16b
!  begin A: Gnn(n)sigmaL(n)
       call zgemm('N','N',m,m,m,zone,&
       &Gnn(1:m,1:m,n),size(Gmid1,dim=1),&
       &sigmaLT(1:m,1:m),size(Gmid2,dim=1),&
       &zzero,Gmid(1:m,1:m),size(Gmid,dim=1))   !first part of 5.12a
       call zgemm('N','N',m,m,m,zone,&
       &sigmaLT(1:m,1:m),size(Gmid2,dim=1),&
       &Gnn(1:m,1:m,n),size(Gmid1,dim=1),&
       &zzero,Gmid1(1:m,1:m),size(Gmid,dim=1))   !second part of 5.12a
       do i = 1, M
          CurxEq(i,n)=2*real(Gmid(i,i)-Gmid1(i,i))          ! 5.12a
       end do

! begin B -- 5.15a
       call zgemm('N','N',m,m,m,zone,&
       &Gnn(1:m,1:m,n),size(Gmid1,dim=1),&
       &GammaRT(1:m,1:m),size(Gmid2,dim=1),&
       &zzero,Gmid(1:m,1:m),size(Gmid,dim=1))   !first half
!-----------------
       sigmaR=transpose(Gnn(1:m,1:m,n))
       gammaR(1:m,1:m) = conjg(sigmaR(1:m,1:m))  ! Gnn*
       do i = 1, M-1
          CuryEq(i,n)=-2*real(Gnn(i+1,i,n)-GammaR(i+1,i))          ! 5.12b знак - перед двойкой
       end do
!----------End of the equilibrium current calculation
!       goto 22
       sigmaL=transpose(sigmaLT(1:m,1:m))
       gammaL(1:m,1:m) = conjg(sigmaL(1:m,1:m)) !this is (sigmaLT)*
!  Gnn(n)GammaR(n)(Gnn(n))*(sigmaL(n))*
       call zgemm('N','N',m,m,m,zone,&
       &gammaR(1:m,1:m),size(Gmid1,dim=1),&
       &gammaL(1:m,1:m),size(Gmid2,dim=1),&
       &zzero,Gmid1(1:m,1:m),size(Gmid,dim=1))   !second half

       call zgemm('N','N',m,m,m,zone,&
       &Gmid(1:m,1:m),size(Gmid1,dim=1),&
       &Gmid1(1:m,1:m),size(Gmid2,dim=1),&
       &zzero,Gmid2(1:m,1:m),size(Gmid,dim=1))       !!5.15a

       do i = 1, M
          CurxN(i,n)=-2*aimag(Gmid2(i,i))          ! 5.14a
       end do

       call zgemm('N','N',m,m,m,zone,&
       &Gmid(1:m,1:m),size(Gmid1,dim=1),&
       &gammaR(1:m,1:m),size(Gmid2,dim=1),&
       &zzero,Gmid2(1:m,1:m),size(Gmid,dim=1))  !5.15b

       do i = 1, M-1
          CuryN(i,n)=2*aimag(Gmid2(i+1,i))
       end do
22  end do  ! end of Nx-loop
  deallocate(gammaL)
  deallocate(gammaR)
  deallocate(sigmaL)
  deallocate(sigmaR)

  deallocate(gnn)
  deallocate(gn1)
  deallocate(g1n)
  deallocate(gNx_n)
  deallocate(gn_Nx)
  deallocate(sigmaRT)
  deallocate(sigmaLT)
  deallocate(gammaRT)
  deallocate(GLnn)
  deallocate(GRnn)
  deallocate(GLn1)
  deallocate(GL1n)
  deallocate(GRnNx)
  deallocate(GRNxn)
  deallocate(Gmid)
  deallocate(Gmid1)
  deallocate(Gmid2)
  deallocate(p)

      write(*,*) 'Current(x,y)'
!      CuryN = CuryEQ
      CuryN = CuryN + CuryEQ
!      CurxN = CurxEQ
      CurxN = CurxN + CurxEQ
!----------------
!      CuryN = -CuryN + CuryEQ
!      CurxN = -CurxN + CurxEQ
  write(*,*) Ef,gtot, Btesla
!  if(gtot.lt.0) then
!     Ef=Ef+0.01*deltaE
!     goto 15
!  end if
!----------------
!      write(*,*) 'Ef_last=',Ef
    goto 27
    CuryN = CuryN*deltaE
    CurxN = CurxN*deltaE
        Mc=M/2
    do n = 2, Nx-3
       sumDown=0
       sumUp=0
       do i = 1, M-1
          jx = curxN(i,n+1)
          if(i.lt.Mc) sumDown=sumDown+ jx
          if(i.ge.Mc) sumUp=sumUp+ jx
       end do
          write(DTOUT1,99) n*b, -(sumDown+sumUP), -sumDown, -sumUP
    end do
24    Jmax=0.
    njmax=0
    mjmax=0
    do n = 2, Nx-2
       do i = 2, M-1
!          sum = sqrt((curxN(i,n))**2+(curyN(i-1,n))**2)
jx = 0
jy = 0
if (curxN(i,n-1).gt.0) jx=jx+curxN(i,n-1)
if (curxN(i,n).lt.0)   jx=jx-curxN(i,n)
if (curyN(i-1,n).gt.0) jy=jy+curyN(i-1,n)
if (curyN(i,n).lt.0)   jy=jy-curyN(i,n)
sum = sqrt(jx**2 + jy**2)
          if(sum.gt.Jmax) then
              Jmax=sum
              njmax = n
              mjmax = i
          endif
       end do
    end do
write(*,*) 'Jmax=',Jmax, 'b=',b,'mJmax=',mjmax,'nJmax=',njmax
write(*,*) 'Jxmax=',curxN(mjmax,njmax),'jymax=',curyN(mjmax,njmax)
27    do n = 2, Nx-3
       do i = 1, M-1
          sum = sqrt((curxN(i,n))**2+(curyN(i,n+1))**2)
jx = curxN(i,n+1)
jy = curyN(i,n)
!if (curxN(i,n+1).gt.0) jx=jx+curxN(i,n+1)
!if (curxN(i,n).lt.0)   jx=jx-curxN(i,n)
!if (curyN(i,n).gt.0) jy=jy+curyN(i,n)
!if (curyN(i-1,n).lt.0)   jy=jy-curyN(i-1,n)
!          write(DTOUT,99) n*b, b*i, jx, jy !!!!!!!!!!!!!!!!!
!         write(DTOUT,99) n*b, b*i, b*jx/Jmax, b*jy/Jmax
!         write(DTOUT,98) b*i, sum
          write(DTOUT2,'(g13.6,1x,$)') sum
       end do
       write (DTOUT2,*)
    end do
26  close(DTOUT2)
  deallocate(curxN)
  deallocate(curyN)
  deallocate(curxEq)
  deallocate(curyEq)
  stop
98  format(g13.6,3x,g13.6,3x,g13.6)
99  format(g13.6,3x,g13.6,3x,g13.6,3x,g13.6)
end program Sheglov

