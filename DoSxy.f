program LDSXY
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

  complex(WP), allocatable :: GRnn(:,:,:), GLnn(:,:,:), Gnn(:,:,:)
  complex(WP), allocatable :: GRNxn(:,:,:)
  complex(WP), allocatable :: GRnNx(:,:,:)
  complex(WP), allocatable :: GLn1(:,:,:), GL1n(:,:,:)
  complex(WP), allocatable :: Gn1(:,:,:), G1n(:,:,:)
  complex(WP), allocatable :: GNx_n(:,:,:), Gn_Nx(:,:,:)
  complex(WP), allocatable :: p(:), pc(:)
  complex(WP), allocatable :: Gmid1(:,:), Gmid2(:,:)
  complex(WP), allocatable :: Gmid(:,:)
  complex(WP), allocatable :: gammaL(:,:),gammaR(:,:)
  complex(WP), allocatable :: sigmaL(:,:),sigmaR(:,:)

  type(p2d_t) :: vxy
  type(p2d_t) :: vxy1
  real(WP) :: vmin, btesla, b, be, ef, hmag, dtmp, lds
  integer :: Nx,m,i,j,k,n,ief,iefmin,iefmax,ibt,ibtmax,ierr
  real(WP) :: sumd, sumUp, sumDown, y, x, Lx
  real(WP) :: sum,sled, gtot
  complex(WP) :: emt, emkt
  real(WP) :: tv, V0, hx, hy, dy, U0, ddy, Umin, Umax, Umean, Utemp
  complex(WP), parameter::zzero=(0.0,0.0)
  complex(WP), parameter::zone = (1.0,0.0)   
  integer Nper, N_x, N_y, isum, ixshift, iyshift
  call p2d_read_n(vxy1,'../MyPoisson3d/uidSQad_Ns15-700246.dat')!ideal case
!     call p2d_mesh_thin_out(vxy1) !!!!!!
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
     sum = 0
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
!   call p2d_write(vxy,'UxyLDS.dat')
  open(DTOUT1,file=DATFILE,status='replace',iostat=ierr)!GatE_B
  open(DTOUT,file=DATAFILE,status='replace',iostat=ierr)!LDSxy
  if(ierr.ne.0) then
     write(*,*) 'Error writing file ',DATAFILE
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
      Ef=Emin
      hmag = btesla * b**2 * 0.006/25.
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
!  time(1) = dsecnd() ; timename(1) = 'start'
  call calc_gamma(m,vxy%u(1,1:m),tv,hmag,ef,Pc,sigmaL,gammaL)
  call calc_gamma(m,vxy%u(Nx,1:m),tv,hmag,ef,P,sigmaR,gammaR)
      write(*,*) 'gamma'
!  time(2) = dsecnd() ; timename(2) = 'gamma'
  call calc_GreenL(M,Nx,vxy%u(1:Nx,1:M),tv,hmag,ef,P,sigmaL,sigmaR,GLnn,GLn1,GL1n)
      write(*,*) 'GreenL'
  call calc_GreenR(M,Nx,vxy%u(1:Nx,1:M),tv,hmag,ef,P,sigmaL,sigmaR,GRnn,GRNxn,GRnNx)
      write(*,*) 'GreenR'
  deallocate(sigmaL)
  deallocate(sigmaR)
  call calc_GreenT(M,Nx,P,GLnn,GLn1,GL1n, GRnn,GRNxn,GRnNx, Gnn,Gn1,G1n,Gn_Nx,GNx_n)
      write(*,*) 'GreenT'
!  time(3) = dsecnd() ; timename(3) = 'greenLRT'
  deallocate(p)
  deallocate(gn1)
  deallocate(g1n)
  deallocate(gNx_n)
  deallocate(gn_Nx)
  deallocate(GLnn)
  deallocate(GRnn)
  deallocate(GL1n)
  deallocate(GRnNx)
  deallocate(GRNxn)

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

       deallocate(gammaL)
       deallocate(gammaR)
       deallocate(GLn1)

       call zgemm('N','N',m,m,m,zone,&
       &Gmid1(1:m,1:m),size(Gmid1,dim=1),&
       &Gmid2(1:m,1:m),size(Gmid2,dim=1),&
       &zzero,Gmid(1:m,1:m),size(Gmid,dim=1))
       deallocate(Gmid1)
       deallocate(Gmid2)

    gtot=0.
    do j = 1, M
          gtot=gtot+(Gmid(j,j))
    end do
    write(*,*) Btesla, gtot, Ef
       deallocate(Gmid)
!    time(4) = dsecnd() ; timename(4) = 'gtot'
    write(DTOUT1,98) Ef, Btesla, gtot
    lds = 0 
    do n = 2, Nx-1
       do i = 1, M
           lds = lds - aimag(Gnn(i,i,n))/PI
          !write(DTOUT1,98) n*b, b*i, -aimag(Gnn(i,i,n))/PI 
          write(DTOUT1,'(g13.6,1x,$)') -aimag(Gnn(i,i,n))/PI 
       end do
       write (DTOUT1,*)
    end do
  deallocate(Gnn)
  write(DTOUT,98) Ef, Btesla, gtot
  write(DTOUT,98) Btesla, Ef, lds/(Nx-2)/M
  close(DTOUT)
  close(DTOUT1)
  stop
98  format(g13.6,3x,g13.6,3x,g13.6)
99  format(g13.6,3x,g13.6,3x,g13.6,3x,g13.6)
end program LDSXY
                         
