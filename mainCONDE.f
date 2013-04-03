! -*- f90 -*-
program Sheglov
  !  include 'mpif.h'
  !  use mpi
  use params
  use potential2d

  implicit none

  interface
     subroutine calc_gamma(m,v,tv,hmag,ef,p,sigma,gamma)
       use params
       integer, intent(in) :: m      !number of points in v
       real(WP), intent(in) :: v(m)  !potential in the channel
       real(WP), intent(in) :: tv    !cell size
       real(WP), intent(in) :: hmag  !magnetic field
       real(WP), intent(in) :: ef    !Fermi energy
       complex(WP), intent(in) :: p(m) !p(k)=exp(2*i*Pi*hmag*k)
       complex(WP), intent(out) :: gamma(1:m,1:m)
       complex(WP), intent(out) :: sigma(1:m,1:m)
     end subroutine calc_gamma
  end interface

  complex(WP), allocatable :: GLn1(:,:,:), GLnn(:,:,:), GL1n(:,:,:)
  complex(WP), allocatable :: p(:), pc(:)
  real(WP), allocatable :: vL(:), vR(:)
  complex(WP), allocatable :: Gn1_N(:,:), Gn1_NTs(:,:)
  complex(WP), allocatable :: Gmid1(:,:), Gmid2(:,:)
  complex(WP), allocatable :: gammaL(:,:),gammaR(:,:)
  complex(WP), allocatable :: sigmaL(:,:),sigmaR(:,:)

  Integer nprocs, myrank, ierr, istart, iend, irank
  type(p2d_t) :: vxy
  type(p2d_t) :: vxy1,vxy2
  real(WP) :: vmin, btesla, b, be, ef, hmag, dtmp, Ewc, Ef0
  integer :: Nx,m,i,j,k,n,ief,iefmin,iefmax,ibt,ibtmax,isim,ixshift,iyshift
  real(WP) :: sled, gtot
  complex(WP) :: emt, emkt
  real(WP) :: tv
  complex(WP), parameter::zzero=(0.0,0.0)
  complex(WP), parameter::zone = (1.0,0.0)

  ! read precomputed 2D potential into array
  call p2d_read(vxy1,'u-new-map.dat')
  write(*,*) vxy1%nx, vxy1%ny, vxy1%hx,vxy1%hy
  !  write(*,*) vxy%nx, vxy%ny, vxy%hx,vxy%hy
  !  vxy%u = (vxy%u+vxy1%u+vxy2%u)/3.
  !scale from eV to meV
  vxy1%u = vxy1%u * 1000.
  call p2d_write_xyz(vxy1,'test1.dat')
  iyshift=30!20!20!60!20
  ixshift=0!30!24!10!60!15
  vxy%xmax=vxy1%xmax-ixshift*vxy1%hx
  vxy%ymax=vxy1%ymax-iyshift*vxy1%hy
  vxy%xmin=vxy1%xmin+ixshift*vxy1%hx
  vxy%ymin=vxy1%ymin+iyshift*vxy1%hy
  vxy%nx=vxy1%nx-2*ixshift!41!81!
  vxy%ny=vxy1%ny-2*iyshift
  vxy%hy=vxy1%hy
  vxy%hx=vxy1%hx
  call p2d_alloc(vxy,vxy%nx,vxy%ny)
  write(*,*) vxy%nx, vxy%ny, vxy%hx,vxy%hy
  !  call p2d_alloc(vxy,vxy1%nx,vxy1%ny)
  do i=1,vxy%nx  ! одно сужение, вход
     do j=1,vxy%ny
        vxy%u(i,j) = vxy1%u(i+ixshift,j+iyshift)
     end do
  end do
  call p2d_remesh(vxy,2*vxy%nx-1,2*vxy%ny-1)
  write(*,*) vxy%nx, vxy%ny, vxy%hx,vxy%hy
  call p2d_write_xyz(vxy,'test2.dat')
  !  stop
5 open(DTOUT,file=DATAFILE,status='replace',iostat=ierr)
  write (DTOUT,'(1x,t1,''G 11 3 176'')')
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
  allocate(vL(1:m))
  allocate(vR(1:m))
  allocate(gLn1(1:m,1:m,1:Nx))
  allocate(gLnn(1:m,1:m,1:Nx))
  allocate(gL1n(1:m,1:m,1:Nx))
  allocate(gn1_NTs(1:m,1:m))
  allocate(gn1_N(1:m,1:m))

  allocate(gammaL(1:m,1:m))
  allocate(gammaR(1:m,1:m))
  allocate(sigmaL(1:m,1:m))
  allocate(sigmaR(1:m,1:m))
  write(*,*) "number of points along y=", m
  allocate(p(1:m))
  allocate(pc(1:m))
  tv=1./be
  allocate(Gmid1(1:m,1:m))
  allocate(Gmid2(1:m,1:m))
  iEfmax = (Emax-Emin)/deltaE+1
  ibtmax = (B_max-B_min)/deltaB+1
  iEfmin=1
  Ef =-0.
  vL(1:m) = vxy%u(1,1:vxy%ny);
  vR(1:m) = vxy%u(Nx,1:vxy%ny);

  do ief=iEfmin,iEfmax
     Ef=Emin+deltaE*(iEf-1)
     !   Ef0=Emin+deltaE*(iEf-1)
15   do ibt=1,ibtmax
        Btesla=B_min+deltaB*(ibt-1)
        !      Ewc=4.*PI*1.602/6.625*1.0e-4*E00*Btesla
        !      Ef=Ef0!+Ewc
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
        !left transmission!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call calc_gamma(vxy%ny,vL(1:m),tv,hmag,ef,Pc,sigmaL,gammaL)
        call calc_gamma(vxy%ny,vR(1:m),tv,hmag,ef,P,sigmaR,gammaR)
        call calc_GreenL(M,Nx,vxy%u(1:Nx,1:M),tv,hmag,ef,P,sigmaL,sigmaR,GLnn,GLn1,GL1n)
        do j = 1, M
           do i = 1, M
              gn1_N(i,j)=GLn1(i,j,Nx)
           end do
        end do
        call zgemm('N','N',m,m,m,zone,&
             &gammaR(1:m,1:m),size(gammaR,dim=1),&
             &gn1_N(1:m,1:m),size(gn1_N,dim=1),&  !???
             &zzero,Gmid1(1:m,1:m),size(Gmid1,dim=1))
        gn1_NTs=transpose(gn1_N)
        do j = 1, M
           do i = 1, M
              gn1_N(i,j)=conjg(Gn1_NTs(i,j))
           end do
        end do
        call zgemm('N','N',m,m,m,zone,&
             &gammaL(1:m,1:m),size(gammaL,dim=1),&
             &gn1_N(1:m,1:m),size(gn1_N,dim=1),&
             &zzero,Gmid2(1:m,1:m),size(Gmid2,dim=1))

        call zgemm('N','N',m,m,m,zone,&
             &Gmid1(1:m,1:m),size(Gmid1,dim=1),&
             &Gmid2(1:m,1:m),size(Gmid2,dim=1),&
             &zzero,Gn1_N(1:m,1:m),size(Gn1_N,dim=1))
        sled=0.
        do j = 1, M
           sled=sled+(Gn1_N(j,j))
        end do
        if(sled.lt.0) then
           Ef=Ef+0.001*deltaE
           goto 15
        end if
        !    write(DTOUT,99) Ef, sled, Btesla!, grtot+gtot
        write(DTOUT,98) Btesla, sled, Ef!, grtot+gtot
        write(*,*) Btesla, sled, Ef
     end do
  end do
  close(DTOUT)
  deallocate(vL)
  deallocate(vR)
  deallocate(gammaL)
  deallocate(gammaR)
  deallocate(sigmaL)
  deallocate(sigmaR)
  deallocate(Gn1_N)
  deallocate(Gn1_NTs)
  deallocate(GLnn)
  deallocate(GLn1)
  deallocate(GL1n)
  deallocate(Gmid1)
  deallocate(Gmid2)
  deallocate(p)
  stop
99 format(g13.6,3x,g13.6,3x,g13.6,3x,g13.6)
98 format(g13.6,3x,g15.8,3x,g13.6)
end program Sheglov
