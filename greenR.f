! version %I% of %G%
!Compute Green's functions 
subroutine calc_GreenR(m,Nx,v,tv,hmag,ef,p,sigmaL,sigmaR,Gnn,GNxn,GnNx)
  use params
  implicit none
  logical  ldummy(1)
  integer, intent(in) :: m      !number of points along y
  integer, intent(in) :: Nx      !number of points along x
  real(WP), intent(in) :: v(1:Nx,1:M)  !potential in the channel
  real(WP), intent(in) :: tv    !cell size
  real(WP), intent(in) :: hmag  !magnetic field
  real(WP), intent(in) :: ef    !Fermi energy
  complex(WP), intent(in) :: p(1:m) !p(k)=exp(2*i*Pi*hmag*k)
  complex(WP), intent(out) :: GnNx(1:m,1:m,1:Nx)
  complex(WP), intent(out) :: GNxn(1:m,1:m,1:Nx)
  complex(WP), intent(out) :: Gnn(1:m,1:m,1:Nx)
  complex(WP), intent(in) :: sigmaL(1:m,1:m)
  complex(WP), intent(in) :: sigmaR(1:m,1:m)

  complex(WP), allocatable :: g1i(:,:)
  complex(WP), allocatable :: gnNxold(:,:)
  complex(WP), allocatable :: gNxnold(:,:)
  complex(WP), allocatable :: gnnold(:,:)
!  complex(WP), allocatable :: gnn(:,:,:)
  complex(WP), allocatable :: g1(:,:,:)
  complex(WP), allocatable :: g1new(:,:)
  complex(WP), allocatable :: gmid(:,:)
  complex(WP), allocatable :: gnmid(:,:)
  complex(WP), allocatable :: gmid1(:,:)
  complex(WP), allocatable :: gmid2(:,:)
  complex(WP), allocatable :: EmHi(:,:)

  integer, allocatable :: IPIV(:)
  complex(WP), allocatable :: WORKI(:)
  
  complex(WP) :: emt, emkt, Pij, Pcij, clj
  integer :: k, i, j, n, l
  real(WP) :: eh, err, rr, sled
  integer :: info,LDAI
  integer, save :: LworkI
  complex(WP), parameter::zzero=(0.0,0.0)
  complex(WP), parameter::zone = (1.0,0.0)   
! definition  of matrixs for Green recursive functions
  if (lworki .eq. 0) lworki = 64*M !TODO: lworki as estimate on m
  allocate(GnNxold(1:m,1:m))
  allocate(GNxnold(1:m,1:m))
  allocate(Gnnold(1:m,1:m))
  allocate(G1new(1:m,1:m))
  allocate(Gmid(1:m,1:m))
  allocate(Gnmid(1:m,1:m))
  allocate(Gmid1(1:m,1:m))
  allocate(Gmid2(1:m,1:m))
  allocate(g1(1:m,1:m,1:Nx))
!  allocate(gnn(1:m,1:m,1:Nx))
  allocate(g1i(1:m,1:m))
  allocate(EmHi(1:m,1:m))
  allocate(WORKI(max(lworki,1)))
  allocate(IPIV(1:m))

  LDAI=m
! EmHi = E - Hi
  Emhi = (0.0_WP,0.0_WP)

! First we calculate Green's function G1 for every isolated column i=1,2,..,N
  do i = 1, Nx
  Emhi = (0.0_WP,0.0_WP)
     do k = 1, M
        EmHi(  k,  k) = -4.0 + ( cmplx(ef,0.0000001) - v( i, k ) )/tv
     end do
  EmHi(1,2) = 1.
  do k = 2, M-1
     EmHi(k,k-1) = 1.
     EmHi(k,k+1) = 1.
  end do
  EmHi(M,M-1) = 1.
  g1i=EmHi
  if(i.eq.1) g1i=g1i-sigmaL
  if(i.eq.Nx) g1i=g1i-sigmaR
        call zgetrf( M, M, g1i, LDAI, IPIV, INFO )
124     call zgetri( M, g1i, LDAI, IPIV, WORKI,LWORKI,INFO )
        if (info .eq. 0) then 
          if (real(worki(1)) .gt. lworki) then
            print *, 'lworki adjusted: ',lworki,'->',int(real(worki(1)))
            lworki = worki(1) 
            deallocate(worki)
            allocate(worki(lworki))
          endif
        else if (info .eq. -6) then 
          print *, 'lworki adjusted: ',lworki,'->',int(real(worki(1)))
          lworki = worki(1) 
          deallocate(worki)
           allocate(worki(lworki))
          goto 124
        else
          print *, 'cannot continue, info=', info, ' worki(1)=', worki(1)
          stop
        endif
     do j = 1, M
        do k = 1, M
           G1(  k,  j, i) = g1i( k, j)
        end do
     end do
  end do
! Then we calculate "right" Green's functions
! Vn,n+1(l)=-t*congj(P(l))=-tv*delta_{l,l'}*exp(cmplx(0,-2.0*PI*l*hmag,WP))
! Vn+1,n(l)=-t*P(l)=-tv*delta_{l,l'}*exp(cmplx(0,2.0*PI*l*hmag,WP))
  do j = 1, M
     do l = 1, M
        clj  = G1(l,j,Nx)
        GNxn(l,j,Nx) = clj    
        GnNx(l,j,Nx) = clj
        Gnn(l,j,Nx) = clj
        Gnnold(l,j)= clj
        GNxnold(l,j)= clj
        GnNxold(l,j)= clj
!        write(*,*) clj
     end do
  end do
!  do n = 1, Nx-1
  do n = Nx-1, 1, -1
    do j = 1, M
       do l = 1, M
          clj  = G1(l,j,n) ! isolated Green's function of n+1 column
          G1new(l,j)  = clj
!        write(*,*) clj
          Gmid1(l,j) = -clj*conjg(P(j))
          Gmid2(l,j) = -Gnnold(l,j)*P(j)
       end do
    end do
!   Gmid = Gmid1*Gmid2
       call zgemm('N','N',m,m,m,zone,&
       &Gmid1(1:m,1:m),size(Gmid1,dim=1),&
       &Gmid2(1:m,1:m),size(Gmid2,dim=1),&
       &zzero,Gmid(1:m,1:m),size(Gmid,dim=1))
!   Gmid => (1-Gmid)^{-1}
       Gmid = -Gmid
    do j = 1, M
          Gmid(j,j) = 1.+Gmid(j,j)
    end do
        call zgetrf( M, M, Gmid, LDAI, IPIV, INFO )
125     call zgetri( M, Gmid, LDAI, IPIV, WORKI,LWORKI,INFO )
        if (info .eq. 0) then 
          if (real(worki(1)) .gt. lworki) then
            print *, 'lworki adjusted: ',lworki,'->',int(real(worki(1)))
            lworki = worki(1) 
            deallocate(worki)
            allocate(worki(lworki))
          endif
        else if (info .eq. -6) then 
          print *, 'lworki adjusted: ',lworki,'->',int(real(worki(1)))
          lworki = worki(1) 
          deallocate(worki)
           allocate(worki(lworki))
          goto 125
        else
          print *, 'cannot continue, info=', info, ' worki(1)=', worki(1)
          stop
        endif

!   Gnn(n+1)= Gmid*G1(n+1)
       call zgemm('N','N',m,m,m,zone,&
       &Gmid(1:m,1:m),size(Gmid,dim=1),&
       &G1new(1:m,1:m),size(G1new,dim=1),&
       &zzero,Gnmid(1:m,1:m),size(Gnmid,dim=1))
    do j = 1, M
       do l = 1, M
          Gnn(l,j,n)  = Gnmid(l,j) ! Green's function G(n+1,n+1) for n+1 strip
!        write(*,*) Gnmid(l,j)
       end do
    end do
    Gnnold = Gnmid
!   GNxn(n)= GNxn(n+1)*Gnmid
    do j = 1, M
       do l = 1, M
          Gmid1(l,j) = -Gnmid(l,j)*P(l)      ! ��������� -, ������ ��� Vn+1,n=-P;Vn,n+1=-P*!!!
          Gmid2(l,j) = -Gnmid(l,j)*conjg(P(j)) 
       end do
    end do
       call zgemm('N','N',m,m,m,zone,&
       &GNxnold(1:m,1:m),size(GNxnold,dim=1),&  
       &Gmid1(1:m,1:m),size(Gmid1,dim=1),&
!       &Gmid1(1:m,1:m),size(Gmid1,dim=1),&
!       &Gn1old(1:m,1:m),size(Gn1old,dim=1),&  
       &zzero,Gnmid(1:m,1:m),size(Gnmid,dim=1))

!   GnNx(n)=Gnmid*GnNx(n+1)
       call zgemm('N','N',m,m,m,zone,&
       &Gmid2(1:m,1:m),size(Gmid2,dim=1),&
       &GnNxold(1:m,1:m),size(GnNxold,dim=1),&
       &zzero,Gmid1(1:m,1:m),size(Gmid1,dim=1))
    do j = 1, M
       do l = 1, M
          GNxn(l,j,n)  = Gnmid(l,j) ! Green's function G(n+1,1) for n+1 strip 
          GnNx(l,j,n)  = Gmid1(l,j) ! Green's function G(1,n+1) for n+1 strip 
       end do
    end do
    GNxnold = Gnmid
    GnNxold = Gmid1
  end do

  deallocate(GnNxold)
  deallocate(GNxnold)
  deallocate(Gnnold)
  deallocate(G1new)
  deallocate(Gmid)
  deallocate(Gnmid)
  deallocate(Gmid1)
  deallocate(Gmid2)
  deallocate(g1)
  deallocate(g1i)
  deallocate(EmHi)
  deallocate(WORKI)
  deallocate(IPIV) 

end subroutine calc_GreenR
