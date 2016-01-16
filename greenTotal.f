! version %I% of %G%
!Compute Green's functions 
subroutine calc_GreenT(m,Nx,P,GLnn,GLn1,GL1n, GRnn,GRNxn,GRnNx, Gnn,Gn1,G1n,GnNx,GNxn)
  use params
  implicit none
  logical  ldummy(1)
  integer, intent(in) :: m      !number of points along y
  integer, intent(in) :: Nx      !number of points along x
  complex(WP), intent(in) :: p(1:m) !p(k)=exp(2*i*Pi*hmag*k)
  complex(WP), intent(in) :: GLnn(1:m,1:m,1:Nx)
  complex(WP), intent(in) :: GRnn(1:m,1:m,1:Nx)
  complex(WP), intent(in) :: GLn1(1:m,1:m,1:Nx)
  complex(WP), intent(in) :: GL1n(1:m,1:m,1:Nx)
  complex(WP), intent(in) :: GRnNx(1:m,1:m,1:Nx)
  complex(WP), intent(in) :: GRNxn(1:m,1:m,1:Nx)

  complex(WP), intent(out) :: Gn1(1:m,1:m,1:Nx)
  complex(WP), intent(out) :: G1n(1:m,1:m,1:Nx)
  complex(WP), intent(out) :: Gnn(1:m,1:m,1:Nx)
  complex(WP), intent(out) :: GnNx(1:m,1:m,1:Nx)
  complex(WP), intent(out) :: GNxn(1:m,1:m,1:Nx)

  complex(WP), allocatable :: gmid(:,:)
  complex(WP), allocatable :: gnmid(:,:)
  complex(WP), allocatable :: gmid1(:,:)
  complex(WP), allocatable :: gmid2(:,:)

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
  allocate(Gmid(1:m,1:m))
  allocate(Gnmid(1:m,1:m))
  allocate(Gmid1(1:m,1:m))
  allocate(Gmid2(1:m,1:m))
  allocate(WORKI(max(lworki,1)))
  allocate(IPIV(1:m))

! calculation "total" Green's functions
! Vn,n+1(l)=-t*congj(P(l))=-tv*delta_{l,l'}*exp(cmplx(0,-2.0*PI*l*hmag,WP))
! Vn+1,n(l)=-t*P(l)=-tv*delta_{l,l'}*exp(cmplx(0,2.0*PI*l*hmag,WP))
  LDAI=m
! 4.22a,b
  do n = 1,Nx-1
    do j = 1, M
       do l = 1, M
          Gmid1(l,j) = -GLnn(l,j,n)*conjg(P(j))
          Gnmid(l,j) = -GL1n(l,j,n)*conjg(P(j))
          Gmid2(l,j) = -GRnn(l,j,n+1)*P(j)
       end do
    end do

       call zgemm('N','N',m,m,m,zone,&
       &Gnmid(1:m,1:m),size(Gmid1,dim=1),&  
       &Gmid2(1:m,1:m),size(Gmid2,dim=1),&
       &zzero,Gmid(1:m,1:m),size(Gmid,dim=1))
       call zgemm('N','N',m,m,m,zone,&
       &Gmid(1:m,1:m),size(Gnmid,dim=1),&  
       &Gnn(1:m,1:m,n),size(Gnmid,dim=1),&
       &zzero,Gnmid(1:m,1:m),size(Gnmid,dim=1))

       G1n(1:m,1:m,n) = GL1n(1:m,1:m,n)+Gnmid(1:m,1:m) !4.22b

!   Gmid = Gnmid*Gmid2
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
       &GLn1(1:m,1:m,n),size(Gmid,dim=1),&
       &zzero,Gn1(1:m,1:m,n),size(Gmid,dim=1)) !4.22a
!
       call zgemm('N','N',m,m,m,zone,&
       &Gmid(1:m,1:m),size(Gmid,dim=1),&
       &GLnn(1:m,1:m,n),size(Gmid,dim=1),&
       &zzero,Gnn(1:m,1:m,n),size(Gmid,dim=1)) !4.22c

    do j = 1, M
       do l = 1, M
          Gmid1(l,j) = -Gnn(l,j,n)*P(l)
          Gmid2(l,j) = -Gnn(l,j,n)*conjg(P(j)) 
       end do
    end do

       call zgemm('N','N',m,m,m,zone,&
       &GRNxn(1:m,1:m,n+1),size(Gmid1,dim=1),&  
       &Gmid1(1:m,1:m),size(Gmid1,dim=1),&
       &zzero,GNxn(1:m,1:m,n),size(Gnmid,dim=1))  !4.22c
       call zgemm('N','N',m,m,m,zone,&
       &Gmid2(1:m,1:m),size(Gnmid,dim=1),&  
       &GRnNx(1:m,1:m,n+1),size(Gnmid,dim=1),&
       &zzero,GnNx(1:m,1:m,n),size(Gnmid,dim=1))  !4.22d
  end do

  deallocate(Gmid)
  deallocate(Gnmid)
  deallocate(Gmid1)
  deallocate(Gmid2)
  deallocate(WORKI)
  deallocate(IPIV) 

end subroutine calc_GreenT
