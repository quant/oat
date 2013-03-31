! version %I% of %G%
!Compute Green function g11 in the semiinfinit lead
!-section is given in array v.
subroutine calc_gamma(m,v,tv,hmag,ef,p,sigma,gamma)
  use params
  implicit none
  logical  ldummy(1)
  integer, intent(in) :: m      !number of points in v
  real(WP), intent(in) :: v(m)  !potential in the channel
  real(WP), intent(in) :: tv    !cell size
  real(WP), intent(in) :: hmag  !magnetic field
  real(WP), intent(in) :: ef    !Fermi energy
  complex(WP), intent(in) :: p(m) !p(k)=exp(2*i*Pi*hmag*k)
!  complex(WP), intent(out) :: g11(1:m,1:m)
  complex(WP), intent(out) :: gamma(1:m,1:m)
  complex(WP), intent(out) :: sigma(1:m,1:m)

  complex(WP), allocatable :: g11old(:,:)
  complex(WP), allocatable :: g11(:,:)
  complex(WP), allocatable :: deltag(:,:)
  complex(WP), allocatable :: ds(:,:)
  complex(WP), allocatable :: db(:,:)
  complex(WP), allocatable :: dbv(:,:)
  complex(WP), allocatable :: Amid(:,:)
  complex(WP), allocatable :: Amidr(:,:)
  complex(WP), allocatable :: adbva(:,:)
  complex(WP), allocatable :: adbvb(:,:)
  complex(WP), allocatable :: bdbvb(:,:)
  complex(WP), allocatable :: bdbva(:,:)
  complex(WP), allocatable :: sigmaT(:,:)

  integer, allocatable :: IPIV(:)
  complex(WP), allocatable :: WORKI(:)
  
  complex(WP) :: emt, emkt, Pij, Pcij
  integer :: k, i, j, n
  real(WP) :: eh, err, rr
  integer :: info,LDAI
  integer, save :: LworkI
  complex(WP), parameter::zzero=(0.0,0.0)
  complex(WP), parameter::zone = (1.0,0.0)   

! definition  of matrixs for Green recursive functions
  if (lworki .eq. 0) lworki = 64*M !TODO: lworki as estimate on m

  allocate(g11(1:m,1:m))
  allocate(g11old(1:m,1:m))
  allocate(deltag(1:m,1:m))
  allocate(adbvb(1:m,1:m))
  allocate(adbva(1:m,1:m))
  allocate(bdbvb(1:m,1:m))
  allocate(bdbva(1:m,1:m))
  allocate(Amid(1:m,1:m))
  allocate(Amidr(1:m,1:m))
  allocate(ds(1:m,1:m))
  allocate(db(1:m,1:m))
  allocate(dbv(1:m,1:m))
  allocate(WORKI(max(lworki,1)))
  allocate(IPIV(1:m))

  LDAI=m
!dimensionless Hamiltonian: (E-H)/t; t=hbar^2/(2m_e*(1nm)^2*b^2)=E0/b^2
! d=ds; D=db; A=-delta_{l,l'}*exp(cmplx(0,- 2.0*PI*l*hmag,WP))
! B=-delta_{l,l'}*exp(cmplx(0,2.0*PI*l*hmag,WP))
  ds = (0.0_WP,0.0_WP)
  ds(1,2) = 1.
  do k = 2, M-1
     ds(k,k-1) = 1.
     ds(k,k+1) = 1.
  end do
  ds(M,M-1) = 1.
  do k = 1, M
     ds(  k,  k) = -4.0 + ( cmplx(ef,0.000001) - v(k) )/tv
!     ds(  k,  k) = -4.0 + ( cmplx(ef,0.0000001) - v(k) )/tv
  end do
  db=ds
  dbv= db
  n=1

1       call zgetrf( M, M, dbv, LDAI, IPIV, INFO )
124     call zgetri( M, dbv, LDAI, IPIV, WORKI,LWORKI,INFO )
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
  if(n.eq.1) then
!!! A=-P*; B=-P if right
!!! A=-P; B=-P* if left
    do i = 1, M
       do j = 1, M
          Pij  = P(i)*P(j)
          Pcij = conjg(P(i))*P(j)
!          Pij  = P(i)*P(j)
!          Pcij = conjg(P(i))*P(j)
          adbvb( i, j) = dbv( i, j)*Pcij
          bdbva( i, j) = dbv( i, j)*conjg(Pcij) !P(i)*conjg(P(j))
          adbva( i, j) = dbv( i, j)*conjg(Pij)  !conjg(P(i))*conjg(P(j))
          bdbvb( i, j) = dbv( i, j)*Pij
       end do
    end do
    ds = ds-adbvb
    db = db-adbvb-bdbva
    dbv = db
    n=2
    goto 1
  else 
!      adbvb=A'*dbv*B'=Amid*B'
       call zgemm('N','N',m,m,m,zone,&
       &adbva(1:m,1:m),size(adbva,dim=1),&
       &dbv(1:m,1:m),size(dbv,dim=1),&
       &zzero,Amid(1:m,1:m),size(Amid,dim=1))
       call zgemm('N','N',m,m,m,zone,&
       &Amid(1:m,1:m),size(Amid,dim=1),&
       &bdbvb(1:m,1:m),size(bdbvb,dim=1),&
       &zzero,Amidr(1:m,1:m),size(Amidr,dim=1))
       adbvb=Amidr
       ds = ds-adbvb
!      bdbva=B'*dbv*A'=Amid*A'
       call zgemm('N','N',m,m,m,zone,&
       &bdbvb(1:m,1:m),size(bdbvb,dim=1),&
       &dbv(1:m,1:m),size(dbv,dim=1),&
       &zzero,Amid(1:m,1:m),size(Amid,dim=1))
       call zgemm('N','N',m,m,m,zone,&
       &Amid(1:m,1:m),size(Amid,dim=1),&
       &adbva(1:m,1:m),size(adbva,dim=1),&
       &zzero,Amidr(1:m,1:m),size(Amidr,dim=1))
       bdbva=Amidr
       db = db-adbvb-bdbva
!      adbva=A'*dbv*A'=Amid*A'
       call zgemm('N','N',m,m,m,zone,&
       &adbva(1:m,1:m),size(adbva,dim=1),&
       &dbv(1:m,1:m),size(dbv,dim=1),&
       &zzero,Amid(1:m,1:m),size(Amid,dim=1))
       call zgemm('N','N',m,m,m,zone,&
       &Amid(1:m,1:m),size(Amid,dim=1),&
       &adbva(1:m,1:m),size(adbva,dim=1),&
       &zzero,Amidr(1:m,1:m),size(Amidr,dim=1))
       adbva=Amidr
!      bdbvb=B'*dbv*B'=Amid*B'
       call zgemm('N','N',m,m,m,zone,&
       &bdbvb(1:m,1:m),size(bdbvb,dim=1),&
       &dbv(1:m,1:m),size(dbv,dim=1),&
       &zzero,Amid(1:m,1:m),size(Amid,dim=1))
       call zgemm('N','N',m,m,m,zone,&
       &Amid(1:m,1:m),size(Amid,dim=1),&
       &bdbvb(1:m,1:m),size(bdbvb,dim=1),&
       &zzero,Amidr(1:m,1:m),size(Amidr,dim=1))
       bdbvb=Amidr
       dbv = db
  end if
  g11old = g11
  g11= ds
        call zgetrf( M, M, g11, LDAI, IPIV, INFO )
123     call zgetri( M, g11, LDAI, IPIV, WORKI,LWORKI,INFO )
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
          goto 123
        else
          print *, 'cannot continue, info=', info, ' worki(1)=', worki(1)
          stop
        endif
  if(n.gt.1) then
  deltag=g11old - g11
    err=0
    do i = 1, M
       do j = 1, M
          rr=abs(g11(i,j))
!          rr=abs(deltag(i,j))
!          if (mod(i,2).and.i.eq.j) write(*,*) i,abs(g11(i,j))!rr
          if(abs(deltag(i,j)).gt.err) err=abs(deltag(i,j))
       end do
    end do
!  write(*,*) err,n
  end if
  n=n+1
  if(err.gt.1.e-12) goto 1
!  write(*,*) err,n

  deallocate(adbvb)
  deallocate(adbva)
  deallocate(bdbvb)
  deallocate(bdbva)
  deallocate(dbv)
  deallocate(db)
  deallocate(ds)
  deallocate(Amid)
  deallocate(Amidr)
  deallocate(WORKI)
  deallocate(IPIV) 
  deallocate(g11old)
  deallocate(deltag)

  allocate(sigmaT(1:m,1:m))
    do i = 1, M
       do j = 1, M
          Pcij = conjg(P(i))*P(j)
          sigma( i, j) = g11( i, j)*Pcij
         end do
    end do
    sigmaT=transpose(sigma)   !!!!!убрала транспонирование
    do i = 1, M
       do j = 1, M
          Pcij = sigma( i, j) - conjg( sigmaT( i, j) ) 
          gamma( i, j) = cmplx(-aimag(Pcij), real(Pcij)) 
       end do
    end do
  deallocate(g11)
  deallocate(sigmaT)

end subroutine calc_gamma
