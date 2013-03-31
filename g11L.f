! version %I% of %G%
!Compute Green function g11 in the semiinfinit lead
!-section is given in array v.
subroutine calc_sigmaL(m,v,tv,hmag,ef,sigma,gamma)
  use params
  implicit none
  logical  ldummy(1)
  integer, intent(in) :: m      !number of points in v
  real(WP), intent(in) :: v(m)  !potential in the channel
  real(WP), intent(in) :: tv    !cell size
  real(WP), intent(in) :: hmag  !magnetic field
  real(WP), intent(in) :: ef    !Fermi energy
!  complex(WP), intent(in) :: p(m) !p(k)=exp(2*i*Pi*hmag*k)
!  complex(WP), intent(out) :: g11(1:m,1:m)
  complex(WP), intent(out) :: gamma(1:m,1:m)
  complex(WP), intent(out) :: sigma(1:m,1:m)
  complex(WP), allocatable :: w(:), vr(:,:), p(:)
  real(WP), allocatable :: qa(:)

  complex(WP), allocatable :: sigmaT(:,:)
  complex(WP), allocatable :: Um(:,:), UVm(:,:), Fm(:,:), FVm(:,:), lambdam(:)

  integer, allocatable :: IPIV(:)
  complex(WP), allocatable :: WORKI(:)
  integer, allocatable :: nume(:)
  integer :: Idown,Iup
  
  complex(WP) :: emt, emkt, Pij, Pcij, csumm
  integer :: k, i, j, n
  real(WP) :: eh, err, rr
  integer :: info,LDAI
  integer, save :: LworkI
  complex(WP), parameter::zzero=(0.0,0.0)
  complex(WP), parameter::zone = (1.0,0.0)   

! definition  of matrixs for Green recursive functions
  if (lworki .eq. 0) lworki = 64*M !TODO: lworki as estimate on m
  allocate(nume(1:2*m))
  allocate(w(1:2*m))
  allocate(p(1:2*m))
  allocate(qa(1:2*m))
  allocate(vr(1:2*m,1:2*m))
  call comp_ev1(M, v, 1./tv, hmag, ef, P, w, vr, size(vr,dim=1))
  ! Now we SORT eigenvalues, doing this in two stages.
  ! <SORT1> Build index array 'nume' so that abs(w(nume(1:2*m))) is
  ! sorted in ascending order.
  do j = 1, 2*m
     nume(j) = 2*m - j + 1      !initialization
  end do
  qa(1:2*m) = abs( w(1:2*m) )
  call mysort(2*m, nume, qa)

  ! <SORT2> Order the eigenvalues with |w| too close to 1.0.
  ! Remember that qa is array of |w|, and 
  ! qa(nume(1)) <= qa(nume(2)) <= ... <= qa(nume(2*m)). 
  ! Find lower and upper indexes where qa(nume(k)) close to 1.0.
    call limits(2*m, nume, qa, Idown, Iup)
! nado uporyadochit' v ubyvauishem poryadke
    call mysort_descending(Iup-Idown+1, nume(Idown:Iup), aimag(w(1:2*m)))
  ! Compute array qa: qa(j) = Im[ ln w(j) ] - alpha
  qa(1:2*m) = aimag( log(w(1:2*m)) ) - PI * hmag * (m+1)
!-----------------------------------
  allocate(Um(1:m,1:m))
  allocate(UVm(1:m,1:m))
  allocate(lambdam(1:m))
  allocate(WORKI(max(lworki,1)))
  allocate(IPIV(1:m))

  LDAI=m
    do j=m+1,2*m
       k=nume(j)
       Um(1:m,j-m) = Vr(1:m,k)
       lambdam(j-m) = w(k)
!       write(*,*) abs(w(k)), w(k), j-m
    end do
       UVm=Um
        call zgetrf( M, M, UVm, LDAI, IPIV, INFO )
        call zgetri( M, Uvm, LDAI, IPIV, WORKI,LWORKI,INFO )
  allocate(Fm(1:m,1:m))
  allocate(FVm(1:m,1:m))
        do k=1,M
        do i=1,M
           csumm=(0.0_WP,0.0_WP)
           do j=1,M
              csumm=csumm+Um(i,j)*lambdam(j)*UVm(j,k)
           end do
        Fm(i,k)=csumm
        end do
        end do
        FVm=Fm
        call zgetrf( M, M, FVm, LDAI, IPIV, INFO )
        call zgetri( M, FVm, LDAI, IPIV, WORKI,LWORKI,INFO )
!        sigma = -P*FVm
        do i=1,M
           do j=1,M
              sigma(j,i) = - P(j)*FVm(j,i)
           end do
        end do
  deallocate(WORKI)
  deallocate(IPIV) 
  deallocate(Fvm)
  deallocate(Fm) 
  deallocate(Um) 
  deallocate(UVm) 
  deallocate(Vr)
  deallocate(P)
  deallocate(lambdam)

  allocate(sigmaT(1:m,1:m))
  sigmaT=transpose(sigma) 
    do i = 1, M
       do j = 1, M
          Pcij = sigma( i, j) - conjg( sigmaT( i, j) ) 
          gamma( i, j) = cmplx(-aimag(Pcij), real(Pcij)) 
       end do
    end do
  deallocate(sigmaT)

end subroutine calc_sigmaL
