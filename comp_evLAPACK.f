! version %I% of %G%
!Compute eigenvalues and eigenvectors for uniform channel whose cross
!-section is given in array v.
subroutine comp_ev1(n,v,be,hmag,ef,p,w,vr,ldvr)
  use params
  implicit none

  logical  ldummy(1)
  integer, intent(in) :: n      !number of points in v
  real(WP), intent(in) :: v(n)  !potential in the channel
  real(WP), intent(in) :: be    !cell size
  real(WP), intent(in) :: hmag  !magnetic field
  real(WP), intent(in) :: ef    !Fermi energy
  complex(WP), intent(out) :: w(2*n) !eigenvalues
  integer, intent(in) :: ldvr   !shall be >= 2*n
  complex(WP), intent(out) :: vr(ldvr,2*n) !eigenvectors
  complex(WP) :: vl(1,1) !eigenvectors
  complex(WP), intent(out) :: p(n) !p(k)=exp(2*i*Pi*hmag*k)
  integer :: info
  
  !TODO alignment and size
  complex(WP), allocatable :: aa(:,:)
  complex(WP) :: emt, emkt
  integer :: k
  integer, save :: lwork
  real(WP) :: eh
!  complex(WP) :: dummy(1)
  complex(WP), allocatable :: work(:)
  real(WP), allocatable :: rwork(:)
   

! first executable statement here !!!!!!!!!!!!!!!!!!!!!
  if (lwork .eq. 0) lwork = 4*n*32


  allocate(aa(1:2*n,1:2*n))
  allocate(work(lwork))
!  allocate(rwork(18*n)) ! for zgeev in ESSL
  allocate(rwork(4*n))
  aa = (0.0_WP,0.0_WP)

  emt = exp(cmplx(0.0,2.0*PI*hmag,WP))
  P(1) = emt
  emkt = emt
  aa(1,2) = -emkt
  do k = 2, n-1
     emkt = emkt * emt
     P(k) = emkt
     aa(k,k-1) = -emkt
     aa(k,k+1) = -emkt
  end do
  P(n)=emkt*emt
  aa(n,n-1) = -emkt * emt

  do k = 1, n
     eh = 4.0 + be * (v(k) - cmplx(ef,0.0000001))
     emt = P(k)!exp(cmplx(0,2.0*PI*k*hmag,WP))
     aa(  k,  k) = eh * emt
     aa(  k,n+k) = - emt * emt
     aa(n+k,  k) = 1.0
  end do
! ESSL:
!  call zgeev(1,aa,2*n,W, VR, 2*n, ldummy, 2*n, RWORK,10*n)
! LAPACK
123  call zgeev( JOBVL, JOBVR, 2*n, aa, 2*n, W, VL, 1, VR, 2*n, &
     &                  WORK, LWORK, RWORK, INFO )
  if (info .eq. 0) then
    if (real(work(1)) .gt. lwork) then
      print *, 'lwork adjusted: ',lwork,'->',int(real(work(1)))
      lwork = work(1)
      deallocate(work)
      allocate(work(lwork))
    endif
  else if (info == -6) then
      print *, 'lwork adjusted: ',lwork,'->',int(real(work(1)))
      lwork = work(1)
      deallocate(work)
      allocate(work(lwork))
      goto 123
  else
    print *, 'cannot continue info=', info, ' work(1)=', work(1)
    stop
  endif
! W(2*M) massiv sobstvennyh znachenii
! VR(2*M,2*M)- sobstvennye vektora 
!	if (INFO .ne. 0) then
!      	  write (0,*) 'Error: zgeev failed, info =', INFO
!	  select case (INFO)
!	  case (-12)
!      	    write (0,*) 'LWORK is too small'
!	  end select
!	end if
!       Lwork=WORK1(1)
!      write(*,*) Lwork
!  allocate(work(lwork))

  deallocate(aa)
  deallocate(RWORK)
  deallocate(WORK)
!  if (info < 0) then
!     write (STDERR, *) 'Wrong arg', -info, 'in call to zgeev'
!     stop 1
!  else if (info > 0) then
!     write (STDERR, *) 'QR failed in zgeev; info=', info
!     stop 1
!  end if

end subroutine comp_ev1
