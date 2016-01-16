! version %I% of %G%
subroutine limits(n, idx, q, i, j)
  use params
  implicit none
  integer, intent(in) :: n    
  integer, intent(in) :: idx(n)
  integer, intent(out) :: i
  integer, intent(out) :: j
  real(WP), intent(in) :: q(*)
  ! <SORT2> Order the eigenvalues with |w| too close to 1.0.
  ! Remember that qa is array of |w|, and 
  ! q(idx(1)) <= q(idx(2)) <= ... <= q(idx(n)). 
  ! Find lower and upper indexes where q(idx(k)) close to 1.0.
  i = n/2 ! lower index
  do while (i >   1 .and. q(idx(i)) > 1.0 - SQRTEPS)
     i = i - 1
  end do
  do while (i < n .and. q(idx(i)) < 1.0 - SQRTEPS)
     i = i + 1
  end do
  j = n/2 ! upper index
  do while (j < n .and. q(idx(j)) < 1.0 + SQRTEPS)
     j = j + 1
  end do
  do while (j >   1 .and. q(idx(j)) > 1.0 + SQRTEPS)
     j = j - 1
  end do
end subroutine limits
