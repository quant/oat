! version %I% of %G%
! Given a permutation idx sort it so that q(idx(i)) is in increasing
! order of q. 
subroutine mysort(n, idx, q)
  use params
  implicit none
  integer, intent(in) :: n    
  integer, intent(inout) :: idx(n)
  real(WP), intent(in) :: q(*)
  
  integer :: i,j,k
  real(WP) :: c

  ! insertion sort, O(n**2)
  do i = 2, n
     k = idx(i)
     c = q(k)                   ! pick an element
     j = i - 1
     do while (j > 0 .and. c < q(idx(j)))
        idx(j+1) = idx(j)
        j = j - 1
     end do
     idx(j+1) = k
  end do
end subroutine mysort

! Given a permutation idx sort it so that q(idx(i)) is in decreasing
! order of q. 
subroutine mysort_descending(n, idx, q)
  use params
  implicit none
  integer, intent(in) :: n    
  integer, intent(inout) :: idx(n)
  real(WP), intent(in) :: q(*)
  
  integer :: i,j,k
  real(WP) :: c

  ! insertion sort, O(n**2)
  do i = 2, n
     k = idx(i)
     c = q(k)                   ! pick an element
     j = i - 1
     do while (j > 0 .and. c > q(idx(j)))  !!! here's the only
                                           !!! change from mysort
        idx(j+1) = idx(j)
        j = j - 1
     end do
     idx(j+1) = k
  end do
end subroutine mysort_descending
