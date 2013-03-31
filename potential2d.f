! version %I% of %G%
! This module implements p2d_t structure to hold the approximation of
! a function (i.e. potential) on a 2d regular rectangular mesh.
module potential2d

  use params
  implicit none

  ! description of a 2d-potential
  type p2d_t
     real(WP) :: xmin = 0, xmax = 0
     real(WP) :: ymin = 0, ymax = 0
     real(WP) :: hx = 0, hy = 0
     integer :: nx = 0, ny = 0
     real(WP), pointer :: u(:,:) => null()
  end type p2d_t

contains

  !==================================================
  !allocate array for potential
  !TODO properly aligned
  !TODO avoid unnecessary reallocations
  !TODO initialization of xmin,...
  subroutine p2d_alloc(u,nx,ny)
    implicit none
    type(p2d_t), intent(inout) :: u
    integer, intent(in) :: nx, ny
    if (associated(u%u)) deallocate(u%u)
    allocate(u%u(1:nx, 1:ny))
    u%nx = nx
    u%ny = ny
    u%hx = (u%xmax - u%xmin) / (u%nx - 1)
    u%hy = (u%ymax - u%ymin) / (u%ny - 1)
  end subroutine p2d_alloc

  !==================================================
  ! read p2d_t structure from a text file
  ! Format of the input file:
  ! line 1: nx xmin xmax
  ! line 2: ny ymin ymax
  ! line 3: u(1,1)
  ! line 4: u(2,1)
  ! ...
  ! TODO handle '-' for stdin?
  subroutine p2d_read(u,fname)
    implicit none
    type(p2d_t), intent(inout) :: u
    character(*), intent(in) :: fname

    !local parameters and variables
    integer, parameter :: UN = 11
    integer :: ierr,i,j

    !..................................................
    open(unit=UN,file=fname,status='unknown',iostat=ierr)
    if (ierr.ne.0) then
       write (STDERR,*) 'Error opening file', fname, '; ierr=', ierr
       stop 1
    end if

    read(UN,*) u%nx, u%xmin, u%xmax
    read(UN,*) u%ny, u%ymin, u%ymax

    call p2d_alloc(u,u%nx,u%ny)
    do i = 1, u%nx
       do j = 1, u%ny
!    do i = 1, u%nx
          read(UN,*) u%u(i,j)
       end do
    end do

    close(UN)
  end subroutine p2d_read

  !==================================================
  ! write p2d_t structure to a text file
  ! Format of the output file:
  ! line 1: nx xmin xmax
  ! line 2: ny ymin ymax
  ! line 3: u(1,1)
  ! line 4: u(2,1)
  ! ...
  !TODO handle '-' for stdout?
  subroutine p2d_write(u,fname)
    implicit none
    type(p2d_t), intent(in) :: u
    character(*), intent(in) :: fname

    !local parameters and variables
    integer, parameter :: UN = 11
    integer :: ierr,i,j

    !..................................................
    open (unit=UN,file=fname,status='replace',iostat=ierr)
    if (ierr.ne.0) then
       write (STDERR,*) 'Error opening file', fname, '; ierr=', ierr
       stop 1
    end if

    write(UN,*) u%nx, u%xmin, u%xmax
    write(UN,*) u%ny, u%ymin, u%ymax

    do i = 1, u%nx
       do j = 1, u%ny
!       do i = 1, u%nx
          write(UN,*) u%u(i,j)
       end do
    end do

    close(UN)
  end subroutine p2d_write

  !==================================================
  ! write p2d_t structure to a text file for surfer
  ! Format of the output file:
  ! x1 y1 u(1,1)
  ! x2 y1 u(2,1)
  ! ...
  !TODO handle '-' for stdout?
  subroutine p2d_write_xyz(u,fname)
    implicit none !u%nx = 5
    type(p2d_t), intent(in) :: u
    character(*), intent(in) :: fname

    !local parameters and variables
    integer, parameter :: UN = 11
    integer :: ierr,i,j

    !..................................................
    open (unit=UN,file=fname,status='replace',iostat=ierr)
    if (ierr.ne.0) then
       write (STDERR,*) 'Error opening file', fname, '; ierr=', ierr
       stop 1
    end if

    do j = 1, u%ny
       do i = 1, u%nx
          write(UN,fmt='(1x,3(e15.7,1x))') u%xmin+(i-1)*u%hx, u%ymin+(j-1)*u%hy, u%u(i,j)
       end do
    end do

    close(UN)
  end subroutine p2d_write_xyz

  !==================================================
  ! write p2d_t structure to a Surfer 6 binary grid file
  subroutine p2d_write_surfer6bin(u,fname)
    implicit none
    type(p2d_t), intent(in) :: u
    character(*), intent(in) :: fname

    !local parameters and variables
    integer, parameter :: UN = 11
    integer :: ierr,i,j

    ! surfer 6 data
    character(4), parameter :: DSBB = 'DSBB'
    integer(2) :: nx,ny
    real(8) :: xmin,xmax,ymin,ymax,zmin,zmax,t
    real(4) :: f

    !..................................................
!    open (unit=UN,file=fname,form='unformatted',status='replace',iostat=ierr)
    open (unit=UN,file=fname,form='binary',status='replace',iostat=ierr)
    if (ierr.ne.0) then
       write (STDERR,*) 'Error opening file', fname, '; ierr=', ierr
       stop 1
    end if
    nx = u%nx
    ny = u%ny
    xmin = u%xmin
    xmax = u%xmax
    ymin = u%ymin
    ymax = u%ymax
    zmin = huge(zmin)
    zmax = -huge(zmax)
    do j = 1, u%ny
       do i = 1, u%nx
          t = u%u(i,j)
          if (t > zmax) zmax = t
          if (t < zmin) zmin = t
       end do
    end do

    write (UN) DSBB, nx, ny
    write (UN) xmin, xmax, ymin, ymax, zmin, zmax
    do j = 1, ny
       do i = 1, nx
          f = u%u(i,j)
          write(UN) f
       end do
    end do
    close(UN)
  end subroutine p2d_write_surfer6bin

  !==================================================
  !interpolate potential to another mesh
  !TODO support arbitrary nx ny
  subroutine p2d_remesh(u,nx,ny)
    implicit none
    type(p2d_t), intent(inout) :: u
    integer, intent(in) :: nx, ny
    !local variables
    type(p2d_t) :: z
    integer :: i,j
    !..................................................
    if (nx .ne. 2*u%nx-1 .or. ny .ne. 2*u%ny-1) then
       write(STDERR,*) 'Unsupported parameters'
       stop 1
    end if
    z = u                       !duplicate pointer
    nullify(u%u)                !nullify ptr to old and
    call p2d_alloc(u,nx,ny)     !allocate new
    do j = 1, z%ny
       do i = 1, z%nx
          u%u(2*i-1,2*j-1) = z%u(i,j)
          if (i < z%nx) u%u(2*i,2*j-1) = 0.5*( z%u(i,j) + z%u(i+1,j) )
          if (j < z%ny) u%u(2*i-1,2*j) = 0.5*( z%u(i,j) + z%u(i,j+1) )
          if (i < z%nx .and. j < z%ny) u%u(2*i,2*j) = &
               & 0.25*( z%u(i,j) + z%u(i+1,j) + z%u(i,j+1) + z%u(i+1,j+1) )
       end do
    end do
    deallocate(z%u)             !forget old
  end subroutine p2d_remesh

  !==================================================
  !transpose potential
  subroutine p2d_transpose(u)
    implicit none
    type(p2d_t), intent(inout) :: u

    !local variables
    type(p2d_t) :: z

    !..................................................
    z = u                       !pointer duplicates
    nullify(u%u)                !nullify it (% to please Visual F90)
    call p2d_alloc(u,z%ny,z%nx) !and allocate new data
    u%u(1:u%nx,1:u%ny) = transpose(z%u(1:z%nx,1:z%ny))
    deallocate(z%u)             !forget old data
    u%hx = z%hy
    u%hy = z%hx
    u%xmin = z%ymin
    u%xmax = z%ymax
    u%ymin = z%xmin
    u%ymax = z%xmax
  end subroutine p2d_transpose

  !==================================================
  !flip potential
  subroutine p2d_flipx(u)
    implicit none
    type(p2d_t), intent(inout) :: u
    integer :: i,j,nx
    real(WP) :: t
    nx = u%nx
    do j = 1, u%ny
       do i = 1, nx/2
          t = u%u(i,j)
          u%u(i,j) = u%u(nx+1-i,j)
          u%u(nx+1-i,j) = t
       end do
    end do
  end subroutine p2d_flipx
  subroutine p2d_flipy(u)
    implicit none
    type(p2d_t), intent(inout) :: u
    integer :: j,nx,ny
    real(WP) :: t(u%nx)
    nx = u%nx
    ny = u%ny
    do j = 1, ny/2
       t(1:nx) = u%u(1:nx,j)
       u%u(1:nx,j) = u%u(1:nx,ny+1-j)
       u%u(1:nx,ny+1-j) = t(1:nx)
    end do
  end subroutine p2d_flipy

end module potential2d
