! version %I% of %G%
module params
  implicit none

  !======================================================================
  ! constants
  !======================================================================

  ! maximum length of strings
  integer, parameter :: MAXSLEN = 512

  ! working precision
  integer, parameter :: WP = kind(1d0)
  !real(WP), parameter :: SQRTEPS = epsilon(1.0_WP)**0.5
  real(WP), parameter :: SQRTEPS = epsilon(1.0_WP)*1e5 
  ! xlf90 does not allow **NONINTEGER in initialization

  ! mathematical and physical constants
  real(WP), parameter :: PI = 3.1415926535897932384626433832795028841972_WP

!  real(WP), parameter :: GAM = 0.067! 0.07766_WP !...??
  real(WP), parameter :: GAM = 0.07! 0.07766_WP !...??
  real(WP), parameter :: E00 = 38.1_WP / GAM !544.28
  real(WP), parameter :: HEK = 0.000025_WP ! step in ef for dE/dk
!  real(WP), parameter :: Emin= -6.3_WP, Emax=4.9_WP!-2.5_WP ! range of Ef
!  real(WP), parameter :: Emin= -5.3_WP, Emax=-4.9_WP!-2.5_WP ! range of Ef
!  real(WP), parameter :: Emin= -2.5_WP, Emax=-2.495_WP ! range of Ef
!  real(WP), parameter :: Emin= 4.95_WP, Emax=4.951_WP ! range of Ef
 real(WP), parameter :: Emin= 0.0_WP, Emax=0.005_WP ! range of Ef
!  real(WP), parameter :: Emin= -4.55_WP, Emax=-4.35_WP ! range of Ef
  real(WP), parameter :: B_min=1.0_WP, B_max=2.001_WP ! range of Btesla
!  real(WP), parameter :: B_min=-2.0_WP, B_max=-2.0_WP ! range of Btesla
  real(WP), parameter :: deltaE= 0.01_WP, deltaB=0.0025_WP ! step by Ef,Btesla

  ! standard units
  integer, parameter :: STDIN = 5 !unit 5 is stdin
  integer, parameter :: STDOUT = 6 !unit 6 is stdout
  integer, parameter :: STDERR = 0 !unit 0 is stderr

  character(*) DATAFILE
  parameter (DATAFILE='GofBE0_imp_a.dat')
!  parameter (DATAFILE='GofB0E0i01_uni_a.dat')
  character(*) DATFILE
  parameter (DATFILE='CurrAbs.dat')
  integer, parameter :: DTOUT = 5 ! logical device
  integer, parameter :: DTOUT1 = 3 ! logical device
  integer, parameter :: LAPACK = 1
  integer, parameter :: ESSL = 0
  integer, parameter :: USED_LIB = LAPACK
  integer, parameter :: LDVL = 1
  character, parameter :: jobvl = 'N'
  character, parameter :: jobvr = 'V'

  !======================================================================
  ! interfaces
  !======================================================================
  interface
     real(kind=kind(1d0)) function dznrm2(n,x,incx)
       integer, intent(in) :: n, incx
       complex(kind=kind(1d0)), intent(in) :: x(*)
     end function dznrm2

      subroutine zgeev( jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info )
      character, intent(in) ::  jobvl, jobvr
      integer, intent(out) ::  info
      integer, intent(in) :: lda, ldvl, ldvr, n, lwork
      real(kind=kind(1d0)), intent(inout) ::  rwork( * )
      complex(kind=kind(1d0)), intent(inout) ::  a( lda, * ), work( * )
      complex(kind=kind(1d0)), intent(out) ::  w( * ), vl( ldvl, * ), vr( ldvr, * )
      end subroutine zgeev


  end interface

end module params
