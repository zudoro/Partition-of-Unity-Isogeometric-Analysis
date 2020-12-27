! NAME : HYUNJU KIM
! DATE : 11 / 8 / 2010

MODULE GLBVAR

	USE NEWTYPE

  IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   GLOBAL PARAMETER

	REAL(8), PARAMETER :: PI=DACOS(-1.0D0)
	REAL(8), PARAMETER :: DEGREE = PI/180.0D0
	REAL(8), PARAMETER:: EPS = 5.0D-15
    
    ! PU-IGA boundary layer parameters
    !-----------------------------------------------------------!
	REAL*8, PARAMETER :: EPSLN = 0.00001d0
	REAL*8, PARAMETER :: SIGMA = 0.650D0
	
	REAL*8, PARAMETER :: ALPHA = 0.0040D0
	REAL*8, PARAMETER :: BETA = 0.50D0	
	
	REAL*8, PARAMETER :: KAPPA = 0.50D0
	!-----------------------------------------------------------!
	
	! PU with modified B-splines with mapping methods for fourth order PDE
	REAL*8, PARAMETER :: RLMBD = 1.0D0 ! lambda = singular intensity factor
	REAL*8, PARAMETER :: ANGL_FTR = 3.0d0/2.0d0 ! angle factor = factor of angle in the circular singular zone 

    !! ELPT = Elliptic, CNVD = Convection Diffusion, FOUR = Fourth order PDE
	CHARACTER(LEN=4), PARAMETER :: PDE='ELPT'
	INTEGER, PARAMETER :: PROBLEM=4
    REAL*8, DIMENSION(1:2), PARAMETER :: RADI_SING = (/0.30D0, 0.50D0/)
    
    REAL*8, PARAMETER :: DELTA = (RADI_SING(2) - RADI_SING(1)) ! delta in the physical domain
    
    REAL*8, DIMENSION(1:2), PARAMETER :: INV_DELTA = (/ (RADI_SING(1)/RADI_SING(2))**(1.0D0/RLMBD), DELTA /)     ! (PDE=='ELPT', PROBLEM==3) OR (PDE=='FOUR', PROBLEM==5)
    
    !! Do you want to use a modified B-spline basis functions?
! 	CHARACTER(LEN=1), PARAMETER  :: MODIFIED_BS='Y'    ! 'Y' = Yes,    'N' = No,    'H' = Sort of...means hybrid of B-splines and enrichment functions
	
    INTEGER, PARAMETER :: MODFD_BSORDER = 1
    
!     REAL*8, DIMENSION(1:MODFD_BSORDER+1), PARAMETER :: EXPO_ARRAY = (/2.0D0*1.5440D0, 2.0D0*1.5440D0 + 4.0D0, 2.0D0*1.5440D0 + 8.0D0, 2.0D0*1.5440D0 + 12.0D0, 2.0D0*1.5440D0 + 16.0d0, 4.0d0, 6.0d0, 8.0d0, 10.0d0, 12.0d0, 14.0d0, 16.0d0/)    
!     REAL*8, DIMENSION(1:MODFD_BSORDER+1), PARAMETER :: EXPO_ARRAY = (/2.0D0*1.5440D0, 2.0D0*1.5440D0 + 4.0D0, 2.0D0*1.5440D0 + 8.0D0, 2.0D0*1.5440D0 + 12.0D0, 2.0D0*1.5440D0 + 16.0d0/)    
    
!     REAL*8, DIMENSION(1:MODFD_BSORDER+1), PARAMETER :: EXPO_ARRAY = (/1.0D0, 3.0D0, 4.0D0, 5.0d0, 6.0d0, 7.0d0, 8.0D0, 9.0D0/)
    REAL*8, DIMENSION(1:MODFD_BSORDER+1), PARAMETER :: EXPO_ARRAY = (/1.0D0, 3.0D0/)    
    
	INTEGER, PARAMETER :: NUMPATCH = 2
	INTEGER, PARAMETER  :: PUORDER = 3
	
	!! Do you solve the boundary layer problem?
	CHARACTER(LEN=1) :: BL

	!! Do you want to use enriching functions?
	CHARACTER(LEN=1) :: USE_ENRICH

	!! Do you want to use Shishkin type refinement?
	CHARACTER(LEN=2) :: MESHTYPE
	
	!! Do you impose a homogeneous boundary condition?
	CHARACTER(LEN=1) :: ZEROBD
	
! 	REAL*8, ALLOCATABLE :: PATCHBDX(:,:), PATCHBDY(:,:)
	
	!! B-SPLINE GLOBAL VARIABLES
	INTEGER :: BS_ORDER(NUMPATCH, 2), EXTRA_KNOTS(NUMPATCH, 2), NURBS_REGULARITY(NUMPATCH, 2)
	
	INTEGER :: NUMBS(2), LOC_NUMBS(NUMPATCH,2), DOF, BD_DOF
	TYPE(BD_INFO) :: BDNDX(MAX_LENGTH), ZERO_BDNDX(MAX_LENGTH)
	
	!!  VARIABLES FOR GEOMETRIC MAP
	TYPE(KNOT_VECTOR) :: BASIS_KVEC(NUMPATCH, 2), OLD_BASIS_KVEC(2), GEO_KVEC(2)
	TYPE(CONTROL_POINTS_2D) :: BASIS_CTL(NUMPATCH), GEO_CTL
	
		!!  CONNECTIVITY ARRAY
! 	INTEGER, ALLOCATABLE :: CARRAY(:,:)
	
	!!  BINOMIAL COEFFICIENTS
	REAL(8) :: BINOM(MAX_DIFF_ORDER,0:MAX_DIFF_ORDER)

	! INTEGRAL GLOBAL VARIABLES
	INTEGER :: NUMGSPT
	REAL*8 :: IR_GRID(NUMPATCH, NUMPATCH, 2, MAX_LENGTH), TMPIR_GRID(NUMPATCH, NUMPATCH, 2, MAX_LENGTH), INVIR_GRID(NUMPATCH, NUMPATCH, 2, MAX_LENGTH)
	REAL*8 :: L2_GRID(2, MAX_LENGTH)
	INTEGER :: NUMIR(NUMPATCH, NUMPATCH, 2), TMP_NUMIR(NUMPATCH, NUMPATCH, 2), L2_NUMIR(2), INV_NUMIR(NUMPATCH, NUMPATCH, 2)

	REAL*8, ALLOCATABLE :: GSX(:), GSXW(:), GSY(:), GSYW(:)
	
	INTEGER, ALLOCATABLE :: NDX(:,:)

	CHARACTER(LEN=2) :: CHAR_ORDER(2)
	CHARACTER(LEN=2) :: CHAR_BASIS(2)
	CHARACTER(LEN=5) :: CHAR_KNOT_INTERVAL(2)
	CHARACTER(LEN=2) :: CHAR_SPE_CONST
	CHARACTER(LEN=2) :: CHAR_GEO_CONST
	CHARACTER(LEN=7) :: CHAR_PROBLEM
	CHARACTER(LEN=100) :: FILENAME
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	MATERIAL COEFFICIENTS

	INTEGER, DIMENSION(1:2), PARAMETER :: ORDER1=(/ 1, 2 /)
	INTEGER, DIMENSION(1:2), PARAMETER :: ORDER2=(/ 2, 1 /)

	TYPE(POINT2D), PARAMETER :: PT0 = POINT2D(0.D0,0.D0)

END MODULE GLBVAR
