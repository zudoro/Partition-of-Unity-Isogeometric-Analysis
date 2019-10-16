MODULE GEOMETRY

	USE PATCH_MAPPING

! 	IMPLICIT INTEGER (I-N)
! 	IMPLICIT REAL(8) (A-H,O-Z)
    IMPLICIT NONE
CONTAINS

SUBROUTINE GET_GEO()
	
	INTEGER :: I, J, II, JJ, K, KK, DIR, PATCH, PATCH_ROW, PATCH_COLUMN, PATCHES(2)
	INTEGER :: BASIS_POLY_ORDER, OLD_BASIS_POLY_ORDER, BASIS_NUM_KNOTS, OLD_BASIS_NUM_KNOTS
	INTEGER :: BASIS_MULTIPLICITIES(MAX_LENGTH), OLD_BASIS_MULTIPLICITIES(MAX_LENGTH)
	REAL(8) :: BASIS_KNOTS(MAX_LENGTH), OLD_BASIS_KNOTS(MAX_LENGTH), NEW_KNOTS(0:MAX_LENGTH), TMP_KNOTS(0:MAX_LENGTH)
	CHARACTER(LEN=1) :: AXIS
	LOGICAL, ALLOCATABLE :: BD_MASK(:)
	INTEGER :: EGM_ETA_ORDER
	REAL*8 :: INNER_RADIUS, INNER_LAYER
	
	TYPE(POINT2D) :: PHYPT, PARPT
	INTEGER :: TMPK
	
! 1ST QUARTER - CIRCLE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	GEO_KVEC(1)%POLY_ORDER = 2
! 	GEO_KVEC(2)%POLY_ORDER = 2
! 
! 	GEO_KVEC(1)%LENGTH = 5
! 	GEO_KVEC(2)%LENGTH = 5
! 
! 	GEO_KVEC(1)%KNOTS(0:GEO_KVEC(1)%LENGTH) = (/0.D0, 0.D0, 0.D0, 1.D0, 1.D0, 1.D0/)
! 	GEO_KVEC(2)%KNOTS(0:GEO_KVEC(2)%LENGTH) = (/0.D0, 0.D0, 0.D0, 1.D0, 1.D0, 1.D0/)
! 	
! ! SET CONTROL POINTS AND WEIGTH
! 	
! 	GEO_CTL%D=2
! 	
! 	GEO_CTL%PTS(0,0,0) = POINT2D(0.D0,0.D0)
! 	GEO_CTL%PTS(1,0,0) = POINT2D(0.5D0,0.D0)
! 	GEO_CTL%PTS(2,0,0) = POINT2D(1.D0,0.D0)
! 	
! 	GEO_CTL%PTS(0,1,0) = POINT2D(0.D0,0.5D0)
! 	GEO_CTL%PTS(1,1,0) = POINT2D(0.5D0*DCOS(45.D0*DEGREE),0.5D0*DSIN(45.D0*DEGREE))
! 	GEO_CTL%PTS(2,1,0) = POINT2D(1.D0,-1.D0 + DCOS(45.D0*DEGREE) + DSIN(45.D0*DEGREE))
! 	
! 	GEO_CTL%PTS(0,2,0) = POINT2D(0.D0,1.D0)
! 	GEO_CTL%PTS(1,2,0) = POINT2D(GEO_CTL%PTS(2,1,0)%Y, GEO_CTL%PTS(2,1,0)%X)
! 	GEO_CTL%PTS(2,2,0) = POINT2D(DCOS(45.D0*DEGREE), DSIN(45.D0*DEGREE))
! 	
! 	GEO_CTL%WGTS(0,0,0) = 1.D0
! 	GEO_CTL%WGTS(1,0,0) = 1.D0
! 	GEO_CTL%WGTS(2,0,0) = 1.D0
! 	
! 	GEO_CTL%WGTS(0,1,0) = 1.D0
! 	GEO_CTL%WGTS(1,1,0) = 1.D0
! 	GEO_CTL%WGTS(2,1,0) = DCOS(22.5D0*DEGREE)
! 	
! 	GEO_CTL%WGTS(0,2,0) = 1.D0
! 	GEO_CTL%WGTS(1,2,0) = DCOS(22.5D0*DEGREE)
! 	GEO_CTL%WGTS(2,2,0) = 1.D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 2ND QUARTER - CIRCLE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	GEO_KVEC(1)%POLY_ORDER = 2
! 	GEO_KVEC(2)%POLY_ORDER = 2
! 
! 	GEO_KVEC(1)%LENGTH = 5
! 	GEO_KVEC(2)%LENGTH = 5
! 
! 	GEO_KVEC(1)%KNOTS(0:GEO_KVEC(1)%LENGTH) = (/0.D0, 0.D0, 0.D0, 1.D0, 1.D0, 1.D0/)
! 	GEO_KVEC(2)%KNOTS(0:GEO_KVEC(2)%LENGTH) = (/0.D0, 0.D0, 0.D0, 1.D0, 1.D0, 1.D0/)
! 	
! ! SET CONTROL POINTS AND WEIGTH
! 	
! 	GEO_CTL%D=2
! 	
! 	GEO_CTL%PTS(0,0,0) = POINT2D(-1.D0,0.D0)
! 	GEO_CTL%PTS(1,0,0) = POINT2D(-0.5D0,0.D0)
! 	GEO_CTL%PTS(2,0,0) = POINT2D(0.D0,0.D0)
! 	
! 	GEO_CTL%PTS(0,1,0) = POINT2D(-1.D0,-1.D0 + DCOS(45.D0*DEGREE) + DSIN(45.D0*DEGREE))
! 	GEO_CTL%PTS(1,1,0) = POINT2D(-0.5D0*DCOS(45.D0*DEGREE),0.5D0*DSIN(45.D0*DEGREE))
! 	GEO_CTL%PTS(2,1,0) = POINT2D(0.D0,0.5D0)
! 	
! 	GEO_CTL%PTS(0,2,0) = POINT2D(-DCOS(45.D0*DEGREE), DSIN(45.D0*DEGREE))
! 	GEO_CTL%PTS(1,2,0) = POINT2D(-(-1.D0 + DCOS(45.D0*DEGREE) + DSIN(45.D0*DEGREE)), 1.D0)
! 	GEO_CTL%PTS(2,2,0) = POINT2D(0.D0,1.D0)
! 	
! 	GEO_CTL%WGTS(0,0,0) = 1.D0
! 	GEO_CTL%WGTS(1,0,0) = 1.D0
! 	GEO_CTL%WGTS(2,0,0) = 1.D0
! 	
! 	GEO_CTL%WGTS(0,1,0) = DCOS(22.5D0*DEGREE)
! 	GEO_CTL%WGTS(1,1,0) = 1.D0
! 	GEO_CTL%WGTS(2,1,0) = 1.D0
! 	
! 	GEO_CTL%WGTS(0,2,0) = 1.D0
! 	GEO_CTL%WGTS(1,2,0) = DCOS(22.5D0*DEGREE)
! 	GEO_CTL%WGTS(2,2,0) = 1.D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! UNIT CIRCLE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	GEO_KVEC(1)%POLY_ORDER = 2
! 	GEO_KVEC(2)%POLY_ORDER = 2
! 
! 	GEO_KVEC(1)%LENGTH = 5
! 	GEO_KVEC(2)%LENGTH = 5
! 
! 	GEO_KVEC(1)%KNOTS(0:GEO_KVEC(1)%LENGTH) = (/0.D0,0.D0,0.D0,1.D0,1.D0,1.D0/)
! 	GEO_KVEC(2)%KNOTS(0:GEO_KVEC(2)%LENGTH) = (/0.D0,0.D0,0.D0,1.D0,1.D0,1.D0/)
! 	
! ! SET CONTROL POINTS AND WEIGTH
! 	GEO_CTL%D=2
! 	GEO_CTL%PTS(0,0,0) = POINT2D(0.D0,-1.D0)
! 	GEO_CTL%PTS(1,0,0) = POINT2D(1.D0,-1.D0)
! 	GEO_CTL%PTS(2,0,0) = POINT2D(1.D0,0.D0)
! 	GEO_CTL%PTS(0,1,0) = POINT2D(-1.D0,-1.D0)
! 	GEO_CTL%PTS(1,1,0) = POINT2D(0.D0,0.D0)
! 	GEO_CTL%PTS(2,1,0) = POINT2D(1.D0,1.D0)
! 	GEO_CTL%PTS(0,2,0) = POINT2D(-1.D0,0.D0)
! 	GEO_CTL%PTS(1,2,0) = POINT2D(-1.D0,1.D0)
! 	GEO_CTL%PTS(2,2,0) = POINT2D(0.D0,1.D0)
! 	
! 	GEO_CTL%WGTS(0,0,0) = 1.D0
! 	GEO_CTL%WGTS(1,0,0) = 1.D0/DSQRT(2.D0)
! 	GEO_CTL%WGTS(2,0,0) = 1.D0
! 	GEO_CTL%WGTS(0,1,0) = 1.D0/DSQRT(2.D0)
! 	GEO_CTL%WGTS(1,1,0) = 1.D0
! 	GEO_CTL%WGTS(2,1,0) = 1.D0/DSQRT(2.D0)
! 	GEO_CTL%WGTS(0,2,0) = 1.D0
! 	GEO_CTL%WGTS(1,2,0) = 1.D0/DSQRT(2.D0)
! 	GEO_CTL%WGTS(2,2,0) = 1.D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! UNIT CIRCLE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! 	GEO_KVEC(1)%POLY_ORDER = 2
! ! 	GEO_KVEC(2)%POLY_ORDER = 2
! ! 
! ! 	GEO_KVEC(1)%LENGTH = 5
! ! 	GEO_KVEC(2)%LENGTH = 7
! ! 
! ! 	GEO_KVEC(1)%KNOTS(0:GEO_KVEC(1)%LENGTH) = (/0.D0, 0.D0, 0.D0, 1.D0, 1.D0, 1.D0/)
! ! 	GEO_KVEC(2)%KNOTS(0:GEO_KVEC(2)%LENGTH) = (/0.D0, 0.D0, 0.D0, 0.5D0, 0.5D0, 1.D0, 1.D0, 1.D0/)
! ! 	
! ! ! SET CONTROL POINTS AND WEIGTH
! ! 	
! ! 	GEO_CTL%D=2
! ! 	
! ! 	GEO_CTL%PTS(0,0,0) = POINT2D(0.D0,-1.D0)
! ! 	GEO_CTL%PTS(1,0,0) = POINT2D(1.D0,-1.D0)
! ! 	GEO_CTL%PTS(2,0,0) = POINT2D(1.D0,0.D0)
! ! 	
! ! 	GEO_CTL%PTS(0,1,0) = POINT2D(-1.D0,-1.D0)
! ! 	GEO_CTL%PTS(1,1,0) = POINT2D(0.5D0,-0.5D0)
! ! 	GEO_CTL%PTS(2,1,0) = POINT2D(0.5D0,0.D0)
! ! 	
! ! 	GEO_CTL%PTS(0,2,0) = POINT2D(-1.D0,0.D0)
! ! 	GEO_CTL%PTS(1,2,0) = POINT2D(-0.5D0,0.D0)
! ! 	GEO_CTL%PTS(2,2,0) = POINT2D(0.D0,0.D0)
! ! 	
! ! 	GEO_CTL%PTS(0,3,0) = POINT2D(-1.D0,1.D0)
! ! 	GEO_CTL%PTS(1,3,0) = POINT2D(0.5D0,0.5D0)
! ! 	GEO_CTL%PTS(2,3,0) = POINT2D(0.5D0,0.D0)
! ! 	
! ! 	GEO_CTL%PTS(0,4,0) = POINT2D(0.D0,1.D0)
! ! 	GEO_CTL%PTS(1,4,0) = POINT2D(1.D0,1.D0)
! ! 	GEO_CTL%PTS(2,4,0) = POINT2D(1.D0,0.D0)
! ! 	
! ! 	GEO_CTL%WGTS(0,0,0) = 1.D0
! ! 	GEO_CTL%WGTS(1,0,0) = DCOS(45.D0*DEGREE)
! ! 	GEO_CTL%WGTS(2,0,0) = 1.D0
! ! 	
! ! 	GEO_CTL%WGTS(0,1,0) = DCOS(45.D0*DEGREE)
! ! 	GEO_CTL%WGTS(1,1,0) = DCOS(45.D0*DEGREE)
! ! 	GEO_CTL%WGTS(2,1,0) = 1.D0
! ! 	
! ! 	GEO_CTL%WGTS(0,2,0) = 1.D0
! ! 	GEO_CTL%WGTS(1,2,0) = 1.D0
! ! 	GEO_CTL%WGTS(2,2,0) = 1.D0
! ! 	
! ! 	GEO_CTL%WGTS(0,3,0) = DCOS(45.D0*DEGREE)
! ! 	GEO_CTL%WGTS(1,3,0) = DCOS(45.D0*DEGREE)
! ! 	GEO_CTL%WGTS(2,3,0) = 1.D0
! ! 	
! ! 	GEO_CTL%WGTS(0,4,0) = 1.D0
! ! 	GEO_CTL%WGTS(1,4,0) = DCOS(45.D0*DEGREE)
! ! 	GEO_CTL%WGTS(2,4,0) = 1.D0
! 	
! 	BASIS_CTL%WGTS(:,:,:) = 1.D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! UNIT CIRCLE - PINCH BOTTOM LINE AT THE ORIGIN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	GEO_KVEC(1)%POLY_ORDER = 2
! 	GEO_KVEC(2)%POLY_ORDER = 2
! 
! 	GEO_KVEC(1)%LENGTH = 11
! 	GEO_KVEC(2)%LENGTH = 5
! 
! 	GEO_KVEC(1)%KNOTS(0:GEO_KVEC(1)%LENGTH) = (/0.D0, 0.D0, 0.D0, 0.25D0, 0.25D0, 0.5D0, 0.5D0, 0.75D0, 0.75D0, 1.D0, 1.D0, 1.D0/)
! 	GEO_KVEC(2)%KNOTS(0:GEO_KVEC(2)%LENGTH) = (/0.D0, 0.D0, 0.D0, 1.D0, 1.D0, 1.D0/)
! 	
! ! SET CONTROL POINTS AND WEIGTH
! 	GEO_CTL%D=2
! 	
! 	GEO_CTL%PTS(0:8,0,0) = POINT2D(0.D0,0.D0)
! 
! 	GEO_CTL%PTS(0,2,0) = POINT2D(1.D0,0.D0)
! 	GEO_CTL%PTS(1,2,0) = POINT2D(1.D0,-1.D0)
! 	GEO_CTL%PTS(2,2,0) = POINT2D(0.D0,-1.D0)
! 	GEO_CTL%PTS(3,2,0) = POINT2D(-1.D0,-1.D0)
! 	GEO_CTL%PTS(4,2,0) = POINT2D(-1.D0,0.D0)
! 	GEO_CTL%PTS(5,2,0) = POINT2D(-1.D0,1.D0)
! 	GEO_CTL%PTS(6,2,0) = POINT2D(0.D0,1.D0)
! 	GEO_CTL%PTS(7,2,0) = POINT2D(1.D0,1.D0)
! 	GEO_CTL%PTS(8,2,0) = POINT2D(1.D0,0.D0)
! 	
! 	DO I = 0, 8
! 		GEO_CTL%PTS(I,1,0) = 0.5D0*(GEO_CTL%PTS(I,0,0) + GEO_CTL%PTS(I,2,0))
! 	ENDDO
! 	
! 	DO I=0, 2
! 		GEO_CTL%WGTS(0:8,I,0) = (/1.D0, DCOS(45.D0*DEGREE), 1.D0, DCOS(45.D0*DEGREE), 1.D0, DCOS(45.D0*DEGREE), 1.D0, DCOS(45.D0*DEGREE), 1.D0/)
! 	ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! FIGURE 2.15 IN THE IGA BOOK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	GEO_KVEC(1)%POLY_ORDER = 2
! 	GEO_KVEC(2)%POLY_ORDER = 2
! 
! 	GEO_KVEC(1)%LENGTH = 6
! 	GEO_KVEC(2)%LENGTH = 5
! 
! 	GEO_KVEC(1)%KNOTS(0:GEO_KVEC(1)%LENGTH) = (/0.D0,0.D0,0.D0,0.5D0,1.D0,1.D0,1.D0/)
! 	GEO_KVEC(2)%KNOTS(0:GEO_KVEC(2)%LENGTH) = (/0.D0,0.D0,0.D0,1.D0,1.D0,1.D0/)
! 	
! ! SET CONTROL POINTS AND WEIGTH
! 	GEO_CTL%D=2
! 	GEO_CTL%PTS(0,0,0) = POINT2D(2.D0/2.D0,0.D0)
! 	GEO_CTL%PTS(0,1,0) = POINT2D(1.D0/2.D0,0.D0)
! 	GEO_CTL%PTS(0,2,0) = POINT2D(0.D0,0.D0)
! 	GEO_CTL%PTS(1,0,0) = POINT2D(2.D0/2.D0,1.D0/5.D0)
! 	GEO_CTL%PTS(1,1,0) = POINT2D(1.D0/2.D0,2.D0/5.D0)
! 	GEO_CTL%PTS(1,2,0) = POINT2D(0.D0/2.D0,2.D0/5.D0)
! 	GEO_CTL%PTS(2,0,0) = POINT2D(6.D0/5.D0,1.5D0/5.D0)
! 	GEO_CTL%PTS(2,1,0) = POINT2D(6.D0/5.D0,4.D0/5.D0)
! 	GEO_CTL%PTS(2,2,0) = POINT2D(6.D0/5.D0,5.D0/5.D0)
! 	GEO_CTL%PTS(3,0,0) = POINT2D(7.D0/5.D0,1.5D0/5.D0)
! 	GEO_CTL%PTS(3,1,0) = POINT2D(7.D0/5.D0,4.D0/5.D0)
! 	GEO_CTL%PTS(3,2,0) = POINT2D(7.D0/5.D0,5.D0/5.D0)
! 	
! 	GEO_CTL%WGTS(:,:,0) = 1.D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! POLYGONAL L-SHAPED DOMAIN WITH ONE PATCH
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	GEO_KVEC(1)%POLY_ORDER = 2
! 	GEO_KVEC(2)%POLY_ORDER = 2
! 
! 	GEO_KVEC(1)%LENGTH = 6
! 	GEO_KVEC(2)%LENGTH = 5
! 
! 	GEO_KVEC(1)%KNOTS(0:GEO_KVEC(1)%LENGTH) = (/0.D0,0.D0,0.D0,0.5D0,1.D0,1.D0,1.D0/)
! 	GEO_KVEC(2)%KNOTS(0:GEO_KVEC(2)%LENGTH) = (/0.D0,0.D0,0.D0,1.D0,1.D0,1.D0/)
! 	
! ! SET CONTROL POINTS AND WEIGTH
! 	GEO_CTL%D=2
! 	GEO_CTL%PTS(0,0,0) = POINT2D(0.D0,-1.D0)
! 	GEO_CTL%PTS(1,0,0) = POINT2D(0.D0,0.D0)
! 	GEO_CTL%PTS(2,0,0) = POINT2D(0.D0,0.D0)
! 	GEO_CTL%PTS(3,0,0) = POINT2D(1.D0,0.D0)
! 	GEO_CTL%PTS(0,1,0) = POINT2D(-0.5D0,-1.D0)
! 	GEO_CTL%PTS(1,1,0) = POINT2D(-2.D0/3.D0,1.D0/3.D0)
! 	GEO_CTL%PTS(2,1,0) = POINT2D(-1.D0/3.D0,2.D0/3.D0)
! 	GEO_CTL%PTS(3,1,0) = POINT2D(1.D0,0.5D0)
! 	GEO_CTL%PTS(0,2,0) = POINT2D(-1.D0,-1.D0)
! 	GEO_CTL%PTS(1,2,0) = POINT2D(-1.D0,1.D0)
! 	GEO_CTL%PTS(2,2,0) = POINT2D(-1.D0,1.D0)
! 	GEO_CTL%PTS(3,2,0) = POINT2D(1.D0,1.D0)
! 	
! 	GEO_CTL%WGTS(:,:,0) = 1.D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! EGM - L-shaped domain with arbitrary order of B-spline along eta-direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 	EGM_ETA_ORDER = 3
! 	INNER_RADIUS = 1.0D0
! 	INNER_LAYER = 0.40D0
! 	
! 	GEO_KVEC(1)%POLY_ORDER = 2
! 	GEO_KVEC(2)%POLY_ORDER = EGM_ETA_ORDER
! 
! 	GEO_KVEC(1)%LENGTH = 9
! 	GEO_KVEC(2)%LENGTH = 3*GEO_KVEC(2)%POLY_ORDER + 1
! 
! 	GEO_KVEC(1)%KNOTS(0:GEO_KVEC(1)%LENGTH) = (/0.0D0, 0.0D0, 0.0D0, 0.4d0, 0.4d0, 0.8d0, 0.8d0, 1.0D0, 1.0D0, 1.0D0/)
! 	
! 	DO I = 0, GEO_KVEC(2)%POLY_ORDER
! 		GEO_KVEC(2)%KNOTS(I) = 0.0D0
! 	ENDDO
! 	DO I = GEO_KVEC(2)%POLY_ORDER + 1, 2*GEO_KVEC(2)%POLY_ORDER
! 		GEO_KVEC(2)%KNOTS(I) = 0.50D0
! 	ENDDO
! 	DO I = 2*GEO_KVEC(2)%POLY_ORDER + 1, 3*GEO_KVEC(2)%POLY_ORDER + 1
! 		GEO_KVEC(2)%KNOTS(I) = 1.0D0
! 	ENDDO
! 	
! ! SET CONTROL POINTS AND WEIGTH
! 	GEO_CTL%D=2
! 	
! 	GEO_CTL%PTS(0:6,0:GEO_KVEC(2)%POLY_ORDER-1,0) = POINT2D(0.0D0,0.0D0)
! 
! 	GEO_CTL%PTS(0,GEO_KVEC(2)%LENGTH - GEO_KVEC(2)%POLY_ORDER - 1,0) = POINT2D(0.0D0,					-INNER_RADIUS)
! 	GEO_CTL%PTS(1,GEO_KVEC(2)%LENGTH - GEO_KVEC(2)%POLY_ORDER - 1,0) = POINT2D(-INNER_RADIUS,-INNER_RADIUS)
! 	GEO_CTL%PTS(2,GEO_KVEC(2)%LENGTH - GEO_KVEC(2)%POLY_ORDER - 1,0) = POINT2D(-INNER_RADIUS, 0.0D0)
! 	GEO_CTL%PTS(3,GEO_KVEC(2)%LENGTH - GEO_KVEC(2)%POLY_ORDER - 1,0) = POINT2D(-INNER_RADIUS, INNER_RADIUS)
! 	GEO_CTL%PTS(4,GEO_KVEC(2)%LENGTH - GEO_KVEC(2)%POLY_ORDER - 1,0) = POINT2D(0.0D0,					 INNER_RADIUS)
! 	GEO_CTL%PTS(5,GEO_KVEC(2)%LENGTH - GEO_KVEC(2)%POLY_ORDER - 1,0) = POINT2D(INNER_RADIUS,	 INNER_RADIUS)
! 	GEO_CTL%PTS(6,GEO_KVEC(2)%LENGTH - GEO_KVEC(2)%POLY_ORDER - 1,0) = POINT2D(INNER_RADIUS,	 0.0D0)
! 
! 	DO I = 0, 6
! 		GEO_CTL%PTS(I,GEO_KVEC(2)%POLY_ORDER,0) = ((INNER_RADIUS - INNER_LAYER)/INNER_RADIUS)*&
! 																										GEO_CTL%PTS(I,GEO_KVEC(2)%LENGTH - GEO_KVEC(2)%POLY_ORDER - 1,0)
! 		DO J = GEO_KVEC(2)%POLY_ORDER + 1, GEO_KVEC(2)%LENGTH - GEO_KVEC(2)%POLY_ORDER - 2
! 			IF (GEO_CTL%PTS(I,GEO_KVEC(2)%LENGTH - GEO_KVEC(2)%POLY_ORDER - 1,0)%X<0.0D0) THEN
! 				GEO_CTL%PTS(I,J,0)%X = GEO_CTL%PTS(I,GEO_KVEC(2)%POLY_ORDER,0)%X - 1.0D0*(J-GEO_KVEC(2)%POLY_ORDER)*&
! 																	DABS(GEO_CTL%PTS(I,GEO_KVEC(2)%LENGTH - GEO_KVEC(2)%POLY_ORDER - 1,0)%X - &
! 																			 GEO_CTL%PTS(I,GEO_KVEC(2)%POLY_ORDER,0)%X)/ &
! 																	(1.0D0*(GEO_KVEC(2)%LENGTH - 2*GEO_KVEC(2)%POLY_ORDER - 1))
! 			ELSE
! 				GEO_CTL%PTS(I,J,0)%X = GEO_CTL%PTS(I,GEO_KVEC(2)%POLY_ORDER,0)%X + 1.0D0*(J-GEO_KVEC(2)%POLY_ORDER)*&
! 																	DABS(GEO_CTL%PTS(I,GEO_KVEC(2)%LENGTH - GEO_KVEC(2)%POLY_ORDER - 1,0)%X - &
! 																			 GEO_CTL%PTS(I,GEO_KVEC(2)%POLY_ORDER,0)%X)/ &
! 																	(1.0D0*(GEO_KVEC(2)%LENGTH - 2*GEO_KVEC(2)%POLY_ORDER - 1))
! 			ENDIF
! 			IF (GEO_CTL%PTS(I,GEO_KVEC(2)%LENGTH - GEO_KVEC(2)%POLY_ORDER - 1,0)%Y<0.0D0) THEN
! 				GEO_CTL%PTS(I,J,0)%Y = GEO_CTL%PTS(I,GEO_KVEC(2)%POLY_ORDER,0)%Y - 1.0D0*(J-GEO_KVEC(2)%POLY_ORDER)*&
! 																	DABS(GEO_CTL%PTS(I,GEO_KVEC(2)%LENGTH - GEO_KVEC(2)%POLY_ORDER - 1,0)%Y - &
! 																			 GEO_CTL%PTS(I,GEO_KVEC(2)%POLY_ORDER,0)%Y)/ &
! 																	(1.0D0*(GEO_KVEC(2)%LENGTH - 2*GEO_KVEC(2)%POLY_ORDER - 1))
! 			ELSE
! 				GEO_CTL%PTS(I,J,0)%Y = GEO_CTL%PTS(I,GEO_KVEC(2)%POLY_ORDER,0)%Y + 1.0D0*(J-GEO_KVEC(2)%POLY_ORDER)*&
! 																	DABS(GEO_CTL%PTS(I,GEO_KVEC(2)%LENGTH - GEO_KVEC(2)%POLY_ORDER - 1,0)%Y - &
! 																			 GEO_CTL%PTS(I,GEO_KVEC(2)%POLY_ORDER,0)%Y)/ &
! 																	(1.0D0*(GEO_KVEC(2)%LENGTH - 2*GEO_KVEC(2)%POLY_ORDER - 1))
! 			ENDIF
! 		ENDDO
! 	ENDDO
! 		
! 	DO I=0, GEO_KVEC(2)%LENGTH - GEO_KVEC(2)%POLY_ORDER - 1
! 		GEO_CTL%WGTS(0:6,I,0) = (/1.0D0, DCOS(45.D0*DEGREE), 1.0D0, DCOS(45.D0*DEGREE), 1.0D0, DCOS(45.D0*DEGREE), 1.0D0/)
! 	ENDDO
	
	!----------------------------------------------------------------------------------------------------------------------------!

	!! Rotate the EGM patches
! 	DO PATCH = NUM_IGA_PATCH + 1, NUM_PATCH
! 		DO I = 0, 6
! 			DO J = 0, GEO_KVEC(2)%LENGTH - GEO_KVEC(2)%POLY_ORDER - 1
! 				GEO_CTL(PATCH)%PTS(I,J,0) = ROTATION(GEO_CTL(PATCH)%PTS(I,J,0), -135.0D0*DEGREE)
! 			ENDDO
! 		ENDDO
! 	ENDDO
	
	!!----------------------------------------------------------------------------------------------------------------------------------------!!
	
	!-----------------------------------------------------------------------------------------------------------
! Elevate degree of B-spline in NURBS mapping
! 	CALL ELEVATE_DEGREE(GEO_KVEC, GEO_CTL, 2, 'X')
! 	CALL ELEVATE_DEGREE(GEO_KVEC, GEO_CTL, 2, 'Y')
!-----------------------------------------------------------------------------------------------------------

! 	CALL REMOVE_KNOT(GEO_KVEC, GEO_CTL, 4, 2, 1, 'X')
! 	CALL REMOVE_KNOT(GEO_KVEC, GEO_CTL, 7, 3, 2, 'X')
	
! 	WRITE(*,*)
! 	WRITE(*,*) '<<< SET GEOMETRIC KNOT, CONTROL POINTS, AND WEIGHT : DONE >>>'
! 	WRITE(*,*)
	
!! ----------------------- [[[[[ SET KNOT INFO OF B-SPLINE BASIS FUNCTIONS IN THE PATCH FOR EACH ]]]]] --------------------------
!-----------------------------------------------------------------------------------------------------------
!! BASIS FUNCTION ON THE OMEGA_HAT1
! 	BASIS_KVEC(1,1,1)%POLY_ORDER = 2
! 	BASIS_KVEC(1,1,2)%POLY_ORDER = 2
! 	
! 	BASIS_KVEC(1,1,1)%LENGTH = 8
! 	BASIS_KVEC(1,1,2)%LENGTH = 5
! 	
! 	BASIS_KVEC(1,1,1)%KNOTS(0:BASIS_KVEC(1,1,1)%LENGTH) = (/0.D0,0.D0,0.D0, 0.50D0, 0.50D0, 0.50D0, 1.D0,1.D0,1.D0/)
! 	BASIS_KVEC(1,1,2)%KNOTS(0:BASIS_KVEC(1,1,2)%LENGTH) = (/0.D0,0.D0,0.D0,1.D0,1.D0,1.D0/)
! 	
! 	BASIS_CTL(1,1)%D = 2
! 	BASIS_CTL(1,1)%PTS(0,0,0) = POINT2D(0.0D0, PATCHBDY(1,1)-DELTA)
! 	BASIS_CTL(1,1)%PTS(1,0,0) = POINT2D(0.20D0, PATCHBDY(1,1)-DELTA)
! 	BASIS_CTL(1,1)%PTS(2,0,0) = POINT2D(0.4D0, PATCHBDY(1,1)-DELTA)
! 	BASIS_CTL(1,1)%PTS(3,0,0) = POINT2D(0.6D0, PATCHBDY(1,1)-DELTA)
! 	BASIS_CTL(1,1)%PTS(4,0,0) = POINT2D(0.8D0, PATCHBDY(1,1)-DELTA)
! 	BASIS_CTL(1,1)%PTS(5,0,0) = POINT2D(1.0D0, PATCHBDY(1,1)-DELTA)
! 	
! 	BASIS_CTL(1,1)%PTS(0,1,0) = POINT2D(0.0D0, 0.50D0*SUM(PATCHBDY(1,:)))
! 	BASIS_CTL(1,1)%PTS(1,1,0) = POINT2D(0.20D0, 0.50D0*SUM(PATCHBDY(1,:)))
! 	BASIS_CTL(1,1)%PTS(2,1,0) = POINT2D(0.40D0, 0.50D0*SUM(PATCHBDY(1,:)))
! 	BASIS_CTL(1,1)%PTS(3,1,0) = POINT2D(0.60D0, 0.50D0*SUM(PATCHBDY(1,:)))
! 	BASIS_CTL(1,1)%PTS(4,1,0) = POINT2D(0.80D0, 0.50D0*SUM(PATCHBDY(1,:)))
! 	BASIS_CTL(1,1)%PTS(5,1,0) = POINT2D(1.0D0, 0.50D0*SUM(PATCHBDY(1,:)))
! 	
! 	BASIS_CTL(1,1)%PTS(0,2,0) = POINT2D(0.0D0, PATCHBDY(1,2))
! 	BASIS_CTL(1,1)%PTS(1,2,0) = POINT2D(0.20D0, PATCHBDY(1,2))
! 	BASIS_CTL(1,1)%PTS(2,2,0) = POINT2D(0.40D0, PATCHBDY(1,2))
! 	BASIS_CTL(1,1)%PTS(3,2,0) = POINT2D(0.60D0, PATCHBDY(1,2))
! 	BASIS_CTL(1,1)%PTS(4,2,0) = POINT2D(0.80D0, PATCHBDY(1,2))
! 	BASIS_CTL(1,1)%PTS(5,2,0) = POINT2D(1.0D0, PATCHBDY(1,2))
! 	
! 
! 	BASIS_CTL(1,1)%WGTS(:,:,0) = 1.D0
! 
! !! BASIS FUNCTION ON THE OMEGA_HAT2
! 	BASIS_KVEC(1,2,1)%POLY_ORDER = 2
! 	BASIS_KVEC(1,2,2)%POLY_ORDER = 2
! 	
! 	BASIS_KVEC(1,2,1)%LENGTH = 5
! 	BASIS_KVEC(1,2,2)%LENGTH = 5
! 	
! 	BASIS_KVEC(1,2,1)%KNOTS(0:BASIS_KVEC(1,2,1)%LENGTH) = (/0.D0,0.D0,0.D0,1.D0,1.D0,1.D0/)
! 	BASIS_KVEC(1,2,2)%KNOTS(0:BASIS_KVEC(1,2,2)%LENGTH) = (/0.D0,0.D0,0.D0,1.D0,1.D0,1.D0/)
! 	
! 	BASIS_CTL(1,2)%D = 2
! 	BASIS_CTL(1,2)%PTS(0,0,0) = POINT2D(0.0D0, PATCHBDY(2,1)-DELTA)
! 	BASIS_CTL(1,2)%PTS(1,0,0) = POINT2D(0.50D0, PATCHBDY(2,1)-DELTA)
! 	BASIS_CTL(1,2)%PTS(2,0,0) = POINT2D(1.0D0, PATCHBDY(2,1)-DELTA)
! 	BASIS_CTL(1,2)%PTS(0,1,0) = POINT2D(0.0D0, 0.50D0*SUM(PATCHBDY(2,:)))
! 	BASIS_CTL(1,2)%PTS(1,1,0) = POINT2D(0.50D0, 0.50D0*SUM(PATCHBDY(2,:)))
! 	BASIS_CTL(1,2)%PTS(2,1,0) = POINT2D(1.0D0, 0.50D0*SUM(PATCHBDY(2,:)))
! 	BASIS_CTL(1,2)%PTS(0,2,0) = POINT2D(0.0D0, PATCHBDY(2,2)+DELTA)
! 	BASIS_CTL(1,2)%PTS(1,2,0) = POINT2D(0.50D0, PATCHBDY(2,2)+DELTA)
! 	BASIS_CTL(1,2)%PTS(2,2,0) = POINT2D(1.0D0, PATCHBDY(2,2)+DELTA)
! 	
! 	BASIS_CTL(1,2)%WGTS(:,:,0) = 1.D0
! 	
! !! BASIS FUNCTION ON THE OMEGA_HAT3
! 	BASIS_KVEC(1,3,1)%POLY_ORDER = 2
! 	BASIS_KVEC(1,3,2)%POLY_ORDER = 2
! 	
! 	BASIS_KVEC(1,3,1)%LENGTH = 7
! ! 	BASIS_KVEC(1,3,1)%LENGTH = 5
! 	BASIS_KVEC(1,3,2)%LENGTH = 5
! 	
! 	BASIS_KVEC(1,3,1)%KNOTS(0:BASIS_KVEC(1,3,1)%LENGTH) = (/0.D0,0.D0,0.D0,0.50D0,0.50D0,1.D0,1.D0,1.D0/)
! ! 	BASIS_KVEC(1,3,1)%KNOTS(0:BASIS_KVEC(1,3,1)%LENGTH) = (/0.D0,0.D0,0.D0,1.D0,1.D0,1.D0/)
! 	BASIS_KVEC(1,3,2)%KNOTS(0:BASIS_KVEC(1,3,2)%LENGTH) = (/0.D0,0.D0,0.D0,1.D0,1.D0,1.D0/)
! 	
! 	BASIS_CTL(1,3)%D = 2
! 	BASIS_CTL(1,3)%PTS(0,0,0) = POINT2D(0.0D0,0.0D0)
! 	BASIS_CTL(1,3)%PTS(1,0,0) = POINT2D(0.250D0,0.0D0)
! 	BASIS_CTL(1,3)%PTS(2,0,0) = POINT2D(0.50D0,0.0D0)
! 	BASIS_CTL(1,3)%PTS(3,0,0) = POINT2D(0.750D0,0.0D0)
! 	BASIS_CTL(1,3)%PTS(4,0,0) = POINT2D(1.0D0,0.0D0)
! 	BASIS_CTL(1,3)%PTS(0,1,0) = POINT2D(0.0D0,0.50D0*SUM(PATCHBDY(3,:)))
! 	BASIS_CTL(1,3)%PTS(1,1,0) = POINT2D(0.250D0,0.50D0*SUM(PATCHBDY(3,:)))
! 	BASIS_CTL(1,3)%PTS(2,1,0) = POINT2D(0.50D0,0.50D0*SUM(PATCHBDY(3,:)))
! 	BASIS_CTL(1,3)%PTS(3,1,0) = POINT2D(0.750D0,0.50D0*SUM(PATCHBDY(3,:)))
! 	BASIS_CTL(1,3)%PTS(4,1,0) = POINT2D(1.0D0,0.50D0*SUM(PATCHBDY(3,:)))
! 	BASIS_CTL(1,3)%PTS(0,2,0) = POINT2D(0.0D0,PATCHBDY(3,2))
! 	BASIS_CTL(1,3)%PTS(1,2,0) = POINT2D(0.250D0,PATCHBDY(3,2))
! 	BASIS_CTL(1,3)%PTS(2,2,0) = POINT2D(0.50D0,PATCHBDY(3,2))
! 	BASIS_CTL(1,3)%PTS(3,2,0) = POINT2D(0.750D0,PATCHBDY(3,2))
! 	BASIS_CTL(1,3)%PTS(4,2,0) = POINT2D(1.0D0,PATCHBDY(3,2))
! 	
! 	BASIS_CTL(1,3)%WGTS(:,:,0) = 1.D0
!-----------------------------------------------------------------------------------------------------------


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! [[[[[ REFINEMENT ]]]]] !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!------------------------------------- [[[[[ B-SPLINE ]]]]] -----------------------------------------------
!-----------------------------------------------------------------------------------------------------------
! ELEVATE DEGREE OF B-SPLINE BASIS FUNCTION
! 	DO PATCH_ROW = 1, NUMPATCH(1)
! 		DO PATCH_COLUMN = 1, NUMPATCH(2)
! 			IF (PATCH_ROW==1 .AND. PATCH_COLUMN==3) THEN ! OMEGA_HAT3 - SHOULD HAVE AN INTERPOLANT AT 0.5 IN THE XI-DIRECTIOIN
! 				DO K = 1, 2
! 					IF (BS_ORDER(PATCH_ROW,PATCH_COLUMN,K)>BASIS_KVEC(PATCH_ROW,PATCH_COLUMN,K)%POLY_ORDER) THEN
! 						DO I = 1, BS_ORDER(PATCH_ROW,PATCH_COLUMN,K) - BASIS_KVEC(PATCH_ROW,PATCH_COLUMN,K)%POLY_ORDER
! 							BASIS_KVEC(PATCH_ROW,PATCH_COLUMN,K) = DEGREE_ELEVATION(BASIS_KVEC(PATCH_ROW,PATCH_COLUMN,K))
! 						ENDDO
! 					ENDIF
! 				ENDDO
! 			ELSE
! 				DO K = 1, 2
! 					IF (BS_ORDER(PATCH_ROW,PATCH_COLUMN,K)>BASIS_KVEC(PATCH_ROW,PATCH_COLUMN,K)%POLY_ORDER) THEN
! 						DO I = 1, BS_ORDER(PATCH_ROW,PATCH_COLUMN,K) - BASIS_KVEC(PATCH_ROW,PATCH_COLUMN,K)%POLY_ORDER
! 							BASIS_KVEC(PATCH_ROW,PATCH_COLUMN,K) = DEGREE_ELEVATION_WITH_CONTROL_REGULARITY(BASIS_KVEC(PATCH_ROW,PATCH_COLUMN,K), BASIS_KVEC(PATCH_ROW,PATCH_COLUMN,K)%POLY_ORDER)
! 						ENDDO
! 					ENDIF
! 				ENDDO
! 			ENDIF
! 		ENDDO
! 	ENDDO
!-----------------------------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------------------------
! INSERT NEW KNOTS IN THE KNOTS CORRESPONDING TO B-SPLINE BASIS FUNCTIONS
! 	DO PATCH_ROW = 1, NUMPATCH(1)
! 		DO PATCH_COLUMN = 1, NUMPATCH(2)
! 			DO K = 1, 2
! 				IF (EXTRA_KNOTS(PATCH_ROW,PATCH_COLUMN,K)>0) THEN
! 					BASIS_KVEC(PATCH_ROW,PATCH_COLUMN,K) = UNIFORM_KNOT_INSERTION(BASIS_KVEC(PATCH_ROW,PATCH_COLUMN,K),EXTRA_KNOTS(PATCH_ROW,PATCH_COLUMN,K),1)
! 				ENDIF
! 			ENDDO
! 		ENDDO
! 	ENDDO
!-----------------------------------------------------------------------------------------------------------


!!------------------------------------- [[[[[ NURBS ]]]]] -----------------------------------------------
!-----------------------------------------------------------------------------------------------------------
! ELEVATE DEGREE OF BERNSTEIN POLY. AND FIND NEW CONTROL POINTS

! 	IF (BS_ORDER(1).GT.GEO_KVEC(1)%POLY_ORDER) THEN
! 		DO I=1, BS_ORDER(1) - GEO_KVEC(1)%POLY_ORDER
! 			CALL P_REFINEMENT(BASIS_CTL, BASIS_KVEC(1:2)%POLY_ORDER, 'X')
! 			BASIS_KVEC(1) = DEGREE_ELEVATION(BASIS_KVEC(1))
! 		ENDDO
! 	ENDIF
! 
! 	IF (BS_ORDER(2).GT.GEO_KVEC(2)%POLY_ORDER) THEN
! 		DO I=1, BS_ORDER(2) - GEO_KVEC(2)%POLY_ORDER
! 			CALL P_REFINEMENT(BASIS_CTL, BASIS_KVEC(1:2)%POLY_ORDER, 'Y')
! 			BASIS_KVEC(2) = DEGREE_ELEVATION(BASIS_KVEC(2))
! 		ENDDO
! 	ENDIF
!-----------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------
! ELEVATE DEGREE OF NURBS BASIS FUNCTIONS
! 	DO PATCH_COLUMN = 1, NUMPATCH(2)
! 		IF (BS_ORDER(1,PATCH_COLUMN,1).GT.BASIS_KVEC(1,PATCH_COLUMN,1)%POLY_ORDER) THEN
! 			CALL ELEVATE_DEGREE(BASIS_KVEC(1,PATCH_COLUMN,:), BASIS_CTL(1,PATCH_COLUMN), BS_ORDER(1,PATCH_COLUMN,1)-BASIS_KVEC(1,PATCH_COLUMN,1)%POLY_ORDER, 'X')
! 		ENDIF
! 		IF (BS_ORDER(1,PATCH_COLUMN,2).GT.BASIS_KVEC(1,PATCH_COLUMN,2)%POLY_ORDER) THEN
! 			CALL ELEVATE_DEGREE(BASIS_KVEC(1,PATCH_COLUMN,:), BASIS_CTL(1,PATCH_COLUMN), BS_ORDER(1,PATCH_COLUMN,2)-BASIS_KVEC(1,PATCH_COLUMN,2)%POLY_ORDER, 'Y')
! 		ENDIF
! 	ENDDO
	
!-----------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------
! REFINE KNOT VECTOR INSERTING NEW KNOTS WITH MULTIPLICITIES ARE GREATER THAN 1
	
! 	DO PATCH_COLUMN = 1, NUMPATCH(2)
! 		OLD_BASIS_KVEC(:) = BASIS_KVEC(1,PATCH_COLUMN,:)
! 		DO DIR = 1, 2
! 	! 		PRINT*, 'DIR', DIR
! 			IF (EXTRA_KNOTS(1,PATCH_COLUMN,DIR).GT.0) THEN
! 				IF (MESHTYPE=='SH' .AND. PATCH_COLUMN==3 .AND. DIR==2) THEN
! 					!! ----------------------- [[[[[ Shishkin type k-refinement ]]]]] ------------------------------
! 					DO I = 1, EXTRA_KNOTS(1, PATCH_COLUMN, DIR)
! 						K = BASIS_KVEC(1, PATCH_COLUMN, DIR)%LENGTH - (BASIS_KVEC(1, PATCH_COLUMN, DIR)%POLY_ORDER + 1) + I
! 						BASIS_KVEC(1, PATCH_COLUMN, DIR)%KNOTS(K) = &
! 						1.0D0 - (EXTRA_KNOTS(1,PATCH_COLUMN,DIR) + 1  - I)*(KAPPA*BASIS_KVEC(1, PATCH_COLUMN, DIR)%POLY_ORDER/(1.0D0*EXTRA_KNOTS(1,PATCH_COLUMN,DIR)))
! 					ENDDO
! 					DO I = 1, BASIS_KVEC(1, PATCH_COLUMN, DIR)%POLY_ORDER + 1
! 						K = BASIS_KVEC(1, PATCH_COLUMN, DIR)%LENGTH - (BASIS_KVEC(1, PATCH_COLUMN, DIR)%POLY_ORDER + 1) + EXTRA_KNOTS(1, PATCH_COLUMN, DIR) + I
! 						BASIS_KVEC(1, PATCH_COLUMN, DIR)%KNOTS(K) = 1.0D0
! 					ENDDO
! 					BASIS_KVEC(1, PATCH_COLUMN, DIR)%LENGTH = BASIS_KVEC(1, PATCH_COLUMN, DIR)%LENGTH + EXTRA_KNOTS(1,PATCH_COLUMN,DIR)
! 					
! 					DO I = 0, BASIS_KVEC(1, PATCH_COLUMN, DIR)%LENGTH
! 						TMP_KNOTS(BASIS_KVEC(1, PATCH_COLUMN, DIR)%LENGTH - I) = 1.0D0 - BASIS_KVEC(1, PATCH_COLUMN, DIR)%KNOTS(I)
! 					ENDDO
! 					BASIS_KVEC(1, PATCH_COLUMN, DIR)%KNOTS(:) = TMP_KNOTS(:)
! 					!-----------------------------------------------------------------------------------------------------
! 				ELSE
! 					BASIS_KVEC(1,PATCH_COLUMN,DIR) = UNIFORM_KNOT_INSERTION(BASIS_KVEC(1,PATCH_COLUMN,DIR), EXTRA_KNOTS(1,PATCH_COLUMN,DIR),NURBS_REGULARITY(1,PATCH_COLUMN,DIR))
! 				ENDIF
! 				
! 				!-----------------------------------------------------------------------------------------------------
! 				!	GAMMA - RADICAL MESH, (Here, gamma = 5)
! ! 				IF (DIR==1 .AND. (PATCH==2 .OR. PATCH==3)) THEN
! ! 					DO I=1, INT(EXTRA_KNOTS(1))
! ! 						IGA_KVEC(PATCH,DIR)%KNOTS(IGA_KVEC(PATCH,DIR)%POLY_ORDER + I) = 1.D0 - ((EXTRA_KNOTS(1)-I+1)*1.D0/((EXTRA_KNOTS(1)+1)*1.D0))**5
! ! 					ENDDO
! ! 				ELSEIF (DIR==1 .AND. (PATCH==1)) THEN
! ! 					DO I=1, INT(EXTRA_KNOTS(1))
! ! 						IGA_KVEC(PATCH,DIR)%KNOTS(IGA_KVEC(PATCH,DIR)%POLY_ORDER + I) = (I*1.D0/((EXTRA_KNOTS(1)+1)*1.D0))**5
! ! 					ENDDO
! ! 				ENDIF
! ! 				
! ! 				IF (DIR==2) THEN
! ! 					DO I=1, INT(EXTRA_KNOTS(2))
! ! 						BASIS_KVEC(DIR)%KNOTS(BASIS_KVEC(DIR)%POLY_ORDER + I) = (I*1.D0/((EXTRA_KNOTS(2)+1)*1.D0))**5
! ! 					ENDDO
! ! 				ELSEIF (DIR==2 .AND. (PATCH==3)) THEN
! ! 					DO I=1, INT(EXTRA_KNOTS(2))
! ! 						IGA_KVEC(PATCH,DIR)%KNOTS(IGA_KVEC(PATCH,DIR)%POLY_ORDER + I) = 1.D0 - ((EXTRA_KNOTS(2)-I+1)*1.D0/((EXTRA_KNOTS(2)+1)*1.D0))**5
! ! 					ENDDO
! ! 				ENDIF
! 				!-----------------------------------------------------------------------------------------------------
! 
! 				CALL KNOT_TO_ARRAY(OLD_BASIS_KVEC(DIR), OLD_BASIS_POLY_ORDER, OLD_BASIS_NUM_KNOTS, OLD_BASIS_KNOTS, OLD_BASIS_MULTIPLICITIES)
! 				CALL KNOT_TO_ARRAY(BASIS_KVEC(1,PATCH_COLUMN,DIR), BASIS_POLY_ORDER, BASIS_NUM_KNOTS, BASIS_KNOTS, BASIS_MULTIPLICITIES)
! 				
! 				JJ = 1
! 				K = -1
! 				DO II = 1, BASIS_NUM_KNOTS
! 					IF (DABS(BASIS_KNOTS(II) - OLD_BASIS_KNOTS(JJ)).LE.EPS .AND. BASIS_MULTIPLICITIES(II).EQ.OLD_BASIS_MULTIPLICITIES(JJ)) THEN
! 						JJ = JJ + 1
! 					ELSEIF (BASIS_KNOTS(II).LT.OLD_BASIS_KNOTS(JJ)) THEN
! 						DO I = 1, BASIS_MULTIPLICITIES(II)
! 							K = K + 1
! 							NEW_KNOTS(K) = BASIS_KNOTS(II)
! 						ENDDO
! 					ELSEIF (DABS(BASIS_KNOTS(II) - OLD_BASIS_KNOTS(JJ)).LE.EPS .AND. BASIS_MULTIPLICITIES(II).NE.OLD_BASIS_MULTIPLICITIES(JJ)) THEN
! 						DO I = 1, BASIS_MULTIPLICITIES(II) - OLD_BASIS_MULTIPLICITIES(JJ)
! 							K = K + 1
! 							NEW_KNOTS(K) = BASIS_KNOTS(II)
! 						ENDDO
! 					ENDIF
! 				ENDDO
! 				IF (DIR==1) THEN
! 					AXIS = 'X'
! 				ELSEIF (DIR==2) THEN
! 					AXIS = 'Y'
! 				ENDIF
! 				CALL REFINE_KNOT(OLD_BASIS_KVEC(:), BASIS_CTL(1,PATCH_COLUMN), NEW_KNOTS(0:K), K, AXIS)
! 			ENDIF	
! 		ENDDO
! 	ENDDO
	
! 	PRINT*, 'AFTER REFINEMENT', BS_ORDER(1,2,2), BASIS_KVEC(1,2,2)%POLY_ORDER
!-----------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------
! REFINE KNOT VECTOR INSERTING NEW KNOTS WITH MULTIPLICITIES 1
! 	IF (EXTRA_KNOTS(1).GT.0) THEN
! 		DO I = 1, EXTRA_KNOTS(1)
! 			NEW_KNOTS(1,I-1) = I*1.D0/(1.D0*(EXTRA_KNOTS(1)+1))
! 		ENDDO
! 		CALL REFINE_KNOT(BASIS_KVEC(:), BASIS_CTL, NEW_KNOTS(1,0:EXTRA_KNOTS(1)-1), EXTRA_KNOTS(1)-1, 'X')
! 	ENDIF
! 	
! 	IF (EXTRA_KNOTS(2).GT.0) THEN
! 		DO J = 1, EXTRA_KNOTS(2)
! 			NEW_KNOTS(2,J-1) = J*1.D0/(1.D0*(EXTRA_KNOTS(2)+1))
! 		ENDDO
! 		CALL REFINE_KNOT(BASIS_KVEC(:), BASIS_CTL, NEW_KNOTS(2,0:EXTRA_KNOTS(2)-1), EXTRA_KNOTS(2)-1, 'Y')
! 	ENDIF
!-----------------------------------------------------------------------------------------------------------

! TEST REFINE KNOT
! 	CALL REFINE_KNOT(BASIS_KVEC(:), BASIS_CTL, NEW_KNOTS(0:K), K, 'X')

!-----------------------------------------------------------------------------------------------------------
! INSERT A NEW KNOT INTO A EXISTING KNOT VECTOR

! 	BASIS_KVEC(1) = UNIFORM_KNOT_INSERTION(BASIS_KVEC(1), EXTRA_KNOTS(1))
! 	CALL H_REFINEMENT(BASIS_CTL, BASIS_KVEC(:), 0.5D0, 'X')
! 
! 	BASIS_KVEC(2) = UNIFORM_KNOT_INSERTION(BASIS_KVEC(2), EXTRA_KNOTS(2))
! 	CALL H_REFINEMENT(BASIS_CTL, BASIS_KVEC(:), 0.5D0, 'Y')
!-----------------------------------------------------------------------------------------------------------
	

	!-----------------------------------------------------------------------------------------------------------
!! BASIS FUNCTION ON THE OMEGA_HAT_F
	BASIS_KVEC(1,1)%POLY_ORDER = 1
	BASIS_KVEC(1,1)%LENGTH = 3
	BASIS_KVEC(1,1)%KNOTS(0:BASIS_KVEC(1,1)%LENGTH) = (/0.D0,0.D0,1.D0,1.D0/)
    
    BASIS_KVEC(1,2)%POLY_ORDER = 1
	BASIS_KVEC(1,2)%LENGTH = 3
	BASIS_KVEC(1,2)%KNOTS(0:BASIS_KVEC(1,2)%LENGTH) = (/0.D0,0.D0,1.D0,1.D0/)
	
!! BASIS FUNCTION ON THE OMEGA_HAT_G
	BASIS_KVEC(2,1)%POLY_ORDER = 1
	BASIS_KVEC(2,1)%LENGTH = 3
	BASIS_KVEC(2,1)%KNOTS(0:BASIS_KVEC(2,1)%LENGTH) = (/0.D0,0.D0,1.D0,1.D0/)
    
    BASIS_KVEC(2,2)%POLY_ORDER = 1
	BASIS_KVEC(2,2)%LENGTH = 3
	BASIS_KVEC(2,2)%KNOTS(0:BASIS_KVEC(2,2)%LENGTH) = (/0.D0,0.D0,1.D0,1.D0/)
!-----------------------------------------------------------------------------------------------------------


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! [[[[[ REFINEMENT ]]]]] !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!------------------------------------- [[[[[ B-SPLINE ]]]]] -----------------------------------------------
!-----------------------------------------------------------------------------------------------------------
! ELEVATE DEGREE OF B-SPLINE BASIS FUNCTION (p-refinement)
	DO K = 1, NUMPATCH
        DO J = 1, 2
            IF (BS_ORDER(K,J)>BASIS_KVEC(K,J)%POLY_ORDER) THEN
                DO I = 1, BS_ORDER(K,J) - BASIS_KVEC(K,J)%POLY_ORDER
                    BASIS_KVEC(K,J) = DEGREE_ELEVATION(BASIS_KVEC(K,J))
                ENDDO
            ENDIF
        ENDDO
	ENDDO
!-----------------------------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------------------------
! INSERT NEW KNOTS IN THE KNOTS CORRESPONDING TO B-SPLINE BASIS FUNCTIONS (h-refinement)
	DO K = 1, NUMPATCH
        DO J = 1, 2
            IF (EXTRA_KNOTS(K,J)>0) THEN
                BASIS_KVEC(K,J) = UNIFORM_KNOT_INSERTION(BASIS_KVEC(K,J), EXTRA_KNOTS(K,J))
            ENDIF
        ENDDO
	ENDDO
!-----------------------------------------------------------------------------------------------------------

	WRITE(*,*)
	WRITE(*,*) '<<< SET BASIS KNOT : DONE >>>'
	WRITE(*,*)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!------------------------------- [[[[[ COMPUTE DOF ]]]]] -----------------------------------	
! 	NUMBS(:) = (/BASIS_KVEC(1)%LENGTH - BS_ORDER(1)-1, BASIS_KVEC(2)%LENGTH - BS_ORDER(2)/) ! IMPOSED PERIODIC CONDTION
	
! 	print*, 'ubound', ubound(EXPO_ARRAY,1)
	
	NUMBS(:) = 0
	DO PATCH = 1, NUMPATCH
        LOC_NUMBS(PATCH, 1) = BASIS_KVEC(PATCH, 1)%LENGTH - BS_ORDER(PATCH, 1)
        IF (PATCH==1) THEN
            LOC_NUMBS(PATCH, 2) = BASIS_KVEC(PATCH, 2)%LENGTH - BS_ORDER(PATCH, 2) + UBOUND(EXPO_ARRAY, 1)
        ELSE   
            LOC_NUMBS(PATCH, 2) = BASIS_KVEC(PATCH, 2)%LENGTH - BS_ORDER(PATCH, 2)
        ENDIF
    ENDDO
	
	DO PATCH = 1, NUMPATCH
        NUMBS(PATCH) = PRODUCT(LOC_NUMBS(PATCH,:))
	ENDDO
	
	DOF = 0
	DO PATCH = 1, NUMPATCH
			DOF = DOF + PRODUCT(LOC_NUMBS(PATCH, :))
		ENDDO
	
!! ------------------------------------ [[[[[ GLOBAL INDEX ]]]]] -------------------------------------------
!-----------------------------------------------------------------------------------------------------------
	! PU-IGA
	ALLOCATE(NDX(3,DOF))
	I = 1
	DO PATCH = 1, NUMPATCH
        DO J = 1, PRODUCT(LOC_NUMBS(PATCH, :))
            NDX(1, I) = PATCH
            NDX(2, I) = MOD(J + (LOC_NUMBS(PATCH, 1)-1), LOC_NUMBS(PATCH, 1))
            NDX(3, I) = MOD(INT((J-1)/LOC_NUMBS(PATCH, 1)) + LOC_NUMBS(PATCH, 2), LOC_NUMBS(PATCH, 2))
            I = I + 1
        ENDDO
	ENDDO
	
! 	! IGA
! 	DO I = 1, DOF
! 		NDX(1,I)=MOD(I+(NUMBS(1)-1), NUMBS(1))
! 		NDX(2,I)=MOD(INT((I-1)/NUMBS(1))+NUMBS(2), NUMBS(2))
! 	ENDDO
! 	
! 	! RPPM
! 	DO I=1, DOF
!         NDX(1,I)=1+INT((I-1)/(IORDER+1)**2)
!         NDX(2,I)=MOD(I+IORDER, IORDER+1)
!         NDX(3,I)=MOD(INT((I-1)/(IORDER+1))+IORDER+1, IORDER+1)
!     ENDDO
    
	WRITE(*,*)
	WRITE(*,*) '<<< SET GLOBAL INDEX : DONE >>>'
	WRITE(*,*)
!-----------------------------------------------------------------------------------------------

!--------------------------------- [[[[[ PATCHES ON THE PARAMETER SPACE ]]]]] -----------------------------------------
	!! FIND END POINTS OF EACH PATCHES ON THE PARAMETER SPACE
	
	
!-----------------------------------------------------------------------------------------------

!--------------------------------- [[[[[ INTEGRAL REGIONS ON THE PARAMETER SPACE ]]]]] --------------------------------
    
    DO PATCH = 1, NUMPATCH
        DO J = 1, 2
            K = 1
            IR_GRID(PATCH, PATCH, J, 1) = BASIS_KVEC(PATCH, J)%KNOTS(0)
            DO I = 1, BASIS_KVEC(PATCH, J)%LENGTH
                IF (DABS(BASIS_KVEC(PATCH, J)%KNOTS(I-1) - BASIS_KVEC(PATCH, J)%KNOTS(I)).LE.EPS) THEN
                ELSE
                    K = K + 1
                    IR_GRID(PATCH, PATCH, J, K) = BASIS_KVEC(PATCH, J)%KNOTS(I)
                ENDIF
            ENDDO
            NUMIR(PATCH, PATCH, J) = K
        ENDDO
        !-----------------------------------------------------------------------------------------------------------------------------------------------!	
        J = 2
		IF (PATCH==1) THEN    ! singular mapping F
			DO I = 1, NUMIR(PATCH, PATCH, J)
				IF (DABS(IR_GRID(PATCH, PATCH, J, I) - INV_DELTA(1))<=EPS) THEN 
					GOTO 115
				ENDIF
			ENDDO
			!--------------------------------------------------------------------------------!
			TMPIR_GRID(:, :, :, :) = IR_GRID(:, :, :, :)
			TMP_NUMIR(:, :, :) = NUMIR(:, :, :)
			!--------------------------------------------------------------------------------!
			DO I = 1, NUMIR(PATCH, PATCH, J) - 1
				IF (INV_DELTA(1) > IR_GRID(PATCH, PATCH, J, I) .AND. INV_DELTA(1) < IR_GRID(PATCH, PATCH, J, I + 1)) THEN 
					IR_GRID(PATCH, PATCH, J, I + 1) = INV_DELTA(1)
					TMPK = I
					GOTO 116
				ENDIF
			ENDDO
			PRINT*, 'ERROR: GEOMETRY - INV_DELTA(1)  IS NOT INCLUDED IN THE IR_GRID'
			STOP
			
			116 CONTINUE
			!--------------------------------------------------------------------------------!
			DO I = TMPK + 2, NUMIR(PATCH, PATCH, J) + 1
				IR_GRID(PATCH, PATCH, J, I) = TMPIR_GRID(PATCH, PATCH, J, I - 1)
			ENDDO
			NUMIR(PATCH, PATCH, J) = NUMIR(PATCH, PATCH, J) + 1
        !--------------------------------------------------------------------------------!
		ELSEIF (PATCH==2) THEN        ! regular mapping G
			DO I = 1, NUMIR(PATCH, PATCH, J)
				IF (DABS(IR_GRID(PATCH, PATCH, J, I) - INV_DELTA(2))<=EPS) THEN 
					GOTO 115
				ENDIF
			ENDDO
			!--------------------------------------------------------------------------------!
			TMPIR_GRID(:, :, :, :) = IR_GRID(:, :, :, :)
			TMP_NUMIR(:, :, :) = NUMIR(:, :, :)
			!--------------------------------------------------------------------------------!
			DO I = 1, NUMIR(PATCH, PATCH, J) - 1
				IF (INV_DELTA(2) > IR_GRID(PATCH, PATCH, J ,I) .AND. INV_DELTA(2) < IR_GRID(PATCH, PATCH, J, I + 1)) THEN 
					IR_GRID(PATCH, PATCH, J, I + 1) = INV_DELTA(2)
					TMPK = I
					GOTO 118
				ENDIF
			ENDDO
			PRINT*, 'ERROR: GEOMETRY - INV_DELTA(2)  IS NOT INCLUDED IN THE IR_GRID'
			STOP
			
			118 CONTINUE
			!--------------------------------------------------------------------------------!
			DO I = TMPK + 2, NUMIR(PATCH, PATCH, J) + 1
				IR_GRID(PATCH, PATCH, J, I) = TMPIR_GRID(PATCH, PATCH, J, I - 1)
			ENDDO
			NUMIR(PATCH, PATCH, J) = NUMIR(PATCH, PATCH, J) + 1
		ENDIF
		115 CONTINUE
	ENDDO

	!--------------------------------------------------------------------------------!
	! intersection of F and G
	! xi direction
	J = 1
    DO I = 1, NUMIR(1, 1, J)
        IR_GRID(1, 2, J, I) = IR_GRID(1, 1, J, I)
    ENDDO   
    NUMIR(1, 2, J) = NUMIR(1, 1, J)
    
    DO I = 1, NUMIR(2, 2, J)
        !--------------------------------------------------------------------------------!
        TMPIR_GRID(1, 2, J, :) = IR_GRID(1, 2, J, :)
        TMP_NUMIR(1, 2, J) = NUMIR(1, 2, J)
        !--------------------------------------------------------------------------------!
        DO II = 1, NUMIR(1, 2, J)-1
            PHYPT = GET_PHY_PT(POINT2D(IR_GRID(2, 2, J, I), 0.0D0), 2)
            PARPT = GET_PAR_PT(PHYPT, 1)
            IF (PARPT%X>IR_GRID(1, 2, J, II) .AND. PARPT%X<IR_GRID(1, 2, J, II+1)) THEN
                IR_GRID(1, 2, J, II + 1) = PARPT%X
                TMPK = II
                GOTO 217
            ENDIF
        ENDDO
        
        GOTO 218
        
        217 CONTINUE
        !--------------------------------------------------------------------------------!
        DO II = TMPK + 2, NUMIR(1, 2, J) + 1
            IR_GRID(1, 2, J, II) = TMPIR_GRID(1, 2, J, II - 1)
        ENDDO
        NUMIR(1, 2, J) = NUMIR(1, 2, J) + 1
        !--------------------------------------------------------------------------------!
        
        218 CONTINUE
        
    ENDDO
    
    ! eta direction
    J = 2
    K = 0
    DO I = 1, NUMIR(1,1, J)
        IF (IR_GRID(1, 1, J, I)>=INV_DELTA(1)) THEN  
            K = K + 1
            IR_GRID(1, 2, J, K) = IR_GRID(1, 1, J, I)
        ENDIF
    ENDDO
    NUMIR(1, 2, J) = K
    !--------------------------------------------------------------------------------!
    K = 0
    ! regular mapping G
    DO II = 2, NUMIR(2, 2, J) - 1
        IF (IR_GRID(2, 2, J, II)<=INV_DELTA(2)) THEN 
!             PRINT*, 'II',II, IR_GRID(2,2,J,II)
            PHYPT = GET_PHY_PT(POINT2D(0.0D0, IR_GRID(2, 2, J, II)), 2)
            PARPT = GET_PAR_PT(PHYPT, 1)
            K = K + 1
            INVIR_GRID(2, 2, J, K) = PARPT%Y
        ENDIF
    ENDDO
    INV_NUMIR(2, 2, J) = K
!     PRINT*, INV_NUMIR(2,2,J), INVIR_GRID(2,2,J,1)
	!--------------------------------------------------------------------------------!
    DO II = 1, INV_NUMIR(2, 2, J)
        DO I = 1, NUMIR(1, 2, J)
            IF (DABS(IR_GRID(1, 2, J, I) - INVIR_GRID(2, 2, J, II))<=EPS) THEN 
                GOTO 135
            ENDIF
        ENDDO
        !--------------------------------------------------------------------------------!
        TMPIR_GRID(1, 2, J, :) = IR_GRID(1, 2, J, :)
        TMP_NUMIR(1, 2, J) = NUMIR(1, 2, J)
        !--------------------------------------------------------------------------------!
        DO I = 1, NUMIR(1, 2, J) - 1
            IF (INVIR_GRID(2, 2, J, II) > IR_GRID(1, 2, J, I) .AND. INVIR_GRID(2, 2, J, II) < IR_GRID(1, 2, J, I + 1)) THEN 
                IR_GRID(1, 2, J, I + 1) = INVIR_GRID(2, 2, J, II)
                TMPK = I
                GOTO 138
            ENDIF
        ENDDO
        PRINT*, 'ERROR: GEOMETRY - INVIR_GRID  IS NOT INCLUDED IN THE IR_GRID'
        STOP
        
        138 CONTINUE
        !--------------------------------------------------------------------------------!
        DO I = TMPK + 2, NUMIR(1, 2, J) + 1
            IR_GRID(1, 2, J, I) = TMPIR_GRID(1, 2, J, I - 1)
        ENDDO
        NUMIR(1, 2, J) = NUMIR(1, 2, J) + 1
        135 CONTINUE
    ENDDO
    !--------------------------------------------------------------------------------!
    DO J = 1, 2
        IR_GRID(2, 1, J, :) = IR_GRID(1, 2, J, :)
        NUMIR(2, 1, J) = NUMIR(1, 2, J)
    ENDDO
    
    !-----------------------------------------------------------------------------------------------------------------------------------------------!	
	
	K = 0
	DO KK = 1, 2
		L2_GRID(KK, :) = IR_GRID(1,1,KK, :)
		K = NUMIR(1,1,KK)
		DO I = 1, NUMIR(1, 2, KK)
			IF (L2_GRID(KK, NUMIR(1, 1, KK))<IR_GRID(1, 2, KK, I)) THEN
				K = K + 1
				L2_GRID(KK, K) = IR_GRID(1, 2, KK, I)
			ENDIF
		ENDDO
! 		DO I = 1, NUMIR(1, 3, KK)
! 			IF (L2_GRID(KK, K)<IR_GRID(1, 3, KK, I)) THEN
! 				K = K + 1
! 				L2_GRID(KK, K) = IR_GRID(1, 3, KK, I)
! 			ENDIF
! 		ENDDO
		L2_NUMIR(KK) = K
	ENDDO
	
	WRITE(*,*)
	WRITE(*,*) '<<< SET INTEGRAL REGION : DONE >>>'
	WRITE(*,*)
!-----------------------------------------------------------------------------------------------

	
!! Zero hohomogeneous boundary conditions
!------------------------------------------------------------------------------
	DO J=1, MAX_LENGTH
		ZERO_BDNDX(J)%LC_NDX(:) = 0
	ENDDO
	ZERO_BDNDX(1:MAX_LENGTH)%GL_NDX = 0
	ZERO_BDNDX(1:MAX_LENGTH)%LC_NUM = 0
    
    BD_DOF = 0
	K = 0
	
! 	IF (ZEROBD=='Y') THEN
		IF (PDE=='FOUR') THEN
			DO I = 1, DOF
                IF (NDX(1, I)==1 .AND. (NDX(2, I)==0 .OR. NDX(2, I)==1 .OR. NDX(2, I)==LOC_NUMBS(1, 1)-2 .OR. NDX(2, I)==LOC_NUMBS(1, 1)-1 .OR. NDX(3, I)==0 .OR. NDX(3, I)==1)) THEN   ! singular mapping F
                    K = K + 1
                    ZERO_BDNDX(K)%LC_NDX(:) = NDX(:,I)
                    ZERO_BDNDX(K)%GL_NDX = I
                    BD_DOF = BD_DOF + 1
                ELSEIF (NDX(1, I)==2 .AND. (NDX(2, I)==0 .OR. NDX(2, I)==1 .OR. NDX(2, I)==LOC_NUMBS(2, 1)-2 .OR. NDX(2, I)==LOC_NUMBS(2, 1)-1 .OR. NDX(3, I)==LOC_NUMBS(2, 2)-2 .OR. NDX(3, I)==LOC_NUMBS(2, 2)-1)) THEN   ! regular mapping G
                    K = K + 1
                    ZERO_BDNDX(K)%LC_NDX(:) = NDX(:,I)
                    ZERO_BDNDX(K)%GL_NDX = I
                    BD_DOF = BD_DOF + 1
                ENDIF
			ENDDO
		ELSEIF (PDE=='ELPT') THEN
			DO I=1, DOF
				IF (NDX(1, I)==1 .AND. (NDX(2, I)==0 .OR. NDX(2, I)==LOC_NUMBS(1, 1)-1 .OR. NDX(3, I)==0)) THEN   ! singular mapping F
					K = K + 1
					ZERO_BDNDX(K)%LC_NDX(:) = NDX(:,I)
					ZERO_BDNDX(K)%GL_NDX = I
					BD_DOF = BD_DOF + 1
                ELSEIF (NDX(1, I)==2 .AND. (NDX(2, I)==0 .OR. NDX(2, I)==LOC_NUMBS(2, 1)-1 .OR. NDX(3, I)==LOC_NUMBS(2, 2)-1)) THEN   ! regular mapping G
                    K = K + 1
					ZERO_BDNDX(K)%LC_NDX(:) = NDX(:,I)
					ZERO_BDNDX(K)%GL_NDX = I
					BD_DOF = BD_DOF + 1
				ENDIF
			ENDDO
		ENDIF
! 	ELSE 
! 	ENDIF
	ZERO_BDNDX(:)%LC_NUM = K
!----------------------------------------------------------------------------------------------

! SET INDEX OF BASIS FUNCTIONS IMPOSED ESSENTIAL BOUNDARY CONDITION
	ALLOCATE(BD_MASK(DOF))
	BD_MASK(:) = .FALSE.
	
	DO J=1, MAX_LENGTH
		BDNDX(J)%LC_NDX(:) = 0
	ENDDO
	BDNDX(1:MAX_LENGTH)%GL_NDX = 0
	BDNDX(1:MAX_LENGTH)%LC_NUM = 0
	
	K = 0
! 	BD_DOF = 0
	
! 	IF (ZEROBD=='N') THEN
! 		DO I=1, DOF
! 			IF (NDX(2,I)==3 .AND. NDX(4,I)==0) THEN
! 				DO II = 1, ZERO_BDNDX(1)%LC_NUM
! 					IF (I==ZERO_BDNDX(II)%GL_NDX) THEN
! 						GOTO 324
! 					ENDIF
! 				ENDDO
! 				K = K + 1
! 				BDNDX(K)%LC_NDX(:) = NDX(:,I)
! 				BDNDX(K)%GL_NDX = I
! 				BD_DOF = BD_DOF + 1
! 			ENDIF
! 			324 CONTINUE
! 		ENDDO
! 	ENDIF
	
	BDNDX(:)%LC_NUM = K
	DEALLOCATE(BD_MASK)
	
	WRITE(*,*)
	WRITE(*,*) '<<< SET LOCAL INDEX CORRESPONDING BASIS FUNCTIONS ON BOUNDARY IMPOSED ESSENTIAL BOUNDARY CONDITION : DONE >>>'
	WRITE(*,*)

!----------------------------------------------------------------------------------------------
!!  GENERATE BINOMIAL COEFFICIENTS

	DO I = 1, MAX_DIFF_ORDER
		DO J = 0, MAX_DIFF_ORDER
			BINOM(I,J) = 1.0D0
		ENDDO
	ENDDO

	BINOM(1,0) = 1.0D0
	BINOM(1,1) = 1.0D0
	DO I = 2, MAX_DIFF_ORDER
		BINOM(I,0) = BINOM(I-1,0)
		DO J = 1, I-1
			BINOM(I,J) = BINOM(I-1,J-1) + BINOM(I-1,J)
		ENDDO
		BINOM(I,I) = BINOM(I-1,I-1)
	ENDDO

	WRITE(*,*)
	WRITE(*,*) '<<< GENERATE BINOMIAL COEFFICIENTS : DONE >>>'
	WRITE(*,*)

END SUBROUTINE GET_GEO

! CHARACTER(LEN=1) FUNCTION GET_GAMMA_HAT(PT)
! 	TYPE(POINT2D), INTENT(IN) :: PT
! 	
! 	IF (DABS(PT%X).LE.EPS .AND. DABS(PT%Y).GT.EPS) THEN
! 		GET_GAMMA_HAT = 'L'
! 	ELSEIF (DABS(PT%X - 1.D0).LE.EPS .AND. DABS(PT%Y).GT.EPS) THEN
! 		GET_GAMMA_HAT = 'R'
! 	ELSEIF (DABS(PT%X).GT.EPS .AND. DABS(PT%Y).LE.EPS) THEN
! 		GET_GAMMA_HAT = 'B'
! 	ELSEIF (DABS(PT%X).GT.EPS .AND. DABS(PT%Y - 1.D0).LE.EPS) THEN
! 		GET_GAMMA_HAT = 'T'
! 	ENDIF
! 	
! END FUNCTION GET_GAMMA_HAT

END MODULE GEOMETRY

