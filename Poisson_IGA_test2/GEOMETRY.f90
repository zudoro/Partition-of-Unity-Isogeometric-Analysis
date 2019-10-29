MODULE GEOMETRY

	USE KNOT_HANDLING
	USE PATCH_MAPPING

! 	IMPLICIT INTEGER (I-N)
! 	IMPLICIT REAL(8) (A-H,O-Z)

CONTAINS

SUBROUTINE GET_GEO()
	
	INTEGER :: I, J, II, JJ, K, KK, DIR
	INTEGER :: BASIS_POLY_ORDER, OLD_BASIS_POLY_ORDER, BASIS_NUM_KNOTS, OLD_BASIS_NUM_KNOTS
	INTEGER :: BASIS_MULTIPLICITIES(MAX_LENGTH), OLD_BASIS_MULTIPLICITIES(MAX_LENGTH)
	REAL(8) :: BASIS_KNOTS(MAX_LENGTH), OLD_BASIS_KNOTS(MAX_LENGTH), NEW_KNOTS(0:MAX_LENGTH)
	CHARACTER(LEN=1) :: AXIS
	
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

! UNIT CIRCLE (EXOTIC GEOMETRIC MAPPING)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	GEO_KVEC(1)%POLY_ORDER = 2
! 	GEO_KVEC(2)%POLY_ORDER = 2
! 
! 	GEO_KVEC(1)%LENGTH = 11
! 	GEO_KVEC(2)%LENGTH = 7
! 
! 	GEO_KVEC(1)%KNOTS(0:GEO_KVEC(1)%LENGTH) = (/0.D0, 0.D0, 0.D0, 0.25D0, 0.25D0, 0.5D0, 0.5D0, 0.75D0, 0.75D0, 1.D0, 1.D0, 1.D0/)
! 	GEO_KVEC(2)%KNOTS(0:GEO_KVEC(2)%LENGTH) = (/0.D0, 0.D0, 0.D0, 0.9D0, 0.9D0, 1.D0, 1.D0, 1.D0/)
! 	
! ! SET CONTROL POINTS AND WEIGTH
! 	GEO_CTL%D=2
! 	
! 	GEO_CTL%PTS(0:8,0:1,0) = POINT2D(0.D0,0.D0)
! 
! 	GEO_CTL%PTS(0,4,0) = POINT2D(1.D0,0.D0)
! 	GEO_CTL%PTS(1,4,0) = POINT2D(1.D0,-1.D0)
! 	GEO_CTL%PTS(2,4,0) = POINT2D(0.D0,-1.D0)
! 	GEO_CTL%PTS(3,4,0) = POINT2D(-1.D0,-1.D0)
! 	GEO_CTL%PTS(4,4,0) = POINT2D(-1.D0,0.D0)
! 	GEO_CTL%PTS(5,4,0) = POINT2D(-1.D0,1.D0)
! 	GEO_CTL%PTS(6,4,0) = POINT2D(0.D0,1.D0)
! 	GEO_CTL%PTS(7,4,0) = POINT2D(1.D0,1.D0)
! 	GEO_CTL%PTS(8,4,0) = POINT2D(1.D0,0.D0)
! 	
! 	DO I = 0, 8
! 		GEO_CTL%PTS(I,4,0) = 0.550D0*GEO_CTL%PTS(I,4,0)
! 		GEO_CTL%PTS(I,2,0) = (0.450D0/0.550D0)*GEO_CTL%PTS(I,4,0)
! 		GEO_CTL%PTS(I,3,0) = 0.5D0*(GEO_CTL%PTS(I,2,0) + GEO_CTL%PTS(I,4,0))
! 	ENDDO
! 	
! 	DO I=0, 4
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
	GEO_KVEC(1)%POLY_ORDER = 2
	GEO_KVEC(2)%POLY_ORDER = 2

	GEO_KVEC(1)%LENGTH = 6
	GEO_KVEC(2)%LENGTH = 5

	GEO_KVEC(1)%KNOTS(0:GEO_KVEC(1)%LENGTH) = (/0.D0,0.D0,0.D0,0.5D0,1.D0,1.D0,1.D0/)
	GEO_KVEC(2)%KNOTS(0:GEO_KVEC(2)%LENGTH) = (/0.D0,0.D0,0.D0,1.D0,1.D0,1.D0/)
	
! SET CONTROL POINTS AND WEIGTH
	GEO_CTL%D=2
	GEO_CTL%PTS(0,0,0) = POINT2D(0.D0,-1.D0)
	GEO_CTL%PTS(1,0,0) = POINT2D(0.D0,0.D0)
	GEO_CTL%PTS(2,0,0) = POINT2D(0.D0,0.D0)
	GEO_CTL%PTS(3,0,0) = POINT2D(1.D0,0.D0)
	GEO_CTL%PTS(0,1,0) = POINT2D(-0.5D0,-1.D0)
	GEO_CTL%PTS(1,1,0) = POINT2D(-2.D0/3.D0,1.D0/3.D0)
	GEO_CTL%PTS(2,1,0) = POINT2D(-1.D0/3.D0,2.D0/3.D0)
	GEO_CTL%PTS(3,1,0) = POINT2D(1.D0,0.5D0)
	GEO_CTL%PTS(0,2,0) = POINT2D(-1.D0,-1.D0)
	GEO_CTL%PTS(1,2,0) = POINT2D(-1.D0,1.D0)
	GEO_CTL%PTS(2,2,0) = POINT2D(-1.D0,1.D0)
	GEO_CTL%PTS(3,2,0) = POINT2D(1.D0,1.D0)
	
	GEO_CTL%WGTS(:,:,0) = 1.D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	WRITE(*,*)
	WRITE(*,*) '<<< SET GEOMETRIC KNOT, CONTROL POINTS, AND WEIGHT : DONE >>>'
	WRITE(*,*)
	
! SET KNOT INFO FOR B-SPLINE BASIS FUNCTIONS
	BASIS_CTL = GEO_CTL
	BASIS_KVEC(:) = GEO_KVEC(:)

!-----------------------------------------------------------------------------------------------------------
! ELEVATE DEGREE OF BERNSTEIN POLY.

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
! ELEVATE DEGREE OF B-SPLINE FUNCTION
	IF (BS_ORDER(1).GT.GEO_KVEC(1)%POLY_ORDER) THEN
		CALL ELEVATE_DEGREE(BASIS_KVEC, BASIS_CTL, BS_ORDER(1)-GEO_KVEC(1)%POLY_ORDER, 'X')
	ENDIF
	IF (BS_ORDER(2).GT.GEO_KVEC(2)%POLY_ORDER) THEN
		CALL ELEVATE_DEGREE(BASIS_KVEC, BASIS_CTL, BS_ORDER(2)-GEO_KVEC(2)%POLY_ORDER, 'Y')
	ENDIF
!-----------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------
! REFINE KNOT VECTOR INSERTING NEW KNOTS WITH MULTIPLICITIES ARE GREATER THAN 1
	
	OLD_BASIS_KVEC(:) = BASIS_KVEC(:)
	DO DIR = 1, 2
! 		PRINT*, 'DIR', DIR
		IF (EXTRA_KNOTS(DIR).GT.0) THEN
			BASIS_KVEC(DIR) = UNIFORM_KNOT_INSERTION(BASIS_KVEC(DIR), EXTRA_KNOTS(DIR),NURBS_REGULARITY(DIR))
			CALL KNOT_TO_ARRAY(OLD_BASIS_KVEC(DIR), OLD_BASIS_POLY_ORDER, OLD_BASIS_NUM_KNOTS, OLD_BASIS_KNOTS, OLD_BASIS_MULTIPLICITIES)
			CALL KNOT_TO_ARRAY(BASIS_KVEC(DIR), BASIS_POLY_ORDER, BASIS_NUM_KNOTS, BASIS_KNOTS, BASIS_MULTIPLICITIES)
			
			JJ = 1
			K = -1
			DO II = 1, BASIS_NUM_KNOTS
				IF (DABS(BASIS_KNOTS(II) - OLD_BASIS_KNOTS(JJ)).LE.EPS .AND. BASIS_MULTIPLICITIES(II).EQ.OLD_BASIS_MULTIPLICITIES(JJ)) THEN
					JJ = JJ + 1
				ELSEIF (BASIS_KNOTS(II).LT.OLD_BASIS_KNOTS(JJ)) THEN
! 					PRINT*, II, BASIS_MULTIPLICITIES(II)
					DO I = 1, BASIS_MULTIPLICITIES(II)
						K = K + 1
						NEW_KNOTS(K) = BASIS_KNOTS(II)
					ENDDO
				ELSEIF (DABS(BASIS_KNOTS(II) - OLD_BASIS_KNOTS(JJ)).LE.EPS .AND. BASIS_MULTIPLICITIES(II).NE.OLD_BASIS_MULTIPLICITIES(JJ)) THEN
					DO I = 1, BASIS_MULTIPLICITIES(II) - OLD_BASIS_MULTIPLICITIES(JJ)
						K = K + 1
						NEW_KNOTS(K) = BASIS_KNOTS(II)
					ENDDO
				ENDIF
			ENDDO
			IF (DIR==1) THEN
				AXIS = 'X'
			ELSEIF (DIR==2) THEN
				AXIS = 'Y'
			ENDIF
			CALL REFINE_KNOT(OLD_BASIS_KVEC(:), BASIS_CTL, NEW_KNOTS(0:K), K, AXIS)
		ENDIF	
	ENDDO
	
! 	PRINT*, 'K', K
! 	PRINT*, 'AXIS', AXIS
! 	
! 	DO I = 0, K
! 		PRINT*, NEW_KNOTS(I)
! 	ENDDO
	
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
	

	WRITE(*,*)
	WRITE(*,*) '<<< SET BASIS KNOT : DONE >>>'
	WRITE(*,*)

	NUMBS(:) = (/BASIS_KVEC(1)%LENGTH - BS_ORDER(1), BASIS_KVEC(2)%LENGTH - BS_ORDER(2)/)
	DOF = NUMBS(1)*NUMBS(2)
	
! GLOBAL INDEX
	ALLOCATE(NDX(2,DOF))
	DO I=1, DOF
		NDX(1,I)=MOD(I+(NUMBS(1)-1), NUMBS(1))
		NDX(2,I)=MOD(INT((I-1)/NUMBS(1))+NUMBS(2), NUMBS(2))
	ENDDO

	WRITE(*,*)
	WRITE(*,*) '<<< SET GLOBAL INDEX : DONE >>>'
	WRITE(*,*)

! CONSTRUCT THE FRAME OF LOCAL INTEGRAL REGIONS
	IR_GRID(1,1) = 0.D0
	IR_GRID(2,1) = 0.D0
	
! 	PRINT*, BASIS_KVEC(1)%LENGTH, BASIS_KVEC(2)%LENGTH
	K = 1
	DO I=1, BASIS_KVEC(1)%LENGTH
		K = K + 1
		IF (DABS(BASIS_KVEC(1)%KNOTS(I-1) - BASIS_KVEC(1)%KNOTS(I)).LE.EPS) THEN
			K = K - 1
		ELSE
			IR_GRID(1,K) = BASIS_KVEC(1)%KNOTS(I)
		ENDIF
! 		IF (K.GE.2) THEN
! 			DO J=1, 
! 		ENDIF
	ENDDO
	NUMIR(1) = K
	
	K = 1
	DO I=1, BASIS_KVEC(2)%LENGTH
		K = K + 1
		IF (DABS(BASIS_KVEC(2)%KNOTS(I-1) - BASIS_KVEC(2)%KNOTS(I)).LE.EPS) THEN
			K = K - 1
		ELSE
			IR_GRID(2,K) = BASIS_KVEC(2)%KNOTS(I)
		ENDIF
	ENDDO
	NUMIR(2) = K

	WRITE(*,*)
	WRITE(*,*) '<<< SET INTEGRAL REGION : DONE >>>'
	WRITE(*,*)

! SET INDEX OF BASIS FUNCTIONS IMPOSED ESSENTIAL BOUNDARY CONDITION
	ALLOCATE(BD_MASK(DOF))
	BD_MASK(:) = .FALSE.
	BD_DOF = 0
	
	K = 0
	DO I=1, DOF
		IF (NDX(1,I)==0 .OR. NDX(1,I)==NUMBS(1)-1 .OR. NDX(2,I)==NUMBS(2)-1 .OR. NDX(2,I)==0) THEN
			K = K + 1
			BDNDX(K)%LC_NDX(:) = NDX(:,I)
			BDNDX(K)%GL_NDX = I
			BD_MASK(I)=.TRUE.
			BD_DOF = BD_DOF + 1
		ENDIF
	ENDDO
	
	BD_DOF = BD_DOF
	BDNDX(:)%LC_NUM = K
		
	WRITE(*,*)
	WRITE(*,*) '<<< SET LOCAL INDEX CORRESPONDING BASIS FUNCTIONS ON BOUNDARY IMPOSED ESSENTIAL BOUNDARY CONDITION : DONE >>>'
	WRITE(*,*)

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

CHARACTER(LEN=1) FUNCTION GET_GAMMA_HAT(PT)
	TYPE(POINT2D), INTENT(IN) :: PT
	
	IF (DABS(PT%X).LE.EPS .AND. DABS(PT%Y).GT.EPS) THEN
		GET_GAMMA_HAT = 'L'
	ELSEIF (DABS(PT%X - 1.D0).LE.EPS .AND. DABS(PT%Y).GT.EPS) THEN
		GET_GAMMA_HAT = 'R'
	ELSEIF (DABS(PT%X).GT.EPS .AND. DABS(PT%Y).LE.EPS) THEN
		GET_GAMMA_HAT = 'B'
	ELSEIF (DABS(PT%X).GT.EPS .AND. DABS(PT%Y - 1.D0).LE.EPS) THEN
		GET_GAMMA_HAT = 'T'
	ENDIF
	
END FUNCTION GET_GAMMA_HAT

END MODULE GEOMETRY

