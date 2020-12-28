MODULE PATCH_MAPPING

	USE NURBS

	implicit integer (i-n)
	implicit real(8) (a-h,o-z)

CONTAINS

	!!  GET THE CORRESPONDING POINT ON THE PHYSICAL SPACE OF GIVE POINT ON THE PARAMETRIC SPACE
	TYPE(POINT2D) FUNCTION GET_PHY_PT(PAR_PT)

		TYPE(POINT2D), INTENT(IN) :: PAR_PT
		REAL*8 :: FACTOR
		
		FACTOR = 1.0D0 - PAR_PT%Y
		GET_PHY_PT%X = FACTOR*DCOS(2.0D0*PI*PAR_PT%X)
		GET_PHY_PT%Y = FACTOR*DSIN(2.0D0*PI*PAR_PT%X)

	END FUNCTION GET_PHY_PT


	!!  GET JACOBIAN OF THE GEOMETRIC MAPPING
	TYPE(MATRIX_22) FUNCTION GET_JACOBIAN_MATRIX(PAR_PT)

		TYPE(POINT2D), INTENT(IN) :: PAR_PT
		TYPE(FUNCTION_2D) :: DIFF_NURBS(2)
		
		GET_JACOBIAN_MATRIX%ENT(1,1) = -2.0D0*PI*(1.0D0 - PAR_PT%Y)*DSIN(2.0D0*PI*PAR_PT%X)
		GET_JACOBIAN_MATRIX%ENT(1,2) = -DCOS(2.0D0*PI*PAR_PT%X)
		GET_JACOBIAN_MATRIX%ENT(2,1) = 2.0D0*PI*(1.0D0 - PAR_PT%Y)*DCOS(2.0D0*PI*PAR_PT%X)
		GET_JACOBIAN_MATRIX%ENT(2,2) = -DSIN(2.0D0*PI*PAR_PT%X)

	END FUNCTION GET_JACOBIAN_MATRIX
	
	!!  GET SECOND-ORDER PARTIAL DERIVATIVES
! 	TYPE(SECOND_PARTIAL_DERIVATIVES) FUNCTION GET_PARTIAL(REF_PT)
! 
! 		TYPE(POINT2D), INTENT(IN) :: REF_PT
! 
! 		TYPE(FUNCTION_2D) :: DIFF_NURBS(2)
! 
! 		CALL GET_DIFF_NURVE_SURFACE_2D(DIFF_NURBS,REF_PT,GEO_KVEC,GEO_CTL,2)
! 
! 		GET_PARTIAL%X(1) = DIFF_NURBS(1)%VAL(2,0)
! 		GET_PARTIAL%X(2) = DIFF_NURBS(1)%VAL(1,1)
! 		GET_PARTIAL%X(3) = DIFF_NURBS(1)%VAL(0,2)
! 
! 		GET_PARTIAL%Y(1) = DIFF_NURBS(2)%VAL(2,0)
! 		GET_PARTIAL%Y(2) = DIFF_NURBS(2)%VAL(1,1)
! 		GET_PARTIAL%Y(3) = DIFF_NURBS(2)%VAL(0,2)
! 
! 	END FUNCTION GET_PARTIAL

! 	REAL*8 FUNCTION GET_DET_LN(REF_PT, GAMMA)
! 	
! 		TYPE(POINT2D), INTENT(IN) :: REF_PT
! 		INTEGER, INTENT(IN) :: GAMMA
! 	
! 		TYPE(FUNCTION_2D) :: DIFF_NURBS(2)
! 	
! 		CALL GET_DIFF_NURVE_SURFACE_2D(DIFF_NURBS,REF_PT,GEO_KVEC(:),GEO_CTL,1)
! 	
! 
! 		IF (GAMMA.EQ.2 .OR. GAMMA.EQ.4 .OR. GAMMA.EQ.6) THEN
! ! 			PRINT*, DIFF_NURBS(1:2)%VAL(0,1)
! 			GET_DET_LN = DSQRT(DIFF_NURBS(1)%VAL(0,1)**2 + DIFF_NURBS(2)%VAL(0,1)**2)
! 		ELSEIF (GAMMA==1 .OR. GAMMA==3 .OR. GAMMA==5) THEN
! 			GET_DET_LN = DSQRT(DIFF_NURBS(1)%VAL(1,0)**2 + DIFF_NURBS(2)%VAL(1,0)**2)
! 		ELSE
! 			PRINT*, 'ERROR - GET_DET_LN : 1'
! 			STOP
! 		ENDIF
! 
! 	END FUNCTION GET_DET_LN
! 	
! 	REAL*8 FUNCTION GET_DET_DR(REF_PT, GAMMA_HAT)
! 	
! 		TYPE(POINT2D), INTENT(IN) :: REF_PT
! 		CHARACTER(LEN=1), INTENT(IN) :: GAMMA_HAT
! 	
! 		TYPE(FUNCTION_2D) :: DIFF_NURBS(2)
! 	
! 		CALL GET_DIFF_NURVE_SURFACE_2D(DIFF_NURBS,REF_PT,GEO_KVEC(:),GEO_CTL,1)
! 	
! 		IF (GAMMA_HAT=='B' .OR. GAMMA_HAT=='T') THEN
! 			GET_DET_DR = DSQRT(DIFF_NURBS(1)%VAL(1,0)**2 + DIFF_NURBS(2)%VAL(1,0)**2)
! 		ELSEIF (GAMMA_HAT=='L' .OR. GAMMA_HAT=='R') THEN
! 			GET_DET_DR = DSQRT(DIFF_NURBS(1)%VAL(0,1)**2 + DIFF_NURBS(2)%VAL(0,1)**2)
! 		ELSE
! 			PRINT*, 'ERROR - GET_DET_DR : 1'
! 				STOP
! 		ENDIF
! 
! 	END FUNCTION GET_DET_DR


!!-------------------[[[[[ PATCH + 2*DELTA TO [-1, 1] OF PU DOMAIN ]]]]] ------------------------
TYPE(FUNCTION_1D) FUNCTION PATCH_TO_REFPU(PHYPT, PATCH, AXIS)

    REAL*8, INTENT(IN) :: PHYPT
		INTEGER, INTENT(IN) :: PATCH
		CHARACTER(LEN=1), INTENT(IN) :: AXIS
		
    REAL*8 :: XPTS(4), YPTS(4)
		
		IF (AXIS=='X') THEN
			IF (PATCH==1) THEN
				XPTS(1) = PATCHBDX(PATCH,1) - 2.0D0*DELTA
				XPTS(2) = PATCHBDX(PATCH,1)
				XPTS(3) = PATCHBDX(PATCH,2) - DELTA
				XPTS(4) = PATCHBDX(PATCH,2) + DELTA
				
				IF (PHYPT>=XPTS(2) .AND. PHYPT<=XPTS(3)) THEN
					PATCH_TO_REFPU%VAL(0) = 0.0D0
					PATCH_TO_REFPU%VAL(1) = 0.0D0
				ELSEIF (PHYPT>=XPTS(3) .AND. PHYPT<=XPTS(4)) THEN
					PATCH_TO_REFPU%VAL(0) = (PHYPT-XPTS(3))/(XPTS(4)-XPTS(3))
					PATCH_TO_REFPU%VAL(1) = 1.0D0/(XPTS(4)-XPTS(3))
				ELSE 
					PATCH_TO_REFPU%VAL(0) = 5.0D0
					PATCH_TO_REFPU%VAL(1) = 0.0D0
				ENDIF
				GOTO 11
			ELSEIF (PATCH==NUMPATCH(1)) THEN
				XPTS(1) = PATCHBDX(PATCH,1) - DELTA
				XPTS(2) = PATCHBDX(PATCH,1) + DELTA
				XPTS(3) = PATCHBDX(PATCH,2)
				XPTS(4) = PATCHBDX(PATCH,2) + 2.0D0*DELTA
			
				IF (PHYPT>=XPTS(1) .AND. PHYPT<=XPTS(2)) THEN
					PATCH_TO_REFPU%VAL(0) = (PHYPT-XPTS(2))/(XPTS(2)-XPTS(1))
					PATCH_TO_REFPU%VAL(1) = 1.0D0/(XPTS(2)-XPTS(1))
				ELSEIF (PHYPT>=XPTS(2) .AND. PHYPT<=XPTS(3)) THEN
					PATCH_TO_REFPU%VAL(0) = 0.0D0
					PATCH_TO_REFPU%VAL(1) = 0.0D0
				ELSE 
					PATCH_TO_REFPU%VAL(0) = 5.0D0
					PATCH_TO_REFPU%VAL(1) = 0.0D0
				ENDIF
				GOTO 11
			ELSE
				XPTS(1) = PATCHBDX(PATCH,1) - DELTA
				XPTS(2) = PATCHBDX(PATCH,1) + DELTA
				XPTS(3) = PATCHBDX(PATCH,2) - DELTA
				XPTS(4) = PATCHBDX(PATCH,2) + DELTA
				
				IF (PHYPT>=XPTS(1) .AND. PHYPT<=XPTS(2)) THEN
					PATCH_TO_REFPU%VAL(0) = (PHYPT-XPTS(2))/(XPTS(2)-XPTS(1))
					PATCH_TO_REFPU%VAL(1) = 1.0D0/(XPTS(2)-XPTS(1))
				ELSEIF (PHYPT>=XPTS(2) .AND. PHYPT<=XPTS(3)) THEN
					PATCH_TO_REFPU%VAL(0) = 0.0D0
					PATCH_TO_REFPU%VAL(1) = 0.0D0
				ELSEIF (PHYPT>=XPTS(3) .AND. PHYPT<=XPTS(4)) THEN
					PATCH_TO_REFPU%VAL(0) = (PHYPT-XPTS(3))/(XPTS(4)-XPTS(3))
					PATCH_TO_REFPU%VAL(1) = 1.0D0/(XPTS(4)-XPTS(3))
				ELSE 
					PATCH_TO_REFPU%VAL(0) = 5.0D0
					PATCH_TO_REFPU%VAL(1) = 0.0D0
				ENDIF
				GOTO 11
			ENDIF
			IF ((DABS(PHYPT-XPTS(2)).LE.EPS) .OR. (DABS(PHYPT-XPTS(3)).LE.EPS)) THEN
				PATCH_TO_REFPU%VAL(0) = 0.0D0
				PATCH_TO_REFPU%VAL(1) = 0.0D0
			ENDIF
			GOTO 11
			IF ((DABS(PHYPT-XPTS(1)).LE.EPS) .OR. (DABS(PHYPT-XPTS(4)).LE.EPS)) THEN
				PATCH_TO_REFPU%VAL(0) = 0.0D0
				PATCH_TO_REFPU%VAL(1) = 0.0D0
			ENDIF
			GOTO 11
		ELSEIF (AXIS=='Y') THEN
			IF (PATCH==NUMPATCH(2)) THEN
				YPTS(1) = PATCHBDY(PATCH,1) - 2.0D0*DELTA
				YPTS(2) = PATCHBDY(PATCH,1)
				YPTS(3) = PATCHBDY(PATCH,2) - DELTA
				YPTS(4) = PATCHBDY(PATCH,2) + DELTA
				
				IF (PHYPT>=YPTS(2) .AND. PHYPT<=YPTS(3)) THEN
					PATCH_TO_REFPU%VAL(0) = 0.0D0
					PATCH_TO_REFPU%VAL(1) = 0.0D0
				ELSEIF (PHYPT>=YPTS(3) .AND. PHYPT<=YPTS(4)) THEN
					PATCH_TO_REFPU%VAL(0) = (PHYPT-YPTS(3))/(YPTS(4)-YPTS(3))
					PATCH_TO_REFPU%VAL(1) = 1.0D0/(YPTS(4)-YPTS(3))
				ELSE 
					PATCH_TO_REFPU%VAL(0) = 5.0D0
					PATCH_TO_REFPU%VAL(1) = 0.0D0
				ENDIF
				GOTO 12
			ELSEIF (PATCH==1) THEN
				YPTS(1) = PATCHBDY(PATCH,1) - DELTA
				YPTS(2) = PATCHBDY(PATCH,1) + DELTA
				YPTS(3) = PATCHBDY(PATCH,2)
				YPTS(4) = PATCHBDY(PATCH,2) + 2.0D0*DELTA
			
				IF (PHYPT>=YPTS(1) .AND. PHYPT<=YPTS(2)) THEN
					PATCH_TO_REFPU%VAL(0) = (PHYPT-YPTS(2))/(YPTS(2)-YPTS(1))
					PATCH_TO_REFPU%VAL(1) = 1.0D0/(YPTS(2)-YPTS(1))
				ELSEIF (PHYPT>=YPTS(2) .AND. PHYPT<=YPTS(3)) THEN
					PATCH_TO_REFPU%VAL(0) = 0.0D0
					PATCH_TO_REFPU%VAL(1) = 0.0D0
				ELSE 
					PATCH_TO_REFPU%VAL(0) = 5.0D0
					PATCH_TO_REFPU%VAL(1) = 0.0D0
				ENDIF
				GOTO 12
			ELSE
				YPTS(1) = PATCHBDY(PATCH,1) - DELTA
				YPTS(2) = PATCHBDY(PATCH,1) + DELTA
				YPTS(3) = PATCHBDY(PATCH,2) - DELTA
				YPTS(4) = PATCHBDY(PATCH,2) + DELTA
				
				IF (PHYPT>=YPTS(1) .AND. PHYPT<=YPTS(2)) THEN
					PATCH_TO_REFPU%VAL(0) = (PHYPT-YPTS(2))/(YPTS(2)-YPTS(1))
					PATCH_TO_REFPU%VAL(1) = 1.0D0/(YPTS(2)-YPTS(1))
				ELSEIF (PHYPT>=YPTS(2) .AND. PHYPT<=YPTS(3)) THEN
					PATCH_TO_REFPU%VAL(0) = 0.0D0
					PATCH_TO_REFPU%VAL(1) = 0.0D0
				ELSEIF (PHYPT>=YPTS(3) .AND. PHYPT<=YPTS(4)) THEN
					PATCH_TO_REFPU%VAL(0) = (PHYPT-YPTS(3))/(YPTS(4)-YPTS(3))
					PATCH_TO_REFPU%VAL(1) = 1.0D0/(YPTS(4)-YPTS(3))
				ELSE 
					PATCH_TO_REFPU%VAL(0) = 5.0D0
					PATCH_TO_REFPU%VAL(1) = 0.0D0
				ENDIF
				GOTO 12
			ENDIF
			IF ((DABS(PHYPT-YPTS(2)).LE.EPS) .OR. (DABS(PHYPT-YPTS(3)).LE.EPS)) THEN
				PATCH_TO_REFPU%VAL(0) = 0.0D0
				PATCH_TO_REFPU%VAL(1) = 0.0D0
			ENDIF
			GOTO 12
			IF ((DABS(PHYPT-YPTS(1)).LE.EPS) .OR. (DABS(PHYPT-YPTS(4)).LE.EPS)) THEN
				PATCH_TO_REFPU%VAL(0) = 0.0D0
				PATCH_TO_REFPU%VAL(1) = 0.0D0
			ENDIF
			GOTO 12
		ENDIF
		11 CONTINUE
		12 CONTINUE
END FUNCTION PATCH_TO_REFPU


!----------------------------------------------------------------------------------------
TYPE(FUNCTION_1D) FUNCTION UNITSQ_TO_PARSP(REFPT, PATCHES, AXIS)
	
	REAL*8, INTENT(IN) :: REFPT
	INTEGER, INTENT(IN) :: PATCHES(2)
	CHARACTER(LEN=1), INTENT(IN) :: AXIS
	
	REAL*8 :: LEFT_PT, RIGHT_PT
	IF (AXIS=='X') THEN
		
		UNITSQ_TO_PARSP%VAL(0) = REFPT
		UNITSQ_TO_PARSP%VAL(1) = 1.0D0
	
	ELSEIF (AXIS=='Y') THEN
		
		IF (PATCHES(2)==1) THEN
			LEFT_PT = BETA - DELTA
			RIGHT_PT = 1.0D0
		ELSEIF (PATCHES(2)==2) THEN
			LEFT_PT = ALPHA - DELTA
			RIGHT_PT = BETA + DELTA
		ELSEIF (PATCHES(2)==3) THEN
			LEFT_PT = 0.0D0
			RIGHT_PT = ALPHA + DELTA
		ENDIF
		
		UNITSQ_TO_PARSP%VAL(0) = (RIGHT_PT - LEFT_PT)*REFPT + LEFT_PT
		UNITSQ_TO_PARSP%VAL(1) = (RIGHT_PT - LEFT_PT)
		
	ENDIF
	
END FUNCTION UNITSQ_TO_PARSP

!----------------------------------------------------------------------------------------
TYPE(FUNCTION_1D) FUNCTION PARSP_TO_UNITSQ(PHYPT, PATCHES, AXIS)

	REAL*8, INTENT(IN) :: PHYPT
	INTEGER, INTENT(IN) :: PATCHES(2)
	CHARACTER(LEN=1), INTENT(IN) :: AXIS
	
	REAL*8 :: LEFT_PT, RIGHT_PT
	IF (AXIS=='X') THEN
		
		PARSP_TO_UNITSQ%VAL(0) = PHYPT
		PARSP_TO_UNITSQ%VAL(1) = 1.0D0
		
	ELSEIF (AXIS=='Y') THEN
		
		IF (PATCHES(2)==1) THEN
			LEFT_PT = BETA - DELTA
			RIGHT_PT = 1.0D0
		ELSEIF (PATCHES(2)==2) THEN
			LEFT_PT = ALPHA - DELTA
			RIGHT_PT = BETA + DELTA
		ELSEIF (PATCHES(2)==3) THEN
			LEFT_PT = 0.0D0
			RIGHT_PT = ALPHA + DELTA
		ENDIF
		
		PARSP_TO_UNITSQ%VAL(0) = (PHYPT - LEFT_PT)/(RIGHT_PT - LEFT_PT)
		PARSP_TO_UNITSQ%VAL(1) = 1.0D0/(RIGHT_PT - LEFT_PT)
		
	ENDIF
	
END FUNCTION PARSP_TO_UNITSQ

TYPE(FUNCTION_1D) FUNCTION PHYPATCH1_TO_UNITSQ(PHYPT, AXIS)

	REAL*8, INTENT(IN) :: PHYPT
	CHARACTER(LEN=1), INTENT(IN) :: AXIS
	
	IF (AXIS=='X') THEN
		PHYPATCH1_TO_UNITSQ%VAL(0) = (1.0D0/(2.0D0*(1.0D0 - BETA + DELTA)))*(PHYPT + (1.0D0 - BETA + DELTA))
		PHYPATCH1_TO_UNITSQ%VAL(1) = (1.0D0/(2.0D0*(1.0D0 - BETA + DELTA)))
	ELSEIF (AXIS=='Y') THEN
		PHYPATCH1_TO_UNITSQ%VAL(0) = (1.0D0/(2.0D0*(1.0D0 - BETA + DELTA)))*(PHYPT + (1.0D0 - BETA + DELTA))
		PHYPATCH1_TO_UNITSQ%VAL(1) = (1.0D0/(2.0D0*(1.0D0 - BETA + DELTA)))
	ENDIF

END FUNCTION PHYPATCH1_TO_UNITSQ

END MODULE PATCH_MAPPING
