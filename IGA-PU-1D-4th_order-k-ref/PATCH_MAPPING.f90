MODULE PATCH_MAPPING

	USE NURBS

	IMPLICIT NONE

CONTAINS

! 	!!  GET THE CORRESPONDING POINT ON THE PHYSICAL SPACE OF GIVE POINT ON THE PARAMETRIC SPACE
! 	TYPE(POINT2D) FUNCTION GET_PHY_PT(REF_PT)
! 
! 		TYPE(POINT2D), INTENT(IN) :: REF_PT
! 		INTEGER :: ITER
! 
! ! 		IF (INI_GEO==1) THEN
! ! 			GET_PHY_PT = GET_POINT_NURVE_SURFACE_2D(REF_PT,GEO_KVEC,GEO_CTL)
! ! 		ELSE
! ! 			GET_PHY_PT = GET_POINT_NURVE_SURFACE_2D(REF_PT,BASIS_KVEC,BASIS_CTL)
! ! 		ENDIF
! 
! 		GET_PHY_PT%X = REF_PT%Y**2*(DCOS(2.0D0*PI*(1.0D0-REF_PT%X)))
! 		GET_PHY_PT%Y = REF_PT%Y**2*(DSIN(2.0D0*PI*(1.0D0-REF_PT%X)))
! 		
! 	END FUNCTION GET_PHY_PT
! 	
! 	
! 	!!  GET THE CORRESPONDING POINT ON THE PARAMETRIC SPACE OF GIVE POINT ON THE PHYSICAL SPACE
! 	TYPE(POINT2D) FUNCTION GET_REF_PT(PHY_PT)
! 
! 		TYPE(POINT2D), INTENT(IN) :: PHY_PT
! 
! 		TYPE(FUNCTION_2D) :: DIFF_NURBS(2)
! 		TYPE(POINT2D) :: NEW_PT, OLD_PT, TMP_PT
! 		REAL(8) :: ERROR, INV_J
! 		INTEGER :: ITER
! 
! 		ERROR = 1.0D0
! 		OLD_PT = POINT2D(0.D0,0.D0)
! 		ITER = 0
! 		DO WHILE (ERROR>TOLERANCE)
! 			ITER = ITER + 1
! ! 			PRINT*, 'OLD_PT : ', OLD_PT
! 			CALL GET_DIFF_NURVE_SURFACE_2D(DIFF_NURBS,OLD_PT,GEO_KVEC,GEO_CTL,1)
! 			TMP_PT%X = DIFF_NURBS(1)%VAL(0,0)-PHY_PT%X
! 			TMP_PT%Y = DIFF_NURBS(2)%VAL(0,0)-PHY_PT%Y
! 			ERROR = MAX(DABS(TMP_PT%X),DABS(TMP_PT%Y))
! 			IF (ERROR<TOLERANCE) THEN
! 				GET_REF_PT = OLD_PT
! 				GOTO 999
! 			ELSE
! 				INV_J = 1.0D0 / (DIFF_NURBS(1)%VAL(1,0)*DIFF_NURBS(2)%VAL(0,1)-DIFF_NURBS(1)%VAL(0,1)*DIFF_NURBS(2)%VAL(1,0))
! 				NEW_PT%X = OLD_PT%X - INV_J*(DIFF_NURBS(2)%VAL(0,1)*TMP_PT%X-DIFF_NURBS(1)%VAL(0,1)*TMP_PT%Y)
! 				NEW_PT%Y = OLD_PT%Y - INV_J*(DIFF_NURBS(1)%VAL(1,0)*TMP_PT%Y-DIFF_NURBS(2)%VAL(1,0)*TMP_PT%X)
! 				OLD_PT = NEW_PT
! 			ENDIF
! 		ENDDO
! 		IF (ITER==MAX_ITERATION) THEN
! 			PRINT *, '[ERROR]  CANNOT FIND REFERENCE POINT !'
! 		ENDIF
! 
! 		999 CONTINUE
! 
! 	END FUNCTION GET_REF_PT
! 
! 
! 	!!  GET JACOBIAN OF THE GEOMETRIC MAPPING
! 	TYPE(MATRIX_22) FUNCTION GET_JACOBIAN_MATRIX(REF_PT)
! 
! 		TYPE(POINT2D), INTENT(IN) :: REF_PT
! 
! 		TYPE(FUNCTION_2D) :: DIFF_NURBS(2)
! 
! 		integer :: ix, iy, k
! 
! ! 		IF (INI_GEO==1) THEN
! ! 			CALL GET_DIFF_NURVE_SURFACE_2D(DIFF_NURBS,REF_PT,GEO_KVEC,GEO_CTL,1)
! ! 		ELSE
! ! 			CALL GET_DIFF_NURVE_SURFACE_2D(DIFF_NURBS,REF_PT,BASIS_KVEC,BASIS_CTL,1)
! ! 		ENDIF
! 
! ! 		GET_JACOBIAN_MATRIX%ENT(1,1) = DIFF_NURBS(1)%VAL(1,0)
! ! 		GET_JACOBIAN_MATRIX%ENT(1,2) = DIFF_NURBS(1)%VAL(0,1)
! ! 		GET_JACOBIAN_MATRIX%ENT(2,1) = DIFF_NURBS(2)%VAL(1,0)
! ! 		GET_JACOBIAN_MATRIX%ENT(2,2) = DIFF_NURBS(2)%VAL(0,1)
! 
! 		GET_JACOBIAN_MATRIX%ENT(1,1) = 2.0D0*PI*REF_PT%Y**2.0D0*DSIN(2.0D0*PI*(1.0D0 - REF_PT%X))
! 		GET_JACOBIAN_MATRIX%ENT(1,2) = 2.0D0*REF_PT%Y*DCOS(2.0D0*PI*(1.0D0 - REF_PT%X))
! 		GET_JACOBIAN_MATRIX%ENT(2,1) = -2.0D0*PI*REF_PT%Y**2.0D0*DCOS(2.0D0*PI*(1.0D0 - REF_PT%X))
! 		GET_JACOBIAN_MATRIX%ENT(2,2) = 2.0D0*REF_PT%Y*DSIN(2.0D0*PI*(1.0D0 - REF_PT%X))
! 		
! 	END FUNCTION GET_JACOBIAN_MATRIX
! 	
! 	!!  GET SECOND-ORDER PARTIAL DERIVATIVES
! 	TYPE(SECOND_PARTIAL_DERIVATIVES) FUNCTION GET_PARTIAL(REF_PT)
! 
! 		TYPE(POINT2D), INTENT(IN) :: REF_PT
! 
! 		TYPE(FUNCTION_2D) :: DIFF_NURBS(2)
! 
! ! 		IF (INI_GEO==1) THEN
! ! 			CALL GET_DIFF_NURVE_SURFACE_2D(DIFF_NURBS,REF_PT,GEO_KVEC,GEO_CTL,2)
! ! 		ELSE
! ! 			CALL GET_DIFF_NURVE_SURFACE_2D(DIFF_NURBS,REF_PT,BASIS_KVEC,BASIS_CTL,2)
! ! 		ENDIF
! 
! ! 		GET_PARTIAL%X(1) = DIFF_NURBS(1)%VAL(2,0)
! ! 		GET_PARTIAL%X(2) = DIFF_NURBS(1)%VAL(1,1)
! ! 		GET_PARTIAL%X(3) = DIFF_NURBS(1)%VAL(0,2)
! ! 
! ! 		GET_PARTIAL%Y(1) = DIFF_NURBS(2)%VAL(2,0)
! ! 		GET_PARTIAL%Y(2) = DIFF_NURBS(2)%VAL(1,1)
! ! 		GET_PARTIAL%Y(3) = DIFF_NURBS(2)%VAL(0,2)
! 
! 		GET_PARTIAL%X(1) = -4.0D0*PI**2*REF_PT%Y**2*DCOS(2.0D0*PI*(1.0D0 - REF_PT%X))
! 		GET_PARTIAL%X(2) = 4.0D0*PI*REF_PT%Y*DSIN(2.0D0*PI*(1.0D0 - REF_PT%X))
! 		GET_PARTIAL%X(3) = 2.0D0*DCOS(2.0D0*PI*(1.0D0 - REF_PT%X))
! 
! 		GET_PARTIAL%Y(1) = -4.0D0*PI**2*REF_PT%Y**2*DSIN(2.0D0*PI*(1.0D0 - REF_PT%X))
! 		GET_PARTIAL%Y(2) = -4.0D0*PI*REF_PT%Y*DCOS(2.0D0*PI*(1.0D0 - REF_PT%X))
! 		GET_PARTIAL%Y(3) = 2.0D0*DSIN(2.0D0*PI*(1.0D0 - REF_PT%X))
! 
! 	END FUNCTION GET_PARTIAL
! 
! 	REAL*8 FUNCTION GET_DET_LN(REF_PT, GAMMA)
! 	
! 		TYPE(POINT2D), INTENT(IN) :: REF_PT
! 		INTEGER, INTENT(IN) :: GAMMA
! 	
! 		TYPE(FUNCTION_2D) :: DIFF_NURBS(2)
! 	
! 		IF (INI_GEO==1) THEN
! 			CALL GET_DIFF_NURVE_SURFACE_2D(DIFF_NURBS,REF_PT,GEO_KVEC(:),GEO_CTL,1)
! 		ELSE
! 			CALL GET_DIFF_NURVE_SURFACE_2D(DIFF_NURBS,REF_PT,BASIS_KVEC(:),BASIS_CTL,1)
! 		ENDIF
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
! 		IF (INI_GEO==1) THEN
! 			CALL GET_DIFF_NURVE_SURFACE_2D(DIFF_NURBS,REF_PT,GEO_KVEC(:),GEO_CTL,1)
! 		ELSE
! 			CALL GET_DIFF_NURVE_SURFACE_2D(DIFF_NURBS,REF_PT,BASIS_KVEC(:),BASIS_CTL,1)
! 		ENDIF
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! NURBS PU WITH FLAT-TOP

! MAPPING LOCAL PATCH TO [-1, 1] WITHOUT expanding a patch with delta
! TYPE(FUNCTION_1D) FUNCTION PATCH_TO_REF1D(X, PATCH)
! 
! 	REAL*8, INTENT(IN) :: X
! 	INTEGER, INTENT(IN) :: PATCH
! 	
! 		PATCH_TO_REF1D%VAL(0) = (2.0D0/(PATCHBDPT(PATCH,2) - PATCHBDPT(PATCH,1)))*(X - PATCHBDPT(PATCH,1)) - 1.0D0
! 		PATCH_TO_REF1D%VAL(1) = (2.0D0/(PATCHBDPT(PATCH,2) - PATCHBDPT(PATCH,1)))
! 		PATCH_TO_REF1D%VAL(2) = 0.0D0
! 	
! END FUNCTION PATCH_TO_REF1D

!-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!

! MAPPING FROM [0, 1] TO [0, 0.5]
TYPE(FUNCTION_1D) FUNCTION MAP_F(REFPT)

	REAL*8, INTENT(IN) :: REFPT
	
! 	MAP_F%VAL(0) = 0.50D0*REFPT**(5.0D0)
! 	MAP_F%VAL(1) = 0.50D0*5.0D0*REFPT**(4.0D0)
! 	MAP_F%VAL(2) = 0.50D0*5.0D0*4.0D0*REFPT**(3.0D0)
	
	MAP_F%VAL(0) = 0.50D0*REFPT
	MAP_F%VAL(1) = 0.50D0
	MAP_F%VAL(2) = 0.0D0
	
	IF (MAP_F%VAL(0)<0.0D0 .OR. MAP_F%VAL(0)>0.50D0) THEN 
		PRINT*, 'ERROR: MAP_F - the physical point is out of the range'
		STOP 
	ENDIF
	
END FUNCTION MAP_F

TYPE(FUNCTION_1D) FUNCTION MAP_INVF(PHYPT)

  REAL*8, INTENT(IN) :: PHYPT
  
!   MAP_INVF%VAL(0) = (2.0D0*PHYPT)**(0.20D0)
!   MAP_INVF%VAL(1) = 0.20D0*(2.0D0*PHYPT)**(-0.80D0)*2.0D0
!   MAP_INVF%VAL(2) = -0.80D0*0.20D0*(2.0D0*PHYPT)**(-1.80D0)*4.0D0
  
  MAP_INVF%VAL(0) = (2.0D0*PHYPT)
  MAP_INVF%VAL(1) = 2.0D0
  MAP_INVF%VAL(2) = 0.0D0
  
END FUNCTION MAP_INVF
!-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!

! MAPPING FROM [0, 1] TO [0.4, 1]
TYPE(FUNCTION_1D) FUNCTION MAP_G(REFPT)
	
	REAL*8, INTENT(IN) :: REFPT
	
	MAP_G%VAL(0) = 0.60D0*REFPT + 0.40D0
	MAP_G%VAL(1) = 0.60D0
	MAP_G%VAL(2) = 0.0D0

	IF (MAP_G%VAL(0)<0.40D0 .OR. MAP_G%VAL(0)>1.0D0) THEN 
		PRINT*, 'ERROR:MAP_G - the physical point is the out of the range'
		STOP
	ENDIF
	
END FUNCTION MAP_G

!-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!

! INVERSE MAPPING FROM [0.4, 1] TO [0, 1]
TYPE(FUNCTION_1D) FUNCTION MAP_INVG(PHYPT)
	
	REAL*8, INTENT(IN) :: PHYPT
	
	MAP_INVG%VAL(0) = (PHYPT - 0.40D0)/0.60D0
	MAP_INVG%VAL(1) = 10.0D0/6.0D0
	MAP_INVG%VAL(2) = 0.0D0
	
	IF (MAP_INVG%VAL(0)<0.0D0 .OR. MAP_INVG%VAL(0)>1.0D0) THEN 
		PRINT*, 'ERROR:MAP_INVG - the reference point is the out of the domain'
		STOP
	ENDIF
	
END FUNCTION MAP_INVG

!-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!

! AFFINE MAPPING FROM [0.8, 1] TO [0, 1]
TYPE(FUNCTION_1D) FUNCTION MAP_AF(REFPT)

	REAL*8, INTENT(IN) :: REFPT
	
	IF (REFPT<0.80D0) THEN 
		PRINT*, "ERROR : the reference point is not in the domain of the affine mapping T_2"
		PRINT*, "Reference point : ", REFPT
		STOP
	ENDIF
	
	MAP_AF%VAL(0) = 5.0D0*REFPT - 4.0D0
	MAP_AF%VAL(1) = 5.0D0
	MAP_AF%VAL(2) = 0.0D0

END FUNCTION MAP_AF

!-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!

! MAPPING LOCAL PATCH TO [-1, 1]
TYPE(FUNCTION_1D) FUNCTION PATCH_TO_REF1D(X, PATCH)

	REAL*8, INTENT(IN) :: X
	INTEGER, INTENT(IN) :: PATCH
		
		IF (PATCH==1) THEN
			PATCH_TO_REF1D%VAL(0) = (2.0D0/((PATCHBDPT(PATCH,2) + DELTA) - PATCHBDPT(PATCH,1)))*(X - PATCHBDPT(PATCH,1)) - 1.0D0
			PATCH_TO_REF1D%VAL(1) = (2.0D0/((PATCHBDPT(PATCH,2) + DELTA) - PATCHBDPT(PATCH,1)))
			PATCH_TO_REF1D%VAL(2) = 0.0D0
		ELSEIF (PATCH==NUMPATCH) THEN
			PATCH_TO_REF1D%VAL(0) = (2.0D0/(PATCHBDPT(PATCH,2) - (PATCHBDPT(PATCH,1) - DELTA)))*(X - (PATCHBDPT(PATCH,1) - DELTA)) - 1.0D0
			PATCH_TO_REF1D%VAL(1) = (2.0D0/(PATCHBDPT(PATCH,2) - (PATCHBDPT(PATCH,1) - DELTA)))
			PATCH_TO_REF1D%VAL(2) = 0.0D0
		ELSE 
			PATCH_TO_REF1D%VAL(0) = (2.0D0/((PATCHBDPT(PATCH,2) + DELTA) - (PATCHBDPT(PATCH,1) - DELTA)))*(X - (PATCHBDPT(PATCH,1) - DELTA)) - 1.0D0
			PATCH_TO_REF1D%VAL(1) = (2.0D0/((PATCHBDPT(PATCH,2) + DELTA) - (PATCHBDPT(PATCH,1) - DELTA)))
			PATCH_TO_REF1D%VAL(2) = 0.0D0
		ENDIF
	
END FUNCTION PATCH_TO_REF1D

!-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!

! MAPPING [-1, 1] TO LOCAL PATCH
TYPE(FUNCTION_1D) FUNCTION REF_TO_LOCPATCH1D(REFX, LOCPATCH)

	REAL*8, INTENT(IN) :: REFX
	INTEGER, INTENT(IN) :: LOCPATCH
	
	REF_TO_LOCPATCH1D%VAL(0) = 0.10D0*(REFX + 2.0D0*LOCPATCH - 1.0D0)
	REF_TO_LOCPATCH1D%VAL(1) = 0.10D0
	REF_TO_LOCPATCH1D%VAL(2) = 0.0D0
	
END FUNCTION REF_TO_LOCPATCH1D

!-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!

!!-------------------[[[[[ PATCH + 2*DELTA TO [-1, 1] OF PU DOMAIN ]]]]] ------------------------
!! CAUTIOn: Specifically modified to adjust PU-IGA 4th order problem 1D !!

TYPE(FUNCTION_1D) FUNCTION PATCH_TO_REFPU(PHYX, PATCH) ! PHYX is a point in the physical domain

    REAL*8, INTENT(IN) :: PHYX
	INTEGER, INTENT(IN) :: PATCH
    
    REAL*8 :: XPTS(4)
    TYPE(FUNCTION_1D) :: PARPT, PHYPT
			
		IF (PATCH==1) THEN
			
			XPTS(1) = 0.0D0 - 2.0D0*DELTA
			XPTS(2) = 0.0D0
			XPTS(3) = 0.40D0
			XPTS(4) = 0.50D0
			
! 			XPTS(1) = PATCHBDPT(1,1) - 2.0D0*DELTA
! 			XPTS(2) = PATCHBDPT(1,1)
! 			XPTS(3) = PATCHBDPT(1,2) - DELTA
! 			XPTS(4) = PATCHBDPT(1,2) + DELTA
			
! 			print*, xpts(3), xpts(4)
! 			stop
			IF (PHYX>=XPTS(2) .AND. PHYX<=XPTS(3)) THEN
				PATCH_TO_REFPU%VAL(0) = 0.0D0
				PATCH_TO_REFPU%VAL(1) = 0.0D0
			ELSEIF (PHYX>=XPTS(3) .AND. PHYX<=XPTS(4)) THEN
				PATCH_TO_REFPU%VAL(0) = (PHYX-XPTS(3))/(XPTS(4)-XPTS(3))
				PATCH_TO_REFPU%VAL(1) = 1.0D0/(XPTS(4)-XPTS(3))
			ELSE 
				PATCH_TO_REFPU%VAL(0) = 5.0D0
				PATCH_TO_REFPU%VAL(1) = 0.0D0
			ENDIF
		ELSEIF (PATCH==2) THEN
		!------------------------------------------------------------------------------!
! 			XPTS(1) = PATCHBDPT(PATCH,1) - DELTA
! 			XPTS(2) = PATCHBDPT(PATCH,1) + DELTA
! 			XPTS(3) = PATCHBDPT(PATCH,2)
! 			XPTS(4) = PATCHBDPT(PATCH,2) + 2.0D0*DELTA
		!------------------------------------------------------------------------------!	

			XPTS(1) = 0.40D0
			XPTS(2) = 0.50D0
			XPTS(3) = 1.0D0
			XPTS(4) = 1.0D0 + 2.0D0*DELTA
			
		!------------------------------------------------------------------------------!				
			IF (PHYX>=XPTS(1) .AND. PHYX<=XPTS(2)) THEN
				PATCH_TO_REFPU%VAL(0) = (PHYX-XPTS(2))/(XPTS(2)-XPTS(1))
				PATCH_TO_REFPU%VAL(1) = 1.0D0/(XPTS(2)-XPTS(1))
			ELSEIF (PHYX>=XPTS(2) .AND. PHYX<=XPTS(3)) THEN
				PATCH_TO_REFPU%VAL(0) = 0.0D0
				PATCH_TO_REFPU%VAL(1) = 0.0D0
			ELSE 
				PATCH_TO_REFPU%VAL(0) = 5.0D0
				PATCH_TO_REFPU%VAL(1) = 0.0D0
			ENDIF
! 		ELSE
! 		
! 			XPTS(1) = PATCHBDPT(PATCH,1) - DELTA
! 			XPTS(2) = PATCHBDPT(PATCH,1) + DELTA
! 			XPTS(3) = PATCHBDPT(PATCH,2) - DELTA
! 			XPTS(4) = PATCHBDPT(PATCH,2) + DELTA
! 			
! 			IF (PHYX>=XPTS(1) .AND. PHYX<=XPTS(2)) THEN
! 				PATCH_TO_REFPU%VAL(0) = (PHYX-XPTS(2))/(XPTS(2)-XPTS(1))
! 				PATCH_TO_REFPU%VAL(1) = 1.0D0/(XPTS(2)-XPTS(1))
! 			ELSEIF (PHYX>=XPTS(2) .AND. PHYX<=XPTS(3)) THEN
! 				PATCH_TO_REFPU%VAL(0) = 0.0D0
! 				PATCH_TO_REFPU%VAL(1) = 0.0D0
! 			ELSEIF (PHYX>=XPTS(3) .AND. PHYX<=XPTS(4)) THEN
! 				PATCH_TO_REFPU%VAL(0) = (PHYX-XPTS(3))/(XPTS(4)-XPTS(3))
! 				PATCH_TO_REFPU%VAL(1) = 1.0D0/(XPTS(4)-XPTS(3))
! 			ELSE 
! 				PATCH_TO_REFPU%VAL(0) = 5.0D0
! 				PATCH_TO_REFPU%VAL(1) = 0.0D0
! 			ENDIF
		ENDIF

!     IF ((DABS(PHYX-XPTS(2)).LE.EPS) .OR. (DABS(PHYX-XPTS(3)).LE.EPS)) THEN
!       PATCH_TO_REFPU%VAL(0) = 0.0D0
!       PATCH_TO_REFPU%VAL(1) = 0.0D0
!     ENDIF
! 
!     IF ((DABS(PHYX-XPTS(1)).LE.EPS) .OR. (DABS(PHYX-XPTS(4)).LE.EPS)) THEN
! 	 PATCH_TO_REFPU%VAL(0) = 0.0D0
!       PATCH_TO_REFPU%VAL(1) = 0.0D0
!     ENDIF
	
	PATCH_TO_REFPU%VAL(2) = 0.0D0
	
END FUNCTION PATCH_TO_REFPU

!-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!

! MAPPING [0, 1] TO PATCH
!----------------------------------------------------------------------------------------
TYPE(FUNCTION_1D) FUNCTION UNITLINE_TO_PATCH(REFPT, PATCH)
	
	REAL*8, INTENT(IN) :: REFPT
	INTEGER, INTENT(IN) :: PATCH
	
	REAL*8 :: LEFT_PT, RIGHT_PT

	IF (NUMPATCH==1) THEN 
		UNITLINE_TO_PATCH%VAL(0) = REFPT
		UNITLINE_TO_PATCH%VAL(1) = 1.0D0
		UNITLINE_TO_PATCH%VAL(2) = 0.0D0
		GOTO 76
	ENDIF
	
	IF (PATCH==NUMPATCH) THEN
		LEFT_PT = PATCHBDPT(PATCH,1) - DELTA
		RIGHT_PT = PATCHBDPT(PATCH,2)
	ELSEIF (PATCH==1) THEN
		LEFT_PT = PATCHBDPT(PATCH,1)
		RIGHT_PT = PATCHBDPT(PATCH,2) + DELTA
	ELSE
		LEFT_PT = PATCHBDPT(PATCH,1) - DELTA
		RIGHT_PT = PATCHBDPT(PATCH,2) + DELTA
	ENDIF
	
! 	LEFT_PT = PATCHBDPT(PATCH, 1)
! 	RIGHT_PT = PATCHBDPT(PATCH, 2)
	
	UNITLINE_TO_PATCH%VAL(0) = (RIGHT_PT - LEFT_PT)*REFPT + LEFT_PT
	UNITLINE_TO_PATCH%VAL(1) = (RIGHT_PT - LEFT_PT)
	UNITLINE_TO_PATCH%VAL(2) = 0.0D0
	
	76 CONTINUE
	
END FUNCTION UNITLINE_TO_PATCH

!-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!

! MAPPING PATCH TO [0, 1]
!----------------------------------------------------------------------------------------
TYPE(FUNCTION_1D) FUNCTION PATCH_TO_UNITLINE(PHYPT, PATCH)

	REAL*8, INTENT(IN) :: PHYPT
	INTEGER, INTENT(IN) :: PATCH
	
	REAL*8 :: LEFT_PT, RIGHT_PT
	
	IF (NUMPATCH==1) THEN 
		PATCH_TO_UNITLINE%VAL(0) = PHYPT
		PATCH_TO_UNITLINE%VAL(1) = 1.0D0
		PATCH_TO_UNITLINE%VAL(2) = 0.0D0
		GOTO 76
	ENDIF
	
	IF (PATCH==NUMPATCH) THEN
		LEFT_PT = PATCHBDPT(PATCH,1) - DELTA
		RIGHT_PT = PATCHBDPT(PATCH,2)
	ELSEIF (PATCH==1) THEN
		LEFT_PT = PATCHBDPT(PATCH,1)
		RIGHT_PT = PATCHBDPT(PATCH,2) + DELTA
	ELSE
		LEFT_PT = PATCHBDPT(PATCH,1) - DELTA
		RIGHT_PT = PATCHBDPT(PATCH,2) + DELTA
	ENDIF
	
! 	LEFT_PT = PATCHBDPT(PATCH, 1)
! 	RIGHT_PT = PATCHBDPT(PATCH, 2)
	
	PATCH_TO_UNITLINE%VAL(0) = (PHYPT - LEFT_PT)/(RIGHT_PT - LEFT_PT)
	PATCH_TO_UNITLINE%VAL(1) = 1.0D0/(RIGHT_PT - LEFT_PT)
	PATCH_TO_UNITLINE%VAL(2) = 0.0D0
	
	76 CONTINUE
	
END FUNCTION PATCH_TO_UNITLINE

END MODULE PATCH_MAPPING

