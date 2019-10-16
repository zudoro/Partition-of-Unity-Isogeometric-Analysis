MODULE GEOMETRY

	USE GSQUAD
	USE PATCH_MAPPING

! 	IMPLICIT INTEGER (I-N)
! 	IMPLICIT REAL(8) (A-H,O-Z)

CONTAINS

SUBROUTINE GET_GEO()
	
	INTEGER :: I, J, II, JJ, K, KK, DIR, PATCH_ROW, PATCH_COLUMN
	INTEGER :: BASIS_POLY_ORDER, OLD_BASIS_POLY_ORDER, BASIS_NUM_KNOTS, OLD_BASIS_NUM_KNOTS
	INTEGER :: BASIS_MULTIPLICITIES(MAX_LENGTH), OLD_BASIS_MULTIPLICITIES(MAX_LENGTH)
	REAL(8) :: BASIS_KNOTS(MAX_LENGTH), OLD_BASIS_KNOTS(MAX_LENGTH), NEW_KNOTS(0:MAX_LENGTH)
	CHARACTER(LEN=1) :: AXIS

	TYPE(FUNCTION_1D) :: PARPT, PHYPT
	INTEGER :: TMPK, EXTRA_GSPT
	REAL*8 :: FX, GX
	
! Unit circle with crack singularity along positive x-axis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	PU_KVEC%POLY_ORDER = PUORDER
	
! 	PU_KVEC%LENGTH = 2*(NUMPATCH - 1)*(PU_KVEC%POLY_ORDER - 1) + 2*(PU_KVEC%POLY_ORDER + 1) - 1
	
! 	PU_KVEC%KNOTS(0:PU_KVEC%LENGTH) = (/0.D0, 0.D0, 0.D0, 0.0D0, 0.20D0, 0.20D0, 0.40D0, 0.40D0, 0.60D0, 0.60D0, 0.80D0, 0.80D0, 1.D0, 1.0D0, 1.0D0, 1.0D0/)
	
! 	PU_KVEC = GET_OPEN_KNOT_VECTOR(PUORDER)
	
! 	PU_KVEC = UNIFORM_KNOT_INSERTION(PU_KVEC,2*(NUMPATCH - 1), 1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 	WRITE(*,*)
! 	WRITE(*,*) '<<< SET PU KNOT : DONE >>>'
! 	WRITE(*,*)

!! ----------------------- [[[[[ SET KNOT INFO OF B-SPLINE BASIS FUNCTIONS IN THE PATCH FOR EACH ]]]]] --------------------------
!-----------------------------------------------------------------------------------------------------------
!! BASIS FUNCTION ON THE OMEGA_HAT_F
	BS_KVEC(1)%POLY_ORDER = 2
	
	BS_KVEC(1)%LENGTH = 5
	
	BS_KVEC(1)%KNOTS(0:BS_KVEC(1)%LENGTH) = (/0.D0,0.D0,0.D0,1.D0,1.D0,1.D0/)

!! BASIS FUNCTION ON THE OMEGA_HAT_G
	BS_KVEC(2)%POLY_ORDER = 2
	
	BS_KVEC(2)%LENGTH = 5
	
	BS_KVEC(2)%KNOTS(0:BS_KVEC(2)%LENGTH) = (/0.D0,0.D0,0.D0,1.D0,1.D0,1.D0/)
!-----------------------------------------------------------------------------------------------------------


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! [[[[[ REFINEMENT ]]]]] !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!------------------------------------- [[[[[ B-SPLINE ]]]]] -----------------------------------------------
!-----------------------------------------------------------------------------------------------------------
! ELEVATE DEGREE OF B-SPLINE BASIS FUNCTION (p-refinement)
	DO K = 1, NUMPATCH
		IF (BSORDER(K)>BS_KVEC(K)%POLY_ORDER) THEN
			DO I = 1, BSORDER(K) - BS_KVEC(K)%POLY_ORDER
				BS_KVEC(K) = DEGREE_ELEVATION(BS_KVEC(K))
			ENDDO
		ENDIF
	ENDDO
!-----------------------------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------------------------
! INSERT NEW KNOTS IN THE KNOTS CORRESPONDING TO B-SPLINE BASIS FUNCTIONS (h-refinement)
	DO K = 1, NUMPATCH
		IF (EXTRA_KNOTS(K)>0) THEN
			BS_KVEC(K) = UNIFORM_KNOT_INSERTION(BS_KVEC(K), EXTRA_KNOTS(K))
		ENDIF
	ENDDO
!-----------------------------------------------------------------------------------------------------------


!!------------------------------------- [[[[[ NURBS ]]]]] -----------------------------------------------
!-----------------------------------------------------------------------------------------------------------
! ELEVATE DEGREE OF BERNSTEIN POLY. AND FIND NEW CONTROL POINTS

! 	IF (BSORDER(1).GT.GEO_KVEC(1)%POLY_ORDER) THEN
! 		DO I=1, BSORDER(1) - GEO_KVEC(1)%POLY_ORDER
! 			CALL P_REFINEMENT(BASIS_CTL, BS_KVEC(1:2)%POLY_ORDER, 'X')
! 			BS_KVEC(1) = DEGREE_ELEVATION(BS_KVEC(1))
! 		ENDDO
! 	ENDIF
! 
! 	IF (BSORDER(2).GT.GEO_KVEC(2)%POLY_ORDER) THEN
! 		DO I=1, BSORDER(2) - GEO_KVEC(2)%POLY_ORDER
! 			CALL P_REFINEMENT(BASIS_CTL, BS_KVEC(1:2)%POLY_ORDER, 'Y')
! 			BS_KVEC(2) = DEGREE_ELEVATION(BS_KVEC(2))
! 		ENDDO
! 	ENDIF
!-----------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------
! ELEVATE DEGREE OF NURBS BASIS FUNCTIONS
! 	DO PATCH_COLUMN = 1, NUMPATCH(2)
! 		IF (BSORDER(1,PATCH_COLUMN,1).GT.BS_KVEC(1,PATCH_COLUMN,1)%POLY_ORDER) THEN
! 			CALL ELEVATE_DEGREE(BS_KVEC(1,PATCH_COLUMN,:), BS_CTL(1,PATCH_COLUMN), BSORDER(1,PATCH_COLUMN,1)-BS_KVEC(1,PATCH_COLUMN,1)%POLY_ORDER, 'X')
! 		ENDIF
! 		IF (BSORDER(1,PATCH_COLUMN,2).GT.BS_KVEC(1,PATCH_COLUMN,2)%POLY_ORDER) THEN
! 			CALL ELEVATE_DEGREE(BS_KVEC(1,PATCH_COLUMN,:), BS_CTL(1,PATCH_COLUMN), BSORDER(1,PATCH_COLUMN,2)-BS_KVEC(1,PATCH_COLUMN,2)%POLY_ORDER, 'Y')
! 		ENDIF
! 	ENDDO
	
!-----------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------
! REFINE KNOT VECTOR INSERTING NEW KNOTS WITH MULTIPLICITIES ARE GREATER THAN 1
	
! 	DO PATCH_COLUMN = 1, NUMPATCH(2)
! 		OLD_BS_KVEC(:) = BS_KVEC(1,PATCH_COLUMN,:)
! 		DO DIR = 1, 2
! 	! 		PRINT*, 'DIR', DIR
! 			IF (EXTRA_KNOTS(1,PATCH_COLUMN,DIR).GT.0) THEN
! 				IF (MESHTYPE=='SH' .AND. PATCH_COLUMN==3 .AND. DIR==2) THEN
! 					!! ----------------------- [[[[[ Shishkin type k-refinement ]]]]] ------------------------------
! 					DO I = 1, EXTRA_KNOTS(1, PATCH_COLUMN, DIR)
! 						K = BS_KVEC(1, PATCH_COLUMN, DIR)%LENGTH - (BS_KVEC(1, PATCH_COLUMN, DIR)%POLY_ORDER + 1) + I
! 						BS_KVEC(1, PATCH_COLUMN, DIR)%KNOTS(K) = &
! 						1.0D0 - (EXTRA_KNOTS(1,PATCH_COLUMN,DIR) + 1  - I)*(KAPPA*BS_KVEC(1, PATCH_COLUMN, DIR)%POLY_ORDER/(1.0D0*EXTRA_KNOTS(1,PATCH_COLUMN,DIR)))
! 					ENDDO
! 					DO I = 1, BS_KVEC(1, PATCH_COLUMN, DIR)%POLY_ORDER + 1
! 						K = BS_KVEC(1, PATCH_COLUMN, DIR)%LENGTH - (BS_KVEC(1, PATCH_COLUMN, DIR)%POLY_ORDER + 1) + EXTRA_KNOTS(1, PATCH_COLUMN, DIR) + I
! 						BS_KVEC(1, PATCH_COLUMN, DIR)%KNOTS(K) = 1.0D0
! 					ENDDO
! 					BS_KVEC(1, PATCH_COLUMN, DIR)%LENGTH = BS_KVEC(1, PATCH_COLUMN, DIR)%LENGTH + EXTRA_KNOTS(1,PATCH_COLUMN,DIR)
! 					
! 					DO I = 0, BS_KVEC(1, PATCH_COLUMN, DIR)%LENGTH
! 						TMP_KNOTS(BS_KVEC(1, PATCH_COLUMN, DIR)%LENGTH - I) = 1.0D0 - BS_KVEC(1, PATCH_COLUMN, DIR)%KNOTS(I)
! 					ENDDO
! 					BS_KVEC(1, PATCH_COLUMN, DIR)%KNOTS(:) = TMP_KNOTS(:)
! 					!-----------------------------------------------------------------------------------------------------
! 				ELSE
! 					BS_KVEC(1,PATCH_COLUMN,DIR) = UNIFORM_KNOT_INSERTION(BS_KVEC(1,PATCH_COLUMN,DIR), EXTRA_KNOTS(1,PATCH_COLUMN,DIR),NURBS_REGULARITY(1,PATCH_COLUMN,DIR))
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
! ! 						BS_KVEC(DIR)%KNOTS(BS_KVEC(DIR)%POLY_ORDER + I) = (I*1.D0/((EXTRA_KNOTS(2)+1)*1.D0))**5
! ! 					ENDDO
! ! 				ELSEIF (DIR==2 .AND. (PATCH==3)) THEN
! ! 					DO I=1, INT(EXTRA_KNOTS(2))
! ! 						IGA_KVEC(PATCH,DIR)%KNOTS(IGA_KVEC(PATCH,DIR)%POLY_ORDER + I) = 1.D0 - ((EXTRA_KNOTS(2)-I+1)*1.D0/((EXTRA_KNOTS(2)+1)*1.D0))**5
! ! 					ENDDO
! ! 				ENDIF
! 				!-----------------------------------------------------------------------------------------------------
! 
! 				CALL KNOT_TO_ARRAY(OLD_BS_KVEC(DIR), OLD_BASIS_POLY_ORDER, OLD_BASIS_NUM_KNOTS, OLD_BASIS_KNOTS, OLD_BASIS_MULTIPLICITIES)
! 				CALL KNOT_TO_ARRAY(BS_KVEC(1,PATCH_COLUMN,DIR), BASIS_POLY_ORDER, BASIS_NUM_KNOTS, BASIS_KNOTS, BASIS_MULTIPLICITIES)
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
! 				CALL REFINE_KNOT(OLD_BS_KVEC(:), BS_CTL(1,PATCH_COLUMN), NEW_KNOTS(0:K), K, AXIS)
! 			ENDIF	
! 		ENDDO
! 	ENDDO
! 	
! ! 	PRINT*, 'AFTER REFINEMENT', BSORDER(1,2,2), BS_KVEC(1,2,2)%POLY_ORDER
!-----------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------
! REFINE KNOT VECTOR INSERTING NEW KNOTS WITH MULTIPLICITIES 1
! 	IF (EXTRA_KNOTS(1).GT.0) THEN
! 		DO I = 1, EXTRA_KNOTS(1)
! 			NEW_KNOTS(1,I-1) = I*1.D0/(1.D0*(EXTRA_KNOTS(1)+1))
! 		ENDDO
! 		CALL REFINE_KNOT(BS_KVEC(:), BASIS_CTL, NEW_KNOTS(1,0:EXTRA_KNOTS(1)-1), EXTRA_KNOTS(1)-1, 'X')
! 	ENDIF
! 	
! 	IF (EXTRA_KNOTS(2).GT.0) THEN
! 		DO J = 1, EXTRA_KNOTS(2)
! 			NEW_KNOTS(2,J-1) = J*1.D0/(1.D0*(EXTRA_KNOTS(2)+1))
! 		ENDDO
! 		CALL REFINE_KNOT(BS_KVEC(:), BASIS_CTL, NEW_KNOTS(2,0:EXTRA_KNOTS(2)-1), EXTRA_KNOTS(2)-1, 'Y')
! 	ENDIF
!-----------------------------------------------------------------------------------------------------------

! TEST REFINE KNOT
! 	CALL REFINE_KNOT(BS_KVEC(:), BASIS_CTL, NEW_KNOTS(0:K), K, 'X')

!-----------------------------------------------------------------------------------------------------------
! INSERT A NEW KNOT INTO A EXISTING KNOT VECTOR

! 	BS_KVEC(1) = UNIFORM_KNOT_INSERTION(BS_KVEC(1), EXTRA_KNOTS(1))
! 	CALL H_REFINEMENT(BASIS_CTL, BS_KVEC(:), 0.5D0, 'X')
! 
! 	BS_KVEC(2) = UNIFORM_KNOT_INSERTION(BS_KVEC(2), EXTRA_KNOTS(2))
! 	CALL H_REFINEMENT(BASIS_CTL, BS_KVEC(:), 0.5D0, 'Y')
!-----------------------------------------------------------------------------------------------------------
	

	WRITE(*,*)
	WRITE(*,*) '<<< SET BASIS KNOT : DONE >>>'
	WRITE(*,*)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
!!------------------------------- [[[[[ COMPUTE DOF ]]]]] -----------------------------------	
! 	NUMBS(:) = (/BS_KVEC(1)%LENGTH - BSORDER(1)-1, BS_KVEC(2)%LENGTH - BSORDER(2)/) ! IMPOSED PERIODIC CONDTION
	
	NUMBS(:) = 0
	DO K = 1, NUMPATCH
		NUMBS(K) = BS_KVEC(K)%LENGTH - BSORDER(K)
	ENDDO
	
! 	numbs(1) = 4

	DOF = SUM(NUMBS(:))
	
! 	ALLOCATE(ENRICH_ORDER(NUMBS(1)))
! 	!-----------------------------------------------------------------!
! 	! Set the degrees of the enrichment functions !
! 	ENRICH_ORDER(:) = (/8.0D0, 9.0D0, 10.0D0, 15.0D0/)
	!-----------------------------------------------------------------!
	
	WRITE(*,*)
	WRITE(*,*) '<<< COMPUTE DOF : DONE >>>'
	WRITE(*,*)

!! ------------------------------------ [[[[[ GLOBAL INDEX ]]]]] -------------------------------------------
!-----------------------------------------------------------------------------------------------------------
	! PU-IGA
	ALLOCATE(NDX(2, DOF))
	
	I = 1
	DO PATCH = 1, NUMPATCH
        DO J = 1, NUMBS(PATCH)
            NDX(1,I) = PATCH
            NDX(2,I) = J-1
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

!!-----------------------------[[[[[ GLOBAL & LOCAL INDEX ON THE BOUNDARY ]]]]]-------------------------
!----------------------------------------------------------------------------------------------------------------------------

	!! SET INDEX OF BASIS FUNCTIONS IMPOSED ESSENTIAL BOUNDARY CONDITION
	K = 0
	BD_DOF = 0
	DO I = 1, DOF
! 		IF (NDX(1,I)==0 .OR. NDX(1,I)==NUMBS(1)-1 .OR. NDX(2,I)==0 .OR. NDX(2,I)==NUMBS(2)-1) THEN
		IF (NDX(1,I)==1) THEN
! 			IF (NDX(2,I)==0 .OR. NDX(2,I)==1) THEN
! 				K = K + 1
! 				BDNDX(K)%LC_NDX(:) = NDX(:,I)
! 				BDNDX(K)%GL_NDX = I
! 				BD_DOF = BD_DOF + 1
! 			ENDIF
		ELSEIF (NDX(1,I)==2) THEN
			IF (NDX(2,I)==NUMBS(2)-2 .OR. NDX(2,I)==NUMBS(2)-1) THEN
				K = K + 1
				BDNDX(K)%LC_NDX(:) = NDX(:,I)
				BDNDX(K)%GL_NDX = I
				BD_DOF = BD_DOF + 1
			ENDIF
		ENDIF
	ENDDO
	
	BDNDX(:)%LC_NUM = K

! 	ALLOCATE(LST_SOL(K), BD_COL_PT(K))

! 	BD_DOF = 2
! 	BDNDX(1)%LC_NDX(:) = NDX(:, 1)
! 	BDNDX(1)%GL_NDX = 1
! 	BDNDX(2)%LC_NDX(:) = NDX(:, DOF)
! 	BDNDX(2)%GL_NDX = DOF
	
	WRITE(*,*)
	WRITE(*,*) '<<< SET LOCAL INDEX CORRESPONDING BASIS FUNCTIONS ON BOUNDARY IMPOSED ESSENTIAL BOUNDARY CONDITION : DONE >>>'
	WRITE(*,*)
	
!--------------------------------- [[[[[ INTEGRAL REGIONS ON THE PARAMETER SPACE ]]]]] --------------------------------
!! CONSTRUCT THE FRAME OF LOCAL INTEGRAL REGIONS ON THE PARAMETER SPACE

!! DIVIDE PARAMETER SPACE INTO SUB-INTERVALS WHOSE END POINTS ARE KNOTS
	IR_GRID(:,:,:) = 0.D0
	TMPIR_GRID(:,:,:) = 0.0D0
	
	PATCHBDPT = RESHAPE((/0.0D0, 0.40D0 + DELTA, 0.40D0 + DELTA, 1.0D0/), (/2, NUMPATCH/))
	PARPT = MAP_INVF(PATCHBDPT(1,2) - DELTA)
	FX = PARPT%VAL(0)
	PARPT = MAP_INVG(PATCHBDPT(2,1) + DELTA)
	GX = PARPT%VAL(0)
	
	DO PATCH = 1, NUMPATCH
		K = 1
! 		PARPT = UNITLINE_TO_PATCH(BS_KVEC(PATCH)%KNOTS(0), PATCH)
! 		IR_GRID(PATCH, PATCH, K) = PARPT%VAL(0)
        IR_GRID(PATCH, PATCH, K) = BS_KVEC(PATCH)%KNOTS(0)
		DO I = 1, BS_KVEC(PATCH)%LENGTH
			K = K + 1
			IF (DABS(BS_KVEC(PATCH)%KNOTS(I-1) - BS_KVEC(PATCH)%KNOTS(I)).LE.EPS) THEN
				K = K - 1
			ELSE
! 				PARPT = UNITLINE_TO_PATCH(BS_KVEC(PATCH)%KNOTS(I), PATCH)
! 				IR_GRID(PATCH, PATCH, K) = PARPT%VAL(0)
                IR_GRID(PATCH, PATCH, K) = BS_KVEC(PATCH)%KNOTS(I)
			ENDIF
		ENDDO
		NUMIR(PATCH, PATCH) = K
	
!-----------------------------------------------------------------------------------------------------------------------------------------------!	
		IF (PATCH==1) THEN 
			DO I = 1, NUMIR(PATCH, PATCH)
				IF (DABS(IR_GRID(PATCH, PATCH, I) - FX)<=EPS) THEN 
					GOTO 115
				ENDIF
			ENDDO
			!--------------------------------------------------------------------------------!
			TMPIR_GRID(:, :, :) = IR_GRID(:, :, :)
			TMP_NUMIR(:, :) = NUMIR(:, :)
			!--------------------------------------------------------------------------------!
			DO I = 1, NUMIR(PATCH, PATCH) - 1
				IF (FX > IR_GRID(PATCH, PATCH, I) .AND. FX < IR_GRID(PATCH, PATCH, I + 1)) THEN 
					IR_GRID(PATCH, PATCH, I + 1) = FX
					TMPK = I
					GOTO 116
				ENDIF
			ENDDO
			PRINT*, 'ERROR: GEOMETRY - (PATCHBDPT(1, 2) - DELTA)  IS NOT INCLUDED IN THE IR_GRID'
			STOP
			
			116 CONTINUE
			!--------------------------------------------------------------------------------!
			DO I = TMPK + 2, NUMIR(PATCH, PATCH) + 1
				IR_GRID(PATCH, PATCH, I) = TMPIR_GRID(PATCH, PATCH, I - 1)
			ENDDO
			NUMIR(PATCH, PATCH) = NUMIR(PATCH, PATCH) + 1
        !--------------------------------------------------------------------------------!
		ELSEIF (PATCH==NUMPATCH) THEN
			DO I = 1, NUMIR(PATCH, PATCH)
				IF (DABS(IR_GRID(PATCH, PATCH, I) - GX)<=EPS) THEN 
					GOTO 115
				ENDIF
			ENDDO
			!--------------------------------------------------------------------------------!
			TMPIR_GRID(:, :, :) = IR_GRID(:, :, :)
			TMP_NUMIR(:, :) = NUMIR(:, :)
			!--------------------------------------------------------------------------------!
			DO I = 1, NUMIR(PATCH, PATCH) - 1
				IF (GX > IR_GRID(PATCH, PATCH, I) .AND. GX < IR_GRID(PATCH, PATCH, I + 1)) THEN 
					IR_GRID(PATCH, PATCH, I + 1) = GX
					TMPK = I
					GOTO 118
				ENDIF
			ENDDO
			PRINT*, 'ERROR: GEOMETRY - 0.87818 AND (PATCHBDPT(2, 1) + DELTA)  IS NOT INCLUDED IN THE IR_GRID'
			STOP
			
			118 CONTINUE
			!--------------------------------------------------------------------------------!
			DO I = TMPK + 2, NUMIR(PATCH, PATCH) + 1
				IR_GRID(PATCH, PATCH, I) = TMPIR_GRID(PATCH, PATCH, I - 1)
			ENDDO
			NUMIR(PATCH, PATCH) = NUMIR(PATCH, PATCH) + 1
		ENDIF
		115 CONTINUE
	ENDDO
	!--------------------------------------------------------------------------------!
    K = 0
    DO I = 1, NUMIR(1,1)
        IF (IR_GRID(1, 1, I)>=FX) THEN  
            K = K + 1
            IR_GRID(1, 2, K) = IR_GRID(1, 1, I)
        ENDIF
    ENDDO
    NUMIR(1, 2) = K
    !--------------------------------------------------------------------------------!
    k = 0
    DO J = 2, NUMIR(2, 2) - 1
        IF (IR_GRID(2, 2, J)<=GX) THEN 
            PHYPT = MAP_G(IR_GRID(2, 2, J))
            PARPT = MAP_INVF(PHYPT%VAL(0))
            k = k + 1
            INVIR_GRID(2, 2, k) = PARPT%VAL(0)
        ENDIF
    ENDDO
    INV_NUMIR(2, 2) = k
	!--------------------------------------------------------------------------------!
    DO J = 1, INV_NUMIR(2, 2)
        DO I = 1, NUMIR(1, 2)
            IF (DABS(IR_GRID(1, 2, I) - INVIR_GRID(2, 2, J))<=EPS) THEN 
                GOTO 135
            ENDIF
        ENDDO
        !--------------------------------------------------------------------------------!
        TMPIR_GRID(1, 2, :) = IR_GRID(1, 2, :)
        TMP_NUMIR(1, 2) = NUMIR(1, 2)
        !--------------------------------------------------------------------------------!
        DO I = 1, NUMIR(1, 2) - 1
            IF (INVIR_GRID(2, 2, J) > IR_GRID(1, 2, I) .AND. INVIR_GRID(2, 2, J) < IR_GRID(1, 2, I + 1)) THEN 
                IR_GRID(1, 2, I + 1) = INVIR_GRID(2, 2, J)
                TMPK = I
                GOTO 138
            ENDIF
        ENDDO
        PRINT*, 'ERROR: GEOMETRY - INVIR_GRID  IS NOT INCLUDED IN THE IR_GRID'
        STOP
        
        138 CONTINUE
        !--------------------------------------------------------------------------------!
        DO I = TMPK + 2, NUMIR(1, 2) + 1
            IR_GRID(1, 2, I) = TMPIR_GRID(1, 2, I - 1)
        ENDDO
        NUMIR(1, 2) = NUMIR(1, 2) + 1
        135 CONTINUE
    ENDDO
        !--------------------------------------------------------------------------------!
        IR_GRID(2, 1, :) = IR_GRID(1, 2, :)
        NUMIR(2, 1) = NUMIR(1, 2)
!-----------------------------------------------------------------------------------------------------------------------------------------------!	

	!! GENERATE GAUSS POINTS
	EXTRA_GSPT = 15
	NUMGSPT = INT(0.5*(BSORDER(1) + PUORDER)) + EXTRA_GSPT
	DO I = 2, NUMPATCH
		IF (( INT(0.5*(BSORDER(I) + PUORDER)) + EXTRA_GSPT) > NUMGSPT) THEN
			NUMGSPT = INT(0.5*(BSORDER(I) + PUORDER)) + EXTRA_GSPT
		ENDIF
	ENDDO
	
	ALLOCATE(GSPT(NUMGSPT), GSW(NUMGSPT))
	
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

SUBROUTINE  FIND_GREVILLE_PT(GREV_PTS, KVEC)
	TYPE(KNOT_VECTOR), INTENT(IN) :: KVEC
	REAL*8, INTENT(OUT) :: GREV_PTS(KVEC%LENGTH - KVEC%POLY_ORDER)
	
	REAL*8 :: SUM1
	INTEGER :: I, J
	
	GREV_PTS(:) = 0.0D0
				
	DO I = 0, UBOUND(GREV_PTS,1) - 1
		SUM1 = 0.0D0
		DO J = 1, KVEC%POLY_ORDER
			SUM1 = SUM1 + KVEC%KNOTS(I+J)
		END DO
		SUM1 = SUM1/ KVEC%POLY_ORDER
		GREV_PTS(I+1) = SUM1
	ENDDO
	
END SUBROUTINE FIND_GREVILLE_PT

END MODULE GEOMETRY
