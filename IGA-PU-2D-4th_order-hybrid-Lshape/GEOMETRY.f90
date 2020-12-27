MODULE GEOMETRY

	USE PATCH_MAPPING

! 	IMPLICIT INTEGER (I-N)
! 	IMPLICIT REAL(8) (A-H,O-Z)
    IMPLICIT NONE
CONTAINS

SUBROUTINE GET_GEO()
	
	INTEGER :: I, J, II, JJ, K, KK, DIR, PATCH, PATCH_ROW, PATCH_COLUMN, PATCHES(2)
	INTEGER :: BASIS_POLY_ORDER, OLD_BASIS_POLY_ORDER, BASIS_NUM_KNOTS, OLD_BASIS_NUM_KNOTS, TMP_BASIS_NUM_KNOTS
	INTEGER :: BASIS_MULTIPLICITIES(MAX_LENGTH), OLD_BASIS_MULTIPLICITIES(MAX_LENGTH), TMP_BASIS_MULTIPLICITIES(MAX_LENGTH)
	REAL(8) :: BASIS_KNOTS(MAX_LENGTH), OLD_BASIS_KNOTS(MAX_LENGTH), TMP_BASIS_KNOTS(MAX_LENGTH), NEW_KNOTS(0:MAX_LENGTH), TMP_KNOTS(0:MAX_LENGTH)
	
	CHARACTER(LEN=1) :: AXIS
	LOGICAL, ALLOCATABLE :: BD_MASK(:)
	INTEGER :: EGM_ETA_ORDER
	REAL*8 :: INNER_RADIUS, INNER_LAYER
	
	TYPE(POINT2D) :: PHYPT, PARPT
	INTEGER :: TMPK
	REAL*8 :: XI, ETA
	
	!! POLYGONAL L-SHAPED DOMAIN WITH ONE PATCH
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	GEO_KVEC(1)%POLY_ORDER = 2
	GEO_KVEC(2)%POLY_ORDER = 2

	GEO_KVEC(1)%LENGTH = 6
	GEO_KVEC(2)%LENGTH = 5

	GEO_KVEC(1)%KNOTS(0:GEO_KVEC(1)%LENGTH) = (/0.D0,0.D0,0.D0,0.50D0,1.0D0,1.0D0,1.0D0/)
	GEO_KVEC(2)%KNOTS(0:GEO_KVEC(2)%LENGTH) = (/0.D0,0.D0,0.D0,1.0D0,1.0D0,1.0D0/)
	
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

!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!
! 	j = 2
!     call KNOT_TO_ARRAY(basis_kvec(1, j), OLD_BASIS_POLY_ORDER, OLD_BASIS_NUM_KNOTS, OLD_BASIS_KNOTS, OLD_BASIS_MULTIPLICITIES)
!     K = 0
!     do i = 1, OLD_BASIS_NUM_KNOTS
!         if (OLD_BASIS_KNOTS(i)>INV_DELTA(1) .and. OLD_BASIS_KNOTS(i)<=1.0d0) then 
!         ELSE 
!             K = K + 1
!             TMP_BASIS_KNOTS(K) = OLD_BASIS_KNOTS(I)
!             TMP_BASIS_MULTIPLICITIES(K) = OLD_BASIS_MULTIPLICITIES(I)
!         endif
!     enddo
!     TMP_BASIS_NUM_KNOTS = K
! 	call KNOT_TO_ARRAY(basis_kvec(2, j), OLD_BASIS_POLY_ORDER, OLD_BASIS_NUM_KNOTS, OLD_BASIS_KNOTS, OLD_BASIS_MULTIPLICITIES)
! 	
! 	DO I = 1, OLD_BASIS_NUM_KNOTS
!         IF (OLD_BASIS_KNOTS(I)>0.0D0 .AND. OLD_BASIS_KNOTS(I)<INV_DELTA(2)) THEN 
!             PHYPT = GET_PHY_PT(POINT2D(1.0D0, OLD_BASIS_KNOTS(I)), 2)
!             PARPT = GET_PAR_PT(PHYPT, 1)
!             TMP_BASIS_NUM_KNOTS = TMP_BASIS_NUM_KNOTS + 1
!             TMP_BASIS_KNOTS(TMP_BASIS_NUM_KNOTS) = PARPT%Y
!             TMP_BASIS_MULTIPLICITIES(TMP_BASIS_NUM_KNOTS) = 1
!         ENDIF
!     ENDDO
!     TMP_BASIS_NUM_KNOTS = TMP_BASIS_NUM_KNOTS + 1
!     TMP_BASIS_KNOTS(TMP_BASIS_NUM_KNOTS) = 1.0D0
!     TMP_BASIS_MULTIPLICITIES(TMP_BASIS_NUM_KNOTS) = BASIS_KVEC(1, J)%POLY_ORDER + 1
!     
!     BASIS_KVEC(1, J) = ARRAY_TO_KNOT(BASIS_KVEC(1, J)%POLY_ORDER, TMP_BASIS_NUM_KNOTS, TMP_BASIS_KNOTS, TMP_BASIS_MULTIPLICITIES)
!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!

!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!
    j = 1
    TMP_BASIS_KNOTS(1) = 0.0D0
    TMP_BASIS_MULTIPLICITIES(1) = BASIS_KVEC(1, J)%POLY_ORDER + 1
    TMP_BASIS_NUM_KNOTS = 1
    
!     call KNOT_TO_ARRAY(basis_kvec(1, j), OLD_BASIS_POLY_ORDER, OLD_BASIS_NUM_KNOTS, OLD_BASIS_KNOTS, OLD_BASIS_MULTIPLICITIES)
!     K = 0
!     do i = 1, OLD_BASIS_NUM_KNOTS
!         if (OLD_BASIS_KNOTS(i)>0.0D0 .and. OLD_BASIS_KNOTS(i)<=1.0d0) then 
!         ELSE 
!             K = K + 1
!             TMP_BASIS_KNOTS(K) = OLD_BASIS_KNOTS(I)
!             TMP_BASIS_MULTIPLICITIES(K) = OLD_BASIS_MULTIPLICITIES(I)
!         endif
!     enddo
!     TMP_BASIS_NUM_KNOTS = K

	call KNOT_TO_ARRAY(basis_kvec(2, j), OLD_BASIS_POLY_ORDER, OLD_BASIS_NUM_KNOTS, OLD_BASIS_KNOTS, OLD_BASIS_MULTIPLICITIES)
	
	DO I = 1, OLD_BASIS_NUM_KNOTS
        IF (OLD_BASIS_KNOTS(I)>0.0D0 .AND. OLD_BASIS_KNOTS(I)<1.0D0) THEN 
            PHYPT = GET_PHY_PT(POINT2D(OLD_BASIS_KNOTS(I), 0.0D0), 2)
            PARPT = GET_PAR_PT(PHYPT, 1)
!             print*, parpt%x
            TMP_BASIS_NUM_KNOTS = TMP_BASIS_NUM_KNOTS + 1
            TMP_BASIS_KNOTS(TMP_BASIS_NUM_KNOTS) = PARPT%X
            TMP_BASIS_MULTIPLICITIES(TMP_BASIS_NUM_KNOTS) = 1
        ENDIF
    ENDDO
    TMP_BASIS_NUM_KNOTS = TMP_BASIS_NUM_KNOTS + 1
    TMP_BASIS_KNOTS(TMP_BASIS_NUM_KNOTS) = 1.0D0
    TMP_BASIS_MULTIPLICITIES(TMP_BASIS_NUM_KNOTS) = BASIS_KVEC(1, J)%POLY_ORDER + 1
    
    BASIS_KVEC(1, J) = ARRAY_TO_KNOT(BASIS_KVEC(1, J)%POLY_ORDER, TMP_BASIS_NUM_KNOTS, TMP_BASIS_KNOTS, TMP_BASIS_MULTIPLICITIES)
!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!

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
        ! xi direction
        J = 1
		IF (PATCH==2) THEN        ! regular mapping G
            DO K = 1, 5
                XI = (1.0D0*K)/6.0D0
                DO I = 1, NUMIR(PATCH, PATCH, J)
                    IF (DABS(IR_GRID(PATCH, PATCH, J, I) - XI)<=EPS) THEN 
                        GOTO 125
                    ENDIF
                ENDDO
                !--------------------------------------------------------------------------------!
                TMPIR_GRID(:, :, :, :) = IR_GRID(:, :, :, :)
                TMP_NUMIR(:, :, :) = NUMIR(:, :, :)
                !--------------------------------------------------------------------------------!
                DO I = 1, NUMIR(PATCH, PATCH, J) - 1
                    IF (XI > IR_GRID(PATCH, PATCH, J ,I) .AND. XI < IR_GRID(PATCH, PATCH, J, I + 1)) THEN 
                        IR_GRID(PATCH, PATCH, J, I + 1) = XI
                        TMPK = I
                        GOTO 128
                    ENDIF
                ENDDO
                PRINT*, 'ERROR: GEOMETRY - XI  IS NOT INCLUDED IN THE IR_GRID'
                STOP
                
                128 CONTINUE
                !--------------------------------------------------------------------------------!
                DO I = TMPK + 2, NUMIR(PATCH, PATCH, J) + 1
                    IR_GRID(PATCH, PATCH, J, I) = TMPIR_GRID(PATCH, PATCH, J, I - 1)
                ENDDO
                NUMIR(PATCH, PATCH, J) = NUMIR(PATCH, PATCH, J) + 1
                125 CONTINUE
            ENDDO
		ENDIF
		! eta direction
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
    DO I = 1, NUMIR(2, 2, J)
        IR_GRID(2, 1, J, I) = IR_GRID(2, 2, J, I)
    ENDDO   
    NUMIR(2, 1, J) = NUMIR(2, 2, J)
    
!     DO I = 1, NUMIR(1, 1, J)
!         !--------------------------------------------------------------------------------!
!         TMPIR_GRID(2, 1, J, :) = IR_GRID(2, 1, J, :)
!         TMP_NUMIR(2, 1, J) = NUMIR(2, 1, J)
!         !--------------------------------------------------------------------------------!
!         DO II = 1, NUMIR(2, 1, J)-1
!             PHYPT = GET_PHY_PT(POINT2D(IR_GRID(1, 1, J, I), 1.0D0), 1)
!             PARPT = GET_PAR_PT(PHYPT, 2)
!             IF (PARPT%X>IR_GRID(2, 1, J, II) .AND. PARPT%X<IR_GRID(2, 1, J, II+1)) THEN
!                 IR_GRID(2, 1, J, II + 1) = PARPT%X
!                 TMPK = II
!                 GOTO 217
!             ENDIF
!         ENDDO
!         
!         GOTO 218
!         
!         217 CONTINUE
!         !--------------------------------------------------------------------------------!
!         DO II = TMPK + 2, NUMIR(2, 1, J) + 1
!             IR_GRID(2, 1, J, II) = TMPIR_GRID(2, 1, J, II - 1)
!         ENDDO
!         NUMIR(2, 1, J) = NUMIR(2, 1, J) + 1
!         !--------------------------------------------------------------------------------!
!         
!         218 CONTINUE
!         
!     ENDDO
    
    ! eta direction
    J = 2
    K = 0
    DO I = 1, NUMIR(2, 2, J)
        IF (IR_GRID(2, 2, J, I)<=INV_DELTA(2)) THEN  
            K = K + 1
            IR_GRID(2, 1, J, K) = IR_GRID(2, 2, J, I)
        ENDIF
    ENDDO
    NUMIR(2, 1, J) = K
    !--------------------------------------------------------------------------------!
!     K = 0
!     ! regular mapping G
!     DO II = 2, NUMIR(1, 1, J) - 1
!         IF (IR_GRID(1, 1, J, II)>=INV_DELTA(1)) THEN 
! !             PRINT*, 'II',II, IR_GRID(2,2,J,II)
!             PHYPT = GET_PHY_PT(POINT2D(1.0d0, IR_GRID(1, 1, J, II)), 1)
!             PARPT = GET_PAR_PT(PHYPT, 2)
! !             print*, phypt, parpt
!             K = K + 1
!             INVIR_GRID(1, 1, J, K) = PARPT%Y
!         ENDIF
!     ENDDO
!     INV_NUMIR(1, 1, J) = K
! !     PRINT*, INV_NUMIR(2,2,J), INVIR_GRID(2,2,J,1)
! 	!--------------------------------------------------------------------------------!
!     DO II = 1, INV_NUMIR(1, 1, J)
!         DO I = 1, NUMIR(2, 1, J)
!             IF (DABS(IR_GRID(2, 1, J, I) - INVIR_GRID(1, 1, J, II))<=EPS) THEN 
!                 GOTO 135
!             ENDIF
!         ENDDO
!         !--------------------------------------------------------------------------------!
!         TMPIR_GRID(2, 1, J, :) = IR_GRID(2, 1, J, :)
!         TMP_NUMIR(2, 1, J) = NUMIR(2, 1, J)
!         !--------------------------------------------------------------------------------!
!         DO I = 1, NUMIR(2, 1, J) - 1
!             IF (INVIR_GRID(1, 1, J, II) > IR_GRID(2, 1, J, I) .AND. INVIR_GRID(1, 1, J, II) < IR_GRID(2, 1, J, I + 1)) THEN 
!                 IR_GRID(2, 1, J, I + 1) = INVIR_GRID(1, 1, J, II)
!                 TMPK = I
!                 GOTO 138
!             ENDIF
!         ENDDO
!         PRINT*, 'ERROR: GEOMETRY - INVIR_GRID  IS NOT INCLUDED IN THE IR_GRID'
!         STOP
!         
!         138 CONTINUE
!         !--------------------------------------------------------------------------------!
!         DO I = TMPK + 2, NUMIR(2, 1, J) + 1
!             IR_GRID(2, 1, J, I) = TMPIR_GRID(2, 1, J, I - 1)
!         ENDDO
!         NUMIR(2, 1, J) = NUMIR(2, 1, J) + 1
!         135 CONTINUE
!     ENDDO
!     !--------------------------------------------------------------------------------!
    DO J = 1, 2
        IR_GRID(1, 2, J, :) = IR_GRID(2, 1, J, :)
        NUMIR(1, 2, J) = NUMIR(2, 1, J)
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
	
	if (problem==1) then ! Fourth order
        DO I = 1, DOF
            IF (NDX(1, I)==1 .AND. (NDX(2, I)==0 .OR. NDX(2, I)==1 .OR. NDX(2, I)==LOC_NUMBS(1, 1)-2 .OR. NDX(2, I)==LOC_NUMBS(1, 1)-1 .OR. NDX(3, I)==0 .OR. NDX(3, I)==1)) THEN   ! singular mapping F
                K = K + 1
                ZERO_BDNDX(K)%LC_NDX(1:3) = NDX(:,I)
                ZERO_BDNDX(K)%GL_NDX = I
                BD_DOF = BD_DOF + 1
            ELSEIF (NDX(1, I)==2 .AND. (NDX(2, I)==0 .OR. NDX(2, I)==1 .OR. NDX(2, I)==LOC_NUMBS(2, 1)-2 .OR. NDX(2, I)==LOC_NUMBS(2, 1)-1 .OR. NDX(3, I)==LOC_NUMBS(2, 2)-2 .OR. NDX(3, I)==LOC_NUMBS(2, 2)-1)) THEN   ! regular mapping G
                K = K + 1
                ZERO_BDNDX(K)%LC_NDX(1:3) = NDX(:,I)
                ZERO_BDNDX(K)%GL_NDX = I
                BD_DOF = BD_DOF + 1
            ENDIF
        ENDDO
    elseif (problem==2 .or. problem==4) then ! Eliptic
        DO I = 1, DOF
            IF (NDX(1, I)==1 .AND. (NDX(2, I)==0 .OR. NDX(2, I)==LOC_NUMBS(1, 1)-1 .OR. NDX(3, I)==0)) THEN   ! singular mapping F
                K = K + 1
                ZERO_BDNDX(K)%LC_NDX(1:3) = NDX(:,I)
                ZERO_BDNDX(K)%GL_NDX = I
                BD_DOF = BD_DOF + 1
            ELSEIF (NDX(1, I)==2 .AND. (NDX(2, I)==0 .OR. NDX(2, I)==LOC_NUMBS(2, 1)-1 .OR. NDX(3, I)==LOC_NUMBS(2, 2)-1)) THEN   ! regular mapping G
                K = K + 1
                ZERO_BDNDX(K)%LC_NDX(1:3) = NDX(:,I)
                ZERO_BDNDX(K)%GL_NDX = I
                BD_DOF = BD_DOF + 1
            ENDIF
        ENDDO
    endif
    
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

