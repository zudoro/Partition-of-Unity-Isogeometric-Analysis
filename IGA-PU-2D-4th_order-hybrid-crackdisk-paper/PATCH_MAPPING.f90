MODULE PATCH_MAPPING

	USE NURBS

	implicit integer (i-n)
	implicit real(8) (a-h,o-z)

CONTAINS

	!!  GET THE CORRESPONDING POINT ON THE PHYSICAL SPACE OF GIVE POINT ON THE PARAMETRIC SPACE
	TYPE(POINT2D) FUNCTION GET_PHY_PT(PAR_PT, PATCH)

		TYPE(POINT2D), INTENT(IN) :: PAR_PT
		INTEGER, INTENT(IN) :: PATCH
		REAL*8 :: FACTOR
		
        IF (PDE=='FOUR' .AND. PROBLEM<=3) THEN 
            IF (PATCH==1) THEN ! singular mapping F(xi, eta)
                GET_PHY_PT%X = 2.0D0*PAR_PT%X - 1.0D0
                GET_PHY_PT%Y = 1.250D0*PAR_PT%Y - 1.0D0
            ELSEIF (PATCH==2) THEN ! regular mapping G(xi, eta)
                GET_PHY_PT%X = 2.0D0*PAR_PT%X - 1.0D0
                GET_PHY_PT%Y = 1.250D0*PAR_PT%Y - 0.250D0
            ENDIF
        ELSEIF ((PDE=='ELPT' .AND. PROBLEM==3) .OR. (PDE=='FOUR' .AND. PROBLEM>=4)) THEN
            IF (PATCH==1) THEN ! singular mapping F(xi, eta)
                !-----------------------------------------------------------!
                FACTOR = RADI_SING(2)*PAR_PT%Y**(2)
!                 FACTOR = RADI_SING(2)*PAR_PT%Y
!                 FACTOR = RADI_SING(2)*PAR_PT%Y**(4)
                !-----------------------------------------------------------!
                GET_PHY_PT%X = FACTOR*DCOS(2.0D0*PI*(1.0D0 - PAR_PT%X))
                GET_PHY_PT%Y = FACTOR*DSIN(2.0D0*PI*(1.0D0 - PAR_PT%X))
            ELSEIF (PATCH==2) THEN ! regular mapping G(xi, eta)
                FACTOR = RADI_SING(1) + (1.0D0 - RADI_SING(1))*PAR_PT%Y
                GET_PHY_PT%X = FACTOR*DCOS(2.0D0*PI*(1.0D0 - PAR_PT%X))
                GET_PHY_PT%Y = FACTOR*DSIN(2.0D0*PI*(1.0D0 - PAR_PT%X))
            ENDIF
        ENDIF
		
	END FUNCTION GET_PHY_PT

    TYPE(POINT2D) FUNCTION GET_PAR_PT(PHY_PT, PATCH)
        
        TYPE(POINT2D), INTENT(IN) :: PHY_PT
        INTEGER, INTENT(IN) :: PATCH
        REAL*8 :: R, ARCOS
        
        R = DSQRT(PHY_PT%X**2 + PHY_PT%Y**2)
        ARCOS = DACOS(PHY_PT%X/R)
        
        IF (PDE=='FOUR' .AND. PROBLEM<=3) THEN 
            IF (PATCH==1) THEN      ! singular mapping F
                GET_PAR_PT%X = 0.50D0*(PHY_PT%X + 1.0D0)
                GET_PAR_PT%Y = (PHY_PT%Y + 1.0D0)/1.250D0
            ELSEIF (PATCH==2) THEN      ! reguar mapping G
                GET_PAR_PT%X = 0.50D0*(PHY_PT%X + 1.0D0)
                GET_PAR_PT%Y = (PHY_PT%Y + 0.250D0)/1.250D0
            ENDIF
        ELSEIF ((PDE=='ELPT' .AND. PROBLEM==3) .OR. (PDE=='FOUR' .AND. PROBLEM>=4)) THEN 
            IF (PATCH==1) THEN      ! singular mapping F
                IF (DABS(R)<=EPS) THEN 
                    GET_PAR_PT%X = 0.0D0
                    GET_PAR_PT%Y = 0.0D0
                    GOTO 111
                ELSEIF (DABS(PHY_PT%Y)<=EPS .AND. PHY_PT%X>EPS) THEN 
                    GET_PAR_PT%X = 1.0D0
                ELSEIF (DABS(PHY_PT%Y)<=EPS .AND. PHY_PT%X<EPS) THEN
                    GET_PAR_PT%X = 0.50D0
                ELSEIF (PHY_PT%Y<EPS) THEN 
                    GET_PAR_PT%X = (1.0D0/(2.0D0*PI))*ARCOS
                ELSEIF (PHY_PT%Y>EPS) THEN
                    GET_PAR_PT%X = 1.0D0 - (1.0D0/(2.0D0*PI))*ARCOS
                ENDIF
                !-----------------------------------------------------------!
                GET_PAR_PT%Y = DSQRT(R/RADI_SING(2))
!                 GET_PAR_PT%Y = R/RADI_SING(2)
!                 GET_PAR_PT%Y = (R/RADI_SING(2))**(0.250D0)
                !-----------------------------------------------------------!
            ELSEIF (PATCH==2) THEN      ! reguar mapping G
                IF (DABS(R)<=EPS) THEN 
                    GET_PAR_PT%X = 0.0D0
                    GET_PAR_PT%Y = 0.0D0
                    GOTO 111
                ELSEIF (DABS(PHY_PT%Y)<=EPS .AND. PHY_PT%X>EPS) THEN 
                    GET_PAR_PT%X = 1.0D0
                ELSEIF (DABS(PHY_PT%Y)<=EPS .AND. PHY_PT%X<EPS) THEN
                    GET_PAR_PT%X = 0.50D0
                ELSEIF (PHY_PT%Y<EPS) THEN 
                    GET_PAR_PT%X = (1.0D0/(2.0D0*PI))*ARCOS
                ELSEIF (PHY_PT%Y>EPS) THEN
                    GET_PAR_PT%X = 1.0D0 - (1.0D0/(2.0D0*PI))*ARCOS
                ENDIF
                GET_PAR_PT%Y = (R - RADI_SING(1))/(1.0D0 - RADI_SING(1))
            ENDIF
        ENDIF
        
        111 continue
        
    END FUNCTION

	!!  GET JACOBIAN OF THE GEOMETRIC MAPPING
	TYPE(MATRIX_22) FUNCTION GET_JACOBIAN_MATRIX(PAR_PT, PATCH)

		TYPE(POINT2D), INTENT(IN) :: PAR_PT
        INTEGER, INTENT(IN) :: PATCH
        
		TYPE(FUNCTION_2D) :: DIFF_NURBS(2)
		
!		TYPE(MATRIX_22) :: TEST_MATRIX1, TEST_MATRIX2, test_result
        
        IF (PDE=='FOUR' .AND. PROBLEM<=3) THEN 
            IF (PATCH==1) THEN ! singular mapping F(xi, eta)
                GET_JACOBIAN_MATRIX%ENT(1,1) = 2.0D0 !DX/DXI
                GET_JACOBIAN_MATRIX%ENT(1,2) = 0.0D0 !DX/DETA
                GET_JACOBIAN_MATRIX%ENT(2,1) = 0.0D0 !DY/DXI
                GET_JACOBIAN_MATRIX%ENT(2,2) = 1.250D0  !DY/DETA
            ELSEIF (PATCH==2) THEN ! regular mapping G(xi, eta)
                GET_JACOBIAN_MATRIX%ENT(1,1) = 2.0D0  !DX/DXI
                GET_JACOBIAN_MATRIX%ENT(1,2) = 0.0D0 !DX/DETA
                GET_JACOBIAN_MATRIX%ENT(2,1) = 0.0D0 !DY/DXI
                GET_JACOBIAN_MATRIX%ENT(2,2) = 1.250D0 !DY/DETA
            ENDIF
        ELSEIF ((PDE=='ELPT' .AND. PROBLEM==3) .OR. (PDE=='FOUR' .AND. PROBLEM>=4)) THEN 
            IF (PATCH==1) THEN ! singular mapping F(xi, eta)
                !----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!
                GET_JACOBIAN_MATRIX%ENT(1,1) = RADI_SING(2)*2.0D0*PI*PAR_PT%Y**(2)*DSIN(2.0D0*PI*(1.0D0 - PAR_PT%X))     !DX/DXI
                GET_JACOBIAN_MATRIX%ENT(1,2) = RADI_SING(2)*2.0D0*PAR_PT%Y*DCOS(2.0D0*PI*(1.0D0 - PAR_PT%X))       !DX/DETA
                GET_JACOBIAN_MATRIX%ENT(2,1) = RADI_SING(2)*(-2.0D0)*PI*PAR_PT%Y**(2)*DCOS(2.0D0*PI*(1.0D0 - PAR_PT%X))       !DY/DXI
                GET_JACOBIAN_MATRIX%ENT(2,2) = RADI_SING(2)*2.0D0*PAR_PT%Y*DSIN(2.0D0*PI*(1.0D0 - PAR_PT%X))        !DY/DETA
                
!                 GET_JACOBIAN_MATRIX%ENT(1,1) = RADI_SING(2)*2.0D0*PI*PAR_PT%Y*DSIN(2.0D0*PI*(1.0D0 - PAR_PT%X))     !DX/DXI
!                 GET_JACOBIAN_MATRIX%ENT(1,2) = RADI_SING(2)*DCOS(2.0D0*PI*(1.0D0 - PAR_PT%X))       !DX/DETA
!                 GET_JACOBIAN_MATRIX%ENT(2,1) = RADI_SING(2)*(-2.0D0)*PI*PAR_PT%Y*DCOS(2.0D0*PI*(1.0D0 - PAR_PT%X))       !DY/DXI
!                 GET_JACOBIAN_MATRIX%ENT(2,2) = RADI_SING(2)*DSIN(2.0D0*PI*(1.0D0 - PAR_PT%X))        !DY/DETA
                
!                 GET_JACOBIAN_MATRIX%ENT(1,1) = RADI_SING(2)*2.0D0*PI*PAR_PT%Y**(4)*DSIN(2.0D0*PI*(1.0D0 - PAR_PT%X))     !DX/DXI
!                 GET_JACOBIAN_MATRIX%ENT(1,2) = RADI_SING(2)*4.0D0*PAR_PT%Y**(3)*DCOS(2.0D0*PI*(1.0D0 - PAR_PT%X))       !DX/DETA
!                 GET_JACOBIAN_MATRIX%ENT(2,1) = RADI_SING(2)*(-2.0D0)*PI*PAR_PT%Y**(4)*DCOS(2.0D0*PI*(1.0D0 - PAR_PT%X))       !DY/DXI
!                 GET_JACOBIAN_MATRIX%ENT(2,2) = RADI_SING(2)*4.0D0*PAR_PT%Y**(3)*DSIN(2.0D0*PI*(1.0D0 - PAR_PT%X))        !DY/DETA
                !----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!
            ELSEIF (PATCH==2) THEN ! regular mapping G(xi, eta)
                GET_JACOBIAN_MATRIX%ENT(1,1) = 2.0D0*PI*(RADI_SING(1) + (1.0D0 - RADI_SING(1))*PAR_PT%Y)*DSIN(2.0D0*PI*(1.0D0 - PAR_PT%X))     !DX/DXI
                GET_JACOBIAN_MATRIX%ENT(1,2) = (1.0D0 - RADI_SING(1))*DCOS(2.0D0*PI*(1.0D0 - PAR_PT%X))       !DX/DETA
                GET_JACOBIAN_MATRIX%ENT(2,1) = -2.0D0*PI*(RADI_SING(1) + (1.0D0 - RADI_SING(1))*PAR_PT%Y)*DCOS(2.0D0*PI*(1.0D0 - PAR_PT%X))     !DY/DXI
                GET_JACOBIAN_MATRIX%ENT(2,2) = (1.0D0 - RADI_SING(1))*DSIN(2.0D0*PI*(1.0D0 - PAR_PT%X))     !DY/DETA
            ENDIF
        ENDIF
        
	END FUNCTION GET_JACOBIAN_MATRIX
	
	!!  GET SECOND-ORDER PARTIAL DERIVATIVES
	TYPE(SECOND_PARTIAL_DERIVATIVES) FUNCTION GET_PARTIAL(REF_PT, PATCH)

		TYPE(POINT2D), INTENT(IN) :: REF_PT
		INTEGER, INTENT(IN) :: PATCH
        
        IF (PDE=='FOUR' .AND. PROBLEM<=3) THEN
            IF (PATCH==1) THEN ! singular mapping F
                GET_PARTIAL%X(1) = 0.0D0
                GET_PARTIAL%X(2) = 0.0D0
                GET_PARTIAL%X(3) = 0.0D0
                
                GET_PARTIAL%Y(1) = 0.0D0
                GET_PARTIAL%Y(2) = 0.0D0
                GET_PARTIAL%Y(3) = 0.0D0
                
            ELSEIF (PATCH==2) THEN ! regular mapping G
                GET_PARTIAL%X(1) = 0.0D0
                GET_PARTIAL%X(2) = 0.0D0
                GET_PARTIAL%X(3) =0.0D0
                
                GET_PARTIAL%Y(1) = 0.0D0
                GET_PARTIAL%Y(2) = 0.0D0
                GET_PARTIAL%Y(3) = 0.0D0
            ENDIF
        ELSEIF ((PDE=='ELPT' .AND. PROBLEM==3) .OR. (PDE=='FOUR' .AND. PROBLEM>=4)) THEN 
            IF (PATCH==1) THEN ! singular mapping F
                !----------------------------------------------------------------------------------------------------------------------------------------------------!
                GET_PARTIAL%X(1) = -4.0D0*PI**(2)*RADI_SING(2)*REF_PT%Y**(2)*DCOS(2.0D0*PI*(1.0D0 - REF_PT%X))
                GET_PARTIAL%X(2) = 4.0D0*PI*RADI_SING(2)*REF_PT%Y*DSIN(2.0D0*PI*(1.0D0 - REF_PT%X))
                GET_PARTIAL%X(3) = 2.0D0*RADI_SING(2)*DCOS(2.0D0*PI*(1.0D0 - REF_PT%X))
                
                GET_PARTIAL%Y(1) = -4.0D0*PI**(2)*RADI_SING(2)*REF_PT%Y**(2)*DSIN(2.0D0*PI*(1.0D0 - REF_PT%X))
                GET_PARTIAL%Y(2) = -4.0D0*PI*RADI_SING(2)*REF_PT%Y*DCOS(2.0D0*PI*(1.0D0 - REF_PT%X))
                GET_PARTIAL%Y(3) = 2.0D0*RADI_SING(2)*DSIN(2.0D0*PI*(1.0D0 - REF_PT%X))
                
!                 GET_PARTIAL%X(1) = -4.0D0*PI**(2.0d0)*RADI_SING(2)*REF_PT%Y*DCOS(2.0D0*PI*(1.0D0 - REF_PT%X))
!                 GET_PARTIAL%X(2) = 2.0D0*PI*RADI_SING(2)*DSIN(2.0D0*PI*(1.0D0 - REF_PT%X))
!                 GET_PARTIAL%X(3) = 0.0d0
!                 
!                 GET_PARTIAL%Y(1) = -4.0D0*PI**(2)*RADI_SING(2)*REF_PT%Y*DSIN(2.0D0*PI*(1.0D0 - REF_PT%X))
!                 GET_PARTIAL%Y(2) = -2.0D0*PI*RADI_SING(2)*DCOS(2.0D0*PI*(1.0D0 - REF_PT%X))
!                 GET_PARTIAL%Y(3) = 0.0d0
                
!                 GET_PARTIAL%X(1) = -4.0D0*PI**(2)*RADI_SING(2)*REF_PT%Y**(4)*DCOS(2.0D0*PI*(1.0D0 - REF_PT%X))
!                 GET_PARTIAL%X(2) = 8.0D0*PI*RADI_SING(2)*REF_PT%Y**(3)*DSIN(2.0D0*PI*(1.0D0 - REF_PT%X))
!                 GET_PARTIAL%X(3) = 12.0D0*RADI_SING(2)*REF_PT%Y**(2)*DCOS(2.0D0*PI*(1.0D0 - REF_PT%X))
!                 
!                 GET_PARTIAL%Y(1) = -4.0D0*PI**(2)*RADI_SING(2)*REF_PT%Y**(4)*DSIN(2.0D0*PI*(1.0D0 - REF_PT%X))
!                 GET_PARTIAL%Y(2) = -8.0D0*PI*RADI_SING(2)*REF_PT%Y**(3)*DCOS(2.0D0*PI*(1.0D0 - REF_PT%X))
!                 GET_PARTIAL%Y(3) = 12.0D0*RADI_SING(2)*REF_PT%Y**(2)*DSIN(2.0D0*PI*(1.0D0 - REF_PT%X))
                !----------------------------------------------------------------------------------------------------------------------------------------------------!
            ELSEIF (PATCH==2) THEN ! regular mapping G
                GET_PARTIAL%X(1) = -4.0D0*PI**(2)*(RADI_SING(1) + (1.0D0 - RADI_SING(1))*REF_PT%Y)*DCOS(2.0D0*PI*(1.0D0 - REF_PT%X))
                GET_PARTIAL%X(2) = 2.0D0*PI*(1.0D0 - RADI_SING(1))*DSIN(2.0D0*PI*(1.0D0 - REF_PT%X))
                GET_PARTIAL%X(3) =0.0D0
                
                GET_PARTIAL%Y(1) = -4.0D0*PI**(2)*(RADI_SING(1) + (1.0D0 - RADI_SING(1))*REF_PT%Y)*DSIN(2.0D0*PI*(1.0D0 - REF_PT%X))
                GET_PARTIAL%Y(2) = -2.0D0*PI*(1.0D0 - RADI_SING(1))*DCOS(2.0D0*PI*(1.0D0 - REF_PT%X))
                GET_PARTIAL%Y(3) = 0.0D0
            ENDIF
        ENDIF

	END FUNCTION GET_PARTIAL

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
TYPE(FUNCTION_1D) FUNCTION PATCH_TO_REFPU(PHYR, PATCH)

    REAL*8, INTENT(IN) :: PHYR
    INTEGER, INTENT(IN) :: PATCH
		
    REAL*8 :: RPTS(4)
		
		
		IF (PDE=='FOUR' .AND. PROBLEM<=3) THEN 
            RPTS(1) = -1.0D0
            RPTS(2) = RADI_SING(1)
            RPTS(3) = RADI_SING(2)
            RPTS(4) = 1.0D0
        ELSEIF ((PDE=='ELPT' .AND. PROBLEM==3) .OR. (PDE=='FOUR' .AND. PROBLEM>=4)) THEN
            RPTS(1) = 0.0D0
            RPTS(2) = RADI_SING(1)
            RPTS(3) = RADI_SING(2)
            RPTS(4) = 1.0D0
        ENDIF
        
		IF (PATCH==1) THEN      ! singular mapping F
            IF (PHYR>=RPTS(1) .AND. PHYR<=RPTS(2)) THEN
                PATCH_TO_REFPU%VAL(0) = 0.0D0
                PATCH_TO_REFPU%VAL(1) = 0.0D0
            ELSEIF (PHYR>=RPTS(2) .AND. PHYR<=RPTS(3)) THEN
                PATCH_TO_REFPU%VAL(0) = (PHYR-RPTS(2))/(RPTS(3)-RPTS(2))
                PATCH_TO_REFPU%VAL(1) = 1.0D0/(RPTS(3)-RPTS(2))
            ELSEIF (PHYR>=RPTS(3) .AND. PHYR<=RPTS(4)) THEN
                PATCH_TO_REFPU%VAL(0) = 5.0D0
                PATCH_TO_REFPU%VAL(1) = 0.0D0
            ENDIF
        ELSEIF (PATCH==2) THEN      ! regular mapping G
            IF (PHYR>=RPTS(1) .AND. PHYR<=RPTS(2)) THEN
                PATCH_TO_REFPU%VAL(0) = 5.0D0
                PATCH_TO_REFPU%VAL(1) = 0.0D0
            ELSEIF (PHYR>=RPTS(2) .AND. PHYR<=RPTS(3)) THEN
                PATCH_TO_REFPU%VAL(0) = (PHYR-RPTS(3))/(RPTS(3)-RPTS(2))
                PATCH_TO_REFPU%VAL(1) = 1.0D0/(RPTS(3)-RPTS(2))
            ELSEIF (PHYR>=RPTS(3) .AND. PHYR<=RPTS(4)) THEN
                PATCH_TO_REFPU%VAL(0) = 0.0D0
                PATCH_TO_REFPU%VAL(1) = 0.0D0
            ENDIF
        ENDIF
        
END FUNCTION PATCH_TO_REFPU


! !----------------------------------------------------------------------------------------
TYPE(FUNCTION_1D) FUNCTION RADIUS_MAP(PARPT, PATCH)
    TYPE(POINT2D), INTENT(IN) :: PARPT
    INTEGER :: PATCH
    
    IF (PDE=='FOUR' .AND. PROBLEM<=3) THEN 
        IF (PATCH==1) THEN ! singular mapping F(xi, eta)
            RADIUS_MAP%VAL(0) = 1.250D0*PARPT%Y - 1.0D0
            RADIUS_MAP%VAL(1) = 1.250D0
            RADIUS_MAP%VAL(2) = 0.0D0
        ELSEIF (PATCH==2) THEN ! regular mapping G(xi, eta)
            RADIUS_MAP%VAL(0) = 1.250D0*PARPT%Y - 0.250D0
            RADIUS_MAP%VAL(1) = 1.250D0
            RADIUS_MAP%VAL(2) = 0.0D0
        ENDIF
    ELSEIF ((PDE=='ELPT' .AND. PROBLEM==3) .OR. (PDE=='FOUR' .AND. PROBLEM>=4)) THEN 
        IF (PATCH==1) THEN  ! singular mapping F
            !------------------------------------------------------------------------------!
            RADIUS_MAP%VAL(0) = RADI_SING(2)*PARPT%Y**(2)
            RADIUS_MAP%VAL(1) = 2.0D0*RADI_SING(2)*PARPT%Y
            RADIUS_MAP%VAL(2) = 2.0D0*RADI_SING(2)

!             RADIUS_MAP%VAL(0) = RADI_SING(2)*PARPT%Y
!             RADIUS_MAP%VAL(1) = RADI_SING(2)
!             RADIUS_MAP%VAL(2) = 0.0d0
            
!             RADIUS_MAP%VAL(0) = RADI_SING(2)*PARPT%Y**(4)
!             RADIUS_MAP%VAL(1) = 4.0D0*RADI_SING(2)*PARPT%Y**(3)
!             RADIUS_MAP%VAL(2) = 12.0D0*RADI_SING(2)*PARPT%Y**(2)
            !------------------------------------------------------------------------------!
        ELSEIF (PATCH==2) THEN   ! regular mapping G
            RADIUS_MAP%VAL(0) = RADI_SING(1) + (1.0D0 - RADI_SING(1))*PARPT%Y
            RADIUS_MAP%VAL(1) = (1.0D0 - RADI_SING(1))
            RADIUS_MAP%VAL(2) = 0.0D0
        ENDIF
    ENDIF
    
END FUNCTION RADIUS_MAP

!----------------------------------------------------------------------------------------
! TYPE(FUNCTION_1D) FUNCTION UNITSQ_TO_PARSP(REFPT, PATCHES, AXIS)
! 	
! 	REAL*8, INTENT(IN) :: REFPT
! 	INTEGER, INTENT(IN) :: PATCHES(2)
! 	CHARACTER(LEN=1), INTENT(IN) :: AXIS
! 	
! 	REAL*8 :: LEFT_PT, RIGHT_PT
! 	IF (AXIS=='X') THEN
! 		
! 		UNITSQ_TO_PARSP%VAL(0) = REFPT
! 		UNITSQ_TO_PARSP%VAL(1) = 1.0D0
! 	
! 	ELSEIF (AXIS=='Y') THEN
! 		
! 		IF (PATCHES(2)==1) THEN
! 			LEFT_PT = BETA - DELTA
! 			RIGHT_PT = 1.0D0
! 		ELSEIF (PATCHES(2)==2) THEN
! 			LEFT_PT = ALPHA - DELTA
! 			RIGHT_PT = BETA + DELTA
! 		ELSEIF (PATCHES(2)==3) THEN
! 			LEFT_PT = 0.0D0
! 			RIGHT_PT = ALPHA + DELTA
! 		ENDIF
! 		
! 		UNITSQ_TO_PARSP%VAL(0) = (RIGHT_PT - LEFT_PT)*REFPT + LEFT_PT
! 		UNITSQ_TO_PARSP%VAL(1) = (RIGHT_PT - LEFT_PT)
! 		
! 	ENDIF
! 	
! END FUNCTION UNITSQ_TO_PARSP

!----------------------------------------------------------------------------------------


! TYPE(FUNCTION_1D) FUNCTION PARSP_TO_UNITSQ(PHYPT, PATCHES, AXIS)
! 
! 	REAL*8, INTENT(IN) :: PHYPT
! 	INTEGER, INTENT(IN) :: PATCHES(2)
! 	CHARACTER(LEN=1), INTENT(IN) :: AXIS
! 	
! 	REAL*8 :: LEFT_PT, RIGHT_PT
! 	IF (AXIS=='X') THEN
! 		
! 		PARSP_TO_UNITSQ%VAL(0) = PHYPT
! 		PARSP_TO_UNITSQ%VAL(1) = 1.0D0
! 		
! 	ELSEIF (AXIS=='Y') THEN
! 		
! 		IF (PATCHES(2)==1) THEN
! 			LEFT_PT = BETA - DELTA
! 			RIGHT_PT = 1.0D0
! 		ELSEIF (PATCHES(2)==2) THEN
! 			LEFT_PT = ALPHA - DELTA
! 			RIGHT_PT = BETA + DELTA
! 		ELSEIF (PATCHES(2)==3) THEN
! 			LEFT_PT = 0.0D0
! 			RIGHT_PT = ALPHA + DELTA
! 		ENDIF
! 		 
! 		PARSP_TO_UNITSQ%VAL(0) = (PHYPT - LEFT_PT)/(RIGHT_PT - LEFT_PT)
! 		PARSP_TO_UNITSQ%VAL(1) = 1.0D0/(RIGHT_PT - LEFT_PT)
! 		
! 	ENDIF
! 	
! END FUNCTION PARSP_TO_UNITSQ
! 
! TYPE(FUNCTION_1D) FUNCTION PHYPATCH1_TO_UNITSQ(PHYPT, AXIS)
! 
! 	REAL*8, INTENT(IN) :: PHYPT
! 	CHARACTER(LEN=1), INTENT(IN) :: AXIS
! 	
! 	IF (AXIS=='X') THEN
! 		PHYPATCH1_TO_UNITSQ%VAL(0) = (1.0D0/(2.0D0*(1.0D0 - BETA + DELTA)))*(PHYPT + (1.0D0 - BETA + DELTA))
! 		PHYPATCH1_TO_UNITSQ%VAL(1) = (1.0D0/(2.0D0*(1.0D0 - BETA + DELTA)))
! 	ELSEIF (AXIS=='Y') THEN
! 		PHYPATCH1_TO_UNITSQ%VAL(0) = (1.0D0/(2.0D0*(1.0D0 - BETA + DELTA)))*(PHYPT + (1.0D0 - BETA + DELTA))
! 		PHYPATCH1_TO_UNITSQ%VAL(1) = (1.0D0/(2.0D0*(1.0D0 - BETA + DELTA)))
! 	ENDIF
! 
! END FUNCTION PHYPATCH1_TO_UNITSQ

END MODULE PATCH_MAPPING
