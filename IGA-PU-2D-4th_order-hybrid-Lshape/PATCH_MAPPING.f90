MODULE PATCH_MAPPING

	USE NURBS

	IMPLICIT NONE

CONTAINS

	!!  GET THE CORRESPONDING POINT ON THE PHYSICAL SPACE OF GIVE POINT ON THE PARAMETRIC SPACE
	TYPE(POINT2D) FUNCTION GET_PHY_PT(PAR_PT, PATCH)

		TYPE(POINT2D), INTENT(IN) :: PAR_PT
		INTEGER, INTENT(IN) :: PATCH
		
		REAL*8 :: FACTOR, ANGLE
		TYPE(FUNCTION_1D) :: SY(3), CL(2), C(2, 2)
		
        IF (PATCH==1) THEN 
            !--------------------------------------------!
            FACTOR = RADI_SING(2)*PAR_PT%Y**(RLMBD)
            ANGLE = ANGL_FTR*PI*(1.0D0 - PAR_PT%X)
            !--------------------------------------------!
            GET_PHY_PT%X = FACTOR*DCOS(ANGLE)
            GET_PHY_PT%Y = FACTOR*DSIN(ANGLE)
        ELSEIF (PATCH==2) THEN 
        !-----------------------------------------------------------!
            CALL SLD_MODF_BEZIER_Y(SY, PAR_PT%Y)
            CALL C1BSCURVE(CL, PAR_PT%X)
            CALL CINFCURVE(C, PAR_PT%X)
            GET_PHY_PT%X = SY(1)%VAL(0)*C(1, 1)%VAL(0) + SY(2)%VAL(0)*C(2, 1)%VAL(0) + SY(3)%VAL(0)*CL(1)%VAL(0)
            GET_PHY_PT%Y = SY(1)%VAL(0)*C(1, 2)%VAL(0) + SY(2)%VAL(0)*C(2, 2)%VAL(0) + SY(3)%VAL(0)*CL(2)%VAL(0)
        !-----------------------------------------------------------!
!             GET_PHY_PT = GET_POINT_NURVE_SURFACE_2D(PAR_PT, GEO_KVEC, GEO_CTL)
        !-----------------------------------------------------------!
        ENDIF
		
	END FUNCTION GET_PHY_PT

    TYPE(POINT2D) FUNCTION GET_PAR_PT(PHY_PT, PATCH)
        
        TYPE(POINT2D), INTENT(IN) :: PHY_PT
        INTEGER, INTENT(IN) :: PATCH
        REAL*8 :: R, ARCOS
        
        R = DSQRT(PHY_PT%X**2 + PHY_PT%Y**2)
        if (dabs(phy_pt%x)<=eps) then 
            arcos = 0.50d0*pi
        else 
            ARCOS = DACOS(PHY_PT%X/R)
        endif
        
        IF (PATCH==1) THEN      ! singular mapping F
            if (dabs(phy_pt%y)<=eps) then
                get_par_pt%x = 1.0d0
                goto 372
            elseIF (PHY_PT%Y<EPS) THEN 
                GET_PAR_PT%X = (2.0D0/(3.0D0*PI))*ARCOS - (1.0D0/3.0D0)
            ELSE
                GET_PAR_PT%X = 1.0D0 - (2.0D0/(3.0D0*PI))*ARCOS
            ENDIF
            372 continue
            !-----------------------------------------------------------!
            GET_PAR_PT%Y = (R/RADI_SING(2))**(1.0D0/RLMBD)
            !-----------------------------------------------------------!
        ELSEIF (PATCH==2) THEN      ! reguar mapping G
            GET_PAR_PT = NEWTON_REFPT(PHY_PT, PATCH)
            print*, phy_pt, patch
        ENDIF
        
    END FUNCTION GET_PAR_PT

    !!  GET THE CORRESPONDING POINT ON THE PARAMETRIC SPACE OF GIVE POINT ON THE PHYSICAL SPACE
	TYPE(POINT2D) FUNCTION NEWTON_REFPT(PHY_PT, PATCH, SIDE)
        
		TYPE(POINT2D), INTENT(IN) :: PHY_PT
		INTEGER, INTENT(IN) :: PATCH
		CHARACTER(LEN=1), OPTIONAL, INTENT(IN) :: SIDE
        
        TYPE(MATRIX_22) :: JACOBIAN
		TYPE(POINT2D) :: NEW_PT, OLD_PT, TMP_PT
		REAL(8) :: ERROR, INV_J
		INTEGER :: ITER
        
		ERROR = 1.0D0
		
		OLD_PT = GOOD_INITIAL(PHY_PT, PATCH)
!         print*, 'OLD_PT', OLD_PT
		ITER = 0
		DO WHILE (ERROR>TOLERANCE)
			ITER = ITER + 1
			TMP_PT =  GET_PHY_PT(OLD_PT, PATCH)
            
			TMP_PT%X = TMP_PT%X - PHY_PT%X
			TMP_PT%Y = TMP_PT%Y - PHY_PT%Y

			ERROR = MAX(DABS(TMP_PT%X), DABS(TMP_PT%Y))
			
			IF (ERROR<TOLERANCE .AND. OLD_PT%X.GE.0.D0 .AND. OLD_PT%X.LE.1.D0 .AND. OLD_PT%Y.GE.0.D0 .AND. OLD_PT%Y.LE.1.D0) THEN
				NEWTON_REFPT = OLD_PT
				GOTO 999
			ELSE
                JACOBIAN = GET_JACOBIAN_MATRIX(OLD_PT, PATCH)
				INV_J = 1.0D0 / (JACOBIAN%ENT(1,1)*JACOBIAN%ENT(2,2) - JACOBIAN%ENT(1,2)*JACOBIAN%ENT(2,1))
				NEW_PT%X = OLD_PT%X - INV_J*(JACOBIAN%ENT(2,2)*TMP_PT%X - JACOBIAN%ENT(1,2)*TMP_PT%Y)
				NEW_PT%Y = OLD_PT%Y - INV_J*(JACOBIAN%ENT(1,1)*TMP_PT%Y - JACOBIAN%ENT(2,1)*TMP_PT%X)
				OLD_PT = NEW_PT
				
! 				PRINT*, PHY_PT, OLD_PT
! 				PRINT*, 'GET_REF_PT-HERE2'
			ENDIF
		ENDDO
		
		IF (ITER==MAX_ITERATION) THEN
			PRINT *, '[ERROR]  CANNOT FIND REFERENCE POINT !'
		ENDIF

		999 CONTINUE
		
! 		PRINT*, 'GET_REF_PT-HERE3'
		
	END FUNCTION NEWTON_REFPT
	
	
	!! HOMOTOPY OR CONTINUATION ALGORITHM TO FIND GOOD INITIAL APPROXIMATION FOR NEWTON'S OR BROYDEN'S METHOD
	!! REFERENCE :: NUMERICAL ANALYSIS (BURDEN & FAIR)
	TYPE(POINT2D) FUNCTION GOOD_INITIAL(PHY_PT, PATCH)
		
		TYPE(POINT2D), INTENT(IN) :: PHY_PT
		INTEGER, INTENT(IN) :: PATCH
		
		TYPE(FUNCTION_2D) :: DIFF_NURBS(2)
		TYPE(MATRIX_22) :: JACOB, INV_JACOB
		TYPE(POINT2D) :: K(4), OLD_PT, INI_B
		INTEGER :: I, J, II, JJ, MAX_IT
		REAL*8 :: STEP_H
		
		MAX_IT = 20
		STEP_H = 1.D0/(1.D0*MAX_IT)
		OLD_PT = POINT2D(0.50D0,0.50D0)
		
! 		CALL GET_DIFF_NURVE_SURFACE_2D(DIFF_NURBS,OLD_PT,GEO_KVEC(PATCH,:),GEO_CTL(PATCH),1)

		INI_B = GET_PHY_PT(OLD_PT, PATCH)
		INI_B%X = -STEP_H*(INI_B%X - PHY_PT%X)
		INI_B%Y = -STEP_H*(INI_B%Y - PHY_PT%Y)
		
		DO I = 1, MAX_IT
			JACOB = GET_JACOBIAN_MATRIX(OLD_PT, PATCH)
			INV_JACOB = .INVERSE.JACOB
			K(1)%X = DOT_PRODUCT((/INV_JACOB%ENT(1,1), INV_JACOB%ENT(1,2)/), (/INI_B%X, INI_B%Y/))
			K(1)%Y = DOT_PRODUCT((/INV_JACOB%ENT(2,1), INV_JACOB%ENT(2,2)/), (/INI_B%X, INI_B%Y/))
			
			JACOB = GET_JACOBIAN_MATRIX(POINT2D(OLD_PT%X + 0.5D0*K(1)%X, OLD_PT%Y + 0.5D0*K(1)%Y), PATCH)
			INV_JACOB = .INVERSE.JACOB
			K(2)%X = DOT_PRODUCT((/INV_JACOB%ENT(1,1), INV_JACOB%ENT(1,2)/), (/INI_B%X, INI_B%Y/))
			K(2)%Y = DOT_PRODUCT((/INV_JACOB%ENT(2,1), INV_JACOB%ENT(2,2)/), (/INI_B%X, INI_B%Y/))
			
			JACOB = GET_JACOBIAN_MATRIX(POINT2D(OLD_PT%X + 0.5D0*K(2)%X, OLD_PT%Y + 0.5D0*K(2)%Y), PATCH)
			INV_JACOB = .INVERSE.JACOB
			K(3)%X = DOT_PRODUCT((/INV_JACOB%ENT(1,1), INV_JACOB%ENT(1,2)/), (/INI_B%X, INI_B%Y/))
			K(3)%Y = DOT_PRODUCT((/INV_JACOB%ENT(2,1), INV_JACOB%ENT(2,2)/), (/INI_B%X, INI_B%Y/))
			
			JACOB = GET_JACOBIAN_MATRIX(POINT2D(OLD_PT%X + 0.5D0*K(3)%X, OLD_PT%Y + 0.5D0*K(3)%Y), PATCH)
			INV_JACOB = .INVERSE.JACOB
			K(4)%X = DOT_PRODUCT((/INV_JACOB%ENT(1,1), INV_JACOB%ENT(1,2)/), (/INI_B%X, INI_B%Y/))
			K(4)%Y = DOT_PRODUCT((/INV_JACOB%ENT(2,1), INV_JACOB%ENT(2,2)/), (/INI_B%X, INI_B%Y/))
			
			OLD_PT%X = OLD_PT%X + (1.D0/6.D0)*(K(1)%X + 2.D0*K(2)%X + 2.D0*K(3)%X + K(4)%X)
			OLD_PT%Y = OLD_PT%Y + (1.D0/6.D0)*(K(1)%Y + 2.D0*K(2)%Y + 2.D0*K(3)%Y + K(4)%Y)
		ENDDO
		
		GOOD_INITIAL = OLD_PT
	END FUNCTION GOOD_INITIAL
	
	
	!!  GET JACOBIAN OF THE GEOMETRIC MAPPING
	TYPE(MATRIX_22) FUNCTION GET_JACOBIAN_MATRIX(PAR_PT, PATCH)

		TYPE(POINT2D), INTENT(IN) :: PAR_PT
        INTEGER, INTENT(IN) :: PATCH
        
		TYPE(FUNCTION_2D) :: DIFF_NURBS(2)
		TYPE(FUNCTION_1D) :: SY(3), CL(2), C(2, 2)
		REAL*8 :: ANGLE
		
!		TYPE(MATRIX_22) :: TEST_MATRIX1, TEST_MATRIX2, test_result
        
        IF (PATCH==1) THEN ! singular mapping F(xi, eta)
            ANGLE = ANGL_FTR*PI*(1.0D0 - PAR_PT%X)
            !----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!
            GET_JACOBIAN_MATRIX%ENT(1,1) = RADI_SING(2)*(ANGL_FTR)*PI*PAR_PT%Y**(RLMBD)*DSIN(ANGLE)     !DX/DXI
            GET_JACOBIAN_MATRIX%ENT(1,2) = RADI_SING(2)*(RLMBD)*PAR_PT%Y**(RLMBD - 1.0D0)*DCOS(ANGLE)       !DX/DETA
            GET_JACOBIAN_MATRIX%ENT(2,1) = RADI_SING(2)*(-ANGL_FTR)*PI*PAR_PT%Y**(RLMBD)*DCOS(ANGLE)       !DY/DXI
            GET_JACOBIAN_MATRIX%ENT(2,2) = RADI_SING(2)*(RLMBD)*PAR_PT%Y**(RLMBD - 1.0D0)*DSIN(ANGLE)        !DY/DETA
            !----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!
        ELSEIF (PATCH==2) THEN ! regular mapping G(xi, eta)
!-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!        
            CALL SLD_MODF_BEZIER_Y(SY, PAR_PT%Y)
            CALL C1BSCURVE(CL, PAR_PT%X)
            CALL CINFCURVE(C, PAR_PT%X)
            GET_JACOBIAN_MATRIX%ENT(1,1) = SY(1)%VAL(0)*C(1, 1)%VAL(1) + SY(2)%VAL(0)*C(2, 1)%VAL(1) + SY(3)%VAL(0)*CL(1)%VAL(1)     !DX/DXI
            GET_JACOBIAN_MATRIX%ENT(1,2) = SY(1)%VAL(1)*C(1, 1)%VAL(0) + SY(2)%VAL(1)*C(2, 1)%VAL(0) + SY(3)%VAL(1)*CL(1)%VAL(0)       !DX/DETA
            GET_JACOBIAN_MATRIX%ENT(2,1) = SY(1)%VAL(0)*C(1, 2)%VAL(1) + SY(2)%VAL(0)*C(2, 2)%VAL(1) + SY(3)%VAL(0)*CL(2)%VAL(1)     !DY/DXI
            GET_JACOBIAN_MATRIX%ENT(2,2) = SY(1)%VAL(1)*C(1, 2)%VAL(0) + SY(2)%VAL(1)*C(2, 2)%VAL(0) + SY(3)%VAL(1)*CL(2)%VAL(0)     !DY/DETA
!-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!
!             CALL GET_DIFF_NURVE_SURFACE_2D(DIFF_NURBS, PAR_PT, GEO_KVEC, GEO_CTL, 1)
!             GET_JACOBIAN_MATRIX%ENT(1,1) = DIFF_NURBS(1)%VAL(1,0)
!             GET_JACOBIAN_MATRIX%ENT(1,2) = DIFF_NURBS(1)%VAL(0,1)
!             GET_JACOBIAN_MATRIX%ENT(2,1) = DIFF_NURBS(2)%VAL(1,0)
!             GET_JACOBIAN_MATRIX%ENT(2,2) = DIFF_NURBS(2)%VAL(0,1)
!-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!            
        ENDIF
        
	END FUNCTION GET_JACOBIAN_MATRIX
	
	!!  GET SECOND-ORDER PARTIAL DERIVATIVES
	TYPE(SECOND_PARTIAL_DERIVATIVES) FUNCTION GET_PARTIAL(REF_PT, PATCH)
        
		TYPE(POINT2D), INTENT(IN) :: REF_PT
		INTEGER, INTENT(IN) :: PATCH
		
		TYPE(FUNCTION_1D) :: SY(3), CL(2), C(2, 2)
		TYPE(FUNCTION_2D) :: DIFF_NURBS(2)
		REAL*8 :: ANGLE
        
        IF (PATCH==1) THEN ! singular mapping F
            ANGLE = ANGL_FTR*PI*(1.0D0 - REF_PT%X)
            !----------------------------------------------------------------------------------------------------------------------------------------------------!
            GET_PARTIAL%X(1) = (-ANGL_FTR**2)*PI**(2)*RADI_SING(2)*REF_PT%Y**(RLMBD)*DCOS(ANGLE)
            GET_PARTIAL%X(2) = (ANGL_FTR)*PI*RADI_SING(2)*(RLMBD)*REF_PT%Y**(RLMBD - 1.0D0)*DSIN(ANGLE)
            GET_PARTIAL%X(3) = (RLMBD)*(RLMBD - 1.0D0)*RADI_SING(2)*REF_PT%Y**(RLMBD - 2.0D0)*DCOS(ANGLE)
            
            GET_PARTIAL%Y(1) = (-ANGL_FTR**2)*PI**(2)*RADI_SING(2)*REF_PT%Y**(RLMBD)*DSIN(ANGLE)
            GET_PARTIAL%Y(2) = (-ANGL_FTR)*PI*RADI_SING(2)*(RLMBD)*REF_PT%Y**(RLMBD - 1.0D0)*DCOS(ANGLE)
            GET_PARTIAL%Y(3) = (RLMBD)*(RLMBD - 1.0D0)*RADI_SING(2)*REF_PT%Y**(RLMBD - 2.0D0)*DSIN(ANGLE)
            !----------------------------------------------------------------------------------------------------------------------------------------------------!
        ELSEIF (PATCH==2) THEN ! regular mapping G
!-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!        
            CALL SLD_MODF_BEZIER_Y(SY, REF_PT%Y)
            CALL C1BSCURVE(CL, REF_PT%X)
            CALL CINFCURVE(C, REF_PT%X)
            GET_PARTIAL%X(1) = SY(1)%VAL(0)*C(1, 1)%VAL(2) + SY(2)%VAL(0)*C(2, 1)%VAL(2) + SY(3)%VAL(0)*CL(1)%VAL(2)
            GET_PARTIAL%X(2) = SY(1)%VAL(1)*C(1, 1)%VAL(1) + SY(2)%VAL(1)*C(2, 1)%VAL(1) + SY(3)%VAL(1)*CL(1)%VAL(1)
            GET_PARTIAL%X(3) = SY(1)%VAL(2)*C(1, 1)%VAL(0) + SY(2)%VAL(2)*C(2, 1)%VAL(0) + SY(3)%VAL(2)*CL(1)%VAL(0)
            
            GET_PARTIAL%Y(1) = SY(1)%VAL(0)*C(1, 2)%VAL(2) + SY(2)%VAL(0)*C(2, 2)%VAL(2) + SY(3)%VAL(0)*CL(2)%VAL(2)
            GET_PARTIAL%Y(2) = SY(1)%VAL(1)*C(1, 2)%VAL(1) + SY(2)%VAL(1)*C(2, 2)%VAL(1) + SY(3)%VAL(1)*CL(2)%VAL(1)
            GET_PARTIAL%Y(3) = SY(1)%VAL(2)*C(1, 2)%VAL(0) + SY(2)%VAL(2)*C(2, 2)%VAL(0) + SY(3)%VAL(2)*CL(2)%VAL(0)
!-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!
!             CALL GET_DIFF_NURVE_SURFACE_2D(DIFF_NURBS,REF_PT,GEO_KVEC,GEO_CTL,2)
! 
!             GET_PARTIAL%X(1) = DIFF_NURBS(1)%VAL(2,0)
!             GET_PARTIAL%X(2) = DIFF_NURBS(1)%VAL(1,1)
!             GET_PARTIAL%X(3) = DIFF_NURBS(1)%VAL(0,2)
! 
!             GET_PARTIAL%Y(1) = DIFF_NURBS(2)%VAL(2,0)
!             GET_PARTIAL%Y(2) = DIFF_NURBS(2)%VAL(1,1)
!             GET_PARTIAL%Y(3) = DIFF_NURBS(2)%VAL(0,2)
!-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!            
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


TYPE(FUNCTION_1D) FUNCTION PATCH_TO_REFPU(PHYR, PATCH)

    REAL*8, INTENT(IN) :: PHYR
    INTEGER, INTENT(IN) :: PATCH
		
    REAL*8 :: RPTS(4)
		
    RPTS(1) = 0.0D0
    RPTS(2) = RADI_SING(1)
    RPTS(3) = RADI_SING(2)
    RPTS(4) = 1.0D0
    
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
    
    PATCH_TO_REFPU%VAL(2) = 0.0D0
        
END FUNCTION PATCH_TO_REFPU


! !----------------------------------------------------------------------------------------
TYPE(FUNCTION_1D) FUNCTION RADIUS_MAP(PARPT, PATCH)
    TYPE(POINT2D), INTENT(IN) :: PARPT
    INTEGER :: PATCH
    
    IF (PATCH==1) THEN  ! singular mapping F
        !------------------------------------------------------------------------------!
        RADIUS_MAP%VAL(0) = RADI_SING(2)*PARPT%Y**(RLMBD)
        RADIUS_MAP%VAL(1) = (RLMBD)*RADI_SING(2)*PARPT%Y**(RLMBD - 1.0D0)
        RADIUS_MAP%VAL(2) = (RLMBD)*(RLMBD - 1.0D0)*RADI_SING(2)*PARPT%Y**(RLMBD - 2.0D0)
        !------------------------------------------------------------------------------!
    ELSEIF (PATCH==2) THEN   ! regular mapping G
!         if (parpt%y<=inv_delta(2)) then 
            RADIUS_MAP%VAL(0) = (3.0D0*DELTA**3 + 6.0D0*DELTA*PARPT%Y**2 - 4.0D0*PARPT%Y**3)/(10.0D0*DELTA**3)
            RADIUS_MAP%VAL(1) = (6.0D0*(DELTA - PARPT%Y)*PARPT%Y)/(5.0D0*DELTA**3)
            RADIUS_MAP%VAL(2) = (6.0D0*(DELTA - 2.0D0*PARPT%Y))/(5.0D0*DELTA**3)
!         else    
!             RADIUS_MAP%val(0) = 0.80d0
!             RADIUS_MAP%val(1) = 0.0d0
!             RADIUS_MAP%val(2) = 0.0d0
!         endif
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


TYPE(FUNCTION_1D) FUNCTION PARSP_TO_UNITSQ(REFX, INTERVAL, AXIS)

	REAL*8, INTENT(IN) :: REFX
	INTEGER, INTENT(IN) :: INTERVAL
	CHARACTER(LEN=1), INTENT(IN) :: AXIS
	
	REAL*8 :: LEFT_PT, RIGHT_PT
    INTEGER :: I, J, K
    
    IF (AXIS=='X') THEN 
        PARSP_TO_UNITSQ%VAL(0) = 6.0D0*(REFX - 1.0D0*(INTERVAL-1)/6.0D0)
        PARSP_TO_UNITSQ%VAL(1) = 6.0D0
        PARSP_TO_UNITSQ%VAL(2) = 0.0D0
    ELSEIF (AXIS=='Y') THEN 
        IF (INTERVAL==1) THEN 
            PARSP_TO_UNITSQ%VAL(0) = REFX/DELTA
            PARSP_TO_UNITSQ%VAL(1) = 1.0D0/DELTA
            PARSP_TO_UNITSQ%VAL(2) = 0.0D0
        ELSEIF (INTERVAL==2) THEN
            PARSP_TO_UNITSQ%VAL(0) = (1.0D0/(1.0D0 - DELTA))*(REFX - DELTA)
            PARSP_TO_UNITSQ%VAL(1) = (1.0D0/(1.0D0 - DELTA))
            PARSP_TO_UNITSQ%VAL(2) = 0.0D0
        ENDIF
    ENDIF

END FUNCTION PARSP_TO_UNITSQ

TYPE(FUNCTION_1D) FUNCTION MODF_BEZIER(REFX, INDX) ! Modified C1 Bezier functions
    
    REAL*8, INTENT(IN) :: REFX
    INTEGER, INTENT(IN) :: INDX
    
    IF (INDX==1) THEN 
        MODF_BEZIER%VAL(0) = (1.0D0 - REFX)**2*(1.0D0 + 2.0D0*REFX)
        MODF_BEZIER%VAL(1) = 6.0D0*(-1.0D0 + REFX)*REFX
        MODF_BEZIER%VAL(2) = -6.0D0 + 12.0D0*REFX
    ELSEIF (INDX==2) THEN
        MODF_BEZIER%VAL(0) = REFX**2*(3.0D0 - 2.0D0*REFX)
        MODF_BEZIER%VAL(1) = -6.0D0*(-1.0D0 + REFX)*REFX
        MODF_BEZIER%VAL(2) = 6.0D0 - 12.0D0*REFX
    ELSEIF (INDX==3) THEN
        MODF_BEZIER%VAL(0) = 3.0D0*(1.0D0 - REFX)**2*REFX**2
        MODF_BEZIER%VAL(1) = 6.0D0*REFX*(1.0D0 - 3.0D0*REFX + 2.0D0*REFX**2)
        MODF_BEZIER%VAL(2) = 6.0D0*(1.0D0 - 6.0D0*REFX + 6.0D0*REFX**2)
    ENDIF
    
END FUNCTION MODF_BEZIER

SUBROUTINE SLD_MODF_BEZIER_Y(S, REFY)
    
    REAL*8, INTENT(IN) :: REFY
    TYPE(FUNCTION_1D), INTENT(OUT) :: S(3)
    
    REAL*8 :: ETA(0:2)
    INTEGER :: I, J, K
    TYPE(FUNCTION_1D) :: N(2), T
    
    ETA(0) = 0.0D0
    ETA(1) = DELTA
    ETA(2) = 1.0D0
    
    ! Define S_1
    IF (REFY>=ETA(0) .AND. REFY<=ETA(1)) THEN 
        T = PARSP_TO_UNITSQ(REFY, 1, 'Y')
        N(1) = MODF_BEZIER(T%VAL(0), 1)
        S(1)%VAL(0) = N(1)%VAL(0)
        S(1)%VAL(1) = N(1)%VAL(1)*T%VAL(1)
        S(1)%VAL(2) = N(1)%VAL(2)*T%VAL(1)**2 + N(1)%VAL(1)*T%VAL(2)
    ELSE 
        S(1)%VAL(:) = 0.0D0
    ENDIF
    
    ! Define S_2
    DO I = 2, 2
        IF (REFY>=ETA(I-2) .AND. REFY<ETA(I-1)) THEN 
            T = PARSP_TO_UNITSQ(REFY, I-1, 'Y')
            N(2) = MODF_BEZIER(T%VAL(0), 2)
            S(I)%VAL(0) = N(2)%VAL(0)
            S(I)%VAL(1) = N(2)%VAL(1)*T%VAL(1)
            S(I)%VAL(2) = N(2)%VAL(2)*T%VAL(1)**2 + N(2)%VAL(1)*T%VAL(2)
        ELSEIF (REFY>=ETA(I-1) .AND. REFY<=ETA(I)) THEN 
            T = PARSP_TO_UNITSQ(REFY, I, 'Y')
            N(1) = MODF_BEZIER(T%VAL(0), 1)
            S(I)%VAL(0) = N(1)%VAL(0)
            S(I)%VAL(1) = N(1)%VAL(1)*T%VAL(1)
            S(I)%VAL(2) = N(1)%VAL(2)*T%VAL(1)**2 + N(1)%VAL(1)*T%VAL(2)
        ELSE 
            S(I)%VAL(:) = 0.0D0
        ENDIF
    ENDDO
    
    ! Define S_3
    IF (REFY>=ETA(1) .AND. REFY<=ETA(2)) THEN 
        T = PARSP_TO_UNITSQ(REFY, 2, 'Y')
        N(2) = MODF_BEZIER(T%VAL(0), 2)
        S(3)%VAL(0) = N(2)%VAL(0)
        S(3)%VAL(1) = N(2)%VAL(1)*T%VAL(1)
        S(3)%VAL(2) = N(2)%VAL(2)*T%VAL(1)**2 + N(2)%VAL(1)*T%VAL(2)
    ELSE 
        S(3)%VAL(:) = 0.0D0
    ENDIF
    
END SUBROUTINE SLD_MODF_BEZIER_Y

SUBROUTINE SLD_MODF_BEZIER_X(S, REFX)
    
    REAL*8, INTENT(IN) :: REFX
    TYPE(FUNCTION_1D), INTENT(OUT) :: S(7)
    
    REAL*8 :: XI(0:6)
    INTEGER :: I, J, K
    TYPE(FUNCTION_1D) :: N(2), T
    
    DO I = 0, 6
        XI(I) = 1.0D0*I/6.0D0
    ENDDO
    
    ! Define S_1
    IF (REFX>=XI(0) .AND. REFX<=XI(1)) THEN 
        T = PARSP_TO_UNITSQ(REFX, 1, 'X')
        N(1) = MODF_BEZIER(T%VAL(0), 1)
        S(1)%VAL(0) = N(1)%VAL(0)
        S(1)%VAL(1) = N(1)%VAL(1)*T%VAL(1)
        S(1)%VAL(2) = N(1)%VAL(2)*T%VAL(1)**2 + N(1)%VAL(1)*T%VAL(2)
    ELSE 
        S(1)%VAL(:) = 0.0D0
    ENDIF
    
    ! Define S_2 ~ S_6
    DO I = 2, 6
        IF (REFX>=XI(I-2) .AND. REFX<XI(I-1)) THEN 
            T = PARSP_TO_UNITSQ(REFX, I-1, 'X')
            N(2) = MODF_BEZIER(T%VAL(0), 2)
            S(I)%VAL(0) = N(2)%VAL(0)
            S(I)%VAL(1) = N(2)%VAL(1)*T%VAL(1)
            S(I)%VAL(2) = N(2)%VAL(2)*T%VAL(1)**2 + N(2)%VAL(1)*T%VAL(2)
        ELSEIF (REFX>=XI(I-1) .AND. REFX<=XI(I)) THEN 
            T = PARSP_TO_UNITSQ(REFX, I, 'X')
            N(1) = MODF_BEZIER(T%VAL(0), 1)
            S(I)%VAL(0) = N(1)%VAL(0)
            S(I)%VAL(1) = N(1)%VAL(1)*T%VAL(1)
            S(I)%VAL(2) = N(1)%VAL(2)*T%VAL(1)**2 + N(1)%VAL(1)*T%VAL(2)
        ELSE 
            S(I)%VAL(:) = 0.0D0
        ENDIF
    ENDDO
    
    ! Define S_7
    IF (REFX>=XI(5) .AND. REFX<=XI(6)) THEN 
        T = PARSP_TO_UNITSQ(REFX, 6, 'X')
        N(2) = MODF_BEZIER(T%VAL(0), 2)
        S(7)%VAL(0) = N(2)%VAL(0)
        S(7)%VAL(1) = N(2)%VAL(1)*T%VAL(1)
        S(7)%VAL(2) = N(2)%VAL(2)*T%VAL(1)**2 + N(2)%VAL(1)*T%VAL(2)
    ELSE 
        S(7)%VAL(:) = 0.0D0
    ENDIF
    
END SUBROUTINE SLD_MODF_BEZIER_X

SUBROUTINE C1BSCURVE(CL, REFX)
    
    REAL*8, INTENT(IN) :: REFX
    TYPE(FUNCTION_1D), INTENT(OUT) :: CL(2)
    
    INTEGER :: I, J, K
    TYPE(POINT2D) :: CTLPT(7)
    TYPE(FUNCTION_1D) :: SX(7)
    
    CTLPT = (/ POINT2D(0.0D0, -1.0D0), POINT2D(-1.0D0, -1.0D0), POINT2D(-1.0D0, 0.0D0), POINT2D(-1.0D0, 1.0D0), POINT2D(0.0D0, 1.0D0), POINT2D(1.0D0, 1.0D0), POINT2D(1.0D0, 0.0D0) /)
    
    CL(:)%VAL(0) = 0.0D0
    CL(:)%VAL(1) = 0.0D0
    CL(:)%VAL(2) = 0.0D0
    
    CALL SLD_MODF_BEZIER_X(SX, REFX)
    
    DO K = 1, 7
        DO J = 0, 2
            CL(1)%VAL(J) = CL(1)%VAL(J) + CTLPT(K)%X*SX(K)%VAL(J)
            CL(2)%VAL(J) = CL(2)%VAL(J) + CTLPT(K)%Y*SX(K)%VAL(J)
        ENDDO
    ENDDO
    
END SUBROUTINE C1BSCURVE


SUBROUTINE CINFCURVE(C, REFX)

    REAL*8, INTENT(IN) :: REFX
    TYPE(FUNCTION_1D), INTENT(OUT) :: C(2, 2)     ! C(1, 1) = C_1(XI), C(1, 2) = C_1(ETA), C(2, 1) = C_2(XI), C(2, 2) = C_2(ETA)
    
    REAL*8 :: ANGLE
    
    ANGLE = (ANGL_FTR)*PI*(1.0D0 - REFX)
    
    C(1, 1)%VAL(0) = RADI_SING(1)*DCOS(ANGLE)
    C(1, 1)%VAL(1) = (ANGL_FTR)*RADI_SING(1)*PI*DSIN(ANGLE)
    C(1, 1)%VAL(2) = -(ANGL_FTR**2)*RADI_SING(1)*PI**2*DCOS(ANGLE)
    
    C(1, 2)%VAL(0) = RADI_SING(1)*DSIN(ANGLE)
    C(1, 2)%VAL(1) = -ANGL_FTR*RADI_SING(1)*PI*DCOS(ANGLE)
    C(1, 2)%VAL(2) = -(ANGL_FTR**2)*RADI_SING(1)*PI**2*DSIN(ANGLE)
    
    C(2, 1)%VAL(0) = RADI_SING(2)*DCOS(ANGLE)
    C(2, 1)%VAL(1) = ANGL_FTR*RADI_SING(2)*PI*DSIN(ANGLE)
    C(2, 1)%VAL(2) = -(ANGL_FTR**2)*RADI_SING(2)*PI**2*DCOS(ANGLE)
    
    C(2, 2)%VAL(0) = RADI_SING(2)*DSIN(ANGLE)
    C(2, 2)%VAL(1) = -ANGL_FTR*RADI_SING(2)*PI*DCOS(ANGLE)
    C(2, 2)%VAL(2) = -(ANGL_FTR**2)*RADI_SING(2)*PI**2*DSIN(ANGLE)
    
END SUBROUTINE CINFCURVE

            
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
