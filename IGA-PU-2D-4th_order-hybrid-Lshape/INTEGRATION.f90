	MODULE INTEGRATION

	USE GSQUAD
	USE NURBS_BASIS
	USE LOADFUNCTION

    IMPLICIT NONE

CONTAINS


!------------------INTEGRAL CODE FOR THE STIFFNESS MATRIX ELEMENT--------------
SUBROUTINE GEN_KF(STIF_K, LOAD_F)

	REAL*8, INTENT(OUT) :: STIF_K(DOF, DOF), LOAD_F(DOF)
	TYPE(RECPATCH) :: IRBOX
	TYPE(POINT2D) :: GSPT, PHYPT, PARPT
	TYPE(INT2D), ALLOCATABLE :: INDX_ROW(:), INDX_COLUMN(:)
	TYPE(MATRIX_22) :: JACOB(2)
	TYPE(FVALUE), ALLOCATABLE :: SF_ROW(:,:), SF_COLUMN(:,:)
	REAL*8, ALLOCATABLE :: WEIGHT(:), LOAD_Q_VEC(:)
	INTEGER :: I,J, II, JJ, KK, GLOBAL_INDX(2), PATCH_ROW, PATCH_COLUMN, PATCHES(2), N, PATCH, KSPAN_NDX(2)
	REAL*8 :: DIFF_XX, DIFF_YY, DIFF_XY, DIFF_YX, DIFF_XN, DIFF_YN, DIFF_NX, DIFF_NY, DIFF_NN, DET_M, DIFF_XXXX, DIFF_XXYY, DIFF_YYXX, DIFF_YYYY
	
	STIF_K(:,:) = 0.D0; LOAD_F(:) = 0.D0
            
	!! ---------------------------- [[[[[ OMEGA_HAT(1) VS. OMEGA_HAT(1) ]]]]] ----------------------------
	!! ---------------------------- [[[[[ OMEGA_HAT(2) VS. OMEGA_HAT(2) ]]]]] ----------------------------
    DO PATCH_ROW = 1, NUMPATCH
        
        IF (PATCH_ROW==1) THEN 
            NUMGSPT = 2*((BASIS_KVEC(PATCH_ROW, 1)%POLY_ORDER+1) + (BASIS_KVEC(PATCH_ROW, 2)%POLY_ORDER+1+EXPO_ARRAY(UBOUND(EXPO_ARRAY,1))) + PUORDER)
            ALLOCATE(INDX_ROW((BASIS_KVEC(PATCH_ROW, 1)%POLY_ORDER+1)*(BASIS_KVEC(PATCH_ROW, 2)%POLY_ORDER+1+UBOUND(EXPO_ARRAY,1))), SF_ROW(NUMGSPT**2, (BASIS_KVEC(PATCH_ROW, 1)%POLY_ORDER+1)*(BASIS_KVEC(PATCH_ROW, 2)%POLY_ORDER+1+UBOUND(EXPO_ARRAY,1))))
        ELSEIF (PATCH_ROW==2) THEN 
            NUMGSPT = 2*((BASIS_KVEC(PATCH_ROW, 1)%POLY_ORDER+1) + (BASIS_KVEC(PATCH_ROW, 2)%POLY_ORDER+1) + PUORDER)
            ALLOCATE(INDX_ROW((BASIS_KVEC(PATCH_ROW, 1)%POLY_ORDER+1)*(BASIS_KVEC(PATCH_ROW, 2)%POLY_ORDER+1)), &
                                SF_ROW(NUMGSPT**2, (BASIS_KVEC(PATCH_ROW, 1)%POLY_ORDER+1)*(BASIS_KVEC(PATCH_ROW, 2)%POLY_ORDER+1)))
        ENDIF
        N = NUMGSPT
! 		print*, 'here-2'
		ALLOCATE(GSX(NUMGSPT), GSXW(NUMGSPT), GSY(NUMGSPT), GSYW(NUMGSPT))
		ALLOCATE(WEIGHT(NUMGSPT**2), LOAD_Q_VEC(NUMGSPT**2))
		
        IF (PATCH_ROW==1) THEN  
            KSPAN_NDX(1) = (BASIS_KVEC(PATCH_ROW, 1)%POLY_ORDER+1)*(BASIS_KVEC(PATCH_ROW, 2)%POLY_ORDER+1+UBOUND(EXPO_ARRAY,1))
            KSPAN_NDX(2) = (BASIS_KVEC(PATCH_ROW, 1)%POLY_ORDER+1)*(BASIS_KVEC(PATCH_ROW, 2)%POLY_ORDER+1+UBOUND(EXPO_ARRAY,1))
        ELSEIF (PATCH_ROW==2) THEN
            KSPAN_NDX(1) = (BASIS_KVEC(PATCH_ROW, 1)%POLY_ORDER+1)*(BASIS_KVEC(PATCH_ROW, 2)%POLY_ORDER+1)
            KSPAN_NDX(2) = (BASIS_KVEC(PATCH_ROW, 1)%POLY_ORDER+1)*(BASIS_KVEC(PATCH_ROW, 2)%POLY_ORDER+1)
        ENDIF
        
		DO JJ=1, NUMIR(PATCH_ROW, PATCH_ROW, 2)-1
			DO II=1, NUMIR(PATCH_ROW, PATCH_ROW, 1)-1
				IRBOX%PT1 = POINT2D(IR_GRID(PATCH_ROW, PATCH_ROW, 1, II), IR_GRID(PATCH_ROW, PATCH_ROW, 2, JJ))
				IRBOX%PT2 = POINT2D(IR_GRID(PATCH_ROW, PATCH_ROW, 1, II+1), IR_GRID(PATCH_ROW, PATCH_ROW, 2, JJ))
				IRBOX%PT3 = POINT2D(IR_GRID(PATCH_ROW, PATCH_ROW, 1, II+1), IR_GRID(PATCH_ROW, PATCH_ROW, 2, JJ+1))
				IRBOX%PT4 = POINT2D(IR_GRID(PATCH_ROW, PATCH_ROW, 1, II), IR_GRID(PATCH_ROW, PATCH_ROW, 2, JJ+1))
				!------------------! VECTORIZE GAUSS POINTS AND WEIGHTS------------------------------
				CALL GAULEG(IRBOX%PT1%X, IRBOX%PT2%X, GSX, GSXW, NUMGSPT)
				CALL GAULEG(IRBOX%PT1%Y, IRBOX%PT4%Y, GSY, GSYW, NUMGSPT)
				KK = 0
				DO I = 1, NUMGSPT
					DO J = 1, NUMGSPT
						KK = KK + 1
						GSPT = POINT2D(GSX(I), GSY(J))
						JACOB(1) = GET_JACOBIAN_MATRIX(GSPT, PATCH_ROW)
						DET_M = .DETERMINANT.JACOB(1)
						WEIGHT(KK) = DET_M*GSXW(I)*GSYW(J)
						LOAD_Q_VEC(KK) = LDFT2D(GSPT, PATCH_ROW)
                        IF (PATCH_ROW==1) THEN 
                            CALL GET_PATCHWISE_BS_HYBRID(SF_ROW(KK, :), INDX_ROW(:), GSPT, PATCH_ROW, JACOB(1))     ! GET_PATCHWISE_BS(SF, INDX, PAR_PT, PATCH, JACOB)
                        ELSEIF (PATCH_ROW==2) THEN 
                            CALL GET_PATCHWISE_BS(SF_ROW(KK, :), INDX_ROW(:), GSPT, PATCH_ROW, JACOB(1))     ! GET_PATCHWISE_BS(SF, INDX, PAR_PT, PATCH, JACOB)
                        ENDIF
					ENDDO
				ENDDO
				!--------------------EVALUATE INTEGRAND AT EACH GAUSS POINT IN SUBMATRIX-----------!
				! --------------------- [[[[[ ROW ]]]]] -------------------------------!
				DO I = 1, KSPAN_NDX(1)
					DO KK=1, DOF
						IF (NDX(1,KK)==PATCH_ROW .AND. NDX(2,KK)==INDX_ROW(I)%A .AND. NDX(3,KK)==INDX_ROW(I)%B) THEN
							GLOBAL_INDX(1) = KK
							GOTO 111
						ENDIF
					ENDDO
					GLOBAL_INDX(1) = 0
					111 CONTINUE
				! ------------------- [[[[[ COLUMN ]]]]] -------------------------------
					DO J = 1, KSPAN_NDX(2)
						DO KK=1, DOF
							IF (NDX(1,KK)==PATCH_ROW .AND. NDX(2,KK)==INDX_ROW(J)%A .AND. NDX(3,KK)==INDX_ROW(J)%B) THEN
								GLOBAL_INDX(2) = KK
								GOTO 222
							ENDIF
						ENDDO
						GLOBAL_INDX(2) = 0
						222 CONTINUE
				!-----------------------------------------------------------------------
						IF (PDE=='ELPT' .AND. GLOBAL_INDX(1).NE.0 .AND. GLOBAL_INDX(2).NE.0) THEN
						! Poisson
						!-----------------------------------------------------------------------------------------------------------------------------
		! 					DIFF_NN = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D00, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_ROW(:,J)%D00)
							DIFF_XX = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D10, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_ROW(:,J)%D10)
							DIFF_YY = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D01, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_ROW(:,J)%D01)
		! 					DIFF_NX = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D00, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_ROW(:,J)%D10)
		! 					DIFF_NY = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D00, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_ROW(:,J)%D01)
		! 					DIFF_YX = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D01, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_ROW(:,J)%D10)
		! 					DIFF_XY = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D10, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_ROW(:,J)%D01)
							STIF_K(GLOBAL_INDX(1),GLOBAL_INDX(2)) = STIF_K(GLOBAL_INDX(1),GLOBAL_INDX(2)) + DIFF_XX + DIFF_YY
						!-----------------------------------------------------------------------------------------------------------------------------
						ELSEIF (PDE=='CNVD' .AND. GLOBAL_INDX(1).NE.0 .AND. GLOBAL_INDX(2).NE.0) THEN
						! Convection-Diffusion
						!-----------------------------------------------------------------------------------------------------------------------------
	! 						DIFF_NN = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D00, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_ROW(:,J)%D00)
							DIFF_XX = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D10, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_ROW(:,J)%D10)
							DIFF_YY = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D01, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_ROW(:,J)%D01)
		! 					DIFF_NX = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D00, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), &
		! 																PRODUCT(RESHAPE((/FT_B(:),   SF_ROW(:,J)%D10/), (/N**2, 2/), ORDER=ORDER1), 2))
! 							DIFF_YN = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D01, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_ROW(:,J)%D00)
							DIFF_NY = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D00, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_ROW(:,J)%D01)
		! 					DIFF_YX = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D01, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_ROW(:,J)%D10)
		! 					DIFF_XY = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D10, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_ROW(:,J)%D01)
							
							STIF_K(GLOBAL_INDX(1),GLOBAL_INDX(2)) = STIF_K(GLOBAL_INDX(1),GLOBAL_INDX(2)) + (EPSLN*(DIFF_XX + DIFF_YY) - DIFF_NY)
		!					STIF_K(GLOBAL_INDX(1),GLOBAL_INDX(2)) = STIF_K(GLOBAL_INDX(1),GLOBAL_INDX(2)) + (EPSLN*(DIFF_XX + DIFF_YY) - DIFF_NX + 1.5D0*DIFF_NN)
		! 					STIF_K(GLOBAL_INDX(2),GLOBAL_INDX(1)) = STIF_K(GLOBAL_INDX(2),GLOBAL_INDX(1)) + DIFF_XX + DIFF_YY
		! 					STIF_K(GLOBAL_INDX(1),GLOBAL_INDX(2)) = STIF_K(GLOBAL_INDX(1),GLOBAL_INDX(2)) + ((DIFF_XX + DIFF_YY) + DIFF_NN)
		! 					STIF_K(GLOBAL_INDX(1),GLOBAL_INDX(2)) = STIF_K(GLOBAL_INDX(1),GLOBAL_INDX(2)) + ((DIFF_XX + DIFF_YY) - DIFF_NX + DIFF_NN)
						!-----------------------------------------------------------------------------------------------------------------------------
                        ELSEIF (PDE=='FOUR' .AND. GLOBAL_INDX(1).NE.0 .AND. GLOBAL_INDX(2).NE.0) THEN
                        ! 4th order PDE
                            DIFF_XXXX = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D20, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_ROW(:,J)%D20)
                            DIFF_XXYY = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D20, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_ROW(:,J)%D02)
                            DIFF_YYXX = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D02, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_ROW(:,J)%D20)
                            DIFF_YYYY = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D02, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_ROW(:,J)%D02)
                            
                            STIF_K(GLOBAL_INDX(1),GLOBAL_INDX(2)) = STIF_K(GLOBAL_INDX(1),GLOBAL_INDX(2)) + (DIFF_XXXX + DIFF_XXYY  + DIFF_YYXX + DIFF_YYYY)
                            
!                             if (GLOBAL_INDX(1)==78 .and. GLOBAL_INDX(2)==78) then 
!                                 DO KK = 1, N**2
!                                     print*, kk, weight(kk), i, sf_row(kk,i)%d20, j, sf_row(kk,j)%d02
!                                 ENDDO
!                             endif
                            
						ENDIF
					ENDDO
					IF (GLOBAL_INDX(1).NE.0) THEN 
                        DIFF_NN = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D00, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), LOAD_Q_VEC(:))
                        LOAD_F(GLOBAL_INDX(1)) = LOAD_F(GLOBAL_INDX(1)) + DIFF_NN
                    ENDIF
				ENDDO
			ENDDO
		ENDDO
		DEALLOCATE(INDX_ROW, SF_ROW, GSX, GSY, GSXW, GSYW, WEIGHT, LOAD_Q_VEC)
	ENDDO
	
	
    DO PATCH = 1, 2
        IF (PATCH==1) THEN
            PATCH_ROW = 1; PATCH_COLUMN = 2
            KSPAN_NDX(PATCH_ROW) = (BASIS_KVEC(PATCH_ROW, 1)%POLY_ORDER+1)*(BASIS_KVEC(PATCH_ROW, 2)%POLY_ORDER+1+UBOUND(EXPO_ARRAY,1))
            KSPAN_NDX(PATCH_COLUMN) = (BASIS_KVEC(PATCH_COLUMN, 1)%POLY_ORDER+1)*(BASIS_KVEC(PATCH_COLUMN, 2)%POLY_ORDER+1)
        ELSEIF (PATCH==2) THEN
            PATCH_ROW = 2; PATCH_COLUMN = 1
            KSPAN_NDX(PATCH_ROW) = (BASIS_KVEC(PATCH_ROW, 1)%POLY_ORDER+1)*(BASIS_KVEC(PATCH_ROW, 2)%POLY_ORDER+1)
            KSPAN_NDX(PATCH_COLUMN) = (BASIS_KVEC(PATCH_COLUMN, 1)%POLY_ORDER+1)*(BASIS_KVEC(PATCH_COLUMN, 2)%POLY_ORDER+1+UBOUND(EXPO_ARRAY,1))
        ENDIF
        NUMGSPT = (BASIS_KVEC(PATCH_ROW, 1)%POLY_ORDER+1) + (BASIS_KVEC(PATCH_ROW, 2)%POLY_ORDER+1) + (BASIS_KVEC(PATCH_COLUMN, 1)%POLY_ORDER+1) + (BASIS_KVEC(PATCH_COLUMN, 2)%POLY_ORDER+1) + EXPO_ARRAY(UBOUND(EXPO_ARRAY,1))
    
		N = NUMGSPT
		
		ALLOCATE(INDX_ROW(KSPAN_NDX(PATCH_ROW)), INDX_COLUMN(KSPAN_NDX(PATCH_COLUMN)),  SF_ROW(NUMGSPT**2, KSPAN_NDX(PATCH_ROW)), &
						SF_COLUMN(NUMGSPT**2, KSPAN_NDX(PATCH_COLUMN)))
		
		ALLOCATE(GSX(NUMGSPT), GSXW(NUMGSPT), GSY(NUMGSPT), GSYW(NUMGSPT))
		
		ALLOCATE(WEIGHT(NUMGSPT**2))
		
		DO JJ=1, NUMIR(PATCH_ROW, PATCH_COLUMN, 2)-1
			DO II=1, NUMIR(PATCH_ROW, PATCH_COLUMN, 1)-1
				IRBOX%PT1 = POINT2D(IR_GRID(PATCH_ROW, PATCH_COLUMN, 1, II), IR_GRID(PATCH_ROW, PATCH_COLUMN, 2, JJ))
				IRBOX%PT2 = POINT2D(IR_GRID(PATCH_ROW, PATCH_COLUMN, 1, II+1), IR_GRID(PATCH_ROW, PATCH_COLUMN, 2, JJ))
				IRBOX%PT3 = POINT2D(IR_GRID(PATCH_ROW, PATCH_COLUMN, 1, II+1), IR_GRID(PATCH_ROW, PATCH_COLUMN, 2, JJ+1))
				IRBOX%PT4 = POINT2D(IR_GRID(PATCH_ROW, PATCH_COLUMN, 1, II), IR_GRID(PATCH_ROW, PATCH_COLUMN, 2, JJ+1))
				!------------------! VECTORIZE GAUSS POINTS AND WEIGHTS------------------------------
				CALL GAULEG(IRBOX%PT1%X, IRBOX%PT2%X, GSX, GSXW, NUMGSPT)
				CALL GAULEG(IRBOX%PT1%Y, IRBOX%PT4%Y, GSY, GSYW, NUMGSPT)
				KK = 0
				DO I = 1, NUMGSPT
					DO J = 1, NUMGSPT
						KK = KK + 1
						GSPT = POINT2D(GSX(I), GSY(J))
                        !---------------------------------------------!
                        PHYPT = GET_PHY_PT(GSPT, 2)
                        PARPT = GET_PAR_PT(PHYPT, 1)
                        !---------------------------------------------!
						JACOB(1) = GET_JACOBIAN_MATRIX(parpt, 1)
						JACOB(2) = GET_JACOBIAN_MATRIX(gspt, 2)
						DET_M = .DETERMINANT.JACOB(2)
						WEIGHT(KK) = DET_M*GSXW(I)*GSYW(J)
                        IF (PATCH==1) THEN
                            CALL GET_PATCHWISE_BS_HYBRID(SF_ROW(KK, :), INDX_ROW(:), parpt, PATCH_ROW, JACOB(1))     ! GET_PATCHWISE_BS(SF, INDX, PAR_PT, PATCH, JACOB)
                            CALL GET_PATCHWISE_BS(SF_COLUMN(KK, :), INDX_COLUMN(:), gspt, PATCH_COLUMN, JACOB(2))     ! GET_PATCHWISE_BS(SF, INDX, PAR_PT, PATCH, JACOB)
                        ELSEIF (PATCH==2) THEN
                            CALL GET_PATCHWISE_BS(SF_ROW(KK, :), INDX_ROW(:), gspt, PATCH_ROW, JACOB(2))     ! GET_PATCHWISE_BS(SF, INDX, PAR_PT, PATCH, JACOB)
                            CALL GET_PATCHWISE_BS_HYBRID(SF_COLUMN(KK, :), INDX_COLUMN(:), parpt, PATCH_COLUMN, JACOB(1))     ! GET_PATCHWISE_BS(SF, INDX, PAR_PT, PATCH, JACOB)
                        ENDIF
					ENDDO
				ENDDO
				!--------------------EVALUATE INTEGRAND AT EACH GAUSS POINT IN SUBMATRIX-----------
				! --------------------- [[[[[ ROW ]]]]] -------------------------------
				DO I = 1, KSPAN_NDX(PATCH_ROW)
					DO KK=1, DOF
						IF (NDX(1,KK)==PATCH_ROW .AND. NDX(2,KK)==INDX_ROW(I)%A .AND. NDX(3,KK)==INDX_ROW(I)%B) THEN
							GLOBAL_INDX(1) = KK
							GOTO 333
						ENDIF
					ENDDO
					GLOBAL_INDX(1) = 0
					333 CONTINUE
				!-----------------------------------------------------------------------
				! ------------------- [[[[[ COLUMN ]]]]] -------------------------------
					DO J = 1, KSPAN_NDX(PATCH_COLUMN)
						DO KK=1, DOF
							IF (NDX(1,KK)==PATCH_COLUMN .AND. NDX(2,KK)==INDX_COLUMN(J)%A .AND. NDX(3,KK)==INDX_COLUMN(J)%B) THEN
								GLOBAL_INDX(2) = KK
								GOTO 444
							ENDIF
						ENDDO
						GLOBAL_INDX(2) = 0
						444 CONTINUE
				!-----------------------------------------------------------------------
						IF (PDE=='ELPT' .AND. GLOBAL_INDX(1).NE.0 .AND. GLOBAL_INDX(2).NE.0) THEN
						! Poisson
						!-----------------------------------------------------------------------------------------------------------------------------
		! 					DIFF_NN = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D00, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_ROW(:,J)%D00)
							DIFF_XX = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D10, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_COLUMN(:,J)%D10)
							DIFF_YY = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D01, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_COLUMN(:,J)%D01)
		! 					DIFF_NX = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D00, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_ROW(:,J)%D10)
		! 					DIFF_NY = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D00, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_ROW(:,J)%D01)
		! 					DIFF_YX = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D01, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_ROW(:,J)%D10)
		! 					DIFF_XY = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D10, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_ROW(:,J)%D01)
							STIF_K(GLOBAL_INDX(1),GLOBAL_INDX(2)) = STIF_K(GLOBAL_INDX(1),GLOBAL_INDX(2)) + DIFF_XX + DIFF_YY
						!-----------------------------------------------------------------------------------------------------------------------------
						ELSEIF (PDE=='CNVD' .AND. GLOBAL_INDX(1).NE.0 .AND. GLOBAL_INDX(2).NE.0) THEN
						! Convection-Diffusion
						!-----------------------------------------------------------------------------------------------------------------------------
	! 						DIFF_NN = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D00, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_ROW(:,J)%D00)
							DIFF_XX = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D10, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_COLUMN(:,J)%D10)
							DIFF_YY = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D01, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_COLUMN(:,J)%D01)
		! 					DIFF_NX = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D00, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), &
		! 																PRODUCT(RESHAPE((/FT_B(:),   SF_ROW(:,J)%D10/), (/N**2, 2/), ORDER=ORDER1), 2))
! 							DIFF_YN = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D01, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_COLUMN(:,J)%D00)
							DIFF_NY = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D00, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_COLUMN(:,J)%D01)
		! 					DIFF_YX = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D01, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_ROW(:,J)%D10)
		! 					DIFF_XY = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D10, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_ROW(:,J)%D01)
							
							STIF_K(GLOBAL_INDX(1),GLOBAL_INDX(2)) = STIF_K(GLOBAL_INDX(1),GLOBAL_INDX(2)) + (EPSLN*(DIFF_XX + DIFF_YY) - DIFF_NY)
		!					STIF_K(GLOBAL_INDX(1),GLOBAL_INDX(2)) = STIF_K(GLOBAL_INDX(1),GLOBAL_INDX(2)) + (EPSLN*(DIFF_XX + DIFF_YY) - DIFF_NX + 1.5D0*DIFF_NN)
		! 					STIF_K(GLOBAL_INDX(2),GLOBAL_INDX(1)) = STIF_K(GLOBAL_INDX(2),GLOBAL_INDX(1)) + DIFF_XX + DIFF_YY
		! 					STIF_K(GLOBAL_INDX(1),GLOBAL_INDX(2)) = STIF_K(GLOBAL_INDX(1),GLOBAL_INDX(2)) + ((DIFF_XX + DIFF_YY) + DIFF_NN)
		! 					STIF_K(GLOBAL_INDX(1),GLOBAL_INDX(2)) = STIF_K(GLOBAL_INDX(1),GLOBAL_INDX(2)) + ((DIFF_XX + DIFF_YY) - DIFF_NX + DIFF_NN)
						!-----------------------------------------------------------------------------------------------------------------------------
                        ELSEIF (PDE=='FOUR' .AND. GLOBAL_INDX(1).NE.0 .AND. GLOBAL_INDX(2).NE.0) THEN
                        ! 4th order PDE
                            DIFF_XXXX = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D20, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_COLUMN(:,J)%D20)
                            DIFF_XXYY = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D20, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_COLUMN(:,J)%D02)
                            DIFF_YYXX = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D02, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_COLUMN(:,J)%D20)
                            DIFF_YYYY = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D02, WEIGHT(:)/), (/N**2, 2/), ORDER=ORDER1), 2), SF_COLUMN(:,J)%D02)
                            
                            STIF_K(GLOBAL_INDX(1),GLOBAL_INDX(2)) = STIF_K(GLOBAL_INDX(1),GLOBAL_INDX(2)) + (DIFF_XXXX + DIFF_XXYY + DIFF_YYXX + DIFF_YYYY)
						ENDIF
					ENDDO
				ENDDO
			ENDDO
		ENDDO
		DEALLOCATE(INDX_ROW, INDX_COLUMN, SF_ROW, SF_COLUMN, GSX, GSY, GSXW, GSYW, WEIGHT)
	ENDDO
	
	WRITE(*,*)
	WRITE(*,*) '<<< ASSEMBLE STIFFNESS MATRIX AND LOAD VECTOR : DONE >>>'
	WRITE(*,*)
	
END SUBROUTINE GEN_KF





! SUBROUTINE GEN_BD_LN(SUB_K, SUB_F)
! 
! 	REAL*8, INTENT(OUT) :: SUB_K(BDNDX(1)%LC_NUM, BDNDX(1)%LC_NUM)
! 	REAL*8, INTENT(OUT) :: SUB_F(BDNDX(1)%LC_NUM)
! 	
! 	REAL*8 :: DIFF_NN, DET_M
! 	REAL*8, ALLOCATABLE :: EXACT_ON_BD(:), WEIGHT(:)
! 	TYPE(FVALUE), ALLOCATABLE :: SF_ROW(:, :)
! 	TYPE(INT2D) :: INDX_ROW((BASIS_KVEC(1,3,1)%POLY_ORDER+1)*(BASIS_KVEC(1,3,2)%POLY_ORDER+1))
! 	TYPE(POINT2D) :: GSPT, DMY_GSPT
! 	INTEGER :: I, J, II, JJ, KK, N, GLOBAL_INDX(2)
! 	TYPE(VEC2D) :: SUB_LN
! 	TYPE(MATRIX_22) :: JACOB
! 	CHARACTER(LEN=1) :: GAMMA_HAT
! 	
! 	IF (USE_ENRICH=='Y') THEN
! ! 		NUMGSPT = 2*(PUORDER + SUM(BS_ORDER(1, 3, :))) + 10
! 		NUMGSPT = 80
! 	ELSE 
! 		NUMGSPT = 2*(PUORDER + SUM(BS_ORDER(1, 3, :)))
! 	ENDIF
! 	
! 	ALLOCATE(EXACT_ON_BD(NUMGSPT), WEIGHT(NUMGSPT))
! 	ALLOCATE(SF_ROW(NUMGSPT, (BASIS_KVEC(1,3,1)%POLY_ORDER+1)*(BASIS_KVEC(1,3,2)%POLY_ORDER+1)))
! 	ALLOCATE(GSX(NUMGSPT), GSXW(NUMGSPT))
! 	
! 	N = NUMGSPT; SUB_K(:,:)=0.D0; SUB_F(:) = 0.D0
! 
! ! !------------------------- LEFT ----------------------------------
! ! 	DO JJ = 1, NUMIR(2)-1
! ! 		SUB_LN = VEC2D(POINT2D(0.D0,IR_GRID(2,JJ)), POINT2D(0.D0,IR_GRID(2,JJ+1)))
! ! 		DIFF_NN = 0.D0
! ! 		CALL GAULEG(SUB_LN%U%Y, SUB_LN%V%Y, GSY, GSYW, NUMGSPT)
! ! 		DO I=1,N
! ! 			GSPT = POINT2D(0.D0, GSY(I))
! ! 			DMY_GSPT = GSPT
! ! 			JACOB = GET_JACOBIAN_MATRIX(GSPT)
! ! 			DET_M = .DETERMINANT.JACOB
! ! ! 			GAMMA_HAT = GET_GAMMA_HAT(GSPT)
! ! ! 			DET_M = GET_DET_DR(GSPT, GAMMA_HAT)
! ! 			WEIGHT(I) = DET_M*GSYW(I)
! ! 			CALL GET_ALL_PHY_BASIS_IN_INT(SF_ROW(I,:), INDX, DMY_GSPT, JACOB)
! ! ! 			CALL GET_ALL_PHY_NURBS_SURFACE_2D(SF_ROW(I,:), INDX, DMY_GSPT, JACOB)
! ! 			EXACT_ON_BD(I) = EX_DISP(DMY_GSPT)
! ! 		ENDDO
! ! ! 		PRINT*, 'HERE-1'
! ! 		DO I = 1, (BASIS_KVEC(1)%POLY_ORDER+1)*(BASIS_KVEC(2)%POLY_ORDER+1)
! ! 			DO KK=1, BDNDX(1)%LC_NUM
! ! 				IF (BDNDX(KK)%LC_NDX(1).EQ.INDX_ROW(I)%A .AND. BDNDX(KK)%LC_NDX(2).EQ.INDX_ROW(I)%B) THEN
! ! 					GLOBAL_INDX(1) = KK
! ! 					GOTO 111
! ! 				ENDIF
! ! 			ENDDO
! ! 			GOTO 112
! ! ! 			print*, 'here-2'
! ! 			111 CONTINUE
! ! 			DO J = 1, (BASIS_KVEC(1)%POLY_ORDER+1)*(BASIS_KVEC(2)%POLY_ORDER+1)
! ! 				DO KK=1, BDNDX(1)%LC_NUM
! ! 					IF (BDNDX(KK)%LC_NDX(1).EQ.INDX_ROW(J)%A .AND. BDNDX(KK)%LC_NDX(2).EQ.INDX_ROW(J)%B) THEN
! ! 						GLOBAL_INDX(2) = KK
! ! 						GOTO 222
! ! 					ENDIF
! ! 				ENDDO
! ! 				GOTO 223
! ! ! 				print*, 'here-3'
! ! 				222 CONTINUE
! ! 				DIFF_NN = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D00, WEIGHT(:)/), (/N, 2/), ORDER=ORDER1), 2), SF_ROW(:,J)%D00)
! ! 				SUB_K(GLOBAL_INDX(1), GLOBAL_INDX(2)) = SUB_K(GLOBAL_INDX(1), GLOBAL_INDX(2)) + DIFF_NN
! ! 				223 CONTINUE
! ! 			ENDDO
! ! 			
! ! 			DIFF_NN = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D00, WEIGHT(:)/), (/N, 2/), ORDER=ORDER1), 2), EXACT_ON_BD(:))
! ! 			SUB_F(GLOBAL_INDX(1)) = SUB_F(GLOBAL_INDX(1)) + DIFF_NN
! ! 			112 CONTINUE
! ! 		ENDDO
! ! ! 		PRINT*, 'HERE-4'
! ! 	ENDDO
! ! 
! ! !------------------------- RIGHT ----------------------------------
! ! 	DO JJ = 1, NUMIR(2)-1
! ! 		SUB_LN = VEC2D(POINT2D(1.D0,IR_GRID(2,JJ)), POINT2D(1.D0,IR_GRID(2,JJ+1)))
! ! 		DIFF_NN = 0.D0
! ! 		CALL GAULEG(SUB_LN%U%Y, SUB_LN%V%Y, GSY, GSYW, NUMGSPT)
! ! 		DO I=1,N
! ! 			GSPT = POINT2D(1.D0, GSY(I))
! ! 			DMY_GSPT = GSPT
! ! 			JACOB = GET_JACOBIAN_MATRIX(GSPT)
! ! 			DET_M = .DETERMINANT.JACOB
! ! ! 			GAMMA_HAT = GET_GAMMA_HAT(GSPT)
! ! ! 			DET_M = GET_DET_DR(GSPT, GAMMA_HAT)
! ! 			WEIGHT(I) = DET_M*GSYW(I)
! ! 			CALL GET_ALL_PHY_BASIS_IN_INT(SF_ROW(I,:), INDX, DMY_GSPT, JACOB)
! ! ! 			CALL GET_ALL_PHY_NURBS_SURFACE_2D(SF_ROW(I,:), INDX, DMY_GSPT, JACOB)
! ! 			EXACT_ON_BD(I) = EX_DISP(DMY_GSPT)
! ! 		ENDDO
! ! 		DO I = 1, (BASIS_KVEC(1)%POLY_ORDER+1)*(BASIS_KVEC(2)%POLY_ORDER+1)
! ! 			DO KK=1, BDNDX(1)%LC_NUM
! ! 				IF (BDNDX(KK)%LC_NDX(1).EQ.INDX_ROW(I)%A .AND. BDNDX(KK)%LC_NDX(2).EQ.INDX_ROW(I)%B) THEN
! ! 					GLOBAL_INDX(1) = KK
! ! 					GOTO 333
! ! 				ENDIF
! ! 			ENDDO
! ! 			GOTO 334
! ! 			333 CONTINUE
! ! 			DO J = 1, (BASIS_KVEC(1)%POLY_ORDER+1)*(BASIS_KVEC(2)%POLY_ORDER+1)
! ! 				DO KK=1, BDNDX(1)%LC_NUM
! ! 					IF (BDNDX(KK)%LC_NDX(1).EQ.INDX_ROW(J)%A .AND. BDNDX(KK)%LC_NDX(2).EQ.INDX_ROW(J)%B) THEN
! ! 						GLOBAL_INDX(2) = KK
! ! 						GOTO 444
! ! 					ENDIF
! ! 				ENDDO
! ! 				GOTO 445
! ! 				444 CONTINUE
! ! 				DIFF_NN = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D00, WEIGHT(:)/), (/N, 2/), ORDER=ORDER1), 2), SF_ROW(:,J)%D00)
! ! 				SUB_K(GLOBAL_INDX(1), GLOBAL_INDX(2)) = SUB_K(GLOBAL_INDX(1), GLOBAL_INDX(2)) + DIFF_NN
! ! 				445 CONTINUE
! ! 			ENDDO
! ! 			DIFF_NN = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D00, WEIGHT(:)/), (/N, 2/), ORDER=ORDER1), 2), EXACT_ON_BD(:))
! ! 			SUB_F(GLOBAL_INDX(1)) = SUB_F(GLOBAL_INDX(1)) + DIFF_NN
! ! 			334 CONTINUE
! ! 		ENDDO
! ! 	ENDDO
! 
! !------------------------- BOTTOM ----------------------------------
! 	DO JJ = 1, NUMIR(1,3,1)-1
! 		SUB_LN = VEC2D(POINT2D(IR_GRID(1,3,1,JJ),0.D0), POINT2D(IR_GRID(1,3,1,JJ+1),0.D0))
! 		DIFF_NN = 0.D0
! 		CALL GAULEG(SUB_LN%U%X, SUB_LN%V%X, GSX, GSXW, NUMGSPT)
! 		DO I=1,N
! 			GSPT = POINT2D(GSX(I),0.D0)
! 			JACOB = GET_JACOBIAN_MATRIX(GSPT)
! 			DET_M = .DETERMINANT.JACOB
! ! 			GAMMA_HAT = GET_GAMMA_HAT(GSPT)
! ! 			DET_M = GET_DET_DR(GSPT, GAMMA_HAT)
! 			WEIGHT(I) = DET_M*GSXW(I)
! 			CALL GET_PATCHWISE_BS(SF_ROW(I, :), INDX_ROW(:), GSPT, (/1, 3/), JACOB)
! ! 			CALL GET_ALL_PHY_BASIS_IN_INT(SF_ROW(I,:), INDX, GSPT, JACOB)
! ! 			CALL GET_ALL_PHY_NURBS_SURFACE_2D(SF_ROW(I,:), INDX, GSPT, JACOB)
! 			EXACT_ON_BD(I) = EX_DISP(GSPT)
! 		ENDDO
! 		DO I = 1, (BASIS_KVEC(1,3,1)%POLY_ORDER+1)*(BASIS_KVEC(1,3,2)%POLY_ORDER+1)
! 			DO KK=1, BDNDX(1)%LC_NUM
! 				IF (BDNDX(KK)%LC_NDX(1)==1 .AND. BDNDX(KK)%LC_NDX(2)==3 .AND. BDNDX(KK)%LC_NDX(3).EQ.INDX_ROW(I)%A .AND. BDNDX(KK)%LC_NDX(4).EQ.INDX_ROW(I)%B) THEN
! 					GLOBAL_INDX(1) = KK
! 					GOTO 555
! 				ENDIF
! 			ENDDO
! 			GOTO 556
! 			555 CONTINUE
! 			DO J = 1, (BASIS_KVEC(1,3,1)%POLY_ORDER+1)*(BASIS_KVEC(1,3,2)%POLY_ORDER+1)
! 				DO KK=1, BDNDX(1)%LC_NUM
! 				IF (BDNDX(KK)%LC_NDX(1)==1 .AND. BDNDX(KK)%LC_NDX(2)==3 .AND. BDNDX(KK)%LC_NDX(3).EQ.INDX_ROW(J)%A .AND. BDNDX(KK)%LC_NDX(4).EQ.INDX_ROW(J)%B) THEN
! 						GLOBAL_INDX(2) = KK
! 						GOTO 666
! 					ENDIF
! 				ENDDO
! 				GOTO 667
! 				666 CONTINUE
! 				DIFF_NN = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D00, WEIGHT(:)/), (/N, 2/), ORDER=ORDER1), 2), SF_ROW(:,J)%D00)
! 				SUB_K(GLOBAL_INDX(1), GLOBAL_INDX(2)) = SUB_K(GLOBAL_INDX(1), GLOBAL_INDX(2)) + DIFF_NN
! 				667 CONTINUE
! 			ENDDO
! 			DIFF_NN = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D00, WEIGHT(:)/), (/N, 2/), ORDER=ORDER1), 2), EXACT_ON_BD(:))
! 			SUB_F(GLOBAL_INDX(1)) = SUB_F(GLOBAL_INDX(1)) + DIFF_NN
! 			556 CONTINUE
! 		ENDDO
! 	ENDDO
! 	
! 	DEALLOCATE(EXACT_ON_BD, WEIGHT)
! 	DEALLOCATE(SF_ROW)
! 	DEALLOCATE(GSX, GSXW)
! 	
! ! !------------------------- TOP ----------------------------------
! ! 	DO JJ = 1, NUMIR(1)-1
! ! 		SUB_LN = VEC2D(POINT2D(IR_GRID(1,JJ),1.D0), POINT2D(IR_GRID(1,JJ+1),1.D0))
! ! 		DIFF_NN = 0.D0
! ! 		CALL GAULEG(SUB_LN%U%X, SUB_LN%V%X, GSX, GSXW, NUMGSPT)
! ! 		DO I=1,N
! ! 			GSPT = POINT2D(GSX(I),1.D0)
! ! 			DMY_GSPT = GSPT
! ! 			JACOB = GET_JACOBIAN_MATRIX(GSPT)
! ! 			DET_M = .DETERMINANT.JACOB
! ! ! 			GAMMA_HAT = GET_GAMMA_HAT(GSPT)
! ! ! 			DET_M = GET_DET_DR(GSPT, GAMMA_HAT)
! ! 			WEIGHT(I) = DET_M*GSXW(I)
! ! 			CALL GET_ALL_PHY_BASIS_IN_INT(SF_ROW(I,:), INDX, DMY_GSPT, JACOB)
! ! ! 			CALL GET_ALL_PHY_NURBS_SURFACE_2D(SF_ROW(I,:), INDX, DMY_GSPT, JACOB)
! ! 			EXACT_ON_BD(I) = EX_DISP(DMY_GSPT)
! ! 		ENDDO
! ! 		DO I = 1, (BASIS_KVEC(1)%POLY_ORDER+1)*(BASIS_KVEC(2)%POLY_ORDER+1)
! ! 			DO KK=1, BDNDX(1)%LC_NUM
! ! 				IF (BDNDX(KK)%LC_NDX(1).EQ.INDX_ROW(I)%A .AND. BDNDX(KK)%LC_NDX(2).EQ.INDX_ROW(I)%B) THEN
! ! 					GLOBAL_INDX(1) = KK
! ! 					GOTO 777
! ! 				ENDIF
! ! 			ENDDO
! ! 			GOTO 778
! ! 			777 CONTINUE
! ! 			DO J = 1, (BASIS_KVEC(1)%POLY_ORDER+1)*(BASIS_KVEC(2)%POLY_ORDER+1)
! ! 				DO KK=1, BDNDX(1)%LC_NUM
! ! 					IF (BDNDX(KK)%LC_NDX(1).EQ.INDX_ROW(J)%A .AND. BDNDX(KK)%LC_NDX(2).EQ.INDX_ROW(J)%B) THEN
! ! 						GLOBAL_INDX(2) = KK
! ! 						GOTO 888
! ! 					ENDIF
! ! 				ENDDO
! ! 				GOTO 889
! ! 				888 CONTINUE
! ! 				DIFF_NN = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D00, WEIGHT(:)/), (/N, 2/), ORDER=ORDER1), 2), SF_ROW(:,J)%D00)
! ! 				SUB_K(GLOBAL_INDX(1), GLOBAL_INDX(2)) = SUB_K(GLOBAL_INDX(1), GLOBAL_INDX(2)) + DIFF_NN
! ! 				889 CONTINUE
! ! 			ENDDO
! ! 			DIFF_NN = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%D00, WEIGHT(:)/), (/N, 2/), ORDER=ORDER1), 2), EXACT_ON_BD(:))
! ! 			SUB_F(GLOBAL_INDX(1)) = SUB_F(GLOBAL_INDX(1)) + DIFF_NN
! ! 			778 CONTINUE
! ! 		ENDDO
! ! 	ENDDO
! 		
! END SUBROUTINE GEN_BD_LN

END MODULE INTEGRATION
