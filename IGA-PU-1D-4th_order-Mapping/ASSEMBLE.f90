	MODULE ASSEMBLE

	USE GSQUAD
!	USE LUDECOMPOSITION
	USE NURBS_BASIS
! 	USE INTERPOLATION
	USE LOADFUNCTION

CONTAINS

!------------------INTEGRAL CODE FOR THE STIFFNESS MATRIX ELEMENT--------------
SUBROUTINE GEN_KF(STIF_K, LOAD_F)

	REAL*8, INTENT(OUT) :: STIF_K(DOF, DOF), LOAD_F(DOF)
	
	INTEGER, ALLOCATABLE :: INDX_ROW(:), INDX_COLUMN(:)
	TYPE(FUNCTION_1D), ALLOCATABLE :: SF_ROW(:,:), SF_COLUMN(:,:)
	TYPE(FUNCTION_1D) :: JACOBI, INVGMAPF
	
	REAL*8 :: LOAD_Q_VEC(NUMGSPT)
	INTEGER :: I,J, K, II, JJ, KK, GLOBAL_INDX(2), PATCH_ROW, PATCH_COLUMN, PATCHES(2), N
	REAL*8 :: DIFF_XXXX, DIFF_NN, DIFF_TS, DIFF_XX
	
	STIF_K(:,:) = 0.D0; LOAD_F(:) = 0.D0
	N = NUMGSPT
	!! ---------------------------- [[[[[ OMEGA_HAT_F (PATCH 1) VS. OMEGA_HAT_G (PATCH 2) ]]]]] ----------------------------
	DO PATCH_ROW = 1, NUMPATCH
		
! 		print*, 'n=', n, 'numir=', numir(patch_row, patch_row)
		
		ALLOCATE(INDX_ROW(BS_KVEC(PATCH_ROW)%POLY_ORDER+1), SF_ROW(N, BS_KVEC(PATCH_ROW)%POLY_ORDER+1))
! 		SF_ROW(:,:)%VAL(0) = 0.0D0
! 		SF_ROW(:,:)%VAL(1) = 0.0D0
! 		SF_ROW(:,:)%VAL(2) = 0.0D0
! 		INDX_ROW(:) = 0
		DO II=1, NUMIR(PATCH_ROW, PATCH_ROW)-1
			!------------------! VECTORIZE GAUSS POINTS AND WEIGHTS------------------------------
			CALL GAULEG(IR_GRID(PATCH_ROW, PATCH_ROW, II), IR_GRID(PATCH_ROW, PATCH_ROW, II + 1), GSPT, GSW, N)
! 			print*, 'ii=', ii
			DIFF_TS = 0.0D0
			DO I = 1, N
! 				PRINT*, 'GSPT', GSPT(I)
				IF (PATCH_ROW==1) THEN
					JACOBI = MAP_F(GSPT(I))
				ELSEIF (PATCH_ROW==2) THEN
					JACOBI = MAP_G(GSPT(I))
				ENDIF
				LOAD_Q_VEC(I) = LDF1D(GSPT(I), PATCH_ROW)
				GSW(I) = GSW(I)*DABS(JACOBI%VAL(1))
				CALL GET_ALL_PHY_BASIS1D(SF_ROW(I,:), INDX_ROW(:), GSPT(I), PATCH_ROW)
				
! 				IF (PATCH_ROW==1) THEN 
!                   PRINT*, PATCH_ROW, DABS(0.50D0 - DABS(JACOBI%VAL(1)))
!                 ELSEIF (PATCH_ROW==2) THEN 
!                   PRINT*, PATCH_ROW, DABS(0.60D0 - DABS(JACOBI%VAL(1)))
!                 ENDIF
! 				PRINT*, 'GSW', GSW(I)
! 				PRINT*, 'SF_ROW', SF_ROW(I, 1)%VAL(2)
! 				DIFF_TS = DIFF_TS + SF_ROW(I,1)%VAL(2)*GSW(I)
! 				PRINT*, I, DIFF_TS
! 				DO J = 1, BS_KVEC(PATCH_ROW)%POLY_ORDER+1
! 					PRINT*, SF_ROW(I, J)%VAL(2), INDX_ROW(J)
! 				ENDDO
			ENDDO
! 			PRINT*, DIFF_TS
! 			STOP
			!--------------------EVALUATE INTEGRAND AT EACH GAUSS POINT IN SUBMATRIX-----------
			! --------------------- [[[[[ ROW ]]]]] -------------------------------
			DO I = 1, BS_KVEC(PATCH_ROW)%POLY_ORDER+1
! 				print*, 'i=', i
				DO KK=1, DOF
					IF (NDX(1,KK)==PATCH_ROW .AND. NDX(2,KK).EQ.INDX_ROW(I)) THEN
						GLOBAL_INDX(1) = KK
						GOTO 111
					ENDIF
				ENDDO
				GLOBAL_INDX(1) = 0
				111 CONTINUE
			!-----------------------------------------------------------------------
			! ------------------- [[[[[ COLUMN ]]]]] -------------------------------
				DO J = 1, BS_KVEC(PATCH_ROW)%POLY_ORDER+1
! 					print*, 'j=', j
					DO KK=1, DOF
						IF (NDX(1,KK)==PATCH_ROW .AND. NDX(2,KK).EQ.INDX_ROW(J)) THEN
							GLOBAL_INDX(2) = KK
							GOTO 222
						ENDIF
					ENDDO
					GLOBAL_INDX(2) = 0
					222 CONTINUE
			!-----------------------------------------------------------------------
! 					PRINT*, GLOBAL_INDX(1), GLOBAL_INDX(2)
					IF (GLOBAL_INDX(1).NE.0 .AND. GLOBAL_INDX(2).NE.0) THEN
						IF (EQN=='FR') THEN 
						! 4th order
							DIFF_XXXX = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%VAL(2), GSW(:)/), (/N, 2/), ORDER=ORDER1), 2), SF_ROW(:,J)%VAL(2))
	! 						PRINT*, GLOBAL_INDX(1), GLOBAL_INDX(2), DIFF_XXXX
							STIF_K(GLOBAL_INDX(1), GLOBAL_INDX(2)) = STIF_K(GLOBAL_INDX(1), GLOBAL_INDX(2)) + DIFF_XXXX
						ELSEIF (EQN=='PS') THEN 
							DIFF_XX = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%VAL(1), GSW(:)/), (/N, 2/), ORDER=ORDER1), 2), SF_ROW(:,J)%VAL(1))
							STIF_K(GLOBAL_INDX(1), GLOBAL_INDX(2)) = STIF_K(GLOBAL_INDX(1), GLOBAL_INDX(2)) + DIFF_XX
						ENDIF
					ENDIF
				ENDDO
				DIFF_NN = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%VAL(0), GSW(:)/), (/N, 2/), ORDER=ORDER1), 2), LOAD_Q_VEC(:))
				LOAD_F(GLOBAL_INDX(1)) = LOAD_F(GLOBAL_INDX(1)) + DIFF_NN
			ENDDO
		ENDDO
		DEALLOCATE(INDX_ROW, SF_ROW)
	ENDDO
	
	
	!! ---------------------------- [[[[[ OMEGA_HAT(1) VS. OMEGA_HAT(2) ]]]]] ----------------------------
	!! ---------------------------- [[[[[ OMEGA_HAT(2) VS. OMEGA_HAT(1) ]]]]] ----------------------------
	
	DO K = 1, 2
		!! CAUTION: Specifically coded to 
		IF (K==1) THEN
			PATCHES = (/1, 2/)
		ELSEIF (K==2) THEN
			PATCHES = (/2, 1/)
		ENDIF
		
		ALLOCATE(INDX_ROW(BS_KVEC(PATCHES(1))%POLY_ORDER+1), &
								INDX_COLUMN(BS_KVEC(PATCHES(2))%POLY_ORDER+1), &
								SF_ROW(N, BS_KVEC(PATCHES(1))%POLY_ORDER+1), &
								SF_COLUMN(N, BS_KVEC(PATCHES(2))%POLY_ORDER+1))
		
		DO II=1, NUMIR(PATCHES(1), PATCHES(2))-1
			!------------------! VECTORIZE GAUSS POINTS AND WEIGHTS------------------------------
			CALL GAULEG(IR_GRID(PATCHES(1), PATCHES(2), II), IR_GRID(PATCHES(1), PATCHES(2), II + 1), GSPT, GSW, N)
			DO I = 1, N
				JACOBI = MAP_F(GSPT(I))
				GSW(I) = GSW(I)*DABS(JACOBI%VAL(1))
				INVGMAPF = MAP_INVG(JACOBI%VAL(0))
				IF (K==1) THEN 
                  CALL GET_ALL_PHY_BASIS1D(SF_ROW(I,:), INDX_ROW(:), GSPT(I), PATCHES(1))
                  CALL GET_ALL_PHY_BASIS1D(SF_COLUMN(I,:), INDX_COLUMN(:), INVGMAPF%VAL(0), PATCHES(2))
                ELSEIF (K==2) THEN 
                  CALL GET_ALL_PHY_BASIS1D(SF_ROW(I,:), INDX_ROW(:), INVGMAPF%VAL(0), PATCHES(1))
                  CALL GET_ALL_PHY_BASIS1D(SF_COLUMN(I,:), INDX_COLUMN(:), GSPT(I), PATCHES(2))
                ENDIF
				
! 				PRINT*, DABS(INVGMAPF%VAL(0) - ((0.50D0*GSPT(I) - 0.40D0)/0.60D0))
			ENDDO
			!--------------------EVALUATE INTEGRAND AT EACH GAUSS POINT IN SUBMATRIX-----------
			! --------------------- [[[[[ ROW ]]]]] -------------------------------
			DO I = 1, BS_KVEC(PATCHES(1))%POLY_ORDER+1
				DO KK=1, DOF
					IF (NDX(1,KK)==PATCHES(1) .AND. NDX(2,KK).EQ.INDX_ROW(I)) THEN
						GLOBAL_INDX(1) = KK
						GOTO 333
					ENDIF
				ENDDO
				GLOBAL_INDX(1) = 0
				333 CONTINUE
			!-----------------------------------------------------------------------
			! ------------------- [[[[[ COLUMN ]]]]] -------------------------------
				DO J = 1, BS_KVEC(PATCHES(2))%POLY_ORDER+1
					DO KK=1, DOF
						IF (NDX(1,KK)==PATCHES(2) .AND. NDX(2,KK).EQ.INDX_COLUMN(J)) THEN
							GLOBAL_INDX(2) = KK
							GOTO 444
						ENDIF
					ENDDO
					GLOBAL_INDX(2) = 0
					444 CONTINUE
			!-----------------------------------------------------------------------
					IF (GLOBAL_INDX(1).NE.0 .AND. GLOBAL_INDX(2).NE.0) THEN
						IF (EQN=='FR') THEN 
						! 4th order
							DIFF_XXXX = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%VAL(2), GSW(:)/), (/N, 2/), ORDER=ORDER1), 2), SF_COLUMN(:,J)%VAL(2))
	! 						PRINT*, GLOBAL_INDX(1), GLOBAL_INDX(2), DIFF_XXXX
							STIF_K(GLOBAL_INDX(1), GLOBAL_INDX(2)) = STIF_K(GLOBAL_INDX(1), GLOBAL_INDX(2)) + DIFF_XXXX
						ELSEIF (EQN=='PS') THEN 
							DIFF_XX = DOT_PRODUCT(PRODUCT(RESHAPE((/SF_ROW(:,I)%VAL(1), GSW(:)/), (/N, 2/), ORDER=ORDER1), 2), SF_COLUMN(:,J)%VAL(1))
							STIF_K(GLOBAL_INDX(1), GLOBAL_INDX(2)) = STIF_K(GLOBAL_INDX(1), GLOBAL_INDX(2)) + DIFF_XX
! 							PRINT*, GLOBAL_INDX(1), GLOBAL_INDX(2), DIFF_XX
						ENDIF
					ENDIF
				ENDDO
			ENDDO
		ENDDO
		DEALLOCATE(INDX_ROW, SF_ROW, INDX_COLUMN, SF_COLUMN)
	ENDDO
	
	
	WRITE(*,*)
	WRITE(*,*) '<<< ASSEMBLE STIFFNESS MATRIX AND LOAD VECTOR : DONE >>>'
	WRITE(*,*)
	
END SUBROUTINE GEN_KF


! SUBROUTINE GEN_BD_LN()
! 
! 	REAL*8 :: SUB_K(BDNDX(1)%LC_NUM, BDNDX(1)%LC_NUM)
! 	REAL*8 :: SUB_F(BDNDX(1)%LC_NUM)
! 	REAL*8 :: DD
! 	INTEGER :: I, J, II, JJ, KK, GLB_NDX, LU_INDX(BDNDX(1)%LC_NUM)
! 	
! 	TYPE(FVALUE) :: SF((BS_KVEC(1)%POLY_ORDER+1)*(BS_KVEC(2)%POLY_ORDER+1))
! 	TYPE(INT2D) :: INDX((BS_KVEC(1)%POLY_ORDER+1)*(BS_KVEC(2)%POLY_ORDER+1))
! 	TYPE(POINT2D) :: COLPT
! 	TYPE(MATRIX_22) :: JACOB
! 	CHARACTER(LEN=1) :: GAMMA_HAT
! 	
! 	IF (BDNDX(1)%LC_NUM.NE.UBOUND(BD_COL_PT,1)) THEN
! 		PRINT*, 'NUM. OF BASIS ON BD. IS NOT EQUAL TO NUM. OF COL. PT. ON BD.'
! 		PRINT*, 'NUM. OF BASIS ON BD. : ', BDNDX(1)%LC_NUM
! 		PRINT*, 'NUM. OF. COL. PT. ON BD. : ', UBOUND(BD_COL_PT,1)
! 		STOP
! 	ENDIF
! 	
! 	SUB_K(:,:)=0.D0; SUB_F(:) = 0.D0
! 
! 	DO II = 1, BDNDX(1)%LC_NUM
! 		JACOB = GET_JACOBIAN_MATRIX(BD_COL_PT(II))
! 		CALL GET_ALL_PHY_BASIS_IN_INT(SF(:), INDX, BD_COL_PT(II), JACOB)
! ! 		CALL GET_ALL_PHY_NURBS_SURFACE_2D(SF(:), INDX(:), BD_COL_PT(II), JACOB)
! 		DO J = 1, (BS_KVEC(1)%POLY_ORDER+1)*(BS_KVEC(2)%POLY_ORDER+1)
! 			DO KK = 1, BDNDX(1)%LC_NUM
! 				IF (BDNDX(KK)%LC_NDX(1).EQ.INDX(J)%A .AND. 
! BDNDX(KK)%LC_NDX(2).EQ.INDX(J)%B) THEN
! 					GLB_NDX = KK
! 					GOTO 222
! 				ELSE
! 					GLB_NDX = 0
! 				ENDIF
! 			ENDDO
! 			222 CONTINUE
! 			IF (GLB_NDX.NE.0) THEN
! 				SUB_K(II, GLB_NDX) = SF(J)%D00
! 			ENDIF
! 		ENDDO
! 		SUB_F(II) = EX_DISP(BD_COL_PT(II))
! 	ENDDO
! 	
! ! 	OPEN(124, FILE = './data/sub_k')
! ! 	OPEN(125, FILE = './data/sub_f')
! ! 	DO I=1, BDNDX(1)%LC_NUM
! ! 		WRITE(124, *) (SUB_K(I,J), J=1,BDNDX(1)%LC_NUM)
! ! 		WRITE(125, *) SUB_F(I)
! ! 	ENDDO
! ! 	CLOSE(124)
! ! 	CLOSE(125)
! 	
! 	CALL LUDCMP(SUB_K, BDNDX(1)%LC_NUM, BDNDX(1)%LC_NUM, LU_INDX, DD)
! 	CALL LUBKSB(SUB_K, BDNDX(1)%LC_NUM, BDNDX(1)%LC_NUM, LU_INDX, SUB_F)
! 	
! ! 	OPEN(126, FILE = './data/sub_s')
! ! 	DO I=1, BDNDX(1)%LC_NUM
! ! 		WRITE(126, *) SUB_F(I)
! ! 	ENDDO
! ! 	CLOSE(126)
! 	
! 	LST_SOL(:) = SUB_F(:)
! 	
! END SUBROUTINE GEN_BD_LN

END MODULE ASSEMBLE
