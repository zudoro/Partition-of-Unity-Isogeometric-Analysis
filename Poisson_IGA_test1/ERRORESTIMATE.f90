MODULE ERRORESTIMATE

	USE GEOMETRY
	USE NURBS_BASIS
	USE LOADFUNCTION
	USE GSQUAD

    IMPLICIT NONE
    
CONTAINS

!----MAX. NORM ESTIMATE----
SUBROUTINE MAXNORM(MAX_NORM, COEFF_SOL)

	REAL*8, INTENT(OUT) :: MAX_NORM(5)
	REAL*8, INTENT(IN) :: COEFF_SOL(DOF)
	TYPE(POINT2D) :: MESH, PHY_PT, REF_PT
	TYPE(INT2D) :: MAXGRID
	INTEGER :: II, I, J
	TYPE(DISPLACEMENT) :: EXDISP, APDISP, ERR_DISP, ERR_LINE_DISP, ERR_POLAR_DISP, TMP_ERR_DISP, MAXDISP
	TYPE(STRESS) :: EXSTR, APSTR, ERR_STR, TMP_ERR_STR, ERR_LINE_STR, ERR_POLAR_STR
	REAL*8 :: R, THETA
	INTEGER :: PATCH

	MAXGRID = INT2D(100,100)
	MESH=POINT2D(1.D0/MAXGRID%A, 1.D0/MAXGRID%B)
	
	ERR_STR = STRESS(0.D0,0.D0,0.D0,0.D0,0.D0)
	ERR_DISP = DISPLACEMENT(0.D0,0.D0,0.D0)
	MAXDISP = DISPLACEMENT(0.D0,0.D0,0.D0)

	OPEN(11, FILE = './data/ext_disp')
	OPEN(21, FILE = './data/app_disp')
	OPEN(31, FILE = './data/err_disp')

	DO PATCH=1, 1
		DO J=1, MAXGRID%B+1
			DO I=1, MAXGRID%A+1
				REF_PT = POINT2D(MESH%X*(I-1), MESH%Y*(J-1))
				PHY_PT = GET_PHY_PT(REF_PT)
				IF(DABS(PHY_PT%Y).LE.EPS .AND. PHY_PT%X.GE.EPS) THEN
					EXDISP%PIX = 0.D0; APDISP%PIX = 0.D0
				ELSE
					CALL APPROXSOL(APDISP%PIX, REF_PT, COEFF_SOL)
					EXDISP%PIX = EX_DISP(REF_PT)
				ENDIF
				TMP_ERR_DISP%PIX = DABS(EXDISP%PIX - APDISP%PIX)
				WRITE(11,*) PHY_PT, EXDISP%PIX
				WRITE(21,*) PHY_PT, APDISP%PIX
				WRITE(31,*) PHY_PT, TMP_ERR_DISP%PIX
				IF (TMP_ERR_DISP%PIX.GE.ERR_DISP%PIX) THEN
					ERR_DISP%PIX = TMP_ERR_DISP%PIX
				ENDIF
				IF (DABS(EXDISP%PIX).GE.MAXDISP%PIX) THEN
					MAXDISP%PIX = DABS(EXDISP%PIX)
				ENDIF
			ENDDO
			WRITE(11,*) ""
			WRITE(21,*) ""
			WRITE(31,*) ""
		ENDDO
		WRITE(11,*) ""
		WRITE(21,*) ""
		WRITE(31,*) ""
	ENDDO
	CLOSE(11)
	CLOSE(21)
	CLOSE(31)

	write(char_order(1),fmt='(i2.2)') BS_ORDER(1)
	write(char_order(2),fmt='(i2.2)') BS_ORDER(2)
	write(char_basis(1),fmt='(i2.2)') NUMBS(1)
	write(char_basis(2),fmt='(i2.2)') NUMBS(2)

! 	filename = trim('./data/pro') // char_problem // trim('_') // trim('result_p') // char_order(1) // trim('X') // char_order(2) // trim('_b') // char_basis(1) // trim('X') // char_basis(2)

! 	open(13,file=filename,status='unknown')
! 	WRITE(13,*) "ORDER OF B-SPLINE : ", BS_ORDER
! 	IF (EXTRA_KNOTS(1).EQ.0 .AND. EXTRA_KNOTS(2).EQ.0) THEN
! 		WRITE(13,*) "KNOT INSERTION : NONE"
! 	ELSE
! 		WRITE(13,*) "KNOT INSERTION : ", EXTRA_KNOTS
! 	ENDIF
! 	WRITE(13,*) "BOUNDARY CONDITION : ", BCTYPE
! 	WRITE(13,*) "NUMBER OF BASIS FTS : ", NUMBS
! 	WRITE(13,*) 'DEGREE OF FREEDOM : ', DOF - BD_DOF
! 	WRITE(13,*) 'DEGREE OF FREEDOM ON BOUNDARY : ', BD_DOF
! 	WRITE(13,*) ""
! 	WRITE(13,*) "MAXIMUM ERROR (%) : ", ERR_DISP%PIX*100.D0
! 	CLOSE(13)

	MAX_NORM(1) = ERR_DISP%PIX*100.D0/MAXDISP%PIX
	
END SUBROUTINE MAXNORM

SUBROUTINE ENRGY(ENRG, DMY_K, SOL)

	REAL*8, INTENT(IN) :: DMY_K(DOF,DOF), SOL(DOF)
	REAL*8, INTENT(OUT) :: ENRG(4)
	REAL*8 :: AP, TR
	INTEGER :: I, J

	AP = 0.50D0*DOT_PRODUCT(SOL, MATMUL(DMY_K, SOL))

	IF (PROBLEM.EQ.1) THEN
		TR = 2.D0/3.D0
	ELSEIF (PROBLEM.EQ.2) THEN
		TR = 0.785398163397448280D0
	ELSEIF (PROBLEM.EQ.3) THEN
		TR = 2.0943951023931948D0
	ENDIF

! 	CALL EXT_ENERGY(TR)

	ENRG(3) = DSQRT(DABS(AP - TR))*100.D0
	ENRG(4) = DSQRT(DABS(AP - TR) / TR)*100.D0

	ENRG(1) = TR; ENRG(2) = AP

!   OPEN(13, FILE=filename, STATUS='unknown', POSITION='append')
! 	WRITE(13,*) ''
! 	WRITE(13,*) 'EXACT ENERGY : ', ENRG(1)
! 	WRITE(13,*) 'APPR. ENERGY : ', ENRG(2)
! 	WRITE(13,*) 'ABS.  ENERGY NORM ERROR (%) : ', ENRG(3)
! 	WRITE(13,*) 'REL.  ENERGY NORM ERROR (%) : ', ENRG(4)
! 	CLOSE(13)
! 

	WRITE(*,*) ""
	WRITE(*,*) 'APPR. ENERGY : ', AP
	PRINT*, ""
	PRINT*, 'REL.  ENERGY NORM ERROR (%) : ', ENRG(4)
	PRINT*, ""


END SUBROUTINE ENRGY

SUBROUTINE L2_NORM_ERROR(L2_NORM, COEFF_SOL)
	
	REAL*8, INTENT(OUT) :: L2_NORM(2)
	REAL*8, INTENT(IN) :: COEFF_SOL(DOF)
	INTEGER :: I, J, II, JJ, SUB_NUMLN(2), N, KK, K, PATCH
	REAL*8 :: SUB_LN_LENGTH(2), WEIGHT, DIFF_NN, DET_M
	TYPE(VEC2D) :: SUB_LN
	TYPE(POINT2D) :: GSPT, DMY_GSPT, PHY_PT
	TYPE(MATRIX_22) :: JACOB
	TYPE(DISPLACEMENT) :: EXDISP, APDISP, INTF(2)
	TYPE(STRESS) :: APSTRESS
	
	N = NUMGSPT
	SUB_NUMLN = (/1,1/)
	INTF(:) = DISPLACEMENT(0.D0,0.D0,0.D0)
	
	DO JJ = 1, NUMIR(1)-1
		SUB_LN_LENGTH(1) = DABS(IR_GRID(1,JJ+1) - IR_GRID(1,JJ))/(1.D0*SUB_NUMLN(1))
		DO J=1, SUB_NUMLN(1)
			CALL GAULEG(IR_GRID(1,JJ)+(J-1)*SUB_LN_LENGTH(1), IR_GRID(1,JJ)+J*SUB_LN_LENGTH(1), GSX, GSXW, N)
			DO II = 1, NUMIR(2)-1
				SUB_LN_LENGTH(2) = DABS(IR_GRID(2,II+1) - IR_GRID(2,II))/(1.D0*SUB_NUMLN(2))
				DO I = 1, SUB_NUMLN(2)
					CALL GAULEG(IR_GRID(2,II)+(I-1)*SUB_LN_LENGTH(2), IR_GRID(2,II)+I*SUB_LN_LENGTH(2), GSY, GSYW, N)
					DO KK = 1, N
						DO K = 1, N
							GSPT = POINT2D(GSX(KK),GSY(K))
							DMY_GSPT = GSPT
							JACOB= GET_JACOBIAN_MATRIX(DMY_GSPT)
							DET_M = .DETERMINANT.JACOB
! 							PHY_PT = GET_PHY_PT(DMY_GSPT)
							EXDISP%PIX = EX_DISP(DMY_GSPT)
							CALL APPROXSOL(APDISP%PIX, DMY_GSPT, COEFF_SOL)
							INTF(1)%PIX = INTF(1)%PIX + (EXDISP%PIX - APDISP%PIX)**2*DET_M*GSXW(KK)*GSYW(K)
! 							INTF(1)%PIY = INTF(1)%PIY + (EXDISP%PIY - APDISP%PIY)**2*DET_M*GSXW(KK)*GSYW(K)
							INTF(2)%PIX = INTF(2)%PIX + EXDISP%PIX**2*DET_M*GSXW(KK)*GSYW(K)
! 							INTF(2)%PIY = INTF(2)%PIY + EXDISP%PIY**2*DET_M*GSXW(KK)*GSYW(K)
						ENDDO
					ENDDO
				ENDDO
			ENDDO
		ENDDO
	ENDDO
	
	L2_NORM(1) = DSQRT(INTF(1)%PIX/INTF(2)%PIX)*100.D0
! 	L2_NORM(2) = DSQRT(INTF(1)%PIY/INTF(2)%PIY)*100.D0
	
! 	OPEN(13, FILE=filename, STATUS='unknown', POSITION='append')
! 	WRITE(13,*) ''
! 	WRITE(13,*) 'RELATIVE ERROR IN L2 NORM (%) OF DISPLACEMENT ALONG X-DIRECTION : ', L2_NORM(1)
! 	WRITE(13,*) 'RELATIVE ERROR IN L2 NORM (%) OF DISPLACEMENT ALONG Y-DIRECTION : ', L2_NORM(2)
! 	CLOSE(1)
	
END SUBROUTINE L2_NORM_ERROR

SUBROUTINE APPROXSOL(APDISP, REF_PT, COEFF_SOL)

	TYPE(POINT2D), INTENT(IN) :: REF_PT
	REAL*8, INTENT(IN) :: COEFF_SOL(DOF)
	REAL*8, INTENT(OUT) :: APDISP
	TYPE(FVALUE) :: SF(DOF), LOCAL_SF((BASIS_KVEC(1)%POLY_ORDER+1)*(BASIS_KVEC(2)%POLY_ORDER+1))
	INTEGER :: NUM_XCORDS, NUM_YCORDS, VEC_N, I, J, II, JJ, KK, GLOBAL_INDX
	TYPE(INT2D) :: INDX((BASIS_KVEC(1)%POLY_ORDER+1)*(BASIS_KVEC(2)%POLY_ORDER+1))
	TYPE(POINT2D) :: PHY_PT

	APDISP = 0.D0
	SF(:) = FVALUE(0.D0,0.D0,0.D0,0.D0,0.D0,0.D0)

	CALL GET_ALL_PHY_NURBS_SURFACE_2D(LOCAL_SF(:), INDX, REF_PT)

	DO JJ=1, UBOUND(LOCAL_SF,1)
		DO II=1, DOF
			IF (NDX(1,II).EQ.INDX(JJ)%A .AND. NDX(2,II).EQ.INDX(JJ)%B) THEN
				GLOBAL_INDX = II
				GOTO 111
			ENDIF
		ENDDO
		111 CONTINUE
		SF(GLOBAL_INDX) = LOCAL_SF(JJ)
	ENDDO
	APDISP = DOT_PRODUCT(COEFF_SOL(1:DOF), SF(:)%D00)
	
END SUBROUTINE APPROXSOL

END MODULE ERRORESTIMATE
