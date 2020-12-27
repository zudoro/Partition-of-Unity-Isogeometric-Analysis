MODULE ERRORESTIMATE

	USE GEOMETRY
	USE NURBS_BASIS
	USE LOADFUNCTION
	USE GSQUAD

    IMPLICIT NONE
    
CONTAINS

!----MAX. NORM ESTIMATE----
SUBROUTINE MAXNORM(MAX_NORM, COEFF_SOL)

    REAL*8, INTENT(OUT) :: MAX_NORM
	REAL*8, INTENT(IN) :: COEFF_SOL(DOF)
	
	TYPE(POINT2D) :: MESH, PHY_PT, REF_PT, PAR_PT
	TYPE(INT2D) :: MAXGRID
	INTEGER :: II, JJ, I, J, PATCH_ROW, PATCH_COLUMN, NUMGSPT, N
	REAL*8 :: R, THETA, ERR, TMP_ERR, MAXDISP, APDISP, EXDISP
	INTEGER :: PATCH
	
	MAX_NORM = 0.0D0
	
	ERR = 0.0D0
	TMP_ERR = 0.0D0
	MAXDISP = 0.0D0
	
! 	print*, 'here-1'
	
	OPEN(11, FILE = './data/ext_disp')
	OPEN(21, FILE = './data/app_disp')
	OPEN(31, FILE = './data/err_disp')
	
	DO PATCH = 1, NUMPATCH
        DO I = 1, 99
            DO J = 0, 100
                PAR_PT = POINT2D(0.010D0*I, 0.010D0*J)
                PHY_PT = GET_PHY_PT(PAR_PT, PATCH)
                APDISP = APSOL(PAR_PT, COEFF_SOL, PATCH)
                EXDISP = EX_DISP(PAR_PT, PATCH)
                TMP_ERR = DABS(EXDISP - APDISP)
                
                WRITE(11,*) PHY_PT, EXDISP
                WRITE(21,*) PHY_PT, APDISP
                WRITE(31,*) PHY_PT, TMP_ERR
                
                IF (TMP_ERR.GE.ERR) THEN
                    ERR = TMP_ERR
                ENDIF
                IF (DABS(EXDISP).GE.MAXDISP) THEN
                    MAXDISP = DABS(EXDISP)
                ENDIF
                123 CONTINUE
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

	PRINT*, 'MAXDISP= ', MAXDISP
	PRINT*, 'ERR= ', ERR
	
	MAX_NORM = ERR*100.D0/MAXDISP
	
END SUBROUTINE MAXNORM

SUBROUTINE ENRGY(ENRG, DMY_K, SOL)

	REAL*8, INTENT(IN) :: DMY_K(DOF,DOF), SOL(DOF)
	REAL*8, INTENT(OUT) :: ENRG(4)
	REAL*8 :: AP, TR
	INTEGER :: I, J

	AP = 0.50D0*DOT_PRODUCT(SOL, MATMUL(DMY_K, SOL))
	
	IF (PDE=='ELPT') THEN
		IF (PROBLEM==0) THEN
			TR = 0.196350D0
		ELSEIF (PROBLEM.EQ.1) THEN
			TR = PI/4.D0
		ELSEIF (PROBLEM.EQ.2) THEN
			TR = 0.785398163397448280D0
		ELSEIF (PROBLEM.EQ.3) THEN
			TR = 2.0943951023931948D0
! 		ELSEIF (PROBLEM==4) THEN
! 			TR = 0.91811333093758096258D0
		ENDIF
	ELSEIF (PDE=='CNVD') THEN
		IF (PROBLEM==1) THEN
			TR = 0.50D0*0.03926990D0
		ELSEIF (PROBLEM.EQ.2) THEN
			TR = 0.50d0*0.001036290d0
		ELSEIF (PROBLEM.EQ.3) THEN
			TR = 0.50d0*0.06509790d0
		ELSEIF (PROBLEM==4) THEN
			IF (DABS(EPSLN - 0.10D0)<=EPS) THEN
				TR = 0.50D0*0.154080d0
			ELSEIF (DABS(EPSLN - 0.0010D0)<=EPS) THEN
				TR = 0.7421950D0
			ELSEIF (DABS(EPSLN - 0.00010D0)<=EPS) THEN
				TR = 0.7391550D0
			ELSEIF (DABS(EPSLN - 0.000010D0)<=EPS) THEN
				TR = 0.7388510D0
			ENDIF
		ENDIF
    ELSEIF (PDE=='FOUR') THEN 
        IF (PROBLEM==1) THEN
            TR = 32768.0D0/1225.0D0
        ELSEIF (PROBLEM==2) THEN
            TR = 156237824.0d0/5780775.0d0
        ELSEIF (PROBLEM==3) THEN
!             TR = (243.0d0*(739860151.0d0-27372582.0d0*dexp(4.0d0)+253151.0d0*dexp(8.0d0)))/(32.0d0*dexp(4.0d0))
            TR = 42.7749926294929279270D0
        ELSEIF (PROBLEM==4) THEN
            TR = 9.359328113819591620D0
        ELSEIF (PROBLEM==5) THEN 
            TR = 16.75516081914556210D0
        ENDIF
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

! SUBROUTINE L2_NORM_ERROR(L2_NORM, COEFF_SOL)
! 	
! 	REAL*8, INTENT(OUT) :: L2_NORM(2)
! 	REAL*8, INTENT(IN) :: COEFF_SOL(DOF)
! 	
! 	INTEGER :: I, J, II, JJ, N, KK, K, PATCH, PATCH_COLUMN, PATCH_ROW
! 	REAL*8 :: WEIGHT, DIFF_NN, DET_M
! 	REAL*8, ALLOCATABLE :: L2_GSX(:), L2_GSXW(:), L2_GSY(:), L2_GSYW(:)
! 	TYPE(POINT2D) :: GSPT, PHY_PT
! 	TYPE(MATRIX_22) :: JACOB
! 	REAL*8 :: EXDISP, APDISP, INTF(2)
! 	TYPE(RECPATCH) :: IRBOX
! 	
! 		N = 2*(PUORDER + SUM(BS_ORDER(:,1)) + SUM(BS_ORDER(:,2)))
! 	
! 	ALLOCATE(L2_GSX(N), L2_GSY(N), L2_GSXW(N), L2_GSYW(N))
! 	
! 	PATCH_ROW = 1
! 	
! 	INTF(:) = DISPLACEMENT(0.D0,0.D0,0.D0)
! 	
! ! 	DO PATCH_COLUMN = 1, NUMPATCH(2)
! 		DO II = 1, L2_NUMIR(1) - 1
! 			DO JJ = 1, L2_NUMIR(2) - 1
! 				IRBOX%PT1 = POINT2D(L2_GRID(1,II), L2_GRID(2,JJ))
! 				IRBOX%PT2 = POINT2D(L2_GRID(1,II+1), L2_GRID(2,JJ))
! 				IRBOX%PT3 = POINT2D(L2_GRID(1,II+1), L2_GRID(2,JJ+1))
! 				IRBOX%PT4 = POINT2D(L2_GRID(1,II), L2_GRID(2,JJ+1))
! 				!------------------! VECTORIZE GAUSS POINTS AND WEIGHTS------------------------------
! 				CALL GAULEG(IRBOX%PT1%X, IRBOX%PT2%X, L2_GSX, L2_GSXW, N)
! 				CALL GAULEG(IRBOX%PT1%Y, IRBOX%PT4%Y, L2_GSY, L2_GSYW, N)
! 				DO KK = 1, N
! 					DO K = 1, N
! 						GSPT = POINT2D(L2_GSX(KK),L2_GSY(K))
! 						JACOB= GET_JACOBIAN_MATRIX(GSPT)
! 						DET_M = .DETERMINANT.JACOB
! ! 							PHY_PT = GET_PHY_PT(DMY_GSPT)
! 						EXDISP%PIX = EX_DISP(GSPT)
! 						CALL APPROXSOL(APDISP%PIX, GSPT, COEFF_SOL)
! 						INTF(1)%PIX = INTF(1)%PIX + (EXDISP%PIX - APDISP%PIX)**2*DET_M*L2_GSXW(KK)*L2_GSYW(K)
! 						INTF(2)%PIX = INTF(2)%PIX + EXDISP%PIX**2*DET_M*L2_GSXW(KK)*L2_GSYW(K)
! 					ENDDO
! 				ENDDO
! 			ENDDO
! 		ENDDO
! ! 	ENDDO
! 	
! 	L2_NORM(1) = DSQRT(INTF(1)%PIX/INTF(2)%PIX)*100.D0
! ! 	L2_NORM(2) = DSQRT(INTF(1)%PIY/INTF(2)%PIY)*100.D0
! 	
! ! 	OPEN(13, FILE=filename, STATUS='unknown', POSITION='append')
! ! 	WRITE(13,*) ''
! ! 	WRITE(13,*) 'RELATIVE ERROR IN L2 NORM (%) OF DISPLACEMENT ALONG X-DIRECTION : ', L2_NORM(1)
! ! 	WRITE(13,*) 'RELATIVE ERROR IN L2 NORM (%) OF DISPLACEMENT ALONG Y-DIRECTION : ', L2_NORM(2)
! ! 	CLOSE(1)
! 	
! END SUBROUTINE L2_NORM_ERROR

! SUBROUTINE RESIDUALNORM(RES_K, SOL, F)
!     REAL*8, INTENT(IN) :: RES_K(DOF, DOF), SOL(DOF), F(DOF)
!     
!     REAL*8 :: APP_F(DOF), ERR
!     INTEGER :: I, J, K
!     
!     ERR = 0.0D0 
!     
!     APP_F = MATMUL(RES_K, SOL)
!     
!     DO I = 1, DOF
!         IF (DABS(APP_F(I) - F(I))>ERR) THEN 
!             ERR = DABS(APP_F(I) - F(I))
!         ENDIF
!     ENDDO
!     
!     PRINT*, ERR
!     
! END SUBROUTINE RESIDUALNORM

REAL*8 FUNCTION APSOL(REFPT, COEFF_SOL, PATCH)

	TYPE(POINT2D), INTENT(IN) :: REFPT
	REAL*8, INTENT(IN) :: COEFF_SOL(DOF)
	INTEGER, INTENT(IN) :: PATCH
	
	TYPE(POINT2D) :: PHYPT, PARPT
	TYPE(FVALUE) :: BS
	INTEGER :: I, II, J, JJ
	
! 	TYPE(POINT2D) :: PHY_PT
! 	TYPE(MATRIX_22) :: JACOB
	
	APSOL = 0.0D0
    
    IF (PATCH==1 .AND. REFPT%Y>INV_DELTA(1)) THEN 
      II = 0
      DO J = 1, LOC_NUMBS(PATCH, 2)
        DO I = 1, LOC_NUMBS(PATCH, 1)
          II = II + 1
          BS = GET_BS(REFPT, (/I-1, J-1/), PATCH)
          APSOL = APSOL + COEFF_SOL(II)*BS%D00
        ENDDO
      ENDDO
      !-------------------------------------------------------!
      PHYPT = GET_PHY_PT(REFPT, PATCH)
      PARPT = GET_PAR_PT(PHYPT, 2)
      !-------------------------------------------------------!
      II = LOC_NUMBS(1, 1)*LOC_NUMBS(1, 2)
      DO J = 1, LOC_NUMBS(2, 2)
        DO I = 1, LOC_NUMBS(2, 1)
          II = II + 1
          BS = GET_BS(PARPT, (/I-1, J-1/), 2)
          APSOL = APSOL + COEFF_SOL(II)*BS%D00
        ENDDO
      ENDDO
    ELSEIF (PATCH==2 .AND. REFPT%Y<INV_DELTA(2)) THEN
      II = LOC_NUMBS(1, 1)*LOC_NUMBS(1, 2)
      DO J = 1, LOC_NUMBS(PATCH, 2)
        DO I = 1, LOC_NUMBS(PATCH, 1)
          II = II + 1
          BS = GET_BS(REFPT, (/I-1, J-1/), PATCH)
          APSOL = APSOL + COEFF_SOL(II)*BS%D00
        ENDDO
      ENDDO
      PHYPT = GET_PHY_PT(REFPT, PATCH)
      PARPT = GET_PAR_PT(PHYPT, 1)
      II = 0
      DO J = 1, LOC_NUMBS(1, 2)
        DO I = 1, LOC_NUMBS(1, 1)
          II = II + 1
          BS = GET_BS(PARPT, (/I-1, J-1/), 1)
          APSOL = APSOL + COEFF_SOL(II)*BS%D00
        ENDDO
      ENDDO
    ELSE
      IF (PATCH==1) THEN
        II = 0
      ELSEIF (PATCH==2) THEN 
        II = LOC_NUMBS(1, 1)*LOC_NUMBS(1, 2)
      ENDIF
      DO J = 1, LOC_NUMBS(PATCH, 2)
        DO I = 1, LOC_NUMBS(PATCH, 1)
          II = II + 1
          BS = GET_BS(REFPT, (/I-1, J-1/), PATCH)
          APSOL = APSOL + COEFF_SOL(II)*BS%D00
        ENDDO
      ENDDO
    ENDIF
    
END FUNCTION APSOL

END MODULE ERRORESTIMATE
