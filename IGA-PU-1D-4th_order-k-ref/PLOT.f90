MODULE PLOT

	USE GEOMETRY
	USE NURBS_BASIS

	implicit integer (i-n)
	implicit real(8) (a-h,o-z)

CONTAINS

SUBROUTINE PLOTGM()
	
	INTEGER :: PATCH, PATCH_ROW, PATCH_COLUMN, I, J, K, II, JJ, KK
	REAL*8 :: PUSUM, BSSUM, parpt
	TYPE(FUNCTION_1D) :: PUTEST, BSTEST, REFVAL, PHYPU
	CHARACTER(LEN=2) :: CHAR_PATCH, CHAR_LOCBS
	
	OPEN(1, FILE = './data/ndx')
	DO I=1, DOF
		WRITE(1,*) NDX(:, I)
	ENDDO
	CLOSE(1)
	
	OPEN(1, FILE = './data/patchbdpt')
	DO I = 1, UBOUND(PATCHBDPT,1)
		WRITE(1,*) PATCHBDPT(I,:)
	ENDDO
	CLOSE(1)
	
	open(1, file='./data/ir_grid')
	DO PATCH_ROW = 1, NUMPATCH
		DO PATCH_COLUMN = 1, NUMPATCH
			DO I = 1, NUMIR(PATCH_ROW, PATCH_COLUMN)
				WRITE(1, *) PATCH_ROW, PATCH_COLUMN, IR_GRID(PATCH_ROW, PATCH_COLUMN, I)
			enddo
		enddo
	enddo
	
	DO PATCH = 1, NUMPATCH
		write(CHAR_PATCH,fmt='(i2)') PATCH
		FILENAME = trim('./data/bspline') // CHAR_PATCH
		OPEN(1, FILE=FILENAME, STATUS='UNKNOWN')
		DO K = 1, NUMBS(PATCH)
			DO I = 1, 1001
				parpt = 0.0010D0*(I-1)
				BSTEST = GET_PHY_BSPLINE1D(K-1, parpt, PATCH)
				WRITE(1, *) parpt, BSTEST%VAL(0), BSTEST%VAL(1), BSTEST%VAL(2)
			ENDDO
			WRITE(1,*)
		ENDDO
		CLOSE(1)
	ENDDO

	DO PATCH = 1, NUMPATCH
		write(CHAR_PATCH,fmt='(i2)') PATCH
		FILENAME = trim('./data/pu') // CHAR_PATCH
		OPEN(1, FILE=FILENAME, STATUS='UNKNOWN')
		DO I = 1, 1001
			parpt = 0.0010D0*(I-1)
			IF (PATCH==1) THEN 
              PHYPU = MAP_F(PARPT)
            ELSEIF (PATCH==2) THEN 
              PHYPU = MAP_G(PARPT)
            ENDIF
			PUTEST = GET_PHYPU1D(PARPT, PATCH)
			WRITE(1, *) PHYPU%VAL(0), PUTEST%VAL(0), PUTEST%VAL(1), PUTEST%VAL(2)
		ENDDO
		CLOSE(1)
	ENDDO
	
! 	OPEN(1, FILE = './data/subirndx')
! 	DO I = 1, NUMPATCH
! 		DO J = 1, NUMPATCH
! 			WRITE(1,*) SUBIRNDX(I,J,:)
! 		ENDDO
! 	ENDDO
! 	CLOSE(1)

! 	OPEN(1, FILE = './data/pu_knot')
! 	WRITE(1,*) 'POLY_ORDER : ', PU_KVEC%POLY_ORDER
! 	WRITE(1,*) 'LENGTH : ', PU_KVEC%LENGTH
! 	DO I=0, PU_KVEC%LENGTH
! 		WRITE(1,102) PU_KVEC%KNOTS(I)
! 	ENDDO
! 	CLOSE(1)

! 	OPEN(1, FILE = './data/col_pt')
! 	DO I = 1, UBOUND(COL_PT,1)
! 		WRITE(1,*) COL_PT(I), 0.0D0
! 	ENDDO
! 	CLOSE(1)

! 	OPEN(1, FILE = './data/jeta')
! 	DO I = 1, UBOUND(JETA,1)
! 		WRITE(1,*) JETA(I), 0.0d0
! 	ENDDO
! 	CLOSE(1)
	
	
	
! 	OPEN(1, FILE = 'data/gspt')
! 	DO I = 1, UBOUND(JETA,1)-1
! 		DO J = 1, NUMGSPT
! 			WRITE(1,*) GSPT(I, J), 0.0d0
! 		ENDDO
! 	ENDDO
	
! 	DO PATCH = 1, NUMPATCH
! 		write(CHAR_PATCH,fmt='(i2)') PATCH
! 		FILENAME = trim('./data/pu') // CHAR_PATCH
! 		OPEN(111, FILE=FILENAME, STATUS='UNKNOWN')
! 		DO I = 1, 101
! 			PUTEST = GET_BSPLINEPU1D(0.010D0*(I-1), PATCH)
! 			WRITE(111,*) 0.010D0*(I-1), PUTEST%VAL(0:2)
! 		ENDDO
! 		CLOSE(111)
! 	ENDDO
	
! 	DO I = 0, 1000
! 		PUSUM = 0.0D0
! 		BSSUM = 0.0D0
! 		DO PATCH = 1, NUMPATCH
! 			PUTEST = GET_BSPLINEPU1D(0.0010D0*I, PATCH)
! 			DO J = 0, IORDER
! 				BSTEST = GET_LOCBASIS(0.0010D0*I, J, PATCH)
! 				BSSUM = BSSUM + PUTEST%VAL(0)*BSTEST%VAL(0)
! 				IF (DABS(0.0010D0*I - 0.20D0)<=EPS) THEN
! 					PRINT*, PATCH, J, BSTEST
! 				ENDIF
! 			ENDDO
! 			PUSUM = PUSUM + PUTEST%VAL(0)
! 		ENDDO
! 		IF (DABS(PUSUM - 1.0D0)>EPS) THEN
! 			PRINT*, 'FAIL PU TEST', PUSUM
! 			STOP
! 		ENDIF
! 		IF (DABS(BSSUM - 1.0D0)>EPS) THEN
! 			PRINT*, 'FAIL BS TEST', BSSUM, 0.0010D0*I, PATCH
! 			STOP
! 		ENDIF
! 	ENDDO
	
! 	DO I = 0, 1000
! 		BSTEST = GET_LOCBASIS(0.010D0*(I-1), J, PATCH)
! 	ENDDO

! 	DO PATCH = 1, NUMPATCH
! 		DO I = 1,101
! 			REFVAL = PATCH_TO_REF1D(X, PATCH)
	
! 	DO PATCH = 1, NUMPATCH
! 	DO J = 0, IORDER
! 		write(CHAR_LOCBS,fmt='(i2)') J
! 		write(CHAR_PATCH, fmt='(i1)') PATCH
! 		FILENAME = trim('./data/bs') // CHAR_PATCH // CHAR_LOCBS
! 		OPEN(111, FILE=FILENAME, STATUS='UNKNOWN')
! 		DO I = 1, 101
! 			BSTEST = GET_LOCBASIS(0.010D0*(I-1), J, PATCH)
! 			WRITE(111,*) 0.010D0*(I-1), BSTEST%VAL(0:2)
! 		ENDDO
! 		CLOSE(111)
! 	ENDDO
! 	ENDDO

	102 FORMAT(1000(f20.8,1x))
END SUBROUTINE PLOTGM

END MODULE PLOT
