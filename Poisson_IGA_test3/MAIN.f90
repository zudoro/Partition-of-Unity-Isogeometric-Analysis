PROGRAM MAIN

  USE GEOMETRY
  USE KNOT_HANDLING
  USE PLOT
  USE BOUNDARY
  USE ERRORESTIMATE
	USE LOADFUNCTION
	
  implicit integer (i-n)
	implicit real(8) (a-h,o-z)

	REAL*8, ALLOCATABLE :: STIF_K(:,:), DMY_K(:,:), STIF_F(:)
  REAL*8 :: ERROR, D, ENRG(4), MAX_NORM(5), L2_NORM(2)
  REAL :: CPTS, CPTE
 	INTEGER :: I, J, II, JJ
 	INTEGER, ALLOCATABLE :: INDX(:)

	PROBLEM = 2

	write(char_problem,fmt='(i1)') PROBLEM
	FILENAME = trim('./data/test-knot-removal') // char_problem
! 	FILENAME = trim('./data/test_pro') // char_problem

	OPEN(131, FILE=FILENAME, STATUS='UNKNOWN')
	DO I = 3, 3
		BS_ORDER = (/I,I/)
		EXTRA_KNOTS = (/I,I/)
		NURBS_REGULARITY = (/bs_order(1)-1, bs_order(2)-1/)

		CALL CPU_TIME(CPTS)

		CALL GET_GEO()

! 		STOP

		ALLOCATE(STIF_K(DOF,DOF), DMY_K(DOF,DOF), STIF_F(DOF), INDX(DOF), CARRAY(DOF,2))

		CALL PLOTGM()

! 		STOP

		CALL GEN_KF(STIF_K, STIF_F)

		DMY_K = STIF_K

! 		OPEN(1, FILE = './data/f1')
! 		DO I=1, DOF
! 			WRITE(1,*) STIF_F(I)
! 		ENDDO
! 		CLOSE(1)
! 	! 
! 		OPEN(1, FILE = './data/k1')
! 		DO I=1, DOF
! 			WRITE(1,*) (STIF_K(I,J), J=1,DOF)
! 		ENDDO
! 		CLOSE(1)
		
		CALL IMPOSEBD(STIF_K,STIF_F)

	! 	OPEN(1, FILE = './data/f2')
	! 	DO I=1, DOF
	! 		WRITE(1,*) F(I)
	! 	ENDDO
	! 	CLOSE(1)
	! 	
	! 	OPEN(1, FILE = './data/k2')
	! 	DO I=1, DOF
	! 		WRITE(1,*) (STIF_K(I,J), J=1,DOF)
	! 	ENDDO
	! 	CLOSE(1)

	! 	RES_K = STIF_K
	! 	DMY_F = F
	! 	
	! 	OPEN(1, FILE = './data/f2')
	! 	DO I=1, DOF
	! 		WRITE(1,*) F(I)
	! 	ENDDO
	! 	CLOSE(1)
		
		CALL LUDCMP(STIF_K, DOF, DOF, INDX, D)
		CALL LUBKSB(STIF_K, DOF, DOF, INDX, STIF_F)
		
	! 	CALL RESIDUALNORM(RES_K, F, DMY_F)

	! 	OPEN(1, FILE = './data/sol')
	! 	DO I=1, DOF
	! 		WRITE(1,*) F(I)
	! 	ENDDO
	! 	CLOSE(1)
		
		CALL MAXNORM(MAX_NORM, STIF_F)
		CALL ENRGY(ENRG, DMY_K, STIF_F)
		CALL L2_NORM_ERROR(L2_NORM, STIF_F)

		PRINT*, BS_ORDER(1), DOF, MAX_NORM(1), L2_NORM(1), ENRG(4)
		WRITE(131, *) DOF, BS_ORDER(1), MAX_NORM(1), L2_NORM(1), ENRG(4), ENRG(2)

		CALL CPU_TIME(CPTE)

		PRINT*, ""
		PRINT*, "ELAPSED CPU TIME : ", CPTE - CPTS
		PRINT*, ""

		DEALLOCATE(STIF_K, DMY_K, STIF_F, INDX, BD_MASK, NDX, CARRAY)

	ENDDO
	CLOSE(131)
END PROGRAM MAIN
