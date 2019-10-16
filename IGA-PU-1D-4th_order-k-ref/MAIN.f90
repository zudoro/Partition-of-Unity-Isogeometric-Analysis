PROGRAM MAIN

	USE GEOMETRY
	USE PLOT
	USE BOUNDARY
	USE ERRORESTIMATE
	USE ASSEMBLE
	USE LUDECOMPOSITION
	
  implicit integer (i-n)
	implicit real(8) (a-h,o-z)

	REAL*8, ALLOCATABLE :: STIF_K(:,:), DMY_K(:,:), STIF_F(:)
	REAL*8 :: ERROR, D, ENRG(4), MAX_NORM, L2_NORM
	REAL :: CPTS, CPTE
 	INTEGER :: I, J, II, JJ
 	INTEGER, ALLOCATABLE :: INDX(:)
 	
 	PUORDER = 3
	
	write(char_problem,fmt='(i1)') PROBLEM
	write(char_puorder,fmt='(i1)') PUORDER
	write(char_iorder,fmt='(i2)') IORDER
	write(char_numpatch,fmt='(i2)') NUMPATCH

! 	FILENAME = trim('./data/test') // char_problem // trim('PU') // char_puorder // trim('NP') // char_numpatch
	FILENAME = trim('./data/test') // char_problem
	PUNAME = trim('./data/pu')
	
	OPEN(131, FILE=FILENAME, STATUS='UNKNOWN')
! 	OPEN(132, FILE=PUNAME, STATUS='UNKNOWN')
	
	DO I = 3, 10
 
		EXTRA_KNOTS = (/I, I+10/)
		BSORDER = (/I, I/)
! 		NURBS_REGULARITY = (/BS_ORDER(1)-1, BS_ORDER(2)-1/)

		CALL CPU_TIME(CPTS) 
		
		CALL GET_GEO()
		
! 		STOP

		ALLOCATE(STIF_K(DOF,DOF), DMY_K(DOF, DOF), STIF_F(DOF), INDX(DOF), CARRAY(DOF,2))

		CALL PLOTGM()

! 		STOP

		CALL GEN_KF(STIF_K, STIF_F)

		DMY_K = STIF_K

! 		OPEN(1, FILE = './data/f1')
! 		DO II=1, DOF
! 			WRITE(1,*) STIF_F(II)
! 		ENDDO
! 		CLOSE(1)
	! 
! 		OPEN(1, FILE = './data/k1')
! 		DO II=1, DOF
! 			WRITE(132,*) (STIF_K(II,JJ), JJ=1,DOF)
! 		ENDDO
! 		CLOSE(1)
! 		
		CALL IMPOSEBD(STIF_K,STIF_F)
! 
! 		OPEN(1, FILE = './data/f2')
! 		DO II=1, DOF
! 			WRITE(1,*) STIF_F(II)
! 		ENDDO
! 		CLOSE(1)
		
! 		OPEN(1, FILE = './data/k2')
! 		DO II=1, DOF
! 			WRITE(132,*) (STIF_K(II,J), J=1,DOF)
! 		ENDDO
! 		WRITE(132,*) ("", J=1,DOF)
! 		CLOSE(1)

! 		RES_K = STIF_K
! 		DMY_F = F
		
		CALL LUDCMP(STIF_K, DOF, DOF, INDX, D)
		CALL LUBKSB(STIF_K, DOF, DOF, INDX, STIF_F)
		
	! 	CALL RESIDUALNORM(RES_K, F, DMY_F)

! 		OPEN(1, FILE = './data/sol')
! 		DO II = 1, DOF
! 			WRITE(1,*) STIF_F(II)
! 		ENDDO
! 		CLOSE(1)
		
		MAX_NORM = MAXNORM(STIF_F)
		CALL ENRGY(ENRG, DMY_K, STIF_F)
! 		CALL L2_NORM_ERROR(L2_NORM, STIF_F)
		
		CALL CPU_TIME(CPTE)

		PRINT*, PUORDER, BSORDER, NUMPATCH, DOF-4, MAX_NORM, ENRG(4), ENRG(2), ENRG(1)
		WRITE(131, *) PUORDER, BSORDER, NUMPATCH, DOF-4, MAX_NORM, ENRG(4), ENRG(2), ENRG(1)

		PRINT*, ""
		PRINT*, "ELAPSED CPU TIME : ", CPTE - CPTS
		PRINT*, ""

		DEALLOCATE(STIF_K, DMY_K, STIF_F, NDX, INDX, CARRAY, GSPT, GSW)
	ENDDO

	CLOSE(131)
! 	CLOSE(132)
END PROGRAM MAIN
