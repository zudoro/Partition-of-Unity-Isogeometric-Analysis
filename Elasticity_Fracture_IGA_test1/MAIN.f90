PROGRAM MAIN

  USE GEOMETRY
  USE PLOT
  USE BOUNDARY
  USE ERRORESTIMATE
	
  implicit integer (i-n)
	implicit real(8) (a-h,o-z)

	REAL*8, ALLOCATABLE :: STIF_K(:,:), DMY_K(:,:), STIF_F(:)
  REAL*8 :: ERROR, D, ENRG(4), MAX_NORM(5), L2_NORM(2)
  REAL :: CPTS, CPTE
 	INTEGER :: I, J, II, JJ
 	INTEGER, ALLOCATABLE :: INDX(:)
	TYPE(POINT2D) :: TSPT, PHYPT

	PROBLEM = 5
	BCTYPE = 'diric'

	PRINT*, '---------------------------------------------------'
	PRINT*, 'Elasticity on L-shaped domain by using Isogeometric Analysis with one patch'
	PRINT*, 'Problem : ', PROBLEM
	PRINT*, 'Boundary Condition : ', BCTYPE
	PRINT*, 'Plane type : ', PLANE_TYPE
	PRINT*, 'Lambda : ', LAMBDA
	PRINT*, 'Q : ', MU
	PRINT*, '---------------------------------------------------'
	
	write(char_problem,fmt='(i1)') PROBLEM
	
	FILENAME = trim('./data/Elasticity_L_shaped_k_ref_IGA_pro') // char_problem // trim('_') // BCTYPE // trim('_') // PLANE_TYPE

! 	OPEN(131, FILE=FILENAME, STATUS='UNKNOWN')
	OPEN(131, FILE='test')

	DO I = 2, 2
		EXTRA_KNOTS = (/I,I/)
		IGA_ORDER = (/I,I/)
		NURBS_REGULARITY = (/IGA_ORDER(1)-1, IGA_ORDER(2)-1/)
		NUMGSPT = 18
		
		CALL CPU_TIME(CPTS)

		CALL GET_GEO()

		ALLOCATE(STIF_K(2*DOF,2*DOF), DMY_K(2*DOF,2*DOF), STIF_F(2*DOF), INDX(2*DOF), CARRAY(DOF,2))
		ALLOCATE(GSX(NUMGSPT), GSXW(NUMGSPT), GSY(NUMGSPT), GSYW(NUMGSPT))

! 		PRINT*, 'MAIN - HERE0'

		CALL PLOTGM()
! 		STOP
! 		PRINT*, 'MAIN - HERE1'

		CALL GEN_KF(STIF_K, STIF_F)

! 		CALL STIF_BUGKILLER(STIF_K, STIF_F)

		DMY_K = STIF_K

! 		OPEN(1, FILE = './data/f1')
! 		DO II=1, 2*DOF
! 			WRITE(1,*) STIF_F(II)
! 		ENDDO
! 		CLOSE(1)
! 	! 
! 		OPEN(1, FILE = './data/k1')
! 		DO II=1, 2*DOF
! 			WRITE(1,*) (STIF_K(II,JJ), JJ=1,2*DOF)
! 		ENDDO
! 		CLOSE(1)
		
		CALL IMPOSEBD(STIF_K,STIF_F)

	! 	OPEN(1, FILE = './data/f2')
	! 	DO I=1, 2*DOF
	! 		WRITE(1,*) F(I)
	! 	ENDDO
	! 	CLOSE(1)
	! 	
	! 	OPEN(1, FILE = './data/k2')
	! 	DO I=1, 2*DOF
	! 		WRITE(1,*) (STIF_K(I,J), J=1, 2*DOF)
	! 	ENDDO
	! 	CLOSE(1)

	! 	RES_K = STIF_K
	! 	DMY_F = F
	! 	
	! 	OPEN(1, FILE = './data/f2')
	! 	DO I=1, 2*DOF
	! 		WRITE(1,*) F(I)
	! 	ENDDO
	! 	CLOSE(1)
		
		CALL LUDCMP(STIF_K, 2*DOF, 2*DOF, INDX, D)
		CALL LUBKSB(STIF_K, 2*DOF, 2*DOF, INDX, STIF_F)
		
	! 	CALL RESIDUALNORM(RES_K, F, DMY_F)

	! 	OPEN(1, FILE = './data/sol')
	! 	DO I=1, 2*DOF
	! 		WRITE(1,*) F(I)
	! 	ENDDO
	! 	CLOSE(1)
		
		CALL MAXNORM(MAX_NORM, STIF_F)
		CALL ENRGY(ENRG, DMY_K, STIF_F)
		CALL L2_NORM_ERROR(L2_NORM, STIF_F)

		PRINT*, IGA_ORDER(1), EXTRA_KNOTS(1), 2*DOF-2*BD_DOF, MAX_NORM(:), L2_NORM(:), ENRG(4), ENRG(2)
		WRITE(131, *) IGA_ORDER(1), EXTRA_KNOTS(1), 2*DOF-2*BD_DOF, MAX_NORM(:), L2_NORM(:), ENRG(4), ENRG(2)
		
! 		PRINT*, BS_ORDER(1), DOF, ENRG(4)
! 		WRITE(131, *) DOF, BS_ORDER(1), ENRG(4), ENRG(2)
		
		
		CALL CPU_TIME(CPTE)

		PRINT*, ""
		PRINT*, "ELAPSED CPU TIME : ", CPTE - CPTS
		PRINT*, ""

		DEALLOCATE(STIF_K, DMY_K, STIF_F, INDX, BD_MASK, NDX, DMY_NDX, INTERFACE_NDX, CARRAY)
		DEALLOCATE(GSX, GSY, GSXW, GSYW)

	ENDDO
	CLOSE(131)
END PROGRAM MAIN
