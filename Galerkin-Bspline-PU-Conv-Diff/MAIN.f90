PROGRAM MAIN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! --------------------------- [[[[[ Title of Paper ]]]]] ------------------------------
! Partition of unity with flat-top constructed by a geometrical mapping for analysis of 
! elliptic boundary value problems containing corner or boundary layer singularities

!! --------------------------- [[[[[ Example Number ]]]]] ------------------------------
! Example 1 - Convection-diffusion equation with boundary layer singularity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  USE GEOMETRY
  USE PLOT
  USE BOUNDARY
  USE ERRORESTIMATE
	USE LOADFUNCTION

	REAL*8, ALLOCATABLE :: STIF_K(:,:), DMY_K(:,:), STIF_F(:)
  REAL*8 :: ERROR, D, ENRG(4), MAX_NORM(5), L2_NORM(2)
  REAL :: CPTS, CPTE
 	INTEGER :: I, J, II, JJ
 	INTEGER, ALLOCATABLE :: INDX(:)

!-----------------!
	PDE = 'CNVD'
! 	PDE = 'ELPT'
!-----------------!

!-----------------!
	BL = 'Y'
! 	BL = 'N'
!-----------------!

!-----------------!
	USE_ENRICH = 'Y'
! 	USE_ENRICH = 'N'
!-----------------!

!-----------------!
! 	MESHTYPE = 'SH'
	MESHTYPE = 'NA'
!-----------------!

!-----------------!
! 	ZEROBD = 'Y'
	ZEROBD = 'N'
!-----------------!

	PROBLEM = 4
	PUORDER = 3
	
	ALLOCATE(PATCHBDX(NUMPATCH(1),2), PATCHBDY(NUMPATCH(2),2))
! 	
	PATCHBDX(1,1) = 0.0D0
	PATCHBDX(1,2) = 1.0D0
	
	PATCHBDY(1,:) = (/BETA, 1.0D0/)
	PATCHBDY(2,:) = (/ALPHA, BETA/)
	PATCHBDY(3,:) = (/0.0D0, ALPHA/)
	
	write(char_problem,fmt='(i1)') PROBLEM
	FILENAME = trim('./data/test_pro') // char_problem
! 	FILENAME = trim('./data/test_pro') // char_problem

	OPEN(131, FILE=FILENAME, STATUS='UNKNOWN')
	
	WRITE(131,*) 'PDE=',PDE
	IF (PDE=='CNVD') THEN
		WRITE(131,*) 'Epsilon=', EPSLN
	ENDIF
	WRITE(131,*) 'Boundary Layer Problem?  ', BL
	WRITE(131,*) 'Do you use enriching functions?  ', USE_ENRICH
	WRITE(131,*) 'Do you Apply any special type of mesh refinement?  ', MESHTYPE
	WRITE(131,*) 'Do you impose zero boundary condition?  ', ZEROBD
	WRITE(131, *) 'DOF - BD_DOF  ', 'ORDER  ', 'MAX.  ', 'L2  ', 'ENRG'
	DO I = 3, 5, 2
		
		BS_ORDER(1,1,:) = (/I,I/)
		BS_ORDER(1,2,:) = (/I,I/)
		BS_ORDER(1,3,:) = (/I,I/)
		
! 		BS_ORDER(1,1,:) = (/3,3/)
! 		BS_ORDER(1,2,:) = (/3,3/)
! 		BS_ORDER(1,3,:) = (/3,3/)

! 		EXTRA_KNOTS(1,1,:) = (/2,3/)
! 		EXTRA_KNOTS(1,2,:) = (/4,5/)
! 		EXTRA_KNOTS(1,3,:) = (/1,6/)
		
! 		EXTRA_KNOTS(1,1,:) = (/int(0.5*(i-1)),int(0.5*(i-1))/)
! 		EXTRA_KNOTS(1,2,:) = (/int(0.5*(i-1)),int(0.5*(i-1))/)
! 		EXTRA_KNOTS(1,3,:) = (/i,i/)
		
			EXTRA_KNOTS(1,1,:) = (/i,i/)
			EXTRA_KNOTS(1,2,:) = (/i,i/)
			EXTRA_KNOTS(1,3,:) = (/i,i/)
			
! 		EXTRA_KNOTS(1,1,:) = (/0,0/)
! 		EXTRA_KNOTS(1,2,:) = (/0,0/)
! 		EXTRA_KNOTS(1,3,:) = (/0,0/)
		
		NURBS_REGULARITY(1,1,:) = (/BS_ORDER(1,1,1)-1, BS_ORDER(1,1,2)-1/)
		NURBS_REGULARITY(1,2,:) = (/BS_ORDER(1,2,1)-1, BS_ORDER(1,2,2)-1/)
		NURBS_REGULARITY(1,3,:) = (/BS_ORDER(1,3,1)-1, BS_ORDER(1,3,2)-1/)
		
		CALL CPU_TIME(CPTS)

		CALL GET_GEO()

! 		STOP
		
		ALLOCATE(STIF_K(DOF,DOF), DMY_K(DOF,DOF), STIF_F(DOF), INDX(DOF))

		CALL PLOTGM()

! 		STOP

		CALL GEN_KF(STIF_K, STIF_F)

		DMY_K = STIF_K

! 		OPEN(1, FILE = './data/f1')
! 		DO II=1, DOF
! 			WRITE(1,*) STIF_F(II)
! 		ENDDO
! 		CLOSE(1)
! 	! 
! 		OPEN(1, FILE = './data/k1')
! 		DO II=1, DOF
! 			WRITE(1,*) (STIF_K(II,JJ), JJ=1,DOF)
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
! 		
! 		OPEN(1, FILE = './data/k2')
! 		DO II=1, DOF
! 			WRITE(1,*) (STIF_K(II,JJ), JJ=1,DOF)
! 		ENDDO
! 		CLOSE(1)

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

! 		OPEN(1, FILE = './data/sol')
! 		DO II = 1, DOF
! 			WRITE(1,*) STIF_F(II)
! 		ENDDO
! 		CLOSE(1)
		
		CALL MAXNORM(MAX_NORM, STIF_F)
		CALL ENRGY(ENRG, DMY_K, STIF_F)
		CALL L2_NORM_ERROR(L2_NORM, STIF_F)

		PRINT*, BS_ORDER(1,3,1), DOF - BD_DOF, MAX_NORM(1), L2_NORM(1), ENRG(2), ENRG(4)
		WRITE(131, *) DOF - BD_DOF, BS_ORDER(1,:,1), MAX_NORM(1), L2_NORM(1), ENRG(2), ENRG(4)

		CALL CPU_TIME(CPTE)

		PRINT*, ""
		PRINT*, "ELAPSED CPU TIME : ", CPTE - CPTS
		PRINT*, ""

		DEALLOCATE(STIF_K, DMY_K, STIF_F, INDX, NDX)
	ENDDO
	DEALLOCATE(PATCHBDX, PATCHBDY)
	CLOSE(131)
END PROGRAM MAIN
