PROGRAM MAIN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! --------------------------- [[[[[ Title of Paper ]]]]] ------------------------------
! Partition of Unity isogeometric analysis (PU-IGA) for numerical solutions of fourth order partial differential equations containing singularities

!! --------------------------- [[[[[ Example Number ]]]]] ------------------------------
! Example 4.1 the fourth order PDE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    USE GEOMETRY
    USE PLOT
    USE BOUNDARY
    USE ERRORESTIMATE
    USE LOADFUNCTION

	REAL*8, ALLOCATABLE :: STIF_K(:,:), DMY_K(:,:), STIF_F(:), RES_K(:,:), DMY_F(:)
  REAL*8 :: ERROR, D, ENRG(4), MAX_NORM, L2_NORM(2)
  REAL :: CPTS, CPTE
 	INTEGER :: I, J, II, JJ
 	INTEGER, ALLOCATABLE :: INDX(:)

!-----------------!
! 	BL = 'Y'
! 	BL = 'N'
!-----------------!

!-----------------!
! 	USE_ENRICH = 'Y'
	USE_ENRICH = 'N'
!-----------------!

!-----------------!
! 	MESHTYPE = 'SH'
! 	MESHTYPE = 'NA'
!-----------------!

!-----------------!
! 	ZEROBD = 'Y'
! 	ZEROBD = 'N'
!-----------------!
	
	write(char_problem,fmt='(i1)') PROBLEM
	FILENAME = trim('./data/test') // PDE // char_problem

	OPEN(131, FILE=FILENAME, STATUS='UNKNOWN')
	
! 	WRITE(131,*) 'PDE=',PDE
! 	IF (PDE=='CNVD') THEN
! 		WRITE(131,*) 'Epsilon=', EPSLN
! 	ENDIF
! 	WRITE(131,*) 'Boundary Layer Problem?  ', BL
! 	WRITE(131,*) 'Do you use enriching functions?  ', USE_ENRICH
! 	WRITE(131,*) 'Do you Apply any special type of mesh refinement?  ', MESHTYPE
! 	WRITE(131,*) 'Do you impose zero boundary condition?  ', ZEROBD
! 	WRITE(131, *) 'DOF - BD_DOF  ', 'ORDER  ', 'MAX.  ', 'L2  ', 'ENRG'
    
	DO I = 3, 3
        
        BS_ORDER(1, :) = (/I, 2/)
		BS_ORDER(2, :) = (/I, I/)
		
!         EXTRA_KNOTS(1, :) = (/I, i/)
!         EXTRA_KNOTS(2, :) = (/i, i/)
		
		EXTRA_KNOTS(1, :) = (/5, 0/)
		EXTRA_KNOTS(2, :) = (/5, 0/)

! 		NURBS_REGULARITY(1, :) = (/BS_ORDER(1, 1)-1, BS_ORDER(1, 2)-1/)
! 		NURBS_REGULARITY(2, :) = (/BS_ORDER(2, 1)-1, BS_ORDER(2, 2)-1/)
		
		CALL CPU_TIME(CPTS)

		CALL GET_GEO()

! 		STOP
		
		ALLOCATE(STIF_K(DOF,DOF), DMY_K(DOF,DOF), RES_K(DOF,DOF), DMY_F(DOF), STIF_F(DOF), INDX(DOF))
		CALL PLOTGM()

! 		STOP
!         PRINT*, 'HERE-1'
        
		CALL GEN_KF(STIF_K, STIF_F)
		DMY_K = STIF_K
        
!         PRINT*, 'HERE-2'
        
		OPEN(1, FILE = './data/f1')
		DO II=1, DOF
			WRITE(1,*) STIF_F(II)
		ENDDO
		CLOSE(1)
	
		OPEN(1, FILE = './data/k1')
		DO II=1, DOF
			WRITE(1,*) (STIF_K(II,JJ), JJ=1,DOF)
		ENDDO
		CLOSE(1)
		
		CALL IMPOSEBD(STIF_K,STIF_F)

		OPEN(1, FILE = './data/f2')
		DO II=1, DOF
			WRITE(1,*) STIF_F(II)
		ENDDO
		CLOSE(1)
! 		
		OPEN(1, FILE = './data/k2')
		DO II=1, DOF
			WRITE(1,*) (STIF_K(II,JJ), JJ=1,DOF)
		ENDDO
		CLOSE(1)

		RES_K = STIF_K
		DMY_F = STIF_F

		CALL LUDCMP(STIF_K, DOF, DOF, INDX, D)
		CALL LUBKSB(STIF_K, DOF, DOF, INDX, STIF_F)
		
		CALL RESIDUALNORM(RES_K, STIF_F, DMY_F)

		OPEN(1, FILE = './data/sol')
		DO II = 1, DOF
			WRITE(1,*) STIF_F(II)
		ENDDO
		CLOSE(1)
		
		CALL MAXNORM(MAX_NORM, STIF_F)
		CALL ENRGY(ENRG, DMY_K, STIF_F)
! 		CALL L2_NORM_ERROR(L2_NORM, STIF_F)

		PRINT*, BS_ORDER(1,1), DOF - BD_DOF, MAX_NORM, ENRG(2), ENRG(4)
		WRITE(131, *) BS_ORDER(1,1), DOF - BD_DOF, MAX_NORM, ENRG(2), ENRG(4)

		CALL CPU_TIME(CPTE)

		PRINT*, ""
		PRINT*, "ELAPSED CPU TIME : ", CPTE - CPTS
		PRINT*, ""

		DEALLOCATE(STIF_K, DMY_K, RES_K, DMY_F, STIF_F, INDX, NDX)
	ENDDO
	CLOSE(131)
END PROGRAM MAIN


