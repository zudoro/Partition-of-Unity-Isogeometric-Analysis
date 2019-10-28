MODULE BOUNDARY

	USE INTEGRATION
	USE LUDECOMPOSITION

    IMPLICIT NONE

CONTAINS

SUBROUTINE IMPOSEBD(K,F)

	REAL*8, INTENT(INOUT) :: K(DOF,DOF), F(DOF)
	REAL*8 :: SUB_K(BDNDX(1)%LC_NUM,BDNDX(1)%LC_NUM), DMY_SUB_K(BDNDX(1)%LC_NUM, BDNDX(1)%LC_NUM)
	REAL*8 :: DD
	REAL*8 :: SUB_F(BDNDX(1)%LC_NUM)
	INTEGER :: I, J, II, JJ, KK, INDX(BDNDX(1)%LC_NUM)

!----------------------------------- [[[[[ IMPOSE NON-HOMOGENEOUS BOUNDARY CONDITION ]]]]] ---------------------------------------
!---------------------------------------------------------------------------------------------------------------------------------
!! FIND COEFFICIENTS OF BASIS FUNCTIONS ON BOUNDARY BY USING LEAST SQUARE METHOD

	SUB_K(:,:) = 0.D0
	SUB_F(:) = 0.D0
! 	OPEN(121, FILE = './data/least_square')
	
! 	IF (ZEROBD=='N') THEN
! 		CALL GEN_BD_LN(SUB_K, SUB_F)

! 		OPEN(124, FILE = './data/sub_k')
! 		OPEN(125, FILE = './data/sub_f')
! 		DO I=1, BDNDX(1)%LC_NUM
! 			WRITE(124, *) (SUB_K(I,J), J=1,BDNDX(1)%LC_NUM)
! 			WRITE(125, *) SUB_F(I)
! 		ENDDO
! 		CLOSE(124)
! 		CLOSE(125)
		
	! 	102 FORMAT(1000(f20.8,1x))
	! 	
! 		CALL LUDCMP(SUB_K, BDNDX(1)%LC_NUM, BDNDX(1)%LC_NUM, INDX, DD)
! 		CALL LUBKSB(SUB_K, BDNDX(1)%LC_NUM, BDNDX(1)%LC_NUM, INDX, SUB_F)
		
! 		OPEN(126, FILE = './data/sub_s')
! 		DO I=1, BDNDX(1)%LC_NUM
! 			WRITE(126, *) SUB_F(I)
! 		ENDDO
! 		CLOSE(126)

	!---------------------------------------------------------------------------------------------------------------------------------

		!! SUBTRACT EXACT SOLUTIONS FROM THE LOAD VECTOR
! 		DO I=1, BDNDX(1)%LC_NUM
! 			DO J=1, DOF
! 				F(J) = F(J) - SUB_F(I)*K(J, BDNDX(I)%GL_NDX)
! 			ENDDO
! 		ENDDO
! 
! 		DO I=1, BDNDX(1)%LC_NUM
! 			F(BDNDX(I)%GL_NDX) = SUB_F(I)
! 			DO J=1, DOF
! 				K(BDNDX(I)%GL_NDX,J) = 0.0D0
! 				K(J,BDNDX(I)%GL_NDX) = 0.0D0
! 			ENDDO
! 			K(BDNDX(I)%GL_NDX, BDNDX(I)%GL_NDX) = 1.0D0
! 		ENDDO

	!---------------------------------------------------------------------------------------------------------------------------------
! 	ENDIF
	
!! CONSTRAINTS EQN. (21)
! 	DO I = DOF - NUMBS(1) + 2, DOF
! 		DO J = DOF - NUMBS(1) + 2, DOF
! 			IF (I==J) THEN
! 				K(I,J) = 1.0D0
! 			ELSE 
! 				K(I,J) = 0.0D0
! 			ENDIF
! 		ENDDO
! 		K(I, DOF - NUMBS(1) + 1) = -1.0D0
! 		F(I) = 0.0D0
! 	ENDDO
		
!!--------------------------------------------------------------------------------------------------------------------!!
!! Impose zero Dirichlet boundary conditions

! 	IF (ZEROBD=='Y') THEN
		DO I=1, ZERO_BDNDX(1)%LC_NUM
			F(ZERO_BDNDX(I)%GL_NDX) = 0.D0
	! 			PRINT*, F(ZERO_BDNDX(II,I)%GL_NDX), F(DOF+BDNDX(II,I)%GL_NDX)
			DO J=1, DOF
				K(ZERO_BDNDX(I)%GL_NDX, J) = 0.0D0
				K(J, ZERO_BDNDX(I)%GL_NDX) = 0.0D0
			ENDDO
			K(ZERO_BDNDX(I)%GL_NDX, ZERO_BDNDX(I)%GL_NDX) = 1.0D0
		ENDDO
! 	ELSE 
! 	ENDIF

	WRITE(*,*)
	WRITE(*,*) '<<< IMPOSING ESSENTIAL BOUNDARY CONDITIONS BY USING LEAST SQUARE METHOD : DONE >>>'
	WRITE(*,*)

END SUBROUTINE IMPOSEBD

END MODULE BOUNDARY