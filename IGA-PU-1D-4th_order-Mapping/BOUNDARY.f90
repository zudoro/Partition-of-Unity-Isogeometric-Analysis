MODULE BOUNDARY

	USE GLBVAR

CONTAINS

SUBROUTINE IMPOSEBD(K,F)

	REAL*8, INTENT(INOUT) :: K(DOF,DOF), F(DOF)
	INTEGER :: I, J, II, JJ, KK


! ! FIND COEFFICIENTS OF BASIS FUNCTIONS ON BOUNDARY BY USING LEAST SQUARE METHOD
! ! 	OPEN(121, FILE = './data/least_square')
! 	
! 	CALL GEN_BD_LN()
! 
! 	102 FORMAT(1000(f20.8,1x))
! 	
! ! 	SUBTRACT EXACT SOLUTIONS FROM THE LOAD VECTOR
! 	DO I=1, BDNDX(1)%LC_NUM
! 		DO J=1, DOF
! 			F(J) = F(J) - LST_SOL(I)*K(J, BDNDX(I)%GL_NDX)
! 		ENDDO
! 	ENDDO
! 
! 	DO I=1, BDNDX(1)%LC_NUM
! 		F(BDNDX(I)%GL_NDX) = LST_SOL(I)
! 		DO J=1, DOF
! 			K(BDNDX(I)%GL_NDX,J) = 0.0D0
! 			K(J,BDNDX(I)%GL_NDX) = 0.0D0
! 		ENDDO
! 		K(BDNDX(I)%GL_NDX, BDNDX(I)%GL_NDX) = 1.0D0
! 	ENDDO
! 	
! 	DO I = 1, NUMBS(1)
! 		F(I) = 0.D0
! 		DO J = 1, DOF
! 			K(I,J) = 0.D0
! 			K(J,I) = 0.D0
! 		ENDDO
! 		K(I,I) = 1.D0
! 	ENDDO
	
! 	F(1) = 0.0D0
! 	F(2) = 0.0D0
	
	F(DOF-1) = 0.0D0
	F(DOF) = 0.0D0
	
	DO I = 1, DOF
! 		K(1,I) = 0.0D0
! 		K(2,I) = 0.0D0
! 		
! 		K(I,1) = 0.0D0
! 		K(I,2) = 0.0D0
		
		K(DOF,I) = 0.0D0
		K(DOF-1, I) = 0.0D0
		
		K(I,DOF) = 0.0D0
		K(I, DOF-1) = 0.0D0
	ENDDO
	
! 	K(1,1) = 1.0D0
! 	K(2,2) = 1.0D0
	
	K(DOF-1, DOF-1) = 1.0D0
	K(DOF,DOF) = 1.0D0
	
	WRITE(*,*)
	WRITE(*,*) '<<< IMPOSING ESSENTIAL BOUNDARY CONDITIONS BY USING LEAST SQUARE METHOD : DONE >>>'
	WRITE(*,*)

END SUBROUTINE IMPOSEBD

END MODULE BOUNDARY
