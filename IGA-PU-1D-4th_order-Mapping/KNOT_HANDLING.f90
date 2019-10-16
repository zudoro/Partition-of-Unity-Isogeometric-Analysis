MODULE KNOT_HANDLING

	USE NEWTYPE

	IMPLICIT INTEGER (I-N)
	IMPLICIT REAL(8) (A-H,O-Z)

	INTERFACE OPERATOR (+)
		MODULE PROCEDURE TRANSLATE_KNOT_VECTOR
	END INTERFACE

	INTERFACE OPERATOR (*)
		MODULE PROCEDURE RESCALE_KNOT_VECTOR
	END INTERFACE

	INTERFACE OPERATOR (.NORMAL.)
		MODULE PROCEDURE NORMALIZE_KNOT_VECTOR
	END INTERFACE

CONTAINS

	!!  FIND KNOT SPAN
	!!  INPUT  : PT - POINT, KVEC - A VALID KNOT VECTOR
	!!  OUTPUT : "INTEGER" KNOT SPAN INDEX
	INTEGER FUNCTION FIND_KNOT_SPAN(PT,KVEC)

		REAL(8), INTENT(IN) :: PT
		TYPE(KNOT_VECTOR), INTENT(IN) :: KVEC

		INTEGER :: LOW, HIGH, MID

		IF ((PT<KVEC%KNOTS(0)-TOLERANCE) .OR. (PT>KVEC%KNOTS(KVEC%LENGTH)+TOLERANCE)) THEN
			PRINT *, '[ERROR]  GIVEN POINT IS OUT OF RANGE !'
			FIND_KNOT_SPAN = -1
		ELSE IF (ABS(PT-KVEC%KNOTS(KVEC%LENGTH-KVEC%POLY_ORDER))<TOLERANCE) THEN
			FIND_KNOT_SPAN = KVEC%LENGTH - KVEC%POLY_ORDER - 1
		ELSE
			LOW = KVEC%POLY_ORDER
			HIGH = KVEC%LENGTH - KVEC%POLY_ORDER
			MID = (LOW+HIGH)/2
			DO WHILE (PT<KVEC%KNOTS(MID) .OR. PT>=KVEC%KNOTS(MID+1))
				IF (PT<KVEC%KNOTS(MID)) THEN
					HIGH = MID
				ELSE
					LOW = MID
				ENDIF
				MID = (LOW+HIGH)/2
			ENDDO
			FIND_KNOT_SPAN = MID
		ENDIF

	END FUNCTION FIND_KNOT_SPAN
	
	!!  GENERATE BASE OPEN KNOT VECTOR
	TYPE(KNOT_VECTOR) FUNCTION GET_BASE_OPEN_KNOT_VECTOR(POLY_ORDER,NUMBER_OF_BASIS)

		INTEGER, INTENT(IN) :: POLY_ORDER, NUMBER_OF_BASIS

		INTEGER :: NADD

		IF (NUMBER_OF_BASIS>0 .AND. NUMBER_OF_BASIS<(POLY_ORDER+1)) THEN
			PRINT *, '[ERROR]  GIVEN NUMBER OF BASIS IS TOO SMALL !'
			NADD = 0
		ENDIF

		IF (NUMBER_OF_BASIS<=0) THEN
			NADD = 0
		ELSE
			NADD = NUMBER_OF_BASIS - (POLY_ORDER+1)
		ENDIF

		GET_BASE_OPEN_KNOT_VECTOR%POLY_ORDER = POLY_ORDER
		GET_BASE_OPEN_KNOT_VECTOR%LENGTH = -1
		DO I = 0, POLY_ORDER
			GET_BASE_OPEN_KNOT_VECTOR%KNOTS(GET_BASE_OPEN_KNOT_VECTOR%LENGTH+1) = 0.0D0
			GET_BASE_OPEN_KNOT_VECTOR%LENGTH = GET_BASE_OPEN_KNOT_VECTOR%LENGTH + 1
		ENDDO
		DO I = 1, NADD
			GET_BASE_OPEN_KNOT_VECTOR%KNOTS(GET_BASE_OPEN_KNOT_VECTOR%LENGTH+1) = (I*1.0D0)/(NADD+1)
			GET_BASE_OPEN_KNOT_VECTOR%LENGTH = GET_BASE_OPEN_KNOT_VECTOR%LENGTH + 1
		ENDDO
		DO I = 0, POLY_ORDER
			GET_BASE_OPEN_KNOT_VECTOR%KNOTS(GET_BASE_OPEN_KNOT_VECTOR%LENGTH+1) = 1.0D0
			GET_BASE_OPEN_KNOT_VECTOR%LENGTH = GET_BASE_OPEN_KNOT_VECTOR%LENGTH + 1
		ENDDO

	END FUNCTION GET_BASE_OPEN_KNOT_VECTOR

	!!  GENERATE BASE OPEN KNOT VECTOR
	TYPE(KNOT_VECTOR) FUNCTION GET_OPEN_KNOT_VECTOR(POLY_ORDER)

		INTEGER, INTENT(IN) :: POLY_ORDER

		GET_OPEN_KNOT_VECTOR%POLY_ORDER = POLY_ORDER
		GET_OPEN_KNOT_VECTOR%LENGTH = -1
		DO I = 0, POLY_ORDER
			GET_OPEN_KNOT_VECTOR%KNOTS(GET_OPEN_KNOT_VECTOR%LENGTH+1) = 0.0D0
			GET_OPEN_KNOT_VECTOR%LENGTH = GET_OPEN_KNOT_VECTOR%LENGTH + 1
		ENDDO
		DO I = 0, POLY_ORDER
			GET_OPEN_KNOT_VECTOR%KNOTS(GET_OPEN_KNOT_VECTOR%LENGTH+1) = 1.0D0
			GET_OPEN_KNOT_VECTOR%LENGTH = GET_OPEN_KNOT_VECTOR%LENGTH + 1
		ENDDO

	END FUNCTION GET_OPEN_KNOT_VECTOR


	!!  TRANSLATE KNOT VECTOR
	TYPE(KNOT_VECTOR) FUNCTION TRANSLATE_KNOT_VECTOR(PT,KVEC)

		REAL(8), INTENT(IN) :: PT
		TYPE(KNOT_VECTOR), INTENT(IN) :: KVEC

		DO K = 0, KVEC%LENGTH
			TRANSLATE_KNOT_VECTOR%KNOTS(K) = KVEC%KNOTS(K) + PT
		ENDDO
		TRANSLATE_KNOT_VECTOR%LENGTH = KVEC%LENGTH
		TRANSLATE_KNOT_VECTOR%POLY_ORDER = KVEC%POLY_ORDER

	END FUNCTION TRANSLATE_KNOT_VECTOR


	!!  RESCALE KNOT VECTOR
	TYPE(KNOT_VECTOR) FUNCTION RESCALE_KNOT_VECTOR(SCALAR,KVEC)

		REAL(8), INTENT(IN) :: SCALAR
		TYPE(KNOT_VECTOR), INTENT(IN) :: KVEC

		DO K = 0, KVEC%LENGTH
			RESCALE_KNOT_VECTOR%KNOTS(K) = SCALAR * KVEC%KNOTS(K)
		ENDDO
		RESCALE_KNOT_VECTOR%LENGTH = KVEC%LENGTH
		RESCALE_KNOT_VECTOR%POLY_ORDER = KVEC%POLY_ORDER

	END FUNCTION RESCALE_KNOT_VECTOR


	!!  NORMALIZE KNOT VECTOR
	TYPE(KNOT_VECTOR) FUNCTION NORMALIZE_KNOT_VECTOR(KVEC)

		TYPE(KNOT_VECTOR), INTENT(IN) :: KVEC

		REAL(8) :: FACTOR

		FACTOR = KVEC%KNOTS(KVEC%LENGTH) - KVEC%KNOTS(0)
		DO K = 0, KVEC%LENGTH
			NORMALIZE_KNOT_VECTOR%KNOTS(K) = (KVEC%KNOTS(K)-KVEC%KNOTS(0)) / FACTOR
		ENDDO
		NORMALIZE_KNOT_VECTOR%LENGTH = KVEC%LENGTH
		NORMALIZE_KNOT_VECTOR%POLY_ORDER = KVEC%POLY_ORDER

	END FUNCTION NORMALIZE_KNOT_VECTOR


	!!  INSERT KNOT VALUES. DOES NOT CHANGE THE POLYNOMIAL ORDER
	TYPE(KNOT_VECTOR) FUNCTION KNOT_INSERTION(KVEC)

		TYPE(KNOT_VECTOR), INTENT(IN) :: KVEC

		REAL(8) :: KNOTS(MAX_LENGTH), NEW_KNOT
		INTEGER :: POLY_ORDER, MULTIPLICITIES(MAX_LENGTH), NUM_KNOTS, MENU
		LOGICAL :: INSERTED

		KNOT_INSERTION = KVEC
		CALL KNOT_TO_ARRAY(KVEC,POLY_ORDER,NUM_KNOTS,KNOTS,MULTIPLICITIES)
		WRITE(*,*)
		WRITE(*,FMT='(A,I2,A)') ' ** THERE ARE ', NUM_KNOTS, ' DISTINCT KNOT VALUES IN THE GIVEN KNOT VECTOR.'
		WRITE(*,*)
		DO K = 1, NUM_KNOTS
			WRITE(*,101) K, KNOTS(K), MULTIPLICITIES(K)
		ENDDO
		WRITE(*,*)
		100 WRITE(*,FMT='(A)',ADVANCE='NO') ' ** SELECT ONE  :  [1] INSERT A KNOT, [0] QUIT  >> '
		READ(*,*,ERR=100) MENU
		IF (MENU/=1 .AND. MENU/=0) GOTO 100
		WRITE(*,*)
		DO WHILE (MENU/=0)
			WRITE(*,FMT='(A)',ADVANCE='NO') '    ENTER A KNOT FROM (0,1) >> '
			READ(*,*) NEW_KNOT
			IF (NEW_KNOT>0.0D0 .AND. NEW_KNOT<1.0D0) THEN
				INSERTED = .FALSE.
				DO K = 2, NUM_KNOTS
					IF ((.NOT.INSERTED) .AND. ABS(NEW_KNOT-KNOTS(K))<TOLERANCE) THEN
						IF (MULTIPLICITIES(K)/=KVEC%POLY_ORDER) THEN
							MULTIPLICITIES(K) = MULTIPLICITIES(K) + 1
							INSERTED = .TRUE.
						ENDIF
					ELSE IF ((.NOT.INSERTED) .AND. (NEW_KNOT>KNOTS(K-1)) .AND. (NEW_KNOT<KNOTS(K)))THEN
						DO J = NUM_KNOTS, K, -1
							KNOTS(J+1) = KNOTS(J)
							MULTIPLICITIES(J+1) = MULTIPLICITIES(J)
						ENDDO
						KNOTS(K) = NEW_KNOT
						MULTIPLICITIES(K) = 1
						NUM_KNOTS = NUM_KNOTS + 1
					ENDIF
				ENDDO
			ENDIF
			WRITE(*,*)
			WRITE(*,FMT='(A,I2,A)') ' ** THERE ARE ', NUM_KNOTS, ' DISTINCT KNOT VALUES IN THE GIVEN KNOT VECTOR.'
			WRITE(*,*)
			DO K = 1, NUM_KNOTS
				WRITE(*,101) K, KNOTS(K), MULTIPLICITIES(K)
			ENDDO
			WRITE(*,*)
			200 WRITE(*,FMT='(A)',ADVANCE='NO') ' ** SELECT ONE  :  [1] INSERT A KNOT, [9] SAVE AND QUIT, [0] QUIT  >> '
			READ(*,*,ERR=200) MENU
			IF (MENU/=1 .AND. MENU/=9 .AND. MENU/=0) GOTO 200
			WRITE(*,*)
			IF (MENU==9) THEN
				KNOT_INSERTION = ARRAY_TO_KNOT(POLY_ORDER,NUM_KNOTS,KNOTS,MULTIPLICITIES)
				MENU = 0
			ENDIF
		ENDDO

		101 FORMAT(5X,'[',I2.2,']',2XF8.5,5X,'X ',I2)

	END FUNCTION KNOT_INSERTION


	!!  INSERT KNOT VALUES. DOES NOT CHANGE THE POLYNOMIAL ORDER
	TYPE(KNOT_VECTOR) FUNCTION UNIFORM_KNOT_INSERTION(KVEC,NUM_INSERTION, REGULARITY)

		TYPE(KNOT_VECTOR), INTENT(IN) :: KVEC
		INTEGER, INTENT(IN) :: NUM_INSERTION
		INTEGER, OPTIONAL, INTENT(IN) :: REGULARITY

		REAL(8) :: KNOTS(MAX_LENGTH), NEW_KNOTS(MAX_LENGTH)
		INTEGER :: POLY_ORDER, MULTIPLICITIES(MAX_LENGTH), NEW_MULTIPLICITIES(MAX_LENGTH), NUM_KNOTS

		CALL KNOT_TO_ARRAY(KVEC,POLY_ORDER,NUM_KNOTS,KNOTS,MULTIPLICITIES)
		NEW_KNOTS(1) = KNOTS(1)
		NEW_MULTIPLICITIES(1) = MULTIPLICITIES(1)
		DO K = 2, NUM_KNOTS
			DO J = 1, NUM_INSERTION
				INDX = (K-1)*(NUM_INSERTION+1)-NUM_INSERTION+J
				NEW_KNOTS(INDX) = (1.0D0*J)*(KNOTS(K)-KNOTS(K-1))/(NUM_INSERTION+1) + KNOTS(K-1)
				IF (PRESENT(REGULARITY)) THEN
					IF (REGULARITY.GE.POLY_ORDER) THEN
						PRINT*, 'ERROR-UNIFORM_KNOT_INSERTION : REGULARITY IS GREATER THAN POLY_ORDER'
						STOP
					ELSE
						NEW_MULTIPLICITIES(INDX) = POLY_ORDER - REGULARITY
					ENDIF
				ELSE
					NEW_MULTIPLICITIES(INDX) = 1
				ENDIF
			ENDDO
			NEW_KNOTS(K*(NUM_INSERTION+1)-NUM_INSERTION) = KNOTS(K)
			NEW_MULTIPLICITIES(K*(NUM_INSERTION+1)-NUM_INSERTION) = MULTIPLICITIES(K)
		ENDDO
		UNIFORM_KNOT_INSERTION = ARRAY_TO_KNOT(POLY_ORDER,NUM_KNOTS*(NUM_INSERTION+1)-NUM_INSERTION,NEW_KNOTS,NEW_MULTIPLICITIES)

	END FUNCTION UNIFORM_KNOT_INSERTION
	
	!!  INSERT KNOT VALUES. DOES NOT CHANGE THE POLYNOMIAL ORDER
	TYPE(KNOT_VECTOR) FUNCTION DEGREE_ELEVATION(KVEC,REGULARITY_OPTION)

		TYPE(KNOT_VECTOR), INTENT(IN) :: KVEC
		LOGICAL, OPTIONAL, INTENT(IN) :: REGULARITY_OPTION

		REAL(8) :: KNOTS(MAX_LENGTH)
		INTEGER :: POLY_ORDER, MULTIPLICITIES(MAX_LENGTH), NUM_KNOTS
		LOGICAL :: INCREASE_REGULARITY

		CALL KNOT_TO_ARRAY(KVEC,POLY_ORDER,NUM_KNOTS,KNOTS,MULTIPLICITIES)
		IF (PRESENT(REGULARITY_OPTION)) THEN
			IF (REGULARITY_OPTION) THEN
				INCREASE_REGULARITY = .TRUE.
			ELSE
				INCREASE_REGULARITY = .FALSE.
			ENDIF
		ELSE
				INCREASE_REGULARITY = .FALSE.
		ENDIF

		IF (INCREASE_REGULARITY) THEN
			MULTIPLICITIES(1) = MULTIPLICITIES(1) + 1
			MULTIPLICITIES(NUM_KNOTS) = MULTIPLICITIES(NUM_KNOTS) + 1
			DEGREE_ELEVATION = ARRAY_TO_KNOT(POLY_ORDER+1,NUM_KNOTS,KNOTS,MULTIPLICITIES)
		ELSE
			DO K = 1, NUM_KNOTS
				MULTIPLICITIES(K) = MULTIPLICITIES(K) + 1
			ENDDO
			DEGREE_ELEVATION = ARRAY_TO_KNOT(POLY_ORDER+1,NUM_KNOTS,KNOTS,MULTIPLICITIES)
		ENDIF

	END FUNCTION DEGREE_ELEVATION

	!!  INSERT KNOT VALUES. DOES NOT CHANGE THE POLYNOMIAL ORDER
	TYPE(KNOT_VECTOR) FUNCTION DEGREE_ELEVATION_WITH_CONTROL_REGULARITY(KVEC,REGULARITY_OPTION)

		TYPE(KNOT_VECTOR), INTENT(IN) :: KVEC
		INTEGER, OPTIONAL, INTENT(IN) :: REGULARITY_OPTION

		REAL(8) :: KNOTS(MAX_LENGTH)
		INTEGER :: POLY_ORDER, MULTIPLICITIES(MAX_LENGTH), NUM_KNOTS

		CALL KNOT_TO_ARRAY(KVEC,POLY_ORDER,NUM_KNOTS,KNOTS,MULTIPLICITIES)
		IF (PRESENT(REGULARITY_OPTION)) THEN
			IF (REGULARITY_OPTION.EQ.-1) THEN
				MULTIPLICITIES(1) = MULTIPLICITIES(1) + 1
				MULTIPLICITIES(NUM_KNOTS) = MULTIPLICITIES(NUM_KNOTS) + 1
				DEGREE_ELEVATION_WITH_CONTROL_REGULARITY = ARRAY_TO_KNOT(POLY_ORDER+1,NUM_KNOTS,KNOTS,MULTIPLICITIES)
			ELSE
				DO K = 1, NUM_KNOTS
					IF (MULTIPLICITIES(K).LE.POLY_ORDER) THEN
						MULTIPLICITIES(K) = POLY_ORDER + 1 - REGULARITY_OPTION
					ELSE
						MULTIPLICITIES(K) = MULTIPLICITIES(K) + 1
					ENDIF
				ENDDO
				DEGREE_ELEVATION_WITH_CONTROL_REGULARITY = ARRAY_TO_KNOT(POLY_ORDER+1,NUM_KNOTS,KNOTS,MULTIPLICITIES)
			ENDIF
		ELSE
			MULTIPLICITIES(1) = MULTIPLICITIES(1) + 1
			MULTIPLICITIES(NUM_KNOTS) = MULTIPLICITIES(NUM_KNOTS) + 1
			DEGREE_ELEVATION_WITH_CONTROL_REGULARITY = ARRAY_TO_KNOT(POLY_ORDER+1,NUM_KNOTS,KNOTS,MULTIPLICITIES)
		ENDIF

	END FUNCTION DEGREE_ELEVATION_WITH_CONTROL_REGULARITY
	
	!!  CONVERT A KNOT VECTOR TO LISTS OF KNOTS AND MULTIPLICITIES
	SUBROUTINE KNOT_TO_ARRAY(KVEC,POLY_ORDER,NUM_KNOTS,KNOTS,MULTIPLICITIES)

		TYPE(KNOT_VECTOR), INTENT(IN) :: KVEC
		INTEGER, INTENT(OUT) :: POLY_ORDER, NUM_KNOTS, MULTIPLICITIES(MAX_LENGTH)
		REAL(8), INTENT(OUT) :: KNOTS(MAX_LENGTH)

		POLY_ORDER = KVEC%POLY_ORDER
		NUM_KNOTS = 1
		KNOTS(NUM_KNOTS) = KVEC%KNOTS(0)
		MULTIPLICITIES(NUM_KNOTS) = 1
		DO K = 1, KVEC%LENGTH
			IF (ABS(KNOTS(NUM_KNOTS)-KVEC%KNOTS(K))<TOLERANCE) THEN
				MULTIPLICITIES(NUM_KNOTS) = MULTIPLICITIES(NUM_KNOTS) + 1
			ELSE
				NUM_KNOTS = NUM_KNOTS + 1
				KNOTS(NUM_KNOTS) = KVEC%KNOTS(K)
				MULTIPLICITIES(NUM_KNOTS) = 1
			ENDIF
		ENDDO

	END SUBROUTINE KNOT_TO_ARRAY


	!!  CONVERT LISTS OF KNOTS AND MULTIPLICITIES TO A KNOT VECTOR
	TYPE(KNOT_VECTOR) FUNCTION ARRAY_TO_KNOT(POLY_ORDER,NUM_KNOTS,KNOTS,MULTIPLICITIES)

		INTEGER, INTENT(IN) :: POLY_ORDER, NUM_KNOTS, MULTIPLICITIES(MAX_LENGTH)
		REAL(8), INTENT(IN) :: KNOTS(MAX_LENGTH)

		ARRAY_TO_KNOT%POLY_ORDER = POLY_ORDER
		ARRAY_TO_KNOT%LENGTH = -1
		DO K = 1, NUM_KNOTS
			DO J = 1, MULTIPLICITIES(K)
				ARRAY_TO_KNOT%LENGTH = ARRAY_TO_KNOT%LENGTH + 1
				ARRAY_TO_KNOT%KNOTS(ARRAY_TO_KNOT%LENGTH) = KNOTS(K)
			ENDDO
		ENDDO

	END FUNCTION ARRAY_TO_KNOT

	SUBROUTINE P_REFINEMENT(OLD_CTL, P, AXIS)
	! P = CURRENT POLYNOMIAL DEGREE
	! P_REFINEMENT
	! FOR THE CASE OF BERNSTEIN POLYNOMIAL

		TYPE(CONTROL_POINTS_2D), INTENT(INOUT) :: OLD_CTL
		INTEGER, INTENT(IN) :: P(2)
		CHARACTER(LEN=1), INTENT(IN) :: AXIS
		
		TYPE(CONTROL_POINTS_2D) :: NEW_CTL
		INTEGER :: I, J, K, II, JJ, KK
		
		IF (OLD_CTL%D.NE.2) THEN
			PRINT*, 'ERROR - FIND_CTL_WGT_BERNSTEIN : DIMENSION OF INPUT CONTROL POINT IS WRONG'
			STOP
		ENDIF

! 		DO I = 0, P(1)
! 			PRINT*, I, OLD_CTL%WGTS(I, 0, 0)
! 		ENDDO
! 		PRINT*, 

		NEW_CTL%D = 2
		
		IF (AXIS=='X') THEN
			DO J = 0, P(2)
				NEW_CTL%PTS(0, J, 0) = OLD_CTL%PTS(0, J, 0)
				NEW_CTL%PTS(P(1)+1, J, 0) = OLD_CTL%PTS(P(1), J, 0)
				
				NEW_CTL%WGTS(0, J, 0) = OLD_CTL%WGTS(0, J, 0)
				NEW_CTL%WGTS(P(1)+1, J, 0) = OLD_CTL%WGTS(P(1), J, 0)
! 				PRINT*, J, OLD_CTL%WGTS(P(1), J, 0), NEW_CTL%WGTS(P(1)+1, J, 0)
			ENDDO
			DO J = 0, P(2)
				DO I = 1, P(1)
					NEW_CTL%WGTS(I,J,0) = ((P(1)+1-I)*OLD_CTL%WGTS(I,J,0) + I*OLD_CTL%WGTS(I-1,J,0))/(1.D0*(P(1)+1))
					
					NEW_CTL%PTS(I,J,0)%X = ((P(1)+1-I)*OLD_CTL%WGTS(I,J,0)*OLD_CTL%PTS(I,J,0)%X + I*OLD_CTL%WGTS(I-1,J,0)*OLD_CTL%PTS(I-1,J,0)%X)/(1.D0*(P(1)+1)*NEW_CTL%WGTS(I,J,0))
					NEW_CTL%PTS(I,J,0)%Y = ((P(1)+1-I)*OLD_CTL%WGTS(I,J,0)*OLD_CTL%PTS(I,J,0)%Y + I*OLD_CTL%WGTS(I-1,J,0)*OLD_CTL%PTS(I-1,J,0)%Y)/(1.D0*(P(1)+1)*NEW_CTL%WGTS(I,J,0))
				ENDDO
			ENDDO
			
		ELSEIF (AXIS=='Y') THEN
			DO I = 0, P(1)
				NEW_CTL%PTS(I, 0, 0) = OLD_CTL%PTS(I, 0, 0)
				NEW_CTL%PTS(I, P(2)+1, 0) = OLD_CTL%PTS(I, P(2), 0)
				
				NEW_CTL%WGTS(I, 0, 0) = OLD_CTL%WGTS(I, 0, 0)
				NEW_CTL%WGTS(I, P(2)+1, 0) = OLD_CTL%WGTS(I, P(2), 0)
			ENDDO
			DO I = 0, P(1)
				DO J = 1, P(2)
					NEW_CTL%WGTS(I,J,0) = ((P(2)+1-J)*OLD_CTL%WGTS(I,J,0) + J*OLD_CTL%WGTS(I,J-1,0))/(1.D0*(P(2)+1))
					
					NEW_CTL%PTS(I,J,0)%X = ((P(2)+1-J)*OLD_CTL%WGTS(I,J,0)*OLD_CTL%PTS(I,J,0)%X + J*OLD_CTL%WGTS(I,J-1,0)*OLD_CTL%PTS(I,J-1,0)%X)/(1.D0*(P(2)+1)*NEW_CTL%WGTS(I,J,0))
					NEW_CTL%PTS(I,J,0)%Y = ((P(2)+1-J)*OLD_CTL%WGTS(I,J,0)*OLD_CTL%PTS(I,J,0)%Y + J*OLD_CTL%WGTS(I,J-1,0)*OLD_CTL%PTS(I,J-1,0)%Y)/(1.D0*(P(2)+1)*NEW_CTL%WGTS(I,J,0))
				ENDDO
			ENDDO
		ENDIF
		
		OLD_CTL = NEW_CTL
		
	END SUBROUTINE P_REFINEMENT
	
	SUBROUTINE H_REFINEMENT(OLD_CTL, NEW_KVEC, NEW_KNOT, AXIS)
	! h_REFINEMENT
	! FOR THE CASE OF BERNSTEIN POLYNOMIAL
		
		TYPE(CONTROL_POINTS_2D), INTENT(INOUT) :: OLD_CTL
		TYPE(KNOT_VECTOR), INTENT(IN) :: NEW_KVEC(2)
		REAL*8, INTENT(IN) :: NEW_KNOT
		CHARACTER(LEN=1), INTENT(IN) :: AXIS
		
		TYPE(CONTROL_POINTS_2D) :: NEW_CTL
		INTEGER :: I, J, K, II, JJ, KK, N(2)
		REAL*8 :: ALPHA, OLD_KNOTS(MAX_LENGTH)
		
		IF (OLD_CTL%D.NE.2) THEN
			PRINT*, 'ERROR - FIND_CTL_WGT_BERNSTEIN : DIMENSION OF INPUT CONTROL POINT IS WRONG'
			STOP
		ENDIF

		NEW_CTL%D = 2
		
		IF (AXIS=='X') THEN
			N(1) = NEW_KVEC(1)%LENGTH - NEW_KVEC(1)%POLY_ORDER - 1
			N(2) = NEW_KVEC(2)%LENGTH - NEW_KVEC(2)%POLY_ORDER
			
			II = 0
			DO I = 0, NEW_KVEC(1)%LENGTH
				IF (DABS(NEW_KVEC(1)%KNOTS(I) - NEW_KNOT).LE.EPS) THEN
					K = I
				ELSE
					II = II + 1
					OLD_KNOTS(II) = NEW_KVEC(1)%KNOTS(I)
				ENDIF
			ENDDO
			
! 			PRINT*, 'K', K, 'N', N(:)
! 			PRINT*, (K - NEW_KVEC(1)%POLY_ORDER), (N(1) + NEW_KVEC(1)%POLY_ORDER + 2)

! 			DO I = 1, NEW_KVEC(1)%LENGTH
! 				PRINT*, OLD_KNOTS(I)
! 			ENDDO
			
			IF (II.NE.NEW_KVEC(1)%LENGTH) THEN
				PRINT*, 'ERROR - H_REFINEMENT (AXIS = X) : NUMBER OF KNOTS IN OLD KNOT VECTOR IS NOT EQUAL TO THE NUMBER OF KNOTS IN NEW KNOT VECTOR - 1'
				STOP
			ENDIF
			
			DO J = 0, N(2) - 1
				DO I = 1, N(1) + 1
					IF (I.GE.1 .AND. I.LE.(K - NEW_KVEC(1)%POLY_ORDER)) THEN
					
						NEW_CTL%PTS(I-1, J, 0)%X = OLD_CTL%PTS(I-1, J, 0)%X
						NEW_CTL%PTS(I-1, J, 0)%Y = OLD_CTL%PTS(I-1, J, 0)%Y
						NEW_CTL%WGTS(I-1, J, 0) = OLD_CTL%WGTS(I-1, J, 0)
					
! 						PRINT*, I, NEW_CTL%PTS(I-1, J, 0), NEW_CTL%WGTS(I-1, J, 0)
					
					ELSEIF (I.GE.(K - NEW_KVEC(1)%POLY_ORDER + 1) .AND. I.LE.K) THEN
						
						ALPHA = (NEW_KNOT - OLD_KNOTS(I))/(OLD_KNOTS(I + NEW_KVEC(1)%POLY_ORDER) - OLD_KNOTS(I))
						
						NEW_CTL%WGTS(I-1, J, 0) = ALPHA*OLD_CTL%WGTS(I-1, J, 0) + (1.D0 - ALPHA)*OLD_CTL%WGTS(I-2, J, 0)
						
						NEW_CTL%PTS(I-1, J, 0)%X = (ALPHA*OLD_CTL%WGTS(I-1, J, 0)/NEW_CTL%WGTS(I-1, J, 0))*OLD_CTL%PTS(I-1, J, 0)%X + &
																			 ((1.D0 - ALPHA)*OLD_CTL%WGTS(I-2, J, 0)/NEW_CTL%WGTS(I-1, J, 0))*OLD_CTL%PTS(I-2, J, 0)%X
																			 
						NEW_CTL%PTS(I-1, J, 0)%Y = (ALPHA*OLD_CTL%WGTS(I-1, J, 0)/NEW_CTL%WGTS(I-1, J, 0))*OLD_CTL%PTS(I-1, J, 0)%Y + &
																			 ((1.D0 - ALPHA)*OLD_CTL%WGTS(I-2, J, 0)/NEW_CTL%WGTS(I-1, J, 0))*OLD_CTL%PTS(I-2, J, 0)%Y
						
! 						NEW_CTL%WGTS(I-1, J, 0) = DCOS(22.5D0*DACOS(-1.D0)/180.D0)
						
! 						PRINT*, I, NEW_CTL%PTS(I-1, J, 0), NEW_CTL%WGTS(I-1, J, 0)
					
					ELSEIF (I.GE.(K + 1) .AND. I.LE.(N(1) + NEW_KVEC(1)%POLY_ORDER + 2)) THEN
						
						NEW_CTL%PTS(I-1, J, 0)%X = OLD_CTL%PTS(I-2, J, 0)%X
						NEW_CTL%PTS(I-1, J, 0)%Y = OLD_CTL%PTS(I-2, J, 0)%Y
						NEW_CTL%WGTS(I-1, J, 0) = OLD_CTL%WGTS(I-2, J, 0)
						
! 						PRINT*, I, NEW_CTL%PTS(I-1, J, 0), NEW_CTL%WGTS(I-1, J, 0)
					
					ENDIF
				ENDDO
! 				STOP
			ENDDO
		ELSEIF (AXIS=='Y') THEN
			N(1) = NEW_KVEC(1)%LENGTH - NEW_KVEC(1)%POLY_ORDER
			N(2) = NEW_KVEC(2)%LENGTH - NEW_KVEC(2)%POLY_ORDER - 1
			
			JJ = 0
			DO J = 0, NEW_KVEC(2)%LENGTH
				IF (DABS(NEW_KVEC(2)%KNOTS(J) - NEW_KNOT).LE.EPS) THEN
					K = J
				ELSE
					JJ = JJ + 1
					OLD_KNOTS(JJ) = NEW_KVEC(2)%KNOTS(J)
				ENDIF
			ENDDO
			
			IF (JJ.NE.NEW_KVEC(2)%LENGTH) THEN
				PRINT*, 'ERROR - H_REFINEMENT (AXIS=Y) : NUMBER OF KNOTS IN OLD KNOT VECTOR IS NOT EQUAL TO THE NUMBER OF KNOTS IN NEW KNOT VECTOR - 1'
				STOP
			ENDIF
			
			DO I = 0, N(1) - 1
				DO J = 1, N(2) + 1
					IF (J.GE.1 .AND. J.LE.(K - NEW_KVEC(2)%POLY_ORDER)) THEN
					
						NEW_CTL%PTS(I, J-1, 0)%X = OLD_CTL%PTS(I, J-1, 0)%X
						NEW_CTL%PTS(I, J-1, 0)%Y = OLD_CTL%PTS(I, J-1, 0)%Y
						NEW_CTL%WGTS(I, J-1, 0) = OLD_CTL%WGTS(I, J-1, 0)
					
					ELSEIF (J.GE.(K - NEW_KVEC(2)%POLY_ORDER + 1) .AND. J.LE.K) THEN
						
						ALPHA = (NEW_KNOT - OLD_KNOTS(J))/(OLD_KNOTS(J + NEW_KVEC(2)%POLY_ORDER) - OLD_KNOTS(J))
						
						NEW_CTL%WGTS(I, J-1, 0) = ALPHA*OLD_CTL%WGTS(I, J-1, 0) + (1.D0 - ALPHA)*OLD_CTL%WGTS(I, J-2, 0)
						
						NEW_CTL%PTS(I, J-1, 0)%X = (ALPHA*OLD_CTL%WGTS(I, J-1, 0)/NEW_CTL%WGTS(I, J-1, 0))*OLD_CTL%PTS(I, J-1, 0)%X + &
																			 ((1.D0 - ALPHA)*OLD_CTL%WGTS(I, J-2, 0)/NEW_CTL%WGTS(I, J-1, 0))*OLD_CTL%PTS(I, J-2, 0)%X
						NEW_CTL%PTS(I, J-1, 0)%Y = (ALPHA*OLD_CTL%WGTS(I, J-1, 0)/NEW_CTL%WGTS(I, J-1, 0))*OLD_CTL%PTS(I, J-1, 0)%Y + &
																			 ((1.D0 - ALPHA)*OLD_CTL%WGTS(I, J-2, 0)/NEW_CTL%WGTS(I, J-1, 0))*OLD_CTL%PTS(I, J-2, 0)%Y
						
					ELSEIF (J.GE.(K + 1) .AND. J.LE.(N(2) + NEW_KVEC(2)%POLY_ORDER + 2)) THEN
						
						NEW_CTL%PTS(I, J-1, 0)%X = OLD_CTL%PTS(I, J-2, 0)%X
						NEW_CTL%PTS(I, J-1, 0)%Y = OLD_CTL%PTS(I, J-2, 0)%Y
						
						NEW_CTL%WGTS(I, J-1, 0) = OLD_CTL%WGTS(I, J-2, 0)
					
					ENDIF
				ENDDO
			ENDDO
		ENDIF 
		
		OLD_CTL = NEW_CTL
		
	END SUBROUTINE H_REFINEMENT
	
	SUBROUTINE REFINE_KNOT(OLD_KVEC, OLD_CTL, NEW_KNOTS, R, AXIS)
		
		TYPE(KNOT_VECTOR), INTENT(INOUT) :: OLD_KVEC(2)
		TYPE(CONTROL_POINTS_2D), INTENT(INOUT) :: OLD_CTL
		INTEGER, INTENT(IN) :: R
		REAL*8, INTENT(IN) :: NEW_KNOTS(0:R)
		CHARACTER(LEN=1), INTENT(IN) :: AXIS
		
		TYPE(KNOT_VECTOR) :: NEW_KVEC(2)
		TYPE(CONTROL_POINTS_2D) :: NEW_CTL
		INTEGER :: I, J, K, L, II, JJ, KK, N(2), M(2), P(2), A, B, IND
		REAL*8 :: ALPHA
		
! 		PRINT*, NEW_KNOTS(0:R)
		
		NEW_CTL%D = 2
		
		M(:) = OLD_KVEC(:)%LENGTH
		P(:) = OLD_KVEC(:)%POLY_ORDER
		N(1) = M(1) - P(1) - 1
		N(2) = M(2) - P(2) - 1
		
! 		PRINT*, 'M', M
! 		PRINT*, 'P', P
! 		PRINT*, 'N', N
		
		DO J = 0, N(2)
			DO I = 0, N(1)
				OLD_CTL%PTS(I, J, 0) = OLD_CTL%WGTS(I, J, 0)*OLD_CTL%PTS(I, J, 0)
			ENDDO
		ENDDO
		
		IF (AXIS=='X') THEN
			PRINT*, 'Refine knot vector along xi - direction'
			A = FIND_KNOT_SPAN(NEW_KNOTS(0),OLD_KVEC(1))
			B = FIND_KNOT_SPAN(NEW_KNOTS(R),OLD_KVEC(1))
			B = B + 1
			
! 			PRINT*, 'A', A
! 			PRINT*, 'B', B
			
			DO J = 0, N(2)
				DO I = 0, A - P(1)
					NEW_CTL%PTS(I, J, 0) = OLD_CTL%PTS(I, J, 0)
					NEW_CTL%WGTS(I, J, 0) = OLD_CTL%WGTS(I, J, 0)
				ENDDO
				DO I = B-1, N(1)
					NEW_CTL%PTS(I + R + 1, J, 0) = OLD_CTL%PTS(I, J, 0)
					NEW_CTL%WGTS(I + R + 1, J, 0) = OLD_CTL%WGTS(I, J, 0)
				ENDDO
			ENDDO
			
			DO I = 0, A
				NEW_KVEC(1)%KNOTS(I) = OLD_KVEC(1)%KNOTS(I)
			ENDDO
			
			DO I = B + P(1), M(1)
				NEW_KVEC(1)%KNOTS(I + R + 1) = OLD_KVEC(1)%KNOTS(I)
			ENDDO
			
			I = B + P(1) - 1
			K = B + P(1) + R
			
			DO J = R, 0, -1
				DO WHILE (NEW_KNOTS(J).LE.OLD_KVEC(1)%KNOTS(I) .AND. I.GT.A)
					DO JJ = 0, N(2)
						NEW_CTL%PTS(K - P(1) - 1, JJ, 0) = OLD_CTL%PTS(I - P(1) - 1, JJ, 0)
						NEW_CTL%WGTS(K - P(1) - 1, JJ, 0) = OLD_CTL%WGTS(I - P(1) - 1, JJ, 0)
					ENDDO
					NEW_KVEC(1)%KNOTS(K) = OLD_KVEC(1)%KNOTS(I)
					K = K - 1
					I = I - 1
				ENDDO
				
				DO JJ = 0, N(2)
					NEW_CTL%PTS(K - P(1) - 1, JJ, 0) = NEW_CTL%PTS(K - P(1), JJ, 0)
					NEW_CTL%WGTS(K - P(1) - 1, JJ, 0) = NEW_CTL%WGTS(K - P(1), JJ, 0)
				ENDDO
				
				DO L = 1, P(1)
					IND = K - P(1) + L
! 					PRINT*, '----------------------'
! 					PRINT*, 'L', L
! 					PRINT*, 'K', K
! 					PRINT*, 'IND', IND
					ALPHA = NEW_KVEC(1)%KNOTS(K+L) - NEW_KNOTS(J)
! 					PRINT*, 'L', L, 'ALPHA', ALPHA
					IF (DABS(ALPHA).LE.EPS) THEN
						DO JJ = 0, N(2)
							NEW_CTL%PTS(IND-1, JJ, 0) = NEW_CTL%PTS(IND, JJ, 0)
							NEW_CTL%WGTS(IND-1, JJ, 0) = NEW_CTL%WGTS(IND, JJ, 0)
						ENDDO
					ELSE
						ALPHA = ALPHA / (NEW_KVEC(1)%KNOTS(K+L) - OLD_KVEC(1)%KNOTS(I - P(1) + L))
! 						PRINT*, 'L', L, 'ALPHA', ALPHA
						DO JJ = 0, N(2)
							NEW_CTL%PTS(IND-1, JJ, 0)%X = ALPHA*NEW_CTL%PTS(IND-1, JJ, 0)%X + (1.D0 - ALPHA)*NEW_CTL%PTS(IND, JJ, 0)%X
							NEW_CTL%PTS(IND-1, JJ, 0)%Y = ALPHA*NEW_CTL%PTS(IND-1, JJ, 0)%Y + (1.D0 - ALPHA)*NEW_CTL%PTS(IND, JJ, 0)%Y
							NEW_CTL%WGTS(IND-1, JJ, 0) = ALPHA*NEW_CTL%WGTS(IND-1, JJ, 0) + (1.D0 - ALPHA)*NEW_CTL%WGTS(IND, JJ, 0)
						ENDDO
					ENDIF
				ENDDO
				NEW_KVEC(1)%KNOTS(K) = NEW_KNOTS(J)
				K = K - 1
			ENDDO
			
			NEW_KVEC(1)%LENGTH = OLD_KVEC(1)%LENGTH + UBOUND(NEW_KNOTS, 1) + 1
			NEW_KVEC(1)%POLY_ORDER = OLD_KVEC(1)%POLY_ORDER
			
			OLD_KVEC(1) = NEW_KVEC(1)
			
		ELSEIF (AXIS=='Y') THEN
			PRINT*, 'Refine knot vector along eta - direction'
			A = FIND_KNOT_SPAN(NEW_KNOTS(0),OLD_KVEC(2))
			B = FIND_KNOT_SPAN(NEW_KNOTS(R),OLD_KVEC(2))
			B = B + 1
			
			DO I = 0, N(1)
				DO J = 0, A - P(2)
					NEW_CTL%PTS(I, J, 0) = OLD_CTL%PTS(I, J, 0)
					NEW_CTL%WGTS(I, J, 0) = OLD_CTL%WGTS(I, J, 0)
				ENDDO
				DO J = B-1, N(2)
					NEW_CTL%PTS(I, J + R + 1, 0) = OLD_CTL%PTS(I, J, 0)
					NEW_CTL%WGTS(I, J + R + 1, 0) = OLD_CTL%WGTS(I, J, 0)
				ENDDO
			ENDDO
			
			DO I = 0, A
				NEW_KVEC(2)%KNOTS(I) = OLD_KVEC(2)%KNOTS(I)
			ENDDO
			
			DO I = B + P(2), M(2)
				NEW_KVEC(2)%KNOTS(I + R + 1) = OLD_KVEC(2)%KNOTS(I)
			ENDDO
			
			I = B + P(2) - 1
			K = B + P(2) + R
			
			DO J = R, 0, -1
				DO WHILE (NEW_KNOTS(J).LE.OLD_KVEC(2)%KNOTS(I) .AND. I.GT.A)
					DO II = 0, N(1)
						NEW_CTL%PTS(II, K - P(2) - 1, 0) = OLD_CTL%PTS(II, I - P(2) - 1, 0)
						NEW_CTL%WGTS(II, K - P(2) - 1, 0) = OLD_CTL%WGTS(II, I - P(2) - 1, 0)
					ENDDO
					NEW_KVEC(2)%KNOTS(K) = OLD_KVEC(2)%KNOTS(I)
					K = K - 1
					I = I - 1
				ENDDO
				
				DO II = 0, N(1)
					NEW_CTL%PTS(II, K - P(2) - 1, 0) = NEW_CTL%PTS(II, K - P(2), 0)
					NEW_CTL%WGTS(II, K - P(2) - 1, 0) = NEW_CTL%WGTS(II, K - P(2), 0)
				ENDDO
				
				DO L = 1, P(2)
					IND = K - P(2) + L
					ALPHA = NEW_KVEC(2)%KNOTS(K+L) - NEW_KNOTS(J)
! 					PRINT*, 'L', L, 'ALPHA', ALPHA
					IF (DABS(ALPHA).LE.EPS) THEN
						DO II = 0, N(1)
							NEW_CTL%PTS(II, IND-1, 0) = NEW_CTL%PTS(II, IND, 0)
							NEW_CTL%WGTS(II, IND-1, 0) = NEW_CTL%WGTS(II, IND, 0)
						ENDDO
					ELSE
						ALPHA = ALPHA / (NEW_KVEC(2)%KNOTS(K+L) - OLD_KVEC(2)%KNOTS(I - P(2) + L))
! 						PRINT*, 'L', L, 'ALPHA', ALPHA
						DO II = 0, N(1)
							NEW_CTL%PTS(II, IND-1, 0)%X = ALPHA*NEW_CTL%PTS(II, IND-1, 0)%X + (1.D0 - ALPHA)*NEW_CTL%PTS(II, IND, 0)%X
							NEW_CTL%PTS(II, IND-1, 0)%Y = ALPHA*NEW_CTL%PTS(II, IND-1, 0)%Y + (1.D0 - ALPHA)*NEW_CTL%PTS(II, IND, 0)%Y
							NEW_CTL%WGTS(II, IND-1, 0) = ALPHA*NEW_CTL%WGTS(II, IND-1, 0) + (1.D0 - ALPHA)*NEW_CTL%WGTS(II, IND, 0)
						ENDDO
					ENDIF
				ENDDO
				NEW_KVEC(2)%KNOTS(K) = NEW_KNOTS(J)
				K = K - 1
			ENDDO
			
			NEW_KVEC(2)%LENGTH = OLD_KVEC(2)%LENGTH + UBOUND(NEW_KNOTS, 1) + 1
			NEW_KVEC(2)%POLY_ORDER = OLD_KVEC(2)%POLY_ORDER
			
			OLD_KVEC(2) = NEW_KVEC(2)
		ENDIF
		
		DO J = 0, MAX_ORDER
			DO I = 0, MAX_ORDER
				NEW_CTL%PTS(I, J, 0)%X = NEW_CTL%PTS(I, J, 0)%X/NEW_CTL%WGTS(I, J, 0)
				NEW_CTL%PTS(I, J, 0)%Y = NEW_CTL%PTS(I, J, 0)%Y/NEW_CTL%WGTS(I, J, 0)
			ENDDO
		ENDDO
		
		OLD_CTL = NEW_CTL
		
	END SUBROUTINE REFINE_KNOT
	
	SUBROUTINE ELEVATE_DEGREE(OLD_KVEC, OLD_CTL, NP, AXIS)
		! NP = NEW POLY ORDER
		TYPE(KNOT_VECTOR), INTENT(INOUT) :: OLD_KVEC(2)
		TYPE(CONTROL_POINTS_2D), INTENT(INOUT) :: OLD_CTL
		INTEGER, INTENT(IN) :: NP
		CHARACTER(LEN=1), INTENT(IN) :: AXIS
		
		REAL*8, ALLOCATABLE :: BEZALFS(:,:), ALPHAS(:,:), BIN(:,:)
		TYPE(KNOT_VECTOR) :: NEW_KVEC(2)
		TYPE(CONTROL_POINTS_2D) :: BPTS, EBPTS, NEXTBPTS, NEW_CTL
		INTEGER :: I, J, K, II, JJ, KK, M(2), N(2), P(2), PH, PH2, BIN_MAX, MPI, MH, KND
		INTEGER :: R, A, B, CND, MUL, OLDR, LBZ, RBZ, SV, S, FIRST, LAST, TR, KJ
		REAL*8 :: UA, UB, INV, NUMER, DEN, BET, ALF, GAM
		
		BIN_MAX = 20
		
		ALLOCATE(BIN(0:BIN_MAX, 0:BIN_MAX))
		!! COMPUTE AND STORE BINOMIAL COEFFICIENTS
		DO I = 1, BIN_MAX
			DO J = 0, BIN_MAX
				BIN(I,J) = 1.0D0
			ENDDO
		ENDDO

		BIN(1,0) = 1.0D0
		BIN(1,1) = 1.0D0
		DO I = 2, BIN_MAX
			BIN(I,0) = BIN(I-1,0)
			DO J = 1, I-1
				BIN(I,J) = BIN(I-1,J-1) + BIN(I-1,J)
			ENDDO
			BIN(I,I) = BIN(I-1,I-1)
		ENDDO
	
		NEW_CTL%D = 2
		
		M(:) = OLD_KVEC(:)%LENGTH
		P(:) = OLD_KVEC(:)%POLY_ORDER
		N(1) = M(1) - P(1) - 1
		N(2) = M(2) - P(2) - 1
		
! 		PRINT*, 'M', M
! 		PRINT*, 'P', P
! 		PRINT*, 'N', N
		
		DO J = 0, N(2)
			DO I = 0, N(1)
				OLD_CTL%PTS(I, J, 0) = OLD_CTL%WGTS(I, J, 0)*OLD_CTL%PTS(I, J, 0)
			ENDDO
		ENDDO
		
		IF (AXIS=='X') THEN
			PRINT*, 'Elevate the degree of B-Spline along xi - direction'
			ALLOCATE(BEZALFS(0:P(1)+NP+1, 0:P(1)+1), ALPHAS(0:P(1)-1, 0:N(2)))
			PH = P(1) + NP
! 			PRINT*, 'PH', PH
			PH2 = INT(0.5*PH)
			!! Compute Bezier degree elevation coefficients
			BEZALFS(0,0) = 1.D0
			BEZALFS(PH,P(1)) = 1.D0
			DO I = 1, PH2
				INV = 1.D0/BIN(PH,I)
				MPI = MIN(P(1), I)
				DO J = MAX(0, I-NP), MPI
					BEZALFS(I,J) = INV*BIN(P(1),J)*BIN(NP,I-J)
				ENDDO
			ENDDO
			DO I = PH2+1, PH-1
				MPI = MIN(P(1),I)
				DO J = MAX(0,I-NP), MPI
					BEZALFS(I,J) = BEZALFS(PH-I,P(1)-J)
				ENDDO
			ENDDO
			MH = PH; KND = PH+1
			R = -1;  A = P(1)
			B = P(1) + 1; CND = 1
			UA = OLD_KVEC(1)%KNOTS(0)
! 			PRINT*, 'HERE-1'
			DO JJ = 0, N(2)
				NEW_CTL%PTS(0, JJ, 0) = OLD_CTL%PTS(0, JJ, 0)
				NEW_CTL%WGTS(0, JJ, 0) = OLD_CTL%WGTS(0, JJ, 0)
! 				PRINT*, 'JJ',JJ,'WEIGHT',NEW_CTL%WGTS(0,JJ,0)
			ENDDO
			DO I = 0, PH
				NEW_KVEC(1)%KNOTS(I) = UA
			ENDDO
				!! Initialize first Bezier seg
			DO I = 0, P(1)
				DO JJ = 0, N(2)
					BPTS%PTS(I,JJ,0) = OLD_CTL%PTS(I,JJ,0)
					BPTS%WGTS(I,JJ,0) = OLD_CTL%WGTS(I,JJ,0)
				ENDDO
			ENDDO
			!! Big loop thru knot vector
			DO WHILE (B<M(1))
				I = B
				DO WHILE (B<M(1) .AND. OLD_KVEC(1)%KNOTS(B)==OLD_KVEC(1)%KNOTS(B+1))
					B = B + 1
				ENDDO
				MUL = B - I + 1
				MH = MH + MUL + NP
				UB = OLD_KVEC(1)%KNOTS(B)
				OLDR = R; R = P(1) - MUL
				!! Insert knot u(b) r times
				IF (OLDR.GT.0) THEN
					LBZ = INT(0.5*(OLDR + 2))
				ELSE
					LBZ = 1
				ENDIF
				IF (R.GT.0) THEN
					RBZ = PH - INT(0.5*(R+1))
				ELSE
					RBZ = PH
				ENDIF
				
				IF (R.GT.0) THEN
				!! Insert knot to get Bezier segment
					NUMER = UB - UA
					DO K = P(1), MUL + 1, -1
						DO JJ = 0, N(2)
							ALPHAS(K-MUL-1, JJ) = NUMER/(OLD_KVEC(1)%KNOTS(A+K) - UA)
						ENDDO
					ENDDO
					DO J = 1, R
						SV = R - J; S = MUL + J
						DO K = P(1), S, -1
							DO JJ = 0, N(2)
! 								PRINT*, 'K', K, 'JJ', JJ
								BPTS%PTS(K, JJ, 0)%X = ALPHAS(K-S, JJ)*BPTS%PTS(K, JJ, 0)%X + (1.D0 - ALPHAS(K-S, JJ))*BPTS%PTS(K-1, JJ, 0)%X
								BPTS%PTS(K, JJ, 0)%Y = ALPHAS(K-S, JJ)*BPTS%PTS(K, JJ, 0)%Y + (1.D0 - ALPHAS(K-S, JJ))*BPTS%PTS(K-1, JJ, 0)%Y
								BPTS%WGTS(K, JJ, 0) = ALPHAS(K-S, JJ)*BPTS%WGTS(K, JJ, 0) + (1.D0 - ALPHAS(K-S, JJ))*BPTS%WGTS(K-1, JJ, 0)
							ENDDO
						ENDDO
						NEXTBPTS%PTS(SV, :, 0) = BPTS%PTS(P(1), :, 0)
						NEXTBPTS%WGTS(SV, :, 0) = BPTS%WGTS(P(1), :, 0)
					ENDDO
				ENDIF !! End of "insert knot"
				
				!! Degree elevate Bezier
				DO I = LBZ, PH	!! Only points lbz, ... , ph are used below
					EBPTS%PTS(I, :, 0) = POINT2D(0.D0, 0.D0)
					EBPTS%WGTS(I, :, 0) = 0.D0
					MPI = MIN(P(1),I)
					DO J = MAX(0,I-NP), MPI
						DO JJ = 0, N(2)
							EBPTS%PTS(I, JJ, 0)%X = EBPTS%PTS(I, JJ, 0)%X + BEZALFS(I,J)*BPTS%PTS(J, JJ, 0)%X
							EBPTS%PTS(I, JJ, 0)%Y = EBPTS%PTS(I, JJ, 0)%Y + BEZALFS(I,J)*BPTS%PTS(J, JJ, 0)%Y
							EBPTS%WGTS(I, JJ, 0) = EBPTS%WGTS(I, JJ, 0) + BEZALFS(I,J)*BPTS%WGTS(J, JJ, 0)
						ENDDO
					ENDDO
				ENDDO	!! End of degree elevating Bezier
				
				IF (OLDR > 1) THEN
					!! Must remove knot u=U[a] oldr times
					FIRST = KND - 2; LAST = KND
					DEN = UB - UA
					BET = (UB - NEW_KVEC(1)%KNOTS(KND-1))/DEN
! 					PRINT*, 'HERE-2'
					DO TR = 1, OLDR - 1
						!!	Knot removal loop
						I = FIRST; J = LAST; KJ = J - KND + 1
						DO WHILE (J-I > TR)	!! Loop and compute the new control points for one removal step
							IF (I < CND) THEN
								ALF = (UB - NEW_KVEC(1)%KNOTS(I)) / (UA - NEW_KVEC(1)%KNOTS(I))
								DO JJ = 0, N(2)
									NEW_CTL%PTS(I, JJ, 0)%X = ALF*NEW_CTL%PTS(I, JJ, 0)%X + (1.D0 - ALF)*NEW_CTL%PTS(I-1, JJ, 0)%X
									NEW_CTL%PTS(I, JJ, 0)%Y = ALF*NEW_CTL%PTS(I, JJ, 0)%Y + (1.D0 - ALF)*NEW_CTL%PTS(I-1, JJ, 0)%Y
									NEW_CTL%WGTS(I, JJ, 0) = ALF*NEW_CTL%WGTS(I, JJ, 0) + (1.D0 - ALF)*NEW_CTL%WGTS(I-1, JJ, 0)
! 									PRINT*, 'TR',TR,'J-I',J-I,'JJ',JJ,'WEIGHT',NEW_CTL%WGTS(I,JJ,0)
								ENDDO
							ENDIF
							IF (J.GE.LBZ) THEN
								IF ((J-TR).LE.(KND - PH + OLDR)) THEN
									GAM = (UB - NEW_KVEC(1)%KNOTS(J-TR))/DEN
									DO JJ = 0, N(2)
										EBPTS%PTS(KJ, JJ, 0)%X = GAM*EBPTS%PTS(KJ, JJ, 0)%X + (1.D0 - GAM)*EBPTS%PTS(KJ+1, JJ, 0)%X
										EBPTS%PTS(KJ, JJ, 0)%Y = GAM*EBPTS%PTS(KJ, JJ, 0)%Y + (1.D0 - GAM)*EBPTS%PTS(KJ+1, JJ, 0)%Y
										EBPTS%WGTS(KJ, JJ, 0) = GAM*EBPTS%WGTS(KJ, JJ, 0) + (1.D0 - GAM)*EBPTS%WGTS(KJ+1, JJ, 0)
									ENDDO
								ELSE
									DO JJ = 0, N(2)
										EBPTS%PTS(KJ, JJ, 0)%X = BET*EBPTS%PTS(KJ, JJ, 0)%X + (1.D0 - BET)*EBPTS%PTS(KJ+1, JJ, 0)%X
										EBPTS%PTS(KJ, JJ, 0)%Y = BET*EBPTS%PTS(KJ, JJ, 0)%Y + (1.D0 - BET)*EBPTS%PTS(KJ+1, JJ, 0)%Y
										EBPTS%WGTS(KJ, JJ, 0) = BET*EBPTS%WGTS(KJ, JJ, 0) + (1.D0 - BET)*EBPTS%WGTS(KJ+1, JJ, 0)
									ENDDO
								ENDIF
							ENDIF
							I = I + 1; J = J - 1; KJ = KJ - 1
						ENDDO
						FIRST = FIRST - 1; LAST = LAST + 1
					ENDDO
				ENDIF	!! End of removing knot, u=U[a]
				
				IF (A.NE.P(1)) THEN	!! Load the knot ua
					DO I = 0, PH - OLDR - 1
						NEW_KVEC(1)%KNOTS(KND) = UA; KND = KND + 1
					ENDDO
				ENDIF
! 				PRINT*, 'HERE-3'
				DO J = LBZ, RBZ	!! Load control points into Qw (NEW_CTL)
					DO JJ = 0, N(2)
						NEW_CTL%PTS(CND, JJ, 0)%X = EBPTS%PTS(J, JJ, 0)%X
						NEW_CTL%PTS(CND, JJ, 0)%Y = EBPTS%PTS(J, JJ, 0)%Y
						NEW_CTL%WGTS(CND, JJ, 0) = EBPTS%WGTS(J, JJ, 0)
! 						PRINT*, 'J',J,'CND',CND,'JJ',JJ,'WEIGHT',NEW_CTL%WGTS(CND,JJ,0)
					ENDDO
					CND = CND + 1
				ENDDO
				IF (B.LT.M(1)) THEN
					!! Set up for next pass thru loop
					DO J = 0, R - 1
						DO JJ = 0, N(2)
							BPTS%PTS(J, JJ, 0)%X = NEXTBPTS%PTS(J, JJ, 0)%X
							BPTS%PTS(J, JJ, 0)%Y = NEXTBPTS%PTS(J, JJ, 0)%Y
							BPTS%WGTS(J, JJ, 0) = NEXTBPTS%WGTS(J, JJ, 0)
						ENDDO
					ENDDO
					DO J = R, P(1)
						DO JJ = 0, N(2)
							BPTS%PTS(J, JJ, 0)%X = OLD_CTL%PTS(B - P(1) + J, JJ, 0)%X
							BPTS%PTS(J, JJ, 0)%Y = OLD_CTL%PTS(B - P(1) + J, JJ, 0)%Y
							BPTS%WGTS(J, JJ, 0) = OLD_CTL%WGTS(B - P(1) + J, JJ, 0)
						ENDDO
					ENDDO
					A = B; B = B + 1; UA = UB
				ELSE
					DO I = 0, PH
						NEW_KVEC(1)%KNOTS(KND + I) = UB
					ENDDO
				ENDIF	!! End knot
			ENDDO
			NEW_KVEC(1)%POLY_ORDER = PH
			NEW_KVEC(1)%LENGTH = MH
			OLD_KVEC(1) = NEW_KVEC(1)
! 			PRINT*, 'MH',MH,'PH',PH
			
		ELSEIF (AXIS=='Y') THEN
			PRINT*, 'Elevate the degree of B-Spline along eta - direction'
			ALLOCATE(BEZALFS(0:P(2)+NP+1, 0:P(2)+1), ALPHAS(0:P(2)-1, 0:N(1)))
			PH = P(2) + NP
! 			PRINT*, 'PH', PH
			PH2 = INT(0.5*PH)
			!! Compute Bezier degree elevation coefficients
			BEZALFS(0,0) = 1.D0
			BEZALFS(PH,P(2)) = 1.D0

			DO I = 1, PH2
				INV = 1.D0/BIN(PH,I)
				MPI = MIN(P(2), I)
				DO J = MAX(0, I-NP), MPI
					BEZALFS(I,J) = INV*BIN(P(2),J)*BIN(NP,I-J)
				ENDDO
			ENDDO

			DO I = PH2+1, PH-1
				MPI = MIN(P(2),I)
				DO J = MAX(0,I-NP), MPI
					BEZALFS(I,J) = BEZALFS(PH-I,P(2)-J)
				ENDDO
			ENDDO

			MH = PH; KND = PH+1
			R = -1;  A = P(2)
			B = P(2) + 1; CND = 1
			UA = OLD_KVEC(2)%KNOTS(0)
! 			PRINT*, 'HERE-1'
			DO II = 0, N(1)
				NEW_CTL%PTS(II, 0, 0) = OLD_CTL%PTS(II, 0, 0)
				NEW_CTL%WGTS(II, 0, 0) = OLD_CTL%WGTS(II, 0, 0)
! 				PRINT*, 'II',II,'WEIGHT',NEW_CTL%WGTS(0,JJ,0)
			ENDDO
			
			DO I = 0, PH
				NEW_KVEC(2)%KNOTS(I) = UA
			ENDDO
			
			!! Initialize first Bezier seg
			DO I = 0, P(2)
				DO II = 0, N(1)
					BPTS%PTS(II,I,0) = OLD_CTL%PTS(II,I,0)
					BPTS%WGTS(II,I,0) = OLD_CTL%WGTS(II,I,0)
				ENDDO
			ENDDO
			!! Big loop thru knot vector
			DO WHILE (B<M(2))
				I = B
				DO WHILE (B<M(2) .AND. OLD_KVEC(2)%KNOTS(B)==OLD_KVEC(2)%KNOTS(B+1))
					B = B + 1
				ENDDO
				MUL = B - I + 1
				MH = MH + MUL + NP
				UB = OLD_KVEC(2)%KNOTS(B)
				OLDR = R; R = P(2) - MUL
				!! Insert knot u(b) r times
				IF (OLDR.GT.0) THEN
					LBZ = INT(0.5*(OLDR + 2))
				ELSE
					LBZ = 1
				ENDIF
				IF (R.GT.0) THEN
					RBZ = PH - INT(0.5*(R+1))
				ELSE
					RBZ = PH
				ENDIF
				
				IF (R.GT.0) THEN
				!! Insert knot to get Bezier segment
					NUMER = UB - UA
					DO K = P(2), MUL + 1, -1
						DO II = 0, N(1)
							ALPHAS(K-MUL-1, II) = NUMER/(OLD_KVEC(2)%KNOTS(A+K) - UA)
						ENDDO
					ENDDO
					DO J = 1, R
						SV = R - J; S = MUL + J
						DO K = P(2), S, -1
							DO II = 0, N(1)
								BPTS%PTS(II, K, 0)%X = ALPHAS(K-S, II)*BPTS%PTS(II, K, 0)%X + (1.D0 - ALPHAS(K-S, II))*BPTS%PTS(II, K-1, 0)%X
								BPTS%PTS(II, K, 0)%Y = ALPHAS(K-S, II)*BPTS%PTS(II, K, 0)%Y + (1.D0 - ALPHAS(K-S, II))*BPTS%PTS(II, K-1, 0)%Y
								BPTS%WGTS(II, K, 0) = ALPHAS(K-S, II)*BPTS%WGTS(II, K, 0) + (1.D0 - ALPHAS(K-S, II))*BPTS%WGTS(II, K-1, 0)
							ENDDO
						ENDDO
						NEXTBPTS%PTS(:, SV, 0) = BPTS%PTS(:, P(2), 0)
						NEXTBPTS%WGTS(:, SV, 0) = BPTS%WGTS(:, P(2), 0)
					ENDDO
				ENDIF !! End of "insert knot"
				
				!! Degree elevate Bezier
				DO I = LBZ, PH	!! Only points lbz, ... , ph are used below
					EBPTS%PTS(:, I, 0) = POINT2D(0.D0, 0.D0)
					EBPTS%WGTS(:, I, 0) = 0.D0
					MPI = MIN(P(2),I)
					DO J = MAX(0,I-NP), MPI
						DO II = 0, N(1)
							EBPTS%PTS(II, I, 0)%X = EBPTS%PTS(II, I, 0)%X + BEZALFS(I,J)*BPTS%PTS(II, J, 0)%X
							EBPTS%PTS(II, I, 0)%Y = EBPTS%PTS(II, I, 0)%Y + BEZALFS(I,J)*BPTS%PTS(II, J, 0)%Y
							EBPTS%WGTS(II, I, 0) = EBPTS%WGTS(II, I, 0) + BEZALFS(I,J)*BPTS%WGTS(II, J, 0)
						ENDDO
					ENDDO
				ENDDO	!! End of degree elevating Bezier
				
				IF (OLDR > 1) THEN
					!! Must remove knot u=U[a] oldr times
					FIRST = KND - 2; LAST = KND
					DEN = UB - UA
					BET = (UB - NEW_KVEC(2)%KNOTS(KND-1))/DEN
! 					PRINT*, 'HERE-2'
					DO TR = 1, OLDR - 1
						!!	Knot removal loop
						I = FIRST; J = LAST; KJ = J - KND + 1
						DO WHILE (J-I > TR)	!! Loop and compute the new control points for one removal step
							IF (I < CND) THEN
								ALF = (UB - NEW_KVEC(2)%KNOTS(I)) / (UA - NEW_KVEC(2)%KNOTS(I))
								DO II = 0, N(1)
									NEW_CTL%PTS(II, I, 0)%X = ALF*NEW_CTL%PTS(II, I, 0)%X + (1.D0 - ALF)*NEW_CTL%PTS(II, I-1, 0)%X
									NEW_CTL%PTS(II, I, 0)%Y = ALF*NEW_CTL%PTS(II, I, 0)%Y + (1.D0 - ALF)*NEW_CTL%PTS(II, I-1, 0)%Y
									NEW_CTL%WGTS(II, I, 0) = ALF*NEW_CTL%WGTS(II, I, 0) + (1.D0 - ALF)*NEW_CTL%WGTS(II, I-1, 0)
! 									PRINT*, 'TR',TR,'J-I',J-I,'II',II,'WEIGHT',NEW_CTL%WGTS(II,I,0)
								ENDDO
							ENDIF
							IF (J.GE.LBZ) THEN
								IF ((J-TR).LE.(KND - PH + OLDR)) THEN
									GAM = (UB - NEW_KVEC(2)%KNOTS(J-TR))/DEN
									DO II = 0, N(1)
										EBPTS%PTS(II, KJ, 0)%X = GAM*EBPTS%PTS(II, KJ, 0)%X + (1.D0 - GAM)*EBPTS%PTS(II, KJ+1, 0)%X
										EBPTS%PTS(II, KJ, 0)%Y = GAM*EBPTS%PTS(II, KJ, 0)%Y + (1.D0 - GAM)*EBPTS%PTS(II, KJ+1, 0)%Y
										EBPTS%WGTS(II, KJ, 0) = GAM*EBPTS%WGTS(II, KJ, 0) + (1.D0 - GAM)*EBPTS%WGTS(II, KJ+1, 0)
									ENDDO
								ELSE
									DO II = 0, N(1)
										EBPTS%PTS(II, KJ, 0)%X = BET*EBPTS%PTS(II, KJ, 0)%X + (1.D0 - BET)*EBPTS%PTS(II, KJ+1, 0)%X
										EBPTS%PTS(II, KJ, 0)%Y = BET*EBPTS%PTS(II, KJ, 0)%Y + (1.D0 - BET)*EBPTS%PTS(II, KJ+1, 0)%Y
										EBPTS%WGTS(II, KJ, 0) = BET*EBPTS%WGTS(II, KJ, 0) + (1.D0 - BET)*EBPTS%WGTS(II, KJ+1, 0)
									ENDDO
								ENDIF
							ENDIF
							I = I + 1; J = J - 1; KJ = KJ - 1
						ENDDO
						FIRST = FIRST - 1; LAST = LAST + 1
					ENDDO
				ENDIF	!! End of removing knot, u=U[a]
				
				IF (A.NE.P(2)) THEN	!! Load the knot ua
					DO I = 0, PH - OLDR - 1
						NEW_KVEC(2)%KNOTS(KND) = UA; KND = KND + 1
					ENDDO
				ENDIF
! 				PRINT*, 'HERE-3'
				DO J = LBZ, RBZ	!! Load control points into Qw (NEW_CTL)
					DO II = 0, N(1)
						NEW_CTL%PTS(II, CND, 0)%X = EBPTS%PTS(II, J, 0)%X
						NEW_CTL%PTS(II, CND, 0)%Y = EBPTS%PTS(II, J, 0)%Y
						NEW_CTL%WGTS(II, CND, 0) = EBPTS%WGTS(II, J, 0)
! 						PRINT*, 'J',J,'CND',CND,'II',II,'WEIGHT',NEW_CTL%WGTS(II,CND,0)
					ENDDO
					CND = CND + 1
				ENDDO
				IF (B.LT.M(2)) THEN
					!! Set up for next pass thru loop
					DO J = 0, R - 1
						DO II = 0, N(1)
							BPTS%PTS(II, J, 0)%X = NEXTBPTS%PTS(II, J, 0)%X
							BPTS%PTS(II, J, 0)%Y = NEXTBPTS%PTS(II, J, 0)%Y
							BPTS%WGTS(II, J, 0) = NEXTBPTS%WGTS(II, J, 0)
						ENDDO
					ENDDO
					DO J = R, P(2)
						DO II = 0, N(1)
							BPTS%PTS(II, J, 0)%X = OLD_CTL%PTS(II, B - P(2) + J, 0)%X
							BPTS%PTS(II, J, 0)%Y = OLD_CTL%PTS(II, B - P(2) + J, 0)%Y
							BPTS%WGTS(II, J, 0) = OLD_CTL%WGTS(II, B - P(2) + J, 0)
						ENDDO
					ENDDO
					A = B; B = B + 1; UA = UB
				ELSE
					DO I = 0, PH
						NEW_KVEC(2)%KNOTS(KND + I) = UB
					ENDDO
				ENDIF	!! End knot
			ENDDO
			NEW_KVEC(2)%POLY_ORDER = PH
			NEW_KVEC(2)%LENGTH = MH
			OLD_KVEC(2) = NEW_KVEC(2)
! 			PRINT*, 'MH',MH,'PH',PH
		ENDIF
		
		DO J = 0, MAX_ORDER
			DO I = 0, MAX_ORDER
				NEW_CTL%PTS(I, J, 0)%X = NEW_CTL%PTS(I, J, 0)%X/NEW_CTL%WGTS(I, J, 0)
				NEW_CTL%PTS(I, J, 0)%Y = NEW_CTL%PTS(I, J, 0)%Y/NEW_CTL%WGTS(I, J, 0)
			ENDDO
		ENDDO
		
		OLD_CTL = NEW_CTL
		
		DEALLOCATE(BEZALFS, ALPHAS)
		
	END SUBROUTINE ELEVATE_DEGREE

END MODULE KNOT_HANDLING