MODULE GEOMETRY

	USE GSQUAD
	USE KNOT_HANDLING
	USE PATCH_MAPPING
	USE NURBS
	USE BESSEL

! 	IMPLICIT INTEGER (I-N)
! 	IMPLICIT REAL(8) (A-H,O-Z)

CONTAINS

SUBROUTINE GET_GEO()
	
	INTEGER :: I, J, II, JJ, K, KK

! Quarter - circle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	GEO_KVEC(1)%POLY_ORDER = 2
! 	GEO_KVEC(2)%POLY_ORDER = 2
! 
! 	GEO_KVEC(1)%LENGTH = 5
! 	GEO_KVEC(2)%LENGTH = 5
! 
! 	GEO_KVEC(1)%KNOTS(0:GEO_KVEC(1)%LENGTH) = (/0.D0, 0.D0, 0.D0, 1.D0, 1.D0, 1.D0/)
! 	GEO_KVEC(2)%KNOTS(0:GEO_KVEC(2)%LENGTH) = (/0.D0, 0.D0, 0.D0, 1.D0, 1.D0, 1.D0/)
! 	
! ! SET CONTROL POINTS AND WEIGTH
! 	
! 	GEO_CTL%D=2
! 	
! 	GEO_CTL%PTS(0,0,0) = POINT2D(0.D0,0.D0)
! 	GEO_CTL%PTS(1,0,0) = POINT2D(0.5D0,0.D0)
! 	GEO_CTL%PTS(2,0,0) = POINT2D(1.D0,0.D0)
! 	
! 	GEO_CTL%PTS(0,1,0) = POINT2D(0.D0,0.5D0)
! 	GEO_CTL%PTS(1,1,0) = POINT2D(0.5D0*DCOS(45.D0*DEGREE),0.5D0*DSIN(45.D0*DEGREE))
! 	GEO_CTL%PTS(2,1,0) = POINT2D(1.D0,-1.D0 + DCOS(45.D0*DEGREE) + DSIN(45.D0*DEGREE))
! 	
! 	GEO_CTL%PTS(0,2,0) = POINT2D(0.D0,1.D0)
! 	GEO_CTL%PTS(1,2,0) = POINT2D(GEO_CTL%PTS(2,1,0)%Y, GEO_CTL%PTS(2,1,0)%X)
! 	GEO_CTL%PTS(2,2,0) = POINT2D(DCOS(45.D0*DEGREE), DSIN(45.D0*DEGREE))
! 	
! 	GEO_CTL%WGTS(0,0,0) = 1.D0
! 	GEO_CTL%WGTS(1,0,0) = 1.D0
! 	GEO_CTL%WGTS(2,0,0) = 1.D0
! 	
! 	GEO_CTL%WGTS(0,1,0) = 1.D0
! 	GEO_CTL%WGTS(1,1,0) = 1.D0
! 	GEO_CTL%WGTS(2,1,0) = DCOS(22.5D0*DEGREE)
! 	
! 	GEO_CTL%WGTS(0,2,0) = 1.D0
! 	GEO_CTL%WGTS(1,2,0) = DCOS(22.5D0*DEGREE)
! 	GEO_CTL%WGTS(2,2,0) = 1.D0
! 	
! 	BASIS_CTL%WGTS(:,:,:) = 1.D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! UNIT CIRCLE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	GEO_KVEC(1)%POLY_ORDER = 2
	GEO_KVEC(2)%POLY_ORDER = 2

	GEO_KVEC(1)%LENGTH = 5
	GEO_KVEC(2)%LENGTH = 7

	GEO_KVEC(1)%KNOTS(0:GEO_KVEC(1)%LENGTH) = (/0.D0, 0.D0, 0.D0, 1.D0, 1.D0, 1.D0/)
	GEO_KVEC(2)%KNOTS(0:GEO_KVEC(2)%LENGTH) = (/0.D0, 0.D0, 0.D0, 0.5D0, 0.5D0, 1.D0, 1.D0, 1.D0/)
	
! SET CONTROL POINTS AND WEIGTH
	
	GEO_CTL%D=2
	
	GEO_CTL%PTS(0,0,0) = POINT2D(0.D0,-1.D0)
	GEO_CTL%PTS(1,0,0) = POINT2D(1.D0,-1.D0)
	GEO_CTL%PTS(2,0,0) = POINT2D(1.D0,0.D0)
	
	GEO_CTL%PTS(0,1,0) = POINT2D(-1.D0,-1.D0)
	GEO_CTL%PTS(1,1,0) = POINT2D(0.5D0,-0.5D0)
	GEO_CTL%PTS(2,1,0) = POINT2D(0.5D0,0.D0)
	
	GEO_CTL%PTS(0,2,0) = POINT2D(-1.D0,0.D0)
	GEO_CTL%PTS(1,2,0) = POINT2D(-0.5D0,0.D0)
	GEO_CTL%PTS(2,2,0) = POINT2D(0.D0,0.D0)
	
	GEO_CTL%PTS(0,3,0) = POINT2D(-1.D0,1.D0)
	GEO_CTL%PTS(1,3,0) = POINT2D(0.5D0,0.5D0)
	GEO_CTL%PTS(2,3,0) = POINT2D(0.5D0,0.D0)
	
	GEO_CTL%PTS(0,4,0) = POINT2D(0.D0,1.D0)
	GEO_CTL%PTS(1,4,0) = POINT2D(1.D0,1.D0)
	GEO_CTL%PTS(2,4,0) = POINT2D(1.D0,0.D0)
	
	GEO_CTL%WGTS(0,0,0) = 1.D0
	GEO_CTL%WGTS(1,0,0) = DCOS(45.D0*DEGREE)
	GEO_CTL%WGTS(2,0,0) = 1.D0
	
	GEO_CTL%WGTS(0,1,0) = DCOS(45.D0*DEGREE)
	GEO_CTL%WGTS(1,1,0) = DCOS(45.D0*DEGREE)
	GEO_CTL%WGTS(2,1,0) = 1.D0
	
	GEO_CTL%WGTS(0,2,0) = 1.D0
	GEO_CTL%WGTS(1,2,0) = 1.D0
	GEO_CTL%WGTS(2,2,0) = 1.D0
	
	GEO_CTL%WGTS(0,3,0) = DCOS(45.D0*DEGREE)
	GEO_CTL%WGTS(1,3,0) = DCOS(45.D0*DEGREE)
	GEO_CTL%WGTS(2,3,0) = 1.D0
	
	GEO_CTL%WGTS(0,4,0) = 1.D0
	GEO_CTL%WGTS(1,4,0) = DCOS(45.D0*DEGREE)
	GEO_CTL%WGTS(2,4,0) = 1.D0
	
	BASIS_CTL%WGTS(:,:,:) = 1.D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	WRITE(*,*)
	WRITE(*,*) '<<< SET GEOMETRIC KNOT, CONTROL POINTS, AND WEIGHT : DONE >>>'
	WRITE(*,*)
	
! SET KNOT INFO FOR B-SPLINE BASIS FUNCTIONS
	BASIS_KVEC(:) = GEO_KVEC(:)

	IF (BS_ORDER(1).GT.GEO_KVEC(1)%POLY_ORDER) THEN
		DO I=1, BS_ORDER(1) - GEO_KVEC(1)%POLY_ORDER
			BASIS_KVEC(1) = DEGREE_ELEVATION(BASIS_KVEC(1))
		ENDDO
	ENDIF
	IF (BS_ORDER(2).GT.GEO_KVEC(2)%POLY_ORDER) THEN
		DO I=1, BS_ORDER(2) - GEO_KVEC(2)%POLY_ORDER
			BASIS_KVEC(2) = DEGREE_ELEVATION(BASIS_KVEC(2))
		ENDDO
	ENDIF
	
	IF (EXTRA_KNOTS(1)>0) THEN
		BASIS_KVEC(1) = UNIFORM_KNOT_INSERTION(BASIS_KVEC(1),EXTRA_KNOTS(1))
	ENDIF
	IF (EXTRA_KNOTS(2)>0) THEN
		BASIS_KVEC(2) = UNIFORM_KNOT_INSERTION(BASIS_KVEC(2),EXTRA_KNOTS(2))
	ENDIF

	WRITE(*,*)
	WRITE(*,*) '<<< SET BASIS KNOT : DONE >>>'
	WRITE(*,*)

	NUMBS(:) = (/BASIS_KVEC(1)%LENGTH - BS_ORDER(1), BASIS_KVEC(2)%LENGTH - BS_ORDER(2)/)
	DOF = NUMBS(1)*NUMBS(2)
	
! GLOBAL INDEX
	ALLOCATE(NDX(2,DOF))
	DO I=1, DOF
		NDX(1,I)=MOD(I+(NUMBS(1)-1), NUMBS(1))
		NDX(2,I)=MOD(INT((I-1)/NUMBS(1))+NUMBS(2), NUMBS(2))
	ENDDO

	WRITE(*,*)
	WRITE(*,*) '<<< SET GLOBAL INDEX : DONE >>>'
	WRITE(*,*)

! CONSTRUCT THE FRAME OF LOCAL INTEGRAL REGIONS
	IR_GRID(1,1) = 0.D0
	IR_GRID(2,1) = 0.D0
	
! 	PRINT*, BASIS_KVEC(1)%LENGTH, BASIS_KVEC(2)%LENGTH
	K = 1
	DO I=1, BASIS_KVEC(1)%LENGTH
		K = K + 1
		IF (DABS(BASIS_KVEC(1)%KNOTS(I-1) - BASIS_KVEC(1)%KNOTS(I)).LE.EPS) THEN
			K = K - 1
		ELSE
			IR_GRID(1,K) = BASIS_KVEC(1)%KNOTS(I)
		ENDIF
! 		IF (K.GE.2) THEN
! 			DO J=1, 
! 		ENDIF
	ENDDO
	NUMIR(1) = K
	
	K = 1
	DO I=1, BASIS_KVEC(2)%LENGTH
		K = K + 1
		IF (DABS(BASIS_KVEC(2)%KNOTS(I-1) - BASIS_KVEC(2)%KNOTS(I)).LE.EPS) THEN
			K = K - 1
		ELSE
			IR_GRID(2,K) = BASIS_KVEC(2)%KNOTS(I)
		ENDIF
	ENDDO
	NUMIR(2) = K

	WRITE(*,*)
	WRITE(*,*) '<<< SET INTEGRAL REGION : DONE >>>'
	WRITE(*,*)

! SET INDEX OF BASIS FUNCTIONS IMPOSED ESSENTIAL BOUNDARY CONDITION
	ALLOCATE(BD_MASK(DOF))
	BD_MASK(:) = .FALSE.
	BD_DOF = 0
	
	K = 0
	DO I=1, DOF
		IF (NDX(1,I)==0 .OR. NDX(1,I)==NUMBS(1)-1 .OR. NDX(2,I)==NUMBS(2)-1 .OR. NDX(2,I)==0) THEN
			K = K + 1
			BDNDX(K)%LC_NDX(:) = NDX(:,I)
			BDNDX(K)%GL_NDX = I
			BD_MASK(I)=.TRUE.
			BD_DOF = BD_DOF + 1
		ENDIF
	ENDDO
	
	BD_DOF = BD_DOF
	BDNDX(:)%LC_NUM = K
		
	WRITE(*,*)
	WRITE(*,*) '<<< SET LOCAL INDEX CORRESPONDING BASIS FUNCTIONS ON BOUNDARY IMPOSED ESSENTIAL BOUNDARY CONDITION : DONE >>>'
	WRITE(*,*)

	!!  GENERATE BINOMIAL COEFFICIENTS

	DO I = 1, MAX_DIFF_ORDER
		DO J = 0, MAX_DIFF_ORDER
			BINOM(I,J) = 1.0D0
		ENDDO
	ENDDO

	BINOM(1,0) = 1.0D0
	BINOM(1,1) = 1.0D0
	DO I = 2, MAX_DIFF_ORDER
		BINOM(I,0) = BINOM(I-1,0)
		DO J = 1, I-1
			BINOM(I,J) = BINOM(I-1,J-1) + BINOM(I-1,J)
		ENDDO
		BINOM(I,I) = BINOM(I-1,I-1)
	ENDDO

END SUBROUTINE GET_GEO

END MODULE GEOMETRY

