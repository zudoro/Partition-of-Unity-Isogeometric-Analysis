MODULE LOADFUNCTION

	USE GLBVAR
	USE PATCH_MAPPING

CONTAINS

! REAL*8 FUNCTION EX_DISP(PT2D)
! 
! 	TYPE(POINT2D), INTENT(IN) :: PT2D
! 	TYPE(POINT2D) :: PHY_PT
! 	REAL*8 :: R, THETA, X, Y
! 
! 	PHY_PT = GET_PHY_PT(PT2D)
! 	
! 	X = PHY_PT%X; Y = PHY_PT%Y
! 	
! 	CALL GET_RTHETA(R, THETA, PHY_PT)
! 
! 	IF (PROBLEM.EQ.1) THEN
! 		EX_DISP = X*Y
! 	ELSEIF (PROBLEM.EQ.2) THEN
! 		EX_DISP = (X**2-1.0D0)*(Y**2-1.0D0) + DSIN(X)*DCOS(Y)
! 	ELSEIF (PROBLEM.EQ.3) THEN
! 		EX_DISP = DSQRT(R)*(1.D0-R)*(DSIN(0.5D0*THETA) + DSIN(1.5D0*THETA))
! ! 	ELSEIF (PROBLEM==4) THEN
! ! 		EX_DISP = R**(2.D0/3.D0)*DSIN(2.D0*THETA/3.D0)
! 	ENDIF
! 	
! END FUNCTION EX_DISP
! 
! 
! REAL*8 FUNCTION LDFT2D(PT2D)
! 
! 	TYPE(POINT2D), INTENT(IN) :: PT2D
! 	TYPE(POINT2D) :: PHY_PT
! 	REAL*8 :: X, Y, R, THETA
! 
! 	PHY_PT = GET_PHY_PT(PT2D)
! 	
! 	X = PHY_PT%X; Y = PHY_PT%Y
! 	
! 	CALL GET_RTHETA(R, THETA, PHY_PT)
! 	
! 	IF (PROBLEM.EQ.1) THEN
! 		LDFT2D = 0.D0
! 	ELSEIF (PROBLEM.EQ.2) THEN
! 		LDFT2D = -2.0D0*(-2.0D0 + X**2 + Y**2 - DCOS(Y)*DSIN(X))
! 	ELSEIF (PROBLEM.EQ.3) THEN
! 		LDFT2D = 2.D0*(1.D0 + R + 2.D0*DCOS(THETA))*DSIN(0.5D0*THETA)/R**(1.5D0)
! ! 	ELSEIF (PROBLEM==4) THEN
! ! 		LDFT2D = 0.D0
! 	ENDIF
! 
! 	IF (DABS(R)<=EPS) THEN
! 		LDFT2D = 0.0D0
! 	ENDIF
! 	
! END FUNCTION LDFT2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! NURBS PU W/ FLAT-TOP 1D

REAL*8 FUNCTION LDF1D(REFPT, PATCH)

	REAL*8, INTENT(IN) :: REFPT
	INTEGER, INTENT(IN) :: PATCH
	
	TYPE(FUNCTION_1D) :: PHYPT
	REAL*8 :: X
	
	IF (PATCH==1) THEN 
		PHYPT = MAP_F(REFPT)
	ELSEIF (PATCH==2) THEN 
		PHYPT = MAP_G(REFPT)
	ENDIF
	
	X = PHYPT%VAL(0)
	
	IF (PROBLEM==0) THEN   ! Fourth order DE
		LDF1D = 48.0D0*PI*(-1.0D0 + 2.0D0*X)*DCOS(PI*X) - 8.0D0*PI**3*(-1.0D0 + 2.0D0*X)*(-X + X**2)*DCOS(PI*X) + 24.0D0*DSIN(PI*X) + PI**4*(-X + X**2)**2.0D0*DSIN(PI*X) - 6.0D0*PI**2*(2.0D0*(-1.0D0 + 2.0D0*X)**2 + 4.0D0*(-X + X**2))*DSIN(PI*X)
	ELSEIF (PROBLEM==1) THEN   ! Fourth order DE
		LDF1D = ((-2.0D0 + ALPHA)*ALPHA*X**(-4.0D0 + ALPHA/2.0D0)*((8.0D0 + 2.0D0*ALPHA - ALPHA**2)*X + 8.0D0*(3.0D0 - 4*ALPHA + ALPHA**2)*X**(ALPHA/2.0D0)))/8.0D0
	ELSEIF (PROBLEM==2) THEN  ! Poisson's equation
		LDF1D = -2.0D0
    ELSEIF (PROBLEM==3) THEN  ! Poisson's equation
        LDF1D = PI**2*DSIN(PI*X)
	ENDIF

END FUNCTION LDF1D

REAL*8 FUNCTION EXSOL(REFPT, PATCH)

	REAL*8, INTENT(IN) :: REFPT
	INTEGER, INTENT(IN) :: PATCH
	
	TYPE(FUNCTION_1D) :: PHYPT
	REAL*8 :: X
    
    IF (PATCH==1) THEN 
		PHYPT = MAP_F(REFPT)
	ELSEIF (PATCH==2) THEN 
		PHYPT = MAP_G(REFPT)
	ENDIF
	
	X = PHYPT%VAL(0)
	
	IF (PROBLEM==0) THEN
		EXSOL = DSIN(PI*X)*(X**2 - X)**2
	ELSEIF (PROBLEM==1) THEN
		EXSOL = (X**(0.50D0*ALPHA) - X)**2
	ELSEIF (PROBLEM==2) THEN 
		EXSOL = X**2 - X
    ELSEIF (PROBLEM==3) THEN 
        EXSOL = DSIN(PI*X)
	ENDIF

END FUNCTION EXSOL

END MODULE LOADFUNCTION
