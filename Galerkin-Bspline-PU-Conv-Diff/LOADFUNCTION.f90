MODULE LOADFUNCTION

	USE GEOMETRY

CONTAINS

REAL*8 FUNCTION EX_DISP(PT2D)

	TYPE(POINT2D), INTENT(IN) :: PT2D
	TYPE(POINT2D) :: PHY_PT
	REAL*8 :: R, THETA, X, Y

	PHY_PT = GET_PHY_PT(PT2D)
	
	X = PHY_PT%X; Y = PHY_PT%Y
	
	CALL GET_RTHETA(R, THETA, PHY_PT)
	
	IF (PDE=='ELPT') THEN
		IF (PROBLEM==0) THEN
			EX_DISP = X**2*Y**2
		ELSEIF (PROBLEM==1) THEN
			EX_DISP = 1.0D0 - (X**2 + Y**2)
		ELSEIF (PROBLEM==2) THEN
			EX_DISP = (1.0D0 - (X**2 + Y**2))*DEXP(X*Y)
		ENDIF
	ELSEIF (PDE=='CNVD') THEN
		IF (PROBLEM==1) THEN
			EX_DISP = X**2*Y**2
		ELSEIF (PROBLEM==2) THEN
			EX_DISP = (1.D0 - DSQRT(X**2 + Y**2))*X**2*Y**2
		ELSEIF (PROBLEM==3) THEN
			EX_DISP = (1.D0 - R)*DCOS(Y)*DSIN(X)
		ELSEIF (PROBLEM.EQ.4) THEN
			IF (DABS(R - 1.0D0)<=EPS) THEN
				EX_DISP = 0.0D0
			ELSE
				EX_DISP = (1.0D0 - X**2)**2*(DSQRT(1.0D0 - X**2) - Y + ((DSQRT(1.0D0 - X**2) + Y)*EPSLN)/(1.0D0 - X**2)**1.50D0)
			ENDIF
		ELSEIF (PROBLEM==5) THEN
			EX_DISP = 1.0D0 - (X**2 + Y**2)
		ENDIF
	ENDIF
	
END FUNCTION EX_DISP


REAL*8 FUNCTION LDFT2D(PT2D)

	TYPE(POINT2D), INTENT(IN) :: PT2D
	TYPE(POINT2D) :: PHY_PT
	REAL*8 :: X, Y, R, THETA

	PHY_PT = GET_PHY_PT(PT2D)
	
	X = PHY_PT%X; Y = PHY_PT%Y
	
	CALL GET_RTHETA(R, THETA, PHY_PT)
	
	IF (PDE=='ELPT') THEN
		IF (PROBLEM==0) THEN
			LDFT2D = -2.0D0*(X**2 + Y**2)
		ELSEIF (PROBLEM==1) THEN
			LDFT2D = 4.0D0
		ELSEIF (PROBLEM==2) THEN
			LDFT2D = DEXP(X*Y)*(4.0D0 + X**4 + 8*X*Y - Y**2 + Y**4 + X**2*(-1.0D0 + 2*Y**2))
		ENDIF
	ELSEIF (PDE=='CNVD') THEN 
		IF (PROBLEM==1) THEN
			LDFT2D = -2.0D0*(Y**2*EPSLN + X**2*(Y + EPSLN))
		ELSEIF (PROBLEM==2) THEN
			LDFT2D = (X**2*Y*(2.0D0*X**2 + 3.0D0*Y**2 - 2.0D0*DSQRT(X**2 + Y**2)) + EPSLN*(2.0D0*X**4 + 2.0D0*Y**4 - 2.0D0*Y**2*DSQRT(X**2 + Y**2) + X**2*(13.0D0*Y**2 - 2.0D0*DSQRT(X**2 + Y**2))))/DSQRT(X**2 + Y**2)
		ELSEIF (PROBLEM.EQ.3) THEN
			LDFT2D = (2.0D0*X*EPSLN*DCOS(X)*DCOS(Y) + DSIN(X)*((Y + EPSLN - 2.0D0*X**2.0D0*EPSLN - 2.0D0*Y**2.0D0*EPSLN + 2.0D0*DSQRT(X**2 + Y**2)*EPSLN)*DCOS(Y) + (-X**2 - Y**2 + DSQRT(X**2 + Y**2) - 2.0D0*Y*EPSLN)*DSIN(Y)))/DSQRT(X**2 + Y**2)
		ELSEIF (PROBLEM.EQ.4) THEN
! 			LDFT2D = ((1.0D0 - X**2)**3.50D0 + EPSLN**2*(2.0D0*DSQRT(1.0D0 - X**2) - 2.0D0*X**2*DSQRT(1.0D0 - X**2) + Y) - 4.0D0*EPSLN*(-1.0D0 + X**2)*(1.0D0 + 5.0D0*X**4 - DSQRT(1.0D0 - X**2)*Y + 3.0D0*X**2*(-2.0D0 + DSQRT(1.0D0 - X**2)*Y)))/(1.0D0 - X**2)**1.50D0
! 			LDFT2D = (1.0D0 + (4.0D0/DSQRT(1.0D0 - X**2) - 4.0D0*Y)*EPSLN + (2.0D0 + Y/DSQRT(1.0D0 - X**2))*EPSLN**2 + X**6*(-1.0D0 - (20.0D0*EPSLN)/DSQRT(1.0D0 - X**2)) + X**4*(3.0D0 + (44.0D0*EPSLN)/DSQRT(1.0D0 - X**2) - 12.0D0*Y*EPSLN) + X**2*(-3.0D0 - (28.0D0*EPSLN)/DSQRT(1.0D0 - X**2) + 16.0D0*Y*EPSLN - 2.0D0*EPSLN**2))/(1.0D0 - X**2)
			LDFT2D = (1.0D0 - X**2)**2
		ELSEIF (PROBLEM==5) THEN
			LDFT2D = 4.0D0*EPSLN + 2.0D0*Y
		ENDIF
	ENDIF
	
END FUNCTION LDFT2D

END MODULE LOADFUNCTION
