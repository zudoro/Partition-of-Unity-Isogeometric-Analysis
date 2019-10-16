MODULE LOADFUNCTION

	USE GEOMETRY

CONTAINS

REAL*8 FUNCTION EX_DISP(PT2D, PATCH)

	TYPE(POINT2D), INTENT(IN) :: PT2D
	INTEGER, INTENT(IN) :: PATCH
	
	TYPE(POINT2D) :: PHY_PT
	REAL*8 :: R, THETA, X, Y

	PHY_PT = GET_PHY_PT(PT2D, PATCH)
	
	X = PHY_PT%X; Y = PHY_PT%Y
	
	CALL GET_RTHETA(R, THETA, PHY_PT)
	
	IF (PDE=='ELPT') THEN
		IF (PROBLEM==0) THEN
			EX_DISP = X**2*Y**2
		ELSEIF (PROBLEM==1) THEN
			EX_DISP = 1.0D0 - (X**2 + Y**2)
		ELSEIF (PROBLEM==2) THEN
			EX_DISP = (1.0D0 - (X**2 + Y**2))*DEXP(X*Y)
        ELSEIF (PROBLEM==3) THEN 
            EX_DISP = DSQRT(R)*(1.D0-R)*(DSIN(0.5D0*THETA) + DSIN(1.5D0*THETA))
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
    ELSEIF (PDE=='FOUR') THEN 
        IF (PROBLEM==1) THEN
            EX_DISP = (1.0D0 - X**2)**2*(1.0D0 - Y**2)**2
        ELSEIF (PROBLEM==2) THEN 
            EX_DISP = (1.0D0 - X**2)**3*(1.0D0 - Y**2)**3
        ELSEIF (PROBLEM==3) THEN
            EX_DISP = DEXP(X+Y)*(1.0D0 - X**2)**3*(1.0D0 - Y**2)**3
        ELSEIF (PROBLEM==4) THEN 
            EX_DISP = (1.0D0 - R)**(2)*R**(2.0D0)*(DSIN(1.50D0*THETA) - 3.0D0*DSIN(0.50D0*THETA) + DCOS(1.50D0*THETA) - DCOS(0.50D0*THETA))
        ELSEIF (PROBLEM==5) THEN 
            EX_DISP = (1.0D0 - R)**(2)*R**(1.50D0)*(DSIN(1.50D0*THETA) - 3.0D0*DSIN(0.50D0*THETA) + DCOS(1.50D0*THETA) - DCOS(0.50D0*THETA))
        ENDIF
	ENDIF
	
END FUNCTION EX_DISP


REAL*8 FUNCTION LDFT2D(PT2D, PATCH)

	TYPE(POINT2D), INTENT(IN) :: PT2D
	INTEGER, INTENT(IN) :: PATCH
	
	TYPE(POINT2D) :: PHY_PT
	REAL*8 :: X, Y, R, THETA

	PHY_PT = GET_PHY_PT(PT2D, PATCH)
	
	X = PHY_PT%X; Y = PHY_PT%Y
	
	CALL GET_RTHETA(R, THETA, PHY_PT)
	
	IF (PDE=='ELPT') THEN
		IF (PROBLEM==0) THEN
			LDFT2D = -2.0D0*(X**2 + Y**2)
		ELSEIF (PROBLEM==1) THEN
			LDFT2D = 4.0D0
		ELSEIF (PROBLEM==2) THEN
			LDFT2D = DEXP(X*Y)*(4.0D0 + X**4 + 8*X*Y - Y**2 + Y**4 + X**2*(-1.0D0 + 2*Y**2))
        ELSEIF (PROBLEM==3) THEN 
            LDFT2D = 2.D0*(1.D0 + R + 2.D0*DCOS(THETA))*DSIN(0.5D0*THETA)/R**(1.5D0)
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
    ELSEIF (PDE=='FOUR') THEN 
        IF (PROBLEM==1) THEN
            LDFT2D = 8.0D0*(10.0D0 + 3.0D0*X**4 - 18.0D0*Y**2 + 3.0D0*Y**4 + 18.0D0*X**2*(2.0D0*Y**2 - 1.0D0))
        ELSEIF (PROBLEM==2) THEN
            LDFT2D = 72.0D0*(3.0D0 - 14.0D0*Y**2 + 8.0D0*Y**4 - Y**6 + X**6*(-1.0D0 + 5.0D0*Y**2) + X**4*(8.0D0 - 45.0D0*Y**2 + 25.0D0*Y**4) + X**2*(-14.0D0 + 66.0D0*Y**2 - 45.0D0*Y**4 + 5.0D0*Y**6))
        ELSEIF (PROBLEM==3) THEN
            LDFT2D = 4.0D0*DEXP(X + Y)*(31.0D0 + 96.0D0*Y - 147.0D0*Y**2 - 168.0D0*Y**3 + 51.0D0*Y**4 + 24.0D0*Y**5 - 7.0D0*Y**6 + 12.0D0*X**5*(2.0D0 + 6.0D0*Y - 15.0D0*Y**2 - 12.0D0*Y**3 + 12.0D0*Y**4 + 6.0D0*Y**5 + Y**6) + X**6*(-7.0D0 - 60.0D0*Y + 21.0D0*Y**2 + 96.0D0*Y**3 + 57.0D0*Y**4 + 12.0D0*Y**5 + Y**6) + 24.0D0*X**3*(-7.0D0 - 6.0D0*Y + 30.0D0*Y**2 + 12.0D0*Y**3 - 27.0D0*Y**4 - 6.0D0*Y**5 + 4.0D0*Y**6) - 12.0D0*X*(-8.0D0 - 6.0D0*Y + 33.0D0*Y**2 + 12.0D0*Y**3 - 30.0D0*Y**4 - 6.0D0*Y**5 + 5.0D0*Y**6) + 3.0D0*X**2*(-49.0D0 - 132.0D0*Y + 255.0D0*Y**2 + 240.0D0*Y**3 - 141.0D0*Y**4 - 60.0D0*Y**5 + 7.0D0*Y**6) + 3.0D0*X**4*(17.0D0 + 120.0D0*Y - 141.0D0*Y**2 - 216.0D0*Y**3 + 33.0D0*Y**4 + 48.0D0*Y**5 + 19.0D0*Y**6))
        ELSEIF (PROBLEM==4) THEN 
            LDFT2D = (15.0D0*(1.0D0 + 7.0D0*(2.0D0 - 9.0D0*R)*R)*DCOS(THETA*0.50D0) + (-63.0D0 + 5.0D0*R*(54.0D0 + 77.0D0*R))*DCOS((3.0D0*THETA)*0.50D0) + 2.0D0*(-9.0D0 + 25.0D0*(18.0D0 - 49.0D0*R)*R + (-63.0D0 + 5.0D0*R*(54.0D0 + 77.0D0*R))*DCOS(THETA))*DSIN(THETA*0.50D0))/(16.0D0*R**2)
        ELSEIF (PROBLEM==5) THEN 
!             LDFT2D = (8.0D0*(-3.0D0*R*DCOS(0.50D0*THETA) + 2.0D0*DCOS(1.50D0*THETA) + (2.0D0 - 9.0D0*R + 4.0D0*DCOS(THETA))*DSIN(0.50D0*THETA)))/R**1.50D0
            LDFT2D = -R**(-1.50D0)*(24.0D0*R*DCOS(0.50D0*THETA) - 16.0D0*DSQRT(2.0D0)*DSIN(3.0D0*0.50D0*THETA + 0.250D0*PI) + 72.0D0*R*DSIN(0.50D0*THETA))
        ENDIF
	ENDIF
	
END FUNCTION LDFT2D

END MODULE LOADFUNCTION
