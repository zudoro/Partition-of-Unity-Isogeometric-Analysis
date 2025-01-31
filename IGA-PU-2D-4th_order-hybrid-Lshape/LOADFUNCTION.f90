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
	
    IF (PROBLEM==1) THEN ! Fourth order
        EX_DISP = (X**2 - 1.0D0)**2*(Y**2 - 1.0D0)**2*R**(1.5440D0)*(DSIN(2.0D0*THETA/3.0D0) - (1.0D0/3.0D0)*DSIN(2.0D0*THETA))
    elseif (PROBLEM==2) then ! Elliptic
        EX_DISP = R**(2.D0/3.D0)*DSIN(2.D0*THETA/3.D0)*((R*DCOS(THETA))**2 - 1.0D0)*((R*DSIN(THETA))**2 - 1.0D0)
    elseif (PROBLEM==3) then ! Fourth order
        EX_DISP = (X**2 - 1.0D0)**2*(Y**2 - 1.0D0)**2*R**(1.50D0)*(DSIN(2.0D0*THETA/3.0D0) - (1.0D0/3.0D0)*DSIN(2.0D0*THETA))
    elseif (problem==4) then ! Elliptic
        EX_DISP = x*y*(x**2 - 1.0d0)*(y**2 - 1.0d0)
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
	
    IF (PROBLEM==1) THEN ! Fourth order
        LDFT2D = (24.0D0*R**8.9120D0*DCOS(THETA)**12.0D0*DSIN((2.0D0*THETA)/3.0D0) - 16.0D0*R**8.9120D0*DCOS(THETA)**13.0D0*DSIN(THETA) + DCOS(THETA)**11.0D0*DSIN(THETA)*(127.0609920D0*R**6.9120D0 + 64.0D0*R**8.9120D0*DCOS((2.0D0*THETA)/3.0D0) - 64.0D0*R**8.9120D0*DCOS(2.0D0*THETA) - 447.9989760D0*R**8.9120D0*DSIN(THETA)**2) + DCOS(THETA)**10*DSIN((2.0D0*THETA)/3.0D0)*(-275.92482133333333330D0*R**6.9120D0 + 927.9984640D0*R**8.9120D0*DSIN(THETA)**2) + DCOS(THETA)**6*DSIN((2.0D0*THETA)/3.0D0)*(-39.440231571253730D0*R**2.9120D0 - (68.0D0/9.0D0)*R**2.9120D0 + 1807.91057554576040D0*R**4.9120D0*DSIN(THETA)**2 - 6399.959763044430D0*R**6.9120D0*DSIN(THETA)**4) + DCOS(THETA)**8*DSIN((2.0D0*THETA)/3.0D0)*(465.51313147984910D0*R**4.9120D0 - 2593.19462323703240D0*R**6.9120D0*DSIN(THETA)**2+ 3409.07930864607170D0*R**8.9120D0*DSIN(THETA)**4) + DCOS(THETA)**2*DSIN((2.0D0*THETA)/3.0D0)*DSIN(THETA)**4*(-118.32069471376123*R**2.9120D0 - (68.0D0/3.0D0)*R**2.9120D0 + 1807.91057554576080D0*R**4.9120D0*DSIN(THETA)**2 - 2593.19462323703240D0*R**6.9120D0*DSIN(THETA)**4 + 927.9984640D0*R**8.9120D0*DSIN(THETA)**6) + DCOS(THETA)**4*DSIN((2.0D0*THETA)/3.0D0)*(-0.45870613081757970D0*R**0.9120D0 + (-118.320694713761230D0*R**2.9120D0 - (68.0D0/3.0D0)*R**2.9120D0)*DSIN(THETA)**2 - 6399.9597630444280D0*R**6.9120D0*DSIN(THETA)**6 + 3409.079308646070D0*R**8.9120D0*DSIN(THETA)**8) + DCOS(THETA)**9*DSIN(THETA)*(-161.360926040064*R**4.9120D0 + 1505.0172033761282*R**6.9120D0*DSIN(THETA)**2 - 2159.99556266939950D0*R**8.9120D0*DSIN(THETA)**4 + DCOS(2.0D0*THETA)*(254.1219840D0*R**6.9120D0 - 255.99931733333330D0*R**8.9120D0*DSIN(THETA)**2) + DCOS((2.0D0*THETA)/3.0D0)*(-292.0479099259259*R**6.9120D0 + 293.925243259259330D0*R**8.9120D0*DSIN(THETA)**2)) + DCOS(THETA)**7*DSIN(THETA)*(-42.444539215872034*R**2.9120D0 + 24.0D0*R**2.9120D0 - 880.2756362403836*R**4.9120D0*DSIN(THETA)**2 + 3879.7466501283850d0*R**6.9120D0*DSIN(THETA)**4 + DCOS(2.0D0*THETA)*(-17.5779840d0*R**4.9120D0 + 508.2439680D0*R**6.9120D0*DSIN(THETA)**2 - (8999968.0d0/4687.0d0)*R**8.9120D0*DSIN(THETA)**4) + DCOS((2.0D0*THETA)/3.0D0)*(55.503909925925890d0*R**4.9120D0 - 584.09581985185170D0*R**6.9120D0*DSIN(THETA)**2 + 229.925243259259280d0*R**8.9120D0*DSIN(THETA)**4)) + DCOS(THETA)**5*DSIN(THETA)*(-4.0854787440640d0*R**0.9120D0 - 127.333617647616010d0*R**2.9120D0*DSIN(THETA)**2 + 72.0D0*R**2.9120D0*DSIN(THETA)**2 + 3879.746650128385*R**6.9120D0*DSIN(THETA)**6 - 2159.9955626693986*R**8.9120D0*DSIN(THETA)**8 + DCOS((2.0D0*THETA)/3.0D0)*(55.5039099259256260D0*R**4.9120D0*DSIN(THETA)**2 - 229.92524325925953*R**8.9120D0*DSIN(THETA)**6) + DCOS(2.0D0*THETA)*(- 17.5779840d0*R**4.9120D0*DSIN(THETA)**2 + (8999968.0d0/4687.0d0)*R**8.9120D0*DSIN(THETA)**6)) + DCOS(THETA)*DSIN(THETA)**3*(-4.0854787440640d0*R**0.9120D0*DSIN(THETA)**2 + (-42.444539215872034*R**2.9120D0 + 24.0D0*R**2.9120D0)*DSIN(THETA)**4 - 161.360926040064*R**4.9120D0*DSIN(THETA)**6 + 127.0609920D0*R**6.9120D0*DSIN(THETA)**8 - 16.0D0*R**8.9120D0*DSIN(THETA)**10 + DCOS((2.0D0*THETA)/3.0D0)*(- 55.50390992592589*R**4.9120D0*DSIN(THETA)**4 + 292.0479099259259*R**6.9120D0*DSIN(THETA)**6 - 64.0D0*R**8.9120D0*DSIN(THETA)**8)) + DCOS(THETA)**3*DSIN(THETA)*(DSIN(THETA)**4*(-127.33361764761605*R**2.9120D0 + 72.*R**2.9120D0 - 880.2756362403836*R**4.9120D0*DSIN(THETA)**2 + 1505.0172033761278*R**6.9120D0*DSIN(THETA)**4 - 447.9989760D0*R**8.9120D0*DSIN(THETA)**6) + DCOS((2.0D0*THETA)/3.0D0)*(- 55.50390992592580D0*R**4.9120D0*DSIN(THETA)**4 + 584.09581985185180D0*R**6.9120D0*DSIN(THETA)**6 - 293.92524325925930D0*R**8.9120D0*DSIN(THETA)**8) + DCOS(2.0D0*THETA)*(17.5779840d0*R**4.9120D0*DSIN(THETA)**4 - 508.2439680D0*R**6.9120D0*DSIN(THETA)**6 + 255.99931733333330D0*R**8.9120D0*DSIN(THETA)**8)) - 1.0213696860159989*R**0.9120D0*DSIN(2.0D0*THETA)**3 - 44.932169387520d0*R**4.9120D0*DSIN(2.0D0*THETA)**5 - 26.9999466667093450d0*R**8.9120D0*DSIN(2.0D0*THETA)**7 + DSIN((2.0D0*THETA)/3.0D0)*(-0.4587061308175797*R**0.9120D0*DSIN(THETA)**4 + (-39.440231571253730D0*R**2.9120D0 - (68.0D0/9.0D0)*R**2.9120D0)*DSIN(THETA)**6 + 465.51313147984910D0*R**4.9120D0*DSIN(THETA)**8 - 275.92482133333333330D0*R**6.9120D0*DSIN(THETA)**10 + 24.*R**8.9120D0*DSIN(THETA)**12 - 0.22935306540879160d0*R**0.9120D0*DSIN(2.0D0*THETA)**2 + 167.799680508239*R**4.9120D0*DSIN(2.0D0*THETA)**4 + 78.28377639518963*R**8.9120D0*DSIN(2.0D0*THETA)**6) + 4.3944960D0*R**4.9120D0*DSIN(THETA)**6*DSIN(4*THETA) - 63.5304960D0*R**6.9120D0*DSIN(THETA)**8*DSIN(4*THETA) + 16.0D0*R**8.9120D0*DSIN(THETA)**10*DSIN(4*THETA))/R**3.3680D0
    elseif (PROBLEM==2) then ! Elliptic
        ldft2d = -dsin((-8.0d0 + 128.0d0*r**2 - 49.0d0*r**4)*theta + 13.0d0*r**4*theta*dcos(4.0d0*theta) - 18.0d0*r**4*dsin(4.0d0*theta))/(27.0d0*r**(4.0d0/3.0d0))
    elseif (PROBLEM==3) then ! Fourth order
        ldft2d = -((421120.0d0 + 8668672.0d0*r**2 - 59568064.0d0*r**4 + 51593408.0d0*r**6 - 13743090.0d0*r**8 + 6.0d0*(120960.0d0 + 532224.0d0*r**2 + 4465312.0d0*r**4 - 5795872.0d0*r**6 + 1914067.0d0*r**8)*dcos((4.0d0*theta)/3.0d0) - 12.0d0*r**4*(511624.0d0 - 1842184.0d0*r**2 + 665329.0d0*r**4)*dcos((8.0d0*theta)/3.0d0) - 5235520.0d0*r**4*dcos(4.0d0*theta) - 13721792.0d0*r**6*dcos(4.0d0*theta) + 10534552.0d0*r**8*dcos(4.0d0*theta) - 1887840.0d0*r**4*dcos((16.0d0*theta)/3.0d0) - 1609632.0d0*r**6*dcos((16.0d0*theta)/3.0d0) - 1701708.0d0*r**8*dcos((16.0d0*theta)/3.0d0) + 671187.0d0*r**8*dcos((20.0d0*theta)/3.0d0) + 554330.0d0*r**8*dcos(8.0d0*theta) + 184275.0d0*r**8*dcos((28.0d0*theta)/3.0d0))*dSin((2.0d0*theta)/3.0d0))/(165888.0d0*r**2.50d0)
    elseif (problem==4) then !Elliptic
        LDFT2D = -6.0D0*X*Y*(X**2 + Y**2 - 2.0D0)
    ENDIF
	
END FUNCTION LDFT2D

END MODULE LOADFUNCTION
