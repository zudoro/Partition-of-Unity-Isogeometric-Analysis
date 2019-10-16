MODULE NURBS_BASIS

	USE GLBVAR
	USE PATCH_MAPPING

    IMPLICIT NONE

CONTAINS

!!  GET BASIS FUNCTION ON PHYSICAL SPACE


!! ----------------------------- [[[[[ MODIFIED B-SPLINE BASIS FUNCTIONS FOR INTEGRATION]]]]] ----------------------------- 
!-----------------------------------------------------------------------------------------------------------

SUBROUTINE GET_PATCHWISE_BS_HYBRID(SF, INDX, PAR_PT, PATCH, JACOB)
	INTEGER, INTENT(IN) :: PATCH
	TYPE(POINT2D), INTENT(IN) :: PAR_PT
	TYPE(MATRIX_22), INTENT(IN) :: JACOB
	
	TYPE(FVALUE), INTENT(OUT) :: SF((BASIS_KVEC(PATCH, 1)%POLY_ORDER+1)*(BASIS_KVEC(PATCH, 2)%POLY_ORDER+1+UBOUND(EXPO_ARRAY, 1)))
	TYPE(INT2D), INTENT(OUT) :: INDX((BASIS_KVEC(PATCH, 1)%POLY_ORDER+1)*(BASIS_KVEC(PATCH, 2)%POLY_ORDER+1+UBOUND(EXPO_ARRAY, 1)))
	
	TYPE(DIFF_BSPLINES) :: NX, NY, BY
	TYPE(MATRIX_22) :: INV_JACOB
	TYPE(FUNCTION_1D) :: PUFT, ENRICH, RFG, PHYPU
	TYPE(FUNCTION_2D) :: ENRICH2D
	TYPE(POINT2D) :: REF_PT, PHY_PT
	REAL*8 :: RADI, THETA
	
	REAL*8 :: DJ, DJXI, DJETA
	TYPE(SECOND_PARTIAL_DERIVATIVES) :: SECOND_PD
	INTEGER :: I, J, K, II, JJ, KK

	REAL*8 :: JXI, JETA, S(2,3)
    INTEGER :: ETA_NUMBS, NY_INDX
    
! ---------------------------------------------------------------------------------------------------------------------------- !
	INV_JACOB = .INVERSE.JACOB
	
	SECOND_PD = GET_PARTIAL(PAR_PT, PATCH)				!! P^2(X)/P(XI^2), P^2(X)/P(XI,ETA), P^2(X)/P(ETA^2), P^2(Y)/P(XI^2), P^2(Y)/P(XI,ETA), P^2(Y)/P(ETA^2)

 	DJ = JACOB%ENT(1,1)*JACOB%ENT(2,2)-JACOB%ENT(1,2)*JACOB%ENT(2,1)      ! determinent of Jacobian matrix
	JXI = SECOND_PD%X(1)*JACOB%ENT(2,2)+JACOB%ENT(1,1)*SECOND_PD%Y(2)-SECOND_PD%X(2)*JACOB%ENT(2,1)-JACOB%ENT(1,2)*SECOND_PD%Y(1)      ! derivative of the determinent Jacobian w.r.t. xi
	JETA = SECOND_PD%X(2)*JACOB%ENT(2,2)+JACOB%ENT(1,1)*SECOND_PD%Y(3)-SECOND_PD%X(3)*JACOB%ENT(2,1)-JACOB%ENT(1,2)*SECOND_PD%Y(2)   ! derivative of the determinent Jacobian w.r.t. eta

	S(1,1) = -(JXI*JACOB%ENT(2,2)*JACOB%ENT(2,2)-JETA*JACOB%ENT(2,1)*JACOB%ENT(2,2))/DJ**3+(SECOND_PD%Y(2)*JACOB%ENT(2,2)-SECOND_PD%Y(3)*JACOB%ENT(2,1))/DJ**2
	S(1,2) = (JXI*JACOB%ENT(1,2)*JACOB%ENT(2,2)-JETA*JACOB%ENT(1,2)*JACOB%ENT(2,1))/DJ**3-(SECOND_PD%X(2)*JACOB%ENT(2,2)-SECOND_PD%X(3)*JACOB%ENT(2,1))/DJ**2
	S(1,3) = -(JXI*JACOB%ENT(1,2)*JACOB%ENT(1,2)-JETA*JACOB%ENT(1,1)*JACOB%ENT(1,2))/DJ**3+(SECOND_PD%X(2)*JACOB%ENT(1,2)-SECOND_PD%X(3)*JACOB%ENT(1,1))/DJ**2
	S(2,1) = (JXI*JACOB%ENT(2,1)*JACOB%ENT(2,2)-JETA*JACOB%ENT(2,1)*JACOB%ENT(2,1))/DJ**3-(SECOND_PD%Y(1)*JACOB%ENT(2,2)-SECOND_PD%Y(2)*JACOB%ENT(2,1))/DJ**2
	S(2,2) = -(JXI*JACOB%ENT(1,1)*JACOB%ENT(2,2)-JETA*JACOB%ENT(1,1)*JACOB%ENT(2,1))/DJ**3+(SECOND_PD%X(1)*JACOB%ENT(2,2)-SECOND_PD%X(2)*JACOB%ENT(2,1))/DJ**2
	S(2,3) = (JXI*JACOB%ENT(1,1)*JACOB%ENT(1,2)-JETA*JACOB%ENT(1,1)*JACOB%ENT(1,1))/DJ**3-(SECOND_PD%X(1)*JACOB%ENT(1,2)-SECOND_PD%X(2)*JACOB%ENT(1,1))/DJ**2
	
	NX = GET_ALL_DIFF_BSPLINES(PAR_PT%X, BASIS_KVEC(PATCH, 1), 2)
    BY = GET_ALL_DIFF_BSPLINES(PAR_PT%Y, BASIS_KVEC(PATCH, 2), 2)
    ETA_NUMBS = BY%POLY_ORDER + UBOUND(EXPO_ARRAY, 1)
    
!     print*, 'bs-here-1'
    
    DO J = 1, UBOUND(EXPO_ARRAY, 1)
        BY%N(BY%POLY_ORDER + J, 0) = (PAR_PT%Y*DSQRT(0.50D0))**EXPO_ARRAY(J)
        BY%N(BY%POLY_ORDER + J, 1) = DSQRT(0.50D0)**(EXPO_ARRAY(J))*EXPO_ARRAY(J)*PAR_PT%Y**(EXPO_ARRAY(J) - 1.0D0)
        BY%N(BY%POLY_ORDER + J, 2) = DSQRT(0.50D0)**(EXPO_ARRAY(J))*EXPO_ARRAY(J)*(EXPO_ARRAY(J) - 1.0D0)*PAR_PT%Y**(EXPO_ARRAY(J) - 2.0D0)
    ENDDO
    
!     print*, 'bs-here-2'
    
	NY%INIT = BY%INIT
	NY%POLY_ORDER = BY%POLY_ORDER
    
! ---------------------------------------------------------------------------------------------------------------------------- !
    RFG = RADIUS_MAP(PAR_PT, PATCH)
	PUFT = PHYPU1D(RFG%VAL(0), PATCH)
	PHYPU%VAL(0) = PUFT%VAL(0)
	PHYPU%VAL(1) = PUFT%VAL(1)*RFG%VAL(1)
	PHYPU%VAL(2) = PUFT%VAL(2)*RFG%VAL(1)**2 + PUFT%VAL(1)*RFG%VAL(2)
! ---------------------------------------------------------------------------------------------------------------------------- !	
	
! 	if (patch==1) then
!         print*, par_pt%y, puft%val(0)
!     endif
    
	K = 0
	DO J=0, ETA_NUMBS
!         print*, 'j=',j
        NY%N(J,0) = BY%N(J, 0)*PHYPU%VAL(0)
        NY%N(J,1) = BY%N(J, 1)*PHYPU%VAL(0) + BY%N(J, 0)*PHYPU%VAL(1)
        NY%N(J,2) = BY%N(J, 2)*PHYPU%VAL(0) + 2.0D0*BY%N(J, 1)*PHYPU%VAL(1) + BY%N(J, 0)*PHYPU%VAL(2)
        IF (J>NY%POLY_ORDER) THEN
            NY_INDX = BASIS_KVEC(PATCH, 2)%LENGTH - NY%POLY_ORDER + (J - NY%POLY_ORDER) - 1
        ELSE 
            NY_INDX = NY%INIT+J
        ENDIF
!         print*, 'bs-here-3'
		DO I=0, NX%POLY_ORDER
            K = K + 1
            INDX(K) = INT2D(NX%INIT+I, NY_INDX)
			!--------------------- periodic condition ----------------------!
! 			IF (INDX(K)%A==LOC_NUMBS(PATCH, 1)-1) THEN
!                 INDX(K)%A = 0
! 			ENDIF
			!----------------------------------------------------------------------!
			SF(K)%D00 = NX%N(I,0)*NY%N(J,0)
			SF(K)%D10 = NX%N(I,1)*NY%N(J,0)*INV_JACOB%ENT(1,1) + NX%N(I,0)*NY%N(J,1)*INV_JACOB%ENT(2,1)
			SF(K)%D01 = NX%N(I,1)*NY%N(J,0)*INV_JACOB%ENT(1,2) + NX%N(I,0)*NY%N(J,1)*INV_JACOB%ENT(2,2)
			SF(K)%D20 = NX%N(I,2)*NY%N(J,0)*INV_JACOB%ENT(1,1)**2 + 2.0D0*NX%N(I,1)*NY%N(J,1)*INV_JACOB%ENT(1,1)*INV_JACOB%ENT(2,1) + NX%N(I,0)*NY%N(J,2)*INV_JACOB%ENT(2,1)**2 + NX%N(I,1)*NY%N(J,0)*S(1,1) + NX%N(I,0)*NY%N(J,1)*S(2,1)
			SF(K)%D11 = NX%N(I,2)*NY%N(J,0)*INV_JACOB%ENT(1,1)*INV_JACOB%ENT(1,2) + NX%N(I,1)*NY%N(J,1)*(INV_JACOB%ENT(1,1)*INV_JACOB%ENT(2,2)+INV_JACOB%ENT(1,2)*INV_JACOB%ENT(2,1)) + NX%N(I,0)*NY%N(J,2)*INV_JACOB%ENT(2,1)*INV_JACOB%ENT(2,2) + NX%N(I,1)*NY%N(J,0)*S(1,2) + NX%N(I,0)*NY%N(J,1)*S(2,2)
			SF(K)%D02 = NX%N(I,2)*NY%N(J,0)*INV_JACOB%ENT(1,2)**2 + 2.0D0*NX%N(I,1)*NY%N(J,1)*INV_JACOB%ENT(1,2)*INV_JACOB%ENT(2,2) + NX%N(I,0)*NY%N(J,2)*INV_JACOB%ENT(2,2)**2 + NX%N(I,1)*NY%N(J,0)*S(1,3) + NX%N(I,0)*NY%N(J,1)*S(2,3)
		ENDDO
	ENDDO
	
END SUBROUTINE GET_PATCHWISE_BS_HYBRID


SUBROUTINE GET_PATCHWISE_BS(SF, INDX, PAR_PT, PATCH, JACOB)
	INTEGER, INTENT(IN) :: PATCH
	TYPE(POINT2D), INTENT(IN) :: PAR_PT
	TYPE(MATRIX_22), INTENT(IN) :: JACOB
	
	TYPE(FVALUE), INTENT(OUT) :: SF((BASIS_KVEC(PATCH, 1)%POLY_ORDER+1)*(BASIS_KVEC(PATCH, 2)%POLY_ORDER+1))
	TYPE(INT2D), INTENT(OUT) :: INDX((BASIS_KVEC(PATCH, 1)%POLY_ORDER+1)*(BASIS_KVEC(PATCH, 2)%POLY_ORDER+1))
	
	TYPE(DIFF_BSPLINES) :: NX, NY, BY
	TYPE(MATRIX_22) :: INV_JACOB
	TYPE(FUNCTION_1D) :: PUFT, ENRICH, RFG, PHYPU
	TYPE(FUNCTION_2D) :: ENRICH2D
	TYPE(POINT2D) :: REF_PT, PHY_PT
	REAL*8 :: RADI, THETA
	
	REAL*8 :: DJ, DJXI, DJETA
	TYPE(SECOND_PARTIAL_DERIVATIVES) :: SECOND_PD
	INTEGER :: I, J, K, II, JJ, KK

	REAL*8 :: JXI, JETA, S(2,3)

! ---------------------------------------------------------------------------------------------------------------------------- !
	INV_JACOB = .INVERSE.JACOB
	
	SECOND_PD = GET_PARTIAL(PAR_PT, PATCH)				!! P^2(X)/P(XI^2), P^2(X)/P(XI,ETA), P^2(X)/P(ETA^2), P^2(Y)/P(XI^2), P^2(Y)/P(XI,ETA), P^2(Y)/P(ETA^2)

 	DJ = JACOB%ENT(1,1)*JACOB%ENT(2,2)-JACOB%ENT(1,2)*JACOB%ENT(2,1)      ! determinent of Jacobian matrix
	JXI = SECOND_PD%X(1)*JACOB%ENT(2,2)+JACOB%ENT(1,1)*SECOND_PD%Y(2)-SECOND_PD%X(2)*JACOB%ENT(2,1)-JACOB%ENT(1,2)*SECOND_PD%Y(1)      ! derivative of the determinent Jacobian w.r.t. xi
	JETA = SECOND_PD%X(2)*JACOB%ENT(2,2)+JACOB%ENT(1,1)*SECOND_PD%Y(3)-SECOND_PD%X(3)*JACOB%ENT(2,1)-JACOB%ENT(1,2)*SECOND_PD%Y(2)   ! derivative of the determinent Jacobian w.r.t. eta

	S(1,1) = -(JXI*JACOB%ENT(2,2)*JACOB%ENT(2,2)-JETA*JACOB%ENT(2,1)*JACOB%ENT(2,2))/DJ**3+(SECOND_PD%Y(2)*JACOB%ENT(2,2)-SECOND_PD%Y(3)*JACOB%ENT(2,1))/DJ**2
	S(1,2) = (JXI*JACOB%ENT(1,2)*JACOB%ENT(2,2)-JETA*JACOB%ENT(1,2)*JACOB%ENT(2,1))/DJ**3-(SECOND_PD%X(2)*JACOB%ENT(2,2)-SECOND_PD%X(3)*JACOB%ENT(2,1))/DJ**2
	S(1,3) = -(JXI*JACOB%ENT(1,2)*JACOB%ENT(1,2)-JETA*JACOB%ENT(1,1)*JACOB%ENT(1,2))/DJ**3+(SECOND_PD%X(2)*JACOB%ENT(1,2)-SECOND_PD%X(3)*JACOB%ENT(1,1))/DJ**2
	S(2,1) = (JXI*JACOB%ENT(2,1)*JACOB%ENT(2,2)-JETA*JACOB%ENT(2,1)*JACOB%ENT(2,1))/DJ**3-(SECOND_PD%Y(1)*JACOB%ENT(2,2)-SECOND_PD%Y(2)*JACOB%ENT(2,1))/DJ**2
	S(2,2) = -(JXI*JACOB%ENT(1,1)*JACOB%ENT(2,2)-JETA*JACOB%ENT(1,1)*JACOB%ENT(2,1))/DJ**3+(SECOND_PD%X(1)*JACOB%ENT(2,2)-SECOND_PD%X(2)*JACOB%ENT(2,1))/DJ**2
	S(2,3) = (JXI*JACOB%ENT(1,1)*JACOB%ENT(1,2)-JETA*JACOB%ENT(1,1)*JACOB%ENT(1,1))/DJ**3-(SECOND_PD%X(1)*JACOB%ENT(1,2)-SECOND_PD%X(2)*JACOB%ENT(1,1))/DJ**2
	
	NX = GET_ALL_DIFF_BSPLINES(PAR_PT%X, BASIS_KVEC(PATCH, 1), 2)
    BY = GET_ALL_DIFF_BSPLINES(PAR_PT%Y, BASIS_KVEC(PATCH, 2), 2)
    
	NY%INIT = BY%INIT
	NY%POLY_ORDER = BY%POLY_ORDER

! ---------------------------------------------------------------------------------------------------------------------------- !
    RFG = RADIUS_MAP(PAR_PT, PATCH)
	PUFT = PHYPU1D(RFG%VAL(0), PATCH)
	PHYPU%VAL(0) = PUFT%VAL(0)
	PHYPU%VAL(1) = PUFT%VAL(1)*RFG%VAL(1)
	PHYPU%VAL(2) = PUFT%VAL(2)*RFG%VAL(1)**2 + PUFT%VAL(1)*RFG%VAL(2)
! ---------------------------------------------------------------------------------------------------------------------------- !	
	
! 	if (patch==1) then
!         print*, par_pt%y, puft%val(0)
!     endif
    
	K = 0
	DO J=0, NY%POLY_ORDER
        NY%N(J,0) = BY%N(J, 0)*PHYPU%VAL(0)
        NY%N(J,1) = BY%N(J, 1)*PHYPU%VAL(0) + BY%N(J, 0)*PHYPU%VAL(1)
        NY%N(J,2) = BY%N(J, 2)*PHYPU%VAL(0) + 2.0D0*BY%N(J, 1)*PHYPU%VAL(1) + BY%N(J, 0)*PHYPU%VAL(2)
		DO I=0, NX%POLY_ORDER
			K = K + 1
			INDX(K) = INT2D(NX%INIT+I, NY%INIT+J)
			!--------------------- periodic condition ----------------------!
! 			IF (INDX(K)%A==LOC_NUMBS(PATCH, 1)-1) THEN
!                 INDX(K)%A = 0
! 			ENDIF
			!----------------------------------------------------------------------!
			SF(K)%D00 = NX%N(I,0)*NY%N(J,0)
			SF(K)%D10 = NX%N(I,1)*NY%N(J,0)*INV_JACOB%ENT(1,1) + NX%N(I,0)*NY%N(J,1)*INV_JACOB%ENT(2,1)
			SF(K)%D01 = NX%N(I,1)*NY%N(J,0)*INV_JACOB%ENT(1,2) + NX%N(I,0)*NY%N(J,1)*INV_JACOB%ENT(2,2)
			SF(K)%D20 = NX%N(I,2)*NY%N(J,0)*INV_JACOB%ENT(1,1)**2 + 2.0D0*NX%N(I,1)*NY%N(J,1)*INV_JACOB%ENT(1,1)*INV_JACOB%ENT(2,1) + NX%N(I,0)*NY%N(J,2)*INV_JACOB%ENT(2,1)**2 + NX%N(I,1)*NY%N(J,0)*S(1,1) + NX%N(I,0)*NY%N(J,1)*S(2,1)
			SF(K)%D11 = NX%N(I,2)*NY%N(J,0)*INV_JACOB%ENT(1,1)*INV_JACOB%ENT(1,2) + NX%N(I,1)*NY%N(J,1)*(INV_JACOB%ENT(1,1)*INV_JACOB%ENT(2,2)+INV_JACOB%ENT(1,2)*INV_JACOB%ENT(2,1)) + NX%N(I,0)*NY%N(J,2)*INV_JACOB%ENT(2,1)*INV_JACOB%ENT(2,2) + NX%N(I,1)*NY%N(J,0)*S(1,2) + NX%N(I,0)*NY%N(J,1)*S(2,2)
			SF(K)%D02 = NX%N(I,2)*NY%N(J,0)*INV_JACOB%ENT(1,2)**2 + 2.0D0*NX%N(I,1)*NY%N(J,1)*INV_JACOB%ENT(1,2)*INV_JACOB%ENT(2,2) + NX%N(I,0)*NY%N(J,2)*INV_JACOB%ENT(2,2)**2 + NX%N(I,1)*NY%N(J,0)*S(1,3) + NX%N(I,0)*NY%N(J,1)*S(2,3)
		ENDDO
	ENDDO
	
END SUBROUTINE GET_PATCHWISE_BS


! TYPE(FVALUE) FUNCTION GET_BS_HYBRID(PAR_PT, INDX, PATCH)
!     TYPE(POINT2D), INTENT(IN) :: PAR_PT
!     INTEGER, INTENT(IN) :: INDX(2), PATCH
!     
!     TYPE(FUNCTION_1D) :: PUFT, ENRICH, RFG, NX, NY
!     INTEGER :: NY_INDX
!     
!     NX = GET_DIFF_BSPLINE(PAR_PT%X, BASIS_KVEC(PATCH, 1), INDX(1), 2)
!     
!     IF (INDX(2)>BASIS_KVEC(PATCH, 2)%LENGTH - BS_ORDER(PATCH, 2)-1) THEN 
!         NY_INDX = INDX(2) - (BASIS_KVEC(PATCH, 2)%LENGTH - BS_ORDER(PATCH, 2)-1)
!         NY%VAL(0) = (PAR_PT%Y*DSQRT(0.50D0))**EXPO_ARRAY(NY_INDX)
!         NY%VAL(1) = DSQRT(0.50D0)**(EXPO_ARRAY(INDX(2)+1))*EXPO_ARRAY(INDX(2)+1)*PAR_PT%Y**(EXPO_ARRAY(NY_INDX) - 1.0D0)
!         NY%VAL(2) = DSQRT(0.50D0)**(EXPO_ARRAY(INDX(2)+1))*EXPO_ARRAY(INDX(2)+1)*(EXPO_ARRAY(INDX(2)+1) - 1.0D0)*PAR_PT%Y**(EXPO_ARRAY(NY_INDX) - 2.0D0)
!     ELSE
!         NY = GET_DIFF_BSPLINE(PAR_PT%Y, BASIS_KVEC(PATCH, 2), INDX(2), 2)
!     ENDIF
!     
!     ! ---------------------------------------------------------------------------------------------------------------------------- !
!     RFG = RADIUS_MAP(PAR_PT, PATCH)
! 	PUFT = PHYPU1D(RFG%VAL(0), PATCH)
! 	PUFT%VAL(1) = PUFT%VAL(1)*RFG%VAL(1)
! 	PUFT%VAL(2) = PUFT%VAL(2)*RFG%VAL(1)**2 + PUFT%VAL(1)*RFG%VAL(2)
!     ! ---------------------------------------------------------------------------------------------------------------------------- !
!     
!     GET_BS%D00 = NX%VAL(0)*NY%VAL(0)*PUFT%VAL(0)
!     GET_BS%D10 = NX%VAL(1)*NY%VAL(0)*PUFT%VAL(0)
!     GET_BS%D01 = NX%VAL(0)*(NY%VAL(1)*PUFT%VAL(0) + NY%VAL(0)*PUFT%VAL(1))
!     GET_BS%D20 = NX%VAL(2)*NY%VAL(0)*PUFT%VAL(0)
!     GET_BS%D11 = NX%VAL(1)*(NY%VAL(1)*PUFT%VAL(0) + NY%VAL(0)*PUFT%VAL(1))
!     GET_BS%D02 = NX%VAL(0)*(NY%VAL(2)*PUFT%VAL(0) + 2.0D0*NY%VAL(1)*PUFT%VAL(1) + NY%VAL(0)*PUFT%VAL(2))
!     
! END FUNCTION GET_BS_HYBRID


TYPE(FVALUE) FUNCTION GET_BS(PAR_PT, INDX, PATCH)
    TYPE(POINT2D), INTENT(IN) :: PAR_PT
    INTEGER, INTENT(IN) :: INDX(2), PATCH
    
    TYPE(FUNCTION_1D) :: PUFT, ENRICH, RFG, NX, NY, PHYPU
    INTEGER :: NY_INDX
    
    NX = GET_DIFF_BSPLINE(PAR_PT%X, BASIS_KVEC(PATCH, 1), INDX(1), 2)
    
    IF (PATCH==1) THEN 
        IF (INDX(2)>BASIS_KVEC(PATCH, 2)%LENGTH - BS_ORDER(PATCH, 2)-1) THEN 
            NY_INDX = INDX(2) - (BASIS_KVEC(PATCH, 2)%LENGTH - BS_ORDER(PATCH, 2)-1)
            NY%VAL(0) = (PAR_PT%Y*DSQRT(0.50D0))**EXPO_ARRAY(NY_INDX)
            NY%VAL(1) = DSQRT(0.50D0)**(EXPO_ARRAY(INDX(2)+1))*EXPO_ARRAY(INDX(2)+1)*PAR_PT%Y**(EXPO_ARRAY(NY_INDX) - 1.0D0)
            NY%VAL(2) = DSQRT(0.50D0)**(EXPO_ARRAY(INDX(2)+1))*EXPO_ARRAY(INDX(2)+1)*(EXPO_ARRAY(INDX(2)+1) - 1.0D0)*PAR_PT%Y**(EXPO_ARRAY(NY_INDX) - 2.0D0)
        ELSE
            NY = GET_DIFF_BSPLINE(PAR_PT%Y, BASIS_KVEC(PATCH, 2), INDX(2), 2)
        ENDIF
	ELSE 
        NY = GET_DIFF_BSPLINE(PAR_PT%Y, BASIS_KVEC(PATCH, 2), INDX(2), 2)
    ENDIF
    
    ! ---------------------------------------------------------------------------------------------------------------------------- !
    RFG = RADIUS_MAP(PAR_PT, PATCH)
	PUFT = PHYPU1D(RFG%VAL(0), PATCH)
	PHYPU%VAL(0) = PUFT%VAL(0)
	PHYPU%VAL(1) = PUFT%VAL(1)*RFG%VAL(1)
	PHYPU%VAL(2) = PUFT%VAL(2)*RFG%VAL(1)**2 + PUFT%VAL(1)*RFG%VAL(2)
    ! ---------------------------------------------------------------------------------------------------------------------------- !
    
    GET_BS%D00 = NX%VAL(0)*NY%VAL(0)*PHYPU%VAL(0)
    GET_BS%D10 = NX%VAL(1)*NY%VAL(0)*PHYPU%VAL(0)
    GET_BS%D01 = NX%VAL(0)*(NY%VAL(1)*PHYPU%VAL(0) + NY%VAL(0)*PHYPU%VAL(1))
    GET_BS%D20 = NX%VAL(2)*NY%VAL(0)*PHYPU%VAL(0)
    GET_BS%D11 = NX%VAL(1)*(NY%VAL(1)*PHYPU%VAL(0) + NY%VAL(0)*PHYPU%VAL(1))
    GET_BS%D02 = NX%VAL(0)*(NY%VAL(2)*PHYPU%VAL(0) + 2.0D0*NY%VAL(1)*PHYPU%VAL(1) + NY%VAL(0)*PHYPU%VAL(2))
    
END FUNCTION GET_BS

!--------------------------PHYSICAL PU-1D-----------------------------

TYPE(FUNCTION_1D) FUNCTION PHYPU1D(PHYR, PATCH)     ! Phi_g \circ T(xi, eta)

    REAL*8, INTENT(IN) :: PHYR
    INTEGER, INTENT(IN) :: PATCH
    
    TYPE(FUNCTION_1D) :: T_MAP

    T_MAP = PATCH_TO_REFPU(PHYR, PATCH)
		
		PHYPU1D%VAL(0) = PU1D(T_MAP%VAL(0), 0)
		PHYPU1D%VAL(1) = PU1D(T_MAP%VAL(0), 1)*T_MAP%VAL(1)
		PHYPU1D%VAL(2) = PU1D(T_MAP%VAL(0), 2)*T_MAP%VAL(1)**2

END FUNCTION PHYPU1D

REAL*8 FUNCTION PU1D(X, FLAGX)

	REAL*8, INTENT(IN) :: X
	INTEGER, INTENT(IN) :: FLAGX

	IF (FLAGX.EQ. 0) THEN													! PU-FUNCTION
		IF ((X + 1.0D0 .GE. EPS) .AND. (X .LT. EPS)) THEN
			IF (PUORDER.EQ.1) THEN
				PU1D=(1.0D0+X)
			ELSEIF (PUORDER.EQ.2) THEN
				PU1D=((1.0D0+X)**2)*(1.0D0-2.0D0*X)
			ELSEIF (PUORDER.EQ.3) THEN
				PU1D=((1.0D0+X)**3)*(1.0D0-3.0D0*X+6.0D0*X**2)
			ELSEIF (PUORDER.EQ.4) THEN
				PU1D=((1.0D0+X)**4)*(1.0D0-4.0D0*X+10.0D0*X**2-20.0D0*X**3)
			ELSEIF (PUORDER.EQ.5) THEN
				PU1D=((1.0D0+X)**5)*(1.0D0-5.0D0*X+15.0D0*X**2-35.0D0*X**3+70.0D0*X**4)
			ELSEIF (PUORDER.EQ.6) THEN
				PU1D=((1.0D0+X)**6.0D0)*(1.0D0-6.0D0*X+21.0D0*X**2.0D0-56.0D0*X**3.0D0+126.0D0*X**4.0D0 &
				-252.0D0*X**5.0D0)
			ELSEIF (PUORDER.EQ.7) THEN
				PU1D=((1.0D0+X)**7)*(1.0D0-7.0D0*X+28.0D0*X**2-84.0D0*X**3+210.0D0*X**4-462.0D0*X**5 &
				+924.0D0*X**6)
			ENDIF
		ELSEIF ((X .GE. EPS) .AND. (X - 1.0D0 .LT. EPS)) THEN
			IF (PUORDER.EQ.1) THEN
				PU1D=(1.0D0-X)
			ELSEIF (PUORDER.EQ.2) THEN
				PU1D=((1.0D0-X)**2)*(1.0D0+2.0D0*X)
			ELSEIF (PUORDER.EQ.3) THEN
				PU1D=((1.0D0-X)**3)*(1.0D0+3.0D0*X+6.0D0*X**2)
			ELSEIF (PUORDER.EQ.4) THEN
				PU1D=((1.0D0-X)**4)*(1.0D0+4.0D0*X+10.0D0*X**2+20.0D0*X**3)
			ELSEIF (PUORDER.EQ.5) THEN
				PU1D=((1.0D0-X)**5)*(1.0D0+5.0D0*X+15.0D0*X**2.0D0+35.0D0*X**3.0D0+70.0D0*X**4)
			ELSEIF (PUORDER.EQ.6) THEN
				PU1D=((1.0D0-X)**6)*(1.0D0+6.0D0*X+21.0D0*X**2+56.0D0*X**3+126.0D0*X**4+252.0D0*X**5)
			ELSEIF (PUORDER.EQ.7) THEN
				PU1D=((1.0D0-X)**7)*(1.0D0+7.0D0*X+28.0D0*X**2+84.0D0*X**3+210.0D0*X**4+462.0D0*X**5 &
				+924.0D0*X**6)
			ENDIF
		ELSE
			PU1D=0.0D0
		ENDIF
	ELSEIF (FLAGX.EQ.1) THEN												! DIFFERENTIATION OF PU-FUNCTION
		IF ((X + 1.0D0 .GE. EPS) .AND. (X .LT. EPS)) THEN
			IF (PUORDER.EQ.1) THEN
				PU1D=1.0D0
			ELSEIF (PUORDER.EQ.2) THEN
				PU1D=-6.0D0*X*(1.0D0 + X)
			ELSEIF (PUORDER.EQ.3) THEN
				PU1D=30.0D0*X**2*(1.0D0 + X)**2
			ELSEIF (PUORDER.EQ.4) THEN
				PU1D=-140.0D0*X**3*(1.0D0 + X)**3
			ELSEIF (PUORDER.EQ.5) THEN
				PU1D=630.0D0*X**4*(1.0D0 + X)**4
			ELSEIF (PUORDER.EQ.6) THEN
				PU1D=-2772.0D0*X**5*(1.0D0 + X)**5
			ELSEIF (PUORDER.EQ.7) THEN
				PU1D=12012.0D0*X**6*(1.0D0 + X)**6
			ENDIF
		ELSEIF ((X .GT. EPS) .AND. (X - 1.0D0 .LE. EPS)) THEN
			IF (PUORDER.EQ.1) THEN
				PU1D=-1.0D0
			ELSEIF (PUORDER.EQ.2) THEN
				PU1D=6.0D0*(-1.0D0 + X)*X
			ELSEIF (PUORDER.EQ.3) THEN
				PU1D=-30.0D0*(-1.0D0 + X)**2*X**2
			ELSEIF (PUORDER.EQ.4) THEN
				PU1D=140.0D0*(-1.0D0 + X)**3*X**3
			ELSEIF (PUORDER.EQ.5) THEN
				PU1D=-630.0D0*(-1.0D0 + X)**4*X**4
			ELSEIF (PUORDER.EQ.6) THEN
				PU1D=2772.0D0*(-1.0D0 + X)**5*X**5
			ELSEIF (PUORDER.EQ.7) THEN
				PU1D=-12012.0D0*(-1.0D0 + X)**6*X**6
			ENDIF
        ELSE
            PU1D=0.0D0
		ENDIF
	ELSEIF (FLAGX.EQ.2) THEN												! DIFFERENTIATION OF PU-FUNCTION
		IF ((X + 1.0D0 .GE. EPS) .AND. (X .LT. EPS)) THEN
			IF (PUORDER.EQ.1) THEN
				PU1D=0.0D0
			ELSEIF (PUORDER.EQ.2) THEN
				PU1D=-6*x - 6*(1.0D0 + x)
			ELSEIF (PUORDER.EQ.3) THEN
				PU1D=60*x**2*(1.0D0 + X) + 60*x*(1.0D0 + X)**2
			ELSEIF (PUORDER.EQ.4) THEN
				PU1D=-420*x**3*(1.0D0 + X)**2 - 420*x**2*(1.0D0 + X)**3
			ELSEIF (PUORDER.EQ.5) THEN
				PU1D=2520*x**4*(1.0D0 + X)**3 + 2520*x**3*(1.0D0 + X)**4
			ELSEIF (PUORDER.EQ.6) THEN
				PU1D=-13860*x**5*(1.0D0 + X)**4 - 13860*x**4*(1.0D0 + X)**5
			ELSEIF (PUORDER.EQ.7) THEN
				PU1D=72072*x**6*(1.0D0 + X)**5 + 72072*x**5*(1.0D0 + X)**6
			ENDIF
		ELSEIF ((X .GT. EPS) .AND. (X - 1.0D0 .LE. EPS)) THEN
			IF (PUORDER.EQ.1) THEN
				PU1D=0.0D0
			ELSEIF (PUORDER.EQ.2) THEN
				PU1D=6*x + 6*(-1.0D0 + X)
			ELSEIF (PUORDER.EQ.3) THEN
				PU1D=-60*x**2*(-1.0D0 + X) - 60*x*(-1.0D0 + X)**2
			ELSEIF (PUORDER.EQ.4) THEN
				PU1D=420*x**3*(-1.0D0 + X)**2 + 420*x**2*(-1.0D0 + X)**3
			ELSEIF (PUORDER.EQ.5) THEN
				PU1D=-2520*x**4*(-1.0D0 + X)**3 - 2520*x**3*(-1.0D0 + X)**4
			ELSEIF (PUORDER.EQ.6) THEN
				PU1D=13860*x**5*(-1.0D0 + X)**4 + 13860*x**4*(-1.0D0 + X)**5
			ELSEIF (PUORDER.EQ.7) THEN
				PU1D=-72072*x**6*(-1.0D0 + X)**5 - 72072*x**5*(-1.0D0 + X)**6
			ENDIF
        ELSE
            PU1D=0.0D0
		ENDIF
	ENDIF
END FUNCTION PU1D

END MODULE NURBS_BASIS
