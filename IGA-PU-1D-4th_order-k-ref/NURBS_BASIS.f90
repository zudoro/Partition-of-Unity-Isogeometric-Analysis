MODULE NURBS_BASIS

	USE NURBS_PU

	IMPLICIT INTEGER (I-N)
	IMPLICIT REAL(8) (A-H,O-Z)

CONTAINS

!!  GET BASIS FUNCTION ON PHYSICAL SPACE

TYPE(FUNCTION_1D) FUNCTION GET_PHY_BSPLINE1D(INDX, PARPT, PATCH)
	
	INTEGER, INTENT(IN) :: INDX, PATCH
	REAL*8, INTENT(IN) :: PARPT
	
	TYPE(FUNCTION_1D) :: UNITPT, PHYPT, PHYPU, PARBS, JACOBI, NX
	INTEGER :: I, J, K
	
	IF (PATCH==1 .AND. (PARPT<PATCHBDPT(PATCH,1) .OR. PARPT>(PATCHBDPT(PATCH,2) + DELTA))) THEN
		GOTO 11
	ELSEIF (PATCH==2 .AND. (PARPT<(PATCHBDPT(PATCH,1) - DELTA) .OR. PARPT>PATCHBDPT(PATCH,2))) THEN
		GOTO 11
	ENDIF
	
	IF (PATCH==1) THEN 
		JACOBI = MAP_F(PARPT)
	ELSEIF (PATCH==2) THEN 
		JACOBI = MAP_G(PARPT)
	ENDIF
	
	UNITPT = PATCH_TO_UNITLINE(PARPT, PATCH)
	
! 	PHYPU = GET_PHYPU1D(PARPT, PATCH)
	
	NX = GET_DIFF_BSPLINE(UNITPT%VAL(0), BS_KVEC(PATCH), INDX, 2)
	
	GET_PHY_BSPLINE1D%VAL(0) = NX%VAL(0)
	GET_PHY_BSPLINE1D%VAL(1) = NX%VAL(1)*UNITPT%VAL(1)
	GET_PHY_BSPLINE1D%VAL(2) = NX%VAL(2)*(UNITPT%VAL(1))**2 + NX%VAL(1)*UNITPT%VAL(2)
	
! 	GET_PHY_BASIS1D%VAL(0) = PARBS%VAL(0)*PHYPU%VAL(0)
! 	GET_PHY_BASIS1D%VAL(1) = (PARBS%VAL(1)*PHYPU%VAL(0) + PARBS%VAL(0)*PHYPU%VAL(1))*(1.0D0/JACOBI%VAL(1))
! 	GET_PHY_BASIS1D%VAL(2) = (PARBS%VAL(2)*PHYPU%VAL(0) + 2.0D0*PARBS%VAL(1)*PHYPU%VAL(1) + PARBS%VAL(0)*PHYPU%VAL(2))*(1.0D0/JACOBI%VAL(1)**2) - (PARBS%VAL(1)*PHYPU%VAL(0) + PARBS%VAL(0)*PHYPU%VAL(1))*JACOBI%VAL(2)*(1.0D0/JACOBI%VAL(1))**3
	
	GOTO 12
	
	11 CONTINUE
	
	GET_PHY_BSPLINE1D%VAL(0) = 0.0D0
	GET_PHY_BSPLINE1D%VAL(1) = 0.0D0
	GET_PHY_BSPLINE1D%VAL(2) = 0.0D0
	
	12 CONTINUE
	
END FUNCTION GET_PHY_BSPLINE1D


TYPE(FUNCTION_1D) FUNCTION GET_PHY_BASIS1D(INDX, PARPT, PATCH)
	
	INTEGER, INTENT(IN) :: INDX, PATCH
	REAL*8, INTENT(IN) :: PARPT
	
	TYPE(DIFF_BSPLINES) :: ENRICH
	TYPE(FUNCTION_1D) :: UNITPT, PHYPT, PHYPU, PARBS, JACOBI, NX
	INTEGER :: I, J, K
	
	PHYPU = GET_PHYPU1D(PARPT, PATCH)
	
	
	IF (PATCH==1) THEN 
		JACOBI = MAP_F(PARPT)
! 		ENRICH = GET_ALL_ENRICH1D(PARPT, BS_KVEC(PATCH))
! 		NX%VAL(0) = ENRICH%N(INDX,0)
! 		NX%VAL(1) = ENRICH%N(INDX,1)
! 		NX%VAL(2) = ENRICH%N(INDX,2)
		
		NX = GET_DIFF_BSPLINE(PARPT, BS_KVEC(PATCH), INDX, 2)
		
		GET_PHY_BASIS1D%VAL(0) = NX%VAL(0)*PHYPU%VAL(0)
		GET_PHY_BASIS1D%VAL(1) = (NX%VAL(1)*PHYPU%VAL(0) + NX%VAL(0)*PHYPU%VAL(1))*(1.0D0/JACOBI%VAL(1))
		GET_PHY_BASIS1D%VAL(2) = (NX%VAL(2)*PHYPU%VAL(0) + 2.0D0*NX%VAL(1)*PHYPU%VAL(1) + NX%VAL(0)*PHYPU%VAL(2))*(1.0D0/JACOBI%VAL(1)**2) - (NX%VAL(1)*PHYPU%VAL(0) + NX%VAL(0)*PHYPU%VAL(1))*JACOBI%VAL(2)*(1.0D0/JACOBI%VAL(1))**3
		
	ELSEIF (PATCH==2) THEN
		JACOBI = MAP_G(PARPT)
        NX = GET_DIFF_BSPLINE(PARPT, BS_KVEC(PATCH), INDX, 2)
        
		GET_PHY_BASIS1D%VAL(0) = NX%VAL(0)*PHYPU%VAL(0)
		GET_PHY_BASIS1D%VAL(1) = (NX%VAL(1)*PHYPU%VAL(0) + NX%VAL(0)*PHYPU%VAL(1))*(1.0D0/JACOBI%VAL(1))
		GET_PHY_BASIS1D%VAL(2) = (NX%VAL(2)*PHYPU%VAL(0) + 2.0D0*NX%VAL(1)*PHYPU%VAL(1) + NX%VAL(0)*PHYPU%VAL(2))*(1.0D0/JACOBI%VAL(1)**2) - (NX%VAL(1)*PHYPU%VAL(0) + NX%VAL(0)*PHYPU%VAL(1))*JACOBI%VAL(2)*(1.0D0/JACOBI%VAL(1))**3
	ENDIF
	
END FUNCTION GET_PHY_BASIS1D


TYPE(DIFF_BSPLINES) FUNCTION GET_ALL_ENRICH1D(PARPT, BS_KVEC)
  
  REAL*8, INTENT(IN) :: PARPT
  TYPE(KNOT_VECTOR), INTENT(IN) :: BS_KVEC
  
  INTEGER :: I, J, K
  
  GET_ALL_ENRICH1D%INIT = 0
  GET_ALL_ENRICH1D%POLY_ORDER = BS_KVEC%POLY_ORDER
  
  GET_ALL_ENRICH1D%N(0,0) = PARPT**8
  GET_ALL_ENRICH1D%N(1,0) = PARPT**9
  GET_ALL_ENRICH1D%N(2,0) = PARPT**10
  GET_ALL_ENRICH1D%N(3,0) = PARPT**15
  GET_ALL_ENRICH1D%N(4,0) = PARPT**20
!   GET_ALL_ENRICH1D%N(5,0) = PARPT**25
!   GET_ALL_ENRICH1D%N(6,0) = PARPT**30
  
  GET_ALL_ENRICH1D%N(0,1) = 8.0D0*PARPT**7
  GET_ALL_ENRICH1D%N(1,1) = 9.0D0*PARPT**8
  GET_ALL_ENRICH1D%N(2,1) = 10.0D0*PARPT**9
  GET_ALL_ENRICH1D%N(3,1) = 15.0D0*PARPT**14
  GET_ALL_ENRICH1D%N(4,1) = 20.0D0*PARPT**19
!   GET_ALL_ENRICH1D%N(5,1) = 25.0D0*PARPT**24
!   GET_ALL_ENRICH1D%N(6,1) = 30.0D0*PARPT**29
  
  GET_ALL_ENRICH1D%N(0,2) = 7.0D0*8.0D0*PARPT**6
  GET_ALL_ENRICH1D%N(1,2) = 8.0D0*9.0D0*PARPT**7
  GET_ALL_ENRICH1D%N(2,2) = 9.0D0*10.0D0*PARPT**8
  GET_ALL_ENRICH1D%N(3,2) = 14.0D0*15.0D0*PARPT**13
  GET_ALL_ENRICH1D%N(4,2) = 19.0D0*20.0D0*PARPT**18
!   GET_ALL_ENRICH1D%N(5,2) = 24.0D0*25.0D0*PARPT**23
!   GET_ALL_ENRICH1D%N(6,2) = 29.0D0*30.0D0*PARPT**28
  
END FUNCTION GET_ALL_ENRICH1D


SUBROUTINE GET_ALL_PHY_BASIS1D(BS, INDX, PARPT, PATCH) ! PARPT is a point in the parameter space
	
	REAL*8, INTENT(IN) :: PARPT
	INTEGER, INTENT(IN) :: PATCH
	
	TYPE(FUNCTION_1D), INTENT(OUT) :: BS(BS_KVEC(PATCH)%POLY_ORDER+1)
	INTEGER, INTENT(OUT) :: INDX(BS_KVEC(PATCH)%POLY_ORDER+1)
	
	TYPE(DIFF_BSPLINES) :: NX
	TYPE(FUNCTION_1D) :: UNITPT, PHYPT, PHYPU, PARBS, JACOBI
	INTEGER :: I, J, K
	
	PHYPU = GET_PHYPU1D(PARPT, PATCH)

! 	print*, patch, parpt, phypu%val(0)
	
	IF (PATCH==1) THEN 
		JACOBI = MAP_F(PARPT)
!         NX = GET_ALL_ENRICH1D(PARPT, BS_KVEC(PATCH))
        NX = GET_ALL_DIFF_BSPLINES(PARPT, BS_KVEC(PATCH), 2)
        
		K = 0
		DO I = 0, NX%POLY_ORDER
			K = K + 1
			INDX(K) = NX%INIT+I
			
			BS(K)%VAL(0) = NX%N(I, 0)*PHYPU%VAL(0)
			BS(K)%VAL(1) = (NX%N(I,1)*PHYPU%VAL(0) + NX%N(I,0)*PHYPU%VAL(1))*(1.0D0/JACOBI%VAL(1))
			BS(K)%VAL(2) = (NX%N(I,2)*PHYPU%VAL(0) + 2.0D0*NX%N(I,1)*PHYPU%VAL(1) + NX%N(I,0)*PHYPU%VAL(2))*(1.0D0/JACOBI%VAL(1)**2) - (NX%N(I,1)*PHYPU%VAL(0) + NX%N(I,0)*PHYPU%VAL(1))*JACOBI%VAL(2)*(1.0D0/JACOBI%VAL(1))**3
		ENDDO
	ELSEIF (PATCH==2) THEN 
		JACOBI = MAP_G(PARPT)
		NX = GET_ALL_DIFF_BSPLINES(PARPT, BS_KVEC(PATCH), 2)
! 		UNITPT = PATCH_TO_UNITLINE(PARPT, PATCH)

		K = 0
		DO I = 0, NX%POLY_ORDER
			K = K + 1
			INDX(K) = NX%INIT+I
			
			BS(K)%VAL(0) = NX%N(I, 0)*PHYPU%VAL(0)
			BS(K)%VAL(1) = (NX%N(I,1)*PHYPU%VAL(0) + NX%N(I,0)*PHYPU%VAL(1))*(1.0D0/JACOBI%VAL(1))
			BS(K)%VAL(2) = (NX%N(I,2)*PHYPU%VAL(0) + 2.0D0*NX%N(I,1)*PHYPU%VAL(1) + NX%N(I,0)*PHYPU%VAL(2))*(1.0D0/JACOBI%VAL(1)**2) - (NX%N(I,1)*PHYPU%VAL(0) + NX%N(I,0)*PHYPU%VAL(1))*JACOBI%VAL(2)*(1.0D0/JACOBI%VAL(1))**3
		ENDDO
	ENDIF
	
! 	GOTO 12
! 	
! 	11 CONTINUE
! 	
! 	DO I = 0, NX%POLY_ORDER
! 		BS(I+1)%VAL(0) = 0.0D0
! 		BS(I+1)%VAL(1) = 0.0D0
! 		BS(I+1)%VAL(2) = 0.0D0
! 	ENDDO
! 	
! 	12 CONTINUE
		
END SUBROUTINE GET_ALL_PHY_BASIS1D


! !!  FIND DERIVATIVES OF NURBS SURFACE ON TWO-DIMENSIONAL SPACE
! !!  INPUT  : PT - POINT, KVEC - A VALID KNOT VECTOR, DIFF_ORDER - ORDER OF DIFFERENTIATION
! !!  OUTPUT : DERIVATIVES OF NURBS SURFACE
! SUBROUTINE GET_ALL_PHY_NURBS_SURFACE_2D(SF, INDX, REF_PT, JACOB)
! 
! 	TYPE(FVALUE), INTENT(OUT) :: SF((BASIS_KVEC(1)%POLY_ORDER+1)*(BASIS_KVEC(2)%POLY_ORDER+1))
! 	TYPE(INT2D), INTENT(OUT) :: INDX((BASIS_KVEC(1)%POLY_ORDER+1)*(BASIS_KVEC(2)%POLY_ORDER+1))
! 	TYPE(POINT2D), INTENT(IN) :: REF_PT
! 	TYPE(MATRIX_22), INTENT(IN) :: JACOB
! 
! 	TYPE(SECOND_PARTIAL_DERIVATIVES) :: SECOND_PD
! 	TYPE(MATRIX_22) :: INV_JACOB, JACOBIAN
! 	TYPE(FUNCTION_2D) :: DIFF_NURBS
! 	TYPE(DIFF_BSPLINES) :: DBSFUN(2)
! 	REAL(8) :: DWSUM(0:MAX_DIFF_ORDER,0:MAX_DIFF_ORDER), J, JXI, JETA, S(2,3)
! 	INTEGER :: INDX_X, INDX_Y, DIFF_ORDER
! 	INTEGER :: I, K, L, J1, J2, I1, I2, LC_INDX
! 
! 	DIFF_ORDER = 2
! 
! 	INV_JACOB = .INVERSE.JACOB					!! P(XI)/P(X), P(XI)/P(Y) // P(ETA)/P(X), P(ETA)/P(Y)
! 
! ! 	ELSE
! ! 		JACOBIAN = GET_JACOBIAN_MATRIX(REF_PT,PATCH)			!! P(X)/P(XI), P(X)/P(ETA) // P(Y)/P(XI), P(Y)/P(ETA)
! ! 		INV_JACOB = .INVERSE.JACOBIAN					!! P(XI)/P(X), P(XI)/P(Y) // P(ETA)/P(X), P(ETA)/P(Y)
! ! 	ENDIF
! 	
! 	
! 	SECOND_PD = GET_PARTIAL(REF_PT)				!! P^2(X)/P(XI^2), P^2(X)/P(XI,ETA), P^2(X)/P(ETA^2), P^2(Y)/P(XI^2), P^2(Y)/P(XI,ETA), P^2(Y)/P(ETA^2)
! 	
!  	J = JACOB%ENT(1,1)*JACOB%ENT(2,2)-JACOB%ENT(1,2)*JACOB%ENT(2,1)
! 	JXI = SECOND_PD%X(1)*JACOB%ENT(2,2)+JACOB%ENT(1,1)*SECOND_PD%Y(2)-SECOND_PD%X(2)*JACOB%ENT(2,1)-JACOB%ENT(1,2)*SECOND_PD%Y(1)
! 	JETA = SECOND_PD%X(2)*JACOB%ENT(2,2)+JACOB%ENT(1,1)*SECOND_PD%Y(3)-SECOND_PD%X(3)*JACOB%ENT(2,1)-JACOB%ENT(1,2)*SECOND_PD%Y(2)
! 	
! 	S(1,1) = -(JXI*JACOB%ENT(2,2)*JACOB%ENT(2,2)-JETA*JACOB%ENT(2,1)*JACOB%ENT(2,2))/J**3+(SECOND_PD%Y(2)*JACOB%ENT(2,2)-SECOND_PD%Y(3)*JACOB%ENT(2,1))/J**2
! 	S(1,2) = (JXI*JACOB%ENT(1,2)*JACOB%ENT(2,2)-JETA*JACOB%ENT(1,2)*JACOB%ENT(2,1))/J**3-(SECOND_PD%X(2)*JACOB%ENT(2,2)-SECOND_PD%X(3)*JACOB%ENT(2,1))/J**2
! 	S(1,3) = -(JXI*JACOB%ENT(1,2)*JACOB%ENT(1,2)-JETA*JACOB%ENT(1,1)*JACOB%ENT(1,2))/J**3+(SECOND_PD%X(2)*JACOB%ENT(1,2)-SECOND_PD%X(3)*JACOB%ENT(1,1))/J**2
! 	S(2,1) = (JXI*JACOB%ENT(2,1)*JACOB%ENT(2,2)-JETA*JACOB%ENT(2,1)*JACOB%ENT(2,1))/J**3-(SECOND_PD%Y(1)*JACOB%ENT(2,2)-SECOND_PD%Y(2)*JACOB%ENT(2,1))/J**2
! 	S(2,2) = -(JXI*JACOB%ENT(1,1)*JACOB%ENT(2,2)-JETA*JACOB%ENT(1,1)*JACOB%ENT(2,1))/J**3+(SECOND_PD%X(1)*JACOB%ENT(2,2)-SECOND_PD%X(2)*JACOB%ENT(2,1))/J**2
! 	S(2,3) = (JXI*JACOB%ENT(1,1)*JACOB%ENT(1,2)-JETA*JACOB%ENT(1,1)*JACOB%ENT(1,1))/J**3-(SECOND_PD%X(1)*JACOB%ENT(1,2)-SECOND_PD%X(2)*JACOB%ENT(1,1))/J**2
! 	
! 	DBSFUN(1) =  GET_ALL_DIFF_BSPLINES(REF_PT%X,BASIS_KVEC(1),DIFF_ORDER)
! 	DBSFUN(2) =  GET_ALL_DIFF_BSPLINES(REF_PT%Y,BASIS_KVEC(2),DIFF_ORDER)
! 	
! 	
! 	DO K = 0, DIFF_ORDER
! 	DO L = 0, DIFF_ORDER-K
! 		DWSUM(K,L) = 0.0D0
! 		DO J1 = 0, BASIS_KVEC(1)%POLY_ORDER
! 			INDX_X = DBSFUN(1)%INIT + J1
! 			IF (INDX_X>=0) THEN
! 				DO J2 = 0, BASIS_KVEC(2)%POLY_ORDER
! 					INDX_Y = DBSFUN(2)%INIT + J2
! 					IF (INDX_Y>=0) THEN
! 						DWSUM(K,L) = DWSUM(K,L) + DBSFUN(1)%N(J1,K)*DBSFUN(2)%N(J2,L)*BASIS_CTL%WGTS(INDX_X,INDX_Y,0)
! 					ENDIF
! 				ENDDO
! 			ENDIF
! 		ENDDO
! 	ENDDO
! 	ENDDO
! 
! 	LC_INDX = 0
! 	
! 	DO I1 = 0, BASIS_KVEC(1)%POLY_ORDER
! 		
! 		INDX_X = DBSFUN(1)%INIT + I1
! 
! 	DO I2 = 0, BASIS_KVEC(2)%POLY_ORDER
! 		
! 		INDX_Y = DBSFUN(2)%INIT + I2
! 		LC_INDX = LC_INDX + 1
! 		INDX(LC_INDX) = INT2D(INDX_X, INDX_Y)
! 		
! 		DO K = 0, DIFF_ORDER
! 		DO L = 0, DIFF_ORDER-K
! 		
! 			DIFF_NURBS%VAL(K,L) = DBSFUN(1)%N(I1,K)*DBSFUN(2)%N(I2,L)*BASIS_CTL%WGTS(INDX_X,INDX_Y,0)
! 			
! 			DO J1 = 1, K
! 				DIFF_NURBS%VAL(K,L) = DIFF_NURBS%VAL(K,L) - BINOM(K,J1)*DWSUM(J1,0)*DIFF_NURBS%VAL(K-J1,L)
! 			ENDDO
! 			DO J2 = 1, L
! 				DIFF_NURBS%VAL(K,L) = DIFF_NURBS%VAL(K,L) - BINOM(L,J2)*DWSUM(0,J2)*DIFF_NURBS%VAL(K,L-J2)
! 			ENDDO
! 			DO J1 = 1, K
! 			DO J2 = 1, L
! 				DIFF_NURBS%VAL(K,L) = DIFF_NURBS%VAL(K,L) - BINOM(K,J1)*BINOM(L,J2)*DWSUM(J1,J2)*DIFF_NURBS%VAL(K-J1,L-J2)
! 			ENDDO
! 			ENDDO
! 			DIFF_NURBS%VAL(K,L) = DIFF_NURBS%VAL(K,L) / DWSUM(0,0)
! 		ENDDO
! 		ENDDO
! 		
! 		SF(LC_INDX)%D00 = DIFF_NURBS%VAL(0,0)
! 		SF(LC_INDX)%D10 = DIFF_NURBS%VAL(1,0)*INV_JACOB%ENT(1,1) + DIFF_NURBS%VAL(0,1)*INV_JACOB%ENT(2,1)
! 		SF(LC_INDX)%D01 = DIFF_NURBS%VAL(1,0)*INV_JACOB%ENT(1,2) + DIFF_NURBS%VAL(0,1)*INV_JACOB%ENT(2,2)
! 		SF(LC_INDX)%D20 = DIFF_NURBS%VAL(2,0)*INV_JACOB%ENT(1,1)**2 + 2.0D0*DIFF_NURBS%VAL(1,1)*INV_JACOB%ENT(1,1)*INV_JACOB%ENT(2,1) + DIFF_NURBS%VAL(0,2)*INV_JACOB%ENT(2,1)**2 + DIFF_NURBS%VAL(1,0)*S(1,1) + DIFF_NURBS%VAL(0,1)*S(2,1)
! 		SF(LC_INDX)%D11 = DIFF_NURBS%VAL(2,0)*INV_JACOB%ENT(1,1)*INV_JACOB%ENT(1,2) + DIFF_NURBS%VAL(1,1)*(INV_JACOB%ENT(1,1)*INV_JACOB%ENT(2,2)+INV_JACOB%ENT(1,2)*INV_JACOB%ENT(2,1)) + DIFF_NURBS%VAL(0,2)*INV_JACOB%ENT(2,1)*INV_JACOB%ENT(2,2) + DIFF_NURBS%VAL(1,0)*S(1,2) + DIFF_NURBS%VAL(0,1)*S(2,2)
! 		SF(LC_INDX)%D02 = DIFF_NURBS%VAL(2,0)*INV_JACOB%ENT(1,2)**2 + 2.0D0*DIFF_NURBS%VAL(1,1)*INV_JACOB%ENT(1,2)*INV_JACOB%ENT(2,2) + DIFF_NURBS%VAL(0,2)*INV_JACOB%ENT(2,2)**2 + DIFF_NURBS%VAL(1,0)*S(1,3) + DIFF_NURBS%VAL(0,1)*S(2,3)
! 		
! 	ENDDO
! 	ENDDO
! 	
! END SUBROUTINE GET_ALL_PHY_NURBS_SURFACE_2D
! 
! 
! 
! !! FIND ALL NONVANISHING BASIS FUNCTION ON PHYSICAL SPACE
! SUBROUTINE GET_ALL_PHY_BASIS_IN_INT(SF, INDX, REF_PT, JACOB)
! 
! 	TYPE(FVALUE), INTENT(OUT) :: SF((BASIS_KVEC(1)%POLY_ORDER+1)*(BASIS_KVEC(2)%POLY_ORDER+1))
! 	TYPE(INT2D), INTENT(OUT) :: INDX((BASIS_KVEC(1)%POLY_ORDER+1)*(BASIS_KVEC(2)%POLY_ORDER+1))
! 	TYPE(POINT2D), INTENT(IN) :: REF_PT
! 	TYPE(MATRIX_22), INTENT(IN) :: JACOB
! 	
! 	REAL*8 :: DJ, DJXI, DJETA, STEPFT(3), INV_JXI(2,2), INV_JETA(2,2), INV_J(2,2), DUF(2), DUXF(2), DUYF(2), DUFXI(2), DUFETA(2)
! 	TYPE(SECOND_PARTIAL_DERIVATIVES) :: SECOND_PD
! 	TYPE(DIFF_BSPLINES) :: NX, NY
! 	TYPE(MATRIX_22) :: INV_JACOB
! 	INTEGER :: I, K, II, JJ, KK
! 
! 	REAL*8 :: JXI, JETA, S(2,3)
! 	
! 	INV_JACOB = .INVERSE.JACOB
! 	
! ! 	INV_J(1,1) = INV_JACOB%ENT(1,1)
! ! 	INV_J(1,2) = INV_JACOB%ENT(1,2)
! ! 	INV_J(2,1) = INV_JACOB%ENT(2,1)
! ! 	INV_J(2,2) = INV_JACOB%ENT(2,2)
! 	
! 	SECOND_PD = GET_PARTIAL(REF_PT)				!! P^2(X)/P(XI^2), P^2(X)/P(XI,ETA), P^2(X)/P(ETA^2), P^2(Y)/P(XI^2), P^2(Y)/P(XI,ETA), P^2(Y)/P(ETA^2)
! 
! !  	DJ = JACOB%ENT(1,1)*JACOB%ENT(2,2)-JACOB%ENT(1,2)*JACOB%ENT(2,1)
! ! 	DJXI = SECOND_PD%X(1)*JACOB%ENT(2,2)+JACOB%ENT(1,1)*SECOND_PD%Y(2)-SECOND_PD%y(1)*JACOB%ENT(2,1)-JACOB%ENT(1,2)*SECOND_PD%x(2)
! ! 	DJETA = SECOND_PD%X(2)*JACOB%ENT(2,2)+JACOB%ENT(1,1)*SECOND_PD%Y(3)-SECOND_PD%y(2)*JACOB%ENT(2,1)-JACOB%ENT(1,2)*SECOND_PD%x(3)
! ! 
! ! 	INV_JXI(1,1) =  (SECOND_PD%Y(2)*DJ - JACOB%ENT(2,2)*DJXI)/DJ**2
! ! 	INV_JXI(1,2) = -(SECOND_PD%Y(1)*DJ - JACOB%ENT(1,2)*DJXI)/DJ**2
! ! 	INV_JXI(2,1) = -(SECOND_PD%X(2)*DJ - JACOB%ENT(2,1)*DJXI)/DJ**2
! ! 	INV_JXI(2,2) =  (SECOND_PD%X(1)*DJ - JACOB%ENT(1,1)*DJXI)/DJ**2
! ! 	
! ! 	INV_JETA(1,1) = (SECOND_PD%Y(3)*DJ - JACOB%ENT(2,2)*DJETA)/DJ**2
! ! 	INV_JETA(1,2) = (SECOND_PD%Y(2)*DJ - JACOB%ENT(1,2)*DJETA)/DJ**2
! ! 	INV_JETA(2,1) = (SECOND_PD%X(3)*DJ - JACOB%ENT(2,1)*DJETA)/DJ**2
! ! 	INV_JETA(2,2) = (SECOND_PD%X(2)*DJ - JACOB%ENT(1,1)*DJETA)/DJ**2
! 
!  	DJ = JACOB%ENT(1,1)*JACOB%ENT(2,2)-JACOB%ENT(1,2)*JACOB%ENT(2,1)
! 	JXI = SECOND_PD%X(1)*JACOB%ENT(2,2)+JACOB%ENT(1,1)*SECOND_PD%Y(2)-SECOND_PD%X(2)*JACOB%ENT(2,1)-JACOB%ENT(1,2)*SECOND_PD%Y(1)
! 	JETA = SECOND_PD%X(2)*JACOB%ENT(2,2)+JACOB%ENT(1,1)*SECOND_PD%Y(3)-SECOND_PD%X(3)*JACOB%ENT(2,1)-JACOB%ENT(1,2)*SECOND_PD%Y(2)
! 
! 	S(1,1) = -(JXI*JACOB%ENT(2,2)*JACOB%ENT(2,2)-JETA*JACOB%ENT(2,1)*JACOB%ENT(2,2))/DJ**3+(SECOND_PD%Y(2)*JACOB%ENT(2,2)-SECOND_PD%Y(3)*JACOB%ENT(2,1))/DJ**2
! 	S(1,2) = (JXI*JACOB%ENT(1,2)*JACOB%ENT(2,2)-JETA*JACOB%ENT(1,2)*JACOB%ENT(2,1))/DJ**3-(SECOND_PD%X(2)*JACOB%ENT(2,2)-SECOND_PD%X(3)*JACOB%ENT(2,1))/DJ**2
! 	S(1,3) = -(JXI*JACOB%ENT(1,2)*JACOB%ENT(1,2)-JETA*JACOB%ENT(1,1)*JACOB%ENT(1,2))/DJ**3+(SECOND_PD%X(2)*JACOB%ENT(1,2)-SECOND_PD%X(3)*JACOB%ENT(1,1))/DJ**2
! 	S(2,1) = (JXI*JACOB%ENT(2,1)*JACOB%ENT(2,2)-JETA*JACOB%ENT(2,1)*JACOB%ENT(2,1))/DJ**3-(SECOND_PD%Y(1)*JACOB%ENT(2,2)-SECOND_PD%Y(2)*JACOB%ENT(2,1))/DJ**2
! 	S(2,2) = -(JXI*JACOB%ENT(1,1)*JACOB%ENT(2,2)-JETA*JACOB%ENT(1,1)*JACOB%ENT(2,1))/DJ**3+(SECOND_PD%X(1)*JACOB%ENT(2,2)-SECOND_PD%X(2)*JACOB%ENT(2,1))/DJ**2
! 	S(2,3) = (JXI*JACOB%ENT(1,1)*JACOB%ENT(1,2)-JETA*JACOB%ENT(1,1)*JACOB%ENT(1,1))/DJ**3-(SECOND_PD%X(1)*JACOB%ENT(1,2)-SECOND_PD%X(2)*JACOB%ENT(1,1))/DJ**2
! 
! 	NX = GET_ALL_DIFF_BSPLINES(REF_PT%X, BASIS_KVEC(1), 2)
! 	NY = GET_ALL_DIFF_BSPLINES(REF_PT%Y, BASIS_KVEC(2), 2)
! 	
! ! 	IF (DABS(REF_PT%X-1.0D0)<=EPS) THEN
! ! 		PRINT*, '---------NX-16TH----------'
! ! 		PRINT*, NX%N(7,:)
! ! 		PRINT*, '---------NX-17TH----------'
! ! 		PRINT*, NX%N(8,:)
! ! 		STOP
! ! 	ENDIF
! 
! ! 	IF (DABS(REF_PT%Y-1.0D0)<=EPS) THEN
! ! 		PRINT*, '---------NY-3TH----------'
! ! 		PRINT*, NY%N(2,:)
! ! 		PRINT*, '---------NY-4TH----------'
! ! 		PRINT*, NY%N(3,:)
! ! 		STOP
! ! 	ENDIF
! 
! 	K = 0
! 	DO J=0, NY%POLY_ORDER
! 		DO I=0, NX%POLY_ORDER
! 			K = K + 1
! 			INDX(K) = INT2D(NX%INIT+I, NY%INIT+J)
! 			
! ! 			DUF = MATMUL(INV_J, (/NX%N(I,1)*NY%N(J,0), NX%N(I,0)*NY%N(J,1)/))
! ! 			
! ! 			DUFXI(1) = (INV_JXI(1,1)*NX%N(I,1)*NY%N(J,0) + INV_J(1,1)*NX%N(I,2)*NY%N(J,0)) + (INV_JXI(1,2)*NX%N(I,0)*NY%N(J,1) + INV_J(1,2)*NX%N(I,1)*NY%N(J,1))
! ! 			DUFXI(2) = (INV_JETA(1,1)*NX%N(I,1)*NY%N(J,0) + INV_J(1,1)*NX%N(I,1)*NY%N(J,1)) + (INV_JETA(1,2)*NX%N(I,0)*NY%N(J,1) + INV_J(1,2)*NX%N(I,0)*NY%N(J,2))
! ! 			
! ! 			DUFETA(1) = (INV_JXI(2,1)*NX%N(I,1)*NY%N(J,0) + INV_J(2,1)*NX%N(I,2)*NY%N(J,0)) + (INV_JXI(2,2)*NX%N(I,0)*NY%N(J,1) + INV_J(2,2)*NX%N(I,1)*NY%N(J,1))
! ! 			DUFETA(2) = (INV_JETA(2,1)*NX%N(I,1)*NY%N(J,0) + INV_J(2,1)*NX%N(I,1)*NY%N(J,1)) + (INV_JETA(2,2)*NX%N(I,0)*NY%N(J,1) + INV_J(2,2)*NX%N(I,0)*NY%N(J,2))
! ! 			
! ! 			DUXF = MATMUL(INV_J, DUFXI)
! ! 			DUYF = MATMUL(INV_J, DUFETA)
! ! 			
! ! 			SF(K)%D00 = NX%N(I,0)*NY%N(J,0)
! ! 			SF(K)%D10 = DUF(1)
! ! 			SF(K)%D01 = DUF(2)
! ! 			SF(K)%D20 = DUXF(1)
! ! 			SF(K)%D11 = DUXF(2)
! ! 			SF(K)%D02 = DUYF(2)
! 			
! 			SF(K)%D00 = NX%N(I,0)*NY%N(J,0)
! 			SF(K)%D10 = NX%N(I,1)*NY%N(J,0)*INV_JACOB%ENT(1,1) + NX%N(I,0)*NY%N(J,1)*INV_JACOB%ENT(2,1)
! 			SF(K)%D01 = NX%N(I,1)*NY%N(J,0)*INV_JACOB%ENT(1,2) + NX%N(I,0)*NY%N(J,1)*INV_JACOB%ENT(2,2)
! 			SF(K)%D20 = NX%N(I,2)*NY%N(J,0)*INV_JACOB%ENT(1,1)**2 + 2.0D0*NX%N(I,1)*NY%N(J,1)*INV_JACOB%ENT(1,1)*INV_JACOB%ENT(2,1) + NX%N(I,0)*NY%N(J,2)*INV_JACOB%ENT(2,1)**2 + NX%N(I,1)*NY%N(J,0)*S(1,1) + NX%N(I,0)*NY%N(J,1)*S(2,1)
! 			SF(K)%D11 = NX%N(I,2)*NY%N(J,0)*INV_JACOB%ENT(1,1)*INV_JACOB%ENT(1,2) + NX%N(I,1)*NY%N(J,1)*(INV_JACOB%ENT(1,1)*INV_JACOB%ENT(2,2)+INV_JACOB%ENT(1,2)*INV_JACOB%ENT(2,1)) + NX%N(I,0)*NY%N(J,2)*INV_JACOB%ENT(2,1)*INV_JACOB%ENT(2,2) + NX%N(I,1)*NY%N(J,0)*S(1,2) + NX%N(I,0)*NY%N(J,1)*S(2,2)
! 			SF(K)%D02 = NX%N(I,2)*NY%N(J,0)*INV_JACOB%ENT(1,2)**2 + 2.0D0*NX%N(I,1)*NY%N(J,1)*INV_JACOB%ENT(1,2)*INV_JACOB%ENT(2,2) + NX%N(I,0)*NY%N(J,2)*INV_JACOB%ENT(2,2)**2 + NX%N(I,1)*NY%N(J,0)*S(1,3) + NX%N(I,0)*NY%N(J,1)*S(2,3)
! 
! 		ENDDO
! 	ENDDO
! ! 	PRINT*, INDX(1),SF(1)
! ! 	STOP
! END SUBROUTINE GET_ALL_PHY_BASIS_IN_INT
! 
! 
! 
! 
! 
! 
! TYPE(FVALUE) FUNCTION GET_PHY_NURBS_SURFACE_2D(REF_PT, IX, IY, JACOB)
! 
! 	INTEGER, INTENT(IN) :: IX, IY
! 	TYPE(POINT2D), INTENT(IN) :: REF_PT
! 	TYPE(MATRIX_22), INTENT(IN) :: JACOB
! 
! 	TYPE(MATRIX_22) :: INV_JACOB
! 	TYPE(FUNCTION_2D) :: REF_SF, DIFF_NURBS
! 	TYPE(DIFF_BSPLINES) :: DBSFUN(2)
! 	REAL(8) :: DWSUM(0:MAX_DIFF_ORDER,0:MAX_DIFF_ORDER)
! 	INTEGER :: INDX_X, INDX_Y, DIFF_ORDER
! 	INTEGER :: I, J, K, L, J1, J2, I1, I2, LC_INDX
! 
! 	DIFF_ORDER = 1
! 
! 	INV_JACOB = .INVERSE.JACOB
! 	
! 	DBSFUN(1) =  GET_ALL_DIFF_BSPLINES(REF_PT%X,BASIS_KVEC(1),DIFF_ORDER)
! 	DBSFUN(2) =  GET_ALL_DIFF_BSPLINES(REF_PT%Y,BASIS_KVEC(2),DIFF_ORDER)
! 
! 	DO K = 0, DIFF_ORDER
! 	DO L = 0, DIFF_ORDER-K
! 		DWSUM(K,L) = 0.0D0
! 		DO J1 = 0, BASIS_KVEC(1)%POLY_ORDER
! 			INDX_X = DBSFUN(1)%INIT + J1
! 			IF (INDX_X>=0) THEN
! 				DO J2 = 0, BASIS_KVEC(2)%POLY_ORDER
! 					INDX_Y = DBSFUN(2)%INIT + J2
! 					IF (INDX_Y>=0) THEN
! 						DWSUM(K,L) = DWSUM(K,L) + DBSFUN(1)%N(J1,K)*DBSFUN(2)%N(J2,L)*BASIS_CTL%WGTS(INDX_X,INDX_Y,0)
! 					ENDIF
! 				ENDDO
! 			ENDIF
! 		ENDDO
! 	ENDDO
! 	ENDDO
! 
! 	REF_SF = GET_BASIS(REF_PT,(/IX,IY/))
! 		
! 	DO K = 0, DIFF_ORDER
! 	DO L = 0, DIFF_ORDER-K
! 	
! 		DIFF_NURBS%VAL(K,L) = REF_SF%VAL(K,L)*BASIS_CTL%WGTS(IX,IY,0)
! 		
! 		DO J1 = 1, K
! 			DIFF_NURBS%VAL(K,L) = DIFF_NURBS%VAL(K,L) - BINOM(K,J1)*DWSUM(J1,0)*DIFF_NURBS%VAL(K-J1,L)
! 		ENDDO
! 		DO J2 = 1, L
! 			DIFF_NURBS%VAL(K,L) = DIFF_NURBS%VAL(K,L) - BINOM(L,J2)*DWSUM(0,J2)*DIFF_NURBS%VAL(K,L-J2)
! 		ENDDO
! 		DO J1 = 1, K
! 		DO J2 = 1, L
! 			DIFF_NURBS%VAL(K,L) = DIFF_NURBS%VAL(K,L) - BINOM(K,J1)*BINOM(L,J2)*DWSUM(J1,J2)*DIFF_NURBS%VAL(K-J1,L-J2)
! 		ENDDO
! 		ENDDO
! 
! 		DIFF_NURBS%VAL(K,L) = DIFF_NURBS%VAL(K,L) / DWSUM(0,0)
! 
! 	ENDDO
! 	ENDDO
! 	
! 	GET_PHY_NURBS_SURFACE_2D%D00 = DIFF_NURBS%VAL(0,0)
! 	GET_PHY_NURBS_SURFACE_2D%D10 = DIFF_NURBS%VAL(1,0)*INV_JACOB%ENT(1,1) + DIFF_NURBS%VAL(0,1)*INV_JACOB%ENT(2,1)
! 	GET_PHY_NURBS_SURFACE_2D%D01 = DIFF_NURBS%VAL(1,0)*INV_JACOB%ENT(1,2) + DIFF_NURBS%VAL(0,1)*INV_JACOB%ENT(2,2)
! 		
! END FUNCTION GET_PHY_NURBS_SURFACE_2D
! 
! !!  FIND BASIS FUNCTION ON PARAMETRIC SPACE
! TYPE(FUNCTION_2D) FUNCTION GET_BASIS(REF_PT,INDX)
! 
! 	TYPE(POINT2D), INTENT(IN) :: REF_PT
! 	INTEGER, INTENT(IN) :: INDX(2)
! 
! 	TYPE(FUNCTION_1D) :: NX, NY
! 
! 	NX = GET_DIFF_BSPLINE(REF_PT%X,BASIS_KVEC(1),INDX(1),2)
! 	NY = GET_DIFF_BSPLINE(REF_PT%Y,BASIS_KVEC(2),INDX(2),2)
! 
! 	GET_BASIS%VAL(0,0) = NX%VAL(0) * NY%VAL(0)
! 	GET_BASIS%VAL(1,0) = NX%VAL(1) * NY%VAL(0)
! 	GET_BASIS%VAL(0,1) = NX%VAL(0) * NY%VAL(1)
! 	GET_BASIS%VAL(2,0) = NX%VAL(2) * NY%VAL(0)
! 	GET_BASIS%VAL(1,1) = NX%VAL(1) * NY%VAL(1)
! 	GET_BASIS%VAL(0,2) = NX%VAL(0) * NY%VAL(2)
! 
! END FUNCTION GET_BASIS

END MODULE NURBS_BASIS