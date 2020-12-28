MODULE PLOT

	USE PATCH_MAPPING
	USE GEOMETRY
	USE GSQUAD

	implicit integer (i-n)
	implicit real(8) (a-h,o-z)

CONTAINS

SUBROUTINE PLOTGM()
	
	TYPE(RECPATCH) :: IRBOX
	TYPE(TRANSFORM2D) :: GSPT
	TYPE(POINT2D) :: TSPT, NURBS_PT, NURBS_TSPT
 	TYPE(FUNCTION_1D) :: TS_BS
	INTEGER :: I, J, II, JJ, KK, NI, NJ, NM, PATCH_ROW, PATCH_COLUMN

	OPEN(1, FILE = './data/irbox')
! 	OPEN(2, FILE = './data/phyirbox')
	DO PATCH_ROW = 1, NUMPATCH(1)
		DO PATCH_COLUMN = 1, NUMPATCH(2)
			DO JJ=1, NUMIR(PATCH_ROW, PATCH_COLUMN, 2)-1
				DO II=1, NUMIR(PATCH_ROW, PATCH_COLUMN, 1)-1
					IRBOX%PT1 = POINT2D(IR_GRID(PATCH_ROW, PATCH_COLUMN, 1,II), IR_GRID(PATCH_ROW, PATCH_COLUMN, 2,JJ))
					IRBOX%PT2 = POINT2D(IR_GRID(PATCH_ROW, PATCH_COLUMN, 1,II+1), IR_GRID(PATCH_ROW, PATCH_COLUMN, 2,JJ))
					IRBOX%PT3 = POINT2D(IR_GRID(PATCH_ROW, PATCH_COLUMN, 1,II+1), IR_GRID(PATCH_ROW, PATCH_COLUMN, 2,JJ+1))
					IRBOX%PT4 = POINT2D(IR_GRID(PATCH_ROW, PATCH_COLUMN, 1,II), IR_GRID(PATCH_ROW, PATCH_COLUMN, 2,JJ+1))
					WRITE(1, *) IRBOX%PT1
					WRITE(1, *) IRBOX%PT2
					WRITE(1, *) IRBOX%PT3
					WRITE(1, *) IRBOX%PT4
					WRITE(1, *) IRBOX%PT1
					WRITE(1,*)
! 					WRITE(2, *) GET_PHY_PT(IRBOX%PT1)
! 					WRITE(2, *) GET_PHY_PT(IRBOX%PT2)
! 					WRITE(2, *) GET_PHY_PT(IRBOX%PT3)
! 					WRITE(2, *) GET_PHY_PT(IRBOX%PT4)
! 					WRITE(2, *) GET_PHY_PT(IRBOX%PT1)
! 					WRITE(2,*)
				ENDDO
			ENDDO
		ENDDO
	ENDDO
	CLOSE(1)
! 	CLOSE(2)
	
	OPEN(1, FILE = './data/l2box')
	DO I = 1, L2_NUMIR(1)
		DO J = 1, L2_NUMIR(2)
			IRBOX%PT1 = POINT2D(L2_GRID(1,II), L2_GRID(2,JJ))
			IRBOX%PT2 = POINT2D(L2_GRID(1,II+1), L2_GRID(2,JJ))
			IRBOX%PT3 = POINT2D(L2_GRID(1,II+1), L2_GRID(2,JJ+1))
			IRBOX%PT4 = POINT2D(L2_GRID(1,II), L2_GRID(2,JJ+1))
			WRITE(1, *) IRBOX%PT1
			WRITE(1, *) IRBOX%PT2
			WRITE(1, *) IRBOX%PT3
			WRITE(1, *) IRBOX%PT4
			WRITE(1, *) IRBOX%PT1
			WRITE(1,*)
		ENDDO
	ENDDO
	CLOSE(1)
	
	OPEN(1, FILE = './data/irboxx')
	DO PATCH_COLUMN = 1, NUMPATCH(2)
		WRITE(1,*) 'PATCH_COLUMN: ', PATCH_COLUMN
		WRITE(1,*) (IR_GRID(1, PATCH_COLUMN, 1, II), II = 1, NUMIR(1, PATCH_COLUMN, 1))
	ENDDO
	CLOSE(1)
	
	OPEN(1, FILE = './data/irboxy')
	DO PATCH_COLUMN = 1, NUMPATCH(2)
		WRITE(1,*) 'PATCH_COLUMN: ', PATCH_COLUMN
		WRITE(1,*) (IR_GRID(1, PATCH_COLUMN, 2, II), II = 1, NUMIR(1, PATCH_COLUMN, 2))
	ENDDO
	CLOSE(1)
	
	OPEN(1, FILE = './data/loc_numbs')
	DO I = 1, NUMPATCH(1)
		DO J = 1, NUMPATCH(2)
			WRITE(1,*) LOC_NUMBS(I,J,1:2)
		ENDDO
	ENDDO
	CLOSE(1)
	
	OPEN(1, FILE = './data/ndx')
	DO I=1, DOF
		WRITE(1,*) NDX(:,I)
	ENDDO
	CLOSE(1)
	
	OPEN(1, FILE = './data/zerobdndx')
	DO I=1, ZERO_BDNDX(1)%LC_NUM
		WRITE(1,*) ZERO_BDNDX(I)%LC_NDX(:), ZERO_BDNDX(I)%LC_NUM, ZERO_BDNDX(I)%GL_NDX
	ENDDO
	CLOSE(1)
	
	OPEN(1, FILE = './data/basis_kvec1')
	DO PATCH_COLUMN = 1, NUMPATCH(2)
		WRITE(1,*) 'PATCH_COLUMN: ', PATCH_COLUMN
		DO I = 0, BASIS_KVEC(1,PATCH_COLUMN,1)%LENGTH
			WRITE(1,*) BASIS_KVEC(1,PATCH_COLUMN,1)%KNOTS(I)
		ENDDO
	ENDDO
	CLOSE(1)
	
	OPEN(1, FILE = './data/basis_kvec2')
	DO PATCH_COLUMN = 1, NUMPATCH(2)
		WRITE(1,*) 'PATCH_COLUMN: ', PATCH_COLUMN
		DO J = 0, BASIS_KVEC(1,PATCH_COLUMN,2)%LENGTH
			WRITE(1,*) BASIS_KVEC(1,PATCH_COLUMN,2)%KNOTS(J)
		ENDDO
	ENDDO
	CLOSE(1)
	
	OPEN(1, FILE = './data/bdndx')
	DO I=1, BDNDX(1)%LC_NUM
		WRITE(1,*) BDNDX(I)%LC_NDX(:), BDNDX(I)%LC_NUM, BDNDX(I)%GL_NDX
	ENDDO
	CLOSE(1)
	
! 	OPEN(1, FILE = './data/basis_knot')
! 	WRITE(1,*) 'POLY_ORDER : ', BASIS_KVEC(1)%POLY_ORDER
! 	WRITE(1,*) 'LENGTH : ', BASIS_KVEC(1)%LENGTH
! 	DO I=0, BASIS_KVEC(1)%LENGTH
! 		WRITE(1,102) BASIS_KVEC(1)%KNOTS(I)
! 	ENDDO
! 	WRITE(1,*)
! 	WRITE(1,*) 'POLY_ORDER : ', BASIS_KVEC(2)%POLY_ORDER
! 	WRITE(1,*) 'LENGTH : ', BASIS_KVEC(2)%LENGTH
! 	DO I=0, BASIS_KVEC(2)%LENGTH
! 		WRITE(1,102) BASIS_KVEC(2)%KNOTS(I)
! 	ENDDO
! 	CLOSE(1)

! 	OPEN(1, FILE = './data/geo_basis_knot')
! 	WRITE(1,*) 'POLY_ORDER : ', GEO_KVEC(1)%POLY_ORDER
! 	WRITE(1,*) 'LENGTH : ', GEO_KVEC(1)%LENGTH
! 	DO I=0, GEO_KVEC(1)%LENGTH
! 		WRITE(1,102) GEO_KVEC(1)%KNOTS(I)
! 	ENDDO
! 	WRITE(1,*)
! 	WRITE(1,*) 'POLY_ORDER : ', GEO_KVEC(2)%POLY_ORDER
! 	WRITE(1,*) 'LENGTH : ', GEO_KVEC(2)%LENGTH
! 	DO I=0, GEO_KVEC(2)%LENGTH
! 		WRITE(1,102) GEO_KVEC(2)%KNOTS(I)
! 	ENDDO
! 	CLOSE(1)

! 	OPEN(11, FILE = './data/ctl')
! 	DO I=0, GEO_KVEC(2)%LENGTH - GEO_KVEC(2)%POLY_ORDER - 1
! 		DO J=0, GEO_KVEC(1)%LENGTH - GEO_KVEC(1)%POLY_ORDER - 1
! 			WRITE(11,*) GEO_CTL%PTS(J,I,0), GEO_CTL%WGTS(J,I,0)
! 		ENDDO
! 		WRITE(11,*)
! 	ENDDO
! 	DO J=0, GEO_KVEC(1)%LENGTH - GEO_KVEC(1)%POLY_ORDER - 1
! 		DO I=0, GEO_KVEC(2)%LENGTH - GEO_KVEC(2)%POLY_ORDER - 1
! 			WRITE(11,*) GEO_CTL%PTS(J,I,0), GEO_CTL%WGTS(J,I,0)
! 		ENDDO
! 		WRITE(11,*)
! 	ENDDO
! 	CLOSE(11)
! 
! 	OPEN(11, FILE = './data/basis_ctl')
! 	DO I=0, NUMBS(2)-1
! 		DO J=0, NUMBS(1)-1
! 			WRITE(11,*) BASIS_CTL%PTS(J,I,0), BASIS_CTL%WGTS(J,I,0)
! 		ENDDO
! 		WRITE(11,*)
! 	ENDDO
! 	CLOSE(11)
	
	OPEN(1, FILE = './data/bsplinex3')
	DO J = 0, LOC_NUMBS(1,3,1)-1
		DO I=1, 601
			TSPT%X = 1.D0*(I-1)/600.D0
			TS_BS = GET_DIFF_BSPLINE(TSPT%X, BASIS_KVEC(1,3,1), J, 1)
			WRITE(1,*) TSPT%X, TS_BS%VAL(:)
		ENDDO
		WRITE(1,*)
	ENDDO
	CLOSE(1)
! 	
! 	OPEN(2, FILE = './data/bspliney')
! 	DO J=0, NUMBS(2)-1
! 		DO I=1, 601
! 			TSPT%Y = 1.D0*(I-1)/600.D0
! 			TS_BS = GET_DIFF_BSPLINE(TSPT%Y,BASIS_KVEC(2),J,1)
! 			WRITE(2,*) TSPT%Y, TS_BS%VAL(:)
! 		ENDDO
! 		WRITE(2,*)
! 	ENDDO
! 	CLOSE(2)
! 	
! 	OPEN(1, FILE = './data/geometry')
! 	DO I=1, 101
!  		DO J=1, 101
!  			TSPT%Y = 1.D0*(I-1)/100.D0
! 			TSPT%X = 1.D0*(J-1)/100.D0
! 			NURBS_PT = GET_POINT_NURVE_SURFACE_2D(TSPT,GEO_KVEC,GEO_CTL)
! 			WRITE(1,*) NURBS_PT
! 		ENDDO
! 		WRITE(1,*) 
! 	ENDDO
! 	CLOSE(1)
! 
! 	OPEN(1, FILE = './data/geometry_knot')
! 	DO J = 1, 2
! 		DO I = 1, NUMIR(J)
! 			DO II = 1, 101
! 				IF (J==1) THEN
! 					TSPT = POINT2D(IR_GRID(1, I), 1.0D0*(II-1)/100.0D0)
! 				ELSE
! 					TSPT = POINT2D(1.0D0*(II-1)/100.0D0, IR_GRID(2, I))
! 				ENDIF
! 				NURBS_PT = GET_POINT_NURVE_SURFACE_2D(TSPT,BASIS_KVEC,BASIS_CTL)
! 				WRITE(1,*) NURBS_PT
! 			ENDDO
! 			WRITE(1,*) ''
! 		ENDDO
! 		WRITE(1,*) ''
! 	ENDDO
! 	CLOSE(1)
	
! 	OPEN(1, FILE = './data/basis_nurbs')
! 	DO I=1, 101
!  		DO J=1, 101
!  			TSPT%Y = 1.D0*(I-1)/100.D0
! 			TSPT%X = 1.D0*(J-1)/100.D0
! 			NURBS_PT = GET_POINT_NURVE_SURFACE_2D(TSPT,BASIS_KVEC,BASIS_CTL)
! 			WRITE(1,*) NURBS_PT
! 		ENDDO
! 		WRITE(1,*) 
! 	ENDDO
! 	CLOSE(1)
! 	
! 	OPEN(1, FILE = './data/geometry_bd')
! 	DO J=1, 101
! 		TSPT = POINT2D(1.D0*(J-1)/100.D0, 0.D0)
! 		NURBS_PT = GET_POINT_NURVE_SURFACE_2D(TSPT,GEO_KVEC,GEO_CTL)
! 		WRITE(1,*) NURBS_PT
! 	ENDDO
! 	WRITE(1,*)
! 	DO J=1, 101
! 		TSPT = POINT2D(1.D0, 1.D0*(J-1)/100.D0)
! 		NURBS_PT = GET_POINT_NURVE_SURFACE_2D(TSPT,GEO_KVEC,GEO_CTL)
! 		WRITE(1,*) NURBS_PT
! 	ENDDO
! 	WRITE(1,*)
! 	DO J=1, 101
! 		TSPT = POINT2D(1.D0*(J-1)/100.D0, 1.D0)
! 		NURBS_PT = GET_POINT_NURVE_SURFACE_2D(TSPT,GEO_KVEC,GEO_CTL)
! 		WRITE(1,*) NURBS_PT
! 	ENDDO
! 	WRITE(1,*)
! 	DO J=1, 101
! 		TSPT = POINT2D(0.D0, 1.D0*(J-1)/100.D0)
! 		NURBS_PT = GET_POINT_NURVE_SURFACE_2D(TSPT,GEO_KVEC,GEO_CTL)
! 		WRITE(1,*) NURBS_PT
! 	ENDDO
! 	CLOSE(1)

	102 FORMAT(1000(f20.8,1x))
END SUBROUTINE PLOTGM

END MODULE PLOT
