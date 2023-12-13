C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
C     This file is part of FEDEM - https://openfedem.org
C
C>    @file beamaux.f
C>    @brief Auxilliary subroutines for beam finite elements.
C
C>    @brief Prints the rectangular matrix @a A(MA,NA) to unit @a IW.
      SUBROUTINE MPRT30(A,MA,NA,IW)
C
C **********************************************************************
C
C   FINITE ELEMENT LIBRARY :  MPRT30                       SINTEF / NTH
C
C
C     MPRT30 PRINTS THE MATRIX  A(MA,NA)
C
C     PROGRAMMED BY : K.BELL
C
C     DATE/VERSION  :  05.12.74 / 03
C
C **********************************************************************
C
      IMPLICIT NONE
C
      INTEGER           MA,NA,IW
      DOUBLE PRECISION  A(MA,NA)
      INTEGER           I,J,K,L
C
      DO 40 J=1,NA,8
      K = J+7
      IF(K-NA)20,20,10
   10 K = NA
   20 WRITE(IW,6000) (L,L=J,K)
      DO 30 I=1,MA
   30 WRITE(IW,6010) I,(A(I,L),L=J,K)
   40 CONTINUE
C
      RETURN
C
 6000 FORMAT(/8I15)
 6010 FORMAT(I5,1P,8E15.7)
C
      END
C>    @brief Computes a local-to-global transformation matrix.
      SUBROUTINE DCOS30(C,X,Y,Z)
C
C **********************************************************************
C
C   FINITE ELEMENT LIBRARY :  DCOS30                       SINTEF / NTH
C
C
C     DCOS30 COMPUTES THE 3*3 MATRIX OF DIRECTION COSINES WHICH
C     TRANSFORMS GLOBAL DISPLACEMENT COMPONENTS TO LOCAL DISPLACEMENT
C     COMPONENTS.
C
C     PROGRAMMED BY : F.STENSRUD
C
C     DATE/VERSION  : 25.05.73 / 01
C
C **********************************************************************
C
      IMPLICIT NONE
C
      DOUBLE PRECISION  C(3,3),X(3),Y(3),Z(3)
      DOUBLE PRECISION  AB,CX,CY,CZ
C
C     THE DIRECTION COSINES OF THE LOCAL X-AXIS
C
      CX = X(2) - X(1)
      CY = Y(2) - Y(1)
      CZ = Z(2) - Z(1)
      AB = CX*CX + CY*CY + CZ*CZ
      AB = DSQRT(AB)
      C(1,1) = CX/AB
      C(1,2) = CY/AB
      C(1,3) = CZ/AB
C
C     THE DIRECTION COSINES OF THE LOCAL Y-AXIS
C
      CX = C(1,3)*(Y(3) - Y(1)) - C(1,2)*(Z(3) - Z(1))
      CY = C(1,1)*(Z(3) - Z(1)) - C(1,3)*(X(3) - X(1))
      CZ = C(1,2)*(X(3) - X(1)) - C(1,1)*(Y(3) - Y(1))
      AB = CX*CX + CY*CY + CZ*CZ
      AB = DSQRT(AB)
      C(2,1) = CX/AB
      C(2,2) = CY/AB
      C(2,3) = CZ/AB
C
C     THE DIRECTION COSINES OF THE LOCAL Z-AXIS
C
      CX = C(1,2)*C(2,3) - C(1,3)*C(2,2)
      CY = C(1,3)*C(2,1) - C(1,1)*C(2,3)
      CZ = C(1,1)*C(2,2) - C(1,2)*C(2,1)
      AB = CX*CX + CY*CY + CZ*CZ
      AB = DSQRT(AB)
      C(3,1) = CX/AB
      C(3,2) = CY/AB
      C(3,3) = CZ/AB
C
      RETURN
      END
C>    @brief Performs a congruence transformation of a symmetric matrix.
      SUBROUTINE MPRO30(EK,C,M)
C
C **********************************************************************
C
C   FINITE ELEMENT LIBRARY :  MPRO30                       SINTEF / NTH
C
C
C     MPRO30 PERFORMS A CONGRUENCE TRANSFORMATION OF A SYMMETRIC
C     MATRIX WITH DIMENSION M*M WHERE M IS A MULTIPLE OF 3.
C     THE TRANSFORMATION MATRIX CONSISTS OF A 3*3 SUBMATRIX WHICH IS
C     REPEATED ALONG THE DIAGONAL.
C
C     PROGRAMMED BY : F.STENSRUD
C
C     DATE/VERSION  : 25.05.73 / 01
C
C **********************************************************************
C
      IMPLICIT NONE
C
      INTEGER           M
      DOUBLE PRECISION  EK(M,M),C(3,3),B(3,3)
      INTEGER           I,IA,IB,J,JA,JB,K,KA,L
C
      DO 500 I=0,M-3,3
      DO 500 J=0,M-3,3
C
      DO 100 K=1,3
      DO 100 L=1,3
  100 B(K,L) = 0.0D0
C
      DO 200 IA=1,3
      IB = I + IA
      DO 200 KA=1,3
      JB = J + KA
      DO 200 JA=1,3
  200 B(JA,IA) = B(JA,IA) + C(KA,JA)*EK(JB,IB)
C
      DO 300 IA=1,3
      IB = I + IA
      DO 300 KA=1,3
      JB = J + KA
  300 EK(JB,IB) = 0.0D0
C
      DO 400 IA=1,3
      IB = I + IA
      DO 400 KA=1,3
      JB = J + KA
      DO 400 JA=1,3
  400 EK(JB,IB) = EK(JB,IB) + B(KA,JA)*C(JA,IA)
  500 CONTINUE
C
      DO 1000 I=1,M
      DO 1000 J=I,M
 1000 EK(I,J) = EK(J,I)
C
      RETURN
      END
C>    @brief Transforms an element matrix from end points to nodal points.
      SUBROUTINE TRIX30(A,X1,Y1,Z1,X2,Y2,Z2)
C
C **********************************************************************
C
C   FINITE ELEMENT LIBRARY :  TRIX30                       SINTEF / NTH
C
C
C     TRIX30 TRANSFORM THE STIFFNESS (OR MASS) MATRIX OF A 12-DEGREE-
C     OF-FREEDOM BEAM ELEMENT FROM END POINTS TO ECCENTRICALLY LOCATED
C     NODAL POINTS
C
C     PROGRAMMED BY :  K.BELL
C
C     DATE/VERSION  :  24.02.77 / 01
C
C **********************************************************************
C
      IMPLICIT NONE
C
      DOUBLE PRECISION  A(12,12)
      DOUBLE PRECISION  X1,X2,Y1,Y2,Z1,Z2
      INTEGER           I
C
      DO 25 I=1,12
      A( 4,I) = A( 4,I) + Z1*A(2,I) - Y1*A(3,I)
      A( 5,I) = A( 5,I) - Z1*A(1,I) + X1*A(3,I)
      A( 6,I) = A( 6,I) + Y1*A(1,I) - X1*A(2,I)
      A(10,I) = A(10,I) + Z2*A(8,I) - Y2*A(9,I)
      A(11,I) = A(11,I) - Z2*A(7,I) + X2*A(9,I)
      A(12,I) = A(12,I) + Y2*A(7,I) - X2*A(8,I)
   25 CONTINUE
C
      DO 50 I=1,12
      A(I, 4) = A(I, 4) + Z1*A(I,2) - Y1*A(I,3)
      A(I, 5) = A(I, 5) - Z1*A(I,1) + X1*A(I,3)
      A(I, 6) = A(I, 6) + Y1*A(I,1) - X1*A(I,2)
      A(I,10) = A(I,10) + Z2*A(I,8) - Y2*A(I,9)
      A(I,11) = A(I,11) - Z2*A(I,7) + X2*A(I,9)
      A(I,12) = A(I,12) + Y2*A(I,7) - X2*A(I,8)
   50 CONTINUE
C
      RETURN
      END
C>    @brief Computes a local-to-global transformation matrix.
      SUBROUTINE DIRC30(C,X,Y,Z)
C
C **********************************************************************
C
C   FINITE ELEMENT LIBRARY :   DIRC30                      SINTEF / NTH
C
C
C     DIRC30 DETERMINES THE (3*3) MATRIX OF DIRECTION COSINES RELATING
C     GLOBAL TO LOCAL COORDINATES
C
C     PROGRAMMED BY :  K.BELL
C
C     DATE/VERSION  :  10.02.75 / 01
C
C **********************************************************************
C
      IMPLICIT NONE
C
      DOUBLE PRECISION  C(3,3),X(3),Y(3),Z(3)
      DOUBLE PRECISION  AB,CX,CY,CZ
C
C     THE DIRECTION COSINES OF THE LOCAL X-AXIS
C
      CX = X(2)-X(1)
      CY = Y(2)-Y(1)
      CZ = Z(2)-Z(1)
      AB = CX*CX + CY*CY + CZ*CZ
      AB = DSQRT(AB)
      C(1,1) = CX/AB
      C(1,2) = CY/AB
      C(1,3) = CZ/AB
C
C     THE DIRECTION COSINES OF THE LOCAL Z-AXIS
C
      CX = C(1,2)*(Z(3)-Z(1)) - C(1,3)*(Y(3)-Y(1))
      CY = C(1,3)*(X(3)-X(1)) - C(1,1)*(Z(3)-Z(1))
      CZ = C(1,1)*(Y(3)-Y(1)) - C(1,2)*(X(3)-X(1))
      AB = CX*CX + CY*CY + CZ*CZ
      AB = DSQRT(AB)
      C(3,1) = CX/AB
      C(3,2) = CY/AB
      C(3,3) = CZ/AB
C
C     THE DIRECTION COSINES OF THE LOCAL Y-AXIS
C
      CX = C(3,2)*C(1,3) - C(3,3)*C(1,2)
      CY = C(3,3)*C(1,1) - C(3,1)*C(1,3)
      CZ = C(3,1)*C(1,2) - C(3,2)*C(1,1)
      AB = CX*CX + CY*CY + CZ*CZ
      AB = DSQRT(AB)
      C(2,1) = CX/AB
      C(2,2) = CY/AB
      C(2,3) = CZ/AB
C
      RETURN
      END
