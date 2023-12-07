C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
C     This file is part of FEDEM - https://openfedem.org
C
      SUBROUTINE TEBA31(EK,AK,E,X,Y,TH,IDN,IW,IPSW,IERR)
C
C***********************************************************************
C
C   FINITE ELEMENT LIBRARY :  TEBA31                       SINTEF / NTH
C
C
C     TEBA31 GENERATES THE STIFFNESS MATRIX  EK  (AND THE MATRIX KAPPA)
C     FOR A HYBRID TRIANGULAR PLATE BENDING ELEMENT WITH  9  DEGREES
C     OF FREEDOM.  THE ELEMENT IS A GENERALIZATION OF THE ONE SUGGESTED
C     BY ALLMAN  -  IT MAY HAVE LINEARLY VARYING THICKNESS AND GENERAL
C     ANISOTROPIC MATERIAL PROPERTIES
C
C
C     PROGRAMMED BY :  K.BELL
C
C     DATE/VERSION  :  25.06.73 / 01
C
C***********************************************************************
C
      INTEGER           I,IDN,IERR,II,IPSW,IW,J,JJ,K,M,MM,N,NN
      INTEGER           IP(3),IWORK(9)
      DOUBLE PRECISION  EK(9,9),AK(9,9),E(3,3),X(3),Y(3),TH(3)
      DOUBLE PRECISION  C(3),S(3),SL(3),EI(3,3),B(3,3),G(9,12)
      DOUBLE PRECISION  T(12,9),Z1(7),Z2(7),Z3(7),W(7)
      DOUBLE PRECISION  F,A,RL11,RL12,RL21,RL22,THK,SS,CC,CS
C
      DATA IP(1),IP(2),IP(3)/2,3,1/
      DATA Z1(1),Z1(2),Z1(3),Z1(4),Z1(5),Z1(6),Z1(7) /0.33333333,
     *     0.05971587,2*0.47014206,0.79742699,2*0.10128651/
      DATA Z2(1),Z2(2),Z2(3),Z2(4),Z2(5),Z2(6),Z2(7) /0.33333333,
     *     0.47014206,0.05971587,0.47014206,0.10128651,0.79742699,
     *     0.10128651/
      DATA Z3(1),Z3(2),Z3(3),Z3(4),Z3(5),Z3(6),Z3(7) /0.33333333,
     *     2*0.47014206,0.05971587,2*0.10128651,0.79742699/
      DATA W(1),W(2),W(3),W(4),W(5),W(6),W(7) /0.225,3*0.13239415,
     *     3*0.12593918/
C
C
C     PRINT ?
C
      IF(IPSW-2) 30,10,10
   10 WRITE(IW,6000) IDN
      WRITE(IW,6010)
      WRITE(IW,6020) TH(1),TH(2),TH(3)
      WRITE(IW,6030)
      WRITE(IW,6040) X(1),X(2),X(3)
      WRITE(IW,6050) Y(1),Y(2),Y(3)
      WRITE(IW,6060)
      DO 20 I=1,3
   20 WRITE(IW,6070) E(I,1),E(I,2),E(I,3)
C
C     ELEMENT SIDE LENGTHS
C
   30 DO 40 I=1,3
      J = IP(I)
   40 SL(I) = DSQRT((X(J)-X(I))**2 + (Y(J)-Y(I))**2)
C
C     MAXIMUM SIDE LENGTH
C
      F = SL(1)
      DO 60 I=2,3
      IF(F-SL(I)) 50,60,60
   50 F = SL(I)
   60 CONTINUE
C
C     ELEMENT AREA  A
C
      A = X(1)*Y(2) + X(2)*Y(3) + X(3)*Y(1)
      A = A - X(1)*Y(3) - X(2)*Y(1) - X(3)*Y(2)
      A = 0.5*A
C
C     CORRECT NODE NUMBERING ?
C
      IF(A) 70,70,80
   70 IERR =-1
      GO TO 1000
C
C     DEGENERATED ELEMENT ?
C
   80 IF(50.*A - F*F) 90,90,100
   90 WRITE(IW,6081) IDN,A/(F*F)
C
C     SINE AND COSINE OF THE ANGLE GAMMA
C
  100 DO 110 I=1,3
      J = IP(I)
      S(I) = (X(I)-X(J))/SL(I)
  110 C(I) = (Y(J)-Y(I))/SL(I)
C
C     THE ELEMENTS OF THE INVERSE JACOBIAN MATRIX
C
      RL11 = 0.5*(Y(2)-Y(3))/A
      RL12 = 0.5*(Y(3)-Y(1))/A
      RL21 = 0.5*(X(3)-X(2))/A
      RL22 = 0.5*(X(1)-X(3))/A
C
C     THE INVERSE OF THE ELASTICITY MATRIX
C
      DO 120 I=1,3
      DO 120 J=1,3
  120 EI(I,J) = E(I,J)
      CALL DINV12(3,EI,B,IWORK,IW,IERR)
      IF (IERR) 130,200,130
  130 IERR =-3
      GO TO 1000
C
C     THE MATRIX  B
C
  200 DO 210 I=1,3
      DO 210 J=I,3
  210 B(I,J) = 0.
C
C     INTEGRATION LOOP
C
      DO 225 K=1,7
C
C     THICKNESS AT THE INTEGRATING POINT AND MULTIPLIER
C
      THK = TH(1)*Z1(K) + TH(2)*Z2(K) + TH(3)*Z3(K)
      IF (THK) 215,215,220
  215 IERR =-2
      GO TO 1000
C
  220 F   = 12.*A*W(K)/(THK*THK*THK)
C
      B(1,1) = B(1,1) + F*Z1(K)*Z1(K)
      B(1,2) = B(1,2) + F*Z1(K)*Z2(K)
      B(1,3) = B(1,3) + F*Z1(K)*Z3(K)
      B(2,2) = B(2,2) + F*Z2(K)*Z2(K)
      B(2,3) = B(2,3) + F*Z2(K)*Z3(K)
      B(3,3) = B(3,3) + F*Z3(K)*Z3(K)
  225 CONTINUE
C
      DO 240 I=1,3
      DO 240 J=I,3
  240 B(J,I) = B(I,J)
C
C     THE FLEXIBILITY MATRIX  (STORED TEMPORARILY IN EK)
C
      DO 325 II=1,3
      DO 325 JJ=II,3
      MM = 3*II-3
      NN = 3*JJ-3
      DO 310 I=1,3
      DO 310 J=1,3
      M = MM+I
      N = NN+J
  310 EK(M,N) = EI(II,JJ)*B(I,J)
  325 CONTINUE
C
      DO 350 I=1,9
      DO 350 J=I,9
  350 EK(J,I) = EK(I,J)
C
C     THE INVERSE FLEXIBILITY MATRIX  (TEMPORARILY STORED IN EK)
C
      CALL DINV12(9,EK,G,IWORK,IW,IERR)
      IF (IERR) 360,400,360
  360 IERR =-4
      GO TO 1000
C
C     ZERO THE TRANSFORMATION MATRICES  G  AND  T
C
  400 DO 410 I=1,9
      DO 410 J=1,12
      G(I,J) = 0.
  410 T(J,I) = 0.
C
C
C
C     THE TRANSFORMATION MATRIX  G
C
C
C     THE FIRST THREE COLUMNS
C
      DO 420 I=1,3
      J = IP(I)
      K = IP(J)
      G(I,I)   = S(K)*C(K) - S(I)*C(I)
      G(I+3,I) =-G(I,I)
      G(I+6,I) = C(I)*C(I) - S(I)*S(I) - C(K)*C(K) + S(K)*S(K)
  420 CONTINUE
C
C     THE 4TH, 5TH AND 6TH COLUMNS
C
      DO 440 I=1,3
      J = I+3
      SS = S(I)*S(I)*C(I)
      CC = C(I)*C(I)*S(I)
      CS = C(I)*C(I) - S(I)*S(I)
C
      G(1,J) = (C(I)+SS)*RL11 - CC*RL21
      G(2,J) = (C(I)+SS)*RL12 - CC*RL22
      G(3,J) =-(C(I)+SS)*(RL11+RL12) + CC*(RL21+RL22)
      G(4,J) = (S(I)+CC)*RL21 - SS*RL11
      G(5,J) = (S(I)+CC)*RL22 - SS*RL12
      G(6,J) =-(S(I)+CC)*(RL21+RL22) + SS*(RL11+RL12)
      G(7,J) = S(I)*RL11 + C(I)*RL21 - CS*(S(I)*RL11-C(I)*RL21)
      G(8,J) = C(I)*RL22 + S(I)*RL12 - CS*(S(I)*RL12-C(I)*RL22)
      G(9,J) =-S(I)*(RL11+RL12) - C(I)*(RL21+RL22)
      G(9,J) = G(9,J) + CS*(S(I)*(RL11+RL12) - C(I)*(RL21+RL22))
  440 CONTINUE
C
C     THE LAST SIX COLUMNS
C
      DO 450 I=1,3
      J = IP(I)
      N = 2*I+5
      G(I,N)     = C(I)*C(I)
      G(I+3,N)   = S(I)*S(I)
      G(I+6,N)   = 2.*S(I)*C(I)
      G(J,N+1)   = C(I)*C(I)
      G(J+3,N+1) = S(I)*S(I)
      G(J+6,N+1) = 2.*S(I)*C(I)
  450 CONTINUE
C
C
C     THE TRANSFORMATION MATRIX  T
C
C     A)  CORRESPONDING TO W, W,X AND W,Y
C
      DO 500 I=1,3
      J = IP(I)
      SS = S(I)*SL(I)
      CC = C(I)*SL(I)
      T(I,3*I-2)     = 1.
      T(I+3,3*I-2)   = SL(I)/2.
      T(I+3,3*I-1)   =-SL(I)*SS/12.
      T(I+3,3*I)     = SL(I)*CC/12.
      T(I+3,3*J-2)   = SL(I)/2.
      T(I+3,3*J-1)   = SL(I)*SS/12.
      T(I+3,3*J)     =-SL(I)*CC/12.
      T(2*I+5,3*I-1) =-CC/3.
      T(2*I+5,3*I)   =-SS/3.
      T(2*I+5,3*J-1) =-CC/6.
      T(2*I+5,3*J)   =-SS/6.
      T(2*I+6,3*I-1) =-CC/6.
      T(2*I+6,3*I)   =-SS/6.
      T(2*I+6,3*J-1) =-CC/3.
      T(2*I+6,3*J)   =-SS/3.
  500 CONTINUE
C
C     B)  CORRESPONDING TO W, THETAX AND THETAY
C
      DO 550 I=1,3
      J = 3*I-1
      K = J+1
      DO 525 M =1,12
      F      = T(M,J)
      T(M,J) = T(M,K)
  525 T(M,K) = -F
  550 CONTINUE
C
C     COMPUTE THE MATRIX  A = G*T   (STORE TEMPORARILY IN AK)
C
      DO 600 I=1,9
      DO 600 J=1,9
      AK(I,J) = 0.
      DO 600 K=1,12
  600 AK(I,J) = AK(I,J) + G(I,K)*T(K,J)
C
C     COMPUTE THE MATRIX  KAPPA   (STORE TEMPORARILY IN  G )
C
      DO 700 I=1,9
      DO 700 J=1,9
      G(I,J) = 0.
      DO 700 K=1,9
  700 G(I,J) = G(I,J) + EK(I,K)*AK(K,J)
C
C     COMPUTE THE FINAL STIFFNESS MATRIX
C
      DO 800 I=1,9
      DO 800 J=I,9
      EK(I,J) = 0.
      DO 800 K=1,9
  800 EK(I,J) = EK(I,J) + AK(K,I)*G(K,J)
C
      DO 825 I=1,9
      DO 825 J=I,9
  825 EK(J,I) = EK(I,J)
C
C     TRANSFER THE MATRIX  KAPPA   (FROM  G  TO  AK )
C
      DO 900 I=1,9
      DO 900 J=1,9
  900 AK(I,J) = G(I,J)
C
      IERR = 0
      GO TO 2000
C
C     ERROR RETURN
C
 1000 WRITE(IW,6080)
      WRITE(IW,6090) IERR
      WRITE(IW,6100) IDN
      IF(IERR+2) 2000,1100,2000
 1100 WRITE(IW,6102) THK
C
C     PRINT ?
C
 2000 IF(IPSW-2) 3000,2020,2010
 2010 WRITE(IW,6110)
      CALL MPRT30(EK,9,9,IW)
 2020 WRITE(IW,6120) IDN
C
 3000 RETURN
C
 6000 FORMAT(///28H ENTERING TEBA31 FOR ELEMENT,I8,2H :)
 6010 FORMAT(/25H     CORNER THICKNESSES :)
 6020 FORMAT(9H     H1 =,1PE11.4,8H    H2 =,1PE11.4,8H    H3 =,1PE11.4)
 6030 FORMAT(/25H     CORNER COORDINATES :)
 6040 FORMAT(9H     X1 =,1PE11.4,8H    X2 =,1PE11.4,8H    X3 =,1PE11.4)
 6050 FORMAT(9H     Y1 =,1PE11.4,8H    Y2 =,1PE11.4,8H    Y3 =,1PE11.4)
 6060 FORMAT(/25H     ELASTICITY MATRIX  :)
 6070 FORMAT(1PE18.7,1P,2E15.7)
 6080 FORMAT(///33H *** ERROR RETURN FROM TEBA31 ***)
 6081 FORMAT(33H  ** WARNING FROM TEBA31: ELEMENT,I7,
     +       22H HAS BAD STRETCH VALUE,1PE11.4)
 6090 FORMAT(/17H     ERROR FLAG =,I3)
 6100 FORMAT(/12H     ELEMENT,I8,11H   IN ERROR)
 6102 FORMAT( 35H     INVALID PLATE THICKNESS, THK =,1PE12.5)
 6110 FORMAT(/31H     ELEMENT STIFFNESS MATRIX :)
 6120 FORMAT(///27H LEAVING TEBA31 FOR ELEMENT,I8)
C
      END
      SUBROUTINE TEBA32(SM,AK,Z)
C
C***********************************************************************
C
C   FINITE ELEMENT LIBRARY :  TEBA32                       SINTEF / NTH
C
C
C     TEBA32 GENERATES A STRESS MATRIX SM FOR A HYBRID TRIANGULAR
C     PLATE BENDING ELEMENT WITH 9 DEGREES OF FREEDOM
C
C
C     PROGRAMMED BY :  K.BELL
C
C     DATE/VERSION  :  07.08.74 / 01
C
C***********************************************************************
C
      INTEGER           I,J,K,L,M
      DOUBLE PRECISION  SM(3,9),AK(9,9),Z(3)
C
      DO 100 I=1,3
      L = 3*I-3
      DO 50 J=1,9
      SM(I,J) = 0.
      DO 25 K=1,3
      M = L+K
   25 SM(I,J) = SM(I,J) + Z(K)*AK(M,J)
   50 CONTINUE
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE TEBA33(F,X,Y,P)
C
C***********************************************************************
C
C   FINITE ELEMENT LIBRARY :  TEBA33                       SINTEF / NTH
C
C
C     TEBA33 DETERMINES A LUMPED LOAD VECTOR  F  DUE TO A LINEARLY
C     VARYING DISTRIBUTED LOAD ON THE TRIANGULAR PLATE BENDING ELEMENT
C     TEBA
C
C
C     PROGRAMMED BY :  K.BELL
C
C     DATE/VERSION  :  16.10.74 / 01
C
C***********************************************************************
C
      INTEGER           I
      DOUBLE PRECISION  F(9),X(3),Y(3),P(3),A
C
      DO 10 I=1,9
      F(I) = 0.
   10 CONTINUE
C
C     ELEMENT AREA
C
      A = X(1)*Y(2) + X(2)*Y(3) + X(3)*Y(1)
      A = A - X(1)*Y(3) - X(2)*Y(1) - X(3)*Y(2)
      A = 0.5*A
C
      IF (A) 50,50,20
C
C     THE NON-ZERO ELEMENTS OF  F
C
   20 A = A/12.
      F(1) = A*(2.*P(1) + P(2) + P(3))
      F(4) = A*(P(1) + 2.*P(2) + P(3))
      F(7) = A*(P(1) + P(2) + 2.*P(3))
C
   50 RETURN
      END
      SUBROUTINE TEBA35(EM,X,Y,Z,THK,RHO)
C
C***********************************************************************
C
C     NYTEBA35 DETERMINES A LUMPED MASS MATRIX  EM  FOR THE TRIANGULAR
C     PLATE BENDING ELEMENT  TEBA
C
C
C     PROGRAMMED BY :  KETIL AAMNES
C
C     DATE/VERSION  :  08.11.85 / 1.0
C
C***********************************************************************
C
      DOUBLE PRECISION  EM(9),X(3),Y(3),Z(3),THK,RHO
      DOUBLE PRECISION  AUX1,COSG,SING,SL21,SL31
      DOUBLE PRECISION  X21,X31,Y21,Y31,Z21,Z31,XL(3),YL(3)
      DOUBLE PRECISION  AREA,XTP,YTP,H,B,B1,B2,S1,S2
      DOUBLE PRECISION  IXX,IYY1,IYY2,IYY,IZZ1,IZZ2,IZZ
      DOUBLE PRECISION  IX,IY,IZ,M,GAMMA,TR(3,3)
C
C ---------------------------- Local coordinates -------------------
C
      X21 = X(2)-X(1)
      Y21 = Y(2)-Y(1)
      Z21 = Z(2)-Z(1)
      X31 = X(3)-X(1)
      Y31 = Y(3)-Y(1)
      Z31 = Z(3)-Z(1)
      SL21 = DSQRT(X21*X21 + Y21*Y21 + Z21*Z21)
      SL31 = DSQRT(X31*X31 + Y31*Y31 + Z31*Z31)
      COSG = (X31*X21 + Y31*Y21 + Z31*Z21)/(SL31*SL21)
      AUX1 = 1.-COSG*COSG
      IF (AUX1) 40,40,50
   40 SING = 0.
      GO TO 60
   50 SING = DSQRT(AUX1)
C
   60 XL(1) = 0.
      XL(2) = SL21
      XL(3) = SL31*COSG
      YL(1) = 0.
      YL(2) = 0.
      YL(3) = SL31*SING
C
C ---------------------------- Element area ------------------------
C
      AREA = (XL(2)*YL(3))/2.
C
      IF (AREA) 80,80,70
C
C ---------------------------- Koordinates for TP ------------------
C
   70 XTP = ((XL(3)+((XL(2)-XL(3))/2.)))*(2./3.)
      YTP = (YL(3))/3.
C
C ------------------------------------------------------------------
C
      H = YL(3)
      B = XL(2)
C
      B1 = XL(3)
      B2 = XL(2)-XL(3)
C
      S1 = DSQRT( (((2./3.)*B1)**2) + (((1./3.)*H)**2) )
      S2 = DSQRT( (((2./3.)*B2)**2) + (((1./3.)*H)**2) )
C
C ---------------------------- Moment of inertia about the X-aksis -
C
      IXX = (RHO*THK*B*H) * (((H**2)/36.) + ((THK**2)/24.))
C
C ---------------------------- Moment of inertia about the Y-aksis -
C
      IYY1 = (RHO*THK*B1*H) * (((B1**2)/36.) + ((THK**2)/24.))
      IYY1 = IYY1 + (((XTP-((2./3.)*B1))**2)*(((B1*H)/2.)*RHO*THK))
C
      IYY2 = (RHO*THK*B2*H) * (((B2**2)/36.) + ((THK**2)/24.))
      IYY2 = IYY2 + (((B1-XTP+((1./3.)*B2))**2)*(((B2*H)/2.)*RHO*THK))
C
      IYY = IYY1 + IYY2
C
C ---------------------------- Moment of inertia about the Z-aksis -
C
      IZZ1 = (RHO*THK*H*B1) * (((B1**2)/4.)+((H**2)/12.))
      IZZ1 = IZZ1 - ((S1**2)*RHO*THK*(B1*H)/2.)
      IZZ1 = IZZ1 + (((XTP-((2./3.)*B1))**2)*((B1*H)/2.)*RHO*THK)
C
      IZZ2 = (RHO*THK*H*B2) * (((B2**2)/4.)+((H**2)/12.))
      IZZ2 = IZZ2 - ((S2**2)*RHO*THK*(B2*H)/2.)
      IZZ2 = IZZ2 + (((B1+((1./3.)*B2)-XTP)**2)*((B2*H)/2.)*RHO*THK)
C
      IZZ = IZZ1 + IZZ2
C
C ---------------------------- The 3 X 3 transformation matrix -----
C ---------------------------- V1 ----------------------------------
      TR(1,1) = X(2)-X(1)
      TR(2,1) = Y(2)-Y(1)
      TR(3,1) = Z(2)-Z(1)
      S1 = DSQRT(TR(1,1)**2+TR(2,1)**2+TR(3,1)**2)
      TR(1,1) = TR(1,1)/S1
      TR(2,1) = TR(2,1)/S1
      TR(3,1) = TR(3,1)/S1
C ---------------------------- V3 ----------------------------------
      TR(1,2) = X(3)-X(1)
      TR(2,2) = Y(3)-Y(1)
      TR(3,2) = Z(3)-Z(1)
      TR(1,3) = (TR(2,1)*TR(3,2)-TR(2,2)*TR(3,1))
      TR(2,3) = (TR(1,2)*TR(3,1)-TR(1,1)*TR(3,2))
      TR(3,3) = (TR(1,1)*TR(2,2)-TR(1,2)*TR(2,1))
C
      S1 = DSQRT(TR(1,3)**2+TR(2,3)**2+TR(3,3)**2)
      TR(1,3) = TR(1,3)/S1
      TR(2,3) = TR(2,3)/S1
      TR(3,3) = TR(3,3)/S1
C ---------------------------- V2 ----------------------------------
      TR(1,2) = TR(2,3)*TR(3,1)-TR(2,1)*TR(3,3)
      TR(2,2) = TR(1,1)*TR(3,3)-TR(1,3)*TR(3,1)
      TR(3,2) = TR(1,3)*TR(2,1)-TR(1,1)*TR(2,3)
C
C ---------------------------- Transform to global
C                              moments of inertia ------------------
      IX = ABS(TR(1,1)*IXX+TR(1,2)*IYY+TR(1,3)*IZZ)
      IY = ABS(TR(2,1)*IXX+TR(2,2)*IYY+TR(2,3)*IZZ)
      IZ = ABS(TR(3,1)*IXX+TR(3,2)*IYY+TR(3,3)*IZZ)
C
C ---------------------------- The elements of EM ------------------
C
      M = (RHO*AREA*THK)/3.
      GAMMA = 1./20.
C
      EM(1) = M
      EM(4) = EM(1)
      EM(7) = EM(1)
C
      EM(2) = GAMMA*IX
      EM(5) = EM(2)
      EM(8) = EM(2)
C
      EM(3) = GAMMA*IY
      EM(6) = EM(3)
      EM(9) = EM(3)
C ---------------------------- Return ------------------------------
   80 RETURN
      END
      SUBROUTINE TEBA36(EKG,X,Y,EN)
C
C***********************************************************************
C
C   FINITE ELEMENT LIBRARY :  TEBA36                       SINTEF / NTH
C
C
C     TEBA36 GENERATES A GEOMETRICAL STIFFNESS MATRIX  EKG  FOR THE
C     TRIANGULAR PLATE BENDING ELEMENT  TEBA
C     EKG  IS BASED ON CONSTANT IN-PLANE FORCES AND A LINEAR DISPLACE-
C     MENT FIELD
C
C
C     PROGRAMMED BY :  K.BELL
C
C     DATE/VERSION  :  17.10.74 / 01
C
C***********************************************************************
C
      INTEGER           I,J
      DOUBLE PRECISION  EKG(9,9),X(3),Y(3),EN(3)
      DOUBLE PRECISION  A,X21,X32,X13,Y12,Y23,Y31,SX,SY,SXY
C
      DO 10 I=1,9
      DO 10 J=I,9
      EKG(I,J) = 0.0
   10 CONTINUE
C
C     SOME GEOMETRICAL CONSTANTS
C
      A = X(1)*Y(2) + X(2)*Y(3) + X(3)*Y(1)
      A = A - X(1)*Y(3) - X(2)*Y(1) - X(3)*Y(2)
      A = 0.5*A
C
      IF (A) 30,30,20
C
   20 X32 = X(3)-X(2)
      X13 = X(1)-X(3)
      X21 = X(2)-X(1)
      Y23 = Y(2)-Y(3)
      Y31 = Y(3)-Y(1)
      Y12 = Y(1)-Y(2)
C
C     THE IN-PLANE FORCES MULTIPLIED BY A FACTOR=1/4A
C
      SX  = EN(1)*0.25/A
      SY  = EN(2)*0.25/A
      SXY = EN(3)*0.25/A
C
C     THE GEOMETRICAL STIFFNESS MATRIX
C
      EKG(1,1) = SX*Y23*Y23 + SY*X32*X32 + 2.*SXY*X32*Y23
      EKG(1,4) = SX*Y23*Y31 + SY*X32*X13 + SXY*(X32*Y31+Y23*X13)
      EKG(1,7) = SX*Y23*Y12 + SY*X32*X21 + SXY*(X32*Y12+Y23*X21)
      EKG(4,4) = SX*Y31*Y31 + SY*X13*X13 + 2.*SXY*X13*Y31
      EKG(4,7) = SX*Y31*Y12 + SY*X13*X21 + SXY*(X13*Y12+Y31*X21)
      EKG(7,7) = SX*Y12*Y12 + SY*X21*X21 + 2.*SXY*X21*Y12
C
   30 DO 40 I=1,9
      DO 40 J=I,9
      EKG(J,I) = EKG(I,J)
   40 CONTINUE
C
      RETURN
      END
