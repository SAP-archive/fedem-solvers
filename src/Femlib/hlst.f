C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
C     This file is part of FEDEM - https://openfedem.org
C
      SUBROUTINE HLST31(EK,AK,E,X,Y,THK,IOP,IDN,IW,IPSW,IERR)
C
C***********************************************************************
C
C   FINITE ELEMENT LIBRARY :  HLST31                       SINTEF / NTH
C
C
C     HLST31 GENERATES THE STIFFNESS MATRIX  EK  (AND THE MATRIX KAPPA)
C     FOR A HYBRID TRIANGULAR MEMBRANE ELEMENT WITH 9 DEGREES OF
C     FREEDOM
C     THE ELEMENT HAS A LINEAR STRESS VARIATION AND TWO VERSIONS (CON-
C     TROLLED BY IOP) OF COMPATIBLE (CUBIC) EDGE DISPLACEMENTS
C     THE THICKNESS IS CONSTANT
C     GENERAL ANISOTROPIC MATERIAL PROPERTIES MAY BE SPECIFIED
C
C
C     PROGRAMMED BY :  K.BELL
C
C     DATE/VERSION  :  11.11.74 / 01
C
C***********************************************************************
C
      INTEGER           I,IDN,IERR,IOP,IPSW,IW,J,K,IP(3),IWORK(7)
      DOUBLE PRECISION  EK(9,9),AK(7,9),E(3,3),X(3),Y(3),C(3),S(3)
      DOUBLE PRECISION  SL(3),EI(3,3),A(7,9),AUX(9),XL(3),YL(3)
      DOUBLE PRECISION  A1,A4,A5,AREA,B1,B2,B3,B4,B5,C2,CX,CY,F
      DOUBLE PRECISION  P02,P11,P20,S2,SC,SX,SY,THK,X0,Y0
C
      DATA IP(1),IP(2),IP(3)/2,3,1/
C
C
C     PRINT ?
C
      IF (IPSW-2) 30,10,10
C
   10 WRITE(IW,6000) IDN
      WRITE(IW,6010) THK
      WRITE(IW,6020)
      WRITE(IW,6030) X(1),X(2),X(3)
      WRITE(IW,6040) Y(1),Y(2),Y(3)
      WRITE(IW,6050)
      DO 20 I=1,3
      WRITE(IW,6060) E(I,1),E(I,2),E(I,3)
   20 CONTINUE
C
C     ELEMENT SIDE LENGTHS
C
   30 DO 50 I=1,3
      J = IP(I)
      SL(I) = DSQRT((X(J)-X(I))**2 + (Y(J)-Y(I))**2)
   50 CONTINUE
C
C     MAXIMUM SIDE LENGTH
C
      F = SL(1)
      DO 70 I=2,3
      IF (F-SL(I)) 60,70,70
   60 F = SL(I)
   70 CONTINUE
C
C     ELEMENT AREA
C
      AREA = X(1)*Y(2) + X(2)*Y(3) + X(3)*Y(1)
      AREA = AREA - X(1)*Y(3) - X(2)*Y(1) - X(3)*Y(2)
      AREA = 0.5*AREA
C
C     CORRECT NODE NUMBERING ?
C
      IF (AREA) 80,100,90
   80 IERR = -1
      GO TO 1000
C
C     DEGENERATED ELEMENT ?
C
   90 IF (50.*AREA - F*F) 100,100,110
  100 WRITE(IW,6071) IDN,AREA/(F*F)
C
C     SINE AND COSINE OF THE ANGLE GAMMA
C
  110 DO 120 I=1,3
      J    = IP(I)
      S(I) = (X(I)-X(J))/SL(I)
      C(I) = (Y(J)-Y(I))/SL(I)
  120 CONTINUE
C
C     LOCAL COORDINATES AND INTEGRATION FORMULAS
C
      X0 = (X(1)+X(2)+X(3))/3.
      Y0 = (Y(1)+Y(2)+Y(3))/3.
C
      DO 150 I=1,3
      XL(I) = X(I)-X0
      YL(I) = Y(I)-Y0
  150 CONTINUE
C
      F   = AREA/(12.*THK)
      P20 = F*(XL(1)*XL(1) + XL(2)*XL(2) + XL(3)*XL(3))
      P11 = F*(XL(1)*YL(1) + XL(2)*YL(2) + XL(3)*YL(3))
      P02 = F*(YL(1)*YL(1) + YL(2)*YL(2) + YL(3)*YL(3))
C
C     THE INVERSE  EI  OF THE ELASTICITY MATRIX  E
C
      DO 170 I=1,3
      DO 170 J=1,3
      EI(I,J) = E(I,J)
  170 CONTINUE
      CALL DINV12(3,EI,AK,IWORK,IW,IERR)
      IF (IERR) 180,200,180
  180 IERR = -3
      GO TO 1000
C
C     THE FLEXIBILITY MATRIX  F  (TEMPORARILY STORED IN  AK )
C
  200 DO 210 I=1,7
      DO 210 J=I,7
      AK(I,J) = 0.
  210 CONTINUE
C
      F       = AREA/THK
C
      AK(1,1) = F*EI(1,1)
      AK(1,2) = F*EI(1,2)
      AK(2,2) = F*EI(2,2)
      AK(1,3) = F*EI(1,3)
      AK(2,3) = F*EI(2,3)
      AK(3,3) = F*EI(3,3)
      AK(4,4) = EI(1,1)*P20 - 2.*EI(1,3)*P11 + EI(3,3)*P02
      AK(4,5) = EI(1,2)*P20 - EI(2,3)*P11
      AK(5,5) = EI(2,2)*P20
      AK(4,6) = EI(1,1)*P11 - EI(1,3)*P02
      AK(5,6) = EI(1,2)*P11
      AK(6,6) = EI(1,1)*P02
      AK(4,7) =-EI(1,3)*P20 + (EI(1,2)+EI(3,3))*P11 - EI(2,3)*P02
      AK(5,7) = EI(2,2)*P11 - EI(2,3)*P20
      AK(6,7) = EI(1,2)*P02 - EI(1,3)*P11
      AK(7,7) = EI(2,2)*P02 - 2.*EI(2,3)*P11 + EI(3,3)*P20
C
      DO 230 I=1,7
      DO 230 J=I,7
      AK(J,I) = AK(I,J)
  230 CONTINUE
C
C     THE INVERSE FLEXIBILITY MATRIX  (TEMPORARILY STORED IN  AK )
C
      CALL DINV12(7,AK,A,IWORK,IW,IERR)
      IF (IERR) 250,300,250
  250 IERR = -4
      GO TO 1000
C
C     THE TRANSFORMATION MATRIX  A
C
  300 DO 310 J=1,9
      DO 310 I=1,7
      A(I,J) = 0.
  310 CONTINUE
C
      DO 500 K=1,3
C
      C2 = C(K)*C(K)
      S2 = S(K)*S(K)
      SC = S(K)*C(K)
      CX = C(K)*XL(K)
      CY = C(K)*YL(K)
      SX = S(K)*XL(K)
      SY = S(K)*YL(K)
      F  = SL(K)*SL(K)
C
      A1 = SL(K)/2.
      A4 =-C(K)*F/12.
      A5 =-S(K)*F/12.
C
      DO 400 I=1,2
C
      IF (I-1) 320,320,350
C
  320 J  = 3*K-2
      B4 = -C(K)*F*SL(K)/30.
      B5 = -S(K)*F*SL(K)/30.
      IF (IOP) 330,330,340
  330 B1 = F*(9.*C2 + 10.*S2)/60.
      B2 = -F*SC/60.
      B3 = F*(10.*C2 + 9.*S2)/60.
      GO TO 380
  340 B1 = F/6.
      B2 = 0.
      B3 = B1
      GO TO 380
C
  350 J  = 3*IP(K)-2
      A4 = -A4
      A5 = -A5
      B4 = C(K)*F*SL(K)/20.
      B5 = S(K)*F*SL(K)/20.
      IF (IOP) 360,360,370
  360 B1 = F*(21.*C2 + 20.*S2)/60.
      B2 = F*SC/60.
      B3 = F*(20.*C2 + 21.*S2)/60.
      GO TO 380
  370 B1 = F/3.
      B2 = 0.
      B3 = B1
C
  380 A(1, J ) = A(1, J ) + A1*C(K)
      A(3, J ) = A(3, J ) + A1*S(K)
      A(4, J ) = A(4, J ) + A1*(CX-SY) - 2.*B1*SC - B2*C2
      A(5, J ) = A(5, J ) - B2*S2
      A(6, J ) = A(6, J ) + A1*CY + B1*C2
      A(7, J ) = A(7, J ) - A1*SX + B1*S2 + 2.*B2*SC
      A(2,J+1) = A(2,J+1) + A1*S(K)
      A(3,J+1) = A(3,J+1) + A1*C(K)
      A(4,J+1) = A(4,J+1) - A1*CY - 2.*B2*SC - B3*C2
      A(5,J+1) = A(5,J+1) + A1*SX - B3*S2
      A(6,J+1) = A(6,J+1) + B2*C2
      A(7,J+1) = A(7,J+1) + A1*(SY-CX) + B2*S2 + 2.*B3*SC
      A(1,J+2) = A(1,J+2) + A4*C(K)
      A(2,J+2) = A(2,J+2) + A5*S(K)
      A(3,J+2) = A(3,J+2) + A4*S(K) + A5*C(K)
      A(4,J+2) = A(4,J+2) + A4*(CX-SY) - A5*CY - 2.*B4*SC - B5*C2
      A(5,J+2) = A(5,J+2) + A5*SX - B5*S2
      A(6,J+2) = A(6,J+2) + A4*CY + B4*C2
      A(7,J+2) = A(7,J+2) - A4*SX + A5*(SY-CX) + B4*S2 + 2.*B5*SC
  400 CONTINUE
  500 CONTINUE
C
C     DETERMINE THE MATRIX  KAPPA  (STORED IN  AK )
C
      DO 600 I=1,7
      DO 550 J=1,9
      AUX(J) = AK(I,1)*A(1,J)
      DO 520 K=2,7
      AUX(J) = AUX(J) + AK(I,K)*A(K,J)
  520 CONTINUE
  550 CONTINUE
      DO 570 J=1,9
      AK(I,J) = AUX(J)
  570 CONTINUE
  600 CONTINUE
C
C     DETERMINE THE STIFFNESS MATRIX  EK
C
      DO 700 I=1,9
      DO 650 J=I,9
      EK(I,J) = A(1,I)*AK(1,J)
      DO 620 K=2,7
      EK(I,J) = EK(I,J) + A(K,I)*AK(K,J)
  620 CONTINUE
  650 CONTINUE
  700 CONTINUE
C
      DO 800 I=1,9
      DO 800 J=I,9
      EK(J,I) = EK(I,J)
  800 CONTINUE
C
      IERR = 0
      GO TO 2000
C
C     ERROR RETURN
C
 1000 WRITE(IW,6070)
      WRITE(IW,6080) IERR
      WRITE(IW,6090) IDN
C
C     PRINT ?
C
 2000 IF (IPSW-2) 3000,2020,2010
 2010 WRITE(IW,6100)
      CALL MPRT30(EK,9,9,IW)
 2020 WRITE(IW,6110) IDN
 3000 RETURN
C
 6000 FORMAT(///29H ENTERING HLST31 FOR ELEMENT ,I7,2H :)
 6010 FORMAT(/24H     ELEMENT THICKNESS :,1PE11.4)
 6020 FORMAT(/25H     CORNER COORDINATES :)
 6030 FORMAT(9H     X1 =,1PE11.4,8H    X2 =,1PE11.4,8H    X3 =,1PE11.4)
 6040 FORMAT(9H     Y1 =,1PE11.4,8H    Y2 =,1PE11.4,8H    Y3 =,1PE11.4)
 6050 FORMAT(/25H     ELASTICITY MATRIX  :)
 6060 FORMAT(1PE18.7,1P,2E15.7)
 6070 FORMAT(///33H *** ERROR RETURN FROM HLST31 ***)
 6071 FORMAT(33H  ** WARNING FROM HLST31: ELEMENT,I7,
     +       22H HAS BAD STRETCH VALUE,1PE11.4)
 6080 FORMAT(/17H     ERROR FLAG =,I3)
 6090 FORMAT(/13H     ELEMENT ,I7,11H   IN ERROR)
 6100 FORMAT(/31H     ELEMENT STIFFNESS MATRIX :)
 6110 FORMAT(///28H LEAVING HLST31 FOR ELEMENT ,I7)
C
      END
      SUBROUTINE HLST32(SM,AK,X,Y,Z)
C
C***********************************************************************
C
C   FINITE ELEMENT LIBRARY :  HLST32                       SINTEF / NTH
C
C
C     HLST32 GENERATES A STRESS MATRIX  SM  FOR THE  HLST  ELEMENT
C
C
C     PROGRAMMED BY :  K.BELL
C
C     DATE/VERSION  :  19.12.74 / 01
C
C***********************************************************************
C
      INTEGER           I,J,K
      DOUBLE PRECISION  SM(3,9),AK(7,9),X(3),Y(3),Z(3)
      DOUBLE PRECISION  RN(3,7),XL(3),YL(3),X0,Y0
C
C     LOCAL COORDINATES
C
      X0 = (X(1) + X(2) + X(3))/3.
      Y0 = (Y(1) + Y(2) + Y(3))/3.
      DO 5 I=1,3
      XL(I) = X(I) - X0
      YL(I) = Y(I) - Y0
    5 CONTINUE
C
C     THE MATRIX N  (STORED IN RN)
C
      DO 10 J=1,7
      DO 10 I=1,3
      RN(I,J) = 0.
   10 CONTINUE
C
      DO 20 I=1,3
      RN(1,4) = RN(1,4) + XL(I)*Z(I)
      RN(1,6) = RN(1,6) + YL(I)*Z(I)
      RN(2,5) = RN(2,5) + XL(I)*Z(I)
      RN(2,7) = RN(2,7) + YL(I)*Z(I)
      RN(3,4) = RN(3,4) - YL(I)*Z(I)
      RN(3,7) = RN(3,7) - XL(I)*Z(I)
   20 CONTINUE
C
C     THE STRESS MATRIX
C
      DO 50 I=1,3
      DO 50 J=1,9
      SM(I,J) = AK(I,J)
   50 CONTINUE
C
      DO 100 I=1,3
      DO  90 J=1,9
      DO  80 K=4,7
      SM(I,J) = SM(I,J) + RN(I,K)*AK(K,J)
   80 CONTINUE
   90 CONTINUE
  100 CONTINUE
C
      RETURN
      END
