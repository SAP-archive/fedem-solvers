C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
C     This file is part of FEDEM - https://openfedem.org
C
      SUBROUTINE FQS31 (EK,E,X,Y,Z,TH,ALPHA,BETA,
     +                  IDN,IW,IPSW,IERR,EMOD,NY)
C
C***********************************************************************
C
C     FQS31 GENERATES THE (24*24) STIFFNESS MATRIX  EK  FOR A FLAT
C     QUADRILATERAL THIN SHELL ELEMENT, COMPOSED OF A MEMBRANE ELEMENT
C     WITH ROTATIONAL DEGREES OF FREEDOM (B.LUND) AND THE QBESH BENDING
C     ELEMENT (WANG WIUXI).
C     THE ELEMENT MAY HAVE LINEARLY VARYING THICKNESS AND GENERAL
C     ANISOTROPIC MATERIAL PROPERTIES.
C
C
C     PROGRAMMED BY :  TERJE ROLVAG
C     DATE/VERSION  :  01.09.87/1
C
C***********************************************************************
C
      IMPLICIT NONE
C
      INTEGER           IDN,IW,IPSW,IERR
      INTEGER           I,II,J,JJ,K,L,L1,LL,M,MM,N,NN
      INTEGER           ISHEAR,IPERM(12),JPERM(12)
      DOUBLE PRECISION  EK(24,24),E(3,3),DM(3,3),NY
      DOUBLE PRECISION  X(4),Y(4),Z(4),TH(4),AELPAT(12),AALF(4)
      DOUBLE PRECISION  PEK(12,12),C(3,3),XL(4),YL(4),ALPHA,BETA,EMOD
      DOUBLE PRECISION  AUX3,AUX4,COSG3,SING3,COSG4,SING4,SL21,SL31
      DOUBLE PRECISION  SL41,THK,X21,X31,X41,Y21,Y31,Y41,Z21,Z31,Z41
      DOUBLE PRECISION  ELAS(2,2),ELAB(3,3)
C
      DATA  IPERM /3,4,5,9,10,11,15,16,17,21,22,23/
      DATA  JPERM /1,2,6,7, 8,12,13,14,18,19,20,24/
C
C ---------------------------- Print ? -----------------------------
C
      IF (IPSW-2) 30,10,10
   10 WRITE(IW,6000) IDN
      WRITE(IW,6010)
      WRITE(IW,6020) TH(1),TH(2),TH(3),TH(4)
      WRITE(IW,6030)
      WRITE(IW,6040) X(1),X(2),X(3),X(4)
      WRITE(IW,6050) Y(1),Y(2),Y(3),Y(4)
      WRITE(IW,6055) Z(1),Z(2),Z(3),Z(4)
      WRITE(IW,6060)
      DO 20 I=1,3
         WRITE(IW,6070) E(I,1),E(I,2),E(I,3)
   20 CONTINUE
C
C ----------- Average element thickness THK og Elastisitetsmodul -------
C
   30 THK = (TH(1)+TH(2)+TH(3)+TH(4)) * 0.25D0
C
C ---------------------------- Local coordinates -------------------
C
      X21 = X(2)-X(1)
      Y21 = Y(2)-Y(1)
      Z21 = Z(2)-Z(1)
      X31 = X(3)-X(1)
      Y31 = Y(3)-Y(1)
      Z31 = Z(3)-Z(1)
      X41 = X(4)-X(1)
      Y41 = Y(4)-Y(1)
      Z41 = Z(4)-Z(1)
C
      SL21 = DSQRT(X21*X21 + Y21*Y21 + Z21*Z21)
      SL31 = DSQRT(X31*X31 + Y31*Y31 + Z31*Z31)
      SL41 = DSQRT(X41*X41 + Y41*Y41 + Z41*Z41)
C
C     BEREGNER COSINUS OG SINUS TIL VINKELEN MELLOM
C  -- DEN LOKAL X-AKSEN OG HENHOLDSVIS KN.PUNKTENE 3 (G3) OG 4 (G4) --
C
      COSG3 = (X31*X21 + Y31*Y21 + Z31*Z21)/(SL31*SL21)
      COSG4 = (X41*X21 + Y41*Y21 + Z41*Z21)/(SL41*SL21)
C
      AUX3 = 1.0D0-COSG3*COSG3
      IF (AUX3) 40,40,50
   40 SING3 = 0.0D0
      GO TO 60
   50 SING3 = DSQRT(AUX3)
C
   60 AUX4 = 1.0D0-COSG4*COSG4
      IF (AUX4) 70,70,80
   70 SING4 = 0.0D0
      GO TO 90
   80 SING4 = DSQRT(AUX4)
C
   90 XL(1) = 0.0D0
      XL(2) = SL21
      XL(3) = SL31*COSG3
      XL(4) = SL41*COSG4
      YL(1) = 0.0D0
      YL(2) = 0.0D0
      YL(3) = SL31*SING3
      YL(4) = SL41*SING4
C
C ----------------------------- Initialisering ----------------------
C
      DO 110 I=1,24
         DO 100 J=1,24
            EK(I,J) = 0.0D0
  100    CONTINUE
  110 CONTINUE
C
C ----------------    Calculating 3*3 matrix of direction cosines -
C
      CALL DIRCOS(X,Y,Z,C)
C
C ---------------------------- Membrane stiffness ------------------
C
      DO 130 I=1,3
         DO 120 J=1,3
            DM(I,J) = THK * E(I,J)
  120    CONTINUE
  130 CONTINUE
C
      DO 150 I=1,12
         DO 140 J=1,12
            PEK(I,J) = 0.0D0
  140    CONTINUE
  150 CONTINUE
C
      CALL QMRF31(PEK,XL,YL,DM,EMOD,NY,THK,ALPHA,BETA,IPSW-2,IW,IERR)
      IF (IERR) 1000,160,160
C
C ---------------------------- Add membrane stiffness --------------
C                              into shell stiffness
  160 CONTINUE
C
      DO 200 I=1,12
         L = JPERM(I)
         DO 190 J=1,12
            K = JPERM(J)
            EK(L,K) = PEK(I,J)
  190    CONTINUE
  200 CONTINUE
C
C ---------------------------- Bending stiffness -------------------
C
C -- Elasticity matrix ELAB corresponds to moments Mx,My and Mxy
C -- Elasticity matrix ELAS corresponds to shearing forces Qx og Qy.
C
C ------------------------------------------------------------------
C
      ELAB(1,1)= EMOD*THK*THK*THK/(12*(1.0D0-NY*NY))
      ELAB(2,2)= ELAB(1,1)
      ELAB(1,2)= ELAB(1,1)*NY
      ELAB(2,1)= ELAB(1,2)
      ELAB(3,3)= ELAB(1,1)*(1.0D0-NY)/2
      ELAB(3,1)= 0.0D0
      ELAB(1,3)= 0.0D0
      ELAB(3,2)= 0.0D0
      ELAB(2,3)= 0.0D0
C
      ELAS(1,1)= 5.0D0*EMOD*THK/(12.0D0*(1.0D0+NY))
      ELAS(2,2)= ELAS(1,1)
      ELAS(2,1)= 0.0D0
      ELAS(1,2)= 0.0D0
C
C ---- Rotation angles of the element node
C
      AALF(1)=0.0D0
      AALF(2)=0.0D0
      AALF(3)=0.0D0
      AALF(4)=0.0D0
C
C----- Shear deformation included
C
      ISHEAR = 1
      CALL QBEK31(PEK,AELPAT,XL,YL,AALF,ELAB,ELAS,THK,NY,
     +            ISHEAR,IDN,IPSW,IW,IERR)
      IF (IERR) 1000,300,300
C
C ---------------------------- Add bending stiffness
C                              into shell stiffness ----------------
  300 CONTINUE
      DO 400 I=1,12
         L = IPERM(I)
         DO 350 J=1,12
            K = IPERM(J)
            EK(L,K) = PEK(I,J)
  350    CONTINUE
  400 CONTINUE
C
C ---------------------------- Transform to global coordinates -----
C
      DO 700 MM=1,4
         M = (MM-1)*6
         DO 600 NN=MM,4
            N = (NN-1)*6
C
            DO 450 JJ=1,6
               J = N+JJ
               DO 430 LL=1,2
                  L1 = (LL-1)*3
                  L  = M+L1
                  DO 410 I=1,3
                     II = L1+I
                     PEK(II,JJ) = C(1,I)*EK(L+1,J)+C(2,I)*EK(L+2,J)
     +                            +C(3,I)*EK(L+3,J)
  410             CONTINUE
  430          CONTINUE
  450       CONTINUE
C
            DO 500 I=1,6
               II = M+I
               DO 480 LL=1,2
                  L  = (LL-1)*3
                  L1 = N+L
                  DO 460 J=1,3
                     JJ = L1+J
                     EK(II,JJ) = PEK(I,L+1)*C(1,J)+PEK(I,L+2)
     +                           *C(2,J)+PEK(I,L+3)*C(3,J)
  460             CONTINUE
  480          CONTINUE
  500       CONTINUE
C
  600    CONTINUE
  700 CONTINUE
C
C ---------------------------- Lower triangle of EK (from symmetry)
C
      DO 860 I=1,24
         DO 850 J=I,24
            EK(J,I) = EK(I,J)
  850    CONTINUE
  860 CONTINUE
C
      GO TO 2000
C
C ---------------------------- Error return ------------------------
C
 1000 WRITE(IW,6080)
      WRITE(IW,6090) IDN
C
C ---------------------------- Print ? -----------------------------
C
 2000 IF (IPSW-2) 3000,2020,2010
 2010 WRITE(IW,6100)
      CALL MPRT30(EK,24,24,IW)
 2020 WRITE(IW,6110) IDN
C
 3000 RETURN
C
 6000 FORMAT(///5X,'ENTERING FQS31 FOR ELEMENT',I8,' :')
 6010 FORMAT(/5X,'CORNER THICKNESSES :')
 6020 FORMAT(5X,'H1 =',1PE11.4,5X,'H2 =',1PE11.4,5X,'H3 =',1PE11.4,
     +       5X,'H4 =',1PE11.4)
 6030 FORMAT(/5X,'CORNER COORDINATES :')
 6040 FORMAT(5X,'X1 =',1PE11.4,5X,'X2 =',1PE11.4,5X,'X3 =',1PE11.4,
     +       5X,'X4 =',1PE11.4)
 6050 FORMAT(5X,'Y1 =',1PE11.4,5X,'Y2 =',1PE11.4,5X,'Y3 =',1PE11.4,
     +       5X,'Y4 =',1PE11.4)
 6055 FORMAT(5X,'Z1 =',1PE11.4,5X,'Z2 =',1PE11.4,5X,'Z3 =',1PE11.4,
     +       5X,'Z4 =',1PE11.4)
 6060 FORMAT(/5X,'ELASTICITY MATRIX  :')
 6070 FORMAT(1PE18.7,1P,2E15.7)
 6080 FORMAT(///' *** ERROR RETURN FROM FQS31 ***')
 6090 FORMAT(5X,'ELEMENT',I8,'  IN ERROR')
 6100 FORMAT(/5X,'ELEMENT STIFFNESS MATRIX :')
 6110 FORMAT(///5X,'LEAVING FQS31 FOR ELEMENT',I8)
C
      END
C
      SUBROUTINE FQS32 (EMOD,RNY,E,TH,EX,EY,EZ,EV,SM,SR,
     +                  IDN,IW,IPSW,IERR)
C
C***********************************************************************
C
C     FQS32 COMPUTES THE TOP AND BOTTOM SURFACE STRESSES AT THE FOUR
C     NODES OF A FLAT QUADRILATERAL THIN SHELL ELEMENT, COMPOSED OF
C     A MEMBRANE ELEMENT WITH ROTATIONAL DEGREES OF FREEDOM (B.LUND)
C     AND THE QBESH BENDING ELEMENT (WANG WIUXI).
C     THE SHELL STRESS RESULTANTS ARE ALSO RETURNED IN SR.
C     THE ELEMENT MAY HAVE LINEARLY VARYING THICKNESS AND GENERAL
C     ANISOTROPIC MATERIAL PROPERTIES.
C
C
C     PROGRAMMED BY :  KNUT M. OKSTAD
C     DATE/VERSION  :  03.01.06/1
C
C***********************************************************************
C
      IMPLICIT NONE
C
      INTEGER           IDN,IW,IPSW,IERR
      DOUBLE PRECISION  EMOD,RNY,E(3,3),EX(4),EY(4),EZ(4),TH(4),EV(24)
      DOUBLE PRECISION  SM(3,8),SR(6,4)
      INTEGER           I,J,K,ISHEAR,NS
      DOUBLE PRECISION  THK,X21,Y21,Z21,X31,Y31,Z31,X41,Y41,Z41
      DOUBLE PRECISION  SL21,SL31,SL41,COSG3,SING3,COSG4,SING4
      DOUBLE PRECISION  ELAB(3,3),ELAS(2,2),AALF(4),XL(4),YL(4)
      DOUBLE PRECISION  VML(12),VBL(12),RBF(5,4),ESMB(20,12)
C
C ---------------------------- Print ? -----------------------------
C
      IF (IPSW-2) 30,10,10
   10 WRITE(IW,6000) IDN
      WRITE(IW,6010)
      WRITE(IW,6020) TH(1),TH(2),TH(3),TH(4)
      WRITE(IW,6030)
      WRITE(IW,6040) EX(1),EX(2),EX(3),EX(4)
      WRITE(IW,6050) EY(1),EY(2),EY(3),EY(4)
      WRITE(IW,6055) EZ(1),EZ(2),EZ(3),EZ(4)
      WRITE(IW,6060)
      DO 20 I=1,3
         WRITE(IW,6070) E(I,1),E(I,2),E(I,3)
   20 CONTINUE
C
C --- Average element thickness THK
   30 THK = (TH(1)+TH(2)+TH(3)+TH(4)) * 0.25D0
C
C --- Compute membrane and bendig displacement vectors,
C     VML and VBL, in local coordinate system
      CALL FQS38 (VML,VBL,EV,EX,EY,EZ)
C
C --- BEREGNER LOKALE KNUTEPUNKTSKOORDINATER
C
      X21 = EX(2) - EX(1)
      Y21 = EY(2) - EY(1)
      Z21 = EZ(2) - EZ(1)
      X31 = EX(3) - EX(1)
      Y31 = EY(3) - EY(1)
      Z31 = EZ(3) - EZ(1)
      X41 = EX(4) - EX(1)
      Y41 = EY(4) - EY(1)
      Z41 = EZ(4) - EZ(1)
C
      SL21 = DSQRT(X21*X21 + Y21*Y21 + Z21*Z21)
      SL31 = DSQRT(X31*X31 + Y31*Y31 + Z31*Z31)
      SL41 = DSQRT(X41*X41 + Y41*Y41 + Z41*Z41)
C
C     BEREGNER COSINUS OG SINUS TIL VINKELEN MELLOM
C  -- DEN LOKALE X-AKSEN OG HENHOLDSVIS KN.PUNKTENE 3 (G3) OG 4 (G4) --
C
      COSG3 = (X31*X21 + Y31*Y21 + Z31*Z21)/(SL31*SL21)
      IF (COSG3*COSG3 .LT. 1.0D0) THEN
         SING3 = DSQRT(1.0D0-COSG3*COSG3)
      ELSE
         SING3 = 0.0D0
      END IF
C
      COSG4 = (X41*X21 + Y41*Y21 + Z41*Z21)/(SL41*SL21)
      IF (COSG4*COSG4 .LT. 1.0D0) THEN
         SING4 = DSQRT(1.0D0-COSG4*COSG4)
      ELSE
         SING4 = 0.0D0
      END IF
C
      XL(1) = 0.0D0
      XL(2) = SL21
      XL(3) = SL31*COSG3
      XL(4) = SL41*COSG4
      YL(1) = 0.0D0
      YL(2) = 0.0D0
      YL(3) = SL31*SING3
      YL(4) = SL41*SING4
C
      ELAB(1,1) = EMOD*THK*THK*THK/(12.0D0*(1.0D0-RNY*RNY))
      ELAB(2,2) = ELAB(1,1)
      ELAB(1,2) = ELAB(1,1)*RNY
      ELAB(2,1) = ELAB(1,2)
      ELAB(3,3) = ELAB(1,1)*(1.0D0-RNY)/2.0D0
      ELAB(3,1) = 0.0D0
      ELAB(1,3) = 0.0D0
      ELAB(3,2) = 0.0D0
      ELAB(2,3) = 0.0D0
C
      ELAS(1,1) = 5.0D0*EMOD*THK/(12.0D0*(1.0D0+RNY))
      ELAS(2,2) = ELAS(1,1)
      ELAS(1,2) = 0.0D0
      ELAS(2,1) = 0.0D0
C
      AALF(1) = 0.0D0
      AALF(2) = 0.0D0
      AALF(3) = 0.0D0
      AALF(4) = 0.0D0
C
C --- BEREGNER MEMBRANSPENNINGER DIREKTE (PLANE STRESS : NS = 3)
C
      NS = 3
      CALL QMRF32 (VML,XL,YL,E,RNY,THK,NS,SM,IW,IERR)
      IF (IERR .NE. 0)                    GO TO 1000
C
C --- BEREGNER PLATE SPENNINGS MATRISEN
C
      ISHEAR = 1
      CALL QBEF2 (ESMB,XL,YL,AALF,ELAB,ELAS,THK,RNY,
     +            ISHEAR,IDN,IPSW,IW,IERR)
      IF (IERR .NE. 0)                    GO TO 1000
C
      CALL MAB (ESMB,VBL,RBF,20,12,1,0)
C
      DO 100 i = 1, 4
         j = 4 + i
         DO 50 k = 1, 3
C
C --------- STORE THE NODAL STRESS RESULTANTS
C
            SR(k,i) = SM(k,i) * THK
            SR(3+k,i) = RBF(k,i)
C
C --------- CALCULATE STRESSES IN TOP AND BOTTOM SURFACES
C
            SM(k,j) = SM(k,i) - 6.0D0*RBF(k,i)/(THK*THK)
            SM(k,i) = SM(k,i) + 6.0D0*RBF(k,i)/(THK*THK)
C
   50    CONTINUE
  100 CONTINUE
C
      GO TO 2000
C
C ---------------------------- Error return ------------------------
C
 1000 WRITE(IW,6080)
      WRITE(IW,6090) IDN,IERR
C
C ---------------------------- Print ? -----------------------------
C
 2000 IF (IPSW-2) 3000,2020,2010
 2010 WRITE(IW,6100)
      CALL MPRT30(SM,3,8,IW)
 2020 WRITE(IW,6110) IDN
C
 3000 RETURN
C
 6000 FORMAT(///5X,'ENTERING FQS32 FOR ELEMENT',I8,' :')
 6010 FORMAT(/5X,'CORNER THICKNESSES :')
 6020 FORMAT(5X,'H1 =',1PE11.4,5X,'H2 =',1PE11.4,5X,'H3 =',1PE11.4,
     +       5X,'H4 =',1PE11.4)
 6030 FORMAT(/5X,'CORNER COORDINATES :')
 6040 FORMAT(5X,'X1 =',1PE11.4,5X,'X2 =',1PE11.4,5X,'X3 =',1PE11.4,
     +       5X,'X4 =',1PE11.4)
 6050 FORMAT(5X,'Y1 =',1PE11.4,5X,'Y2 =',1PE11.4,5X,'Y3 =',1PE11.4,
     +       5X,'Y4 =',1PE11.4)
 6055 FORMAT(5X,'Z1 =',1PE11.4,5X,'Z2 =',1PE11.4,5X,'Z3 =',1PE11.4,
     +       5X,'Z4 =',1PE11.4)
 6060 FORMAT(/5X,'ELASTICITY MATRIX  :')
 6070 FORMAT(1PE18.7,1P,2E15.7)
 6080 FORMAT(///' *** ERROR RETURN FROM FQS32 ***')
 6090 FORMAT(5X,'ELEMENT',I8,'  IN ERROR',I4)
 6100 FORMAT(/5X,'ELEMENT STRESSES :')
 6110 FORMAT(///5X,'LEAVING FQS32 FOR ELEMENT',I8)
C
      END
C
      SUBROUTINE FQS33 (F,X,Y,Z,PX,PY,PZ,LOP)
C
C***********************************************************************
C
C     FQS33 DETERMINES A LUMPED LOAD VECTOR  F  (REFERRED TO GLOBAL
C     COORDINATES) DUE TO LINEARLY VARYING DISTRIBUTED SURFACE LOAD
C     (DEFINED IN LOCAL OR GLOBAL COORDINATES) ON THE  FQS  ELEMENT.
C
C
C     PROGRAMMED BY :  KNUT M. OKSTAD
C     DATE/VERSION  :  04.04.08/1
C
C***********************************************************************
C
      IMPLICIT NONE
C
      INTEGER           I,J,K,LOP
      DOUBLE PRECISION  F(24),X(4),Y(4),Z(4),PX(4),PY(4),PZ(4)
      DOUBLE PRECISION  C(3,3),PW(3),S1(3),S2(3),XI,ETA,N1,N2,N3,N4,DS
C
      IF (LOP .LE. 0) THEN
C
C     FIND THE 3*3 MATRIX OF DIRECTION COSINES
C
         CALL DIRCOS (X,Y,Z,C)
C
      END IF
C
C     INITIALIZE THE FORCE VECTOR
C
      DO 20 I = 1, 24
         F(I) = 0.0D0
   20 CONTINUE
C
C     NUMERICAL INTEGRATION OF THE NODAL FORCE VECTOR BY 2X2 GAUSS
C
      DO 200 I = -1, 1, 2
         DO 100 J = -1, 1, 2
            XI  = DBLE(I)/SQRT(3.0D0)
            ETA = DBLE(J)/SQRT(3.0D0)
C
C     FIND THE SURFACE DILATATION, DS = || S1 X S2 ||, WHERE
C     S1 AND S2 ARE THE TWO TANGENT VECTORS d{X}/dXI AND d{X}/dETA
C
            N1 = -0.25D0 * (1.0D0 - ETA)
            N2 =  0.25D0 * (1.0D0 - ETA)
            N3 =  0.25D0 * (1.0D0 + ETA)
            N4 = -0.25D0 * (1.0D0 + ETA)
            S1(1) = N1*X(1) + N2*X(2) + N3*X(3) + N4*X(4)
            S1(2) = N1*Y(1) + N2*Y(2) + N3*Y(3) + N4*Y(4)
            S1(3) = N1*Z(1) + N2*Z(2) + N3*Z(3) + N4*Z(4)
            N1 = -0.25D0 * (1.0D0 - XI)
            N2 = -0.25D0 * (1.0D0 + XI)
            N3 =  0.25D0 * (1.0D0 + XI)
            N4 =  0.25D0 * (1.0D0 - XI)
            S2(1) = N1*X(1) + N2*X(2) + N3*X(3) + N4*X(4)
            S2(2) = N1*Y(1) + N2*Y(2) + N3*Y(3) + N4*Y(4)
            S2(3) = N1*Z(1) + N2*Z(2) + N3*Z(3) + N4*Z(4)
            CALL DVPROD (S1,S2,PW,DS)
C
C     EVALUATE THE SURFACE PRESSURE AT CURRENT POINT
C
            N1 = 0.25D0 * (1.0D0 - XI) * (1.0D0 - ETA)
            N2 = 0.25D0 * (1.0D0 + XI) * (1.0D0 - ETA)
            N3 = 0.25D0 * (1.0D0 + XI) * (1.0D0 + ETA)
            N4 = 0.25D0 * (1.0D0 - XI) * (1.0D0 + ETA)
            PW(1) = (PX(1)*N1 + PX(2)*N2 + PX(3)*N3 + PX(4)*N4) * DS
            PW(2) = (PY(1)*N1 + PY(2)*N2 + PY(3)*N3 + PY(4)*N4) * DS
            PW(3) = (PZ(1)*N1 + PZ(2)*N2 + PZ(3)*N3 + PZ(4)*N4) * DS
C
            IF (LOP .LE. 0) THEN
C
C     TRANSFORM TO GLOBAL COORDINATES
C
               S1(1) = PW(1)
               S1(2) = PW(2)
               S1(3) = PW(3)
               DO 50 K = 1, 3
                  PW(K) = C(1,K)*S1(1) + C(2,K)*S1(2) + C(3,K)*S1(3)
   50          CONTINUE
C
            END IF
C
C     ACCUMULATE NODAL FORCES
C
            DO 60 K = 1, 3
               F(   K) = F(   K) + PW(K)*N1
               F( 6+K) = F( 6+K) + PW(K)*N2
               F(12+K) = F(12+K) + PW(K)*N3
               F(18+K) = F(18+K) + PW(K)*N4
   60       CONTINUE
C
  100    CONTINUE
  200 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE FQS35 (EM,X,Y,Z,TH,RHO)
C
C***********************************************************************
C
C     FQS35 DETERMINES A LUMPED MASS MATRIX  EM  FOR THE
C     QUADRILATERAL SHELL ELEMENT FQS.
C
C
C     PROGRAMMED BY :  KETIL AAMNES
C     DATE/VERSION  :
C
C***********************************************************************
C
      IMPLICIT NONE
C
      DOUBLE PRECISION  EM(24),X(4),Y(4),Z(4)
      DOUBLE PRECISION  TH(4),THK,RHO
C
C
C ---------------------------- Average element thickness THK --------
C
      THK = (TH(1)+TH(2)+TH(3)+TH(4)) * 0.25D0
C
C ---------------------------- Membrane mass matrix -----------------
C
      CALL QMRF35(EM,X,Y,Z,THK,RHO)
C
C ---------------------------- Bending mass matrix ------------------
C
      CALL QBEK35(EM,X,Y,Z,THK,RHO)
C
C ---------------------------- Return -------------------------------
      RETURN
      END
C
      SUBROUTINE FQS38 (VML,VBL,V,X,Y,Z)
C
C***********************************************************************
C
C     FQS38 DETERMINES TWO SUBVECTORS OF LOCAL NODAL POINT PARAMETERS
C     FOR THE  FQS  ELEMENT.  ONE VECTOR (VML) CONTAINS THE LOCAL
C     MEMBRANE PARAMETERS AND THE OTHER (VBL) CONTAINS THE LOCAL BENDING
C     PARAMETERS.  THE TOTAL VECTOR OF GLOBAL NODAL POINT PARAMETERS
C     (V) IS INPUT.
C
C
C     PROGRAMMED BY :  TERJE ROLVAG
C     DATE/VERSION  :  17.08.88 / 1.0
C
C***********************************************************************
C
      IMPLICIT NONE
C
      INTEGER           I,K,L,N
      DOUBLE PRECISION  VML(12),VBL(12),V(24),X(4),Y(4),Z(4)
      DOUBLE PRECISION  A(3),B(3),C(3,3)
C
C     THE 3*3 MATRIX OF DIRECTION COSINES
C
      CALL DIRCOS(X,Y,Z,C)
C
C
C     TRANSFORM V AND TRANSFER ITS ELEMENTS INTO APPROPRIATE LOCATIONS
C     IN VML AND VBL
C
      DO 100 I=1,4
C
      A(1) = V(6*I-5)
      A(2) = V(6*I-4)
      A(3) = V(6*I-3)
      N = -1
C
   20 DO 40 L=1,3
         B(L) = C(L,1)*A(1)
         DO 30 K=2,3
            B(L) = B(L) + C(L,K)*A(K)
   30    CONTINUE
   40 CONTINUE
C
      IF (N) 50,50,80
C
   50 VML(3*I-2) = B(1)
      VML(3*I-1) = B(2)
      VBL(3*I-2) = B(3)
C
      A(1) = V(6*I-2)
      A(2) = V(6*I-1)
      A(3) = V(6*I)
      N = 1
      GO TO 20
C
   80 VBL(3*I-1) = B(1)
      VBL(3*I)   = B(2)
      VML(3*I)   = B(3)
C
  100 CONTINUE
C
      RETURN
      END
C
      subroutine DIRCOS (X,Y,Z,T)
C
C***********************************************************************
C
C     DIRCOS FORMS THE TRANFORMATION MATRIX FROM THE GLOBAL TO THE LOCAL
C     COORDINATE SYSTEM FOR A 4-NODED SHELL ELEMENT.
C
C
C     PROGRAMMED BY :  BJORN HAUGEN
C     DATE/VERSION  :
C
C***********************************************************************
C
      implicit none
C
      double precision X(4),Y(4),Z(4),T(3,3)
      double precision vec1(3),vec2(3),vec3(3),len
      integer          i
C
C     Diagonal 1-3
      vec1(1) = X(3) - X(1)
      vec1(2) = Y(3) - Y(1)
      vec1(3) = Z(3) - Z(1)
C
C     Diagonal 2-4
      vec2(1) = X(4) - X(2)
      vec2(2) = Y(4) - Y(2)
      vec2(3) = Z(4) - Z(2)
C
C     Z-axis as normal to the diagonals 1-3 and 2-4
      call DVPROD (vec1,vec2,vec3,len)
      do i = 1, 3
         T(3,i) = vec3(i)/len
      end do
C
C     Temporary X-axis along side 1-2
      vec1(1) = X(2) - X(1)
      vec1(2) = Y(2) - Y(1)
      vec1(3) = Z(2) - Z(1)
C
C     Y-axis as cross product of Z- and X-axis
      call DVPROD (vec3,vec1,vec2,len)
      do i = 1, 3
         T(2,i) = vec2(i)/len
      end do
C
C     X-axis as cross product of Y- and Z-axis
      call DVPROD (vec2,vec3,vec1,len)
      do i = 1, 3
         T(1,i) = vec1(i)/len
      end do
C
      return
      end
C
      subroutine DVPROD (A,B,C,lenC)
      implicit none
      double precision A(3),B(3),C(3),lenC
C     Form the vector product C = AxB
      C(1) = A(2)*B(3) - A(3)*B(2)
      C(2) = A(3)*B(1) - A(1)*B(3)
      C(3) = A(1)*B(2) - A(2)*B(1)
C     Calculate the length of the vector C
      lenC = sqrt(C(1)*C(1) + C(2)*C(2) + C(3)*C(3))
      return
      end
