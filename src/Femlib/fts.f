C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
C     This file is part of FEDEM - https://openfedem.org
C
      SUBROUTINE FTS31(EK,HH,AKB,E,X,Y,Z,TH,ALPHA,BETA,IDN,IW,IPSW,IERR)
C
C***********************************************************************
C
C     FTS31 GENERATES THE (18*18) STIFFNESS MATRIX  EK  FOR A FLAT
C     TRIANGULAR THIN SHELL ELEMENT COMPOSED OF P.G.BERGAN / C.A.
C     FELIPPA'S MEMBRAN ELEMENT WITH ROTATIONAL DEGREES OF FREEDOM
C     AND THE  TEBA  BENDING ELEMENT. THE ELEMENT MAY HAVE LINEARLY
C     VARYING THICKNESS AND GENERAL ANISOTROPIC MATERIAL PROPERTIES.
C
C
C     PROGRAMMED BY :  KETIL AAMNES
C     DATE/VERSION  :  01.10.85 / 1.0
C
C***********************************************************************
C
      INTEGER           IDN,IW,IPSW,IERR
      INTEGER           I,II,J,JJ,K,L,L1,LL,M,MM,N,NN
      INTEGER           IPERM(9),JPERM(9),LS(12)
      CHARACTER*12      STATUS
      DOUBLE PRECISION  EK(18,18),HH(3,9),AKB(9,9),E(3,3),DM(3,3)
      DOUBLE PRECISION  X(3),Y(3),Z(3),TH(3),ALPHA,BETA
      DOUBLE PRECISION  PEK(9,9),C(3,3),XL(3),YL(3)
      DOUBLE PRECISION  AUX1,COSG,SING,SL21,SL31,THK
      DOUBLE PRECISION  X21,X31,Y21,Y31,Z21,Z31
C
      DATA  IPERM(1),IPERM(2),IPERM(3),IPERM(4),IPERM(5),IPERM(6),
     +      IPERM(7),IPERM(8),IPERM(9) /3,4,5,9,10,11,15,16,17/
      DATA  JPERM(1),JPERM(2),JPERM(3),JPERM(4),JPERM(5),JPERM(6),
     +      JPERM(7),JPERM(8),JPERM(9) /1,2,6,7,8,12,13,14,18/
      DATA  LS(1),LS(2),LS(3),LS(4),LS(5),LS(6),LS(7),LS(8),LS(9)
     +      /1,2,3,4,5,6,7,8,9/
C
C ---------------------------- Print ? -----------------------------
C
      IF (IPSW-2) 30,10,10
   10 WRITE(IW,6000) IDN
      WRITE(IW,6010)
      WRITE(IW,6020) TH(1),TH(2),TH(3)
      WRITE(IW,6030)
      WRITE(IW,6040) X(1),X(2),X(3)
      WRITE(IW,6050) Y(1),Y(2),Y(3)
      WRITE(IW,6055) Z(1),Z(2),Z(3)
      WRITE(IW,6060)
      DO 20 I=1,3
         WRITE(IW,6070) E(I,1),E(I,2),E(I,3)
   20 CONTINUE
C
C ---------------------------- Average element thickness THK -------
C
   30 THK = (TH(1)+TH(2)+TH(3))/3.0D0
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
      AUX1 = 1.0D0 - COSG*COSG
      IF (AUX1) 40,40,50
   40 SING = 0.0D0
      GO TO 60
   50 SING = DSQRT(AUX1)
C
   60 XL(1) = 0.0D0
      XL(2) = SL21
      XL(3) = SL31*COSG
      YL(1) = 0.0D0
      YL(2) = 0.0D0
      YL(3) = SL31*SING
C
      DO 80 I=1,18
         DO 70 J=1,18
            EK(I,J) = 0.0D0
   70    CONTINUE
   80 CONTINUE
C
C ---------------------------- The 3*3 matrix of direction cosines -
C
      CALL DIRC30(C,X,Y,Z)
C
C ---------------------------- Membrane stiffness ------------------
C
      DO 100 I=1,3
         DO 90 J=1,3
            DM(I,J) = THK * E(I,J)
   90    CONTINUE
  100 CONTINUE
C
      DO 120 I=1,9
         DO 110 J=1,9
            PEK(I,J) = 0.0D0
  110    CONTINUE
  120 CONTINUE
C
      CALL TMRF31(XL,YL,DM,ALPHA,BETA,0,LS,HH,PEK,9,STATUS)
      IF (STATUS(1:1).NE.' ') GOTO 1000
C
C ---------------------------- Add membrane stiffness
C                              into shell stiffness ----------------
      DO 200 I=1,9
         L = JPERM(I)
         DO 150 J=1,9
            K = JPERM(J)
            EK(L,K) = PEK(I,J)
  150    CONTINUE
  200 CONTINUE
C
C ---------------------------- Bending stiffness -------------------
C
      CALL TEBA31(PEK,AKB,E,XL,YL,TH,IDN,IW,0,IERR)
      IF (IERR) 1000,300,300
C
C ---------------------------- Add bending stiffness
C                              into shell stiffness ----------------
  300 DO 400 I=1,9
         L = IPERM(I)
         DO 350 J=1,9
            K = IPERM(J)
            EK(L,K) = PEK(I,J)
  350    CONTINUE
  400 CONTINUE
C
C ---------------------------- Transform to global coordinates -----
C
      DO 700 MM=1,3
         M = (MM-1)*6
         DO 600 NN=MM,3
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
      DO 860 I=1,18
         DO 850 J=I,18
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
      IF (IERR) 2000,1100,1100
 1100 IERR = -1
      WRITE(IW,6095) STATUS
C
C ---------------------------- Print ? -----------------------------
C
 2000 IF (IPSW-2) 3000,2020,2010
 2010 WRITE(IW,6100)
      CALL MPRT30(EK,18,18,IW)
 2020 WRITE(IW,6110) IDN
C
 3000 RETURN
C
 6000 FORMAT(///28H ENTERING FTS31 FOR ELEMENT ,I7,2H :)
 6010 FORMAT(/25H     CORNER THICKNESSES :)
 6020 FORMAT(9H     H1 =,1PE11.4,8H    H2 =,1PE11.4,8H    H3 =,1PE11.4)
 6030 FORMAT(/25H     CORNER COORDINATES :)
 6040 FORMAT(9H     X1 =,1PE11.4,8H    X2 =,1PE11.4,8H    X3 =,1PE11.4)
 6050 FORMAT(9H     Y1 =,1PE11.4,8H    Y2 =,1PE11.4,8H    Y3 =,1PE11.4)
 6055 FORMAT(9H     Z1 =,1PE11.4,8H    Z2 =,1PE11.4,8H    Z3 =,1PE11.4)
 6060 FORMAT(/25H     ELASTICITY MATRIX  :)
 6070 FORMAT(1PE18.7,1P,2E15.7)
 6080 FORMAT(///32H *** ERROR RETURN FROM FTS31 ***)
 6090 FORMAT(/13H     ELEMENT ,I7,11H   IN ERROR)
 6095 FORMAT(5X,A)
 6100 FORMAT(/31H     ELEMENT STIFFNESS MATRIX :)
 6110 FORMAT(///27H LEAVING FTS31 FOR ELEMENT ,I7)
C
      END
      SUBROUTINE FTS32(SMM,SMB,HH,AKB,X,Y,Z,DM,ALPHA,IDN,IW,IERR)
C
C***********************************************************************
C
C     FTS32 DETERMINES TWO STRESS MATRICES  SMM  AND  SMB  YIELDING THE
C     LOCAL MEMBRANE FORCES AND THE LOCAL MOMENTS, RESPECTIVELY, AT THE
C     CENTROID OF THE  FTS  ELEMENT
C
C
C     PROGRAMMED BY :  KETIL AAMNES
C     DATE/VERSION  :  151085 / 1.0
C
C***********************************************************************
C
      INTEGER           I,IDN,IW,IERR
      DOUBLE PRECISION  AUX1,COSG,SING,SL21,SL31,X21,X31
      DOUBLE PRECISION  Y21,Y31,Z21,Z31
      DOUBLE PRECISION  SMM(3,9),SMB(3,9),HH(3,9),AKB(9,9)
      DOUBLE PRECISION  X(3),Y(3),Z(3),XL(3),YL(3),ZZ(3)
      DOUBLE PRECISION  DM(3,3),ALPHA
      CHARACTER*12      STATUS
C
      IERR = 0
C
C ---------------------------- Natural (area) coordinates
C                              of the centroid ---------------------
      DO 20 I=1,3
         ZZ(I) = 1.0D0/3.0D0
   20 CONTINUE
C
C ---------------------------- Local coordinates -------------------
C
      X21 = X(2)-X(1)
      Y21 = Y(2)-Y(1)
      Z21 = Z(2)-Z(1)
      X31 = X(3)-X(1)
      Y31 = Y(3)-Y(1)
      Z31 = Z(3)-Z(1)
      SL21 =  DSQRT(X21*X21 + Y21*Y21 + Z21*Z21)
      SL31 =  DSQRT(X31*X31 + Y31*Y31 + Z31*Z31)
      COSG = (X31*X21 + Y31*Y21 + Z31*Z21)/(SL31*SL21)
      AUX1 = 1.0D0 - COSG*COSG
      IF (AUX1) 40,40,50
   40 SING = 0.0D0
      GO TO 60
   50 SING =  DSQRT(AUX1)
C
   60 XL(1) = 0.0D0
      XL(2) = SL21
      XL(3) = SL31*COSG
      YL(1) = 0.0D0
      YL(2) = 0.0D0
      YL(3) = SL31*SING
C
C ---------------------------- The membrane stress matrix  SMM -----
C
      CALL TMRF32(HH,XL,YL,DM,ALPHA,SMM,STATUS)
      IF (STATUS(1:1).NE.' ') GOTO 1000
C
C ---------------------------- The bending stress matrix  SMB ------
C
      CALL TEBA32(SMB,AKB,ZZ)
C
      GO TO 2000
C
C ---------------------------- Error return ------------------------
C
 1000 IERR = -1
      WRITE(IW,6080)
      WRITE(IW,6090) IDN
      WRITE(IW,6095) STATUS
C
 2000 RETURN
C
 6080 FORMAT(///32H *** ERROR RETURN FROM FTS32 ***)
 6090 FORMAT(/13H     ELEMENT ,I7,11H   IN ERROR)
 6095 FORMAT(5X,A)
C
      END
      SUBROUTINE FTS33(F,X,Y,Z,PX,PY,PZ,LOP)
C
C***********************************************************************
C
C     FTS33 DETERMINES A LUMPED LOAD VECTOR  F  (REFERRED TO GLOBAL
C     COORDINATES) DUE TO LINEARLY VARYING DISTRIBUTED SURFACE LOAD
C     (DEFINED IN LOCAL OR GLOBAL COORDINATES) ON THE  FTS  ELEMENT
C
C
C     PROGRAMMED BY :  KETIL AAMNES
C     DATE/VERSION  :  151085 / 1.0
C
C***********************************************************************
C
      INTEGER           I,J,K,L,M,N,LOP,IP(3)
      DOUBLE PRECISION  A,AUX,AUX1,COSG,SING,SL21,SL31
      DOUBLE PRECISION  X21,X31,Y21,Y31,Z21,Z31
      DOUBLE PRECISION  F(18),X(3),Y(3),Z(3),PX(3),PY(3),PZ(3)
      DOUBLE PRECISION  XL(3),YL(3),C(3,3),P1(3),P2(3),P3(3),B(3)
C
      DATA IP(1),IP(2),IP(3)/2,3,1/
C
      DO 20 I=1,18
      F(I) = 0.0D0
   20 CONTINUE
C
C     LOCAL COORDINATES
C
      X21 = X(2)-X(1)
      Y21 = Y(2)-Y(1)
      Z21 = Z(2)-Z(1)
      X31 = X(3)-X(1)
      Y31 = Y(3)-Y(1)
      Z31 = Z(3)-Z(1)
      SL21 =  DSQRT(X21*X21 + Y21*Y21 + Z21*Z21)
      SL31 =  DSQRT(X31*X31 + Y31*Y31 + Z31*Z31)
      COSG = (X31*X21 + Y31*Y21 + Z31*Z21)/(SL31*SL21)
      AUX1 = 1.0D0 - COSG*COSG
      IF (AUX1) 40,40,50
   40 SING = 0.0D0
      GO TO 60
   50 SING =  DSQRT(AUX1)
C
   60 XL(1) = 0.0D0
      XL(2) = SL21
      XL(3) = SL31*COSG
      YL(1) = 0.0D0
      YL(2) = 0.0D0
      YL(3) = SL31*SING
C
C     ELEMENT AREA (DIVIDED BY 12)
C
      A = XL(1)*YL(2) + XL(2)*YL(3) + XL(3)*YL(1)
      A = A - XL(1)*YL(3) - XL(2)*YL(1) - XL(3)*YL(2)
      A = A/24.0D0
C
C     THE 3*3 MATRIX OF DIRECTION COSINES
C
      CALL DIRC30(C,X,Y,Z)
C
      IF (LOP) 80,80,120
C
C     TRANSFER LOCAL LOAD INTENSITIES
C
   80 DO 100 I=1,3
      P1(I) = PX(I)
      P2(I) = PY(I)
      P3(I) = PZ(I)
  100 CONTINUE
      GO TO 150
C
C     TRANSFORM FROM GLOBAL TO LOCAL LOAD INTENSITIES
C
  120 DO 140 I=1,3
      P1(I) = C(1,1)*PX(I) + C(1,2)*PY(I) + C(1,3)*PZ(I)
      P2(I) = C(2,1)*PX(I) + C(2,2)*PY(I) + C(2,3)*PZ(I)
      P3(I) = C(3,1)*PX(I) + C(3,2)*PY(I) + C(3,3)*PZ(I)
  140 CONTINUE
C
C     DETERMINE NON-ZERO ELEMENTS OF THE LOCAL LOAD VECTOR, TRANSFORM
C     AND ADD INTO APPROPRIATE LOCATIONS IN THE GLOBAL LOAD VECTOR
C
  150 DO 200 I=1,3
      J = IP(I)
      K = IP(J)
      B(1) = A*(2.0D0*P1(I) + P1(J) + P1(K))
      B(2) = A*(2.0D0*P2(I) + P2(J) + P2(K))
      B(3) = A*(2.0D0*P3(I) + P3(J) + P3(K))
C
      DO 180 M=1,3
      AUX = C(1,M)*B(1)
      DO 170 N=2,3
      AUX = AUX + C(N,M)*B(N)
  170 CONTINUE
      L = 6*I-6+M
      F(L) = AUX
  180 CONTINUE
  200 CONTINUE
C
      RETURN
      END
      SUBROUTINE FTS35(EM,X,Y,Z,TH,RHO)
C
C***********************************************************************
C
C     FTS35 DETERMINES A LUMPED MASS MATRIX  EM  FOR THE TRIANGULAR
C     SHELL ELEMENT  FTS
C
C
C     PROGRAMMED BY :  KETIL AAMNES
C     DATE/VERSION  :  081085 / 1.0
C
C***********************************************************************
C
      DOUBLE PRECISION  EM(18),EMM(9),EMP(9),X(3),Y(3),Z(3)
      DOUBLE PRECISION  TH(3),THK,RHO
C
C ---------------------------- Average element thickness THK --------
C
      THK = (TH(1)+TH(2)+TH(3))/3.0D0
C
C ---------------------------- Membrane mass matrix -----------------
C
      CALL TMRF35(EMM,X,Y,Z,THK,RHO)
C
C ---------------------------- Bending mass matrix ------------------
C
      CALL TEBA35(EMP,X,Y,Z,THK,RHO)
C
C ---------------------------- Add into SHELL mass matrix -----------
C
      EM( 1) = EMM(1)
      EM( 2) = EMM(2)
      EM( 3) = EMP(1)
      EM( 4) = EMP(2)
      EM( 5) = EMP(3)
      EM( 6) = EMM(3)
C
      EM( 7) = EM(1)
      EM(13) = EM(1)
      EM( 8) = EM(2)
      EM(14) = EM(2)
      EM( 9) = EM(3)
      EM(15) = EM(3)
      EM(10) = EM(4)
      EM(16) = EM(4)
      EM(11) = EM(5)
      EM(17) = EM(5)
      EM(12) = EM(6)
      EM(18) = EM(6)
C
      RETURN
      END
      SUBROUTINE FTS38(VML,VBL,V,X,Y,Z)
C
C***********************************************************************
C
C     FTS38 DETERMIBES TWO SUBVECTORS OF LOCAL NODAL POINT PARAMETERS
C     FOR THE  FTS  ELEMENT.  ONE VECTOR (VML) CONTAINS THE LOCAL
C     MEMBRANE PARAMETERS AND THE OTHER (VBL) CONTAINS THE LOCAL BENDING
C     PARAMETERS.  THE TOTAL VECTOR OF GLOBAL NODAL POINT PARAMETERS
C     (V) IS INPUT
C
C
C     PROGRAMMED BY :  KETIL AAMNES
C     DATE/VERSION  :  151085 / 1.0
C
C***********************************************************************
C
      INTEGER           I,K,L,N
      DOUBLE PRECISION  VML(9),VBL(9),V(18),X(3),Y(3),Z(3)
      DOUBLE PRECISION  A(3),B(3),C(3,3)
C
C     THE 3*3 MATRIX OF DIRECTION COSINES
C
      CALL DIRC30(C,X,Y,Z)
C
C     TRANSFORM V AND TRANSFER ITS ELEMENTS INTO APPROPRIATE LOCATIONS
C     IN VML AND VBL
C
      DO 100 I=1,3
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
   30 CONTINUE
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
