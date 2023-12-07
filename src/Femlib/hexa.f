C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
C     This file is part of FEDEM - https://openfedem.org
C
      SUBROUTINE HEXA31(EK,A1,A2,E,X,Y,Z,NIP,IOP,ID,IW,IPSW,IERR)
C
C***********************************************************************
C
C  FINITE ELEMENT LIBRARY : HEXA31                           SINTEF/NTH
C
C
C     HEXA31 GENERATES THE STIFFNESS MATRIX FOR AN INCOMPATIBLE THREE
C     DIMENSIONAL SOLID HEXAHEDRON ELEMENT WITH 8 NODAL POINTS AND 24
C     DEGREES OF FREEDOM (3 AT EACH NODAL POINT).
C     THE ELEMENT HAS 9 INTERNAL (NODELESS) DEGREES OF FREEDOM (IOP=1).
C     IF  IOP=0  THE INTERNAL DEGREES OF FREEDOM ARE NOT INTRODUCED AND
C     THE ELEMENT REDUCES TO THE STANDARD COMPATIBLE 24-DEGREE-OF-
C     FREEDOM ELEMENT
C     THE ELEMENT MAY HAVE GENERAL, ANISOTROPIC MATERIAL PROPERTIES.
C
C     PROGRAMMED BY: S.L\SETH
C
C     DATE/VERSION : 01.05.74 / 01
C     DATE/VERSION : 83-07-06 / 1.1  REVISION BY S. L\SETH
C     DATE/VERSION : 88-02-05 / 1.2  REVISION BY K. AAMNES
C     DATE/VERSION : 01-03-17 / 1.3  REVISION BY K. OKSTAD
C     DATE/VERSION : 10-05-21 / 1.4  REVISION BY K. OKSTAD
C
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION EK(24,24),A1(24,9),A2(9,9),E(6,6),X(8),Y(8),Z(8),ABC(3),
     +          W(3),D(3,11),BJ(3,3),EP(6,33),EKT(33,33),CXI(8),CETA(8),
     +          CZETA(8),IWORK(9)
C
      DATA (CXI(I),I=1,8)  /-1., 1., 1.,-1.,-1., 1., 1.,-1./
      DATA (CETA(I),I=1,8) /-1.,-1., 1., 1.,-1.,-1., 1., 1./
      DATA (CZETA(I),I=1,8)/-1.,-1.,-1.,-1., 1., 1., 1., 1./
C
C     PRINT ?
C
      IF (IPSW-2) 30,10,10
C
   10 WRITE(IW,6000)ID
      WRITE(IW,6010)
      WRITE(IW,6020)
      WRITE(IW,6030) X
      WRITE(IW,6040) Y
      WRITE(IW,6050) Z
      WRITE(IW,6060)
      DO 15 I=1,6
   15 WRITE(IW,6070)(E(I,J),J=1,6)
C
   30 IERR=0
C
C     PARAMETER NIP OUT OF RANGE ?
C
      IF (NIP-2) 1000,75,50
   50 IF (NIP-3) 1000,75,1000
C
C     CHECK THE JACOBIAN DETERMINANT AT THE NODES
C
   75 DO 80 I=1,8
      XI=CXI(I)
      ETA=CETA(I)
      ZETA=CZETA(I)
      CALL JABN30(BJ,X,Y,Z,XI,ETA,ZETA,DETJ,0)
      IF (DETJ) 1100,80,80
   80 CONTINUE
C
C     ABSCISSAS AND WEIGHT COEFFICIENTS OF THE INTEGRATION FORMULA
C
      IF (NIP-2) 1000,85,90
C
   85 ABC(1)=-SQRT(1.0/3.0)
      ABC(2)=-ABC(1)
C
      W(1)=1.0
      W(2)=1.0
      GO TO 100
C
   90 ABC(1)=-SQRT(0.6)
      ABC(2)= 0.0
      ABC(3)=-ABC(1)
C
      W(1)=5.0/9.0
      W(2)=8.0/9.0
      W(3)=W(1)
C
C     INTERNAL DEGREES OF FREEDOM ?
C
  100 IF (IOP) 105,105,110
C
C     24 DOF
C
  105 NDOF=24
      NODE=8
      GO TO 120
C
C     33 DOF
C
  110 NDOF=33
      NODE=11
C
  120 DO 130 I=1,NDOF
      DO 130 J=1,NDOF
  130 EKT(I,J)=0.
C
C
C     START OF THE INTEGRATION LOOP
C
      DO 500 K=1,NIP
      ZETA=ABC(K)
      ZETP=1.+ZETA
      ZETM=1.-ZETA
      WZ=W(K)
C
      DO 400 J=1,NIP
      ETA=ABC(J)
      ETAP=1.+ETA
      ETAM=1.-ETA
      WY=W(J)
C
      DO 300 I=1,NIP
      XI=ABC(I)
      XIP=1.+XI
      XIM=1.-XI
      WI=WZ*WY*W(I)
C
C     THE INVERSE OF THE JACOBIAN MATRIX
C
      CALL JABN30(BJ,X,Y,Z,XI,ETA,ZETA,DETJ,1)
C
C     THE DERIVATIVES OF THE SHAPE FUNCTION WITH RESPECT TO GLOBAL
C     COORDINATES
C
      DO 150 N=1,3
C
      D(N,1)=.125*(-BJ(N,1)*ETAM*ZETM-BJ(N,2)*XIM*ZETM-BJ(N,3)*XIM*ETAM)
      D(N,2)=.125*( BJ(N,1)*ETAM*ZETM-BJ(N,2)*XIP*ZETM-BJ(N,3)*XIP*ETAM)
      D(N,3)=.125*( BJ(N,1)*ETAP*ZETM+BJ(N,2)*XIP*ZETM-BJ(N,3)*XIP*ETAP)
      D(N,4)=.125*(-BJ(N,1)*ETAP*ZETM+BJ(N,2)*XIM*ZETM-BJ(N,3)*XIM*ETAP)
      D(N,5)=.125*(-BJ(N,1)*ETAM*ZETP-BJ(N,2)*XIM*ZETP+BJ(N,3)*XIM*ETAM)
      D(N,6)=.125*( BJ(N,1)*ETAM*ZETP-BJ(N,2)*XIP*ZETP+BJ(N,3)*XIP*ETAM)
      D(N,7)=.125*( BJ(N,1)*ETAP*ZETP+BJ(N,2)*XIP*ZETP+BJ(N,3)*XIP*ETAP)
      D(N,8)=.125*(-BJ(N,1)*ETAP*ZETP+BJ(N,2)*XIM*ZETP+BJ(N,3)*XIM*ETAP)
C
      IF (IOP) 150,150,140
C
  140 D(N,9)=-2.*BJ(N,1)*XI
      D(N,10)=-2.*BJ(N,2)*ETA
      D(N,11)=-2.*BJ(N,3)*ZETA
  150 CONTINUE
C
C     THE MATRIX E*P MULTIPLIED BY DETJ
C
      DO 200 M=1,NODE
      DO 200 N=1,6
      EP(N,3*M-2)=DETJ*(E(N,1)*D(1,M)+E(N,4)*D(2,M)+E(N,5)*D(3,M))
      EP(N,3*M-1)=DETJ*(E(N,2)*D(2,M)+E(N,4)*D(1,M)+E(N,6)*D(3,M))
  200 EP(N,3*M)=DETJ*(E(N,3)*D(3,M)+E(N,5)*D(1,M)+E(N,6)*D(2,M))
C
C     UPPER TRIANGULAR PART OF THE TOTAL 33*33 (OR 24*24) STIFFNESS
C     MATRIX
C
      DO 250 L=1,NODE
      M=3*L-2
      DO 250 N=M,NDOF
      EKT(M,N)=EKT(M,N)+WI*(D(1,L)*EP(1,N)+D(2,L)*EP(4,N)
     ++D(3,L)*EP(5,N))
      EKT(M+1,N)=EKT(M+1,N)+WI*(D(2,L)*EP(2,N)+D(1,L)*EP(4,N)
     ++D(3,L)*EP(6,N))
      EKT(M+2,N)=EKT(M+2,N)+WI*(D(3,L)*EP(3,N)+D(1,L)*EP(5,N)
     ++D(2,L)*EP(6,N))
  250 CONTINUE
C
C     END OF INTEGRATION LOOP
C
  300 CONTINUE
  400 CONTINUE
  500 CONTINUE
C
C     LOWER PART OF TOTAL STIFFNESS MATRIX FROM SYMMETRY
C
      DO 600 I=1,NDOF
      DO 600 J=I,NDOF
  600 EKT(J,I)=EKT(I,J)
C
      IF (IOP) 605,605,620
C
C     NO INTERNAL DEGREES OF FREEDOM
C
  605 DO 610 I=1,24
      DO 610 J=1,24
  610 EK(I,J)=EKT(I,J)
      GO TO 2200
C
C     COMPUTE A2=INVERSE OF KII
C
  620 DO 700 I=1,9
      DO 700 J=1,9
  700 A2(I,J)=EKT(I+24,J+24)
C
      CALL DINV12(9,A2,A1,IWORK,IW,IERR)
      IF (IERR) 1200,750,1200
C
C     COMPUTE A1=KEI*A2
C
  750 DO 800 I=1,9
      DO 800 J=1,24
      A1(J,I)=0.0
      DO 800 K=1,9
  800 A1(J,I)=A1(J,I)+EKT(J,K+24)*A2(K,I)
C
C     THE STIFFNESS MATRIX EK OBTAINED BY STATIC CONDENSATION
C
      DO 900 I=1,24
      DO 900 J=1,I
      EK(J,I)=EKT(J,I)
      DO 900 K=1,9
  900 EK(J,I)=EK(J,I)-A1(J,K)*EKT(24+K,I)
      DO 950 I=1,24
      DO 950 J=I,24
      EK(J,I)=EK(I,J)
  950 CONTINUE
C
      GO TO 2200
C
C     ERROR MESSAGE
C
 1000 IERR=-1
      GO TO 2000
C
 1100 IERR=-2
      WRITE(IW,6101) I,DETJ
      WRITE(IW,6030) X
      WRITE(IW,6040) Y
      WRITE(IW,6050) Z
      GO TO 2000
C
 1200 IERR=-3
C
 2000 WRITE(IW,6080)
      WRITE(IW,6090) IERR
      WRITE(IW,6100) ID
      IF (IPSW-2) 3000,2500,2500
C
C     PRINT ?
C
 2200 IF (IPSW-2) 3000,2500,2300
 2300 WRITE(IW,6110)
      CALL MPRT30(EK,24,24,IW)
 2500 WRITE(IW,6120)ID
C
 3000 RETURN
C
C
 6000 FORMAT(///33X,29H ENTERING HEXA31 FOR ELEMENT ,I7,2H :)
 6010 FORMAT(/39X,25H NODAL POINT COORDINATES:)
 6020 FORMAT(/10H NODE NO: ,11X,1H1,10X,1H2,10X,1H3,10X,1H4,10X,1H5,10X,
     +1H6,10X,1H7,10X,1H8)
 6030 FORMAT(15H X-COORDINATES:,1P,8E11.4)
 6040 FORMAT(15H Y-COORDINATES:,1P,8E11.4)
 6050 FORMAT(15H Z-COORDINATES:,1P,8E11.4)
 6060 FORMAT(//42X,20H ELASTICITY MATRIX: /)
 6070 FORMAT(16X,1P,6E11.4)
 6080 FORMAT(///33H *** ERROR RETURN FROM HEXA31 ***)
 6090 FORMAT(/16H     ERROR FLAG=,I3)
 6100 FORMAT(/13H     ELEMENT ,I7,11H   IN ERROR)
 6101 FORMAT(/38H     DEGENERATED ELEMENT, NODE, DETJ =,I3,1PE13.5 )
 6110 FORMAT(/35X,27H ELEMENT STIFFNESS MATRIX :)
 6120 FORMAT(///35X,28H LEAVING HEXA31 FOR ELEMENT ,I7,2H :)
C
      END
      SUBROUTINE HEXA32(SI,B,A1,A2,DETJ,E,X,Y,Z,XI,ETA,ZETA,IOP,LOP,IER)
C
C***********************************************************************
C
C   FINITE ELEMENT LIBRARY : HEXA32                           NTH/SINTEF
C
C
C     HEXA32 COMPUTES A STRESS MATRIX 'SI' AND (IF IOP=1) A STRESS
C     CORRECTION MATRIX 'B' FOR A COMPATIBLE (IOP=0) OR INCOMPATIBLE
C     (IOP=1) SOLID HEXAHEDRON ELEMENT WITH 8 NODAL POINTS AND 24
C     DEGREES OF FREEDOM (PLUS, IF IOP=1, 9 INTERNAL DOFS)
C
C     PROGRAMMED BY: S.L\SETH
C
C     DATE/VERSION : 13.05.74 / 01
C     DATE/VERSION : 83-07-06 / 1.1  REVISION BY S. L\SETH
C     DATE/VERSION : 88-02-05 / 1.2  REVISION BY K. AAMNES
C     DATE/VERSION : 01-03-17 / 1.3  REVISION BY K. OKSTAD
C     DATE/VERSION : 10-05-21 / 1.4  REVISION BY K. OKSTAD
C
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION SI(6,24),B(6,9),A1(24,9),A2(9,9),E(6,6),X(8),Y(8),Z(8),
     +          BJ(3,3),P(6,24),PA(6,9),CXI(8),CETA(8),CZETA(8)
C
      DATA (CXI(I),I=1,8)  /-1., 1., 1.,-1.,-1., 1., 1.,-1./
      DATA (CETA(I),I=1,8) /-1.,-1., 1., 1.,-1.,-1., 1., 1./
      DATA (CZETA(I),I=1,8)/-1.,-1.,-1.,-1., 1., 1., 1., 1./
C
      IER=0
C
C     COORDINATES OF THE 'STRESS POINT'
C
      IF(LOP-2) 20,30,10
   10 IF(LOP-3) 20,40,20
C
   20 ABC = 1.0
      GO TO 50
C
   30 ABC = SQRT(1.0/3.0)
      GO TO 50
C
   40 ABC = SQRT(0.6)
C
   50 XXI   = ABC*XI
      XETA  = ABC*ETA
      XZETA = ABC*ZETA
C
      DO 60 I=1,6
      DO 60 J=1,24
      P(I,J)  = 0.
   60 SI(I,J) = 0.
C
      IF(IOP) 75,75,65
C
   65 DO 70 I=1,6
      DO 70 J=1,9
   70 B(I,J) = 0.
C
C     THE INVERSE OF THE JACOBIAN MATRIX
C
   75 CALL JABN30(BJ,X,Y,Z,XXI,XETA,XZETA,DETJ,1)
      IF (DETJ) 80,100,100
   80 IER=1
      GO TO 200
C
C     THE MATRIX 'P'
C
  100 DO 120 I=1,8
      DO 110 J=1,3
      K=3*I+J-3
      P(J,K)=0.125*(BJ(J,1)*CXI(I)*(1.+XETA*CETA(I))*(1.+XZETA*CZETA(I))
     +             +BJ(J,2)*CETA(I)*(1.+XXI*CXI(I))*(1.+XZETA*CZETA(I))
     +             +BJ(J,3)*CZETA(I)*(1.+XXI*CXI(I))*(1.+XETA*CETA(I)))
  110 CONTINUE
      L=3*I-2
      P(4,L)=P(2,L+1)
      P(4,L+1)=P(1,L)
      P(5,L)=P(3,L+2)
      P(5,L+2)=P(1,L)
      P(6,L+1)=P(3,L+2)
      P(6,L+2)=P(2,L+1)
  120 CONTINUE
C
      IF (IOP) 180,180,125
C
C     THE MATRIX 'PA'
C
  125 DO 130 J=1,3
      PA(J,J)=-2.*BJ(J,1)*XXI
      PA(J,J+3)=-2.*BJ(J,2)*XETA
  130 PA(J,J+6)=-2.*BJ(J,3)*XZETA
C
      DO 140 I=1,7,3
      PA(4,I)=PA(2,I+1)
      PA(4,I+1)=PA(1,I)
      PA(5,I)=PA(3,I+2)
      PA(5,I+2)=PA(1,I)
      PA(6,I+1)=PA(3,I+2)
  140 PA(6,I+2)=PA(2,I+1)
C
C     STRESS CORRECTION MATRIX 'B'
C
      DO 150 I=1,6
      DO 150 J=1,9
      DO 150 K=1,9
  150 SI(I,J) = SI(I,J) + PA(I,K)*A2(K,J)
C
      DO 160 I=1,6
      DO 160 J=1,9
      DO 160 K=1,6
  160 B(I,J) = B(I,J) + E(I,K)*SI(K,J)
C
C     MODIFICATION OF 'P' DUE TO STATIC CONDENSATION
C
      DO 170 J=1,24
      DO 170 K=1,6
      DO 170 L=1,9
  170 P(K,J)=P(K,J)-PA(K,L)*A1(J,L)
C
      DO 175 I=1,6
      DO 175 J=1,9
  175 SI(I,J) = 0.
C
C     STRESS MATRIX 'SI'
C
  180 DO 190 I=1,6
      DO 190 J=1,24
      DO 190 K=1,6
  190 SI(I,J)=SI(I,J)+E(I,K)*P(K,J)
C
  200 RETURN
C
      END
      SUBROUTINE HEXA33(F,FI,A1,X,Y,Z,PX,PY,PZ,ISD,IOP)
C
C***********************************************************************
C
C
C   FINITE ELEMENT LIBRARY  :  HEXA33                         SINTEF/NTH
C
C
C     HEXA33 COMPUTES A CONSISTENT LOAD VECTOR FOR A COMPATIBLE (IOP=0)
C     OR AN INCOMPATIBLE (IOP=1) THREE-DIMENSIONAL, SOLID HEXAHEDRON
C     ELEMENT WITH 8 NODAL POINTS AND 24 DEGREES OF FREEDOM
C     THE LOAD VECTOR IS DUE TO CONSTANT VOLUME FORCES PX,PY,PZ OR
C     CONSTANT SURFACE FORCES PX,PY,PZ ACTING ON ONE SURFACE OF THE
C     ELEMENT
C
C     PROGRAMMED BY : S.L\SETH
C
C     DATE/VERSION : 27.08.74 / 01
C     DATE/VERSION : 83-07-06 / 1.1  REVISION BY S. L\SETH
C     DATE/VERSION : 88-02-05 / 1.2  REVISION BY K. AAMNES
C     DATE/VERSION : 01-03-17 / 1.3  REVISION BY K. OKSTAD
C     DATE/VERSION : 10-05-21 / 1.4  REVISION BY K. OKSTAD
C
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION F(24),FI(9),A1(24,9),X(8),Y(8),Z(8),ABC(2),D(11),CXI(8),
     +          CETA(8),CZETA(8),NR(6,4),DT(4),DS(4),BJ(3,3)
C
      DATA (CXI(I),I=1,8)  /-1., 1., 1.,-1.,-1., 1., 1.,-1./
      DATA (CETA(I),I=1,8) /-1.,-1., 1., 1.,-1.,-1., 1., 1./
      DATA (CZETA(I),I=1,8)/-1.,-1.,-1.,-1., 1., 1., 1., 1./
C
      DATA ((NR(I,J),J=1,4),I=1,6)/
     1   5,6,7,8,
     2   3,4,7,8,
     3   1,4,5,8,
     4   1,2,5,6,
     5   2,3,6,7,
     6   1,2,3,4/
C
      IERR=0
C
   10 DO 20 I=1,24
   20 F(I)=0.
C
      IF (IOP) 35,35,25
C
   25 DO 30 I=1,9
   30 FI(I)=0.
C
C
   35 IF (IERR) 3000,40,40
C
C     ABSCISSAS OF THE GAUSSIAN QUADRATURE FORMULA
C
   40 ABC(1)=-SQRT(1.0/3.0)
      ABC(2)=-ABC(1)
C
C
C     VOLUME FORCES ?
C
      IF (ISD) 50,50,55
   55 IF (ISD-6) 1000,1000,3000
C
C     START OF THE INTEGRATION LOOP
C
   50 DO 500 I=1,2
      XI=ABC(I)
      DO 500 J=1,2
      ETA=ABC(J)
      DO 500 K=1,2
      ZETA=ABC(K)
C
C     THE DETERMINANT OF THE JACOBIAN MATRIX
C
      CALL JABN30(BJ,X,Y,Z,XI,ETA,ZETA,DETJ,0)
      IF (DETJ) 60,60,70
   60 IERR=-1
      GO TO 10
C
C     THE SHAPE FUNCTIONS
C
   70 DO 80 N=1,8
   80 D(N)=0.125*(1.+CXI(N)*XI)*(1.+CETA(N)*ETA)*(1.+CZETA(N)*ZETA)
      D(9)=1.-XI*XI
      D(10)=1.-ETA*ETA
      D(11)=1.-ZETA*ZETA
C
C     EXTERNAL AND INTERNAL PART OF THE TOTAL LOAD VECTOR
C
      DO 100 N=1,8
  100 F(3*N)=F(3*N)+D(N)*DETJ
C
      IF (IOP) 500,500,110
C
  110 DO 120 N=1,3
      M=N+8
  120 FI(3*N)=FI(3*N)+D(M)*DETJ
C
  500 CONTINUE
C
C
      DO 510 N=1,8
      F(3*N-2)=F(3*N)*PX
      F(3*N-1)=F(3*N)*PY
  510 F(3*N)=F(3*N)*PZ
C
      IF (IOP) 3000,3000,520
C
  520 DO 530 N=1,3
      FI(3*N-2)=FI(3*N)*PX
      FI(3*N-1)=FI(3*N)*PY
  530 FI(3*N)=FI(3*N)*PZ
      GO TO 2500
C
C
C     SURFACE FORCES
C
C
 1000 DO 1010 I=1,11
 1010 D(I)=0.
C
C     START OF INTEGRATION LOOP
C
      DO 2000 I=1,2
      S=ABC(I)
      DO 2000 J=1,2
      T=ABC(J)
      DTX=0.
      DTY=0.
      DTZ=0.
      DSX=0.
      DSY=0.
      DSZ=0.
C
C     WHICH SURFACE IS LOADED ?
C
      GO TO(1050,1100,1150,1100,1150,1050),ISD
C
C     SURFACE NO.1 OR NO.6.DERIVATIVES AND SHAPE FUNCTIONS
C
 1050 DO 1060 K=1,4
      N=NR(ISD,K)
      DT(K)=0.25*CXI(N)*(1.+CETA(N)*S)
      DS(K)=0.25*CETA(N)*(1.+CXI(N)*T)
      DTX=DTX+DT(K)*X(N)
      DTY=DTY+DT(K)*Y(N)
      DTZ=DTZ+DT(K)*Z(N)
      DSX=DSX+DS(K)*X(N)
      DSY=DSY+DS(K)*Y(N)
      DSZ=DSZ+DS(K)*Z(N)
 1060 D(N)=0.25*         (1.+CXI(N)*T)*(1.+CETA(N)*S)
      D(9)=1.-T*T
      D(10)=1.-S*S
      GO TO 1200
C
C     SURFACE NO.2 OR NO.4. DERIVATIVES AND SHAPE FUNCTIONS
C
 1100 DO 1110 K=1,4
      N=NR(ISD,K)
      DT(K)=0.25*CXI(N)*(1.+CZETA(N)*S)
      DS(K)=0.25*CZETA(N)*(1.+CXI(N)*T)
      DTX=DTX+DT(K)*X(N)
      DTY=DTY+DT(K)*Y(N)
      DTZ=DTZ+DT(K)*Z(N)
      DSX=DSX+DS(K)*X(N)
      DSY=DSY+DS(K)*Y(N)
      DSZ=DSZ+DS(K)*Z(N)
 1110 D(N)=0.25*        (1.+CXI(N)*T)*(1.+CZETA(N)*S)
      D(9)=1.-T*T
      D(11)=1.-S*S
      GO TO 1200
C
C     SURFACE NO.3 OR NO.5. DERIVATIVES AND SHAPE FUNCTIONS
C
 1150 DO 1160 K=1,4
      N=NR(ISD,K)
      DT(K)=0.25*CETA(N)*(1.+CZETA(N)*S)
      DS(K)=0.25*CZETA(N)*(1.+CETA(N)*T)
      DTX=DTX+DT(K)*X(N)
      DTY=DTY+DT(K)*Y(N)
      DTZ=DTZ+DT(K)*Z(N)
      DSX=DSX+DS(K)*X(N)
      DSY=DSY+DS(K)*Y(N)
      DSZ=DSZ+DS(K)*Z(N)
 1160 D(N)=0.25*       (1.+CETA(N)*T)*(1.+CZETA(N)*S)
      D(10)=1.-T*T
      D(11)=1.-S*S
C
 1200 G1=DTY*DSZ-DSY*DTZ
      G2=DTX*DSZ-DSX*DTZ
      G3=DTX*DSY-DSX*DTY
      G=SQRT(G1*G1+G2*G2+G3*G3)
C
C     EXTERNAL AND INTERNAL PART OF THE TOTAL LOAD VECTOR
C
      DO 1210 K=1,8
 1210 F(3*K)=F(3*K)+G*D(K)
C
      IF (IOP) 2000,2000,1215
C
 1215 DO 1220 K=1,3
 1220 FI(3*K)=FI(3*K)+G*D(8+K)
C
 2000 CONTINUE
C
C
      DO 2010 K=1,8
      F(3*K-2)=F(3*K)*PX
      F(3*K-1)=F(3*K)*PY
 2010 F(3*K)=F(3*K)*PZ
C
      IF (IOP) 3000,3000,2020
C
 2020 DO 2030 K=1,3
      FI(3*K-2)=FI(3*K)*PX
      FI(3*K-1)=FI(3*K)*PY
 2030 FI(3*K)=FI(3*K)*PZ
C
C     STATIC CONDENSATION
C
 2500 DO 2510 I=1,24
      DO 2510 J=1,9
 2510 F(I)=F(I)-A1(I,J)*FI(J)
C
 3000 RETURN
C
      END
      SUBROUTINE HEXA34(F,FI,A1,E,X,Y,Z,EC,NIP,IOP)
C
C***********************************************************************
C
C   FINITE ELEMENT LIBRARY : HEXA34                           SINTEF/NTH
C
C
C     HEXA34 COMPUTES A CONSISTENT VECTOR OF NODAL FORCES DUE TO
C     'TRILINEARLY' VARYING INITIAL STRAINS FOR A (COMPATIBLE OR IN-
C     COMPATIBLE) SOLID HEXAHEDRON ELEMENT WITH 8 NODAL POINTS AND
C     24 DEGREES OF FREEDOM
C
C
C     PROGRAMMED BY : S.L\SETH
C
C     DATE/VERSION : 30.07.74 / 01
C     DATE/VERSION : 88-02-05 / 1.2  REVISION BY K. AAMNES
C     DATE/VERSION : 01-03-17 / 1.3  REVISION BY K. OKSTAD
C     DATE/VERSION : 10-05-21 / 1.4  REVISION BY K. OKSTAD
C
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DOUBLE PRECISION  NC
C
      DIMENSION F(24),FI(9),A1(24,9),E(6,6),X(8),Y(8),Z(8),EC(6,8),W(3),
     +          ABC(3),BJ(3,3),CXI(8),CETA(8),CZETA(8),NC(8),D(3,11),
     +          E0(6),EE0(6)
C
      DATA (CXI(I),I=1,8)  /-1., 1., 1.,-1.,-1., 1., 1.,-1./
      DATA (CETA(I),I=1,8) /-1.,-1., 1., 1.,-1.,-1., 1., 1./
      DATA (CZETA(I),I=1,8)/-1.,-1.,-1.,-1., 1., 1., 1., 1./
C
      IERR=0
C
      IF (NIP-2) 20,25,10
   10 IF (NIP-3) 20,30,20
C
   20 IERR = -1
      GO TO 40
C
C     ABCISSAS AND WEIGHT FACTORS OF THE INTEGRATION FORMULAS
C
   25 ABC(1)=-SQRT(1.0/3.0)
      ABC(2)=-ABS(1)
C
      W(1)=1.0
      W(2)=1.0
      GO TO 40
C
   30 ABC(1)=-SQRT(0.6)
      ABC(2)=0.0
      ABC(3)=-ABC(1)
C
      W(1)=5.0/9.0
      W(2)=8.0/9.0
      W(3)=W(1)
C
   40 DO 50 I=1,24
   50 F(I)=0.
C
      IF (IOP) 65,65,55
C
   55 DO 60 I=1,9
   60 FI(I)=0.
C
   65 IF (IERR) 1000,70,70
C
C     START OF INTEGRATION LOOP
C
   70 DO 500 K=1,NIP
      ZETA=ABC(K)
      WZ=W(K)
      DO 400 J=1,NIP
      ETA=ABC(J)
      WY=W(J)
      DO 300 I=1,NIP
      XI=ABC(I)
      WI=W(I)*WY*WZ
C
C     SHAPE FUNCTIONS ASSOCIATED WITH THE NODES
C
      DO 80 L=1,8
   80 NC(L)=0.125*(1.+CXI(L)*XI)*(1.+CETA(L)*ETA)*(1.+CZETA(L)*ZETA)
C
C     THE INVERSE AND THE DETERMINANT OF THE JACOBIAN MATRIX
C
      CALL JABN30(BJ,X,Y,Z,XI,ETA,ZETA,DETJ,1)
      IF (DETJ) 90,90,100
   90 IERR=-1
      GO TO 40
C
C     THE DERIVATIVES OF THE SHAPE FUNCTIONS WITH RESPECT TO GLOBAL
C     COORDINATES
C
  100 DO 110 N=1,8
      DO 110 M=1,3
  110 D(M,N)=0.125*(BJ(M,1)*CXI(N)*(1.+CETA(N)*ETA)*(1.+CZETA(N)*ZETA)
     +             +BJ(M,2)*CETA(N)*(1.+CXI(N)*XI)*(1.+CZETA(N)*ZETA)
     +             +BJ(M,3)*CZETA(N)*(1.+CXI(N)*XI)*(1.+CETA(N)*ETA))
C
      IF(IOP) 140,140,120
C
  120 DO 130 M=1,3
      D(M,9)=-2.0*BJ(M,1)*XI
      D(M,10)=-2.0*BJ(M,2)*ETA
  130 D(M,11)=-2.0*BJ(M,3)*ZETA
  140 CONTINUE
C
C     INITIAL STRAINS E0 MULTIPLIED BY DETJ (AND WI)
C
      DO 160 M=1,6
      E0(M)=0.
      DO 150 N=1,8
  150 E0(M)=E0(M)+EC(M,N)*NC(N)
  160 E0(M)=DETJ*E0(M)*WI
C
C     THE MATRIX E*E0 (MULTIPLIED BY DETJ)
C
      DO 170 M=1,6
      EE0(M)=0.
      DO 170 N=1,6
  170 EE0(M)=EE0(M)+E(M,N)*E0(N)
C
C     CONSISTENT NODAL FORCES
C
      DO 180 N=1,8
      F(3*N-2)=F(3*N-2)+D(1,N)*EE0(1)+D(2,N)*EE0(4)+D(3,N)*EE0(5)
      F(3*N-1)=F(3*N-1)+D(2,N)*EE0(2)+D(1,N)*EE0(4)+D(3,N)*EE0(6)
  180 F(3*N)  =F(3*N)  +D(3,N)*EE0(3)+D(1,N)*EE0(5)+D(2,N)*EE0(6)
C
      IF (IOP) 300,300,190
C
  190 DO 200 N=1,3
      M=8+N
      FI(3*N-2)=FI(3*N-2)+D(1,M)*EE0(1)+D(2,M)*EE0(4)+D(3,M)*EE0(5)
      FI(3*N-1)=FI(3*N-1)+D(2,M)*EE0(2)+D(1,M)*EE0(4)+D(3,M)*EE0(6)
  200 FI(3*N)  =FI(3*N)  +D(3,M)*EE0(3)+D(1,M)*EE0(5)+D(2,M)*EE0(6)
C
C     END OF INTEGRATION LOOP
C
  300 CONTINUE
  400 CONTINUE
  500 CONTINUE
C
      IF (IOP) 1000,1000,550
C
C     STATIC CONDENSATION
C
  550 DO 600 I=1,24
      DO 600 J=1,9
  600 F(I)=F(I)-A1(I,J)*FI(J)
C
C
 1000 RETURN
C
      END
      SUBROUTINE JABN30(BJ,X,Y,Z,XI,ETA,ZETA,DETJ,IOP)
C
C***********************************************************************
C
C     FINITE ELEMENT LIBRARY : JABN30                         SINTEF/NTH
C
C
C     JABN30 COMPUTES THE INVERSE (BJ) AND/OR THE DETERMINANT (DETJ) OF
C     THE JACOBIAN MATRIX AT AN ARBITRARY POINT WITHIN A HEXAHEDRON
C
C
C     PROGRAMMED BY: S.L\SETH
C
C     DATE/VERSION : 01.06.74 / 01
C     DATE/VERSION : 88-02-05 / 1.2  REVISION BY K. AAMNES
C
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DOUBLE PRECISION J11,J12,J13,J21,J22,J23,J31,J32,J33
C
      DIMENSION BJ(3,3),X(8),Y(8),Z(8)
C
      A1=1.+XI
      A2=1.-XI
      B1=1.+ETA
      B2=1.-ETA
      C1=1.+ZETA
      C2=1.-ZETA
C
      IF (IOP) 75,75,25
C
   25 DO 50 I=1,3
      DO 50 J=1,3
   50 BJ(I,J)=0.
C
C     THE ELEMENTS OF THE JACOBIAN MATRIX
C
   75 J11=0.125*( B1*C1*(X(7)-X(8))-B2*C1*(X(5)-X(6))
     +           +B1*C2*(X(3)-X(4))-B2*C2*(X(1)-X(2)))
      J12=0.125*( B1*C1*(Y(7)-Y(8))-B2*C1*(Y(5)-Y(6))
     +           +B1*C2*(Y(3)-Y(4))-B2*C2*(Y(1)-Y(2)))
      J13=0.125*( B1*C1*(Z(7)-Z(8))-B2*C1*(Z(5)-Z(6))
     +           +B1*C2*(Z(3)-Z(4))-B2*C2*(Z(1)-Z(2)))
      J21=0.125*( A1*C1*(X(7)-X(6))+A2*C1*(X(8)-X(5))
     +           +A1*C2*(X(3)-X(2))+A2*C2*(X(4)-X(1)))
      J22=0.125*( A1*C1*(Y(7)-Y(6))+A2*C1*(Y(8)-Y(5))
     +           +A1*C2*(Y(3)-Y(2))+A2*C2*(Y(4)-Y(1)))
      J23=0.125*( A1*C1*(Z(7)-Z(6))+A2*C1*(Z(8)-Z(5))
     +           +A1*C2*(Z(3)-Z(2))+A2*C2*(Z(4)-Z(1)))
      J31=0.125*( A1*B1*(X(7)-X(3))+A2*B1*(X(8)-X(4))
     +           +A2*B2*(X(5)-X(1))+A1*B2*(X(6)-X(2)))
      J32=0.125*( A1*B1*(Y(7)-Y(3))+A2*B1*(Y(8)-Y(4))
     +           +A2*B2*(Y(5)-Y(1))+A1*B2*(Y(6)-Y(2)))
      J33=0.125*( A1*B1*(Z(7)-Z(3))+A2*B1*(Z(8)-Z(4))
     +           +A2*B2*(Z(5)-Z(1))+A1*B2*(Z(6)-Z(2)))
C
C     THE DETERMINANT
C
      DETJ=J11*(J22*J33-J23*J32)-J12*(J21*J33-J23*J31)
     ++J13*(J21*J32-J22*J31)
C
      IF (DETJ) 2000,2000,100
  100 IF (IOP) 2000,2000,200
C
C     THE INVERSE
C
  200 DET=1./DETJ
      BJ(1,1)=DET*(J22*J33-J32*J23)
      BJ(1,2)=-DET*(J12*J33-J13*J32)
      BJ(1,3)=DET*(J12*J23-J13*J22)
      BJ(2,1)=-DET*(J21*J33-J23*J31)
      BJ(2,2)=DET*(J11*J33-J13*J31)
      BJ(2,3)=-DET*(J11*J23-J13*J21)
      BJ(3,1)=DET*(J21*J32-J22*J31)
      BJ(3,2)=-DET*(J11*J32-J12*J31)
      BJ(3,3)=DET*(J11*J22-J12*J21)
C
C
 2000 RETURN
C
      END
C **********************************************************************
C ********      H   H   EEEEE   X   X    AAA   3333   55555     ********
C ******        H   H   E        X X    A   A     33  5           ******
C *****         HHHHH   EEE       X     AAAAA   333   5555         *****
C ******        H   H   E        X X    A   A     33      5       ******
C ********      H   H   EEEEE   X   X   A   A  3333   5555      ********
C **********************************************************************
C **                                                                   *
C **  This file contains the routines :                                *
C **                                                                   *
C **  - HEXA35                                                         *
C **  - PRIMAS                                                         *
C **  - ISOJAC                                                         *
C **  - LINHEX                                                         *
C **                                                                   *
C **********************************************************************
C
C ######################################################################
C ######################################################  HEXA35     ###
C ######################################################################
C
      SUBROUTINE HEXA35(EM,X,Y,Z,RHO,LPU,IERR)
C
C***********************************************************************
C
C     HEXA35 DETERMINES A CONSISTENT MASS MATRIX  EM  FOR THE HEXAHEDRAL
C     8-NODED SOLID ELEMENT  HEXA.
C
C
C     PROGRAMMED BY :  SJUR LASSESEN
C     DATE/VERSION  :  020488 / 1.0
C
C***********************************************************************
C
      DOUBLE PRECISION  EM(24,24),X(8),Y(8),Z(8),RHO,GC,COORDS(3,8)
      INTEGER           LPU,IERR,I
      DOUBLE PRECISION  WORK(216)
      EXTERNAL          LINHEX
C
C -------------  Compute consistent mass matrix using STYBER - routines
C -------------  GC - CONVERSIONFACTOR UNITS: F = M*A/GC
C -------------  IERR   - ERRORCONDITION; = 0 OK         (OUTPUT)
C -------------                           < 0 UNABLE TO COMPUTE EM
C -------------  WORK   - WORKVECTOR; DIMENSION AT LEAST 9*NDOF
C
      GC = -1.0D0
      DO 10 I = 1, 8
         COORDS(1,I) = X(I)
         COORDS(2,I) = Y(I)
         COORDS(3,I) = Z(I)
   10 CONTINUE
C
      CALL PRIMAS(LPU,RHO,GC,24,EM,LINHEX,8,COORDS,3,3,3,3,WORK,IERR)
      IF (IERR.LT.0) WRITE(LPU,6900)
C
      RETURN
 6900 FORMAT(///4X,'***ERROR FROM A PREPRO ROUTINE'
     +         /4X,'   ROUTINE NAME : HEXA35'
     +         /4X,'   UNABLE TO COMPUTE EM')
C
      END
C ######################################################################
C ######################################################  PRIMAS     ###
C ######################################################################
      SUBROUTINE PRIMAS(IOUT,RO,GC,NDOF,CMA,ZSHAPE,NPOI,X,NRX,
     .                  NPX,NPY,NPZ,WORK,IERR)
C
C ... PURPOSE: ASSEMBLE MASSMATRIX FOR ANY ISOPARAMETRIC ELEMENT WITH
C              PRISMATIC BASIS SHAPE
C
C ... METHOD: MAPPING OF NODAL COORDINATES INTO A NORMALIZED CUBE WITH
C             CENTER IN (0.,0.,0.) AND EDGES = 2.
C             SHAPEFUNCTION SUBROUTINE MUST BE PROVIDED BY THE USER
C
C ... PARAMETERS: IOUT   - LOGICAL UNIT# OUTPUTFILE              (INPUT)
C                 RO     - DENSITY OF MATERIAL IN ELEMENT        (INPUT)
C                 GC     - CONVERSIONFACTOR UNITS: F = M*A/GC    (INPUT)
C                 NDOF   - DEGREES OF FREEDOM OF ELEMENT         (INPUT)
C                 CMA    - CONSISTENT MASS MATRIX               (OUTPUT)
C                 ZSHAPE - SHAPEFUNCTIONGENERATING ROUTINE       (INPUT)
C                 NPOI   - # OF NODES OF ELEMENT                 (INPUT)
C                 X      - CARTESIAN COORDS OF NODES             (INPUT)
C                 NRX    - COLUMNSIZE OF X                       (INPUT)
C                 NPX    - # OF INTEGRATIONPOINTS IN X-DIRECTION (INPUT)
C                 NPY    - # OF INTEGRATIONPOINTS IN Y-DIRECTION (INPUT)
C                 NPZ    - # OF INTEGRATIONPOINTS IN Z-DIRECTION (INPUT)
C                 WORK   - WORKVECTOR; DIMENSION AT LEAST 4*NDOF+9
C                 IERR   - ERRORCONDITION; = 0 OK               (OUTPUT)
C                                          < 0 UNABLE TO COMPUTE CMA
C
C ... CODED BY: H STENBERG MARCH 1980
C
C
      INTEGER           IOUT,NDOF,NPOI,NRX,NPX,NPY,NPZ,IERR
      DOUBLE PRECISION  RO,GC,CMA(NDOF,NDOF),X(NRX,NPOI),WORK(*)
      EXTERNAL          ZSHAPE
C
      INTEGER           NSH,NN,NJAC,J,J1,J2,JNX,K,IX,IY,IZ
      DOUBLE PRECISION  XSI,ETA,ZETA,XW,YW,ZW,C, DET3
C
      IERR = 0
C
C ... WORK MUST PROVIDE SPACE FOR A VECTOR OF LENGTH 3*NPOI,
C     A MATRIX OF LENGTH 3*NDOF, AND THE JACOBI-MATRIX (3*3)
C
      NSH  = 1
      NN   = NSH + NPOI*3
      NJAC = NN + NDOF*3
C
C ... CLEAR MASS MATRIX
C
      CALL DCOPY (NDOF*NDOF,0.0D0,0,CMA,1)
      CALL DCOPY (3*NDOF,0.0D0,0,WORK(NN),1)
C
C ... ESTABLISH MASS MATRIX
C
      DO 30 IZ = 1,NPZ
      CALL GAQDTA(IZ,ZETA,ZW,NPZ)
C
      DO 30 IY = 1,NPY
      CALL GAQDTA(IY,ETA,YW,NPY)
C
      DO 30 IX = 1,NPX
      CALL GAQDTA(IX,XSI,XW,NPX)
C
      CALL ISOJAC(XSI,ETA,ZETA,ZSHAPE,3,WORK(NJAC),3,NPOI,
     .            WORK(NSH),3,0,X,NRX,IERR)
      IF (IERR .LT. 0) GO TO 800
C
      C = XW*YW*ZW*RO*DET3(WORK(NJAC))/GC
C
C ... SHAPE FUNCTION MATRIX (N) (3*NDOF)
C
      DO 10 K = 1,3
      DO 10 J = 1,NPOI
      J1      = K + (J-1)*3
      JNX     = (J1-1)*3 + K - 1
      J2      = (J-1)*3 + 1
   10 WORK(NN+JNX) = WORK(J2)
C
C ... INTEGRATION OF N-T * N FOR CURRENT POINT
C
      CALL DGEMM('T','N',NDOF,NDOF,3,C,WORK(NN),3,WORK(NN),3,
     .           1.0D0,CMA,NDOF)
C
   30 CONTINUE
C
      GO TO 999
C
C ... ERROR DETECTED
C
  800 WRITE(IOUT,9010) -IERR
  999 RETURN
 9010 FORMAT('   .. PRIMAS .. COORDINATE NO',I2,' OUTSIDE BOUNDS')
      END
      DOUBLE PRECISION FUNCTION DET3(J)
      DOUBLE PRECISION J(3,3),F1,F2,F3
      F1=J(2,2)*J(3,3)-J(2,3)*J(3,2)
      F2=J(2,3)*J(3,1)-J(2,1)*J(3,3)
      F3=J(2,1)*J(3,2)-J(2,2)*J(3,1)
      DET3=J(1,1)*F1+J(1,2)*F2+J(1,3)*F3
      RETURN
      END
C ######################################################################
C ######################################################  ISOJAC     ###
C ######################################################################
      SUBROUTINE ISOJAC(X,Y,Z,ZSHAPE,NAJ,AJAC,NRA,NPOI,SH,NRSH,IDER,
     .                  XP,NRXP,ICOND)
C
C ... PURPOSE: COMPUTE JACOBIAN MATRIX FOR AN ISOPARAMETRIC ELEMENT
C
C ... METHOD: THE VALUES OF EACH OF THE PARTIAL DERIVATIVES WITH RESPECT
C             TO NORMALIZED COORDINATES ARE COMPUTED IN A SUBROUTINE
C             WITH FORMAL NAME ZSHAPE, THIS ROUTINE IS PROVIDED BY THE
C             USER FOR HIS PARTICULAR FUNCTION
C
C
C ... PARAMETERS: X,Y,Z  - CURVILINIAR COORDINATES OF CURRENT POSITION
C                 ZSHAPE - ROUTINE THAT COMPUTES VALUES ASSOCIATED WITH
C                          THE SHAPEFUNCTION
C                 NAJ    - # OF ROWS AND COLUMNS IN JACOBI MATRIX
C                 AJC    - JACOBI MATRIX                        (OUTPUT)
C                 NRA    - COLUMNSIZE OF AJAC
C                 NPOI   - # OF POINTS WITH ELEMENT
C                 SH     - WORKMATRIX, WILL CONTAIN SHAPEFUNCTIONS AS
C                          ENTRIES OF 1ST ROW (IDER=0), OR PARTIAL
C                          DERIVATIVES AS COLUMNS (IDER=1)      (OUTPUT)
C                 NRSH   - COLUMNSIZE OF SH  (>=3)
C                 IDER   - ORDER OF DERIVATIVES OF SHAPEFUNCTIONS
C                 XP     - MATRIX; COLUMNS CONTAIN CARTESIAN COORDINATES
C                          OF NODES                             (INPUT)
C                 NRXP   - COLUMNSIZE OF XP
C                 ICOND  - < 0 : COMPUTATION FAILED. = 0 COMPUTATION OK
C
C ... CODED BY: H STENBERG MARCH - 1980
C
C
      INTEGER            NAJ,NRA,NPOI,NRSH,IDER,NRXP,ICOND, J
      DOUBLE PRECISION   X,Y,Z,AJAC(NRA,NAJ),SH(NRSH,NPOI),XP(NRXP,NPOI)
      EXTERNAL           ZSHAPE
C
C ... VALID NORMALIZED COORDS?
C
      ICOND = 0
      IF (ABS(X) .GT. 1.0D0) ICOND = -1
      IF (ABS(Y) .GT. 1.0D0) ICOND = -2
      IF (ABS(Z) .GT. 1.0D0) ICOND = -3
      IF (ICOND .LT. 0) GO TO 999
C
C ... CREATE MATRIX WHOSE COLUMNS ARE PARTIAL DERIVATIVES OF
C     EACH SHAPEFUNCTION WITH RESPECT TO XSI, ETA AND ZETA
C
      DO 10 J = 1,NPOI
   10 CALL ZSHAPE(J,X,Y,Z,SH(1,J),1)
C
C ... COMPUTE JACOBIAN MATRIX: AJAC = SH * XP^T
C
      CALL DGEMM ('N','T',NAJ,NAJ,NPOI,1.0D0,SH,NRSH,XP,NRXP,
     .            0.0D0,AJAC,NRA)
C
      IF (IDER .GT. 0) GO TO 999
C
      DO 20 J = 1,NPOI
   20 CALL ZSHAPE (J,X,Y,Z,SH(1,J),0)
C
  999 RETURN
      END

C ######################################################################
C ######################################################  LINHEX     ###
C ######################################################################
      SUBROUTINE LINHEX(NGRID,XSI,ETA,ZETA,W,IDER)
C
C ... PURPOSE: COMPUTE VALUES OF SHAPEFUNCTIONS AND FIRST ORDER PARTIAL
C              DERIVATIVES FOR A LINEAR HEXAHEDRAL ELEMENT IN NORMALIZED
C              COORDINATES.
C
C ... METHOD: SERENDIPITY FAMILY OF PRISMATIC ELEMENTS.
C
C ... PARAMETERS: NGRID - GRIDPOINT NUMBER (SEE TABLES BELOW) (INPUT)
C                 XSI,ETA,ZETA - COORDINATES OF CURRENT POINT (INPUT)
C                 W     - COMPUTED VALUES AS REQUESTED       (OUTPUT)
C                 IDER  - ORDER OF DERIVATIVE I E             (INPUT)
C                         = 0 FOR FUNCTIONVALUE, = 1 FOR 1ST ORDER DERIV
C
C ... CODED BY: H STENBERG MARCH - 80
C
C
      INTEGER            NGRID,IDER
      DOUBLE PRECISION   XSI,ETA,ZETA,W(*)
      DOUBLE PRECISION   A,B,C,X(8),Y(8),Z(8)
C
      DATA  X /  1. , -1. , -1. ,  1. ,  1. , -1. , -1. ,  1. /
      DATA  Y /  1. ,  1. , -1. , -1. ,  1. ,  1. , -1. , -1. /
      DATA  Z /  1. ,  1. ,  1. ,  1. , -1. , -1. , -1. , -1. /
C
      IF (NGRID .LT. 1 .OR. NGRID .GT. 8) GO TO 999
C
      A    = 1.0D0 + XSI*X(NGRID)
      B    = 1.0D0 + ETA*Y(NGRID)
      C    = 1.0D0 + ZETA*Z(NGRID)
C
      IF (IDER .GT. 0) GO TO 10
C
C .. VALUE OF SHAPEFUNCTION AT CURRENT POINT
C
      W(1) = 0.125D0*A*B*C
      GO TO 999
C
C ... VALUE OF FIRST DERIVATIVES AT CURRENT POINT
C
   10 W(1) = 0.125D0*X(NGRID)*B*C
      W(2) = 0.125D0*Y(NGRID)*A*C
      W(3) = 0.125D0*Z(NGRID)*A*B
C
  999 RETURN
      END
      SUBROUTINE GAQDTA(I,A,W,NPOINT)
C
C ... PURPOSE: FIND ABSCISSA AND WEIGHT COEFFICIENT OF THE GAUSSIAN
C              QUADRATURE FORMULA
C
C ... METHOD:  THE ACTUAL VALUES ARE LISTED IN ATABLE FOR A DEFINITE
C              INTEGRAL IN <-1., 1.> , THE REQUESTED VALUES WILL BE
C              FETCHED FROM THESE TABLES
C
C ... PARAMETERS: I      - POINT NUMBER          (INPUT)
C                 A      - ABSCISSA              (OUTPUT)
C                 W      - WEGHT                 (OUTPUT)
C                 NPOINT - NO OF QUADRATUREPOINT (INPUT)
C
C ... REMARK: ORDER OF ACCURACY: (NPOINT*2 - 1).  1 < NPOINT < 11
C
C ... REF: ABRAMOWITZ AND STEGUN: HANDBOOK OF MATHEMATICAL FUNCTIONS
C                                 DOVER
C
C
      INTEGER          I,J,IP,INX,NP,NPOINT,ITP(10)
      DOUBLE PRECISION A,W,ATAB(29),WTAB(29)
C
C ... # OF ENTRIES OF TABLES PRIOR TO NPOINT QUADRATURE
C
      DATA  ITP/ 0 , 0 , 1 , 3 , 5 , 8 , 11 , 15 , 19 , 24 /
C
C ... ABSCISSAE
C
      DATA ATAB/
C                  NPOINT = 2
     .  0.57735026918963  ,
C                  NPOINT = 3
     .  0.77459666924148  ,  0.                ,
C                  NPOINT = 4
     .  0.86113631159405  ,  0.33998104358486  ,
C                  NPOINT = 5
     .  0.90617984593866  ,  0.53846931010568  ,  0.                ,
C                  NPOINT = 6
     .  0.93246951420315  ,  0.66120938646627  ,  0.2386191860832   ,
C                  NPOINT = 7
     .  0.94910791234276  ,  0.74153118559939  ,  0.4058451513774   , 0.
C                  NPOINT = 8
     ., 0.96028985649754  ,  0.79666647741363  ,  0.52553240991633  ,
     .  0.18343464249565  ,
C                  NPOINT = 9
     .  0.96816023950763  ,  0.83603110732664  ,  0.61337143270059  ,
     .  0.32425342340381  ,  0.                ,
C                  NPOINT = 10
     .  0.97390652851717  ,  0.86506336668899  ,  0.67940956829902  ,
     .  0.43339539412925  ,  0.14887433898163  /
C
C ... WEIGHTS
C
      DATA WTAB/
C                  NPOINT = 2
     .  1.                ,
C                  NPOINT = 3
     .  0.55555555555556  ,  0.88888888888889  ,
C                  NPOINT = 4
     .  0.34785484513745  ,  0.65214515486255  ,
C                  NPOINT = 5
     .  0.23692688505619  ,  0.47862867049937  ,  0.56888888888889  ,
C                  NPOINT = 6
     .  0.17132449237917  ,  0.36076157304814  ,  0.46791393457269  ,
C                  NPOINT = 7
     .  0.129484966168870 ,  0.27970539148928  ,  0.38183005050512  ,
     .  0.41795918367347  ,
C                  NPOINT = 8
     .  0.10122853629038  ,  0.22238103445337  ,  0.31370664587789  ,
     .  0.36268378337836  ,
C                  NPOINT = 9
     .  0.081274388361574 ,  0.18064816069486  ,  0.26061069640294  ,
     .  0.31234707704000  ,  0.33023935500126  ,
C                  NPOINT = 10
     .  0.066671344308688 ,  0.14945134915058  ,  0.21908636251598  ,
     .  0.26926671931     ,  0.29552422471475  /
C
C ... SET INITIALVALUES FOR LOCAL VARIABLES
C
      J       = I
      IP      = 1
      NP      = MAX0(NPOINT,2)
      NP      = MIN0(NP,10)
      IF (I .GT. NP) GO TO 999
C
C ... POINT TO THE RIGHT OF ORIGIN?
C
      IF  ((2*I).GT.NP)   J       = NP+1 - I
      IF  ((2*I).GT.NP)   IP      = 2
C
      INX     = ITP(NP) + J
C
C ... SET VALUES
C
      W       = WTAB(INX)
      A       = ATAB(INX)*(-1.)**IP
C
  999 RETURN
      END
