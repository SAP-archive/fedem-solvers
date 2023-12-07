C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
C     This file is part of FEDEM - https://openfedem.org
C
      SUBROUTINE  TMRF31
     +       (X,Y,DM,ALPHA,BETA,IAT,LS,HH,SM,M,STATUS)
C
C             T Y P E   &   D I M E N S I O N
C
      CHARACTER*(*)      STATUS
      INTEGER            IAT,M,LS(12)
      DOUBLE PRECISION   X(3),Y(3),DM(3,3),ALPHA,BETA,SM(M,M)
      DOUBLE PRECISION   XT(3),YT(3),DMT(3,3),F,HH(3,9)
      INTEGER            I,IAS,J,N,T
      INTEGER            NTRIGS(0:3),TNODES(3,4,0:3),LST(9)
      DATA               NTRIGS /1,2,2,4/
      DATA               TNODES /1,2,3,9*0,
     +                           1,2,3, 3,4,1, 6*0,
     +                           1,2,4, 2,3,4, 6*0,
     +                           1,2,3, 3,4,1, 1,2,4, 2,3,4/
C
C             L O G I C
C
      STATUS= ' '
C
      DO 2200 I=1,3
         DO 2200 J=1,3
            DMT(I,J)=DM(I,J)
 2200 CONTINUE
C
      IAS=MAX(0,MIN(IAT,3))
      F=1.0
      IF (IAS.EQ.3) F=0.5
      DO 3000 T=1,NTRIGS(IAS)
         DO 2500 I=1,3
            N=TNODES(I,T,IAS)
            XT(I)=X(N)
            YT(I)=Y(N)
            LST(2*I-1)=LS(3*N-2)
            LST(2*I  )=LS(3*N-1)
            LST(  I+6)=LS(3*N  )
 2500    CONTINUE
C
         CALL SM3MB (XT,YT,DMT,ALPHA,F,LST,SM,M,STATUS)
         IF (BETA.GT.0.0) THEN
            CALL SM3MH (XT,YT,DMT,F*BETA,LST,HH,SM,M,STATUS)
         ENDIF
         IF (STATUS(1:1).NE.' ') RETURN
 3000 CONTINUE
C
      RETURN
C
      END
      SUBROUTINE   SM3MB
     +      (X,Y,DM,ALPHA,F,LS,SM,M,STATUS)
C
C        T Y P E   &   D I M E N S I O N
C
      CHARACTER*(*)     STATUS
      INTEGER           M,LS(9)
      DOUBLE PRECISION  X(3),Y(3),DM(3,3),ALPHA,F,P(9,3),SM(M,M)
      DOUBLE PRECISION  AREA2,C
      DOUBLE PRECISION  D11,D12,D13,D22,D23,D33
      DOUBLE PRECISION  X21,X32,X13,Y21,Y32,Y13
      DOUBLE PRECISION  X12,X23,X31,Y12,Y23,Y31
      DOUBLE PRECISION  S1,S2,S3
      INTEGER           I,J,K,L,N
C
C        L O G I C
C
      STATUS = ' '
      X21 =  X(2)-X(1)
      X12 = -X21
      X32 =  X(3)-X(2)
      X23 = -X32
      X13 =  X(1)-X(3)
      X31 = -X13
      Y21 =  Y(2)-Y(1)
      Y12 = -Y21
      Y32 =  Y(3)-Y(2)
      Y23 = -Y32
      Y13 =  Y(1)-Y(3)
      Y31 = -Y13
C
      AREA2 = Y21*X13 - X21*Y13
      IF (AREA2.LT.-1.0D-16) THEN
         STATUS = 'NEGA_AREA'
         RETURN
      ELSE IF (AREA2.LE.1.0D-16) THEN
         STATUS = 'ZERO_AREA'
         RETURN
      ENDIF
C
      P(1,1) = Y23
      P(2,1) = 0.0
      P(3,1) = Y31
      P(4,1) = 0.0
      P(5,1) = Y12
      P(6,1) = 0.0
      P(1,2) = 0.0
      P(2,2) = X32
      P(3,2) = 0.0
      P(4,2) = X13
      P(5,2) = 0.0
      P(6,2) = X21
      P(1,3) = X32
      P(2,3) = Y23
      P(3,3) = X13
      P(4,3) = Y31
      P(5,3) = X21
      P(6,3) = Y12
      N      = 6
C
      IF (ALPHA.GT.0.0) THEN
         P(7,1) = Y23*(Y13-Y21)*ALPHA/6
         P(7,2) = X32*(X31-X12)*ALPHA/6
         P(7,3) = (X31*Y13-X12*Y21)*ALPHA/3
         P(8,1) = Y31*(Y21-Y32)*ALPHA/6
         P(8,2) = X13*(X12-X23)*ALPHA/6
         P(8,3) = (X12*Y21-X23*Y32)*ALPHA/3
         P(9,1) = Y12*(Y32-Y13)*ALPHA/6
         P(9,2) = X21*(X23-X31)*ALPHA/6
         P(9,3) = (X23*Y32-X31*Y13)*ALPHA/3
         N      = 9
      ENDIF
C
      C   = 0.5*F/AREA2
      D11 = C * DM(1,1)
      D22 = C * DM(2,2)
      D33 = C * DM(3,3)
      D12 = C * DM(1,2)
      D13 = C * DM(1,3)
      D23 = C * DM(2,3)
C
      DO 3000 J=1,N
         L  = LS(J)
         S1 = D11*P(J,1) + D12*P(J,2) + D13*P(J,3)
         S2 = D12*P(J,1) + D22*P(J,2) + D23*P(J,3)
         S3 = D13*P(J,1) + D23*P(J,2) + D33*P(J,3)
         DO 2500 I=1,J
            K       = LS(I)
            SM(K,L) = SM(K,L) + (S1*P(I,1)+S2*P(I,2)+S3*P(I,3))
            SM(L,K) = SM(K,L)
 2500    CONTINUE
 3000 CONTINUE
C
      RETURN
      END
      SUBROUTINE  SM3MH
     +      (X,Y,DM,F,LS,HH,SM,M,STATUS)
C
C        T Y P E   &   D I M E N S I O N
C
      CHARACTER*(*)      STATUS
      INTEGER            LS(9),M
      DOUBLE PRECISION   X(3),Y(3),DM(3,3),F,SM(M,M)
      DOUBLE PRECISION   XC(3),YC(3),XM(3),YM(3)
      DOUBLE PRECISION   GT(9,9),HH(3,9),SQH(3,3),T(9)
      DOUBLE PRECISION   QX(3,3),QY(3,3)
      DOUBLE PRECISION   AREA,AREA2
      DOUBLE PRECISION   A1J,A2J,A3J, B1J,B2J,B3J
      DOUBLE PRECISION   C,CJ,SJ,DL,DX,DY
      DOUBLE PRECISION   D11,D12,D13,D22,D23,D33,JXX,JXY,JYY
      DOUBLE PRECISION   S1,S2,S3,S4,S5,S6
      DOUBLE PRECISION   X0,Y0,XI,YI
      INTEGER            I,J,K,L,ISING,IPERM(9)
C
C        L O G I C
C
      STATUS = ' '
      AREA2 = (Y(2)-Y(1))*(X(1)-X(3)) - (X(2)-X(1))*(Y(1)-Y(3))
C
      IF (AREA2.LT.-1.0D-16) THEN
         STATUS = 'NEGA_AREA'
         RETURN
      ELSE IF (AREA2.LE.1.0D-16) THEN
         STATUS = 'ZERO_AREA'
         RETURN
      ENDIF
C
      X0 = (X(1)+X(2)+X(3))/3.0
      Y0 = (Y(1)+Y(2)+Y(3))/3.0
      AREA = 0.5*AREA2
      C = 1./SQRT(AREA)
C
      XC(1) = C * (X(1)-X0)
      XC(2) = C * (X(2)-X0)
      XC(3) = C * (X(3)-X0)
      YC(1) = C * (Y(1)-Y0)
      YC(2) = C * (Y(2)-Y0)
      YC(3) = C * (Y(3)-Y0)
      XM(1) = 0.5 * (XC(2)+XC(3))
      XM(2) = 0.5 * (XC(3)+XC(1))
      XM(3) = 0.5 * (XC(1)+XC(2))
      YM(1) = 0.5 * (YC(2)+YC(3))
      YM(2) = 0.5 * (YC(3)+YC(1))
      YM(3) = 0.5 * (YC(1)+YC(2))
C
C        Form G' (G transposed) in GT and initialize H in HH
C                                                     h
      DO 1300 I=1,9
         DO 1200 J=1,6
            GT(J,I) = 0.
 1200    CONTINUE
         HH(1,I) = 0.
         HH(2,I) = 0.
         HH(3,I) = 0.
 1300 CONTINUE
C
      D11 = F * DM(1,1)
      D22 = F * DM(2,2)
      D33 = F * DM(3,3)
      D12 = F * DM(1,2)
      D13 = F * DM(1,3)
      D23 = F * DM(2,3)
C
      JXX = -2. * (XC(1)*XC(2) + XC(2)*XC(3) + XC(3)*XC(1))/3.
      JXY =       (XC(1)*YC(1) + XC(2)*YC(2) + XC(3)*YC(3))/3.
      JYY = -2. * (YC(1)*YC(2) + YC(2)*YC(3) + YC(3)*YC(1))/3.
C
      DO 2500 J=1,3
         DX = XM(J)-XC(J)
         DY = YM(J)-YC(J)
         DL = SQRT (DX**2 + DY**2)
         CJ = DX/DL
         SJ = DY/DL
         A1J = -0.5*SJ*CJ**2
         A2J =  0.5*CJ**3
         B2J = -0.5*SJ**3
         B3J =  0.5*SJ**2*CJ
         A3J = -(B2J + A1J + A1J)
         B1J = -(B3J + B3J + A2J)
C
         GT(1,2*J-1) =  1.
         GT(2,2*J  ) =  1.
         GT(3,2*J-1) = -YC(J)
         GT(3,2*J  ) =  XC(J)
         GT(3,  J+6) =  C
         GT(4,2*J-1) =  XC(J)
         GT(6,2*J-1) =  YC(J)
         GT(5,2*J  ) =  YC(J)
         GT(6,2*J  ) =  XC(J)
C
         HH(J,J+6) = 1.
C
         QX(J,1) =  A1J
         QX(J,2) =  B2J
         QX(J,3) = -2.*B3J
         QY(J,1) =  A2J
         QY(J,2) =  B3J
         QY(J,3) = -2.*A1J
C
         S1 = D11*QX(J,1) + D12*QX(J,2) + D13*QX(J,3)
         S2 = D12*QX(J,1) + D22*QX(J,2) + D23*QX(J,3)
         S3 = D13*QX(J,1) + D23*QX(J,2) + D33*QX(J,3)
         S4 = D11*QY(J,1) + D12*QY(J,2) + D13*QY(J,3)
         S5 = D12*QY(J,1) + D22*QY(J,2) + D23*QY(J,3)
         S6 = D13*QY(J,1) + D23*QY(J,2) + D33*QY(J,3)
C
         DO 2200 I=1,3
            XI = XC(I)
            YI = YC(I)
            GT(J+6,2*I-1) =  A1J*XI*XI + 2.*A2J*XI*YI + A3J*YI*YI
            GT(J+6,2*I  ) =  B1J*XI*XI + 2.*B2J*XI*YI + B3J*YI*YI
            GT(J+6,I+6  ) = -C*(CJ*XI+SJ*YI)
 2200    CONTINUE
C
         DO 2400 I=1,J
            SQH(I,J) = JXX * (QX(I,1)*S1+QX(I,2)*S2+QX(I,3)*S3)
     +               + JXY * (QX(I,1)*S4+QX(I,2)*S5+QX(I,3)*S6
     +                       +QY(I,1)*S1+QY(I,2)*S2+QY(I,3)*S3)
     +               + JYY * (QY(I,1)*S4+QY(I,2)*S5+QY(I,3)*S6)
 2400    CONTINUE
 2500 CONTINUE
C
C        Factor G' and backsolve to obtain H
C                                           h
      CALL LUFACT (GT,9,9,IPERM,T,ISING)
      IF (ISING.NE.0) THEN
         STATUS = 'SINGULAR_G'
         RETURN
      ENDIF
C
      CALL LUSOLV (GT,9,9,IPERM,HH,3,-3)
C
C        Form physical stiffness and add to incoming SM
C
      DO 4000 J=1,9
         L = LS(J)
         S1 = SQH(1,1)*HH(1,J) + SQH(1,2)*HH(2,J) + SQH(1,3)*HH(3,J)
         S2 = SQH(1,2)*HH(1,J) + SQH(2,2)*HH(2,J) + SQH(2,3)*HH(3,J)
         S3 = SQH(1,3)*HH(1,J) + SQH(2,3)*HH(2,J) + SQH(3,3)*HH(3,J)
         DO 3500 I=1,J
            K = LS(I)
            SM(K,L) = SM(K,L) + (S1*HH(1,I)+S2*HH(2,I)+S3*HH(3,I))
            SM(L,K) = SM(K,L)
 3500    CONTINUE
 4000 CONTINUE
C
      RETURN
      END
      SUBROUTINE  LUFACT (A,NA,NN,IPERM,V,ISING)
C
C        T Y P E   A N D   D I M E N S I O N
C
      INTEGER            NA,NN,IPERM(*),ISING
      DOUBLE PRECISION   A(NA,*),V(*),MACTOL,X,Y
      INTEGER            I,J,K,L,N
      PARAMETER        ( MACTOL = 2.0D-16 )
C
C        L O G I C
C
C        Initialization
C
      ISING  = 0
      N      = IABS(NN)
      IF (N.LE.0) RETURN
C
C        Calculate the Euclidean norm of each row of A
C
      DO 1500 I=1,N
         Y = 0.0
         DO 1200 J=1,N
            Y = Y + A(I,J)*A(I,J)
 1200    CONTINUE
         IF (Y.GT.0.0) THEN
            V(I) = SQRT(1.0/Y)
         ELSE
            V(I) = 0.0
         ENDIF
 1500 CONTINUE
C
C        Main loop
C
      DO 4000 K=1,N
         IPERM(K) = K
         IF (V(K).LE.0.0) GOTO 4000
         L = K
         X = 0.0
         DO 2000 I=K,N
            Y = 0.0
            DO 1600 J=1,K-1
               Y = Y + A(I,J) * A(J,K)
 1600       CONTINUE
            A(I,K) = A(I,K) - Y
            Y = ABS(V(I)*A(I,K))
            IF (Y.GT.X) THEN
               X = Y
               L = I
            ENDIF
 2000    CONTINUE
C
C        Interchange rows
C
         IF (L.NE.K) THEN
            DO 2100 J=1,N
               Y      = A(K,J)
               A(K,J) = A(L,J)
               A(L,J) = Y
 2100       CONTINUE
            V(L) = V(K)
            IPERM(K) = L
         ENDIF
C
C        Test for singularity
C
         IF (X.LE.MACTOL) THEN
            ISING = N-(K-1)
            RETURN
         ENDIF
         X = 1.0/A(K,K)
         A(K,K) = X
C
C        Calculate elements of strict upper triangle
C
         DO 3600 J=K+1,N
            Y = 0.0
            DO 3400 I=1,K-1
               Y = Y + A(K,I) * A(I,J)
 3400       CONTINUE
            A(K,J) = (A(K,J)-Y) * X
 3600    CONTINUE
 4000 CONTINUE
C
      RETURN
      END
      SUBROUTINE  LUSOLV (A,NA,NN,IPERM,B,NB,M)
C
C        T Y P E   A N D   D I M E N S I O N
C
      INTEGER            NA,NN,IPERM(*),NB,M
      DOUBLE PRECISION   A(NA,*),B(NB,*),T,SUM
      LOGICAL            TRANSP
      INTEGER            I,J,K,N,NRHS
C
C        L O G I C
C
C        Initialization
C
      NRHS = IABS(M)
      TRANSP = .FALSE.
      IF (M.LT.0) TRANSP = .TRUE.
         N = IABS(NN)
C
C        Intercharge elements of rhs matrix b
C
      DO 1800 I=1,N
         K = IPERM(I)
         IF (K.NE.I) THEN
            IF (.NOT.TRANSP) THEN
               DO 1400 J=1,NRHS
                  T      = B(I,J)
                  B(I,J) = B(K,J)
                  B(K,J) = T
 1400          CONTINUE
            ELSE
               DO 1600 J=1,NRHS
                  T      = B(J,I)
                  B(J,I) = B(J,K)
                  B(J,K) = T
 1600          CONTINUE
            ENDIF
         ENDIF
 1800 CONTINUE
C
C        Cycle over right hand sides
C
      DO 4000 J=1,NRHS
         IF(.NOT.TRANSP) THEN
C
C        Solve for rhs: column storage
C
            B(1,J) = B(1,J) * A(1,1)
            DO 2400 I=1,N-1
               SUM = 0.0
               DO 2300 K=1,I
                  SUM = SUM - A(I+1,K) * B(K,J)
 2300          CONTINUE
               B(I+1,J) = (B(I+1,J)+SUM) * A(I+1,I+1)
 2400       CONTINUE
            DO 2800 I=N-1,1,-1
               SUM = 0.0
               DO 2600 K=1,N-I
                  SUM = SUM - A(I,I+K) * B(I+K,J)
 2600          CONTINUE
               B(I,J) = B(I,J) + SUM
 2800       CONTINUE
         ELSE
C
C        Solve for rhs: row storage
C
            B(J,1) = B(J,1) * A(1,1)
            DO 3400 I=1,N-1
               SUM = 0.0
               DO 3200 K=1,I
                  SUM = SUM - A(I+1,K) * B(J,K)
 3200          CONTINUE
               B(J,I+1) = (B(J,I+1)+SUM) * A(I+1,I+1)
 3400       CONTINUE
            DO 3800 I=N-1,1,-1
               SUM = 0.0
               DO 3600 K=1,N-I
                  SUM = SUM - A(I,I+K) * B(J,I+K)
 3600          CONTINUE
               B(J,I) =B(J,I) + SUM
 3800       CONTINUE
         ENDIF
 4000 CONTINUE
C
      RETURN
      END
      SUBROUTINE TMRF32(HH,X,Y,DM,ALPHA,SM,STATUS)
C
C ******************************************************************
C
C     TMRF32 GENERATES A STRESS MATRIX  SM  FOR THE  TMRF  ELEMENT
C
C
C     PROGRAMMED BY :  KETIL AAMNES
C     DATE/VERSION  :  101085 / 1.0
C
C ******************************************************************
C
      INTEGER            I,J,K
      DOUBLE PRECISION   X(3),Y(3),X21,X12,X32,X23,X13,X31
      DOUBLE PRECISION   Y21,Y12,Y32,Y23,Y13,Y31,AREA,X0,Y0
      DOUBLE PRECISION   C,XC(3),YC(3),XM(3),YM(3),DX,DY,DL
      DOUBLE PRECISION   CJ,SJ,A1J,A2J,A3J,B1J,B2J,B3J,S,ALPHA
      DOUBLE PRECISION   L(9,3),TMB(3,9),TMH(3,9),HH(3,9),BH(3,3)
      DOUBLE PRECISION   TM(3,9),SM(3,9),DM(3,3)
      CHARACTER*(*)      STATUS
C
C ---------------------------- Generates the lumping matrix --------
C
      STATUS = ' '
      X21 =  X(2)-X(1)
      X12 = -X21
      X32 =  X(3)-X(2)
      X23 = -X32
      X13 =  X(1)-X(3)
      X31 = -X13
      Y21 =  Y(2)-Y(1)
      Y12 = -Y21
      Y32 =  Y(3)-Y(2)
      Y23 = -Y32
      Y13 =  Y(1)-Y(3)
      Y31 = -Y13
C
      AREA = 0.5 * (Y21*X13 - X21*Y13)
      IF (AREA.LT.-1.0D-16) THEN
         STATUS = 'NEGA_AREA'
         RETURN
      ELSE IF (AREA.LE.1.0D-16) THEN
         STATUS = 'ZERO_AREA'
         RETURN
      ENDIF
C
      L(1,1) = 0.5 * Y23
      L(2,1) = 0.5 * 0.0
      L(4,1) = 0.5 * Y31
      L(5,1) = 0.5 * 0.0
      L(7,1) = 0.5 * Y12
      L(8,1) = 0.5 * 0.0
      L(1,2) = 0.5 * 0.0
      L(2,2) = 0.5 * X32
      L(4,2) = 0.5 * 0.0
      L(5,2) = 0.5 * X13
      L(7,2) = 0.5 * 0.0
      L(8,2) = 0.5 * X21
      L(1,3) = 0.5 * X32
      L(2,3) = 0.5 * Y23
      L(4,3) = 0.5 * X13
      L(5,3) = 0.5 * Y31
      L(7,3) = 0.5 * X21
      L(8,3) = 0.5 * Y12
C
      IF (ALPHA.GT.0.0) THEN
         L(3,1) = 0.5 * (Y23*(Y13-Y21)*ALPHA/6)
         L(3,2) = 0.5 * (X32*(X31-X12)*ALPHA/6)
         L(3,3) = 0.5 * ((X31*Y13-X12*Y21)*ALPHA/3)
         L(6,1) = 0.5 * (Y31*(Y21-Y32)*ALPHA/6)
         L(6,2) = 0.5 * (X13*(X12-X23)*ALPHA/6)
         L(6,3) = 0.5 * ((X12*Y21-X23*Y32)*ALPHA/3)
         L(9,1) = 0.5 * (Y12*(Y32-Y13)*ALPHA/6)
         L(9,2) = 0.5 * (X21*(X23-X31)*ALPHA/6)
         L(9,3) = 0.5 * ((X23*Y32-X31*Y13)*ALPHA/3)
      ENDIF
C                                                         T
C ---------------------------- Generates  TMB = 1/AREA * L  --------
C
      DO 20 I=1,9
         DO 10 J=1,3
            TMB(J,I) = (1./AREA) * L(I,J)
   10    CONTINUE
   20 CONTINUE
C
C ---------------------------- Generates  BH(3,3) ------------------
C
      X0 = (X(1)+X(2)+X(3))/3.0
      Y0 = (Y(1)+Y(2)+Y(3))/3.0
      C = 1./DSQRT(AREA)
C
      XC(1) = C * (X(1)-X0)
      XC(2) = C * (X(2)-X0)
      XC(3) = C * (X(3)-X0)
      YC(1) = C * (Y(1)-Y0)
      YC(2) = C * (Y(2)-Y0)
      YC(3) = C * (Y(3)-Y0)
      XM(1) = 0.5 * (XC(2)+XC(3))
      XM(2) = 0.5 * (XC(3)+XC(1))
      XM(3) = 0.5 * (XC(1)+XC(2))
      YM(1) = 0.5 * (YC(2)+YC(3))
      YM(2) = 0.5 * (YC(3)+YC(1))
      YM(3) = 0.5 * (YC(1)+YC(2))
C
C
      DO 30 J=1,3
         DX = XM(J)-XC(J)
         DY = YM(J)-YC(J)
         DL = DSQRT (DX**2 + DY**2)
         CJ = DX/DL
         SJ = DY/DL
         A1J = -0.5*SJ*CJ**2
         A2J =  0.5*CJ**3
         B2J = -0.5*SJ**3
         B3J =  0.5*SJ**2*CJ
         A3J = -(B2J + A1J + A1J)
         B1J = -(B3J + B3J + A2J)
C
         BH(1,J) = C * ( 2*A1J*XC(J) +   A2J*YC(J) )
         BH(2,J) = C * (   B2J*XC(J) + 2*B3J*YC(J) )
         BH(3,J) = C * (-4*B3J*XC(J) - 4*A1J*YC(J) )
C
   30 CONTINUE
C
C ---------------------------- Generates  TMH = BH * HH ------------
C
      DO 60 I=1,9
         DO 50 J=1,3
            S = 0.0
            DO 40 K=1,3
               S = S + BH(J,K)*HH(K,I)
   40       CONTINUE
            TMH(J,I) = S
   50    CONTINUE
   60 CONTINUE
C
C ---------------------------- Generates  TM = TMB + TMH -----------
C
      DO 80 I=1,3
         DO 70 J=1,9
            TM(I,J) = TMB(I,J) + TMH(I,J)
   70    CONTINUE
   80 CONTINUE
C
C ---------------------------- Generates  SM = DM * TM -------------
C
      DO 110 I=1,9
         DO 100 J=1,3
            S = 0.0
            DO 90 K=1,3
               S = S + DM(J,K)*TM(K,I)
   90       CONTINUE
            SM(J,I) = S
  100    CONTINUE
  110 CONTINUE
C
C ---------------------------- Return ------------------------------
C
      RETURN
C
      END
      SUBROUTINE TMRF35(EM,X,Y,Z,THK,RHO)
C
C***********************************************************************
C
C     TMRF35 DETERMINES A LUMPED MASS MATRIX  EM  FOR THE TRIANGULAR
C     MEMBRANE ELEMENT WITH ROTATIONAL DEGREES OF FREEDOM OF P.G.BERGAN
C     AND C.A.FELIPPA.
C
C     PROGRAMMED BY :  KETIL AAMNES
C     DATE/VERSION  :  081085 / 1.0
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
      EM(2) = EM(1)
      EM(4) = EM(1)
      EM(5) = EM(1)
      EM(7) = EM(1)
      EM(8) = EM(1)
C
      EM(3) = GAMMA*IZ
      EM(6) = EM(3)
      EM(9) = EM(3)
C ---------------------------- Return ------------------------------
   80 RETURN
C
      END

