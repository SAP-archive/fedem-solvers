C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
C     This file is part of FEDEM - https://openfedem.org
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C Navn   : QMRF31
C
C Form}l: ] danne stivhetsmatrisen for membranelement med firkantet
C          geometri,4 knutepunkter,3 frihetsgrader pr. knutepunkt.
C          Utviklingen av stivheten er basert p} fri formulerings-
C          teori. QMRF31 henter basisstivhet fra SM4B og h|yere
C          ordens stivhet fra SM4H. Parameteren BETA er innf|rt
C          for } studere de h|yere ordens polynomers innflytelse
C          p} den totale stivheten av elementet.
C
C Input  : X,Y,DM,EMOD,NY,THK,ALFA,BETA
C
C Output : PEK
C
C Kaller : Subrutinen SM4B  (DANNKB)
C          Subrutinen SM4H  (DANNKH)
C          Subrutinen SKFEIL
C          Subrutinen MPRT30
C
C Versjon:01          Dato:85.10.03          Progr. av:B.Lund
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      SUBROUTINE QMRF31(PEK,X,Y,DM,EMOD,NY,THK,ALFA,BETA,IPSW,LPU,IERR)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION PEK(12,12),X(4),Y(4),DM(3,3),EMOD,NY,THK,
     $                 ALFA,BETA,GAMMA,KB(12,12),KH(12,12),
     $                 D(4),N(4),A1,A2,TETAX,DIMXL,DIMYL
C
      INTEGER I,J,IPSW,LPU,IERR
C
C--------------------------------------------------------------------
C
      IF (IPSW .GT. 2) THEN
         WRITE(LPU,6000) X,Y,((DM(I,J),J=1,3),I=1,3)
         WRITE(LPU,6010) EMOD,NY,THK,ALFA,BETA
      ENDIF
C
      IERR=0
C
C---------------------------------------------------------------------
C     Finn basisstivhet.
C---------------------------------------------------------------------
C
      CALL SM4B(KB,X,Y,DM,ALFA,IERR)
      IF(IERR .LT. 0)THEN
         CALL SKFEIL(IERR,LPU)
         GOTO 50
      ENDIF
      IF (IPSW .GT. 2) THEN
         WRITE(LPU,6099) 'Basic membrane stiffness:'
         CALL MPRT30(KB,12,12,LPU)
      ENDIF
C
C----------------------------------------------------------------------
C     Finn h|yere ordens stivhet.
C----------------------------------------------------------------------
C
      GAMMA=EMOD*THK/(1.0-NY*NY)
      CALL SM4H(KH,X,Y,DIMXL,DIMYL,GAMMA,NY,D,N,TETAX,A1,A2,IERR)
      IF(IERR .LT. 0)THEN
         CALL SKFEIL(IERR,LPU)
         GOTO 50
      ENDIF
      IF (IPSW .GT. 2) THEN
         WRITE(LPU,6099) 'Higher-order membrane stiffness:'
         CALL MPRT30(KH,12,12,LPU)
      ENDIF
C
C--------------------------------------------------------------------
C      Adder sammen bidrag til total elementstivhet.
C--------------------------------------------------------------------
C
      DO 30 I=1,12
         DO 40 J=1,12
            PEK(I,J)=KB(I,J)+BETA*KH(I,J)
 40      CONTINUE
 30   CONTINUE
C
      IF (IPSW .GT. 1) THEN
         WRITE(LPU,6099) 'Total membrane stiffness:'
         CALL MPRT30(PEK,12,12,LPU)
         WRITE(LPU,6100)
      ENDIF
C
 50   R E T U R N
C
 6000 FORMAT(// 5X,'ENTERING QMRF31:'
     +        / 5X,'X    =',1P4E13.5
     +        / 5X,'Y    =',1P4E13.5
     +        / 5X,'DM   =',1P3E13.5 / (11X,3E13.5) )
 6010 FORMAT(   5X,'EMOD =',1PE13.5
     +        / 5X,'NY   =',1PE13.5
     +        / 5X,'THK  =',E13.5
     +        / 5X,'ALFA =',E13.5
     +        / 5X,'BETA =',E13.5 )
 6099 FORMAT(/,1X,A)
 6100 FORMAT(// 5X,'LEAVING QMRF31' )
C
C  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **
C
C      Returnerer beregnet total elementstivhet.
C
C  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **
      E N D
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C Navn   : SM4B                               
C
C Form}l:  ] danne basisstivhet for generelt firkantet membranelement .
C          Elementet har 4 hj|rneknutepunkter med 3 frihetsgrader pr.
C          knutepunkt.Utviklingen av elementet er gjort vhja. s}kalt
C          fri formulerings teori.For } studere sammenhengen mellom
C          mekanisk hj|rnerotasjon og rotasjon av sidekant er faktoren
C          ALFA innf|rt.
C
C Input  : X,Y,DM,ALFA,IFLAG
C
C Output : KB,IFLAG
C
C
C Kaller : Ingen
C
C
C Versjon:02                Dato:85.09.25               Progr. av:B.Lund
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      SUBROUTINE SM4B(KB,X,Y,DM,ALFA,IFLAG)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION KB(12,12),X(4),Y(4),DM(3,3),ALFA,AREAL,LUMP(12,3)
      INTEGER          I,J,K,L,IFLAG
C
C-----------------------------------------------------------------------
C
      IFLAG=0
      AREAL=DABS((Y(1)*(X(2)-X(4))+Y(2)*(X(3)-X(1))+Y(3)*(X(4)-X(2))
     $      +Y(4)*(X(1)-X(3)))*0.5)
      IF(AREAL .LE. 0.)THEN
        IFLAG= IFLAG-1
        RETURN
      ENDIF
C
C-----------------------------------------------------------------------
C      -Etabler LUMP,angir sammenheng mellom kn.pkt.krefter og fordelte
C       krefter.
C-----------------------------------------------------------------------
C
      LUMP(1,1)=0.5*(Y(2)-Y(4))
      LUMP(1,2)=0.
      LUMP(1,3)=0.5*(X(4)-X(2))
      LUMP(2,1)=0.
      LUMP(2,2)=LUMP(1,3)
      LUMP(2,3)=LUMP(1,1)
      LUMP(3,1)=0.5*ALFA*((Y(1)-Y(4))**2-(Y(2)-Y(1))**2)/6.0
      LUMP(3,2)=0.5*ALFA*((X(4)-X(1))**2-(X(1)-X(2))**2)/6.0
      LUMP(3,3)=0.5*ALFA*((X(4)-X(1))*(Y(1)-Y(4))-
     $                    (X(1)-X(2))*(Y(2)-Y(1)))/3.0
      LUMP(4,1)=0.5*(Y(3)-Y(1))
      LUMP(4,2)=0.
      LUMP(4,3)=0.5*(X(1)-X(3))
      LUMP(5,1)=0.
      LUMP(5,2)=LUMP(4,3)
      LUMP(5,3)=LUMP(4,1)
      LUMP(6,1)=0.5*ALFA*((Y(2)-Y(1))**2-(Y(3)-Y(2))**2)/6.0
      LUMP(6,2)=0.5*ALFA*((X(1)-X(2))**2-(X(2)-X(3))**2)/6.0
      LUMP(6,3)=0.5*ALFA*((X(1)-X(2))*(Y(2)-Y(1))-
     $                    (X(2)-X(3))*(Y(3)-Y(2)))/3.0
      LUMP(7,1)=0.5*(Y(4)-Y(2))
      LUMP(7,2)=0.
      LUMP(7,3)=0.5*(X(2)-X(4))
      LUMP(8,1)=0.
      LUMP(8,2)=LUMP(7,3)
      LUMP(8,3)=LUMP(7,1)
      LUMP(9,1)=0.5*ALFA*((Y(3)-Y(2))**2-(Y(4)-Y(3))**2)/6.0
      LUMP(9,2)=0.5*ALFA*((X(2)-X(3))**2-(X(3)-X(4))**2)/6.0
      LUMP(9,3)=0.5*ALFA*((X(2)-X(3))*(Y(3)-Y(2))-
     $                    (X(3)-X(4))*(Y(4)-Y(3)))/3.0
      LUMP(10,1)=0.5*(Y(1)-Y(3))
      LUMP(10,2)=0.
      LUMP(10,3)=0.5*(X(3)-X(1))
      LUMP(11,1)=0.
      LUMP(11,2)=LUMP(10,3)
      LUMP(11,3)=LUMP(10,1)
      LUMP(12,1)=0.5*ALFA*((Y(4)-Y(3))**2-(Y(1)-Y(4))**2)/6.0
      LUMP(12,2)=0.5*ALFA*((X(3)-X(4))**2-(X(4)-X(1))**2)/6.0
      LUMP(12,3)=0.5*ALFA*((X(3)-X(4))*(Y(4)-Y(3))-
     $                     (X(4)-X(1))*(Y(1)-Y(4)))/3.0
C
C-----------------------------------------------------------------------
C     Adderer inn stivhetsbidrag.
C-----------------------------------------------------------------------
C
C
      DO 60 I=1,12
       DO 50 J=1,12
        KB(I,J)=0.
        DO 40 K=1,3
         DO 30 L=1,3
          KB(I,J)=KB(I,J)+(LUMP(J,K)*DM(L,K)*LUMP(I,L))/AREAL
 30      CONTINUE
 40     CONTINUE
 50    CONTINUE
 60   CONTINUE
C
      R E T U R N
C
C  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **
C
C     Returnerer beregnet basisstivhet.
C
C  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **
      E N D
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C Navn   : SM4H (DANNKH)                      
C
C Form}l : ] danne den h|yere ordens stivheten for membranelementet
C          med 4 kn.pkter/ 12 frh.gr. basert p} fri formuleringsteori.
C          Ved } innf|re skaleringsfaktorene DIMX,DIMY tas det 
C          hensyn til at elementet kan v{re langstrakt.
C
C
C
C Input  : X,Y,DIMXL,DIMYL,GAMMA,NY,IFLAG     DIMXL - beregnes autom.
C                                             DIMYL - beregnes autom.
C Output : KH,IFLAG
C
C
C Kaller : Subrutine DANKQH
C          Subrutine DANNHH
C
C
C
C Versjon:01                 Dato:85.09.25         Progr. av:B.Lund   
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      SUBROUTINE SM4H(KH,X,Y,DIMXL,DIMYL,GAMMA,NY,E,N,TETAX,A1,A2,IFLAG)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION   KH(12,12),X(4),Y(4),GAMMA,NY,A1,A2,
     $     AREAL,THETA,TETAX,
     $     XC,YC,XGLOB1,XGLOB2,XGLOB3,XGLOB4,YGLOB1,YGLOB2,
     $     HRC(6,12),YGLOB3,YGLOB4,E(4),N(4),KQH(6,6),HH(6,12),
     $     T(12,12),TEMP(12),DIMXL,DIMYL,DLETA,DLKSI
      INTEGER I,J,IFLAG,K,L
C
      DOUBLE PRECISION   PI
      PARAMETER        ( PI = 3.141592653589793D0 )
C
C-------------------------------------------------------------------
C
      IFLAG=0
      A1=0.5D0*DABS(X(2)*Y(3)-X(3)*Y(2)+
     $              X(3)*Y(1)-X(1)*Y(3)+
     $              X(1)*Y(2)-X(2)*Y(1))
      A2=0.5D0*DABS(X(3)*Y(4)-X(4)*Y(3)+
     $              X(4)*Y(1)-X(1)*Y(4)+
     $              X(1)*Y(3)-X(3)*Y(1))
      IF (A1 .LE. 0.0D0 .OR. A2 .LE. 0.0D0) THEN
         IFLAG=-2
         RETURN
      ENDIF
C
C--------------------------------------------------------------------
C     Orienterer lokalt aksesystem v.hj.a. minste kvadr. metode.
C     TETAX angir vinkel mellom globale retning og den lokale.
C--------------------------------------------------------------------
C
      TETAX=0.0D0
      DO 10 I=1,4
         J=1+MOD(I,4)
         THETA=DATAN2(Y(J)-Y(I),X(J)-X(I))
         IF (THETA.LT.0.0D0) THEN
            TETAX=TETAX+THETA+PI+PI
         ELSE
            TETAX=TETAX+THETA
         ENDIF
   10 CONTINUE
C
      TETAX=0.25D0*TETAX-0.75D0*PI
      IF (TETAX.GT.PI) THEN
         TETAX=TETAX-PI
      ELSE IF (TETAX.GT.0.5D0*PI) THEN
         TETAX=TETAX-0.5D0*PI
      ELSE IF (TETAX.LT.-0.5D0*PI) THEN
         TETAX=TETAX+PI
      ENDIF
C
C--------------------------------------------------------------------
C     Beregner flatesenterets koordinater XC og YC.
C--------------------------------------------------------------------
C
      AREAL=A1+A2
      XC=(A1*(X(1)+X(2)+X(3))+A2*(X(3)+X(4)+X(1)))/(3.0*AREAL)
      YC=(A1*(Y(1)+Y(2)+Y(3))+A2*(Y(3)+Y(4)+Y(1)))/(3.0*AREAL)
C
C--------------------------------------------------------------------
C     Orienterer aksesystem og innf|rer dimensjonsl|se koord.
C     Bestemmer skaleringsfaktorene i de to retningene.
C--------------------------------------------------------------------
C
      XGLOB1=X(1)-XC
      XGLOB2=X(2)-XC
      XGLOB3=X(3)-XC
      XGLOB4=X(4)-XC
      YGLOB1=Y(1)-YC
      YGLOB2=Y(2)-YC
      YGLOB3=Y(3)-YC
      YGLOB4=Y(4)-YC
C
C      IF(XGLOB2 .GT. 0. .AND. XGLOB1 .LT. 0.)THEN
C       DLKSI=(XGLOB2+XGLOB3-XGLOB1-XGLOB4)/4.
C       DLETA=(YGLOB4+YGLOB3-YGLOB1-YGLOB2)/4.
C      ENDIF
C      IF(XGLOB2 .GT. 0. .AND. XGLOB1 .GT. 0.)THEN
C       DLKSI=(XGLOB1+XGLOB2-XGLOB3-XGLOB4)/4.
C       DLETA=(YGLOB3+YGLOB2-YGLOB4-YGLOB1)/4.
C      ENDIF
C      IF(XGLOB2 .LT. 0. .AND. XGLOB1 .GT. 0.)THEN
C       DLKSI=(XGLOB1+XGLOB4-XGLOB2-XGLOB3)/4.
C       DLETA=(YGLOB1+YGLOB2-YGLOB3-YGLOB4)/4.
C      ENDIF
C      IF(XGLOB2 .LT. 0. .AND. XGLOB1 .LT. 0.)THEN
C       DLKSI=(XGLOB4+XGLOB3-XGLOB1-XGLOB2)/4.
C       DLETA=(YGLOB1+YGLOB4-YGLOB2-YGLOB3)/4.     
C      ENDIF
C
C ** NEW TR 3.9.90 **
C
      DLKSI = ABS(XGLOB2+XGLOB3-XGLOB1-XGLOB4) / 4.0D0
      DLETA = ABS(YGLOB3+YGLOB4-YGLOB1-YGLOB2) / 4.0D0
C            
C ** END NEW
C
      DIMXL=1.0D0/DLKSI
      DIMYL=1.0D0/DLETA
C
      E(1)=DIMXL*(XGLOB1*DCOS(TETAX)+YGLOB1*DSIN(TETAX))
      E(2)=DIMXL*(XGLOB2*DCOS(TETAX)+YGLOB2*DSIN(TETAX))
      E(3)=DIMXL*(XGLOB3*DCOS(TETAX)+YGLOB3*DSIN(TETAX))
      E(4)=DIMXL*(XGLOB4*DCOS(TETAX)+YGLOB4*DSIN(TETAX))
      N(1)=DIMYL*(-XGLOB1*DSIN(TETAX)+YGLOB1*DCOS(TETAX))
      N(2)=DIMYL*(-XGLOB2*DSIN(TETAX)+YGLOB2*DCOS(TETAX))
      N(3)=DIMYL*(-XGLOB3*DSIN(TETAX)+YGLOB3*DCOS(TETAX))
      N(4)=DIMYL*(-XGLOB4*DSIN(TETAX)+YGLOB4*DCOS(TETAX))
C
C--------------------------------------------------------------------
C     Finn generalisert h|yere ordens stivhet og HH.(Evnt.HRC).
C     Beregn s} den h|yere ordens stivheten referert til det
C     globale aksesystemet.
C--------------------------------------------------------------------
C
      CALL DANKQH(KQH,E,N,NY,DIMXL,DIMYL,A1,A2,GAMMA)
      CALL DANNHH(HH,HRC,E,N,DIMXL,DIMYL,TETAX,T,IFLAG)
      IF (IFLAG .LT. 0) RETURN
C
      DO 60 I=1,12
       DO 50 J=1,12
        KH(I,J)=0.0D0
        DO 40 K=1,6
         DO 30 L=1,6
          KH(I,J)=KH(I,J)+HH(K,J)*KQH(L,K)*HH(L,I)
 30      CONTINUE                       
 40     CONTINUE                       
 50    CONTINUE                       
 60   CONTINUE 
C
C---------------------------------------------------------------
C     Transformer stivhet til globalt system.
C---------------------------------------------------------------
C
      CALL PREMUL(T,KH,TEMP,12,12,2)
      CALL PSTMUL(T,KH,TEMP,12,12,1)
C
      R E T U R N
C
C  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **
C
C     Returnerer beregnet h|yere ordens stivhet.
C
C  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **
      E N D
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C Navn   : DANNHH               
C
C Form}l:  ] danne sammenhengen mellom generaliserte frihetsgr.
C          og fysiske frihetsgr.Denne relasjonen lagres i HH,
C          som benyttes til } danne den h|yere ordens stivheten.
C          HH framkommer ved invertering av G-matrisen. I tillegg
C          etableres ogs} transformasjonsmatrisen T mellom akse-
C          systemene.
C
C Input  : E,N,DIMX,DIMY,TETAX
C
C Output : HH,HRC,T,IERR
C
C
C Kaller : Subrutinen DANNG
C          Subrutinen DINV12
C          Subrutinen INSUB / S A M
C
C
C
C Versjon: 01        Dato:85.10.03             Progr. av:B.Lund       
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      SUBROUTINE DANNHH (HH,HRC,E,N,DIMX,DIMY,TETAX,T,IERR)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION   HH(6,12),HRC(6,12),E(4),N(4),DIMX,DIMY,
     $                   TETAX,T(12,12)
      DOUBLE PRECISION   T0(3,3),GINV(12,12),WA(12)
      INTEGER            IERR,I,J,IWA(12)
C
C-------------------------------------------------------------------
C     Etablerer transformasjonsmatrisen T.
C-------------------------------------------------------------------
C
      T0(1,1)=DCOS(TETAX)
      T0(1,2)=DSIN(TETAX)
      T0(1,3)=0.
      T0(2,1)=-T0(1,2)
      T0(2,2)=T0(1,1)
      T0(2,3)=0.
      T0(3,1)=0.
      T0(3,2)=0.
      T0(3,3)=1.
C
      DO 20 I=1,12
       DO 10 J=1,12
        T(I,J)=0.
 10    CONTINUE
 20   CONTINUE
C
      CALL INSUB(T,T0,12,12,3,3,1,1,1)
      CALL INSUB(T,T0,12,12,3,3,4,4,1)
      CALL INSUB(T,T0,12,12,3,3,7,7,1)
      CALL INSUB(T,T0,12,12,3,3,10,10,1)
C                                                                  
C---------------------------------------------------------------
C     Dann G-matrisen.
C---------------------------------------------------------------
C
      CALL DANNG(GINV,E,N,DIMX,DIMY)
C
C---------------------------------------------------------------
C     Inverter G og del opp den inverterte matrisen i HH og HRC.
C     Sjekk at G ikke er singul{r eller at G befinner seg sv{rt
C     n{r singularitet.
C---------------------------------------------------------------
C
      CALL DINV12(12,GINV,WA,IWA,0,IERR)
      IF (IERR .NE. 0) THEN
         IERR = -3
         RETURN
      ENDIF
C
C----------------------------------------------------------------
C     HRC utgj|r |vre submtrise,HH nedre.
C----------------------------------------------------------------
C
      DO 50 I=1,6
       DO 40 J=1,12
        HRC(I,J)=GINV(I,J)
 40    CONTINUE
 50   CONTINUE
C
      DO 70 I=7,12
       DO 60 J=1,12
        HH(I-6,J)=GINV(I,J)
 60    CONTINUE
 70   CONTINUE
C
      R E T U R N
C
C
C  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **
C
C     Returnerer HH og HRC.
C
C  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **
      E N D  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C Navn   : DANKQH               
C
C Form}l : ] danne den lokale, generaliserte stivhetsmatrisen for
C          fri formuleringselementet.Samtlige aktuelle polynomledd
C          integreres eksakt av DANINT.De integrerte leddene til-
C          h|rende koeffisientene i KQH beregnes vhja. VEKT.
C
C
C Input  : E,N NY,DIMX,DIMY,A1,A2,GAMMA
C
C Output : KQH
C
C
C Kaller : Subrutinen DANINT
C          Subrutinen VEKT
C
C
C Versjon:01                 Dato:85.10.02         Progr. av:B.Lund   
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      SUBROUTINE DANKQH(KQH,E,N,NY,DIMX,DIMY,A1,A2,GAMMA)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION KQH(6,6),E(4),N(4),NY,DIMX,DIMY,A1,A2,GAMMA,
     $                 N201,N202,E1N101,E1N102,E201,E202,N301,N302,
     $                 E1N201,E1N202,E2N101,E2N102,E301,E302,N401,
     $                 N402,E1N301,E1N302,E2N201,E2N202,E3N101,
     $                 E3N102,E401,E402,COEFF(105),BHL1(6),
     $                 BHL2(6),BHL3(6),FI
      INTEGER          IPEK(6,6),I,J
C
C---------------------------------------------------------------------
C    Danner B-matrisen som 3 vektorer.
C---------------------------------------------------------------------
C
      BHL1(1)=DIMX
      BHL1(2)=0.
      BHL1(3)=0.
      BHL1(4)=0.
      BHL1(5)=2.*DIMX
      BHL1(6)=-2.*DIMX
      BHL2(1)=0.
      BHL2(2)=0.
      BHL2(3)=DIMY
      BHL2(4)=0.
      BHL2(5)=2.*DIMY
      BHL2(6)=2.*DIMY
      BHL3(1)=DIMY
      BHL3(2)=2.*DIMY
      BHL3(3)=DIMX
      BHL3(4)=2.*DIMX
      BHL3(5)=DIMY
      BHL3(6)=DIMY
C
C------------------------------------------------------------------
C     Beregn integrerte polynomledd.
C------------------------------------------------------------------
C
      CALL DANINT(E,N,N201,N202,E1N101,E1N102,E201,E202,N301,N302,
     $                 E1N201,E1N202,E2N101,E2N102,E301,E302,N401,
     $                 N402,E1N301,E1N302,E2N201,E2N202,E3N101,
     $                 E3N102,E401,E402)
C
      FI=DIMX/DIMY
C 
      CALL VEKT(COEFF,IPEK,FI,A1,A2,
     $                 N201,N202,E1N101,E1N102,E201,E202,N301,N302,
     $                 E1N201,E1N202,E2N101,E2N102,E301,E302,N401,
     $                 N402,E1N301,E1N302,E2N201,E2N202,E3N101,
     $                 E3N102,E401,E402)
C
C--------------------------------------------------------------------
C     Dann stivhetsledd.
C--------------------------------------------------------------------
C
      DO 40 I=1,6
       DO 30 J=1,I
        KQH(I,J)=GAMMA*(BHL1(I)*BHL1(J)*COEFF(IPEK(I,J))+NY*BHL2(I)*
     $                  BHL1(J)*COEFF(IPEK(I,J)+1)+NY*BHL1(I)*BHL2(J)*
     $                  COEFF(IPEK(I,J)+2)+BHL2(I)*BHL2(J)*
     $                 COEFF(IPEK(I,J)+3)+0.5*(1.-NY)*BHL3(I)*BHL3(J)*
     $                  COEFF(IPEK(I,J)+4))
 30    CONTINUE
 40   CONTINUE
C
C------------------------------------------------------------------
C     Utnytt symmetri.
C------------------------------------------------------------------
C
      DO 60 I=1,6
       DO 50 J=1,I
        KQH(J,I)=KQH(I,J)
 50    CONTINUE
 60   CONTINUE
C
      R E T U R N
C
C  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **
C
C     Returnerer beregnet generalisert h|yere ordens stivhet.
C
C  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **
      E N D    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C Navn   : DANNG                 
C
C Form}l:  ] danne G-matrisen for b}de lavere og h|yere ordens
C          forskyvningspolynomer for det generelle firkant-
C          elementet med rotasjonsfrihetsgrader. G-matrisen
C          angir sammenhengen mellom knutepunktsforskyvningene
C          og de generaliserte forskyvningsm|nsterne.
C
C
C Input  : E,N,DIMX,DIMY
C
C Output : G
C
C
C Kaller : Ingen.
C
C
C Versjon:  01           Dato:85.10.03         Progr. av:B.Lund       
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      SUBROUTINE DANNG(G,E,N,DIMX,DIMY)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION G(12,12),E(4),N(4),DIMX,DIMY,GHL(12,6),GRC(12,6)
      INTEGER          I,J,K
C
C---------------------------------------------------------------------
C     Utnytter "symmetri" i GHL og GRC og setter disse sammen til G.
C---------------------------------------------------------------------
C
      K=1
C
      DO 10 I=1,4
       GHL(K,1)  =E(I)*N(I)     
       GHL(K+1,1)=0.     
       GHL(K+2,1)=-0.5*DIMY*E(I)     
       GHL(K,2)  =N(I)**2     
       GHL(K+1,2)=0.     
       GHL(K+2,2)=-DIMY*N(I)     
       GHL(K,3)  =0.     
       GHL(K+1,3)=E(I)*N(I)     
       GHL(K+2,3)=0.5*DIMX*N(I)     
       GHL(K,4)  =0.     
       GHL(K+1,4)=E(I)**2     
       GHL(K+2,4)=DIMX*E(I)     
       GHL(K,5)  =N(I)*((E(I)-N(I))**2)
       GHL(K+1,5)=-E(I)*((N(I)-E(I))**2)     
       GHL(K+2,5)=0.5*(DIMX*(-(N(I)**2)+4.*E(I)*N(I)-(3.*E(I)
     $                **2))-DIMY*((E(I)**2)-4.*E(I)*N(I)+(3.*
     $                N(I)**2)))         
       GHL(K,6)  =-N(I)*((E(I)+N(I))**2)     
       GHL(K+1,6)=E(I)*((N(I)+E(I))**2)     
       GHL(K+2,6)=0.5*(DIMX*((N(I)**2)+4.*E(I)*N(I)+(3.*E(I)
     $                 **2))+DIMY*((E(I)**2)+4.*E(I)*N(I)+(3.*
     $                 N(I)**2)))     
       GRC(K,1)  =1.     
       GRC(K+1,1)=0.     
       GRC(K+2,1)=0.     
       GRC(K,2)  =0.     
       GRC(K+1,2)=1.     
       GRC(K+2,2)=0.     
       GRC(K,3)  =-N(I)
       GRC(K+1,3)=E(I)     
       GRC(K+2,3)=0.5*(DIMX+DIMY)     
       GRC(K,4)  =E(I)
       GRC(K+1,4)=0.
       GRC(K+2,4)=0.            
       GRC(K,5)  =0.     
       GRC(K+1,5)=N(I)     
       GRC(K+2,5)=0.     
       GRC(K,6)  =N(I)     
       GRC(K+1,6)=E(I)     
       GRC(K+2,6)=0.5*(DIMX-DIMY)     
       K=K+3
 10   CONTINUE
C
C------------------------------------------------------------------
C     GRC danner venstre halvdel av G , og GHL danner h|yre.
C------------------------------------------------------------------
C
      DO 30 I=1,12
       DO 20 J=1,6
        G(I,J)=GRC(I,J)
 20    CONTINUE
 30   CONTINUE
C
      DO 50 I=1,12
       DO 40 J=7,12
        G(I,J)=GHL(I,J-6)
 40    CONTINUE
 50   CONTINUE
C
      R E T U R N
C
C  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **
C
C     Returnerer G-matrisen.
C
C  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **
      E N D
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C Navn   : DANINT              
C
C Form}l:  ] integrere forskyvningspolynomene for den h|yere ordens
C          stivheten eksakt.Integralene beregnes ved hjelp av en
C          binomialformel.Integrasjonen over det generelle firkant-  
C          elementet deles opp i integrasjon over to trekanter.
C          Hhv. (1,2,3) og (3,4,1) n}r kn.pktene er nummerert mot
C          urviseren.Delintegralene lagres i POLYN1 og POLYN2.
C
C Input  : E,N
C
C Output : Enkeltvariable som inneholder samtlige aktuelle polynom-
C          ledd integrert.
C
C Kaller : Subrutinen FAKUL
C
C
C Versjon: 01                 Dato: 85.09.26          Progr. av:B.Lund
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      SUBROUTINE DANINT(E,N,N201,N202,E1N101,E1N102,E201,E202,N301,N302
     $                 ,E1N201,E1N202,E2N101,E2N102,E301,E302,N401,N402
     $                 ,E1N301,E1N302,E2N201,E2N202,E3N101,E3N102,E401,
     $                  E402)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION E(4),N(4),N201,N202,E1N101,E1N102,E201,E202,N301
     $           ,N302,E1N201,E1N202,E2N101,E2N102,E301,E302,N401,N402,
     $             E1N301,E1N302,E2N201,E2N202,E3N101,E3N102,E401,E402,
     $        EKSP1,EKSP2,EKSP3,EKSP4,EKSP5,EKSP6,POLYN1(15),POLYN2(15)
      INTEGER IA,IB,IC,ID,IEN,IF,IG,IH,IK,IL,IM,INE,IP,IQ,IFAK(14),I
C
C----------------------------------------------------------------------
C
      DO 10 I=1,15
       POLYN1(I)=0.
       POLYN2(I)=0.
 10   CONTINUE
C
C----------------------------------------------------------------------
C     Beregn binom.form. over de to trekantene.Sjekker eksponenter. 
C----------------------------------------------------------------------
C
      I=1
      DO 70 INE=0,4
       DO 60 IM=0,4-INE
        DO 50 IK=0,INE
         DO 40 IP=0,IM
          DO 30 IL=0,IK
           DO 20 IQ=0,IP
            IA=INE-IK
            IB=IK-IL
            IC=IM-IP
            ID=IP-IQ
            IEN=IL+IQ
            IF=IK-IL+IP-IQ
            IG=INE-IK+IM-IP
            IH=IM+INE+2
             IF(INE.LT.IK)IA=0
             IF(IK.LT.IL)IB=0
             IF(IM.LT.IP)IC=0
             IF(IP.LT.IQ)ID=0
             IF(IF.LT.0)IF=0
             IF(IG.LT.0)IG=0
            CALL FAKUL(IA,IB,IC,ID,IEN,IF,IG,IH,IK,IL,IM,INE,IP,
     $                 IQ,IFAK)
             EKSP1=1.
             EKSP2=1.
             EKSP3=1.
             EKSP4=1.
             EKSP5=1.
             EKSP6=1.
             IF(IL.NE.0)EKSP1=E(1)**IL
             IF(IB.NE.0)EKSP2=E(2)**IB
             IF(IA.NE.0)EKSP3=E(3)**IA
             IF(IQ.NE.0)EKSP4=N(1)**IQ
             IF(ID.NE.0)EKSP5=N(2)**ID
             IF(IC.NE.0)EKSP6=N(3)**IC
           POLYN1(I)=POLYN1(I)+(2.*IFAK(12)*IFAK(9)*IFAK(11)*
     $                           IFAK(13)*IFAK(5)
     $                          *IFAK(6)*IFAK(7)*(EKSP1)*(EKSP2)
     $                          *(EKSP3)*(EKSP4)*(EKSP5)
     $                          *(EKSP6))/(IFAK(9)*IFAK(1)*IFAK(10)
     $                          *IFAK(2)*IFAK(13)*IFAK(3)*IFAK(14)*
     $                           IFAK(4)*IFAK(8))
 20        CONTINUE
 30       CONTINUE
 40      CONTINUE
 50     CONTINUE
        I=I+1
 60    CONTINUE
 70   CONTINUE     
C
C
      I=1
      DO 130 INE=0,4
       DO 120 IM=0,4-INE
        DO 110 IK=0,INE
         DO 100 IP=0,IM
          DO 90 IL=0,IK
           DO 80 IQ=0,IP
            IA=INE-IK
            IB=IK-IL
            IC=IM-IP
            ID=IP-IQ
            IEN=IL+IQ
            IF=IK-IL+IP-IQ
            IG=INE-IK+IM-IP
            IH=IM+INE+2
             IF(INE.LT.IK)IA=0
             IF(IK.LT.IL)IB=0
             IF(IM.LT.IP)IC=0
             IF(IP.LT.IQ)ID=0
             IF(IF.LT.0)IF=0
             IF(IG.LT.0)IG=0
            CALL FAKUL(IA,IB,IC,ID,IEN,IF,IG,IH,IK,IL,IM,INE,IP
     $                ,IQ,IFAK)
             EKSP1=1.
             EKSP2=1.
             EKSP3=1.
             EKSP4=1.
             EKSP5=1.
             EKSP6=1.
             IF(IL.NE.0)EKSP1=E(3)**IL
             IF(IB.NE.0)EKSP2=E(4)**IB
             IF(IA.NE.0)EKSP3=E(1)**IA
             IF(IQ.NE.0)EKSP4=N(3)**IQ
             IF(ID.NE.0)EKSP5=N(4)**ID
             IF(IC.NE.0)EKSP6=N(1)**IC
           POLYN2(I)=POLYN2(I)+(2.*IFAK(12)*IFAK(9)*IFAK(11)*
     $                          IFAK(13)*IFAK(5)
     $                          *(EKSP3)*(EKSP4)*(EKSP5)
     $                          *IFAK(6)*IFAK(7)*(EKSP1)*(EKSP2)
     $                          *(EKSP6))/(IFAK(9)*IFAK(1)*IFAK(10)
     $                          *IFAK(2)*IFAK(13)*IFAK(3)*
     $                           IFAK(14)*IFAK(4)*IFAK(8))
 80         CONTINUE
 90        CONTINUE
 100      CONTINUE
 110     CONTINUE
         I=I+1
 120    CONTINUE
 130   CONTINUE
C
C----------------------------------------------------------------------
C      Legg integralene i enkeltvariable.
C----------------------------------------------------------------------
C
       N201=POLYN1(3)
       N202=POLYN2(3)
       E1N101=POLYN1(7)
       E1N102=POLYN2(7)
       E201=POLYN1(10)
       E202=POLYN2(10)
       N301=POLYN1(4)
       N302=POLYN2(4)
       E1N201=POLYN1(8)
       E1N202=POLYN2(8)
       E2N101=POLYN1(11)
       E2N102=POLYN2(11)
       E301=POLYN1(13)
       E302=POLYN2(13)
       N401=POLYN1(5)
       N402=POLYN2(5)
       E1N301=POLYN1(9)
       E1N302=POLYN2(9)
       E2N201=POLYN1(12)
       E2N202=POLYN2(12)
       E3N101=POLYN1(14)
       E3N102=POLYN2(14)
       E401=POLYN1(15)
       E402=POLYN2(15)
C
       R E T U R N
C
C  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **
C
C      Returnerer integrerte polynomledd.
C
C  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **
       E N D
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C Navn   : VEKT                
C
C Form}l : ] finne integralverdier til koeffisientene i KQH.
C          Pekertabellen IPEK etableres ogs}.
C
C
C Input  : Integrerte polynomledd,A1,A2
C
C Output : KOEFF
C
C
C Kaller : Ingen
C
C
C Versjon: 01                 Dato: 85.10.02       Progr. av:B.Lund    
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      SUBROUTINE VEKT(COEFF,IPEK,FI,A1,A2,N201,N202,E1N101,E1N102,E201
     $          ,E202,N301,N302,E1N201,E1N202,E2N101,E2N102,E301,E302,
     $                N401,N402,E1N301,E1N302,E2N201,E2N202,E3N101,
     $                E3N102,E401,E402)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION COEFF(105),A1,A2,N201,N202,E1N101,E1N102,E201,
     $                 E202,N301,N302,E1N201,E1N202,E2N101,E2N102,
     $                 E301,E302,N401,N402,E1N301,E1N302,E2N201,
     $                 E2N202,E3N101,E3N102,E401,E402,FI
      INTEGER          IPEK(6,6),I,J,ISUM
C
C---------------------------------------------------------------------
C     Dann pekertabellen IPEK, og etabler COEFF.
C---------------------------------------------------------------------
C
      ISUM=-4
      DO 20 I=1,6
       DO 10 J=1,I
        ISUM=ISUM+5
        IPEK(I,J)=ISUM
 10    CONTINUE
 20   CONTINUE      
C
      DO 30 I=1,105
       COEFF(I)=0.
30    CONTINUE
C
      COEFF(1)=A1*N201+A2*N202
      COEFF(5)=A1*E201+A2*E202
      COEFF(10)=A1*E1N101+A2*E1N102
      COEFF(15)=COEFF(1)
      COEFF(17)=COEFF(10)
      COEFF(20)=COEFF(10)
      COEFF(25)=COEFF(1)
      COEFF(29)=COEFF(5)
      COEFF(30)=COEFF(1)
      COEFF(35)=COEFF(5)
      COEFF(40)=COEFF(10)
      COEFF(45)=COEFF(10)
      COEFF(50)=COEFF(5)
      COEFF(51)=A1*E1N201+A2*E1N202-(A1*N301+A2*N302)
      COEFF(52)=-1.0*(A1*E1N201+A2*E1N202)+1.*(A1*E2N101+A2*E2N102)
      COEFF(55)=(A1*E301+A2*E302)*(1.-3.*FI)-4.*(1.-FI)*(A1*E2N101+
     $            A2*E2N102)+(3.-FI)*(A1*E1N201+A2*E1N202)
      COEFF(60)=(1.-3.*FI)*(A1*E2N101+A2*E2N102)+(4.*FI-4.)*
     $            (A1*E1N201+A2*E1N202)+(3.-FI)*(A1*N301+A2*N302)
      COEFF(63)=A1*E2N101+A2*E2N102-(A1*E1N201+A2*E1N202)
      COEFF(64)=-1.*(A1*E2N101+A2*E2N102)+1.*(A1*E301+A2*E302)
      COEFF(65)=COEFF(60)
      COEFF(70)=COEFF(55) 
      COEFF(71)=A1*E2N201+A2*E2N202-2.*(A1*E1N301+A2*E1N302)+(A1*
     $            N401+A2*N402) 
      COEFF(72)=-2.*(A1*E2N201+A2*E2N202)+1.*(A1*E3N101+A2*E3N102)
     $           +1.*(A1*E1N301+A2*E1N302) 
      COEFF(73)=COEFF(72) 
      COEFF(74)=1.*(A1*E401+A2*E402)-2.*(A1*E3N101+A2*E3N102)+
     $           1.*(A1*E2N201+A2*E2N202) 
      COEFF(75)=((1.-3.*FI)**2)*(A1*E401+A2*E402)+((1.-3.*FI)*
     $ (4.*FI-4.)+(4.*FI-4.)*(1.-3.*FI))*(A1*E3N101+A2*E3N102)+
     $            ((1.-3.*FI)*(3.-FI)+(3.-FI)*(1.-3.*FI)+((4.*FI-4.)**
     $            2))*(A1*E2N201+A2*E2N202)+((4.*FI-4.)*(3.-FI)+
     $            (3.-FI)*(4.*FI-4.))*(A1*E1N301+A2*E1N302)+((3.-FI)
     $            **2)*(A1*N401+A2*N402) 
      COEFF(76)=A1*E1N201+A2*E1N202+A1*N301+A2*N302 
      COEFF(77)=1.*(A1*E1N201+A2*E1N202)+1.*(A1*E2N101+A2*E2N102) 
      COEFF(80)=(3.*FI-1.)*(A1*E301+A2*E302)+(4.*FI-4.)*(A1*E2N101+
     $            A2*E2N102)+(FI-3.)*(A1*E1N201+A2*E1N202) 
      COEFF(85)=(3.*FI-1.)*(A1*E2N101+A2*E2N102)+(4.*FI-4.)*(A1*
     $            E1N201+A2*E1N202)+(FI-3.)*(A1*N301+A2*N302) 
      COEFF(88)=A1*E2N101+A2*E2N102+A1*E1N201+A2*E1N202 
      COEFF(89)=1.*(A1*E2N101+A2*E2N102)+1.*(A1*E301+A2*E302) 
      COEFF(90)=COEFF(85) 
      COEFF(95)=COEFF(80) 
      COEFF(96)=A1*E2N201+A2*E2N202-(A1*N401+A2*N402) 
      COEFF(97)=1.*(A1*E3N101+A2*E3N102)-
     $             1.*(A1*E1N301+A2*E1N302) 
      COEFF(98)= 1.*(A1*E3N101+A2*E3N102)-
     $             1.*(A1*E1N301+A2*E1N302)
      COEFF(99)=1.*(A1*E401+A2*E402)-1.*(A1*E2N201+A2*E2N202) 
      COEFF(100)=(1.-3.*FI)*(3.*FI-1)*(A1*E401+A2*E402)+((1.-3.*FI)*
     $             (4.*FI-4.)+(4.*FI-4.)*(3.*FI-1.))*(A1*E3N101+A2*
     $             E3N102)+((1.-3.*FI)*(FI-3.)+((4.*FI-4.)**2)+
     $             (3.-FI)*(3.*FI-1.))*(A1*E2N201+A2*E2N202)+
     $             ((4.-4.*FI)*(FI-3.)+(3.-FI)*(4.*FI-4.))*(A1*
     $             E1N301+A2*E1N302)-((3.-FI)**2)*(A1*N401+A2*N402) 
      COEFF(101)=A1*E2N201+A2*E2N202+2.*(A1*E1N301+A2*E1N302)+
     $            (A1*N401+A2*N402) 
      COEFF(102)=2.*(A1*E2N201+A2*E2N202)+1.*(A1*E3N101+A2*E3N102)+
     $            1.*(A1*E1N301+A2*E1N302)
      COEFF(103)=COEFF(102) 
      COEFF(104)=1.*(A1*E2N201+A2*E2N202)+2.*(A1*E3N101+A2*E3N102)+
     $            1.*(A1*E401+A2*E402) 
      COEFF(105)=((3.*FI-1.)**2)*(A1*E401+A2*E402)+((3.*FI-1.)*(4.*FI-
     $             4.)+(4.*FI-4.)*(3.*FI-1.))*(A1*E3N101+A2*E3N102)+
     $             ((3.*FI-1.)*(FI-3.)+((4.*FI-4.)**2)+(FI-3.)*(3.*
     $             FI-1.))*(A1*E2N201+A2*E2N202)+((4.*FI-4.)*(FI-3.)+
     $             (FI-3.)*(4.*FI-4.))*(A1*E1N301+A2*E1N302)+((FI-3.)
     $             **2)*(A1*N401+A2*N402)
C
      R E T U R N
C
C  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **
C
C     Returnerer integralverdier til koeffisientene i KQH.
C
C  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **
      E N D 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C Navn   : FAKUL              
C
C Form}l:  ] beregne fakultetkoeffisientene i binomialformelen
C          som benyttes til } integrere forskyvningspolynomene
C          eksakt over elementet.Fakultetsverdiene lagres i
C          tabellen IFAK.
C
C
C Input  : IA,IB,IC,ID,IE,IF,IG,IH,IK,IL,IM,IN,IP,IQ
C
C Output : IFAK
C
C
C Kaller : Ingen
C
C
C Versjon: 01                Dato: 85.09.26        Progr. av:B.Lund
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      SUBROUTINE FAKUL(IA,IB,IC,ID,IEN,IF,IG,IH,IK,IL,IM,INE,IP,IQ,IFAK)
C
      IMPLICIT NONE
C
      INTEGER IA,IB,IC,ID,IEN,IF,IG,IH,IK,IL,IM,INE,IP,IQ,
     $        IFAK(14),ITALL(14),IGRENS,I,J
C
C----------------------------------------------------------------------
C
      ITALL(1)=IA
      ITALL(2)=IB
      ITALL(3)=IC
      ITALL(4)=ID
      ITALL(5)=IEN
      ITALL(6)=IF
      ITALL(7)=IG
      ITALL(8)=IH
      ITALL(9)=IK
      ITALL(10)=IL
      ITALL(11)=IM
      ITALL(12)=INE
      ITALL(13)=IP
      ITALL(14)=IQ
C
      DO 20 I=1,14
       IF(ITALL(I) .EQ. 0)THEN
        IFAK(I)=1
        GOTO 15
       ENDIF
       IFAK(I)=1
       IGRENS=ITALL(I)
       DO 10 J=1,IGRENS
        IFAK(I)=IFAK(I)*J
 10    CONTINUE
 15    CONTINUE
 20   CONTINUE
C
C
      R E T U R N
C
C  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **
C
C     Returnerer fakultetsverdier til binomialformelen.
C
C  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **
      E N D          
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C Navn   : SKFEIL                    
C
C Form}l : ] skrive feilmeldinger fra QMFR31 (DANNK).
C
C    *********  Denne versjonen skriver feilmeldinger  *********
C    *********  ogs} fra LAGK2. (Altn. II).            *********
C
C
C Input  : IFLAG
C
C Output : Ingen
C
C
C Kaller : Ingen
C
C
C Versjon: 02         Dato:85.10.04               Progr. av:B.Lund   
C Versjon: 03         Dato:02.10.02               Progr. av:K.Okstad
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      SUBROUTINE SKFEIL(IFLAG,LPU)
C
      IMPLICIT NONE
C
      INTEGER IFLAG,LPU
C
C---------------------------------------------------------------------
C     Lag overskrift og spesifiser hva slags feil som er funnet.
C---------------------------------------------------------------------
C
      IF (IFLAG .GE. 0) THEN
         RETURN
      ELSE IF (IFLAG .GT. -11) THEN
         WRITE(LPU,1000)
         WRITE(LPU,1100)
      ELSE IF (IFLAG .GT. -21) THEN
         WRITE(LPU,2000)
         WRITE(LPU,2100)
      ELSE IF (IFLAG .GT. -31) THEN
         WRITE(LPU,3000)
         WRITE(LPU,3100)
      ELSE IF (IFLAG .GT. -41) THEN
         WRITE(LPU,4000)
         WRITE(LPU,4100)
      ELSE IF (IFLAG .GT. -51) THEN
         WRITE(LPU,5000)
         WRITE(LPU,5100)
      ELSE
         RETURN
      ENDIF
C
      GOTO( 1, 2, 3, 4, 5,99,99,99,99,99,
     +     11,12,13,14,99,99,99,99,99,99,
     +     21,22,23,99,99,99,99,99,99,99,
     +     31,32,33,34,99,99,99,99,99,99,
     +     41,42,99,99,99,99,99,99,99,99) -IFLAG
C
C---------------------------------------------------------------------
C     Feil i stivhetsberegningen / Alt. 1
C---------------------------------------------------------------------
    1 WRITE(LPU,1200)
      RETURN
    2 WRITE(LPU,1300)
      RETURN
    3 WRITE(LPU,1400)
      RETURN
    4 WRITE(LPU,1500)
      RETURN
    5 WRITE(LPU,1600)
      RETURN
C
C----------------------------------------------------------------
C     Feilmelding for fordelt belastning
C----------------------------------------------------------------
   11 WRITE(LPU,2200)
      RETURN
   12 WRITE(LPU,2300)
      RETURN
   13 WRITE(LPU,2400)
      RETURN
   14 WRITE(LPU,2400)
      RETURN
C
C-----------------------------------------------------------------
C     Feil i spenningsberegningen Alt. 1
C-----------------------------------------------------------------
   21 WRITE(LPU,3200)
      RETURN
   22 WRITE(LPU,3300)
      RETURN
   23 WRITE(LPU,3400)
      RETURN
C
C-----------------------------------------------------------------
C     Feil i stivhetsberegningen / Alt. 2
C-----------------------------------------------------------------
   31 WRITE(LPU,4200)
      RETURN
   32 WRITE(LPU,4300)
      RETURN
   33 WRITE(LPU,4400)
      RETURN
   34 WRITE(LPU,4400)
      RETURN
C
C------------------------------------------------------------------
C     Feil i spenningsberegningen / Alt. 2
C------------------------------------------------------------------
   41 WRITE(LPU,5200)
      RETURN
   42 WRITE(LPU,5300)
C
   99 RETURN
C
 1000 FORMAT(/15X,'** ERROR IN STIFFNESS, QUAD w/rotational DOFs **')
 1100 FORMAT(15X,'----------------------------------------------------')
 1200 FORMAT(15X,'AREA less than or equal to zero -SM4B')
 1300 FORMAT(15X,'Degenerated quad not acceptable, coord. error? -SM4H')
 1400 FORMAT(15X,'Singular G-matrix -DANNHH')
 1500 FORMAT(15X,'WARNING: The G-matrix is nearly singular -DANNHH')
 1600 FORMAT(15X,'Error establishing transformation matrix -DANNHH')
C
 2000 FORMAT(/15X,'** ERROR IN LOADVECTOR, QUAD w/rotational DOFs **')
 2100 FORMAT(15X,'-------------------------------------------------')
 2200 FORMAT(15X,'Error computing dimensionless coords. -ELLAST')
 2300 FORMAT(15X,'Error constructing the H-matrix -ELLAST')
 2400 FORMAT(15X,'Error during matrix multiplication -ELLAST')
C
 3000 FORMAT(/15X,'** ERROR IN STRESS CALCULATION - STRS08 **')
 3100 FORMAT(15X,'------------------------------------------')
 3200 FORMAT(15X,'Error in call to DANNKH')
 3300 FORMAT(15X,'Error in call to DANNHH')
 3400 FORMAT(15X,'Error in call to SPKNP')
C
 4000 FORMAT(/15X,'** ERROR IN STIFFNESS -- Alt. II **')
 4100 FORMAT(15X,'----------------------------------------------')
 4200 FORMAT(15X,'Error computing higher order stiffness -LAGHOY')
 4300 FORMAT(15X,'Error constructing the H-matrix -LAGHH')
 4400 FORMAT(15X,'Error multiplying KQH and HH -SM4H2')
C
 5000 FORMAT(/15X,'** ERROR - STRS08 -- Alt. II **')
 5100 FORMAT(15X,'--------------------------------------------')
 5200 FORMAT(15X,'Error during construction of HH / HRC -LAGHH')
 5300 FORMAT(15X,'Error during stress calculation -KNPSP2')
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C Navn   : LAGG / F I E S T A
C
C Form}l : ] danne G-matrisen for elementet formulert i skjevt
C          lokalt system.
C
C
C Input  : X,Y,E,N,IFLAG
C
C Output : G,IFLAG
C
C
C Kaller : Ingen
C
C
C Versjon:      01      Dato: 85.10.28            Progr. av:B.Lund   
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      SUBROUTINE LAGG(G,X,Y,E,N,IFLAG)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION G(12,12),X(4),Y(4),GRC(12,6),GH(12,6),
     $                 E(4),N(4),DETJG(4),NGETA(4),NGKSI(4),NGEY,NGKY,
     $                 NGEX,NGKX,A1,A2,AREAL,XC
      INTEGER IFLAG,I,J,K
C
C
C---------------------------------------------------------------------
C     Sjekker orientering av elementet for } f} rett G-matrise.
C---------------------------------------------------------------------
C
C
      IFLAG=0
      A1=DABS(0.5*(X(2)*Y(3)-X(3)*Y(2)+X(3)*Y(1)-X(1)*Y(3)+X(1)*Y(2)-
     $        X(2)*Y(1)))
      A2=DABS(0.5*(X(3)*Y(4)-X(4)*Y(3)+X(4)*Y(1)-X(1)*Y(4)+X(1)*Y(3)-
     $        X(3)*Y(1)))
      AREAL=A1+A2
      IF(A1 .LE. 0. .OR. A2 .LE. 0.)THEN
       IFLAG=IFLAG-1
       GOTO 1000
      ENDIF
C
      XC=(A1*(X(1)+X(2)+X(3))+A2*(X(3)+X(4)+X(1)))/(3.*AREAL)
C
      IF((X(2)-XC).GT. 0. .AND. (X(1)-XC).LT. 0.)THEN
       E(1)=-1.
       E(2)=1.
       E(3)=1.
       E(4)=-1.
       N(1)=-1.
       N(2)=-1.
       N(3)=1.
       N(4)=1.
      ENDIF
      IF((X(2)-XC).GT. 0. .AND. (X(1)-XC).GT. 0.)THEN
       E(1)=1.
       E(2)=1.
       E(3)=-1.
       E(4)=-1.
       N(1)=-1.
       N(2)=1.
       N(3)=1.
       N(4)=-1.
      ENDIF
      IF((X(2)-XC).LT. 0. .AND. (X(1)-XC).GT. 0.)THEN
       E(1)=1.
       E(2)=-1.
       E(3)=-1.
       E(4)=1.
       N(1)=1.
       N(2)=1.
       N(3)=-1.
       N(4)=-1.
      ENDIF
      IF((X(2)-XC).LT. 0. .AND. (X(1)-XC).LT. 0.)THEN
       E(1)=-1.
       E(2)=-1.
       E(3)=1.
       E(4)=1.
       N(1)=1.
       N(2)=-1.
       N(3)=-1.
       N(4)=1.
      ENDIF
C
C
C---------------------------------------------------------------------
C     Gir verdier til GH og GRC hver for seg og setter disse sammen
C     til G.
C---------------------------------------------------------------------
C
C
      DO 30 I=1,12
       DO 20 J=1,6
        GH(I,J)=0.
        GRC(I,J)=0.
 20    CONTINUE
 30   CONTINUE
C
C
      K=1
C
C
C-----------------------------------------------------------------
C     Danner GRC og GH for de 4 knutepunktene.
C-----------------------------------------------------------------
C
C
      DO 50 I=1,4
C
       DO 35 J=1,4
        NGETA(J)=0.
        NGKSI(J)=0.
 35    CONTINUE
C
       NGETA(1)=0.25*N(1)*(1.+(E(1)*E(I)))
       NGETA(2)=0.25*N(2)*(1.+(E(2)*E(I)))
       NGETA(3)=0.25*N(3)*(1.+(E(3)*E(I)))
       NGETA(4)=0.25*N(4)*(1.+(E(4)*E(I)))
       NGKSI(1)=0.25*E(1)*(1.+(N(1)*N(I)))
       NGKSI(2)=0.25*E(2)*(1.+(N(2)*N(I)))
       NGKSI(3)=0.25*E(3)*(1.+(N(3)*N(I)))
       NGKSI(4)=0.25*E(4)*(1.+(N(4)*N(I)))
       NGEY=0.
       NGKY=0.
       NGEX=0.
       NGKX=0.
C
       DO 40 J=1,4
        NGEY=NGEY+NGETA(J)*Y(J)
        NGKY=NGKY+NGKSI(J)*Y(J)
        NGEX=NGEX+NGETA(J)*X(J)
        NGKX=NGKX+NGKSI(J)*X(J)
 40    CONTINUE
C
C------------------------------------------------------------------
C      Jakobideterminanten beregnes.
C------------------------------------------------------------------
C
       DETJG(I)=NGKX*NGEY-NGEX*NGKY
       IF(ABS(DETJG(I)) .LE. 1.0D-16)IFLAG=IFLAG-2
C
       GH(K,1)=E(I)*N(I)
       GH(K,2)=N(I)**2
       GH(K,5)=N(I)*((E(I)-N(I))**2)  
       GH(K,6)=-N(I)*((E(I)+N(I))**2)  
       GH(K+1,3)=E(I)*N(I)  
       GH(K+1,4)=E(I)**2  
       GH(K+1,5)=-E(I)*((N(I)-E(I))**2)
       GH(K+1,6)=E(I)*((N(I)+E(I))**2)  
       GH(K+2,1)=0.5*(NGEX*N(I)-NGKX*E(I))/DETJG(I)
       GH(K+2,2)=0.5*(-2.*NGKX*N(I))/DETJG(I)
       GH(K+2,3)=0.5*(NGEY*N(I)-NGKY*E(I))/DETJG(I)
       GH(K+2,4)=NGEY*E(I)/DETJG(I)
       GH(K+2,5)= 0.5*(NGEY*(-(N(I)**2)+4.*E(I)*N(I)-3.*(E(I)**2))
     $                -NGKY*(2.*(E(I)**2)-2.*E(I)*N(I))
     $                +NGEX*(2.*E(I)*N(I)-2.*(N(I)**2))
     $         -NGKX*((E(I)**2)-4.*E(I)*N(I)+3.*(N(I)**2)))/DETJG(I)
       GH(K+2,6)=0.5*(NGEY*((N(I)**2)+4.*E(I)*N(I)+3.*(E(I)**2))
     $               -NGKY*(2.*E(I)*N(I)+2.*(E(I)**2))
     $               +NGEX*(-2.*E(I)*N(I)-2.*(N(I)**2))
     $         -NGKX*(-(E(I)**2)-4.*E(I)*N(I)-3.*(N(I)**2)))/DETJG(I)
       GRC(K,1)=1.
       GRC(K+1,1)=0.
       GRC(K+2,1)=0.
       GRC(K,2)=0.
       GRC(K+1,2)=1.
       GRC(K+2,2)=0.
       GRC(K,3)=-N(I)
       GRC(K+1,3)=E(I)
       GRC(K+2,3)=0.5*(NGEY+NGKX)/DETJG(I)
       GRC(K,4)=E(I)
       GRC(K+1,4)=0.
       GRC(K+2,4)=0.5*NGEX/DETJG(I)
       GRC(K,5)=0.
       GRC(K+1,5)=N(I)
       GRC(K+2,5)=-0.5*NGKY/DETJG(I)
       GRC(K,6)=N(I)
       GRC(K+1,6)=E(I)
       GRC(K+2,6)=0.5*(NGEY-NGKX)/DETJG(I)
       K=K+3
 50   CONTINUE
C
C
C-----------------------------------------------------------------
C     Setter inn submatriser p} rett plass.
C-----------------------------------------------------------------
C
C
      DO 60 I=1,12
       DO 55 J=1,6
        G(I,J)=GRC(I,J)
 55    CONTINUE
 60   CONTINUE
C
      DO 80 I=1,12
       DO 70 J=7,12
        G(I,J)=GH(I,J-6)
 70    CONTINUE
 80   CONTINUE
C
 1000  R E T U R N
C
C  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **
C
C     Returnerer beregnet G-matrise og orientering av elementet.
C
C  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **
      E N D
C
      SUBROUTINE QMRF35(EM,X,Y,Z,THK,RHO)
C
C***********************************************************************
C
C     QMRF35 DETERMINES A LUMPED MASS MATRIX  EM  FOR THE QUADRILATERAL
C     MEMBRANE ELEMENT WITH ROTATIONAL DEGREES OF FREEDOM OF P.G.BERGAN
C     AND C.A.FELIPPA.
C
C     PROGRAMMED BY :  KETIL AAMNES
C     DATE/VERSION  :  081085 / 1.0
C
C***********************************************************************
C
      IMPLICIT NONE
C
      DOUBLE PRECISION  EM(24),EMM(2,9),X(4),Y(4),Z(4),THK,RHO
      DOUBLE PRECISION  AUX1,COSG,SING,SL21,SL31
      DOUBLE PRECISION  X1,X2,X3,Y1,Y2,Y3,Z1,Z2,Z3,XL(4),YL(4)
      DOUBLE PRECISION  AREA,XTP,YTP,H,B,B1,B2,S1,S2,SL
      DOUBLE PRECISION  IXX,IYY1,IYY2,IYY,IZZ1,IZZ2,IZZ
      DOUBLE PRECISION  IX,IY,IZ,M,GAMMA,TR(3,3)
      INTEGER           I
C
C ---------------------------- Local coordinates -------------------
C
      DO 10 I=1,2
C
C
      X1 = 0.
      Y1 = 0.
      Z1 = 0.
C
      X2 = X(I+1)-X(1)
      Y2 = Y(I+1)-Y(1)
      Z2 = Z(I+1)-Z(1)
C
      X3 = X(I+2)-X(1)
      Y3 = Y(I+2)-Y(1)
      Z3 = Z(I+2)-Z(1)
C
      SL21 = DSQRT(X2*X2 + Y2*Y2 + Z2*Z2)
      SL31 = DSQRT(X3*X3 + Y3*Y3 + Z3*Z3)
C
      COSG = (X3*X2 + Y3*Y2 + Z3*Z2)/(SL31*SL21)
      AUX1 = 1.-COSG*COSG
C
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
C
C ---------------------------- V1 ----------------------------------
      TR(1,1) = X2-X1
      TR(2,1) = Y2-Y1
      TR(3,1) = Z2-Z1
      SL = DSQRT(TR(1,1)**2+TR(2,1)**2+TR(3,1)**2)
      TR(1,1) = TR(1,1)/SL
      TR(2,1) = TR(2,1)/SL
      TR(3,1) = TR(3,1)/SL
C ---------------------------- V3 ----------------------------------
      TR(1,2) = X3-X1
      TR(2,2) = Y3-Y1
      TR(3,2) = Z3-Z1
      TR(1,3) = (TR(2,1)*TR(3,2)-TR(2,2)*TR(3,1))
      TR(2,3) = (TR(1,2)*TR(3,1)-TR(1,1)*TR(3,2))
      TR(3,3) = (TR(1,1)*TR(2,2)-TR(1,2)*TR(2,1))
C
      SL = DSQRT(TR(1,3)**2+TR(2,3)**2+TR(3,3)**2)
      TR(1,3) = TR(1,3)/SL
      TR(2,3) = TR(2,3)/SL
      TR(3,3) = TR(3,3)/SL
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
      EMM(I,1) = M
      EMM(I,2) = M
      EMM(I,4) = M
      EMM(I,5) = M
      EMM(I,7) = M
      EMM(I,8) = M
C
      EMM(I,3) = GAMMA*IZ
      EMM(I,6) = EMM(I,3)
      EMM(I,9) = EMM(I,3)
C
   10 CONTINUE
C
      EM( 1) = EMM(1,1) + EMM(2,1)
      EM( 2) = EMM(1,2) + EMM(2,2)
      EM( 6) = EMM(1,3) + EMM(2,3)
C
      EM( 7) = EMM(1,4)
      EM( 8) = EMM(1,5)
      EM(12) = EMM(1,6)
C
      EM(13) = EMM(1,7) + EMM(2,4)
      EM(14) = EMM(1,8) + EMM(2,5)
      EM(18) = EMM(1,9) + EMM(2,6)
C
      EM(19) = EMM(2,7)
      EM(20) = EMM(2,8)
      EM(24) = EMM(2,9)
C
C---- Fordeler masser og treghetsmomenter likt pa de 4 knutepunktene
C
      EM(1) = (EM(1) + EM(7) + EM(13) + EM(19))/4.0
      EM(2) = (EM(2) + EM(8) + EM(14) + EM(20))/4.0
      EM(6) = (EM(6) + EM(12) + EM(18) + EM(24))/4.0
C
      EM(7)  = EM(1)
      EM(13) = EM(1)
      EM(19) = EM(1)
C
      EM(8)  = EM(2)
      EM(14) = EM(2)
      EM(20) = EM(2)
C
      EM(12) = EM(6)
      EM(18) = EM(6)
      EM(24) = EM(6)
C
C ---------------------------- Return ------------------------------
   80 RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Navn   : QMRF32 / F E D E M
C
C Form}l : ] beregne de opptredene spenninger i fri formulerings-
C          elementet. Elementet er firkantet og har 3 frihetsgrd.
C          pr. knutepunkt. Spenningene beregnes i de 4 knutepktene
C          samt i flatesenteret. Den 4.de spenningskomponenten
C          beregnes desom elementet er i plan t|yningstilstand.
C
C
C Input  : V,XL,YL,E,RNY,TYKK,NS
C
C Output : S
C
C
C Kaller : Subrutinen DANNHH
C          Subrutinen SM4H
C          Subrutinen SPKNP
C          Subrutinen SKFEIL
C
C Versjon:   02     Dato: 85.11.04               Progr. av:B.Lund
C Versjon:   03     Dato: 88.08.20               Progr. av:T.R.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      SUBROUTINE QMRF32(V,XL,YL,EM,RNY,TYKK,NS,S,LPU,IFLAG)
C
      IMPLICIT NONE
C
      INTEGER          NS,LPU,IFLAG,K
      DOUBLE PRECISION V(12),XL(4),YL(4),EM(3,3),RNY,TYKK,S(NS,5)
      DOUBLE PRECISION GAMMA,A1,A2,KH(12,12),E(4),N(4),HH(6,12)
      DOUBLE PRECISION HRC(6,12),TETAX,T(12,12),SPPUNK(2,5),DIMXL,DIMYL
C
C----------------------------------------------------------------------
C     Beregn HH og HRC. Legg spenningspunktenes koord. i SPPUNK.
C----------------------------------------------------------------------
C
      IFLAG=0
      GAMMA=TYKK*EM(1,1)
      CALL SM4H(KH,XL,YL,DIMXL,DIMYL,GAMMA,RNY,E,N,TETAX,A1,A2,IFLAG)
      IF(IFLAG .LT. 0)THEN
         IFLAG=-10
         GOTO 915
      ENDIF
C
      CALL DANNHH(HH,HRC,E,N,DIMXL,DIMYL,TETAX,T,IFLAG)
      IF(IFLAG .LT. 0)THEN
         IFLAG=-20
         GOTO 915
      ENDIF
C
      DO 10 K=1,4
         SPPUNK(1,K)=E(K)
         SPPUNK(2,K)=N(K)
 10   CONTINUE
      SPPUNK(1,5)=0.0D0
      SPPUNK(2,5)=0.0D0
C
C--------------------------------------------------------------------
C     Beregn spenningene i ett punkt av gangen.
C     Adder spenninger inn i spenningsvektoren.
C--------------------------------------------------------------------
C
      DO 30 K=1,5
C
         CALL SPKNP(S(1,K),SPPUNK(1,K),SPPUNK(2,K),HH,HRC,V,EM,
     +              DIMXL,DIMYL,TETAX)
C
         IF (NS .EQ. 4) S(4,K) = RNY*(S(1,K)+S(2,K))
C
 30   CONTINUE
C
      GOTO 999
C
 915  CALL SKFEIL(IFLAG,LPU)
C
 999  R E T U R N
C
C  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **
C
C     Returnerer elementets spenningsvektor.
C
C  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **
C
      E N D
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C Navn   : SPKNP / F E D E M
C
C Form}l : ] beregne spenningene i knutepunktene og i flatesenteret.
C
C
C Input  : A,B,HH,HRC,V,E,DIMX,DIMY
C
C Output : SL
C
C
C Kaller : Subrutinen MAB / S A M
C
C
C Versjon:  01        Dato: 85.10.22              Progr. av:B.Lund
C Versjon:  02        Dato: 88.08.20              Progr. av:T.R.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      SUBROUTINE SPKNP(SL,A,B,HH,HRC,V,E,DIMX,DIMY,TETAX)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION SL(4),A,B,HH(6,12),HRC(6,12),V(12),E(3,3),
     $                 DIMX,DIMY,TETAX
      DOUBLE PRECISION KNKOOR(2,1),BRC(3,6),BH(3,6),
     $                 BRCHRC(3,12),BHHH(3,12),ERC(3,1),EH(3,1),
     $                 EPS(3,1),TREPS(3,1),FACT2,VL(12)
      INTEGER          I,J
C
C-------------------------------------------------------------------
C     Gir verdier til BRC og BH.
C-------------------------------------------------------------------
C
      KNKOOR(1,1)=A
      KNKOOR(2,1)=B
C
      DO 20 I=1,3
         DO 10 J=1,6
            BRC(I,J)=0.
            BH(I,J)=0.
 10      CONTINUE
 20   CONTINUE
C
C----------------------------------------------------------
C     Forskyvn. i lokalt system.
C----------------------------------------------------------
C
      VL(1)=V(1)*COS(TETAX)+V(2)*SIN(TETAX)
      VL(2)=-V(1)*SIN(TETAX)+V(2)*COS(TETAX)
      VL(3)=V(3)
      VL(4)=V(4)*COS(TETAX)+V(5)*SIN(TETAX)
      VL(5)=-V(4)*SIN(TETAX)+V(5)*COS(TETAX)
      VL(6)=V(6)
      VL(7)=V(7)*COS(TETAX)+V(8)*SIN(TETAX)
      VL(8)=-V(7)*SIN(TETAX)+V(8)*COS(TETAX)
      VL(9)=V(9)
      VL(10)=V(10)*COS(TETAX)+V(11)*SIN(TETAX)
      VL(11)=-V(10)*SIN(TETAX)+V(11)*COS(TETAX)
      VL(12)=V(12)
C
      BRC(1,4)=DIMX
      BRC(2,5)=DIMY
      BRC(3,3)=DIMX-DIMY
      BRC(3,6)=DIMX+DIMY
      BH(1,1)=DIMX*KNKOOR(2,1)
      BH(1,5)=DIMX*(2.*KNKOOR(1,1)*KNKOOR(2,1)-2.*(KNKOOR(2,1)**2))
      BH(1,6)=-2.*DIMX*(KNKOOR(1,1)*KNKOOR(2,1)+(KNKOOR(2,1)**2))
      BH(2,3)=DIMY*KNKOOR(1,1)
      BH(2,5)=DIMY*(2.*(KNKOOR(1,1)**2)-2.*KNKOOR(1,1)*KNKOOR(2,1))
      BH(2,6)=DIMY*(2.*KNKOOR(1,1)*KNKOOR(2,1)+2.*(KNKOOR(1,1)**2))
      BH(3,1)=DIMY*KNKOOR(1,1)
      BH(3,2)=2.*DIMY*KNKOOR(2,1)
      BH(3,3)=DIMX*KNKOOR(2,1)
      BH(3,4)=2.*DIMX*KNKOOR(1,1)
      BH(3,5)=DIMY*((KNKOOR(1,1)**2)-4.*KNKOOR(1,1)*KNKOOR(2,1)+3.*
     $              (KNKOOR(2,1)**2))+DIMX*(-(KNKOOR(2,1)**2)+4.*
     $               KNKOOR(1,1)*KNKOOR(2,1)-3.*(KNKOOR(1,1)**2))
      BH(3,6)=-DIMY*((KNKOOR(1,1)**2)+4.*KNKOOR(1,1)*KNKOOR(2,1)+3.*
     $               (KNKOOR(2,1)**2))+DIMX*((KNKOOR(2,1)**2)+4.*
     $                KNKOOR(1,1)*KNKOOR(2,1)+3.*(KNKOOR(1,1)**2))
C
C-----------------------------------------------------------------
C     Danner t|yninger og spenninger ved matrisemult.
C-----------------------------------------------------------------
C
      CALL MAB(BRC,HRC,BRCHRC,3,6,12,0)
      CALL MAB(BH,HH,BHHH,3,6,12,0)
      CALL MAB(BRCHRC,VL,ERC,3,12,1,0)
      CALL MAB(BHHH,VL,EH,3,12,1,0)
C
      DO 30 I=1,3
         EPS(I,1)=ERC(I,1)+EH(I,1)
 30   CONTINUE
C
C----------------------------------------------------------------------
C     Transformerer spenninger - "Mohrs"-sirkel.
C---------------------------------------------------------------------
C
      FACT2=COS(2.*TETAX)
C
      IF(ABS(FACT2) .GT. 1.0E-16)THEN
         TREPS(1,1)=(EPS(1,1)*(1.+COS(2.*TETAX))+EPS(2,1)*
     $              (1.-COS(2.*TETAX))-EPS(3,1)*SIN(2.*TETAX))*0.5
         TREPS(2,1)=EPS(1,1)+EPS(2,1)-TREPS(1,1)
         TREPS(3,1)=((EPS(1,1)-EPS(2,1))*SIN(2.*TETAX)+EPS(3,1))/FACT2
      ENDIF
C
      CALL MAB(E,TREPS,SL,3,3,1,0)
C
      RETURN
C
C **  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **
C
C     Returnerer spenningsvektor i aktuelt punkt.
C
C **  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **  **
C
      END
