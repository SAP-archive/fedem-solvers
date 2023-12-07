C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
C     This file is part of FEDEM - https://openfedem.org
C
      SUBROUTINE JACO31 (JA,DNXI,DNET,DNZE,XG,YG,ZG,MEK)
C
C***********************************************************************
C
C   F I N I T E  E L E M E N T  L I B R A R Y  SUBROUTINE :  JACO31
C
C     JACO31 CALCULATES THE JACOBIAN MATRIX
C     FOR AN ISOPARAMETRIC 3-DIMENSIONAL ELEMENT WITH MEK NODES.
C
C     PROGRAMMED BY: B. AAMODT AND H.F. KLEM
C
C     DATE/VERSION : 22.1.74/01
C
C***********************************************************************
C
C   INPUT:
C     DNXI  \
C     DNET  - THE DERIVATIVES OF THE SHAPE FUNCTION IN THREE DIRECTIONS.
C     DNZE  /
C     XG    \
C     YG    -  GLOBAL NODAL COORDINATES DEFINING THE GEOMETRY.
C     ZG    /
C     MEK   -  NUMBER OF NODES ON THE THREE DIMENSIONAL ELEMENT.
C
C   OUTPUT:
C     JA    -  THE JACOBIAN MATRIX
C***********************************************************************
C
      IMPLICIT NONE
C
      INTEGER            MEK, I, J
      DOUBLE PRECISION   JA(3,3), DNXI(MEK), DNET(MEK), DNZE(MEK),
     1                   XG(MEK), YG(MEK), ZG(MEK)
C
      DO 100 I=1,3
      DO 100 J=1,3
  100 JA(I,J) = 0.0D0
C
      DO 130 I=1,MEK
      JA(1,1) = JA(1,1) + DNXI(I)*XG(I)
      JA(1,2) = JA(1,2) + DNXI(I)*YG(I)
      JA(1,3) = JA(1,3) + DNXI(I)*ZG(I)
      JA(2,1) = JA(2,1) + DNET(I)*XG(I)
      JA(2,2) = JA(2,2) + DNET(I)*YG(I)
      JA(2,3) = JA(2,3) + DNET(I)*ZG(I)
      JA(3,1) = JA(3,1) + DNZE(I)*XG(I)
      JA(3,2) = JA(3,2) + DNZE(I)*YG(I)
  130 JA(3,3) = JA(3,3) + DNZE(I)*ZG(I)
C
      RETURN
      END
      SUBROUTINE JACI31 (JI,DETJ,DNXI,DNET,DNZE,XG,YG,ZG,IEL,MEK,
     1                   IW,IPSW,IERR)
C
C***********************************************************************
C
C   F I N I T E  E L E M E N T  L I B R A R Y  SUBROUTINE :  JACI31
C
C     JACI31 CALCULATES THE DETERMINANT AND INVERSE OF THE JACOBIAN
C     MATRIX FOR AN ISOPARAMETRIC 3-DIMENSIONAL ELEMENT WITH MEK NODES.
C
C     PROGRAMMED BY: B. AAMODT AND H.F. KLEM
C
C     DATE/VERSION : 22.1.74/01
C
C***********************************************************************
C
C   INPUT:
C     DNXI  \
C     DNET  - THE DERIVATIVES OF THE SHAPE FUNCTION IN THREE DIRECTIONS.
C     DNZE  /
C     XG    \
C     YG    -  GLOBAL NODAL COORDINATES DEFINING THE GEOMETRY.
C     ZG    /
C     IEL   -  ELEMENT IDENTIFICATION (E.G., ELEMENT NUMBER).
C     MEK   -  NUMBER OF NODES ON THE THREE DIMENSIONAL ELEMENT.
C     IW    -  WRITE UNIT.
C     IPSW  -  PRINT KEY. IF > 2 PRINT OF CONTROL DATA IS GIVEN.
C              IF > 3 PRINT OF RESULTS IS GIVEN IN ADDITION.
C
C   OUTPUT:
C     IERR  -  ERROR FLAG. IERR=0: NO ERRORS DETECTED. IERR=-1: THE
C              DETERMINANT OF THE JACOBIAN MATRIX IS NUMERICALLY ZERO.
C     JI    -  THE INVERTED JACOBIAN MATRIX
C     DETJ  -  THE DETERMINANT OF THE JACOBIAN MATRIX
C***********************************************************************
C
      IMPLICIT NONE
C
      INTEGER            IEL, MEK, IW, IPSW, IERR
      DOUBLE PRECISION   JI(3,3), DETJ, DNXI(MEK), DNET(MEK), DNZE(MEK),
     1                   XG(MEK), YG(MEK), ZG(MEK)
C
      INTEGER            I1, J1
      DOUBLE PRECISION   EPS, J(3,3)
      PARAMETER        ( EPS = tiny(1.0D0)*100.0D0 )
C
      IERR=0
C
C     PRINT OF INPUT ?
C
      IF(IPSW-3) 110,50,50
   50 CONTINUE
      WRITE(IW,6000) IEL
      WRITE(IW,6010)
      DO 70 I1=1,MEK
      WRITE(IW,6020) I1,XG(I1),YG(I1),ZG(I1)
   70 CONTINUE
      WRITE(IW,6030)
      I1=1
      CALL MPRT30(DNXI,MEK,I1,IW)
      WRITE(IW,6040)
      CALL MPRT30(DNET,MEK,I1,IW)
      WRITE(IW,6050)
      CALL MPRT30(DNZE,MEK,I1,IW)
      IF (IERR) 190,110,110
C
C     THE JACOBIAN MATRIX
C
  110 CONTINUE
      CALL JACO31(J,DNXI,DNET,DNZE,XG,YG,ZG,MEK)
C
C     THE DETERMINANT OF THE JACOBIAN MATRIX
C
      DETJ=J(1,1)*(J(2,2)*J(3,3)-J(2,3)*J(3,2))
     1    +J(1,2)*(J(2,3)*J(3,1)-J(2,1)*J(3,3))
     2    +J(1,3)*(J(2,1)*J(3,2)-J(2,2)*J(3,1))
      IF (ABS(DETJ)-EPS) 170,170,150
C
C     THE INVERSE OF THE JACOBIAN MATRIX
C
  150 CONTINUE
      JI(1,1)=(J(2,2)*J(3,3)-J(2,3)*J(3,2))/DETJ
      JI(1,2)=(J(3,2)*J(1,3)-J(3,3)*J(1,2))/DETJ
      JI(1,3)=(J(1,2)*J(2,3)-J(1,3)*J(2,2))/DETJ
      JI(2,1)=(J(3,1)*J(2,3)-J(3,3)*J(2,1))/DETJ
      JI(2,2)=(J(1,1)*J(3,3)-J(1,3)*J(3,1))/DETJ
      JI(2,3)=(J(2,1)*J(1,3)-J(2,3)*J(1,1))/DETJ
      JI(3,1)=(J(2,1)*J(3,2)-J(2,2)*J(3,1))/DETJ
      JI(3,2)=(J(3,1)*J(1,2)-J(3,2)*J(1,1))/DETJ
      JI(3,3)=(J(1,1)*J(2,2)-J(1,2)*J(2,1))/DETJ
C
      DETJ=ABS(DETJ)
C
      GO TO 210
C
  170 CONTINUE
      IERR=-1
      IF(IPSW-3) 50,190,190
C
C     PRINT OF OUTPUT IF IN ERROR
C
  190 CONTINUE
      WRITE(IW,6500)
      WRITE(IW,6505) ((J(I1,J1),J1=1,3),I1=1,3)
      WRITE(IW,6510) IERR
      WRITE(IW,6520) IEL
      GO TO 240
C
C     PRINT OF OUTPUT ?
C
  210 CONTINUE
      IF (IPSW-3) 250,240,230
  230 WRITE(IW,6550)
      I1=3
      CALL MPRT30(J,I1,I1,IW)
      WRITE(IW,6560)
      CALL MPRT30(JI,I1,I1,IW)
  240 WRITE(IW,6530) DETJ
      WRITE(IW,6540) IEL
  250 RETURN
C
 6000 FORMAT(///2X,30H*ENTERING*JACI31  FOR ELEMENT ,I7,2H :)
 6010 FORMAT(//,'     NODAL POINT COORDINATES:',
     2        /,'     NODE        XG          YG          ZG')
 6020 FORMAT(I9,6X,1P3E12.4)
 6030 FORMAT(/5X,75HTHE DERIVATIVES OF THE SHAPE FUNCTION WITH RESPECT T
     2O XI (DNXI) :                     )
 6040 FORMAT(/5X,75HTHE DERIVATIVES OF THE SHAPE FUNCTION WITH RESPECT T
     2O ETTA (DNET) :                   )
 6050 FORMAT(/5X,75HTHE DERIVATIVES OF THE SHAPE FUNCTION WITH RESPECT T
     2O ZETA (DNZE) :                   )
 6500 FORMAT(///71H *** ERROR RETURN FROM JACI31 *** - THE DETERMINANT I
     2S NUMERICALLY ZERO )
 6505 FORMAT(/25H     THE JACOBIAN MATRIX:/(5X,1P3E12.4))
 6510 FORMAT(/16H     ERROR FLAG=,I3)
 6520 FORMAT(/13H     ELEMENT ,I7,11H   IN ERROR)
 6530 FORMAT(/44H     THE DETERMINANT HAS THE VALUE (DETJ):  ,1PE12.4)
 6540 FORMAT(//2X,29H*LEAVING*JACI31* FOR ELEMENT:,I7)
 6550 FORMAT(/33H     PRINT OF THE JACOBI MATRIX :)
 6560 FORMAT(/42H     PRINT OF THE INVERTED JACOBI MATRIX: )
C
      END
