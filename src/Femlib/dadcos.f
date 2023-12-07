C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
C     This file is part of FEDEM - https://openfedem.org
C
      SUBROUTINE DADCOS (V1,V2,SVXV,DCOSI)
C
C***********************************************************************
C
C  F I N I T E  E L E M E N T  L I B R A R Y  SUBROUTINE :  DADCOS
C
C     DADCOS CROSS MULTIPLIES A VECTOR V1 WITH A VECTOR V2 (V1 X V2) AND
C     STORES THE RESULT IN VXV. THE DIMENSION OF V1, V2 AND VXV IS 3.
C     FURTHER ON IT CALCULATES THE LENGTH (SVXV) OF THE VECTOR VXV AND
C     IT'S DIRECTION COSINES (DCOSI) BY NORMALIZING VXV.
C
C     PROGRAMMED BY: H.F. KLEM
C
C     DATE/VERSION : 15.10.74/01
C
C***********************************************************************
C
C   INPUT:
C     V1    - VECTOR WITH TERMS IN THE X, Y AND Z-DIRECTION.
C     V2    - VECTOR WITH TERMS IN THE X, Y AND Z-DIRECTION.
C
C   OUTPUT:
C     SVXV  - THE LENGTH OF THE VECTOR VXV = V1 x V2.
C     DCOSI - THE VECTOR VXV NORMALIZED, I.E., THE DIRECTION COSINES OF
C             THE VECTOR IN SPACE.
C***********************************************************************
C
      IMPLICIT NONE
C
      DOUBLE PRECISION   V1(3), V2(3), SVXV, DCOSI(3), VXV(3)
C
      VXV(1) = V1(2)*V2(3) - V1(3)*V2(2)
      VXV(2) = V1(3)*V2(1) - V1(1)*V2(3)
      VXV(3) = V1(1)*V2(2) - V1(2)*V2(1)
C
      SVXV = SQRT(VXV(1)*VXV(1) + VXV(2)*VXV(2) + VXV(3)*VXV(3))
C
      DCOSI(1) = VXV(1) / SVXV
      DCOSI(2) = VXV(2) / SVXV
      DCOSI(3) = VXV(3) / SVXV
C
      RETURN
      END
