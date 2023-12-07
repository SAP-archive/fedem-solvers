C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
C     This file is part of FEDEM - https://openfedem.org
C
      SUBROUTINE XERBLA (SRNAME,INFO)
      CHARACTER*6        SRNAME
      INTEGER            INFO
      WRITE(*,9999) SRNAME,INFO
      RETURN
 9999 FORMAT('*** ',A,': Parameter number',I3,' had an illegal value')
      END
