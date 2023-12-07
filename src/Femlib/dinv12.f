C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
C     This file is part of FEDEM - https://openfedem.org
C
C>    @brief Inverts a NxN square matrix using LAPACK.
C
      SUBROUTINE DINV12 (N,A,WORK,IWORK,IWR,IERR)
      INTEGER            N,IWR,IERR
      DOUBLE PRECISION   A(N,N),WORK(N),IWORK(N)
C
      CALL DGETRF (N,N,A,N,IWORK,IERR)
      IF (IERR .NE. 0) THEN
         WRITE(IWR,"(' *** LAPACK::DGETRF, INFO =',I8)") IERR
         RETURN
      END IF
C
      CALL DGETRI (N,A,N,IWORK,WORK,N,IERR)
      IF (IERR .NE. 0) THEN
         WRITE(IWR,"(' *** LAPACK::DGETRI, INFO =',I8)") IERR
         RETURN
      END IF
C
      RETURN
      END
