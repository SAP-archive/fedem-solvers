C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
C     This file is part of FEDEM - https://openfedem.org
C
      subroutine DGETOL (A,LDA,N,IPVT,RELTOL,CNAME,LPRINT,IO,IERR)
      integer            LDA,N,IPVT(N),IO,IERR(2)
      logical            LPRINT
      double precision   A(LDA,N),RELTOL
      character*(*)      CNAME
      double precision   ABSTOL,U,UMIN
C
C     DGETOL is a wrapper of the LAPack subroutine DGETRF,
C     which in addition checks for small pivot elements (instead of
C     only exactly zero), indicating a ill-conditioned matrix.
C     RELTOL, the relative tolerance factor, is added in the call,
C     along with parameters for debug and error message print out.
C     The actual tolerance is ABSTOL = RELTOL*minDiagElm,
C     where minDiagElm is the largest diagonal element of U,
C     the upper triangular factorized matrix.
C
      if (LPRINT) call WRIMAT (A,LDA,N,IO,CNAME)
C
      I = IDAMAX(N,A,LDA+1)
      ABSTOL = RELTOL*A(I,I)
      call DGETRF (LDA,N,A,LDA,IPVT,IERR(2))
      if (IERR(2) .lt. 0) then
         IERR(1) = -2
         write(IO,6000) IERR(2)
      else if (IERR(2) .gt. 0) then
         IERR(1) = -2
         write(IO,6100) CNAME,IERR(2)
      else if (RELTOL .gt. 0.0D0) then
         UMIN = ABSTOL + 1.0D0
         NBAD = 0
         do I = 1, N
            U = DABS(A(I,I))
            if (U .le. ABSTOL) then
               if (U .lt. UMIN) THEN
                  UMIN = U
                  IERR(2) = I
               end if
               NBAD = NBAD + 1
            end if
         end do
         if (NBAD .GT. 0) then
            IERR(1) = -2
            write(IO,6200) CNAME,RELTOL,ABSTOL/RELTOL,ABSTOL
            write(IO,6210) UMIN,IERR(2),NBAD
         end if
      end if
C
      return
C
 6000 format(' *** LAPACK::DGETRF, INFO =',I3)
 6100 format(//' >>>>> Singular ',A
     +       /7X,'First zero pivot occured in equation',I8)
 6200 format(//' >>>>> Ill-conditioned ',A
     +       /7X,'Relative tolerance factor       =',1PE10.3
     +       /7X,'Minimum diagonal value of A     =',1PE10.3
     +       /7X,'Resulting tolerance             =',1PE10.3)
 6210 format( 7X,'Value of smallest element found =',1PE10.3,
     +        2X,'(equation',I8,' )'
     +       /7X,'Number of small diagonal elements: ',I8)
C
      end
