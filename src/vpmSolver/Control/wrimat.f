C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
C     This file is part of FEDEM - https://openfedem.org
C
      subroutine WRIMAT (A,M,N,IO,CNAME)
      integer            M,N,IO
      double precision   A(M,N)
      character*(*)      CNAME
      integer            NELL
      parameter        ( NELL = 6 )
C
      write(IO,6000) CNAME,min(M,N)
      do J = 1, N, NELL
         L = min(J+NELL-1,N)
         write(IO,6010) (K,K=J,L)
         do I = 1, M
            write(IO,6020) I,(A(I,K),K=J,L)
         end do
      end do
C
      return
C
 6000 format(//' >>>>> Dump of Ctrl. Sys. ',A/7X,'NDIM =',I8)
 6010 format(/20I13)
 6020 format(I6,1P20E13.5)
C
      end
