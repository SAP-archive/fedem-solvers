C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
C     This file is part of FEDEM - https://openfedem.org
C
      SUBROUTINE READMAT (MATNAS,NCOL,NTOT,IUN,LPU,IERR)
C
C***********************************************************************
C
C     THIS SUBROUTINE READS A NASTRAN OUTPUT2 MATRIX FILE.
C     DATA BLOCK NAME IS 2 WORDS, NASTRAN TRAILER IS 7 WORDS,
C     AND THE OUTPUT2 LABEL IS 2 WORDS.
C
      IMPLICIT NONE
C
      INTEGER NCOL, NTOT, IUN, LPU, IERR
C
C     DIMENSION ALLOWS AS MANY AS 10000 TERMS, IT CAN BE INCREASED
      INTEGER     MAXCOL
      PARAMETER ( MAXCOL = 10000 )
C
      DOUBLE PRECISION COLUMN(MAXCOL), MATNAS(NTOT,*)
C
      INTEGER KEY, DATE(3), T(7), LABEL(2), NAME(2), H(2)
      INTEGER I, INO, ICOL, IROW
C
C***********************************************************************
C
C     PROCESS THE LABEL
C
      READ(IUN,ERR=910,END=920) KEY
      IF (KEY .EQ. 2) THEN
C        The header is skipped
         REWIND(IUN)
      ELSE
         READ(IUN,ERR=910,END=920) KEY
         IF (KEY .NE. 3) GO TO 930
C        Read the date
         READ(IUN,ERR=910,END=920) DATE
         READ(IUN,ERR=910,END=920) KEY
         IF (KEY .NE. 7) GO TO 930
C        Read the title
         READ(IUN,ERR=910,END=920) T
         READ(IUN,ERR=910,END=920) KEY
         IF (KEY .NE. 2) GO TO 930
C        Read the label
         READ(IUN,ERR=910,END=920) LABEL
         READ(IUN,ERR=910,END=920) KEY
         IF (KEY .NE.-1) GO TO 930
         READ(IUN,ERR=910,END=920) KEY
         IF (KEY .NE. 0) GO TO 930
      END IF
C
C     PROCESS THE HEADER RECORD
C
      READ(IUN,ERR=910,END=920) KEY
      IF (KEY .NE. 2) GO TO 930
C     Read data block name
      READ(IUN,ERR=910,END=920) NAME
      READ(IUN,ERR=910,END=920) KEY
      IF (KEY .NE.-1) GO TO 930
      READ(IUN,ERR=910,END=920) KEY
      IF (KEY .NE. 7) GO TO 930
C     Read trailer
      READ(IUN,ERR=910,END=920) T
      READ(IUN,ERR=910,END=920) KEY
      IF (KEY .NE.-2) GO TO 930
      READ(IUN,ERR=910,END=920) KEY
      IF (KEY .NE. 1) GO TO 930
      READ(IUN,ERR=910,END=920) KEY
      IF (KEY .NE. 0) GO TO 930
      READ(IUN,ERR=910,END=920) KEY
      IF (KEY .LT. 2) GO TO 930
C     Read data from data block header
      READ(IUN,ERR=910,END=920) H
      READ(IUN,ERR=910,END=920) KEY
      IF (KEY .NE.-3) GO TO 930
C
      WRITE(LPU,6) NAME,T
    6 FORMAT(' DATABLOCK  ',2A4 /' WITH TRAILER,',7I8,', STARTED')
C
C     GET THE DATA BLOCK
C
c read record 20 = key
      read(iun,ERR=910,END=920) key
c read record 21 = key
      read(iun,ERR=910,END=920) key
      if (key.eq.0) go to 810
c
c loop on columns
c
      do 300 icol = 1, ncol
         read(IUN,ERR=910,END=920) key
c
c key = number of terms to be read
c divide by 2 for double precision(this example is for double precision)
c
  100    ino = key/2
         if (ino .gt. MAXCOL) then
            write(LPU,*) '*** Too many data records:',ino,MAXCOL
            ierr = -ino
            return
         end if
         read(IUN,ERR=910,END=920) irow,(column(i),i=1,ino)
         if (irow+ino-1 .gt. NTOT) then
            write(LPU,*) '*** Too many matrix rows:',irow+ino-1,NTOT
            ierr = 1-IROW-ino
            return
         end if
c
c store the record data in a matrix
         do 200 i = 1, ino
            matnas(irow,icol) = column(i)
            irow = irow + 1
  200    continue
c
c read key again - if <0, then => eor
         read(IUN,ERR=910,END=920) key
c if key>0, then it is the number of words to be read
         if (key .gt. 0) go to 100
         if (key .eq. 0) go to 800
c end of column
c
c read key = 1 => start new column
         read(IUN,ERR=910,END=920) key
         if (key .ne. 1) go to 900
c
c read column number
         read(IUN,ERR=910,END=920) key
c if key = 0 then done
         if (key.eq.0) go to 800
  300 continue
      write(LPU,*) '*** Too many matrix columns:',ICOL,NCOL
      ierr = -ICOL
      return
c
  800 continue
      write(LPU,*) 'Matrix read'
      return
  810 continue
      write(LPU,*) '*** Data type (record 21) indicates this is a table'
      return
  900 continue
      write(LPU,*) '*** READMAT: KEY =',KEY,' when it should be 1'
      return
  910 WRITE(LPU,*) '*** READMAT: ERROR READING OP2-FILE, KEY =',KEY
      IERR = -1
      return
  920 WRITE(LPU,*) '*** READMAT: UNEXPECTED END-OF-FILE, KEY =',KEY
      IERR = -2
      return
  930 WRITE(LPU,*) '*** READMAT: BAD KEY =',KEY
      IERR = -3
      return
      end
