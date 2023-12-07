C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
C     This file is part of FEDEM - https://openfedem.org
C
      SUBROUTINE CONVERT_CMS (IERR)
C
C***********************************************************************
C  *
C  * MODUL            : CONVERT_CMS
C  *
C***********************************************************************
C  *
C  * PROGRAM CATEGORY : FORTRAN SUBROUTINE
C  *
C  * PURPOSE          : THE PURPOSE OF THIS PROGRAM IS TO READ OUTPUT2
C  *                    FILES CONTAINING MATRIX M, K AND G OUTPUT FROM
C  *                    NASTRAN CMS ANALYSIS AND TRANSFORM THEM TO FEDEM
C  *                    COMPATIBLE REDUCER FILES.
C  *
C  * EXTERNAL CALLS   : NONE
C  *
C  * AUTHOR           : TERJE ROLVAG
C  *
C  * DATE / VERSION   : 2001-06-28 / 1.0
C  *                    2005-06-13 / 2.0  K. M. Okstad (using F90 utils)
C  *
C***********************************************************************
C
      use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getint
      use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getstring
      use FFaFilePathInterface  , only : ffa_getbasename, ffa_checkpath
C
      IMPLICIT NONE
C
C     DEFINE VARIABLES
C
      INTEGER, intent(out) :: IERR
C
C     IOP2 IS THE FORTRAN UNIT OUTPUT2 MASS, STIFFNESS AND GRAVITY
C     DATA WILL BE READ FROM.
C
      INTEGER    STDOUT  , IOP2   , ILOG   , NFILES
      PARAMETER (STDOUT=6, IOP2=71, ILOG=12, NFILES=3)
C
C     File names: OP2 INPUT file and FEDEM reducer files
      CHARACTER(len=256) :: PARTNAME, OP2F
      CHARACTER(len=72)  :: FMXNAME, OUTF
      CHARACTER(len=32)  :: LOGF, FILETAG
C
      INTEGER                       :: I, NCOL, NTOT
      DOUBLE PRECISION, allocatable :: MATNAS(:,:)
C
C     Get input from command-line arguments
      call ffa_cmdlinearg_getint ('ndim',NTOT)
      call ffa_cmdlinearg_getstring ('partName',PARTNAME)
      call ffa_getbasename (PARTNAME,FMXNAME)
      call ffa_checkpath (PARTNAME)
C
      allocate(MATNAS(NTOT,NTOT),STAT=IERR)
      IF (IERR .NE. 0) GOTO 999
C
      LOGF = 'fedem_op2fmx.log'
      OPEN(ILOG,FILE=LOGF,STATUS='UNKNOWN',ERR=999,IOSTAT=IERR)
C
      WRITE(STDOUT,6020) '       CONVERT_CMS STARTING       '
 6020 FORMAT(/11X,10('='),'> ',A,' <',10('='))
C
C --- 3 DIFFERENT MATRICES/FILES ARE CONVERTED (M, K AND G)
C
      DO 10 I = 1, NFILES
         MATNAS = 0.0D0
         IF (I .EQ. 1) THEN
            OP2F = trim(PARTNAME)//'_M.op2'
            OUTF = trim(FMXNAME)//'_M.fmx'
            NCOL = NTOT
            FILETAG = '#FEDEM mass matrix'
         ELSE IF (I .EQ. 2) THEN
            OP2F = trim(PARTNAME)//'_S.op2'
            OUTF = trim(FMXNAME)//'_S.fmx'
            NCOL = NTOT
            FILETAG = '#FEDEM stiffness matrix'
         ELSE IF (I .EQ. 3) THEN
            OP2F = trim(PARTNAME)//'_G.op2'
            OUTF = trim(FMXNAME)//'_G.fmx'
            NCOL = 3
            FILETAG = '#FEDEM gravity force vectors'
         END IF
C
         WRITE(STDOUT,6000) OUTF
 6000    FORMAT(/2X,'Convert and write CMS reduced matrices to : ',A40)
C
C ------ OPEN UNFORMATTED NASTRAN OP2 FILE
C
         OPEN(IOP2,FILE=OP2F,STATUS='OLD',FORM='UNFORMATTED',ERR=915)
C
C ------ READ MATRIX DATA FROM NASTRAN OUTPUT 2 (OP2) FILE
C
         CALL READMAT(MATNAS,NCOL,NTOT,IOP2,ILOG,IERR)
         IF (IERR .NE. 0) GOTO 915
C
         CLOSE(IOP2)
C
C ------ MODIFY EQUATION NUMBERS (PUT COMPONENT MODES LAST IN EQ SET)
C        AND WRITE MATRIX DATA TO FEDEM FILE
C
         CALL WRIMAT(MATNAS,NCOL,NTOT,OUTF,FILETAG,ILOG,IERR)
         IF (IERR .NE. 0) GOTO 916
C
   10 CONTINUE
C
C --- SUCCESSFULL COMPLETION
C
      WRITE(STDOUT,6020) 'CONVERT_CMS SUCCESSFULLY COMPLETED'
      deallocate(MATNAS)
      CLOSE(ILOG)
      RETURN
C
C --- UNSUCCESSFULL COMPLETION
C
  915 WRITE(STDOUT,"(2X,2A)") 'Failed to open/read OP2 file : ',OP2F
      GOTO 999
  916 WRITE(STDOUT,"(2X,2A)") 'Failed to write FMX file : ',OUTF
      GOTO 999
  999 WRITE(STDOUT,6020) '    CONVERT_CMS NOT COMPLETED     '
C
      RETURN
      END
C
C **********************************************************************
C **********************************************************************
C **********************************************************************
C
      SUBROUTINE WRIMAT (MATNAS,NCOL,NTOT,OUTF,FILETAG,ILOG,IERR)
C
C***********************************************************************
C  *
C  * MODUL            : WRIMAT
C  *
C***********************************************************************
C  *
C  * PROGRAM CATEGORY : FORTRAN ROUTINE
C  *
C  * PURPOSE          : THE PURPOSE OF THIS PROGRAM IS TO WRITE FILES
C  *                    CONTAINING MATRIX M, K AND G OUTPUT FROM NASTRAN
C  *                    CMS ANALYSIS.
C  *
C  * EXTERNAL CALLS   : NONE
C  *
C  * AUTHOR           : TERJE ROLVAG
C  *
C  * DATE / VERSION   : 2001-06-28 / 1.0
C  *                    2005-06-13 / 2.0  K. M. Okstad (using F90 utils)
C  *
C***********************************************************************
C
      use kindModule       , only : dp
      use manipMatrixModule, only : writeObject
      use binaryDBInterface, only : openBinaryDB, closeBinaryDB
      use binaryDBInterface, only : writeTagDB, writeDoubleDB, write_p
C
      implicit none
C
      INTEGER       NCOL, NTOT, ILOG, IERR
      REAL(kind=dp) MATNAS(NTOT,NCOL)
      CHARACTER*(*) OUTF, FILETAG
      INTEGER       ifile, chksum
C
      WRITE(ILOG,*) 'Reduced matrix written to file: ',OUTF
      call writeObject (MATNAS,ILOG,nellIn=10)
C
      call openBinaryDB (OUTF,write_p,ifile,ierr)
      if (ierr < 0) return
C
      chksum = 0
      call writeTagDB (ifile,FILETAG,chksum,ierr)
      if (ierr < 0) return
C
      call writeDoubleDB (ifile,MATNAS(1,1),NTOT*NCOL,ierr)
      if (ierr < 0) return
C
      call closeBinaryDB (ifile,ierr)
C
      RETURN
      END
