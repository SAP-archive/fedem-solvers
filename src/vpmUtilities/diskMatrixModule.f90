!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module DiskMatrixModule

  !!============================================================================
  !! This module contains a data type and associated utility routines for
  !! representing a general rectangular matrix stored on a disk file.
  !! The matrix is stored column-wise on disk and it is therefore adviced to
  !! access the matrix in a column-wise fashion to obtain optimal speed.
  !! Note that all data members of the DiskMatrixType are private and can thus
  !! be accessed only through the subroutines/functions defined in this module.
  !!============================================================================

  use KindModule, only : sp, dp, i8

  implicit none

  integer, parameter :: dmOld_p = 0 !! Read-only matrix
  integer, parameter :: dmNsp_p = 1 !! Write matrix in single precision
  integer, parameter :: dmNdp_p = 2 !! Write matrix in double precision

  character(len=30), private, parameter :: dmTag_p = '#FEDEM disk matrix'

  type DiskMatrixType
     private
     logical             :: isReadOnly !! Write protected when .true.
     integer             :: nRows      !! Number of rows in the entire matrix
     integer             :: nCols      !! Number of columns in the entire matrix
     integer             :: nColSwap   !! Number of columns in one swap section
     integer             :: nElSwap    !! Number of elements in one swap section
     integer             :: nElastSwap !! Number of elements in the last section
     integer             :: nSwap      !! Number of swap sections
     integer             :: curSwap    !! Current swap section in core
     integer             :: swapFile   !! Swap file number
     integer             :: swapCnt(2) !! Number of swap operations performed
     integer(i8),pointer :: fileKey(:) !! File address for all swap sections
     real(dp)   ,pointer :: vala(:,:)  !! Matrix elements for one swap section
     real(sp)   ,pointer :: vals(:)    !! Single precision buffer
  end type DiskMatrixType

  private :: dmSwap, dmGetAddress, dmSingleDoubleCast


contains

  function dmGetSwapSize (n1,n2,autoRamSizeRatio) result(nw)

    !!==========================================================================
    !! Returns the requested size of one swap section from command-line option.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 31 Jan 2002/1.0
    !!==========================================================================

    use reportErrorModule     , only : reportError, note_p
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getint
    use FFaProfilerInterface  , only : ffa_getPhysMem

    integer , optional, intent(in) :: n1, n2
    real(sp), optional, intent(in) :: autoRamSizeRatio

    !! Local variables
    logical            :: printNote
    integer            :: nw, swapMB, swapW
    integer, parameter :: MByte = 131072 ! Number of dp words in one megabyte
    character(len=128) :: msg

    !! --- Logic section ---

    call ffa_cmdlinearg_getint('Bramsize',swapMB) ! Swap size in MegaBytes
    call ffa_cmdlinearg_getint('dmramsize',swapW) ! Swap size in dp Words
    if (int(swapMB,i8)*int(MByte,i8) > int(huge(1),i8)) then
       nw = huge(1)        ! = 2 GigaWords
    else if (swapMB >= 0) then
       nw = swapMB * MByte ! = swapMB * 1024^2 / 8
    else
       nw = swapW
    end if

    if (nw >= 0 .or. .not.present(autoRamSizeRatio)) return

    !! Compute max swap size based on the current available memory
    swapMB = ffa_getPhysMem(.false.)
    if (swapMB <= 0) return

    nw = int(real(swapMB,sp)*autoRamSizeRatio)
    if (present(n1) .and. present(n2)) then
       !! Print note only if we actually need to go out-of-core
       printNote = int(nw,i8) <= 1_i8 + int(n1,i8)*int(n2,i8)/int(MByte,i8)
    else
       printNote = .true.
    end if

    if (printNote) then
       write(msg,650) nw, 100*nw/swapMB
650    format('Setting recovery matrix RAM size to max',I5,'MB (',I3, &
            & '% of available memory).')
       call reportError (note_p,msg,'Use option -Bramsize to override.')
    end if

    nw = nw * MByte

  end function dmGetSwapSize


  subroutine dmNullify (this)

    !!==========================================================================
    !! Initializes the DiskMatrixType object.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 14 Feb 2002/1.0
    !!==========================================================================

    type(DiskMatrixType), intent(out) :: this

    !! --- Logic section ---

    this%nRows      =  0
    this%nCols      =  0
    this%nColSwap   =  0
    this%nElSwap    =  0
    this%nElastSwap =  0
    this%nSwap      =  0
    this%curSwap    =  0
    this%swapFile   = -1
    this%swapCnt    =  0

    nullify(this%fileKey)
    nullify(this%vala)
    nullify(this%vals)

  end subroutine dmNullify


  subroutine dmWrite (this,lpu,name)

    !!==========================================================================
    !! Writes out the DiskMatrixType object to unit lpu.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 18 Mar 2003/1.0
    !!==========================================================================

    use manipMatrixModule, only : writeObject

    type(DiskMatrixType), intent(inout) :: this
    integer             , intent(in)    :: lpu
    character(len=*)    , intent(in)    :: name

    !! Local variables
    integer , parameter :: nelLin = 10
    integer             :: i, ierr
    real(dp), pointer   :: column(:)
    character(len=16)   :: chCol

    !! --- Logic section ---

    if (.not.associated(this%vala) .or. this%nSwap < 1) return

#if defined(win32) || defined(win64)
    if (this%nSwap == 1 .and. this%nelSwap < 64000) then
#else
    if (this%nSwap == 1) then
#endif

       !! Write the whole matrix in one go
       call writeObject (this%vala,lpu,name,nelLin)

    else

       !! Write out the matrix column by column
       do i = 1, this%nCols
          column => dmGetColPtr(this,i,ierr)
          if (ierr < 0) exit
          write(chCol,"('  column',I8)") i
          call writeObject (column,lpu,name//chCol,nelLin)
       end do

    end if

  end subroutine dmWrite


  subroutine dmOpen (this,fileName,fileChkSum,nRows,nCols,status,err, &
       &             nColSwap,swapSize,wantTag)

    !!==========================================================================
    !! Opens a disk matrix file for Read/Write or ReadOnly operations.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 1 Jun 2001/1.0
    !!==========================================================================

    use kindModule       , only : nbd_p, nbs_p, nbi_p, maxInt_p
    use allocationModule , only : reAllocate, StrBytes
    use reportErrorModule, only : getErrorFile, reportError, error_p
    use binaryDBInterface, only : openBinaryDB, closeBinaryDB, read_p, update_p
    use binaryDBInterface, only : readTagDB, writeTagDB, getPositionDB
    use binaryDBInterface, only : readFloatDB, writeFloatDB
    use binaryDBInterface, only : readDoubleDB, writeDoubleDB

    type(DiskMatrixType), intent(out)   :: this
    character(len=*)    , intent(in)    :: fileName
    integer             , intent(inout) :: fileChkSum
    integer             , intent(in)    :: nRows, nCols, status
    integer             , intent(out)   :: err
    integer, optional   , intent(in)    :: nColSwap ! Columns per swap section
    integer, optional   , intent(in)    :: swapSize ! Swap section size (words)
    character(len=*), optional, intent(in) :: wantTag

    !! Local variables
    logical           :: storeAsFloat
    integer,parameter :: maxFloatBuf_p = 4096
    integer           :: i, j, lpu, nElSwap, nFloat, nRead
    integer(i8)       :: nBytes
    character(len=32) :: fileTag
    character(len=64) :: errMsg

    !! --- Logic section ---

    this%isReadOnly = status == dmOld_p
    this%nRows      = nRows
    this%nCols      = nCols
    this%nColSwap   = nCols ! Default: store full matrix in memory
    this%curSwap    = 0
    this%swapCnt    = 0
    storeAsFloat    = .false.

    !! Set number of columns in swap section
    if (present(nColSwap)) then
       if (nColSwap > 0) then
          this%nColSwap = adjustSwapSize(nCols,nColSwap)
       end if
    else if (present(swapSize) .and. nRows > 0 .and. nCols > 0) then
       if (swapSize > 0) then
          this%nColSwap = adjustSwapSize(nCols,swapSize/nRows)
       end if
    end if

    !! Check that the swap section size is within 2 Giga words
    nBytes = int(nRows,i8)*int(this%nColSwap,i8)
    if (nBytes > maxInt_p) then
       this%nColSwap = adjustSwapSize(nCols,maxInt_p/nRows)
    end if

    !! Set number of swap sections and number of elements in swap section
    this%nSwap      = 1 + (nCols-1)/max(1,this%nColSwap)
    this%nElSwap    = nRows * this%nColSwap
    this%nElastSwap = int(int(nRows,i8)*int(nCols,i8) - &
         &                int(this%nElSwap,i8)*int(this%nSwap-1,i8))

    !! Report matrix size
    lpu = getErrorFile()
    write(lpu,600) trim(fileName), nRows, nCols, this%nSwap, this%nElSwap
    if (this%nElastSwap < this%nElSwap) write(lpu,601) this%nElastSwap

    if (status == dmOld_p) then
       !! Open swap file and check the file tag
       call openBinaryDB (fileName,read_p,this%swapFile,err)
       if (err == 0) call readTagDB (this%swapFile,fileTag,fileChkSum,err)
    else
       !! Open swap file and write file tag and checksum record
       call openBinaryDB (fileName,update_p,this%swapFile,err)
       if (err == 0) then
          if (status == dmNsp_p) then
             storeAsFloat = .true.
             call writeTagDB (this%swapFile,trim(dmTag_p)//' SP',fileChkSum,err)
          else
             call writeTagDB (this%swapFile,dmTag_p,fileChkSum,err)
          end if
       end if
    end if

    if (err < 0) then
       call reportError (error_p,'Cannot open disk matrix file '//fileName, &
            &            addString='dmOpen')
       goto 110
    else if (status == dmOld_p) then
       !! Check whether the matrix was saved in single precision
       i = len_trim(fileTag)
       if (fileTag(i-2:i) == ' SP') then
          storeAsFloat = .true.
          i = i - 3
       end if
       if (present(wantTag)) then
          if (fileTag(1:i) /= wantTag) err = -1
       else
          if (fileTag(1:i) /= dmTag_p) err = -1
       end if
       if (err < 0) then
          call reportError (error_p,'Invalid disk matrix file '//fileName, &
               &            'File tag: '//fileTag,addString='dmOpen')
          goto 100
       end if
    end if

    !! Report allocated storage in core
#ifdef FT_HAS_INT4ONLY
    nBytes = 4_i8*int(this%nSwap,i8)
#else
    nBytes = 8_i8*int(this%nSwap,i8)
#endif
    nBytes = nBytes + int(11*nbi_p,i8) + int(nbd_p,i8)*int(this%nElSwap,i8)
    if (storeAsFloat) then
       nBytes = nBytes + int(nbs_p,i8)*int(min(maxFloatBuf_p,nRows),i8)
    end if
    errMsg = StrBytes(nBytes)
    write(lpu,602) trim(errMsg)

    !! Report allocated storage on file
    if (storeAsFloat) then
       nBytes = 46_i8 + int(nbs_p,i8)*dmSize(this)
    else
       nBytes = 46_i8 + int(nbd_p,i8)*dmSize(this)
    end if
    errMsg = StrBytes(nBytes)
    write(lpu,603) trim(errMsg)
    if (nRows < 1) return

    !! Allocate the in-core swap array
    call reAllocate ('dmOpen',this%fileKey,this%nSwap,err)
    call reAllocate ('dmOpen',this%vala,this%nRows,this%nColSwap,err)
    if (err < 0) goto 110
    if (status /= dmOld_p) this%vala = 0.0_dp

    if (storeAsFloat) then
       !! Allocate single precision buffer (max size: one matrix column)
       call reAllocate ('dmOpen',this%vals,min(maxFloatBuf_p,this%nRows),err)
       if (err < 0) goto 100
       if (status == dmNsp_p) this%vals = 0.0_sp
    end if

    !! Read or write through the file to initialize the file keys
    this%curSwap = 1
    nElSwap = this%nElSwap
    do i = 1, this%nSwap
       if (i == this%nSwap) nElSwap = this%nElastSwap

       call getPositionDB (this%swapFile,this%fileKey(i))
       if (this%fileKey(i) == 0_i8) then
          err = -i
          exit
       end if

       if (status == dmOld_p) then
          if (storeAsFloat) then
             nFloat = min(size(this%vals),nElSwap)
             do j = 1, nElSwap, max(1,nFloat)
                nRead = min(nFloat,nElSwap+1-j)
                call readFloatDB (this%swapFile,this%vals(1),nRead,err)
                if (err < 0) exit
                if (i == this%nSwap) then ! Last swap section, cast to double
                   call dmSingleDoubleCast (this,j,nRead,.false.)
                end if
             end do
          else
             call readDoubleDB (this%swapFile,this%vala(1,1),nElSwap,err)
          end if
          this%curSwap = i
       else if (storeAsFloat) then
          nFloat = min(size(this%vals),nElSwap)
          do j = 1, nElSwap, max(1,nFloat)
             nRead = min(nFloat,nElSwap+1-j)
             call writeFloatDB (this%swapFile,this%vals(1),nRead,err)
             if (err < 0) exit
          end do
       else
          call writeDoubleDB (this%swapFile,this%vala(1,1),nElSwap,err)
       end if
       if (err < 0) exit

    end do
    if (err >= 0) return

    !! Failure to read/write through the disk matrix file
    if (status == dmOld_p) then
       if (storeAsFloat) then
          nBytes = 46_i8 + dmSize(this)*int(nbs_p,i8)
       else
          nBytes = 46_i8 + dmSize(this)*int(nbd_p,i8)
       end if
       errMsg = 'The file is smaller than the anticipated '//StrBytes(nBytes)
    else
       errMsg = 'Check that you have sufficient disk space available.'
    end if
    call reportError (error_p,'Failed to initialize disk matrix file.', &
         &            errMsg,'File: '//fileName,addString='dmOpen')

100 call closeBinaryDB (this%swapFile,lpu)
110 call reAllocate ('dmOpen',this%fileKey)
    call reAllocate ('dmOpen',this%vala)
    call reAllocate ('dmOpen',this%vals)
    call dmNullify (this)

600 format(/4X,'Opening disk matrix file: ',A, &
         & /4X,'   Number of matrix rows     =',I10, &
         & /4X,'   Number of matrix columns  =',I10, &
         & /4X,'   Number of swap sections   =',I10, &
         & /4X,'   Size of swap sections     =',I10,' words')
601 format( 4X,'   Size of last swap section =',I10,' words')
602 format( 4X,'   Allocated storage         =      ',A)
603 format( 4X,'   Total file size           =      ',A)

  contains

    function adjustSwapSize (nTot,nSwap) result(n)
      integer, intent(in) :: nTot,nSwap
      integer :: n
      !!========================================================================
      !! Adjust the number of columns such that the total number of columns is
      !! "nearly" divisable by the number of swap columns, such that the last
      !! swap section is as close in size to the other ones as possible.
      !!========================================================================
      if (nSwap >= nTot) then
         n = nTot
      else if (nSwap <= 1) then
         n = 1
      else
         n = nSwap
         do while (n > 1 .and. mod(nTot,n) > 0)
            if (nTot/(n-1) > nTot/n .and. mod(nTot,n-1) > 0) exit
            n = n - 1
         end do
      end if
    end function adjustSwapSize

  end subroutine dmOpen


  subroutine dmSetReadOnly (this,err,isReadOnly)

    !!==========================================================================
    !! Toggles between Read/Write and ReadOnly. If switching to ReadOnly,
    !! the present in-core section is first written to file.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 1 Dec 2000/2.0
    !!==========================================================================

    use reportErrorModule, only : reportError, debugFileOnly_p

    type(DiskMatrixType), intent(inout) :: this
    integer             , intent(out)   :: err
    logical, optional   , intent(in)    :: isReadOnly

    !! --- Logic section ---

    if (present(isReadOnly)) then
       if (isReadOnly) then
          call dmSwap (this,0,err)
       else
          err = 0
       end if
       this%isReadOnly = isReadOnly
    else
       call dmSwap (this,0,err)
       this%isReadOnly = .true.
    end if

    if (err < 0) call reportError (debugFileOnly_p,'dmSetReadOnly')

  end subroutine dmSetReadOnly


  function dmSize (this,iDim)

    !!==========================================================================
    !! Returns the size of the given disk matrix, either the total size, in the
    !! the given matrix dimension, or the number of columns in the swap array.
    !!
    !! Programmer : Bjorn Haugen
    !! date/rev   : 29 Jan 1999/1.0
    !!==========================================================================

    type(DiskMatrixType), intent(in) :: this
    integer, optional   , intent(in) :: iDim
    integer(i8)                      :: dmSize

    !! --- Logic section ---

    if (.not. present(idim)) then
       dmSize = int(this%nRows,i8)*int(this%nCols,i8)
    else if (iDim == 1) then
       dmSize = int(this%nRows,i8)
    else if (iDim == 2) then
       dmSize = int(this%nCols,i8)
    else
       dmSize = int(this%nColSwap,i8) ! Number of columns in swap array
    end if

  end function dmSize


  subroutine dmClose (this,err)

    !!==========================================================================
    !! Closes the disk matrix file and releases the associated swap array.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 1 Dec 2000/2.0
    !!==========================================================================

    use allocationModule , only : reAllocate
    use reportErrorModule, only : getErrorFile, reportError, debugFileOnly_p
    use binaryDBInterface, only : closeBinaryDB

    type(DiskMatrixType), intent(inout) :: this
    integer             , intent(out)   :: err

    !! Local variables
    integer :: lpu

    !! --- Logic section ---

    if (this%nElSwap > 0) then
       !! Write current swap section to file
       call dmSwap (this,0,err)
       if (err < 0) goto 110
    else
       err = 0
    end if

    if (this%swapCnt(1) > 0 .or. this%swapCnt(2) > 0) then
       lpu = getErrorFile()
       write(lpu,600) this%swapCnt
    end if

    if (this%nElSwap > 0) then
       !! Close the swap file
       call closeBinaryDB (this%swapFile,err)
       if (err < 0) goto 110
    end if

100 call reAllocate ('dmClose',this%fileKey)
    call reAllocate ('dmClose',this%vala)
    call reAllocate ('dmClose',this%vals)
    call dmNullify (this)
    return

110 call reportError (debugFileOnly_p,'dmClose')
    goto 100

600 format(/4X,'Closing disk matrix file:', &
         & /4X,'   Number write operations performed =',I6, &
         & /4X,'   Number read operations performed  =',I6)

  end subroutine dmClose


  subroutine dmGetAddress (this,rInd,cInd,section,index,err)

    !!==========================================================================
    !! Gets the swap section number and the array index within that swap section
    !! corresponding to the given matrix indices (rInd,cInd).
    !!
    !! Programmer : Bjorn Haugen
    !! date/rev   : 10 Feb 1999/1.0
    !!==========================================================================

    use reportErrorModule, only : internalError

    type(DiskMatrixType), intent(inout) :: this
    integer             , intent(in)    :: rInd, cInd
    integer             , intent(out)   :: section, index, err

    !! --- Logic section ---

    if ( rInd <= 0 .or. rInd > this%nRows .or. &
         cInd <= 0 .or. cInd > this%nCols ) then
       section = 0
       index = 0
       err = internalError('dmGetAddress: Invalid matrix indices')
    else
       section = (cInd-1)/max(1,this%nColSwap) + 1
       index = cInd - this%nColSwap*(section-1)
       err = 0
    end if

  end subroutine dmGetAddress


  subroutine dmSwap (this,newSection,err)

    !!==========================================================================
    !! Writes current section to disk and reads section newSection into core.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 1 Dec 2000/2.0
    !!==========================================================================

    use reportErrorModule   , only : internalError, reportError, error_p
    use binaryDBInterface   , only : setPositionDB
    use binaryDBInterface   , only : readFloatDB, writeFloatDB
    use binaryDBInterface   , only : readDoubleDB, writeDoubleDB
    use FFaProfilerInterface, only : ffa_starttimer, ffa_stoptimer

    type(DiskMatrixType), intent(inout) :: this
    integer             , intent(in)    :: newSection
    integer             , intent(out)   :: err

    !! Local variables
    integer :: j, nElSwap, nFloat, nRead

    !! --- Logic section ---

    if (newSection == this%curSwap) then
       err = 0
       return
    else if (newSection > this%nSwap) then
       err = internalError('dmSwap: Invalid swap section number')
       return
    end if

    if (.not. this%isReadOnly .and. this%curSwap > 0) then
       if (this%curSwap == this%nSwap) then
          nElSwap = this%nElastSwap
       else
          nElSwap = this%nElSwap
       end if

       !! Write current swap section to file
       call ffa_starttimer('dmSwap:write')
       call setPositionDB (this%swapFile,this%fileKey(this%curSwap),err)
       if (err < 0) goto 100
       if (associated(this%vals)) then
          nFloat = min(size(this%vals),nElSwap)
          do j = 1, nElSwap, max(1,nFloat)
             nRead = min(nFloat,nElSwap+1-j)
             call dmSingleDoubleCast (this,j,nRead,.true.) ! Cast to single
             call writeFloatDB (this%swapFile,this%vals(1),nRead,err)
             if (err < 0) goto 110
          end do
       else
          call writeDoubleDB (this%swapFile,this%vala(1,1),nElSwap,err)
          if (err < 0) goto 110
       end if
       call ffa_stoptimer('dmSwap:write')
       this%swapCnt(1) = this%swapCnt(1) + 1
    end if

    if (newSection > 0) then
       if (newSection == this%nSwap) then
          nElSwap = this%nElastSwap
       else
          nElSwap = this%nElSwap
       end if

       !! Read specified section into memory
       call ffa_starttimer('dmSwap:read')
       call setPositionDB (this%swapFile,this%fileKey(newSection),err)
       if (err < 0) goto 100
       if (associated(this%vals)) then
          nFloat = min(size(this%vals),nElSwap)
          do j = 1, nElSwap, max(1,nFloat)
             nRead = min(nFloat,nElSwap+1-j)
             call readFloatDB (this%swapFile,this%vals(1),nRead,err)
             if (err < 0) goto 120
             call dmSingleDoubleCast (this,j,nRead,.false.) ! Cast to double
          end do
       else
          call readDoubleDB (this%swapFile,this%vala(1,1),nElSwap,err)
          if (err < 0) goto 120
       end if
       call ffa_stoptimer('dmSwap:read')
       this%curSwap = newSection
       this%swapCnt(2) = this%swapCnt(2) + 1
    end if

    err = 0
    return

100 call ffa_stoptimer('dmSwap:read')
    call ffa_stoptimer('dmSwap:write')
    call reportError (error_p,'Could not set swap file pointer', &
         &            addString='dmSwap')
    return

110 call ffa_stoptimer('dmSwap:write')
    call reportError (error_p,'Could not write to swap file', &
         &            addString='dmSwap')
    return

120 call ffa_stoptimer('dmSwap:read')
    call reportError (error_p,'Could not read from swap file', &
         &            addString='dmSwap')
    return

  end subroutine dmSwap


  subroutine dmSingleDoubleCast (this,ielSwap,nelSwap,toSingle)

    !!==========================================================================
    !! Casts the in-core matrix elements to/from single precision.
    !! It is assumed that the sp buffer array is not larger than one column.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 19 May 2005/1.0
    !!==========================================================================

    type(DiskMatrixType), intent(inout) :: this
    integer             , intent(in)    :: ielSwap, nelSwap
    logical             , intent(in)    :: toSingle

    !! Local variables
    integer :: i, nr, ic = 1, lastRow = 0

    !! --- Logic section ---

    if (ielSwap == 1) then
       lastRow = 0
       ic = 1
    end if

    nr = min(nelSwap,this%nRows-lastRow)
    if (toSingle) then
       do i = 1, nr
          this%vals(i) = real(this%vala(lastRow+i,ic),sp)
       end do
    else ! toDouble
       do i = 1, nr
          this%vala(lastRow+i,ic) = real(this%vals(i),dp)
       end do
    end if

    lastRow = mod(lastRow+nr,this%nRows)
    if (lastRow == 0) ic = ic + 1
    if (nr == nelSwap) return

    nr = nelSwap - nr
    if (toSingle) then
       do i = 1, nr
          this%vals(nelSwap-nr+i) = real(this%vala(i,ic),sp)
       end do
    else ! to Double
       do i = 1, nr
          this%vala(i,ic) = real(this%vals(nelSwap-nr+i),dp)
       end do
    end if

    lastRow = nr

  end subroutine dmSingleDoubleCast


  subroutine dmSetValue (this,value,rInd,cInd,err)

    !!==========================================================================
    !! Assigns a value to matrix element (rInd,cInd).
    !!
    !! Programmer : Bjorn Haugen
    !! date/rev   : 10 Feb 1999/1.0
    !!==========================================================================

    use reportErrorModule, only : reportError, debugFileOnly_p

    type(DiskMatrixType), intent(inout) :: this
    real(dp)            , intent(in)    :: value
    integer             , intent(in)    :: rInd, cInd
    integer             , intent(out)   :: err

    !! Local variables
    integer :: section, index

    !! --- Logic section ---

    err = -1
    if (this%isReadOnly) goto 100

    call dmGetAddress (this,rInd,cInd,section,index,err)
    if (err < 0) goto 100

    call dmSwap (this,section,err)
    if (err < 0) goto 100

    this%vala(rInd,index) = value
    return

100 call reportError (debugFileOnly_p,'dmSetValue')

  end subroutine dmSetValue


  function dmGetValue (this,rInd,cInd,err)

    !!==========================================================================
    !! Returns the value of matrix element (rInd,cInd).
    !!
    !! Programmer : Bjorn Haugen
    !! date/rev   : 10 Feb 1999/1.0
    !!==========================================================================

    use reportErrorModule, only : reportError, debugFileOnly_p

    type(DiskMatrixType), intent(inout) :: this
    integer             , intent(in)    :: rInd, cInd
    integer             , intent(out)   :: err
    real(dp)                            :: dmGetValue

    !! Local variables
    integer :: index, section

    !! --- Logic section ---

    call dmGetAddress (this,rInd,cInd,section,index,err)
    if (err < 0) goto 100

    call dmSwap (this,section,err)
    if (err < 0) goto 100

    dmGetValue = this%vala(rInd,index)
    return

100 call reportError (debugFileOnly_p,'dmGetValue')
    dmGetValue = 0.0_dp

  end function dmGetValue


  function dmGetSwapSecPtr (this,cInd,err)

    !!==========================================================================
    !! Returns a pointer to a matrix containing the entire swap section
    !! starting at column cInd of the given disk matrix.
    !! The function can be used for speedy random access of the matrix elements.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 13 Sep 2004/1.0
    !!==========================================================================

    use reportErrorModule, only : reportError, debugFileOnly_p

    type(DiskMatrixType), intent(inout) :: this
    integer             , intent(in)    :: cInd
    integer             , intent(out)   :: err
    real(dp)            , pointer       :: dmGetSwapSecPtr(:,:)

    !! Local variables
    integer :: index, section

    !! --- Logic section ---

    call dmGetAddress (this,1,cInd,section,index,err)
    if (err < 0) goto 100

    call dmSwap (this,section,err)
    if (err < 0) goto 100

    dmGetSwapSecPtr => this%vala
    return

100 call reportError (debugFileOnly_p,'dmGetSwapSecPtr')
    nullify(dmGetSwapSecPtr)

  end function dmGetSwapSecPtr


  function dmGetColPtr (this,cInd,err)

    !!==========================================================================
    !! Returns a pointer to column cInd of the given disk matrix.
    !!
    !! Caution: This function circumvents all protection and should be used with
    !! extreme caution and only for read-only purposes when speed is essential.
    !! Other calls to dm... routines can also change the contents of the memory
    !! being pointed to.
    !!
    !! Programmer : Bjorn Haugen and Karl Erik Thoresen
    !! date/rev   : 5 May 1999/1.0
    !!==========================================================================

    use reportErrorModule, only : reportError, debugFileOnly_p

    type(DiskMatrixType), intent(inout) :: this
    integer             , intent(in)    :: cInd
    integer             , intent(out)   :: err
    real(dp)            , pointer       :: dmGetColPtr(:)

    !! Local variables
    integer :: index, section

    !! --- Logic section ---

    call dmGetAddress (this,1,cInd,section,index,err)
    if (err < 0) goto 100

    call dmSwap (this,section,err)
    if (err < 0) goto 100

    dmGetColPtr => this%vala(:,index)
    return

100 call reportError (debugFileOnly_p,'dmGetColPtr')
    nullify(dmGetColPtr)

  end function dmGetColPtr


  subroutine dmGetCol (this,cInd,col,err)

    !!==========================================================================
    !! Gets column cInd of the given disk matrix.
    !!
    !! Programmer : Bjorn Haugen
    !! date/rev   : 10 Feb 1999/1.0
    !!==========================================================================

    use reportErrorModule, only : internalError, reportError, debugFileOnly_p

    type(DiskMatrixType), intent(inout) :: this
    integer             , intent(in)    :: cInd
    real(dp)            , intent(out)   :: col(:)
    integer             , intent(out)   :: err

    !! Local variables
    integer :: index, section

    !! --- Logic section ---

    if (size(col) /= this%nRows) then
       err = internalError('dmGetCol: Matrix dimension mismatch')
       return
    end if

    call dmGetAddress (this,1,cInd,section,index,err)
    if (err < 0) goto 100

    call dmSwap (this,section,err)
    if (err < 0) goto 100

    col = this%vala(:,index)
    return

100 call reportError (debugFileOnly_p,'dmGetCol')

  end subroutine dmGetCol


  subroutine dmGetRow (this,rInd,row,err)

    !!==========================================================================
    !! Gets row rInd of the given disk matrix.
    !!
    !! Programmer : Bjorn Haugen
    !! date/rev   : 10 Feb 1999/1.0
    !!==========================================================================

    use reportErrorModule, only : internalError, reportError, debugFileOnly_p

    type(DiskMatrixType), intent(inout) :: this
    integer             , intent(in)    :: rInd
    real(dp)            , intent(out)   :: row(:)
    integer             , intent(out)   :: err

    !! Local variables
    integer :: i, iSect, cInd, index, section

    !! --- Logic section ---

    if (size(row) < this%nCols) then
       err = internalError('dmGetRow: Matrix dimension mismatch')
       return
    end if

    cInd = 1
    do iSect = 1, this%nSwap
       call dmGetAddress (this,rInd,cInd,section,index,err)
       if (err < 0) exit

       call dmSwap (this,section,err)
       if (err < 0) exit

       do i = 1, this%nColSwap
          if (cInd > this%nCols) exit ! Last section may be less than nColSwap
          row(cInd) = this%vala(rInd,index)
          cInd = cInd + 1
          index = index + 1
       end do

    end do

    if (err < 0) call reportError (debugFileOnly_p,'dmGetRow')

  end subroutine dmGetRow


  subroutine dmMatTimesVec (this,xVec,yVec,err,doInitialize)

    !!==========================================================================
    !! Multiplies a disk matrix with a given vector, y = A*x  or  y = y + A*x.
    !!
    !! Programmer : Bjorn Haugen
    !! date/rev   : 10 Feb 1999/1.0
    !!==========================================================================

    use reportErrorModule, only : internalError, reportError, debugFileOnly_p

    type(DiskMatrixType), intent(inout) :: this
    real(dp)            , intent(in)    :: xVec(:)
    real(dp)            , intent(inout) :: yVec(:)
    integer             , intent(out)   :: err
    logical, optional   , intent(in)    :: doInitialize

    !! Local variables
    integer :: i, nRows, section, index

    !! --- Logic section ---

    err = 0
    if (this%nRows < 1) return
    if (this%nCols < 1) return

    nRows = this%nRows
    if (size(xVec) < this%nCols .or. size(yVec) /= nRows) then
       err = internalError('dmMatTimesVec: Matrix dimension mismatch')
       return
    end if

    if (.not. present(doInitialize)) then
       yVec = 0.0_dp
    else if (doInitialize) then
       yVec = 0.0_dp
    end if

    !! Perform the multiplication using SAXPY in order to finish the columns

    do i = 1, this%nCols
       call dmGetAddress (this,1,i,section,index,err)
       if (err < 0) exit

       call dmSwap (this,section,err)
       if (err < 0) exit

       yVec = yVec + this%vala(:,index)*xVec(i)
    end do

    if (err < 0) call reportError (debugFileOnly_p,'dmMatTimesVec')

  end subroutine dmMatTimesVec


  subroutine dmMatTransTimesVec (this,xVec,yVec,err)

    !!==========================================================================
    !! Multiplies the transpose of a disk matrix with a given vector, y = A^t*x.
    !!
    !! Programmer : Bjorn Haugen
    !! date/rev   : 10 Feb 1999/1.0
    !!==========================================================================

    use reportErrorModule, only : internalError, reportError, debugFileOnly_p

    type(DiskMatrixType), intent(inout) :: this
    real(dp)            , intent(in)    :: xVec(:)
    real(dp)            , intent(out)   :: yVec(:)
    integer             , intent(out)   :: err

    !! Local variables
    integer :: i, nRows, section, index

    !! --- Logic section ---

    err = 0
    if (this%nRows < 1) return
    if (this%nCols < 1) return

    nRows = this%nRows
    if (size(xVec) /= nRows .or. size(yVec) < this%nCols) then
       err = internalError('dmMatTransTimesVec: Matrix dimension mismatch')
       return
    end if

    !! Perform the multiplication using dot-products
    !! in order to work in a column-like fashion

    do i = 1, this%nCols
       call dmGetAddress (this,1,i,section,index,err)
       if (err < 0) exit

       call dmSwap (this,section,err)
       if (err < 0) exit

       yVec(i) = dot_product(this%vala(:,index),xVec)
    end do

    if (err < 0) call reportError (debugFileOnly_p,'dmMatTransTimesVec')

  end subroutine dmMatTransTimesVec


  subroutine dmFindAbsMaxRowValues (this,maxRowVal,err)

    !!==========================================================================
    !! Finds the absolute maximum value for each row in the given disk matrix.
    !!
    !! Programmer : Bjorn Haugen
    !! date/rev   : 25 Mar 2003/1.0
    !!==========================================================================

    use reportErrorModule, only : reportError, debugFileOnly_p

    type(DiskMatrixType), intent(inout) :: this
    real(dp)            , intent(out)   :: maxRowVal(:)
    integer             , intent(out)   :: err

    !! Local variables
    integer  :: i, j, nRows, section, index
    real(dp) :: colVal

    !! --- Logic section ---

    err = 0
    nRows = min(size(maxRowVal),this%nRows)
    if (nRows < 1) return

    !! Search for largest row value column by column

    do i = 1, this%nCols
       call dmGetAddress (this,1,i,section,index,err)
       if (err < 0) exit

       call dmSwap (this,section,err)
       if (err < 0) exit

       if (i == 1) then
          maxRowVal(1:nRows) = abs(this%vala(1:nRows,index))
       else
          do j = 1, nRows
             colVal = abs(this%vala(j,index))
             if (colVal > maxRowVal(j)) maxRowVal(j) = colVal
          end do
       end if
    end do

    if (err < 0) call reportError (debugFileOnly_p,'dmFindAbsMaxRowValues')

  end subroutine dmFindAbsMaxRowValues


  subroutine dmFindMinConnections (this,minIdx,minVal,err)

    !!==========================================================================
    !! Find the lowest connections represented by the given disk matrix.
    !!
    !! Programmer : Knut Morten Okstad                date/rev : 25 Feb 2019/1.0
    !!==========================================================================

    use kindModule        , only : dp, hugeVal_p
    use scratchArrayModule, only : getRealScratchArray
    use reportErrorModule , only : reportError, debugFileOnly_p

    type(DiskMatrixType), intent(inout) :: this
    integer             , intent(out)   :: minIdx(:)
    real(dp)            , intent(out)   :: minVal(:)
    integer             , intent(out)   :: err

    !! Local variables
    integer           :: i
    real(dp), pointer :: absMaxRowValues(:)

    !! --- Logic section ---

    absMaxRowValues => getRealScratchArray(this%nRows,err)
    if (err < 0) goto 100

    call dmFindAbsMaxRowValues (this,absMaxRowValues,err)
    if (err < 0) goto 100

    do i = 1, size(minIdx)
       minIdx(i) = minloc(absMaxRowValues,1)
       if (i <= size(minVal)) then
          minVal(i) = absMaxRowValues(minIdx(i))
       end if
       absMaxRowValues(minIdx(i)) = hugeVal_p
    end do

    return

100 call reportError (debugFileOnly_p,'dmFindMinConnections')

  end subroutine dmFindMinConnections

end module DiskMatrixModule
