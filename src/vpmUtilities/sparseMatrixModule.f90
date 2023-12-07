!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module SparseMatrixModule

  !!============================================================================
  !! This module contains a data type and associated utility routines for
  !! representing a general sparse matrix stored on classic coordinate form.
  !! Note that all data members of the SparseMatrixType are private and can
  !! thus be accessed only through the subroutines contained in this module.
  !!============================================================================

  use KindModule, only : dp

  implicit none

  integer , parameter, private :: initSize_p   = 10000
  real(dp), parameter, private :: growFactor_p = 1.5_dp
  real(dp), parameter, private :: epsValue_p   = 1.0e-30_dp

  type SparseMatrixType
     private
     integer           :: nNonZero     ! Number of non-zero elements
     integer           :: nRows, nCols ! Number of rows and columns
     real(dp), pointer :: vala(:)      ! Values of non-zero elements
     integer , pointer :: iRow(:)      ! Row indices of the non-zero entries
     integer , pointer :: jCol(:)      ! Column indices of the non-zero entries
  end type SparseMatrixType

  private :: smGrow


contains

  subroutine smNullify (this)

    !!==========================================================================
    !! Initialize the SparseMatrixType object.
    !!
    !! Programmer : Knut Morten Okstad                   date/rev : Mar 2003/1.0
    !!==========================================================================

    type(SparseMatrixType), intent(out) :: this

    !! --- Logic section ---

    this%nRows    = 0
    this%nCols    = 0
    this%nNonZero = 0

    nullify(this%vala)
    nullify(this%iRow)
    nullify(this%jCol)

  end subroutine smNullify


  subroutine smAllocate (this,nRows,nCols,useInitSize,err)

    !!==========================================================================
    !! Allocate a SparseMatrixType object.
    !!
    !! Programmer : Bjorn Haugen                         date/rev : Jan 1999/1.0
    !!              Knut Morten Okstad                              Jan 2003/2.0
    !!==========================================================================

    use allocationModule, only : reAllocate

    type(SparseMatrixType), intent(out) :: this
    integer               , intent(in)  :: nRows, nCols
    logical               , intent(in)  :: useInitSize
    integer               , intent(out) :: err

    !! --- Logic section ---

    err = 0
    call smNullify (this)

    this%nRows = nRows
    this%nCols = nCols

    if (useInitSize) then
       call reAllocate ('smAllocate',this%vala,initSize_p,err)
       call reAllocate ('smAllocate',this%iRow,initSize_p,err)
       call reAllocate ('smAllocate',this%jCol,initSize_p,err)
    end if

  end subroutine smAllocate


  subroutine smSetMatrix (this,mRow,mCol,values,err)

    !!==========================================================================
    !! Initiate the SparseMatrixType object with the given data.
    !! Note: The input arrays are referenced directly by the SparseMatrixType
    !! object (pointer-associated) and not copied. They must therefore not be
    !! accessed from outside the SparsMatrixType object after the object has
    !! been destroyed (by calling smDeallocate), as these arrays then will be
    !! deallocated as well.
    !!
    !! Programmer : Knut Morten Okstad                date/rev : 16 Sep 2004/1.0
    !!==========================================================================

    use KindModule      , only : nbi_p, nbd_p
    use allocationModule, only : doLogMem, logAllocMem, reAllocate

    type(SparseMatrixType), intent(inout) :: this
    integer               , pointer       :: mRow(:), mCol(:)
    real(dp)              , pointer       :: values(:)
    integer               , intent(out)   :: err

    !! --- Logic section ---

    err = 0
    if (.not.associated(values)) err = err - 1
    if (.not.associated(mRow))   err = err - 1
    if (.not.associated(mCol))   err = err - 1
    if (err < 0) return

    this%nNonZero = size(values)
    if (this%nNonZero /= size(mRow)) err = err - 1
    if (this%nNonZero /= size(mCol)) err = err - 1
    if (err < 0) return

    call reAllocate ('smSetMatrix',this%vala)
    call reAllocate ('smSetMatrix',this%iRow)
    call reAllocate ('smSetMatrix',this%jCol)
    this%vala => values
    this%iRow => mRow
    this%jCol => mCol

    if (doLogMem) then
       call logAllocMem ('smSetMatrix',0,size(this%vala),nbd_p)
       call logAllocMem ('smSetMatrix',0,size(this%iRow)+size(this%jCol),nbi_p)
    end if

  end subroutine smSetMatrix


  function smSize (this,iDim)

    !!==========================================================================
    !! Return the dimension of a SparseMatrixType object.
    !!
    !! Programmer : Bjorn Haugen                         date/rev : Jan 1999/1.0
    !!==========================================================================

    type(SparseMatrixType), intent(in) :: this
    integer, optional     , intent(in) :: iDim
    integer                            :: smSize

    !! --- Logic section ---

    if (.not. present(idim)) then
       smSize = this%nRows*this%nCols
    else if (iDim == 1) then
       smSize = this%nRows
    else if (iDim == 2) then
       smSize = this%nCols
    else
       smSize = this%nNonZero
    end if

  end function smSize


  subroutine smDeallocate (this)

    !!==========================================================================
    !! Deallocate a SparseMatrixType object.
    !!
    !! Programmer : Bjorn Haugen                         date/rev : Jan 1999/1.0
    !!              Knut Morten Okstad                              Jan 2003/2.0
    !!==========================================================================

    use allocationModule, only : reAllocate

    type(SparseMatrixType), intent(inout) :: this

    !! --- Logic section ---

    call reAllocate ('smDeallocate',this%vala)
    call reAllocate ('smDeallocate',this%iRow)
    call reAllocate ('smDeallocate',this%jCol)
    call smNullify (this)

  end subroutine smDeallocate


  subroutine smTranspose (this)

    !!==========================================================================
    !! Transpose a SparseMatrixType object.
    !!
    !! Programmer : Bjorn Haugen                         date/rev : Jan 1999/1.0
    !!==========================================================================

    type(SparseMatrixType), intent(inout) :: this

    !! Local variables
    integer :: tmp, i

    !! --- Logic section ---

    tmp        = this%nRows
    this%nRows = this%nCols
    this%nCols = tmp

    do i = 1, this%nNonZero
       tmp          = this%iRow(i)
       this%iRow(i) = this%jCol(i)
       this%jCol(i) = tmp
    end do

  end subroutine smTranspose


  subroutine smSetValue (this,value,rInd,cInd,err,lpu)

    !!==========================================================================
    !! Assign a value to the (i,j)-element of a SparseMatrixType object.
    !!
    !! Programmer : Bjorn Haugen                         date/rev : Jan 1999/1.0
    !!==========================================================================

    type(SparseMatrixType), intent(inout) :: this
    real(dp)              , intent(in)    :: value
    integer               , intent(in)    :: rInd, cInd
    integer               , intent(out)   :: err
    integer, optional     , intent(in)    :: lpu

    !! Local variables
    integer :: i

    !! --- Logic section ---

    err = 0
    if (abs(value) < epsValue_p) return

    !! Check if the indices already exist, and set the associated value if so
    do i = 1, this%nNonZero
       if (this%iRow(i) == rInd .and. this%jCol(i) == cInd) then
          this%vala(i) = value
          return
       end if
    end do

    !! New indices, increment number of nonZeroes and insert new value
    this%nNonZero = this%nNonZero + 1
    if (this%nNonZero > size(this%vala)) then
       call smGrow (this,err,lpu)
       if (err < 0) return
    end if

    i = this%nNonZero
    this%vala(i) = value
    this%iRow(i) = rInd
    this%jCol(i) = cInd

  end subroutine smSetValue


  function smGetValue (this,rInd,cInd)

    !!==========================================================================
    !! Return the (i,j)-element of a SparseMatrixType object.
    !!
    !! Programmer : Bjorn Haugen                         date/rev : Jan 1999/1.0
    !!==========================================================================

    type(SparseMatrixType), intent(in) :: this
    integer               , intent(in) :: rInd, cInd
    real(dp)                           :: smGetValue

    !! Local variables
    integer :: i

    !! --- Logic section ---

    smGetValue = 0.0_dp

    !! Check if indices represent a nonZero value
    do i = 1, this%nNonZero
       if (this%iRow(i) == rInd .and. this%jCol(i) == cInd) then
          smGetValue = this%vala(i)
          return
       end if
    end do

  end function smGetValue


  function smGetCol (this,cInd)

    !!==========================================================================
    !! Return the specified column of a SparseMatrixType object.
    !!
    !! Programmer : Bjorn Haugen                         date/rev : Jan 1999/1.0
    !!==========================================================================

    type(SparseMatrixType), intent(in) :: this
    integer               , intent(in) :: cInd
    real(dp)                           :: smGetCol(this%nRows)

    !! Local variables
    integer :: i

    !! --- Logic section ---

    smGetCol = 0.0_dp

    !! Insert the nonZero elements in the column vector
    do i = 1, this%nNonZero
       if (this%jCol(i) == cInd) smGetCol(this%iRow(i)) = this%vala(i)
    end do

  end function smGetCol


  function smGetRow (this,rInd)

    !!==========================================================================
    !! Return the specified row of a SparseMatrixType object.
    !!
    !! Programmer : Bjorn Haugen                         date/rev : Jan 1999/1.0
    !!==========================================================================

    type(SparseMatrixType), intent(in) :: this
    integer               , intent(in) :: rInd
    real(dp)                           :: smGetRow(this%nCols)

    !! Local variables
    integer :: i

    !! --- Logic section ---

    smGetRow = 0.0_dp

    !! Insert the nonZero elements in the row vector
    do i = 1, this%nNonZero
       if (this%iRow(i) == rInd) smGetRow(this%jCol(i)) = this%vala(i)
    end do

  end function smGetRow


  subroutine smMatTimesVec (this,xVec,yVec,err)

    !!==========================================================================
    !! Multiply a sparse matrix with a given vector, y = A*x.
    !!
    !! Programmer : Bjorn Haugen                         date/rev : Jan 1999/1.0
    !!==========================================================================

    use reportErrorModule, only : internalError

    type(SparseMatrixType), intent(in)  :: this
    real(dp)              , intent(in)  :: xVec(:)
    real(dp)              , intent(out) :: yVec(:)
    integer               , intent(out) :: err

    !! Local variables
    integer :: i, j, k

    !! --- Logic section ---

    if (size(xVec) == this%nCols .and. size(yVec) == this%nRows) then
       err = 0
    else
       err = internalError('smMatTimesVec: Dimension mis-match')
       return
    end if

    yVec = 0.0_dp

    !! Perform the multiplication
    do k = 1, this%nNonZero
       i = this%iRow(k)
       j = this%jCol(k)
       yVec(i) = yVec(i) + this%vala(k)*xVec(j)
    end do

  end subroutine smMatTimesVec


  subroutine smMatTransTimesVec (this,xVec,yVec,err)

    !!==========================================================================
    !! Multiply the transpose of a sparse matrix with a given vector, y = A^t*x.
    !!
    !! Programmer : Bjorn Haugen                         date/rev : Jan 1999/1.0
    !!==========================================================================

    use reportErrorModule, only : internalError

    type(SparseMatrixType), intent(in)  :: this
    real(dp)              , intent(in)  :: xVec(:)
    real(dp)              , intent(out) :: yVec(:)
    integer               , intent(out) :: err

    !! Local variables
    integer :: i, j, k

    !! --- Logic section ---

    if (size(xVec) == this%nRows .and. size(yVec) == this%nCols) then
       err = 0
    else
       err = internalError('smMatTransTimesVec: Dimension mis-match')
       return
    end if

    yVec = 0.0_dp

    !! Perform the multiplication
    do k = 1, this%nNonZero
       i = this%jCol(k)
       j = this%iRow(k)
       yVec(i) = yVec(i) + this%vala(k)*xVec(j)
    end do

  end subroutine smMatTransTimesVec


  subroutine smShrinkToFit (this,err)

    !!==========================================================================
    !! Reallocate a SparseMatrixType object to exactly fit its dimensions.
    !!
    !! Programmer : Bjorn Haugen                         date/rev : Jan 1999/1.0
    !!              Knut Morten Okstad                              Jan 2003/2.0
    !!==========================================================================

    use allocationModule, only : reAllocate

    type(SparseMatrixType), intent(inout) :: this
    integer               , intent(out)   :: err

    !! Local variables
    integer            :: newSize, oldSize
    logical, parameter :: preserveContent_p = .true.

    !! --- Logic section ---

    err     = 0
    oldSize = size(this%vala)
    newSize = this%nNonZero
    call reAllocate ('smShrinkToFit',this%vala,newSize,err,preserveContent_p)
    call reAllocate ('smShrinkToFit',this%iRow,newSize,err,preserveContent_p)
    call reAllocate ('smShrinkToFit',this%jCol,newSize,err,preserveContent_p)

  end subroutine smShrinkToFit


  subroutine smWrite (this,name,lpu,writeMatrixElms)

    !!==========================================================================
    !! Write out size parameters for the SparseMatrixType object.
    !! The matrix elements themselves may optionally be written out as well.
    !!
    !! Programmer : Knut Morten Okstad                   date/rev : Mar 2003/1.0
    !!==========================================================================

    use KindModule         , only : nbd_p, nbi_p, i8
    use AllocationModule   , only : StrBytes
    use ScratchArrayModule , only : getIntegerScratchArray
    use SearchAndSortModule, only : quickSortInt

    type(SparseMatrixType), intent(in) :: this
    character(len=*)      , intent(in) :: name
    integer               , intent(in) :: lpu
    logical, optional     , intent(in) :: writeMatrixElms

    !! Local variables
    integer, parameter :: nelLin = 10
    integer            :: i, j, lastRow, lenLine, nnz
    integer, pointer   :: indices(:)
    character(len=16)  :: cBytes

    !! --- Logic section ---

    write(lpu,600) name, this%nRows, this%nCols
    if (this%nNonZero > 0) then
       cBytes = StrBytes(int(this%nNonZero,i8)*int(nbd_p+2*nbi_p,i8))
       write(lpu,610) this%nNonZero,trim(adjustl(cBytes))
    else
       write(lpu,611)
       return
    end if
    if (this%nRows*this%nCols > 0) then
       write(lpu,620) real(this%nNonZero)/real(this%nRows*this%nCols)
    end if

    if (.not. present(writeMatrixElms)) return
    if (.not. writeMatrixElms) return

    !! Sort the non-zeroes in the printing order (row by row)
    nnz = this%nNonZero
    indices => getIntegerScratchArray(2*nnz,i)
    if (i < 0) return
    do i = 1, nnz
       indices(nnz+i) = this%iRow(i)*this%nCols + this%jCol(i)
    end do
    call quickSortInt (indices(nnz+1:2*nnz),indices(1:nnz))

    !! Write out the matrix elements, up to nelLin elements per row
    lastRow = 0
    lenLine = 0
    do i = 1, nnz
       if (this%iRow(indices(i)) /= lastRow) then
          lenLine = 0
          lastRow = this%iRow(indices(i))
          write(lpu,650) lastRow
       else if (lenLine == nelLin) then
          lenLine = 0
          write(lpu,*)
       end if
       if (lenLine == 0) then
          j = i ! Write the column numbers for current row
          do while (this%iRow(indices(j)) == lastRow .and. lenLine < nelLin)
             lenLine = lenLine + 1
             write(lpu,660,ADVANCE='no') this%jCol(indices(j))
             j = j + 1
             if (j > nnz) exit
          end do
          lenLine = 0
          write(lpu,*)
       end if
       lenLine = lenLine + 1
       write(lpu,670,ADVANCE='no') this%vala(indices(i))
    end do
    write(lpu,*)

600 format(/4X,'Storage requirement for sparse matrix ', A, &
         & /7X,'nRows    = ',I10,'  nCols   = ',I10 )
610 format( 7X,'nNonZero = ',I10,'  (',A,')')
611 format( 7X,'There are no non-zero elements in this matrix')
620 format( 7X,'Relative density = ',1PE12.5 )
650 format(/' +++++ Row',I8,' +++++')
660 format(I9,4X)
670 format(1PE13.5)

  end subroutine smWrite


  subroutine smGrow (this,err,lpu)

    !!==========================================================================
    !! Increase the allocated space for a SparseMatrixType object
    !! with a fixed amount.
    !!
    !! Comments by KET:
    !!    Number of times to call this routine to grow a matrix to size x;
    !!
    !!       n = log(x) / log(growfactor)
    !!
    !! Programmer : Bjorn Haugen                         date/rev : Jan 1999/1.0
    !!              Knut Morten Okstad                              Jan 2003/2.0
    !!==========================================================================

    use allocationModule, only : reAllocate

    type(SparseMatrixType), intent(inout) :: this
    integer               , intent(out)   :: err
    integer, optional     , intent(in)    :: lpu

    !! Local variables
    integer            :: newSize, oldSize
    logical, parameter :: preserveContent_p = .true.

    !! --- Logic section ---

    err     = 0
    oldSize = size(this%vala)
    newSize = max(int(oldSize*growFactor_p),initSize_p)
    if (present(lpu)) write(lpu,*) 'SparseMatrixType::smGrow new size :',newSize
    call reAllocate ('smGrow',this%vala,newSize,err,preserveContent_p)
    call reAllocate ('smGrow',this%iRow,newSize,err,preserveContent_p)
    call reAllocate ('smGrow',this%jCol,newSize,err,preserveContent_p)

  end subroutine smGrow

end module SparseMatrixModule
