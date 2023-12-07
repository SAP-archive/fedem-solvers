!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file allocationModule.f90
!> @brief Subroutines for dynamic allocation of integer and real arrays.

!!==============================================================================
!> @brief Module with subroutines for dynamic allocation of arrays.
!>
!> @details This module contains utility subroutines to facilitate allocation,
!> reallocation, and deallocation of integer, real and double precision arrays,
!> with optional log of accumulated memory usage. They work much in the similar
!> way as the C-function realloc(). The allocation subroutines are accessed via
!> the common interface allocationmodule::reallocate().
!>
!> @author Knut Morten Okstad
!> @date 25 Oct 2001

module allocationModule

  use KindModule, only : i8

  implicit none

  private

  !> @brief Allocates, reallocates or deallocates a generic array.
  interface reAllocate
     module procedure reAllocateInt
     module procedure reAllocateFloat1D
     module procedure reAllocateFloat2D
     module procedure reAllocateDouble1D
     module procedure reAllocateDouble2D
#if !defined(FT_HAS_INT4_ONLY) && !defined(FT_USE_INT8)
     module procedure reAllocateInt8
     module procedure reAllocateBigDouble1D
#endif
  end interface

  !> @brief Prints the accumulated memory usage to log file.
  interface logAllocMem
     module procedure logAllocMemI8
#if !defined(FT_HAS_INT4_ONLY) && !defined(FT_USE_INT8)
     module procedure logAllocMemI4
#endif
  end interface

  logical    , save :: doLogMem = .false. !< Toggle for print of memory usage
  integer    , save :: logFile = 0        !< File unit number for memory logging
  integer(i8), save :: totMem = 0_i8      !< Total number of bytes allocatedy
  integer(i8), save :: peakMem = 0_i8     !< Maximum number of bytes allocated

  public :: reAllocate, logAllocMem, closeMemLogFile, writePeakMem, writeStorage
  public :: StrBytes, doLogMem


contains

  !!============================================================================
  !> @brief Allocate, reallocate or deallocate a one-dimensional integer array.
  !>
  !> @param[in] label Name of calling subroutine (for memory logging)
  !> @param array The integer array to (re)allocate
  !> @param[in] nw Number of words in array after reallocation
  !> @param[out] ierr Error flag
  !> @param[in] preserveContent If @e .true., the existing contents is preserved
  !> when the array is reallocated

  subroutine reAllocateInt (label,array,nw,ierr,preserveContent)

    use KindModule       , only : nbi_p
    use reportErrorModule, only : allocationError

    character(len=*), intent(in)    :: label
    integer         , pointer       :: array(:)
    integer,optional, intent(in)    :: nw
    integer,optional, intent(inout) :: ierr
    logical,optional, intent(in)    :: preserveContent

    !! Local variables
    integer          :: newSize, oldSize, minSize, err
    integer, pointer :: tmp(:)

    !! --- Logic section ---

    if (present(nw)) then
       newSize = nw
    else
       newSize = 0
    end if
    oldSize = 0
    minSize = 0

    if (associated(array)) then
       oldSize = size(array)
       if (newSize == oldSize) return
       if (present(preserveContent)) then
          if (preserveContent) then
             minSize = min(oldSize,newSize)
             if (minSize > 0) then
                !! The array already contains data which we want to preserve.
                !! Copy into a temporary array and swap pointers after the
                !! old array has been deallocated.
                allocate(tmp(newSize),STAT=err)
                if (err /= 0) then
                   if (present(ierr)) ierr = ierr + allocationError(label)
                   newSize = 0
                   minSize = 0
                else
                   !! Use BLAS-routine to avoid stack problem for large arrays
                   call ICOPY (minSize,array(1),1,tmp(1),1)
                end if
             end if
          end if
       end if
       !! Release old array
       deallocate(array)
       nullify(array)
    end if

    if (minSize > 0) then ! Reallocation with preserving old content
       array => tmp
    else if (newSize > 0) then ! Allocation of new array
       allocate(array(newSize),STAT=err)
       if (err /= 0) then
          if (present(ierr)) ierr = ierr + allocationError(label)
          newSize = -newSize
          nullify(array)
       end if
    end if

    if (doLogMem) call logAllocMem (label,oldSize,newSize,nbi_p)

  end subroutine reAllocateInt


#if !defined(FT_HAS_INT4_ONLY) && !defined(FT_USE_INT8)
  !!============================================================================
  !> @brief Allocate, reallocate or deallocate a one-dimensional integer array.
  !>
  !> @param[in] label Name of calling subroutine (for memory logging)
  !> @param array The integer*8 array to (re)allocate
  !> @param[in] nw Number of words in array after reallocation
  !> @param[out] ierr Error flag
  !> @param[in] preserveContent If @e .true., the existing contents is preserved
  !> when the array is reallocated

  subroutine reAllocateInt8 (label,array,nw,ierr,preserveContent)

    use reportErrorModule, only : allocationError

    character(len=*), intent(in)    :: label
    integer(i8)     , pointer       :: array(:)
    integer,optional, intent(in)    :: nw
    integer,optional, intent(inout) :: ierr
    logical,optional, intent(in)    :: preserveContent

    !! Local variables
    integer              :: i, newSize, oldSize, minSize, err
    integer(i8), pointer :: tmp(:)

    !! --- Logic section ---

    if (present(nw)) then
       newSize = nw
    else
       newSize = 0
    end if
    oldSize = 0
    minSize = 0

    if (associated(array)) then
       oldSize = size(array)
       if (newSize == oldSize) return
       if (present(preserveContent)) then
          if (preserveContent) then
             minSize = min(oldSize,newSize)
             if (minSize > 0) then
                !! The array already contains data which we want to preserve.
                !! Copy into a temporary array and swap pointers after the
                !! old array has been deallocated.
                allocate(tmp(newSize),STAT=err)
                if (err /= 0) then
                   if (present(ierr)) ierr = ierr + allocationError(label)
                   newSize = 0
                   minSize = 0
                else
                   !! No BLAS for 64-bit integer
                   do i = 1, minSize
                      tmp(i) = array(i)
                   end do
                end if
             end if
          end if
       end if
       !! Release old array
       deallocate(array)
       nullify(array)
    end if

    if (minSize > 0) then ! Reallocation with preserving old content
       array => tmp
    else if (newSize > 0) then
       !! Allocation of new array
       allocate(array(newSize),STAT=err)
       if (err /= 0) then
          if (present(ierr)) ierr = ierr + allocationError(label)
          newSize = -newSize
          nullify(array)
       end if
    end if

    if (doLogMem) call logAllocMem (label,oldSize,newSize,8)

  end subroutine reAllocateInt8
#endif


  !!============================================================================
  !> @brief Allocate, reallocate or deallocate a one-dimensional real array.
  !>
  !> @param[in] label Name of calling subroutine (for memory logging)
  !> @param array The single precision real array to (re)allocate
  !> @param[in] nw Number of words in array after reallocation
  !> @param[out] ierr Error flag
  !>
  !> @author Knut Morten Okstad
  !> @date 30 Dec 2004

  subroutine reAllocateFloat1D (label,array,nw,ierr)

    use KindModule       , only : sp, nbs_p
    use reportErrorModule, only : allocationError

    character(len=*), intent(in)    :: label
    real(sp)        , pointer       :: array(:)
    integer,optional, intent(in)    :: nw
    integer,optional, intent(inout) :: ierr

    !! Local variables
    integer :: newSize, oldSize, err

    !! --- Logic section ---

    if (present(nw)) then
       newSize = nw
    else
       newSize = 0
    end if
    oldSize = 0

    if (associated(array)) then
       oldSize = size(array)
       if (newSize == oldSize) return
       !! Release old array
       deallocate(array)
       nullify(array)
    end if

    if (newSize > 0) then
       !! Allocation of new array
       allocate(array(newSize),STAT=err)
       if (err /= 0) then
          if (present(ierr)) ierr = ierr + allocationError(label)
          newSize = -newSize
          nullify(array)
       end if
    end if

    if (doLogMem) call logAllocMem (label,oldSize,newSize,nbs_p)

  end subroutine reAllocateFloat1D


  !!============================================================================
  !> @brief Allocate, reallocate or deallocate a two-dimensional real array.
  !>
  !> @param[in] label Name of calling subroutine (for memory logging)
  !> @param array The single precision real array to (re)allocate
  !> @param[in] n1 Number of rows in array after reallocation
  !> @param[in] n2 Number of columns in array after reallocation
  !> @param[out] ierr Error flag

  subroutine reAllocateFloat2D (label,array,n1,n2,ierr)

    use KindModule       , only : sp, nbs_p
    use reportErrorModule, only : allocationError

    character(len=*), intent(in)    :: label
    real(sp)        , pointer       :: array(:,:)
    integer,optional, intent(in)    :: n1, n2
    integer,optional, intent(inout) :: ierr

    !! Local variables
    integer :: newSize, oldSize, err

    !! --- Logic section ---

    if (present(n1) .and. present(n2)) then
       newSize = n1*n2
    else
       newSize = 0
    end if
    oldSize = 0

    if (associated(array)) then
       oldSize = size(array)
       if (newSize == oldSize) then
          if (.not.present(n1) .or. .not.present(n2)) return
          if (n1 == size(array,1) .and. n2 == size(array,2)) return
       end if
       !! Release old array
       deallocate(array)
       nullify(array)
    end if

    if (newSize > 0) then
       !! Allocation of new array
       allocate(array(n1,n2),STAT=err)
       if (err /= 0) then
          if (present(ierr)) ierr = ierr + allocationError(label)
          newSize = -newSize
          nullify(array)
       end if
    end if

    if (doLogMem) call logAllocMem (label,oldSize,newSize,nbs_p)

  end subroutine reAllocateFloat2D


  !!============================================================================
  !> @brief Allocate, reallocate or deallocate a one-dimensional real array.
  !>
  !> @param[in] label Name of calling subroutine (for memory logging)
  !> @param array The double precision real array to (re)allocate
  !> @param[in] nw Number of words in array after reallocation
  !> @param[out] ierr Error flag
  !> @param[in] preserveContent If @e .true., the existing contents is preserved
  !> when the array is reallocated

  subroutine reAllocateDouble1D (label,array,nw,ierr,preserveContent)

    use KindModule       , only : dp, nbd_p
    use reportErrorModule, only : allocationError

    character(len=*), intent(in)    :: label
    real(dp)        , pointer       :: array(:)
    integer,optional, intent(in)    :: nw
    integer,optional, intent(inout) :: ierr
    logical,optional, intent(in)    :: preserveContent

    !! Local variables
    integer           :: newSize, oldSize, minSize, err
    real(dp), pointer :: tmp(:)

    !! --- Logic section ---

    if (present(nw)) then
       newSize = nw
    else
       newSize = 0
    end if
    oldSize = 0
    minSize = 0

    if (associated(array)) then
       oldSize = size(array)
       if (newSize == oldSize) return
       if (present(preserveContent)) then
          if (preserveContent) then
             minSize = min(oldSize,newSize)
             if (minSize > 0) then
                !! The array already contains data which we want to preserve.
                !! Copy into a temporary array and swap pointers after the
                !! old array has been deallocated.
                allocate(tmp(newSize),STAT=err)
                if (err /= 0) then
                   if (present(ierr)) ierr = ierr + allocationError(label)
                   newSize = 0
                   minSize = 0
                else
                   !! Use BLAS-routine to avoid stack problem for large arrays
                   call DCOPY (minSize,array(1),1,tmp(1),1)
                end if
             end if
          end if
       end if
       !! Release old array
       deallocate(array)
       nullify(array)
    end if

    if (minSize > 0) then ! Reallocation with preserving old content
       array => tmp
    else if (newSize > 0) then ! Allocation of new array
       allocate(array(newSize),STAT=err)
       if (err /= 0) then
          if (present(ierr)) ierr = ierr + allocationError(label)
          newSize = -newSize
          nullify(array)
       end if
    end if

    if (doLogMem) call logAllocMem (label,oldSize,newSize,nbd_p)

  end subroutine reAllocateDouble1D


#if !defined(FT_HAS_INT4_ONLY) && !defined(FT_USE_INT8)
  !!============================================================================
  !> @brief Allocate, reallocate or deallocate a big one-dimensional real array.
  !>
  !> @param[in] label Name of calling subroutine (for memory logging)
  !> @param array The double precision real array to (re)allocate
  !> @param[in] nw Number of words in array after reallocation
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine is used for arrays that may contain more than
  !> 2^31 words, which is the largest number that can be represented by normal
  !> 4-byte integer variables.
  !>
  !> @author Knut Morten Okstad
  !> @date 12 May 2017

  subroutine reAllocateBigDouble1D (label,array,nw,ierr)

    use KindModule       , only : dp, nbd_p
    use reportErrorModule, only : allocationError

    character(len=*), intent(in)    :: label
    real(dp)        , pointer       :: array(:)
    integer(i8)     , intent(in)    :: nw
    integer         , intent(inout) :: ierr

    !! Local variables
    integer(i8) :: newSize, oldSize
    integer     :: err

    !! --- Logic section ---

    newSize = nw
    oldSize = 0_i8

    if (associated(array)) then
       oldSize = size(array,kind=i8)
       if (newSize == oldSize) return
       !! Release old array
       deallocate(array)
       nullify(array)
    end if

    if (newSize > 0_i8) then ! Allocation of new array
       allocate(array(newSize),STAT=err)
       if (err /= 0) then
          ierr = ierr + allocationError(label)
          newSize = -newSize
          nullify(array)
       end if
    end if

    if (doLogMem) call logAllocMem (label,oldSize,newSize,nbd_p)

  end subroutine reAllocateBigDouble1D
#endif


  !!============================================================================
  !> @brief Allocate, reallocate or deallocate a two-dimensional real array.
  !>
  !> @param[in] label Name of calling subroutine (for memory logging)
  !> @param array The double precision real array to (re)allocate
  !> @param[in] n1 Number of rows in array after reallocation
  !> @param[in] n2 Number of columns in array after reallocation
  !> @param[out] ierr Error flag

  subroutine reAllocateDouble2D (label,array,n1,n2,ierr)

    use KindModule       , only : dp, nbd_p
    use reportErrorModule, only : allocationError

    character(len=*), intent(in)    :: label
    real(dp)        , pointer       :: array(:,:)
    integer,optional, intent(in)    :: n1, n2
    integer,optional, intent(inout) :: ierr

    !! Local variables
    integer :: newSize, oldSize, err

    !! --- Logic section ---

    if (present(n1) .and. present(n2)) then
       newSize = n1*n2
    else
       newSize = 0
    end if
    oldSize = 0

    if (associated(array)) then
       oldSize = size(array)
       if (newSize == oldSize) then
          if (.not.present(n1) .or. .not.present(n2)) return
          if (n1 == size(array,1) .and. n2 == size(array,2)) return
       end if
       !! Release old array
       deallocate(array)
       nullify(array)
    end if

    if (newSize > 0) then
       !! Allocation of new array
       allocate(array(n1,n2),STAT=err)
       if (err /= 0) then
          if (present(ierr)) ierr = ierr + allocationError(label)
          newSize = -newSize
          nullify(array)
       end if
    end if

    if (doLogMem) call logAllocMem (label,oldSize,newSize,nbd_p)

  end subroutine reAllocateDouble2D


  !!============================================================================
  !> @brief Prints the current accumulated memory usage to the log file.
  !>
  !> @param[in] label Name of calling subroutine
  !> @param[in] oldSize Number of words deallocated
  !> @param[in] newSize Number of words allocated
  !> @param[in] nBytes Number of bytes per word

  subroutine logAllocMemI8 (label,oldSize,newSize,nBytes)

    use FileUtilitiesModule, only : findUnitNumber

    character(len=*), intent(in) :: label
    integer(i8)     , intent(in) :: oldSize, newSize
    integer         , intent(in) :: nBytes

    !! Local variables
    integer :: ierr, ip

    !! --- Logic section ---

    if (peakMem == 0_i8) then
       logFile = findUnitNumber(10)
       open(logFile,FILE='memoryUsage.log',IOSTAT=ierr)
       if (ierr /= 0) then
          logFile = 0
          return
       end if
    end if

    totMem  = totMem + (max(0_i8,newSize)-oldSize)*int(nBytes,i8)
    peakMem = max(peakMem,totMem)
    ip = index(label,'::')
    if (ip > 0) then
       ip = ip + 2
    else
       ip = 1
    end if
    if (oldSize > 0 .and. newSize > 0) then
       write(logFile,600) label(ip:), newSize-oldSize, nBytes, totMem, peakMem
    else if (newSize > 0) then
       write(logFile,601) label(ip:), newSize, nBytes, totMem, peakMem
    else if (newSize < 0) then
       write(logFile,603) label(ip:),-newSize, nBytes, totMem
    else if (oldSize > 0) then
       write(logFile,602) label(ip:), oldSize, nBytes, totMem, peakMem
    end if

600 format(A24,': Reallocating array:',I8,' item size:',I4, &
         & ' total:',I12,'  peak:',I12)
601 format(A24,': Allocating array:', I10,' item size:',I4, &
         & ' total:',I12,'  peak:', I12)
602 format(A24,': Deallocating array:',I8,' item size:',I4, &
         & ' total:',I12,'  peak:',I12)
603 format(A24,': Allocation failed:', I9,' item size:',I4,' total:',I12)

  end subroutine logAllocMemI8

#if !defined(FT_HAS_INT4_ONLY) && !defined(FT_USE_INT8)
  !!============================================================================
  !> @brief Prints the current accumulated memory usage to the log file.
  !> @param[in] label Name of calling subroutine
  !> @param[in] oldSize Number of words deallocated
  !> @param[in] newSize Number of words allocated
  !> @param[in] nBytes Number of bytes per word
  subroutine logAllocMemI4 (label,oldSize,newSize,nBytes)
    character(len=*), intent(in) :: label
    integer         , intent(in) :: oldSize, newSize, nBytes
    call logAllocMemI8 (label,int(oldSize,i8),int(newSize,i8),nBytes)
  end subroutine logAllocMemI4
#endif


  !!============================================================================
  !> @brief Closes the memory logging file.
  subroutine closeMemLogFile ()
    if (logFile > 0) then
       close(logFile)
       logFile = 0
    end if
  end subroutine closeMemLogFile


  !!============================================================================
  !> @brief Writes the peak memory usage to file unit @a lpu.
  subroutine writePeakMem (lpu)
    use KindModule, only : dp
    integer, intent(in) :: lpu
    write(lpu,600) real(peakMem,dp) / 1048576.0_dp
600 format(4X,'Peak work memory usage:',F9.2, &
         & ' MB (dynamic allocation in Fortran)')
  end subroutine writePeakMem


  !!============================================================================
  !> @brief Writes storage requirement information for an item to unit lpu.
  !>
  !> @param[in] name Name of the item to write storage info for
  !> @param[in] nwI Number of integer words
  !> @param[in] nwR Number of real words
  !> @param[in] lpu File unit number to write to
  !>
  !> @author Knut Morten Okstad
  !> @date 7 Mar 2003

  subroutine writeStorage (name,nwI,nwR,lpu)

    use KindModule, only : nbi_p, nbd_p

    character(len=*), intent(in) :: name
    integer(i8)     , intent(in) :: nwI, nwR
    integer         , intent(in) :: lpu

    !! Local variables
    character(len=16) :: chBytes

    !! --- Logic section ---

    write(lpu,600) name

    chBytes = StrBytes(nwI*int(nbi_p,i8))
    write(lpu,610) nwI,trim(chBytes)

    chBytes = StrBytes(nwR*int(nbd_p,i8))
    write(lpu,620) nwR,trim(chBytes)

600 format(/4X,'Storage requirement for ',A)
610 format( 7X,'Integer words          :',I12,'  (',A,')')
620 format( 7X,'Double precision words :',I12,'  (',A,')')

  end subroutine writeStorage


  !!============================================================================
  !> @brief Converts a number of bytes to a user-friendly string.
  !>
  !> @param[in] nBytes Number of bytes
  !> @return A string with the number of bytes with a user-friendly byte unit.
  !>
  !> @author Knut Morten Okstad
  !> @date 6 Sep 2006

  function StrBytes (nBytes) result(text)

    integer(i8), intent(in) :: nBytes
    character(len=12)       :: text

    !! --- Logic section ---

#ifdef FT_HAS_INT4_ONLY
    if (nBytes > 1073741824_i8) then
#else
    if (nBytes > 10737418240_i8) then
       write(text,600) (nBytes+536870912_i8)/1073741824_i8,'GBytes'
    else if (nBytes > 1073741824_i8) then
#endif
       write(text,601) real(nBytes)/1073741824.0,'GBytes'
    else if (nBytes > 10485760_i8) then
       write(text,600) (nBytes+524288_i8)/1048576_i8,'MBytes'
    else if (nBytes > 1048576_i8) then
       write(text,601) real(nBytes)/1048576.0,'MBytes'
    else if (nBytes > 10480_i8) then
       write(text,600) (nBytes+512_i8)/1024_i8,'KBytes'
    else if (nBytes > 1048_i8) then
       write(text,601) real(nBytes)/1024.0,'KBytes'
    else
       write(text,600) nBytes,'Bytes'
    end if

600 format(I4,1X,A)
601 format(F4.1,1X,A)

  end function StrBytes

end module allocationModule
