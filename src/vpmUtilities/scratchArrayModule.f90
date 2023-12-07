!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file scratchArrayModule.f90
!> @brief Utilities for dynamic allocation of scratch arrays.

!!==============================================================================
!> @brief Module with subroutines for dynamic allocation of scratch arrays.
!> @details This module is provided to avoid frequent allocation and
!> deallocation of temporary scratch arrays during the computations.
!> It should NOT be used to transfer data between subroutines.
!> That will be dangerous and may sooner or later lead to overwriting errors.
!>
!> @brief Knut Morten Okstad
!> @date 12 July 2002

module scratchArrayModule

  use kindModule, only : sp, dp

  implicit none

  integer , save, private, allocatable, target :: iscr(:) !< Integer array
  logical , save, private, allocatable, target :: lscr(:) !< Logical array
  real(sp), save, private, allocatable, target :: sscr(:) !< Single precision
  real(dp), save, private, allocatable, target :: dscr(:) !< Double precision

  !> @brief Returns a pointer to a real scratch array.
  interface realScratchArray
     module procedure singleScratchArray
     module procedure doubleScratchArray
  end interface

  private :: allocateIscr, allocateLscr, allocateSscr, allocateDscr
  private :: singleScratchArray, doubleScratchArray

contains

  !!============================================================================
  !> @brief Reallocates the integer scratch array to fit at least the size.
  subroutine allocateIscr (nw,ierr)

    use kindModule       , only : nbi_p
    use allocationModule , only : doLogMem, logAllocMem
    use reportErrorModule, only : allocationError

    integer, intent(in)  :: nw
    integer, intent(out) :: ierr

    !! Local variables
    integer :: oldSize

    !! --- Logic section ---

    ierr = 0
    if (.not. allocated(iscr)) then
       oldSize = 0
       allocate(iscr(nw),STAT=ierr)
    else if (size(iscr) < nw) then
       oldSize = size(iscr)
       deallocate(iscr)
       allocate(iscr(nw),STAT=ierr)
    else
       return
    end if

    if (ierr /= 0) then
       ierr = allocationError('ScratchArrayModule::allocateIscr')
       if (doLogMem) call logAllocMem ('allocateIscr',oldSize,-nw,nbi_p)
    else if (doLogMem) then
       call logAllocMem ('allocateIscr',oldSize,nw,nbi_p)
    end if

  end subroutine allocateIscr

  !!============================================================================
  !> @brief Reallocates the logical scratch array to fit at least given size.
  !>
  !> @author Knut Morten Okstad
  !> @date 29 Jan 2004

  subroutine allocateLscr (nw,ierr)

    use kindModule       , only : nbi_p
    use allocationModule , only : doLogMem, logAllocMem
    use reportErrorModule, only : allocationError

    integer, intent(in)  :: nw
    integer, intent(out) :: ierr

    !! Local variables
    integer :: oldSize

    !! --- Logic section ---

    ierr = 0
    if (.not. allocated(lscr)) then
       oldSize = 0
       allocate(lscr(nw),STAT=ierr)
    else if (size(lscr) < nw) then
       oldSize = size(lscr)
       deallocate(lscr)
       allocate(lscr(nw),STAT=ierr)
    else
       return
    end if

    if (ierr /= 0) then
       ierr = allocationError('ScratchArrayModule::allocateLscr')
       if (doLogMem) call logAllocMem ('allocateLscr',oldSize,-nw,nbi_p)
    else if (doLogMem) then
       call logAllocMem ('allocateLscr',oldSize,nw,nbi_p)
    end if

  end subroutine allocateLscr


  !!============================================================================
  !> @brief Reallocates the real(sp) scratch array to fit at least given size.
  !>
  !> @author Knut Morten Okstad
  !> @date 20 Jan 2006

  subroutine allocateSscr (nw,ierr)

    use kindModule       , only : nbs_p
    use allocationModule , only : doLogMem, logAllocMem
    use reportErrorModule, only : allocationError

    integer, intent(in)  :: nw
    integer, intent(out) :: ierr

    !! Local variables
    integer :: oldSize

    !! --- Logic section ---

    if (.not. allocated(sscr)) then
       oldSize = 0
       allocate(sscr(nw),STAT=ierr)
    else if (size(sscr) < nw .or. nw == 0) then
       oldSize = size(sscr)
       deallocate(sscr)
       allocate(sscr(nw),STAT=ierr)
    else
       ierr = 0
       return
    end if

    if (ierr /= 0) then
       ierr = allocationError('ScratchArrayModule::allocateSscr')
       if (doLogMem) call logAllocMem ('allocateSscr',oldSize,-nw,nbs_p)
    else if (doLogMem) then
       call logAllocMem ('allocateSscr',oldSize,nw,nbs_p)
    end if

  end subroutine allocateSscr


  !!============================================================================
  !> @brief Reallocates the real(dp) scratch array to fit at least given size.

  subroutine allocateDscr (nw,ierr)

    use kindModule       , only : nbd_p
    use allocationModule , only : doLogMem, logAllocMem
    use reportErrorModule, only : allocationError

    integer, intent(in)  :: nw
    integer, intent(out) :: ierr

    !! Local variables
    integer :: oldSize

    !! --- Logic section ---

    if (.not. allocated(dscr)) then
       oldSize = 0
       allocate(dscr(nw),STAT=ierr)
    else if (size(dscr) < nw .or. nw == 0) then
       oldSize = size(dscr)
       deallocate(dscr)
       allocate(dscr(nw),STAT=ierr)
    else
       ierr = 0
       return
    end if

    if (ierr /= 0) then
       ierr = allocationError('ScratchArrayModule::allocateDscr')
       if (doLogMem) call logAllocMem ('allocateDscr',oldSize,-nw,nbd_p)
    else if (doLogMem) then
       call logAllocMem ('allocateDscr',oldSize,nw,nbd_p)
    end if

  end subroutine allocateDscr


  !!==========================================================================
  !> @brief Returns a pointer to an integer scratch array.

  function getIntegerScratchArray (nsize,ierr)

    integer, pointer     :: getIntegerScratchArray(:)
    integer, intent(in)  :: nsize
    integer, intent(out) :: ierr

    !! --- Logic section ---

    call allocateIscr (nsize,ierr)
    getIntegerScratchArray => iscr(1:nsize)

  end function getIntegerScratchArray


  !!==========================================================================
  !> @brief Returns a pointer to a logical scratch array.
  !>
  !> @author Knut Morten Okstad
  !> @date 29 Jan 2004

  function getLogicalScratchArray (nsize,ierr)

    logical, pointer     :: getLogicalScratchArray(:)
    integer, intent(in)  :: nsize
    integer, intent(out) :: ierr

    !! --- Logic section ---

    call allocateLscr (nsize,ierr)
    getLogicalScratchArray => lscr(1:nsize)

  end function getLogicalScratchArray


  !!============================================================================
  !> @brief Returns a pointer to a real(sp) scratch array.

  function getSingleScratchArray (nsize,ierr)

    real(sp), pointer     :: getSingleScratchArray(:)
    integer , intent(in)  :: nsize
    integer , intent(out) :: ierr

    !! --- Logic section ---

    call allocateSscr (nsize,ierr)
    getSingleScratchArray => sscr(1:nsize)

  end function getSingleScratchArray


  !!==========================================================================
  !> @brief Returns a pointer to a real(dp) scratch array.

  function getRealScratchArray (nsize,ierr)

    real(dp), pointer     :: getRealScratchArray(:)
    integer , intent(in)  :: nsize
    integer , intent(out) :: ierr

    !! --- Logic section ---

    call allocateDscr (nsize,ierr)
    getRealScratchArray => dscr(1:nsize)

  end function getRealScratchArray


  !!==========================================================================
  !> @brief Returns a pointer to a `nrow`&times;`ncol` real(dp) scratch array.

  function getRealScratchMatrix (nrow,ncol,ierr)

    use manipMatrixModule, only : MatrixPointToArray_real

    real(dp), pointer     :: getRealScratchMatrix(:,:)
    integer , intent(in)  :: nrow, ncol
    integer , intent(out) :: ierr

    !! --- Logic section ---

    call allocateDscr (nrow*ncol,ierr)
    getRealScratchMatrix => MatrixPointToArray_real(dscr,nrow,ncol)

  end function getRealScratchMatrix


  !!============================================================================
  !> @brief Returns a pointer to a real(sp) scratch array.
  !>
  !> @author Knut Morten Okstad
  !> @date 10 Feb 2017

  subroutine singleScratchArray (rscr,nsize,ierr)

    real(sp), pointer     :: rscr(:)
    integer , intent(in)  :: nsize
    integer , intent(out) :: ierr

    !! --- Logic section ---

    call allocateSscr (nsize,ierr)
    rscr => sscr(1:nsize)

  end subroutine singleScratchArray


  !!============================================================================
  !> @brief Returns a pointer to a real(dp) scratch array.
  !>
  !> @author Knut Morten Okstad
  !> @date 10 Feb 2017

  subroutine doubleScratchArray (rscr,nsize,ierr)

    real(dp), pointer     :: rscr(:)
    integer , intent(in)  :: nsize
    integer , intent(out) :: ierr

    !! --- Logic section ---

    call allocateDscr (nsize,ierr)
    rscr => dscr(1:nsize)

  end subroutine doubleScratchArray


  !!============================================================================
  !> @brief Deallocates the scratch arrays, if allocated.
  !>
  !> @author Knut Morten Okstad
  !> @date 15 Jan 2003

  subroutine releaseScratchArrays ()

    use kindModule      , only : nbi_p, nbs_p, nbd_p
    use allocationModule, only : doLogMem, logAllocMem

    !! Local variables
    integer :: oldSize

    !! --- Logic section ---

    if (allocated(iscr)) then
       oldSize = size(iscr)
       deallocate(iscr)
       if (doLogMem) call logAllocMem ('releaseScratchArrays',oldSize,0,nbi_p)
    end if

    if (allocated(lscr)) then
       oldSize = size(lscr)
       deallocate(lscr)
       if (doLogMem) call logAllocMem ('releaseScratchArrays',oldSize,0,nbi_p)
    end if

    if (allocated(sscr)) then
       oldSize = size(sscr)
       deallocate(sscr)
       if (doLogMem) call logAllocMem ('releaseScratchArrays',oldSize,0,nbs_p)
    end if

    if (allocated(dscr)) then
       oldSize = size(dscr)
       deallocate(dscr)
       if (doLogMem) call logAllocMem ('releaseScratchArrays',oldSize,0,nbd_p)
    end if

  end subroutine releaseScratchArrays

end module scratchArrayModule
