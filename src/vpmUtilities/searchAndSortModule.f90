!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file searchAndSortModule.f90
!> @brief Utilities for searching and sorting of arrays.

!!==============================================================================
!> @brief Module with subroutines for searching and sorting of arrays.

module SearchAndSortModule

  implicit none

  private :: quickSortIntArray


contains

  !!============================================================================
  !> @brief Finds indices of the largest values in the given real array.
  !>
  !> @author Knut Morten Okstad
  !> @date Jan 2004

  subroutine findMaxRealValues (values,indMax,valMax,ierr)

    use KindModule        , only : dp
    use ScratchArrayModule, only : getLogicalScratchArray
    use ReportErrorModule , only : reportError, debugFileOnly_p

    real(dp)         , intent(in)  :: values(:)
    integer          , intent(out) :: indMax(:)
    real(dp),optional, intent(out) :: valMax(:)
    integer, optional, intent(out) :: ierr

    !! Local variables
    integer          :: i, loc(1)
    logical, pointer :: mask(:)

    !! --- Logic section ---

    mask => getLogicalScratchArray(size(values),i)
    if (i < 0) then
       if (present(ierr)) ierr = -1
       call reportError (debugFileOnly_p,'findMaxValues')
       return
    end if

    mask = .true.
    do i = 1, size(indMax)
       loc = maxloc(values,mask)
       mask(loc(1)) = .false.
       indMax(i) = loc(1)
       if (present(valMax)) then
          if (i <= size(valMax)) valMax(i) = values(loc(1))
       end if
    end do

  end subroutine findMaxRealValues


  !!============================================================================
  !> @brief Sorts an integer array in increasing order.
  !> @details The actual permutation is returned.
  !>
  !> @author Knut Morten Okstad
  !> @date Mar 2003

  subroutine quickSortInt (values,indval,ierr)

    use ScratchArrayModule, only : getIntegerScratchArray
    use ReportErrorModule , only : reportError, debugFileOnly_p

    integer                  , intent(inout) :: values(:)
    integer, optional, target, intent(out)   :: indval(:)
    integer, optional        , intent(out)   :: ierr

    !! Local variables
    integer          :: i, n
    integer, pointer :: indices(:)

    !! --- Logic section ---

    n = size(values)
    if (present(indval)) then
       indices => indval
       if (size(indices) < n) return
    else
       indices => getIntegerScratchArray(n,i)
       if (i < 0) then
          if (present(ierr)) ierr = -1
          call reportError (debugFileOnly_p,'quickSortInt')
          return
       end if
    end if

    do i = 1, n
       indices(i) = i
    end do

    call quickSortIntArray (indices,values,1,n)

  end subroutine quickSortInt


  !!============================================================================
  !> @brief Sorts an integer array in increasing order.
  !>
  !> @author Bjorn Haugen
  !> @date Nov 1998

  recursive subroutine quickSortIntArray (indexArray,valueArray,lIndex,rIndex)

    integer, intent(inout) :: indexArray(:), valueArray(:)
    integer, intent(in)    :: lIndex, rIndex

    !! Local variables
    integer :: l1, r1, tmp

    !! --- Logic section ---

    if (lIndex >= rIndex) return

    l1 = lIndex
    r1 = rIndex

    do while (l1 < r1)
       do while (l1 < rIndex .and. valueArray(l1) <= valueArray(lIndex))
          l1 = l1 + 1
       end do

       do while (r1 > lIndex .and. valueArray(r1) >= valueArray(lIndex))
          r1 = r1 - 1
       end do

       if (l1 < r1) then
          call swapInt (valueArray(l1), valueArray(r1))
          call swapInt (indexArray(l1), indexArray(r1))
       end if
    end do

    call swapInt (valueArray(lIndex), valueArray(r1))
    call swapInt (indexArray(lIndex), indexArray(r1))

    call quickSortIntArray (indexArray, valueArray, lIndex, r1-1)
    call quickSortIntArray (indexArray, valueArray, r1+1, rIndex)

  contains

    !> @brief Swaps the contents of two integer variables.
    subroutine swapInt (a,b)
      integer, intent(inout) :: a, b
      tmp = a
      a = b
      b = tmp
    end subroutine swapInt

  end subroutine quickSortIntArray

end module SearchAndSortModule

