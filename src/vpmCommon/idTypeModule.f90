!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file idTypeModule.f90
!> @brief Data type for object identification with associated subroutines.

!!==============================================================================
!> @brief Module with a data type for object identification.
!> @details This module contains a data type that is used by all mechanism
!> objects to store the unique identification of each object. It also has
!> associated subroutines and functions for manipulation of the id data.
!> You have to configure doxygen with the option ENABLED_SECTIONS = FULL_DOC
!> to extract detailed documentation of those subroutines and functions.

module IdTypeModule

  implicit none

  integer, parameter :: ldesc_p = 64       !< Max. length of description tag
  integer, parameter :: lId_p = ldesc_p+13 !< Length of string returned by GetId

  !> @brief Data type for unique object identification
  type IdType
     integer                :: baseId   !< The unique base ID of the object
     integer                :: userId   !< The non-unique user ID of the object
     integer, pointer       :: assId(:) !< Assembly path (user ID) of the object
     character(len=ldesc_p) :: descr    !< Description tag of the object
  end type IdType


contains
  !> @cond FULL_DOC

  !> @brief Nullifies the @a id object.
  pure subroutine NullifyId(id)
    type(IdType), intent(out) :: id

    id%baseId = 0
    id%userId = 0
    id%descr  = ''
    nullify(id%assId)

  end subroutine NullifyId


  !> @brief Deallocates the @a id object.
  pure subroutine DeallocateId(id)
    type(IdType), intent(inout) :: id

    if (associated(id%assId)) then
       deallocate(id%assId)
       nullify(id%assId)
    end if

  end subroutine DeallocateId


  !> @brief Initializes the @a id object.
  pure subroutine InitId(id,baseId,userId,descr,stat)
    type(IdType)    , intent(out) :: id
    integer         , intent(in)  :: baseId, userId(:)
    character(len=*), intent(in)  :: descr
    integer         , intent(out) :: stat
    integer :: n

    id%baseId = baseId
    id%userId = userId(1)
    id%descr  = descr

    n = 1
    do while (n < size(userId) .and. userId(n+1) /= 0)
       n = n + 1
    end do
    if (n > 1) then
       allocate(id%assId(n-1),STAT=stat)
       if (stat == 0) id%assId = userId(2:n)
    else
       nullify(id%assId)
       stat = 0
    end if

  end subroutine InitId


  !> @brief Writes the @a id object to file unit @a lpu.
  subroutine WriteId(id,lpu)
    type(IdType), intent(in) :: id
    integer     , intent(in) :: lpu

    write(lpu,*) 'baseId      =', id%baseId
    write(lpu,*) 'userId      =', id%userId
    if (associated(id%assId)) write(lpu,*) 'assId       =', id%assId
    if (id%descr /= '') write(lpu,*) 'description = ',id%descr

  end subroutine WriteId


  !> @brief Writes the @a id object to a result file header.
  subroutine WriteIdHeader(type,id,lpu,doAdvance)
    character(len=*), intent(in) :: type
    type(IdType)    , intent(in) :: id
    integer         , intent(in) :: lpu
    logical,optional, intent(in) :: doAdvance
    character(len=10) :: txtId
    character(len=3)  :: adv

    if (.not.present(doAdvance)) then
       adv = 'NO'
    else if (doAdvance) then
       adv = 'YES'
    else
       adv = 'NO'
    end if

    write(lpu,600,advance='NO') trim(type)
    if (id%baseId > 0) then
       write(txtId,"(i10)") id%baseId
       write(lpu,"(a,';')",advance='NO') trim(adjustl(txtId))
    else
       write(lpu,"(';')",advance='NO')
    end if

    if (id%userId > 0) then
       write(txtId,"(i10)") id%userId
       write(lpu,"(a,';')",advance='NO') trim(adjustl(txtId))
    else
       write(lpu,"(';')",advance='NO')
    end if

    if (id%descr == '') then
       write(lpu,"(';')",advance=adv)
    else
       write(lpu,601) trim(id%descr)
    end if

600 format('{"',a,'";')
601 format('"',a,'";')

  end subroutine WriteIdHeader


  !> @brief Converts the @a id number into a string (negative sign is ignored).
  function GetId(id) result(text)
    type(IdType), intent(in) :: id
    character(len=lId_p) :: text
    integer :: myId, i, j, k

    !! We cannot use formatted write into a string here because that
    !! will result in a recursive IO crash on IBM aix when this function
    !! is used in a write statement in the calling subroutine/function.
    if (id%descr /= '' .or. id%userId /= 0) then
       myId = id%userId
    else
       myId = id%baseId
    end if
    i = 10
    text(1:10) = ' '
    do while (myId /= 0 .and. i > 0)
       text(i:i) = char(ichar('0')+mod(abs(myId),10))
       myId = myId/10
       i = i - 1
    end do

    if (associated(id%assId)) then
       !! Append the assembly Id path as a comma-separated list
       j = lId_p
       do k = 1, size(id%assId)
          myId = id%assId(k)
          do while (myId /= 0 .and. j > 11)
             text(j:j) = char(ichar('0')+mod(abs(myId),10))
             myId = myId/10
             j = j - 1
          end do
          text(j:j) = ','
          j = j - 1
       end do
       if (j > 10) text = text(1:10)//text(j+1:)
       j = 10 + lId_p - j
    else
       j = 10
    end if

    !! Add the description, if any
    if (i < 1) i = 1
    if (id%descr /= '') then
       text = text(i:j)//' "'//trim(id%descr)//'"'
    else if (id%userId /= 0) then
       text = text(i:j)
    else
       text = text(i:j)//' (baseId)'
    end if

  end function GetId


  !> @brief Converts the assembly path of an @a id object to a string
  function GetAssId(id) result(text)
    type(IdType), intent(in) :: id
    character(len=lId_p) :: text
    integer :: myId, i, j

    i = lId_p
    if (associated(id%assId)) then
       !! Build the assembly Id path as a comma-separated list
       do j = 1, size(id%assId)
          myId = id%assId(j)
          do while (myId /= 0)
             text(i:i) = char(ichar('0')+mod(abs(myId),10))
             if (i < 2) return
             myId = myId/10
             i = i - 1
          end do
          text(i:i) = ','
          i = i - 1
       end do
    end if
    if (i+1 < lId_p) then
       text = text(i+2:)
    else
       text = ' '
    end if

  end function GetAssId


  !> @brief Converts the @a id number into a string (negative sign is ignored).
  function StrId(id) result(text)
    integer, intent(in) :: id
    character(len=10) :: text
    integer :: myId, indx

    !! We cannot use formatted write into a string here because that
    !! will result in a recursive IO crash on IBM aix if this function
    !! is used in a write statement in the calling subroutine/function.
    text = ' '
    myId = id
    indx = 10
    do while (myId /= 0 .and. indx > 0)
       text(indx:indx) = char(ichar('0')+mod(abs(myId),10))
       myId = myId/10
       indx = indx - 1
    end do
    if (id == 0) then
       text = ' 0'
    else if (indx > 1) then
       text = text(indx:)
    end if

  end function StrId


  !> @brief Prints an error message related to solver input file parsing.
  subroutine ReportInputError(name,idIn,id,msg)
    use reportErrorModule, only : reportError, error_p
    character(*)          , intent(in) :: name
    integer     , optional, intent(in) :: idIn
    type(IdType), optional, intent(in) :: id
    character(*), optional, intent(in) :: msg
    character(len=128) :: errMsg

    if (present(id)) then
       write(errMsg,"('Invalid specification for ',A,A)") name,getId(id)
    else if (present(idIn)) then
       write(errMsg,"('Could not read ',A,' namelist number',I6)") name,idIn
    else
       write(errMsg,"('Could not read ',A,' namelist')") name
    end if

    if (present(msg)) then
       call reportError (error_p,errMsg,msg)
    else
       call reportError (error_p,errMsg)
    end if

  end subroutine ReportInputError

  !> @endcond
end module IdTypeModule
