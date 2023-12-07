!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module RoadTypeModule

  use KindModule        , only : dp, lfnam_p
  use FunctionTypeModule, only : FunctionType, EngineType, IdType

  implicit none

  character(len=13), parameter:: roadTypes_p(2) =(/ 'FUNCTION_ROAD', &
       &                                            'SURFACE_ROAD ' /)

  integer, parameter ::                              FUNCTION_ROAD_p = 1, &
       &                                             SURFACE_ROAD_P  = 2

  type RoadType
     type(IdType)                :: id       ! General identification data
     integer                     :: idIn     ! Index into the global road array
     integer                     :: type     ! 1: 2D Func road, 2: 3D Surf road
     type(FunctionType), pointer :: roadFunc ! Function to define the 2D road
     type(EngineType)  , pointer :: xTrEng, yTrEng, zTrEng ! Translation engines
     integer                     :: nRopar
     real(dp), pointer           :: roparr(:)
     integer                     :: nCharRoadDataFileName
     character(len=lfnam_p)      :: RoadDataFileName
  end type RoadType

  interface GetPtrToId
     module procedure GetPtrToIdRoad
  end interface

  interface WriteObject
     module procedure WriteRoadType
  end interface


contains

  subroutine WriteRoadType (road,io)

    !!==========================================================================
    !! Standard routine for writing an object to io.
    !!
    !! Programmer : Bjorn Haugen
    !! date/rev   : Aug 2002/1.0
    !!==========================================================================

    use IdTypeModule, only : writeId

    type(RoadType)  , intent(in) :: road
    integer         , intent(in) :: io

    !! --- Logic section ---

    write(io,'(A)') 'Road','{'
    call writeId (road%id,io)

    write(io,*) 'type        = ', road%type

    if (associated(road%roadFunc)) then
       write(io,*) 'roadFunc(id) = ', road%roadFunc%id%baseId
    end if

    if (associated(road%xTrEng)) then
       write(io,*) 'xTrEng(id)  = ', road%xTrEng%id%baseId
    end if

    if (associated(road%yTrEng)) then
       write(io,*) 'yTrEng(id)  = ', road%yTrEng%id%baseId
    end if

    if (associated(road%zTrEng)) then
       write(io,*) 'zTrEng(id)  = ', road%zTrEng%id%baseId
    end if

    if (associated(road%roparr) .and. road%nRopar > 0) then
       write(io,*) 'roparr      = ', road%roparr(1:road%nRopar)
    end if

    if (road%nCharRoadDataFileName > 0) then
       write(io,*) 'dataFileName = ', trim(road%roadDataFileName)
    end if

    write(io,'(A)') '}'

  end subroutine WriteRoadType


  function GetPtrToIdRoad (array,id) result(ptr)
    type(RoadType), pointer :: ptr

    !!==========================================================================
    !! Return pointer to (first) element with specified id or NULL if not found.
    !!
    !! Programmer : Bjorn Haugen
    !! date/rev   : Aug 2002/1.0
    !!==========================================================================

    integer       , intent(in)        :: id
    type(RoadType), intent(in),target :: array(:)

    !! Local variables
    integer :: i

    !! --- Logic section ---

    do i = 1, size(array)
       if (array(i)%id%baseId == id) then
          ptr => array(i)
          return
       end if
    end do

    nullify(ptr)
    write(*,*) '*** GetPtrToIdRoad returned nullified, baseId =',id

  end function GetPtrToIdRoad


  subroutine NullifyRoad (road)

    !!==========================================================================
    !! Initialize the road type object.
    !!
    !! Programmer : Trond Arne Svidal
    !! date/rev   : Aug 2002/1.0
    !!==========================================================================

    use IdTypeModule, only : nullifyId

    type(RoadType), intent(out) :: road

    !! --- Logic section ---

    call nullifyId (road%id)

    road%idIn = 0
    road%type = 0
    nullify(road%roadFunc)
    nullify(road%xTrEng)
    nullify(road%yTrEng)
    nullify(road%zTrEng)
    nullify(road%roparr)
    road%nRopar = 0
    road%nCharRoadDataFileName = 0
    road%RoadDataFileName = ''

  end subroutine NullifyRoad


  subroutine DeallocateRoad (road)

    !!==========================================================================
    !! Deallocate the RoadType object.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 23 Jan 2017/1.0
    !!==========================================================================

    use IdTypeModule, only : deallocateId

    type(RoadType), intent(inout) :: road

    !! --- Logic section ---

    call deallocateId (road%id)

    if (associated(road%roparr)) deallocate(road%roparr)

    call nullifyRoad (road)

  end subroutine DeallocateRoad


  subroutine deallocateRoads (roads)
    type(RoadType), pointer :: roads(:)
    integer :: i
    do i = 1, size(roads)
       call DeallocateRoad (roads(i))
    end do
    deallocate(roads)
    nullify(roads)
  end subroutine deallocateRoads

end module RoadTypeModule
