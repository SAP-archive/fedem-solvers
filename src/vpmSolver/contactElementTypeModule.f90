!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file contactElementTypeModule.f90
!> @brief Contact element object data container.

!!==============================================================================
!> @brief Module with data types representing contact element objects.
!>
!> @details The module also contains subroutines for accessing or manipulating
!> the contact element data.

module ContactElementTypeModule

  use IdTypeModule        , only : IdType, getId
  use TriadTypeModule     , only : TriadType, dp
  use SpringTypeModule    , only : SpringPtrType
  use DamperTypeModule    , only : DamperPtrType
  use FrictionTypeModule  , only : FrictionType
  use ContactSurfaceModule, only : GliderCurveType

  implicit none

  !> @brief Data type representing a contact element object.
  type ContactElementType

     type(IdType) :: id !< General identification data

     integer :: samElNum !< Element number for SAM reference (MPMNPC)
     integer :: nDOFs    !< Number of element DOFs (sum of triad DOFs)
     integer :: nWarning !< Number of warnings detected for this element

     type(TriadType)      , pointer :: triad1 !< Follower triad
     type(GliderCurveType), pointer :: csurf  !< Contact surface definition

     type(SpringPtrType)         :: springs(6) !< Contact springs
     type(DamperPtrType)         :: dampers(6) !< Contact dampers
     type(FrictionType), pointer :: friction   !< Slider friction

     real(dp) :: thickness      !< Contact surface domain thickness
     real(dp) :: width          !< Contact surface domain width
     real(dp) :: radius         !< Contact surface radius (radial contact)
     real(dp) :: FPosInT(3,4)   !< Relative follower position
     real(dp) :: cVar(6,2)      !< Movement of follower w.r.t. contact curve
     real(dp) :: cVarPrev(3)    !< Contact variables at previous time step
     real(dp) :: accContactDist !< Accumulated sliding while in contact
     logical  :: saveVar(2)     !< Flags indicating variables to be saved

     logical  :: wasActive  !< Indicates whether follower was in contact
     logical  :: isActive   !< Indicates whether follower is in contact
     logical  :: isAtCorner !< Indicates whether contact is at a corner point
     integer  :: lastDomain !< Last curve domain follower was in contact with

     real(dp) :: coeff(2)     !< Scaling for the active curve segment
     real(dp) :: CPosInG(3,4) !< Global position of contact point
     real(dp) :: eVec(3,2)    !< Vectors from active triads to contact point

     logical  :: nonRadialDamping      !< Radial contact but rectangular dampers
     real(dp), pointer :: DPosInG(:,:) !< Global position matrix for dampers

  end type ContactElementType


  logical, save :: ignoreEccF = .false. !< Ignore eccentricity Force terms
  logical, save :: ignoreEccN = .false. !< Ignore eccentricity Newton terms


  !> @brief Returns pointer to object with specified ID.
  interface GetPtrToId
     module procedure GetPtrToIdContactElm
  end interface

  !> @brief Standard routine for writing an object to file.
  interface WriteObject
     module procedure WriteContactElementType
  end interface

  !> @brief Updates contact element quantities at convergence.
  interface UpdateAtConvergence
     module procedure UpdateDistanceInContact
  end interface


contains

  !!============================================================================
  !> @brief Returns pointer to (first) contact element with specified ID.
  !>
  !> @param[in] array Array of contactelementtypemodule::contactelementtype
  !>                  objects to search within
  !> @param[in] id Base ID of the object to search for
  !> @param[out] index The array index of the found object
  !>
  !> @details If the contact element is not found, NULL is returned.
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 1 Nov 1999

  function GetPtrToIdContactElm (array,id,index) result(ptr)

    integer                 , intent(in)            :: id
    type(ContactElementType), intent(in) , target   :: array(:)
    integer                 , intent(out), optional :: index
    type(ContactElementType), pointer               :: ptr

    !! Local variables
    integer :: i

    !! --- Logic section ---

    do i = 1, size(array)
       if (array(i)%id%baseId == id) then
          ptr => array(i)
          if (present(index)) index = i
          return
       end if
    end do

    nullify(ptr)
    write(*,*) '*** GetPtrToIdContactElm returned nullified, baseId =',id

  end function GetPtrToIdContactElm


  !!============================================================================
  !> @brief Standard routine for writing an object to io.
  !>
  !> @param[in] contElem The contactelementtypemodule::contactelementtype
  !>                     object to write
  !> @param[in] io File unit number to write to
  !> @param[in] complexity If present, the value indicates the amount of print
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 1 Sep 2001

  subroutine WriteContactElementType (contElem,io,complexity)

    use IdTypeModule        , only : writeId
    use SpringTypeModule    , only : writeObject
    use DamperTypeModule    , only : writeObject
    use ContactSurfaceModule, only : writeObject

    type(ContactElementType), intent(in) :: contElem
    integer                 , intent(in) :: io
    integer,optional        , intent(in) :: complexity

    !! Local variables
    integer :: i

    !! --- Logic section ---

    write(io,'(A)') 'ContactElement','{'
    call writeId (contElem%id,io)

    write(io,*) 'samElNum    = ', contElem%samElNum
    write(io,*) 'nDOFs       = ', contElem%nDOFs
    write(io,*) 'nWarning    = ', contElem%nWarning

    if (associated(contElem%triad1)) then
       write(io,*) 'followerId  = ', contElem%triad1%id%baseId
    end if
    if (associated(contElem%csurf)) then
       call writeObject (contElem%csurf,io)
    end if

    write(io,*) 'thickness   = ', contElem%thickness
    write(io,*) 'width       = ', contElem%width
    write(io,*) 'radius      = ', contElem%radius
    write(io,*) 'lastDomain  = ', contElem%lastDomain

    if (associated(contElem%friction)) then
       if (associated(contElem%friction%param)) then
          write(io,*) 'frictionId  = ',contElem%friction%param%id%baseId
       end if
    end if

    if (present(complexity)) then
       if (complexity >= 1) then

          do i = 1, min(size(contElem%springs),size(contElem%dampers))
             if (associated(contElem%springs(i)%p)) then
                call writeObject (contElem%springs(i)%p,io,complexity)
             end if
             if (associated(contElem%dampers(i)%p)) then
                call writeObject (contElem%dampers(i)%p,io,complexity)
             end if
          end do

       end if
       if (complexity >= 2) then

          write(io,*) 'deflection  = ', contElem%cVar(:,1)
          write(io,*) 'defl_prev   = ', contElem%cVarPrev
          write(io,*) 'velocity    = ', contElem%cVar(:,2)
          write(io,*) 'accContact  = ', contElem%accContactDist

       end if
    end if

    write(io,'(A)') '}'

  end subroutine WriteContactElementType


  !!============================================================================
  !> @brief Initializes a contact element type object.
  !>
  !> @param contElem The contactelementtypemodule::contactelementtype
  !>                 object to initialize
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> date 24 Jun 2002

  pure subroutine NullifyContactElement (contElem)

    use IdTypeModule, only : nullifyId

    type(ContactElementType), intent(out) :: contElem

    !! Local variables
    integer :: i

    !! --- Logic section ---

    call nullifyId (contElem%id)

    contElem%samElNum = 0
    contElem%nDOFs = 0
    contElem%nWarning = 0

    nullify(contElem%triad1)
    nullify(contElem%csurf)
    do i = 1, size(contElem%springs)
       nullify(contElem%springs(i)%p)
    end do
    do i = 1, size(contElem%dampers)
       nullify(contElem%dampers(i)%p)
    end do
    nullify(contElem%friction)

    contElem%thickness = 0.0_dp
    contElem%width = 0.0_dp
    contElem%radius = 0.0_dp
    contElem%cVar = 0.0_dp
    contElem%cVarPrev = 0.0_dp
    contElem%accContactDist = 0.0_dp
    contElem%saveVar = .false.

    contElem%wasActive = .false.
    contElem%isActive = .false.
    contElem%isAtCorner = .false.
    contElem%lastDomain = 0

    contElem%coeff = 0.0_dp
    contElem%FPosInT = 0.0_dp
    contElem%CPosInG = 0.0_dp
    contElem%eVec = 0.0_dp

    contElem%nonRadialDamping = .false.
    nullify(contElem%DPosInG)

  end subroutine NullifyContactElement


  !!============================================================================
  !> @brief Deallocates a contact element object.
  !>
  !> @param contElem The contactelementtypemodule::contactelementtype
  !>                 object to deallocate
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> date 23 Jan 2017

  elemental subroutine DeallocateContactElement (contElem)

    use IdTypeModule, only : deallocateId

    type(ContactElementType), intent(inout) :: contElem

    !! --- Logic section ---

    call deallocateId (contElem%id)

    if (associated(contElem%friction)) deallocate(contElem%friction)
    if (contElem%nonRadialDamping .and. associated(contElem%DPosInG)) then
       deallocate(contElem%DPosInG)
    end if

    call nullifyContactElement (contElem)

  end subroutine DeallocateContactElement


  !!============================================================================
  !> @brief Deallocates an array of contact element objects.
  !>
  !> @param contElms The contactelementtypemodule::contactelementtype
  !>                 objects to deallocate
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> date 23 Jan 2017

  subroutine deallocateContactElements (contElms)

    type(ContactElementType), pointer :: contElms(:)

    !! --- Logic section ---

    call DeallocateContactElement (contElms)

    deallocate(contElms)
    nullify(contElms)

  end subroutine deallocateContactElements


  !!============================================================================
  !> @brief Updates contact element quantities at convergence.
  !>
  !> @param contElem The contactelementtypemodule::contactelementtype
  !>                 object to update
  !>
  !> @details This subroutine updates the accumulated sliding distance while
  !> contact is established. It should be invoked only once for each time step,
  !> after convergence has been achieved.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 3 Oct 2005

  subroutine UpdateDistanceInContact (contElem)

    type(ContactElementType), intent(out) :: contElem

    !! Local variables
    integer  :: i
    logical  :: wasInContact(2), isInContact(2)
    real(dp) :: deltaDis
    real(dp), parameter :: eps_p = 1.0e-16_dp

    !! --- Logic section ---

    if (contElem%isActive .and. contElem%wasActive) then

       isInContact = .false.
       wasInContact = .false.
       do i = 1, 2
          if (associated(contElem%springs(i)%p)) then
             isInContact(i) = abs(contElem%springs(i)%p%force) > eps_p
             wasInContact(i) = abs(contElem%springs(i)%p%forcePrev) > eps_p
          end if
       end do

       deltaDis = abs(contElem%cVar(3,1) - contElem%cVarPrev(3))
       if (any(isInContact) .and. any(wasInContact)) then
          !! Contact is established both in current and previous time step
          contElem%accContactDist = contElem%accContactDist + deltaDis
       else if (any(isInContact) .or. any(wasInContact)) then
          !! Contact status changed in this time step, assume halfways contact
          contElem%accContactDist = contElem%accContactDist + 0.5_dp*deltaDis
       end if

    end if

  end subroutine UpdateDistanceInContact

end module ContactElementTypeModule
