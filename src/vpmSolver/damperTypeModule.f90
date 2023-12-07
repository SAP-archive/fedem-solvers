!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file damperTypeModule.f90
!> @brief Damper object data containers.

!!==============================================================================
!> @brief Module with data types representing damper objects.
!>
!> @details The module also contains subroutines for accessing the damper data.

module DamperTypeModule

  use TriadTypeModule   , only : TriadType, IdType, dp
  use FunctionTypeModule, only : FunctionType, EngineType
  use SpringTypeModule  , only : SpringBaseType

  implicit none

  !> @brief Data type representing a base damper object.
  type DamperBaseType

     type(IdType) :: id !< General identification data

     logical :: isActive !< Turn damper ON and OFF
     integer :: dof      !< Joint variable the damper acts on, 0 if axial damper

     type(FunctionType), pointer :: coeffFunction !< Direct function of velocity
     type(FunctionType), pointer :: forceFunction !< Force-velocity function

     real(dp) :: coeff !< Current damper coefficient
     real(dp) :: dmp0  !< Constant damping coefficient value
     real(dp) :: dmp1  !< Scaling factor for damping coefficient/force function

     !> Optional scaling of the damping coefficient
     type(EngineType), pointer :: coeffScaleEngine
     real(dp) :: scale1 !< current scaling coefficient = scale1*coeffScaleEngine

     !> Associated spring for deformational damper
     type(SpringBaseType), pointer :: spr

     !> @brief Damper length.
     !> @details For joint dampers, this pointer points to the
     !> corresponding joint variable. For axial dampers it is allocated.
     real(dp), pointer :: length

     !> @brief Damper velocity.
     !> @details For regular joint dampers, this pointer points to the velocity
     !> of the corresponding joint variable. For deformational joint dampers
     !> and axial dampers it is allocated.
     real(dp), pointer :: velocity

     real(dp) :: lengthPrev !< Damper length at previous time step
     real(dp) :: force      !< Current damper force
     real(dp) :: forcePrev  !< Previous damper force
     real(dp) :: Edmp       !< Current energy loss

     logical :: saveVar(5) !< Flags indicating which variables should be saved

  end type DamperBaseType


  !> @brief Data type representing a damper element.
  type DamperType

     type(DamperBaseType), pointer :: dmp !< Base damper data

     integer :: samElNum !< Element number for SAM reference (index into MPMNPC)

     !> @brief Number of DOFs for the damper element.
     !> @details Equals @a nJointDOFs for joint dampers,
     !> 6 or 12 for axial dampers, and 3 or 6 for dampers to ground.
     integer :: nDofs

     type(TriadType), pointer :: triad1 !< First triad for axial damper
     type(TriadType), pointer :: triad2 !< Second triad for axial damper

     real(dp), pointer :: forceDir(:) !< Force direction for axial dampers

  end type DamperType


  !> @brief Data type representing a base damper pointer.
  !> @details This data type is used to construct arrays of dampers where each
  !> element is a pointer to a damper object, and not the objects themselves.
  type DamperPtrType
     type(DamperBaseType), pointer :: p !< Pointer to a base damper object
  end type DamperPtrType


  !> @brief Returns pointer to object with specified ID.
  interface GetPtrToId
     module procedure GetPtrToIdDamperBase
  end interface

  !> @brief Standard routine for writing an object to file.
  interface WriteObject
     module procedure WriteDamperType
     module procedure WriteDamperBaseType
  end interface

  !> @brief Initializes a damper object.
  interface NullifyDamper
    module procedure NullifyDamperElement
    module procedure NullifyDamperBase
  end interface

  !> @brief Deallocates an array of damper objects.
  interface DeallocateDampers
    module procedure DeallocateDamperElements
    module procedure DeallocateDamperBases
  end interface


  private :: GetPtrToIdDamperBase, WriteDamperType, WriteDamperBaseType
  private :: NullifyDamperElement, NullifyDamperBase
  private :: DeallocateDamperElement, DeallocateDamperBase
  private :: DeallocateDamperElements, DeallocateDamperBases


contains

  !!============================================================================
  !> @brief Returns pointer to (first) base damper with specified id.
  !>
  !> @param[in] array Array of dampertypemodule::damperbasetype objects
  !>                  to search within
  !> @param[in] id Base ID of the object to search for
  !> @param[out] index The array index of the found object
  !>
  !> @details If the damper is not found, NULL is returned.
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 5 Oct 2001

  function GetPtrToIdDamperBase (array,id,index) result(ptr)

    integer             , intent(in)            :: id
    type(DamperBaseType), intent(in) , target   :: array(:)
    integer             , intent(out), optional :: index
    type(DamperBaseType), pointer               :: ptr

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
    write(*,*) '*** GetPtrToIdDamperBase returned nullified, baseId =',id

  end function GetPtrToIdDamperBase


  !!============================================================================
  !> @brief Standard routine for writing an object to io.
  !>
  !> @param[in] dmp The dampertypemodule::damperbasetype object to write
  !> @param[in] io File unit number to write to
  !> @param[in] complexity If present, the value indicates the amount of print
  !>
  !> @author Karl Erik Thoresen
  !>
  !> @date Dec 1998

  subroutine WriteDamperBaseType (dmp,io,complexity)

    use IdTypeModule, only : writeId

    type(DamperBaseType), intent(in) :: dmp
    integer             , intent(in) :: io
    integer,optional    , intent(in) :: complexity

    !! --- Logic section ---

    write(io,'(A)') 'DamperBase','{'
    call writeId (dmp%id,io)

    write(io,*) 'dof         =', dmp%dof

    if (associated(dmp%coeffFunction)) then
       write(io,*) 'coeffFunction =', dmp%coeffFunction%id%baseid
    end if
    if (associated(dmp%forceFunction)) then
       write(io,*) 'forceFunction =', dmp%forceFunction%id%baseid
    end if
    write(io,*) 'dmp0, dmp1  =', dmp%dmp0, dmp%dmp1

    if (associated(dmp%coeffScaleEngine)) then
       write(io,*) 'coeffScaleEngine =', dmp%coeffScaleEngine%id%baseid
    end if
    write(io,*) 'scale1      =', dmp%scale1

    write(io,*) 'saveVar     =', dmp%saveVar

    if (present(complexity)) then
       if (complexity >= 2) then
          write(io,*) 'coeff       =', dmp%coeff
          if (associated(dmp%length))   write(io,*) 'length      =',dmp%length
          if (associated(dmp%velocity)) write(io,*) 'velocity    =',dmp%velocity
          write(io,*) 'force       =', dmp%force
          write(io,*) 'Edmp        =', dmp%Edmp
       end if
    end if

    write(io,'(A)') '}'

  end subroutine WriteDamperBaseType


  !!============================================================================
  !> @brief Standard routine for writing an object to io.
  !>
  !> @param[in] damper The dampertypemodule::dampertype object to write
  !> @param[in] io File unit number to write to
  !> @param[in] complexity If present, the value indicates the amount of print
  !>
  !> @author Karl Erik Thoresen
  !>
  !> @date Dec 1998

  subroutine WriteDamperType (damper,io,complexity)

    type(DamperType), intent(in) :: damper
    integer         , intent(in) :: io
    integer,optional, intent(in) :: complexity

    !! --- Logic section ---

    write(io,'(A)') 'Damper','{'
    if (associated(damper%dmp)) then
       call writeDamperBaseType (damper%dmp,io,complexity)
    end if

    write(io,*) 'samElNum    =', damper%samElNum
    write(io,*) 'nDofs       =', damper%nDofs

    if (associated(damper%triad1)) then
       write(io,*) 'triad1(id)  =',damper%triad1%id%baseId
    end if
    if (associated(damper%triad2)) then
       write(io,*) 'triad2(id)  =',damper%triad2%id%baseId
    end if
    if (present(complexity)) then
       if (associated(damper%forceDir) .and. complexity >= 2) then
          write(io,*) 'forceDir    =', damper%forceDir
       end if
    end if

    write(io,'(A)') '}'

  end subroutine WriteDamperType


  !!============================================================================
  !> @brief Initializes a base damper object.
  !>
  !> @param dmp The dampertypemodule::damperbasetype object to initialize
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 1 Mar 2001

  subroutine NullifyDamperBase (dmp)

    use IdTypeModule, only : nullifyId

    type(DamperBaseType), intent(out) :: dmp

    !! --- Logic section ---

    call nullifyId (dmp%id)

    dmp%dof = 0
    dmp%isActive = .false.
    dmp%saveVar  = .false.

    nullify(dmp%coeffFunction)
    nullify(dmp%forceFunction)
    nullify(dmp%coeffScaleEngine)
    nullify(dmp%spr)

    dmp%coeff      = 0.0_dp
    dmp%dmp0       = 0.0_dp
    dmp%dmp1       = 0.0_dp
    dmp%scale1     = 1.0_dp ! Changed for bushing elements
    nullify(dmp%length)
    dmp%lengthPrev = 0.0_dp
    nullify(dmp%velocity)
    dmp%force      = 0.0_dp
    dmp%forcePrev  = 0.0_dp
    dmp%Edmp       = 0.0_dp

  end subroutine NullifyDamperBase


  !!============================================================================
  !> @brief Initializes a damper element.
  !>
  !> @param[out] damper The dampertypemodule::dampertype object to initialize
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 1 Mar 2001

  subroutine NullifyDamperElement (damper)

    type(DamperType), intent(out) :: damper

    !! --- Logic section ---

    damper%samElNum = 0
    damper%nDOFs = 0
    nullify(damper%dmp)
    nullify(damper%triad1)
    nullify(damper%triad2)
    nullify(damper%forceDir)

  end subroutine NullifyDamperElement


  !!============================================================================
  !> @brief Deallocates a base damper object.
  !>
  !> @param dmp The dampertypemodule::damperbasetype object to deallocate
  !> @param[in] hasTriad1 If .true., this is an axial damper or damper to ground
  !> @param[in] hasTriad2 If .true., this is an axial damper
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Jan 2017

  subroutine DeallocateDamperBase (dmp,hasTriad1,hasTriad2)

    use IdTypeModule, only : deallocateId

    type(DamperBaseType), intent(inout) :: dmp
    logical, optional   , intent(in)    :: hasTriad1, hasTriad2

    !! --- Logic section ---

    call deallocateId (dmp%id)

    if (present(hasTriad1) .and. present(hasTriad2)) then
       if (hasTriad1 .and. associated(dmp%length)) then
          deallocate(dmp%length)
       end if
       if (hasTriad1 .and. hasTriad2 .and. associated(dmp%velocity)) then
          deallocate(dmp%velocity)
       end if
    else if (associated(dmp%spr) .and. associated(dmp%velocity)) then
       deallocate(dmp%velocity)
    end if

    call nullifyDamperBase (dmp)

  end subroutine DeallocateDamperBase


  !!============================================================================
  !> @brief Deallocates a damper element.
  !>
  !> @param damper The dampertypemodule::dampertype object to deallocate
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Jan 2017

  subroutine DeallocateDamperElement (damper)

    type(DamperType), intent(inout) :: damper

    !! --- Logic section ---

    if (associated(damper%dmp)) then

       call DeallocateDamperBase (damper%dmp, &
            &                     associated(damper%triad1), &
            &                     associated(damper%triad2))
       nullify(damper%dmp)

    end if
    if (associated(damper%forceDir)) deallocate(damper%forceDir)

    nullify(damper%triad1)
    nullify(damper%triad2)
    nullify(damper%forceDir)

  end subroutine DeallocateDamperElement


  !!============================================================================
  !> @brief Deallocates an array of base damper objects.
  !>
  !> @param dampers The dampertypemodule::damperbasetype objects to deallocate
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Jan 2017

  subroutine deallocateDamperBases (dampers)

    type(DamperBaseType), pointer :: dampers(:)

    !! Local variables
    integer :: i

    !! --- Logic section ---

    do i = 1, size(dampers)
       call DeallocateDamperBase (dampers(i))
    end do

    deallocate(dampers)
    nullify(dampers)

  end subroutine deallocateDamperBases


  !!============================================================================
  !> @brief Deallocates an array of damper elements
  !>
  !> @param dampers The dampertypemodule::dampertype objects to deallocate
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Jan 2017

  subroutine deallocateDamperElements (dampers)

    type(DamperType), pointer :: dampers(:)

    !! Local variables
    integer :: i

    !! --- Logic section ---

    do i = 1, size(dampers)
       call DeallocateDamperElement (dampers(i))
    end do

    deallocate(dampers)
    nullify(dampers)

  end subroutine deallocateDamperElements


  !!============================================================================
  !> @brief Evaluates the current force for a damper.
  !>
  !> @param[in] dmp The dampertypemodule::damperbasetype object to evaluate for
  !> @param ierr Error flag
  !> @return Current damping force
  !>
  !> @details The current damping force is evaluated as
  !> @verbatim
  !> force = dmp0*velocity + dmp1*int_0^velocity{coeffFunction}dv, or
  !>       = dmp0*velocity + dmp1*forceFunction(velocity)
  !> @endverbatim
  !>
  !> The variable @a dmp0 represents a constant damping coefficient term both
  !> if the damper is defined as a coefficient-velocity function or a
  !> force-velocity function.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 7 Jun 2002

  function damperForce (dmp,ierr)

    use IdTypeModule      , only : getId
    use FunctionTypeModule, only : FunctionIntegral, FunctionValue
    use reportErrorModule , only : reportError, error_p

    type(DamperBaseType), intent(in)    :: dmp
    integer             , intent(inout) :: ierr

    !! Local variables
    real(dp) :: damperForce, fVal
    integer  :: lerr

    !! --- Logic section ---

    if (associated(dmp%coeffFunction)) then

       lerr = ierr
       fVal = FunctionIntegral(dmp%coeffFunction,0.0_dp,dmp%velocity,1,ierr)
       if (ierr < lerr) then
          damperForce = 0.0_dp
          call reportError (error_p,'Damper Id'//getId(dmp%id), &
               &            'Evaluating damper coefficient function failed.')
          return
       end if

    else if (associated(dmp%forceFunction)) then

       lerr = ierr
       fVal = FunctionValue(dmp%forceFunction,dmp%velocity,ierr)
       if (ierr < lerr) then
          damperForce = 0.0_dp
          call reportError (error_p,'Damper Id'//getId(dmp%id), &
               &            'Evaluating gradient of force function failed.')
          return
       end if

    else
       fVal = 0.0_dp
    end if

    damperForce = dmp%dmp0*dmp%velocity + dmp%dmp1*fVal

  end function damperForce


  !!============================================================================
  !> @brief Evaluates the current damping coefficient for a damper.
  !>
  !> @param[in] dmp The dampertypemodule::damperbasetype object to evaluate for
  !> @param ierr Error flag
  !> @return Current damping coefficient
  !>
  !> @details The current damping coefficient, @a coeff, is evaluated as
  !> @verbatim
  !> coeff = dmp0 + dmp1*coeffFunction, or
  !>       = dmp0 + dmp1*d(forceFunction)/dv
  !> @endverbatim
  !>
  !> The variable @a dmp0 represents a constant damping coefficient term both
  !> if the damper is defined as a coefficient-velocity function or a
  !> force-velocity function.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 7 Jun 2002

  function damperCoeff (dmp,ierr)

    use IdTypeModule      , only : getId
    use FunctionTypeModule, only : FunctionValue, FunctionDerivative
    use reportErrorModule , only : reportError, error_p

    type(DamperBaseType), intent(in)    :: dmp
    integer             , intent(inout) :: ierr

    !! Local variables
    real(dp) :: damperCoeff, fVal
    integer  :: lerr

    !! --- Logic section ---

    if (associated(dmp%coeffFunction)) then

       lerr = ierr
       fVal = FunctionValue(dmp%coeffFunction,dmp%velocity,ierr)
       if (ierr < lerr) then
          damperCoeff = 0.0_dp
          call reportError (error_p,'Damper Id'//getId(dmp%id), &
               &            'Evaluating damper coefficient function failed.')
          return
       end if

    else if (associated(dmp%forceFunction)) then

       lerr = ierr
       fVal = FunctionDerivative(dmp%forceFunction,dmp%velocity,1,ierr)
       if (ierr < lerr) then
          damperCoeff = 0.0_dp
          call reportError (error_p,'Damper Id'//getId(dmp%id), &
               &            'Evaluating damper force function failed.')
          return
       end if

    else
       fVal = 0.0_dp
    end if

    damperCoeff = dmp%dmp0 + dmp%dmp1*fVal

  end function damperCoeff

  !!============================================================================
  !> @brief Updates the deformational velocity in a damper.
  !>
  !> @param dmp The dampertypemodule::damperbasetype object to update for
  !> @param[in] dt Time increment size
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 10 Dec 2020

  subroutine updateDamperVelocity (dmp,dt)

    type(DamperBaseType), intent(inout) :: dmp
    real(dp)            , intent(in)    :: dt

    !! --- Logic section ---

    if (associated(dmp%spr) .and. dt > 0.0_dp) then
       !! Compute damper velocity from the associated spring deflection change
       dmp%velocity = (dmp%spr%deflection - dmp%spr%deflectionPrev) / dt
    end if

  end subroutine updateDamperVelocity

end module DamperTypeModule
