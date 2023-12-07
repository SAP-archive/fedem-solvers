!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file springTypeModule.f90
!> @brief Spring object data containers.

!!==============================================================================
!> @brief Module with data types representing spring objects.
!>
!> @details The module also contains subroutines for accessing the spring data.

module SpringTypeModule

  use KindModule        , only : dp
  use IdTypeModule      , only : IdType
  use TriadTypeModule   , only : TriadPtrType
  use FunctionTypeModule, only : FunctionType, EngineType

  implicit none

  real(dp), private, parameter :: eps_p = 1.0e-8_dp !< Zero deflection tolerance


  !> @brief Data type representing a spring failure object.
  type SpringFailureType
     type(IdType) :: id             !< General identification data
     logical      :: deflActive(2)  !< Activation flags for deflection failure
     real(dp)     :: deflRange(2)   !< Deflection range for failure
     logical      :: forceActive(2) !< Activation flags for force failure
     real(dp)     :: forceRange(2)  !< Force range for failure
     logical      :: compFailure    !< Flag for complete joint failure
  end type SpringFailureType


  !> @brief Data type representing a yield force limit.
  type YieldLimitType
     logical  :: active !< Yield limit activation flag
     real(dp) :: ysf(2) !< Yield limit scaling factors
     type(EngineType), pointer :: engine !< Function for varying yield limit
     real(dp) :: yieldLimit !< Current yield limit, = ysf(1) + ysf(2)*engine
  end type YieldLimitType


  !> @brief Data type representing a spring yield object.
  type SpringYieldType
     type(IdType)         :: id  !< General identification data
     type(YieldLimitType) :: pos !< Yield force limit on tension
     type(YieldLimitType) :: neg !< Yield force limit on compression
     logical :: hasDeflectionMax !< Max yield deflection activation flag
     real(dp) :: deflectionMax   !< Max yield deflection before failure
  end type SpringYieldType


  !> @brief Data type representing a single-DOF base spring object.
  type SpringBaseType

     type(IdType) :: id !< General identification data

     !> @brief Turns the spring ON/OFF.
     !> @details Used by cam joints, contact elements and spring failures
     logical :: isActive

     !> @brief Which DOF(s) the spring is associated with
     !> @details
     !>   - = 0 : axial spring or contact element spring
     !>   - > 0 : the joint variable the spring is acting on
     !>   - < 0 : unattached spring (model inconsistency)
     integer :: dof

     !> @brief Flag for cyclic-plastic springs.
     !> @details The value has the following interpretation:
     !>   - = 0 : Unloading along the main force/deflection curve
     !>   - = 1 : Unloading with tangent stiffness, positive deflections
     !>   - = 2 : As 1, but also for negative deflections
     !>   - = 3 : Unloading with secant stiffness, both sides
     integer :: unLoadType

     !> Stress-free length function
     type(EngineType), pointer :: length0Engine
     real(dp) :: l0 !< Constant stress-free length value
     real(dp) :: l1 !< Scaling factor for stress-free length function

     !> @brief Current stress-free length
     !> @details The stress-free length equals l0 + l1*length0Function
     real(dp) :: length0

     !> Stiffness-deflection function
     type(FunctionType), pointer :: stiffnessFunction
     !> Force-deflection function
     type(FunctionType), pointer :: forceFunction
     real(dp) :: s0 !< Constant spring stiffness value
     real(dp) :: s1 !< Scaling factor for spring stiffness/force function

     !> @brief Current spring stiffness.
     !> @details The spring stiffness equals either s0 + s1*stiffnessFunction,
     !> or s0 + s1*d(forceFunction)/dx where @a s0 is a stiffness coefficient
     !> both for @a stiffnessFunction and @a forceFunction, i.e.,
     !> if @a forceFunction is constant the @a s0 value represents a stiffness
     real(dp) :: stiffness

     !> Optional scaling of stiffness for positive deflection
     type(EngineType), pointer :: stiffScaleEnginePos
     !> Optional scaling of stiffness for negative deflection
     type(EngineType), pointer :: stiffScaleEngineNeg
     real(dp) :: scale1 !< Current scaling coefficient = scale1*stiffScaleEngine

     type(SpringFailureType), pointer :: failure !< Data for failure calculation
     type(SpringYieldType)  , pointer :: yield   !< Data for yield calculation
     real(dp) :: yieldDeflection     !< Yield deflection at current time step
     real(dp) :: yieldDeflectionPrev !< Yield deflection at previous time step
     logical  :: hasFailed !< Set to .true. when the spring has failed
     logical  :: isInYield !< Set to .true. when the spring is yielding

     !> @brief Spring length.
     !> @details For joint springs, this pointer points to the
     !> corresponding joint variable. For axial springs it is allocated.
     real(dp), pointer :: length
     logical :: allocatedLength !< Set to .true. if the @a length is allocated

     real(dp) :: lengthPrev     !< Length at previous time step
     real(dp) :: deflection     !< Current spring deflection
     real(dp) :: deflectionPrev !< Previous spring deflection
     real(dp) :: force          !< Current spring force
     real(dp) :: forcePrev      !< Previous spring force
     real(dp) :: Estr           !< Current strain energy
     real(dp) :: Einp           !< Accumulated input energy
     real(dp) :: Eyield         !< Accumulated yield energy loss

     logical :: saveVar(5) !< Flags indicating which variables should be saved

  end type SpringBaseType


  !> @brief Data type representing a base spring pointer.
  !> @details This data type is used to construct arrays of springs where each
  !> element is a pointer to a spring object, and not the objects themselves.
  type SpringPtrType
     type(SpringBaseType), pointer :: p !< Pointer to a base spring object
  end type SpringPtrType


  !> @brief Data type representing multi-DOF spring element.
  type SpringType

     type(IdType), pointer :: id !< General identification data

     !> Pointers to axial, global or joint DOF springs
     type(SpringPtrType), pointer :: spr(:)

     logical, pointer :: lConn(:,:) !< Spring inter-connectivity table

     integer :: samElNum !< Element number for SAM reference (MPMNPC)

     !> @brief Number of degrees of freedom in the spring element.
     !> @details Equals @a nJointDOFs for joint springs, 2*3 or 2*6
     !> for axial- and global DOF springs, and 3 or 6 for springs to ground
     integer :: nDOFs

     real(dp) :: alpha2 !< Stiffness-proportional damping factor
     logical, pointer :: hasFailed !< .true. if the spring is failed entirely

     type(TriadPtrType), pointer :: triads(:) !< Nullified for joint springs

     real(dp), pointer :: forceDir(:,:) !< For axial- and global DOF springs
     real(dp), pointer :: stiffMat(:,:) !< For inter-connected springs
     real(dp), pointer :: couplStiff(:) !< For explicit coupling stiffness

  end type SpringType


  !> @brief Returns pointer to object with specified ID.
  interface GetPtrToId
     module procedure GetPtrToIdSpringElement
     module procedure GetPtrToIdSpringBase
     module procedure GetPtrToIdSpringFailure
     module procedure GetPtrToIdSpringYield
  end interface

  !> @brief Standard routine for writing an object to file.
  interface WriteObject
    module procedure WriteSpringElementType
    module procedure WriteSpringBaseType
  end interface

  !> @brief Initializes a spring object.
  interface NullifySpring
    module procedure NullifySpringElement
    module procedure NullifySpringBase
  end interface

  !> @brief Deallocates an array of spring objects.
  interface DeallocateSprings
    module procedure DeallocateSpringElements
    module procedure DeallocateSpringBases
    module procedure DeallocateFailures
    module procedure DeallocateYields
  end interface

  !> @brief Updates the state variables pertaining to previous time step.
  interface updateAtConvergence
     module procedure updatePreviousValues
  end interface

  !> @brief Restores the state variables from the last converged time step.
  interface restoreFromLastStep
     module procedure restorePreviousValues
  end interface

  !> Checks if a spring has passed failure state.
  interface checkFailure
     module procedure checkSpringFailure
     module procedure checkDofFailure
  end interface


  private :: GetPtrToIdSpringElement, GetPtrToIdSpringBase
  private :: GetPtrToIdSpringFailure,  GetPtrToIdSpringYield
  private :: WriteSpringElementType, WriteSpringBaseType
  private :: NullifySpringElement, NullifySpringBase, DeallocateSpringBase
  private :: DeallocateSpringElements, DeallocateSpringBases
  private :: updatePreviousValues, restorePreviousValues
  private :: checkSpringFailure, checkDofFailure


contains

  !!============================================================================
  !> @brief Returns pointer to (first) spring element with specified id.
  !>
  !> @param[in] array Array of springtypemodule::springtype objects
  !>                  to search within
  !> @param[in] id Base ID of the object to search for
  !>
  !> @details If the spring is not found, NULL is returned.
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 1 Nov 1999

  function GetPtrToIdSpringElement (array,id) result(ptr)

    integer         , intent(in)         :: id
    type(SpringType), intent(in), target :: array(:)
    type(SpringType), pointer            :: ptr

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
    write(*,*) '*** GetPtrToIdSpringElement returned nullified, BaseId =',id

  end function GetPtrToIdSpringElement


  !!============================================================================
  !> @brief Returns pointer to (first) base spring with specified id.
  !>
  !> @param[in] array Array of springtypemodule::springbasetype objects
  !>                  to search within
  !> @param[in] id Base ID of the object to search for
  !> @param[out] index The array index of the found object
  !>
  !> @details If the spring is not found, NULL is returned.
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 6 Jun 2002

  function GetPtrToIdSpringBase (array,id,index) result(ptr)

    integer             , intent(in)            :: id
    type(SpringBaseType), intent(in) , target   :: array(:)
    integer             , intent(out), optional :: index
    type(SpringBaseType), pointer               :: ptr

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
    write(*,*) '*** GetPtrToIdSpringBase returned nullified, BaseId =',id

  end function GetPtrToIdSpringBase


  !!============================================================================
  !> @brief Returns pointer to (first) spring failure with specified id.
  !>
  !> @param[in] array Array of springtypemodule::springfailuretype objects
  !>                  to search within
  !> @param[in] id Base ID of the object to search for
  !> @param[out] index The array index of the found object
  !>
  !> @details If the spring failure is not found, NULL is returned.
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 24 Apr 2005

  function GetPtrToIdSpringFailure (array,id,index) result(ptr)

    integer                , intent(in)            :: id
    type(SpringFailureType), intent(in) , target   :: array(:)
    integer                , intent(out), optional :: index
    type(SpringFailureType), pointer               :: ptr

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
    write(*,*) '*** GetPtrToIdSpringFailure returned nullified, BaseId =',id

  end function GetPtrToIdSpringFailure


  !!============================================================================
  !> @brief Returns pointer to (first) spring yield with specified id.
  !>
  !> @param[in] array Array of springtypemodule::springyieldtype objects
  !>                  to search within
  !> @param[in] id Base ID of the object to search for
  !> @param[out] index The array index of the found object
  !>
  !> @details If the spring yield is not found, NULL is returned.
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 20 Jun 2005

  function GetPtrToIdSpringYield (array,id,index) result(ptr)

    integer              , intent(in)            :: id
    type(SpringYieldType), intent(in) , target   :: array(:)
    integer              , intent(out), optional :: index
    type(SpringYieldType), pointer               :: ptr

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
    write(*,*) '*** GetPtrToIdSpringYield returned nullified, BaseId =',id

  end function GetPtrToIdSpringYield


  !!============================================================================
  !> @brief Standard routine for writing an object to io.
  !>
  !> @param[in] spr The springtypemodule::springbasetype object to write
  !> @param[in] io File unit number to write to
  !> @param[in] complexity If present, the value indicates the amount of print
  !>
  !> @author Karl Erik Thoresen
  !>
  !> @date Dec 1998

  subroutine WriteSpringBaseType (spr,io,complexity)

    use IdTypeModule, only : writeId

    type(SpringBaseType), intent(in) :: spr
    integer             , intent(in) :: io
    integer, optional   , intent(in) :: complexity

    !! --- Logic section ---

    write(io,'(A)') 'SpringBase','{'
    call writeId (spr%id,io)

    write(io,*) 'dof         =', spr%dof
    write(io,*) 'unLoadType  =', spr%unLoadType

    if (associated(spr%length0Engine)) then
       write(io,*) 'length0Engine =', spr%length0Engine%id%baseId
    end if
    write(io,*) 'l0, l1      =', spr%l0, spr%l1

    if (associated(spr%stiffnessFunction)) then
       write(io,*) 'stiffnessFunction =', spr%stiffnessFunction%id%baseId
    end if
    if (associated(spr%forceFunction)) then
       write(io,*) 'forceFunction =', spr%forceFunction%id%baseId
    end if
    write(io,*) 's0, s1      =', spr%s0, spr%s1

    if (associated(spr%stiffScaleEnginePos)) then
       write(io,*) 'stiffScaleEnginePos =', spr%stiffScaleEnginePos%id%baseId
    end if
    if (associated(spr%stiffScaleEngineNeg)) then
       write(io,*) 'stiffScaleEngineNeg =', spr%stiffScaleEngineNeg%id%baseId
    end if
    write(io,*) 'scale1      =', spr%scale1

    if (associated(spr%failure)) then
       write(io,*) 'compFailure =', spr%failure%compFailure
       if (spr%failure%deflActive(2)) then
          write(io,*) 'deflMax     =', spr%failure%deflRange(2)
       end if
       if (spr%failure%deflActive(1)) then
          write(io,*) 'deflMin     =', spr%failure%deflRange(1)
       end if
       if (spr%failure%forceActive(2)) then
          write(io,*) 'forceMax    =', spr%failure%forceRange(2)
       end if
       if (spr%failure%forceActive(1)) then
          write(io,*) 'forceMin    =', spr%failure%forceRange(1)
       end if
       write(io,*) 'hasFailed   =', spr%hasFailed
    end if

    if (associated(spr%yield)) then
       if (spr%yield%pos%active) then
          write(io,*) 'yieldMax    =', spr%yield%pos%yieldLimit
          write(io,*) 'yMax0,yMax1 =', spr%yield%pos%ysf
          if (associated(spr%yield%pos%engine)) then
             write(io,*) 'yieldMaxEngine =', spr%yield%pos%engine%id%baseId
          end if
       end if
       if (spr%yield%neg%active) then
          write(io,*) 'yieldMin    =', spr%yield%neg%yieldLimit
          write(io,*) 'yMin0,yMin1 =', spr%yield%neg%ysf
          if (associated(spr%yield%neg%engine)) then
             write(io,*) 'yieldMinEngine =', spr%yield%neg%engine%id%baseId
          end if
       end if
       if (spr%yield%hasDeflectionMax) then
          write(io,*) 'yieldDefMax =', spr%yield%deflectionMax
       end if
       write(io,*) 'isInYield   =', spr%isInYield
    end if

    write(io,*) 'saveVar     =', spr%saveVar

    if (present(complexity)) then
       if (complexity >= 2) then
          write(io,*) 'force       =', spr%force
          write(io,*) 'stiffness   =', spr%stiffness
          if (associated(spr%length)) write(io,*) 'length      =', spr%length
          write(io,*) 'length0     =', spr%length0
          write(io,*) 'deflection  =', spr%deflection
          write(io,*) 'yDeflection =', spr%yieldDeflection
          write(io,*) 'Estr        =', spr%Estr
          write(io,*) 'Einp        =', spr%Einp
          write(io,*) 'Eyield      =', spr%Eyield
       end if
    end if

    write(io,'(A)') '}'

  end subroutine WriteSpringBaseType


  !!============================================================================
  !> @brief Standard routine for writing an object to io.
  !>
  !> @param[in] spring The springtypemodule::springtype object to write
  !> @param[in] io File unit number to write to
  !> @param[in] complexity If present, the value indicates the amount of print
  !>
  !> @author Karl Erik Thoresen
  !>
  !> @date Dec 1998

  subroutine WriteSpringElementType (spring,io,complexity)

    use IdTypeModule, only : writeId

    type(SpringType) , intent(in) :: spring
    integer          , intent(in) :: io
    integer, optional, intent(in) :: complexity

    !! Local variables
    integer :: i

    !! --- Logic section ---

    write(io,'(A)') 'Spring','{'
    call writeId (spring%id,io)
    do i = 1, size(spring%spr)
       if (associated(spring%spr(i)%p)) then
          call writeSpringBaseType (spring%spr(i)%p,io,complexity)
       end if
    end do

    write(io,*) 'samElNum    =', spring%samElnum
    write(io,*) 'nDOFs       =', spring%nDOFs
    write(io,*) 'alpha2      =', spring%alpha2
    if (associated(spring%hasFailed)) then
       write(io,*) 'hasFailed   =', spring%hasFailed
    end if
    if (associated(spring%triads)) then
       write(io,*) 'triads(id)  =',(spring%triads(i)%p%id%baseId, &
            &                       i=1,size(spring%triads))
    end if
    if (present(complexity)) then
       if (associated(spring%forceDir) .and. complexity >= 2) then
          do i = 1, size(spring%forceDir,2)
             write(io,*) 'forceDir    =', spring%forceDir(:,i)
          end do
       end if
    end if

    write(io,'(A)') '}'

  end subroutine WriteSpringElementType


  !!============================================================================
  !> @brief Initializes a spring failure object.
  !>
  !> @param failure The springtypemodule::springfailuretype object to initialize
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 14 Mar 2006

  subroutine NullifyFailure (failure)

    use IdTypeModule, only : nullifyId

    type(SpringFailureType), intent(out) :: failure

    !! --- Logic section ---

    call nullifyId (failure%id)

    failure%deflActive = .false.
    failure%deflRange = 0.0_dp

    failure%forceActive = .false.
    failure%forceRange = 0.0_dp

    failure%compFailure = .false.

  end subroutine NullifyFailure


  !!============================================================================
  !> @brief Deallocates an array of spring failure objects.
  !>
  !> @param fails The springtypemodule::springfailuretype objects to deallocate
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Jan 2017

  subroutine DeallocateFailures (fails)

    use IdTypeModule, only : deallocateId

    type(SpringFailureType), pointer :: fails(:)

    !! Local variables
    integer :: i

    !! --- Logic section ---

    do i = 1, size(fails)
       call deallocateId (fails(i)%id)
    end do
    deallocate(fails)
    nullify(fails)

  end subroutine DeallocateFailures


  !!============================================================================
  !> @brief Initializes a spring yield object.
  !>
  !> @param yield The springtypemodule::springyieldtype object to initialize
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 14 Mar 2006

  subroutine NullifyYield (yield)

    use IdTypeModule, only : nullifyId

    type(SpringYieldType), intent(out) :: yield

    !! --- Logic section ---

    call nullifyId (yield%id)

    yield%pos%active = .false.
    nullify(yield%pos%engine)
    yield%pos%yieldLimit = 0.0_dp
    yield%pos%ysf = 0.0_dp

    yield%neg%active = .false.
    nullify(yield%neg%engine)
    yield%neg%yieldLimit = 0.0_dp
    yield%neg%ysf = 0.0_dp

    yield%hasDeflectionMax = .false.
    yield%deflectionMax = 0.0_dp

  end subroutine NullifyYield


  !!============================================================================
  !> @brief Deallocates an array of spring yield objects.
  !>
  !> @param yields The springtypemodule::springyieldtype objects to deallocate
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Jan 2017

  subroutine DeallocateYields (yields)

    use IdTypeModule, only : deallocateId

    type(SpringYieldType), pointer :: yields(:)

    !! Local variables
    integer :: i

    !! --- Logic section ---

    do i = 1, size(yields)
       call deallocateId (yields(i)%id)
    end do
    deallocate(yields)
    nullify(yields)

  end subroutine DeallocateYields


  !!============================================================================
  !> @brief Initializes a spring base object.
  !>
  !> @param spr The springtypemodule::springbasetype object to initialize
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 1 Mar 2001

  subroutine NullifySpringBase (spr)

    use IdTypeModule, only : nullifyId

    type(SpringBaseType), intent(out) :: spr

    !! --- Logic section ---

    call nullifyId (spr%id)

    spr%dof = -1
    spr%unLoadType = 0
    spr%isActive = .true.
    spr%hasFailed = .false.
    spr%isInYield = .false.
    spr%saveVar = .false.
    spr%allocatedLength = .false.

    nullify(spr%length0Engine)
    nullify(spr%stiffnessFunction)
    nullify(spr%forceFunction)
    nullify(spr%stiffScaleEnginePos)
    nullify(spr%stiffScaleEngineNeg)
    nullify(spr%failure)
    nullify(spr%yield)

    spr%stiffness      = 0.0_dp
    spr%length0        = 0.0_dp
    spr%l0             = 0.0_dp
    spr%l1             = 0.0_dp
    spr%s0             = 0.0_dp
    spr%s1             = 0.0_dp
    spr%scale1         = 1.0_dp ! Changed for bushing elements
    nullify(spr%length)
    spr%lengthPrev     = 0.0_dp
    spr%deflection     = 0.0_dp
    spr%deflectionPrev = 0.0_dp
    spr%yieldDeflection = 0.0_dp
    spr%yieldDeflectionPrev = 0.0_dp
    spr%force          = 0.0_dp
    spr%forcePrev      = 0.0_dp
    spr%Estr           = 0.0_dp
    spr%Einp           = 0.0_dp
    spr%Eyield         = 0.0_dp

  end subroutine NullifySpringBase


  !!============================================================================
  !> @brief Initializes a spring element.
  !>
  !> @param spring The springtypemodule::springtype object to initialize
  !> @param[in] sprId Id of the owning joint (if joint spring)
  !> @param[in] nDofs Number of degrees of freedom
  !>
  !> @author  Knut Morten Okstad
  !>
  !> @date 1 Mar 2001

  subroutine NullifySpringElement (spring,sprId,nDofs)

    use IdTypeModule, only : IdType

    type(SpringType)              , intent(out) :: spring
    type(IdType), optional, target, intent(in)  :: sprId
    integer     , optional        , intent(in)  :: nDofs

    !! --- Logic section ---

    if (present(sprId)) then
       spring%id => sprId
    else
       nullify(spring%id)
    end if

    nullify(spring%spr)
    nullify(spring%lConn)

    spring%samElNum = 0
    if (present(nDofs)) then
       spring%nDOFs = nDOFs
    else
       spring%nDOFs = 0
    end if
    spring%alpha2 = 0.0_dp

    nullify(spring%hasFailed)
    nullify(spring%triads)
    nullify(spring%forceDir)
    nullify(spring%stiffMat)
    nullify(spring%couplStiff)

  end subroutine NullifySpringElement


  !!============================================================================
  !> @brief Deallocates a base spring object.
  !>
  !> @param spr The springtypemodule::springbasetype object to deallocate
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Jan 2017

  subroutine DeallocateSpringBase (spr)

    use IdTypeModule, only : deallocateId

    type(SpringBaseType), intent(inout) :: spr

    !! --- Logic section ---

    call deallocateId (spr%id)

    if (spr%allocatedLength) deallocate(spr%length)

    call nullifySpringBase (spr)

  end subroutine DeallocateSpringBase


  !!============================================================================
  !> @brief Deallocates an array of base spring object.
  !>
  !> @param springs The springtypemodule::springbasetype objects to deallocate
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Jan 2017

  subroutine deallocateSpringBases (springs)

    type(SpringBaseType), pointer :: springs(:)

    !! Local variables
    integer :: i

    !! --- Logic section ---

    do i = 1, size(springs)
       call DeallocateSpringBase (springs(i))
    end do

    deallocate(springs)
    nullify(springs)

  end subroutine deallocateSpringBases


  !!============================================================================
  !> @brief Deallocates a spring element object.
  !>
  !> @param spring The springtypemodule::springtype object to deallocate
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Jan 2017

  subroutine DeallocateSpring (spring)

    use IdTypeModule, only : deallocateId

    type(SpringType), intent(inout) :: spring

    !! --- Logic section ---

    if (associated(spring%id)) then
       call deallocateId (spring%id)
       deallocate(spring%id)
    end if
    if (associated(spring%spr))        deallocate(spring%spr)
    if (associated(spring%lConn))      deallocate(spring%lConn)
    if (associated(spring%hasFailed))  deallocate(spring%hasFailed)
    if (associated(spring%triads))     deallocate(spring%triads)
    if (associated(spring%forceDir))   deallocate(spring%forceDir)
    if (associated(spring%stiffMat))   deallocate(spring%stiffMat)
    if (associated(spring%couplStiff)) deallocate(spring%couplStiff)

    call nullifySpringElement (spring)

  end subroutine DeallocateSpring


  !!============================================================================
  !> @brief Deallocates an array of spring element objects.
  !>
  !> @param springs The springtypemodule::springtype objects to deallocate
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Jan 2017

  subroutine deallocateSpringElements (springs)

    type(SpringType), pointer :: springs(:)

    !! Local variables
    integer :: i

    !! --- Logic section ---

    do i = 1, size(springs)
       call DeallocateSpring (springs(i))
    end do

    deallocate(springs)
    nullify(springs)

  end subroutine deallocateSpringElements


  !!============================================================================
  !> @brief Checks whether a spring element is active or not.
  !>
  !> @param[in] spring The springtypemodule::springtype object to check
  !> @return .true. if the spring is still active, otherwise .false.
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 26 Jan 2011

  function isActive (spring)

    type(SpringType), intent(in) :: spring

    !! Local variables
    logical :: isActive
    integer :: i

    !! --- Logic section ---

    isActive = .true.

    do i = 1, size(spring%spr)
       if (associated(spring%spr(i)%p)) then
          if (spring%spr(i)%p%isActive) return
       end if
    end do

    isActive = .false.

  end function isActive


  !!============================================================================
  !> @brief Evaluates the current force in a spring.
  !>
  !> @param[in] spr The base spring to evaluate the force in
  !> @param[in] deflection Current spring deflection
  !> @param[out] ierr Error flag
  !> @return Current spring force
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 7 Jun 2002

  function springForce (spr,deflection,ierr)

    use IdTypeModule      , only : getId
    use FunctionTypeModule, only : FunctionValue, FunctionIntegral
    use reportErrorModule , only : reportError, error_p

    type(SpringBaseType), intent(in)    :: spr
    real(dp)            , intent(in)    :: deflection
    integer             , intent(inout) :: ierr

    !! Local variables
    real(dp) :: springForce, fVal
    integer  :: lerr

    !! --- Logic section ---

    if (spr%hasFailed) then
       springForce = 0.0_dp
       return
    end if

    fVal = 0.0_dp
    if (associated(spr%stiffnessFunction)) then

       lerr = ierr
       fVal = FunctionIntegral(spr%stiffnessFunction,0.0_dp,deflection,1,ierr)
       if (ierr < lerr) then
          call reportError (error_p,'Spring'//getId(spr%id), &
               &            'Integrating elastic stiffness function failed.')
       end if

    else if (associated(spr%forceFunction)) then
       if (spr%unLoadType /= 1 .or. deflection >= -eps_p) then

          lerr = ierr
          if (spr%unLoadType == 1) then ! use the total deflection here
             fVal = FunctionValue(spr%forceFunction,spr%deflection,ierr)
          else
             fVal = FunctionValue(spr%forceFunction,deflection,ierr)
          end if
          if (ierr < lerr) then
             call reportError (error_p,'Spring'//getId(spr%id), &
                  &            'Evaluating elastic force function failed.')
          end if

       end if
    end if

    if (spr%unLoadType > 0 .and. associated(spr%forceFunction)) then
       !! This is a cyclic-plastic spring
       if (spr%unLoadType == 1 .and. deflection < -eps_p) then
          !! Elastic compression ==> no force
          springForce = 0.0_dp
       else if (spr%s0*abs(deflection) <= spr%s1*abs(fVal)) then
          !! We are in the linear unloading range ==> use tangent stiffness, s0
          springForce = spr%s0*deflection
       else
          !! Beyond the linear range ==> use the force/deflection function value
          springForce = spr%s1*fVal
       end if
    else
       springForce = spr%s0*deflection + spr%s1*fVal
    end if

  end function springForce


  !!============================================================================
  !> @brief Evaluates the current stiffness in a spring.
  !>
  !> @param[in] spr The base spring to evaluate the stiffness in
  !> @param[in] deflection Current spring deflection
  !> @param[out] ierr Error flag
  !> @return Current spring stiffness
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 7 Jun 2002

  function springStiff (spr,deflection,ierr)

    use IdTypeModule      , only : getId
    use FunctionTypeModule, only : FunctionValue, FunctionDerivative
    use reportErrorModule , only : reportError, error_p

    type(SpringBaseType), intent(in)    :: spr
    real(dp)            , intent(in)    :: deflection
    integer             , intent(inout) :: ierr

    !! Local variables
    real(dp) :: springStiff, fVal, fLin
    integer  :: lerr

    !! --- Logic section ---

    if (spr%hasFailed) then
       springStiff = 0.0_dp
       return
    end if

    fVal = 0.0_dp
    fLin = 0.0_dp
    if (associated(spr%stiffnessFunction)) then

       lerr = ierr
       fVal = FunctionValue(spr%stiffnessFunction,deflection,ierr)
       if (ierr < lerr) then
          call reportError (error_p,'Spring'//getId(spr%id), &
               &            'Evaluating elastic stiffness function failed.')
       end if

    else if (associated(spr%forceFunction)) then
       fLin = spr%s0*deflection - spr%force
       if ( (spr%unLoadType < 1 .or. abs(fLin) > spr%s0*eps_p) .and. &
            (spr%unLoadType /= 1 .or. deflection >= -eps_p) ) then

          lerr = ierr
          if (spr%unLoadType == 1) then ! use the total deflection here
             fVal = FunctionDerivative(spr%forceFunction,spr%deflection,1,ierr)
          else
             fVal = FunctionDerivative(spr%forceFunction,deflection,1,ierr)
          end if
          if (ierr < lerr) then
             call reportError (error_p,'Spring'//getId(spr%id), &
                  &            'Evaluating elastic tangent stiffness failed.')
          end if

       end if
    end if

    if (spr%unLoadType > 0 .and. associated(spr%forceFunction)) then
       !! This is a cyclic-plastic spring
       if (spr%unLoadType == 1 .and. deflection < -eps_p) then
          !! Elastic compression ==> no stiffness
          springStiff = 0.0_dp
       else if (abs(fLin) <= spr%s0*eps_p) then
          !! We are in the linear unloading range ==> use tangent stiffness, s0
          springStiff = spr%s0
       else
          !! Beyond the linear range ==> use the force/deflection function value
          springStiff = spr%s1*fVal
       end if
    else
       springStiff = spr%s0 + spr%s1*fVal
    end if

  end function springStiff


  !!============================================================================
  !> @brief Evaluates the current strain energy in a spring.
  !>
  !> @param[in] spr The base spring to evaluate the energy in
  !> @param[in] deflection Current spring deflection
  !> @param[out] ierr Error flag
  !> @return Current strain energy
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 20 Jan 2011

  function springEnergy (spr,deflection,ierr)

    use FunctionTypeModule     , only : CanIntegrate, FunctionIntegral
    use ExplicitFunctionsModule, only : funcType_p
    use reportErrorModule      , only : reportError, error_p

    type(SpringBaseType), intent(in)    :: spr
    real(dp)            , intent(in)    :: deflection
    integer             , intent(inout) :: ierr

    !! Local variables
    real(dp) :: springEnergy, fInt

    !! --- Logic section ---

    if (spr%hasFailed) then
       springEnergy = 0.0_dp
       return
    end if

    fInt = 0.0_dp
    if (associated(spr%stiffnessFunction)) then

       if (CanIntegrate(spr%stiffnessFunction,2)) then
          fInt = FunctionIntegral(spr%stiffnessFunction, &
               &                  0.0_dp,deflection,2,ierr)
       else
          ierr = ierr - 1
          call reportError (error_p,'Strain energy integration is not imple'// &
               &            'mented for stiffness functions of type '// &
               &            funcType_p(spr%stiffnessFunction%type))
       end if

    else if (associated(spr%forceFunction)) then

       !!TODO,kmo: Modify for cyclic-plastic springs (unLoadType > 0)
       if (CanIntegrate(spr%forceFunction,1)) then
          fInt = FunctionIntegral(spr%forceFunction,0.0_dp,deflection,1,ierr)
       else
          ierr = ierr - 1
          call reportError (error_p,'Strain energy integration is not imple'// &
               &            'mented for force-deflection functions of type '// &
               &            funcType_p(spr%forceFunction%type))
       end if

    end if

    !! Note that the constant value s0 is interpreted as a stiffness
    !! regardless of whether stiffnessFunction or forceFunction is specified
    springEnergy = 0.5_dp*spr%s0*deflection*deflection + spr%s1*fInt

  end function springEnergy


  !!============================================================================
  !> @brief Updates the state variables pertaining to the previous time step.
  !>
  !> @param spr The springtypemodule::springbasetype object to update state for
  !>
  !> @details This subroutine is invoked once after convergence as been reached.
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 4 Jul 2022

  subroutine updatePreviousValues (spr)

    type(SpringBaseType), intent(inout) :: spr

    !! --- Logic section ---

    spr%lengthPrev          = spr%length
    spr%deflectionPrev      = spr%deflection
    spr%yieldDeflectionPrev = spr%yieldDeflection
    spr%forcePrev           = spr%force

  end subroutine updatePreviousValues


  !!============================================================================
  !> @brief Restores the state variables from the last converged time step.
  !>
  !> @param spr The springtypemodule::springbasetype object to restore state for
  !>
  !> @details This subroutine is invoked when doing iteration cut-back.
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 4 Jul 2022

  subroutine restorePreviousValues (spr)

    type(SpringBaseType), intent(inout) :: spr

    !! --- Logic section ---

    spr%length          = spr%lengthPrev
    spr%deflection      = spr%deflectionPrev
    spr%yieldDeflection = spr%yieldDeflectionPrev
    spr%force           = spr%forcePrev

  end subroutine restorePreviousValues


  !!============================================================================
  !> @brief Checks if a spring has passed failure state.
  !>
  !> @param spr The springtypemodule::springbasetype object to check for
  !>
  !> @author Bjorn Haugen
  !> @date 25 Apr 2005
  !>
  !> @author Knut Morten Okstad
  !> @date 14 Mar 2006

  subroutine checkDofFailure (spr)

    use IdTypeModule     , only : getId
    use reportErrorModule, only : getErrorFile

    type(SpringBaseType), intent(inout), target :: spr

    !! Local variables
    integer                          :: lpu
    type(SpringFailureType), pointer :: f
    type(SpringYieldType)  , pointer :: y

    !! --- Logic section ---

    if (spr%hasFailed) return

    if (associated(spr%failure)) then
       f => spr%failure
       if (f%deflActive(2) .and. spr%deflection > f%deflRange(2)) then
          spr%hasFailed = .true.
       else if (f%deflActive(1) .and. spr%deflection < f%deflRange(1)) then
          spr%hasFailed = .true.
       else if (f%forceActive(2) .and. spr%force > f%forceRange(2)) then
          spr%hasFailed = .true.
       else if (f%forceActive(1) .and. spr%force < f%forceRange(1)) then
          spr%hasFailed = .true.
       end if
    end if

    if (associated(spr%yield)) then
       y => spr%yield
       if (y%hasDeflectionMax) then
          !! The yield failure criterion is identical on tension and compression
          if (spr%yieldDeflection > y%deflectionMax) then
             spr%hasFailed = .true.
          else if (spr%yieldDeflection < -y%deflectionMax) then
             spr%hasFailed = .true.
          end if
       end if
    end if

    if (spr%hasFailed) then
       lpu = getErrorFile()
       write(lpu,600,advance='NO') trim(getId(spr%id))
       if (associated(spr%yield)) then
          write(lpu,601) spr%force,spr%deflection,spr%yieldDeflection
       else
          write(lpu,602) spr%force,spr%deflection
       end if
    end if

600 format(5X,'Spring',A,' has failed')
601 format(', force =',1PE13.5, ', deflection =',E13.5, &
         & ', yield deflection =',E13.5)
602 format(', force =',1PE13.5, ', deflection =',E13.5)

  end subroutine checkDofFailure


  !!============================================================================
  !> @brief Makes the entire spring element fail if one of its DOFs has failed.
  !>
  !> @param spring The springtypemodule::springtype object to check for
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 7 Apr 2010

  subroutine checkSpringFailure (spring)

    use IdTypeModule     , only : getId
    use reportErrorModule, only : getErrorFile

    type(SpringType), intent(out) :: spring

    !! Local variables
    integer :: i, j

    !! --- Logic section ---

    if (.not. associated(spring%hasFailed)) return
    if (spring%hasFailed) return

    do i = 1, size(spring%spr)
       if (associated(spring%spr(i)%p)) then
          if (spring%spr(i)%p%hasFailed) then
             do j = 1, size(spring%spr)
                if (associated(spring%spr(j)%p)) then
                   spring%spr(j)%p%hasFailed = .true.
                end if
             end do
             spring%hasFailed = .true.
             exit
          end if
       end if
    end do

    if (spring%hasFailed) then
       j = getErrorFile()
       write(j,600) trim(getId(spring%id))
    end if

600 format(5X,'Spring element',A,' has failed entirely.')

  end subroutine checkSpringFailure


  !!============================================================================
  !> @brief Enforces the yield limits on current stiffness and force values.
  !>
  !> @param spr The springtypemodule::springbasetype object to check yield for
  !> @param[in] yPos Yield limit on tension
  !> @param[in] yNeg Yield limit on compression
  !>
  !> @author Bjorn Haugen
  !> @date 25 Apr 2005
  !>
  !> @author Knut Morten Okstad
  !> @date 14 Mar 2006

  subroutine checkYieldLimits (spr,yPos,yNeg)

    use IdTypeModule       , only : getId
    use reportErrorModule  , only : getErrorFile
#ifdef FT_DEBUG
    use fileUtilitiesModule, only : getDBGfile
#endif

    type(SpringBaseType), intent(inout) :: spr
    type(YieldLimitType), intent(in)    :: yPos, yNeg

    !! Local variables
    integer  :: lpu
    logical  :: wasInYield
    real(dp) :: elasticForce, elasticStiff, tol

    !! --- Logic section ---

    if (spr%hasFailed) return

    wasInYield   = spr%isInYield
    elasticStiff = spr%stiffness
    elasticForce = spr%force

    if (yPos%active .and. elasticForce > yPos%yieldLimit) then

       !! Yield in tension
       spr%isInYield = .true.
       spr%force     = yPos%yieldLimit

    else if (yNeg%active .and. elasticForce < yNeg%yieldLimit) then

       !! Yield in compression
       spr%isInYield = .true.
       spr%force     = yNeg%yieldLimit

    else if (wasInYield) then

       !! Spring was already yielding, check for re-entry to elastic state
       if (yPos%active .and. yNeg%active) then
          tol = eps_p*max(1.0_dp,yPos%yieldLimit-yNeg%yieldLimit)
          if ( elasticForce < yPos%yieldLimit-tol .and. &
               elasticForce > yNeg%yieldLimit+tol ) spr%isInYield = .false.
       else if (yPos%active) then
          tol = eps_p*max(1.0_dp,abs(yPos%yieldLimit))
          if (elasticForce < yPos%yieldLimit-tol) spr%isInYield = .false.
       else if (yNeg%active) then
          tol = eps_p*max(1.0_dp,abs(yNeg%yieldLimit))
          if (elasticForce > yNeg%yieldLimit+tol) spr%isInYield = .false.
       end if

    end if

    if (spr%isInYield) then

       !! Spring is yielding, set zero stiffness and accumulate yield deflection
       spr%stiffness = 0.0_dp
       if (elasticStiff < -eps_p .or. elasticStiff > eps_p) then
          spr%yieldDeflection = spr%yieldDeflection + &
               &               (elasticForce-spr%force)/elasticStiff
       end if

    end if

#ifdef FT_DEBUG
    if (spr%isInYield .or. wasInYield) then
       lpu = getDBGFile(55,'yield.dbg')
       if (spr%isInYield) then
          write(lpu,600) trim(getId(spr%id)),'is yielding',elasticForce
       else
          write(lpu,600) trim(getId(spr%id)),'enters elastic state',elasticForce
       end if
       write(lpu,*) 'yieldDeflectPrev =', spr%yieldDeflectionPrev
       write(lpu,*) 'yieldDeflection  =', spr%yieldDeflection
       write(lpu,*) 'totalDeflection  =', spr%deflection
       write(lpu,*) 'currentForce     =', spr%force
    end if
#endif

    if (spr%isInYield .eqv. wasInYield) return

    lpu = getErrorFile()
    if (spr%isInYield) then
       write(lpu,600) trim(getId(spr%id)),'is yielding',elasticForce
    else
       write(lpu,600) trim(getId(spr%id)),'re-enters elastic state',elasticForce
    end if

600 format(5X,'Spring',A,1X,A,', elastic force =',1PE13.5)

  end subroutine checkYieldLimits


  !!============================================================================
  !> @brief Updates the constant yield force range for a spring.
  !>
  !> @param spr The springtypemodule::springbasetype object to update
  !> @param[in] fLimit New yield force limit, for tension and compression
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 6 Jul 2022

  subroutine setYieldLimit (spr,fLimit)

    type(SpringBaseType), intent(inout) :: spr
    real(dp)            , intent(in)    :: fLimit

    !! --- Logic section ---

    if (associated(spr%yield)) then
       spr%yield%pos%ysf(1) =  fLimit
       spr%yield%neg%ysf(1) = -fLimit
    end if

  end subroutine setYieldLimit


  !!============================================================================
  !> @brief Checks for elastic unloading of cyclic plastic springs.
  !>
  !> @param spr The springtypemodule::springbasetype object to check for
  !> @param[out] ierr Error flag
  !>
  !> @details The associated plastic (yield) deformation or secant stiffness
  !> are updated if still in plastic range.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 20 Jan 2011

  subroutine checkUnloading (spr,ierr)

    use IdTypeModule      , only : getId
    use FunctionTypeModule, only : FunctionValue
    use reportErrorModule , only : reportError, error_p

    type(SpringBaseType), intent(inout) :: spr
    integer             , intent(inout) :: ierr

    !! Local variables
    integer  :: lerr
    real(dp) :: elasticDeflection, fVal

    !! --- Logic section ---

    if (spr%hasFailed .or. spr%unLoadType < 1) return
    if (.not. associated(spr%forceFunction)) return

    lerr = ierr
    if (spr%unLoadType == 1) then

       elasticDeflection = spr%deflection - spr%yieldDeflection
       fVal = spr%s1*FunctionValue(spr%forceFunction,spr%deflection,ierr)
       if (spr%s0*elasticDeflection > fVal .and. elasticDeflection > eps_p) then
          spr%yieldDeflection = spr%deflection - fVal/spr%s0
       end if

    else if (spr%unLoadType == 3) then

       fVal = spr%s1*FunctionValue(spr%forceFunction,spr%deflection,ierr)
       if ( (spr%s0*spr%deflection > fVal .and. spr%deflection > eps_p) .or. &
            (spr%s0*spr%deflection < fVal .and. spr%deflection < -eps_p) ) then
          !! Beyond the linear range ==> update the secant stiffness
          spr%s0 = fVal/spr%deflection
       end if

    end if

    if (ierr < lerr) then
       call reportError (error_p,'Spring'//getId(spr%id), &
            &            'Evaluating plastic force function failed.', &
            &            addString='checkUnloading')
    end if

  end subroutine checkUnloading

end module SpringTypeModule
