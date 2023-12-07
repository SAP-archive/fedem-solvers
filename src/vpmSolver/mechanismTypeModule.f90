!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file mechanismTypeModule.f90
!> @brief Mechanism data container.

!!==============================================================================
!> @brief Module with mechanism data containers.
!>
!> @details This module contains data containers for all components that
!> together define the mechanism model that is to be analyzed. It also contains
!> some subroutines for accessing or manipulating the mechanism data.

module MechanismTypeModule

  use EnvironmentTypeModule     , only : EnvironmentType
  use TriadTypeModule           , only : TriadType, IdType, dp
  use SupElTypeModule           , only : SupElType
  use SupElLoadTypeModule       , only : SupElLoadType
  use MassTypeModule            , only : MassType
  use ForceTypeModule           , only : ForceType
  use MotionTypeModule          , only : MotionType
  use SpringTypeModule          , only : SpringBaseType, SpringType
  use SpringTypeModule          , only : SpringFailureType, SpringYieldType
  use DamperTypeModule          , only : DamperBaseType, DamperType
  use BushingElementTypeModule  , only : BushingElementType
  use ContactElementTypeModule  , only : ContactElementType
  use ContactSurfaceModule      , only : GliderCurveType
  use MasterSlaveJointTypeModule, only : MasterSlaveJointType, HigherPairType
  use FrictionTypeModule        , only : FrictionParameterType, FrictionPtrType
  use FunctionTypeModule        , only : FunctionType, EngineType
  use SensorTypeModule          , only : SensorType
  use RoadTypeModule            , only : RoadType
  use TireTypeModule            , only : TireType
  use UserdefElTypeModule       , only : UserdefElType
  use WindTurbineTypeModule     , only : TurbineConfig

  implicit none

  !> @brief Data type containing all components of the mechanism model.
  type MechanismType

     type(IdType) :: id !< General identification data

     type(EnvironmentType) :: env        !< Data for environmental loads
     real(dp), pointer     :: gravity(:) !< Should point to @a env%gravity

     real(dp) :: Etot !< Energy checksum
     real(dp) :: Epot !< Potential energy
     real(dp) :: Ekin !< Kinetic energy
     real(dp) :: Estr !< Strain energy
     real(dp) :: Einp !< Input energy
     real(dp) :: Edmp !< Damping energy
     real(dp) :: Eext !< Mixed energies from external "objects"

     logical :: hasEinp !< If .true., the mechanism have input energy
     logical :: hasEext !< If .true., the mechanism have external energy

     logical :: saveVar(6) !< Flags indicating variables to be saved

     !! All mechanism components
     type(TriadType)            , pointer :: triads(:)
     type(SupElType)            , pointer :: sups(:)
     type(SupElLoadType)        , pointer :: supLoads(:)
     type(UserdefElType)        , pointer :: elms(:)
     type(MassType)             , pointer :: masses(:)
     type(ForceType)            , pointer :: forces(:)
     type(MotionType)           , pointer :: motions(:)
     type(SpringBaseType)       , pointer :: baseSprings(:)
     type(SpringType)           , pointer :: axialSprings(:)
     type(SpringFailureType)    , pointer :: springFailures(:)
     type(SpringYieldType)      , pointer :: springYields(:)
     type(DamperBaseType)       , pointer :: baseDampers(:)
     type(DamperType)           , pointer :: dampers(:)
     type(BushingElementType)   , pointer :: bElems(:)
     type(ContactElementType)   , pointer :: cElems(:)
     type(GliderCurveType)      , pointer :: cSurfs(:)
     type(MasterSlaveJointType) , pointer :: joints(:)
     type(HigherPairType)       , pointer :: higherPairs(:)
     type(FrictionParameterType), pointer :: frictionSets(:)
     type(FrictionPtrType)      , pointer :: frictions(:)
     type(FunctionType)         , pointer :: functions(:)
     type(EngineType)           , pointer :: engines(:)
     type(SensorType)           , pointer :: sensors(:)
     type(RoadType)             , pointer :: roads(:)
     type(TireType)             , pointer :: tires(:)
     type(TurbineConfig)        , pointer :: turbine

  end type MechanismType

  integer, save :: mDmpFuncId = 0 !< Base ID of modal damping function


  !> @brief Standard routine for writing an object to file.
  interface WriteObject
     module procedure WriteMechanismType
  end interface


contains

  !!==========================================================================
  !> @brief Standard routine for writing an object to io.
  !>
  !> @param[in] mech Mechanism components of the model
  !> @param[in] io File unit number to write to
  !> @param[in] complexity Option controlling the amount of output
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 27 May 2002

  subroutine WriteMechanismType (mech,io,complexity)

    use IdTypeModule              , only : writeId
    use EnvironmentTypeModule     , only : writeObject
    use TriadTypeModule           , only : writeObject
    use SupElTypeModule           , only : writeObject
    use SupElLoadTypeModule       , only : writeObject
    use UserdefElTypeModule       , only : writeObject
    use MassTypeModule            , only : writeObject
    use ForceTypeModule           , only : writeObject
    use MotionTypeModule          , only : writeObject
    use SpringTypeModule          , only : writeObject
    use DamperTypeModule          , only : writeObject
    use BushingElementTypeModule  , only : writeObject
    use ContactElementTypeModule  , only : writeObject
    use MasterSlaveJointTypeModule, only : writeObject
    use FrictionTypeModule        , only : writeObject
    use FunctionTypeModule        , only : writeObject
    use SensorTypeModule          , only : writeObject
    use RoadTypeModule            , only : writeObject
    use TireTypeModule            , only : writeObject
    use WindTurbineTypeModule     , only : writeObject

    type(MechanismType), intent(in) :: mech
    integer            , intent(in) :: io
    integer, optional  , intent(in) :: complexity

    !! Local variables
    integer :: i

    !! --- Logic section ---

    write(io,'(A)') 'Mechanism','{'
    call writeId (mech%id,io)

    write(io,'(A,1P,3E13.5)') ' gravity     =', mech%gravity
    write(io,*) 'saveVar     =', mech%saveVar

    if (present(complexity)) then
       if (complexity >= 2) then
          write(io,'(A,1P,7E13.5)') ' energies    =', mech%Etot,mech%Epot, &
               &                                      mech%Ekin,mech%Estr, &
               &                                      mech%Einp,mech%Edmp, &
               &                                      mech%Eext
       end if
       if (complexity >= 3 .or. complexity < 0) then
          write(io,*) 'size(triads)       =', size(mech%triads)
          write(io,*) 'size(sups)         =', size(mech%sups)
          write(io,*) 'size(supLoads)     =', size(mech%supLoads)
          write(io,*) 'size(elms)         =', size(mech%elms)
          write(io,*) 'size(masses)       =', size(mech%masses)
          write(io,*) 'size(forces)       =', size(mech%forces)
          write(io,*) 'size(motions)      =', size(mech%motions)
          write(io,*) 'size(baseSprings)  =', size(mech%baseSprings)
          write(io,*) 'size(axialSprings) =', size(mech%axialSprings)
          write(io,*) 'size(springFailurs)=', size(mech%springFailures)
          write(io,*) 'size(springYields) =', size(mech%springYields)
          write(io,*) 'size(baseDampers)  =', size(mech%baseDampers)
          write(io,*) 'size(dampers)      =', size(mech%dampers)
          write(io,*) 'size(bElems)       =', size(mech%bElems)
          write(io,*) 'size(cElems)       =', size(mech%cElems)
          write(io,*) 'size(cSurfs)       =', size(mech%cSurfs)
          write(io,*) 'size(joints)       =', size(mech%joints)
          write(io,*) 'size(higherPairs)  =', size(mech%higherPairs)
          write(io,*) 'size(frictionSets) =', size(mech%frictionSets)
          write(io,*) 'size(frictions)    =', size(mech%frictions)
          write(io,*) 'size(functions)    =', size(mech%functions)
          write(io,*) 'size(engines)      =', size(mech%engines)
          write(io,*) 'size(sensors)      =', size(mech%sensors)
          write(io,*) 'size(roads)        =', size(mech%roads)
          write(io,*) 'size(tires)        =', size(mech%tires)
       end if
    end if

    write(io,'(A/)') '}'

    call writeObject (mech%env,io)

    do i = 1, size(mech%triads)
       call WriteObject (mech%triads(i),io,complexity)
       write (io,'(/)')
    end do

    do i = 1, size(mech%sups)
       call WriteObject (mech%sups(i),io,complexity)
       write (io,'(/)')
    end do

    do i = 1, size(mech%supLoads)
       call WriteObject (mech%supLoads(i),io,complexity)
       write (io,'(/)')
    end do

    do i = 1, size(mech%elms)
       call WriteObject (mech%elms(i),io,complexity)
       write (io,'(/)')
    end do

    do i = 1, size(mech%masses)
       call writeObject (mech%masses(i),io,complexity)
       write (io,'(/)')
    end do

    do i = 1, size(mech%forces)
       call writeObject (mech%forces(i),io,complexity)
       write (io,'(/)')
    end do

    do i = 1, size(mech%motions)
       call writeObject (mech%motions(i),io,complexity)
       write (io,'(/)')
    end do

    do i = 1, size(mech%axialSprings)
       call writeObject (mech%axialSprings(i),io,complexity)
       write (io,'(/)')
    end do

    do i = 1, size(mech%dampers)
       call writeObject (mech%dampers(i),io,complexity)
       write (io,'(/)')
    end do

    do i = 1, size(mech%bElems)
       call writeObject (mech%bElems(i),io,complexity)
       write (io,'(/)')
    end do

    do i = 1, size(mech%cElems)
       call writeObject (mech%cElems(i),io,complexity)
       write (io,'(/)')
    end do

    do i = 1, size(mech%joints)
       call writeObject (mech%joints(i),io,complexity)
       write (io,'(/)')
    end do

    do i = 1, size(mech%frictions)
       call writeObject (mech%frictions(i)%p,io,complexity)
       write (io,'(/)')
    end do

    do i = 1, size(mech%functions)
       call writeObject (mech%functions(i),io)
       write (io,'(/)')
    end do

    do i = 1, size(mech%engines)
       call writeObject (mech%engines(i),io)
       write (io,'(/)')
    end do

    do i = 1, size(mech%sensors)
       call writeObject (mech%sensors(i),io,complexity)
       write (io,'(/)')
    end do

    do i = 1, size(mech%roads)
       call writeObject (mech%roads(i),io)
       write (io,'(/)')
    end do

    do i = 1, size(mech%tires)
       call writeObject (mech%tires(i),io,complexity)
       write (io,'(/)')
    end do

    if (associated(mech%turbine)) then
       call writeObject (mech%turbine,io)
       write (io,'(/)')
    end if

  end subroutine WriteMechanismType


  !!============================================================================
  !> @brief Initializes the MechanismType object.
  !>
  !> @param[out] mech Mechanism components of the model
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 27 May 2002

  subroutine NullifyMechanism (mech)

    use IdTypeModule         , only : nullifyId
    use EnvironmentTypeModule, only : nullifyEnvironment

    type(MechanismType), target, intent(out) :: mech

    !! --- Logic section ---

    call nullifyId (mech%id)
    call nullifyEnvironment (mech%env)

    mech%gravity => mech%env%gravity

    mech%Etot = 0.0_dp
    mech%Epot = 0.0_dp
    mech%Ekin = 0.0_dp
    mech%Estr = 0.0_dp
    mech%Einp = 0.0_dp
    mech%Edmp = 0.0_dp
    mech%Eext = 0.0_dp
    mech%hasEinp = .false.
    mech%hasEext = .false.
    mech%saveVar = .false.

    !! Allocate initially to zero length, such that the size() operator works
    allocate(mech%triads(0), mech%sups(0), mech%supLoads(0), mech%elms(0), &
         &   mech%masses(0), mech%forces(0), mech%motions(0), &
         &   mech%baseSprings(0), mech%axialSprings(0), &
         &   mech%springFailures(0), mech%springYields(0), &
         &   mech%baseDampers(0), mech%dampers(0), &
         &   mech%bElems(0), mech%cElems(0), mech%cSurfs(0), &
         &   mech%joints(0), mech%higherPairs(0), &
         &   mech%frictionSets(0), mech%frictions(0), &
         &   mech%functions(0), mech%engines(0), mech%sensors(0), &
         &   mech%roads(0), mech%tires(0))

    nullify(mech%turbine)

  end subroutine NullifyMechanism


  !!============================================================================
  !> @brief Deallocates the MechanismType object.
  !>
  !> @param mech Mechanism components of the model
  !>
  !> @callergraph
  !>
  !> @details The various mechanism object types are deallocated in the opposite
  !> order as they were allocated, in order to resolve any interconnectivities
  !> safely.
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Jan 2017

  subroutine DeallocateMechanism (mech)

    use IdTypeModule              , only : DeallocateId
    use BushingElementTypeModule  , only : DeallocateBushingElements
    use ContactElementTypeModule  , only : DeallocateContactElements
    use ContactSurfaceModule      , only : DeallocateGliderCurves
    use DamperTypeModule          , only : DeallocateDampers
    use ForceTypeModule           , only : DeallocateForces
    use FrictionTypeModule        , only : DeallocateFrictionPrms
    use FunctionTypeModule        , only : DeallocateFunctions
    use MassTypeModule            , only : DeallocateMasses
    use MasterSlaveJointTypeModule, only : DeallocateJoints
    use MotionTypeModule          , only : DeallocateMotions
    use RoadTypeModule            , only : DeallocateRoads
    use SensorTypeModule          , only : DeallocateSensors
    use SpringTypeModule          , only : DeallocateSprings
    use SupElLoadTypeModule       , only : DeallocateLoads
    use SupElTypeModule           , only : DeallocateSupEls
    use TireTypeModule            , only : DeallocateTires
    use TriadTypeModule           , only : DeallocateTriads
    use UserdefElTypeModule       , only : DeallocateUdElms
    use WindTurbineTypeModule     , only : DeallocateTurbine

    type(MechanismType), intent(inout) :: mech

    !! --- Logic section ---

    if (associated(mech%turbine)) then
       call deallocateTurbine (mech%turbine)
       deallocate(mech%turbine)
       nullify(mech%turbine)
    end if

    call DeallocateSensors (mech%sensors)

    deallocate(mech%frictions)
    nullify(mech%frictions)

    call DeallocateTires (mech%tires)
    call DeallocateRoads (mech%roads)

    call DeallocateBushingElements (mech%bElems)
    call DeallocateContactElements (mech%cElems)
    call DeallocateGliderCurves (mech%cSurfs)

    call DeallocateMotions (mech%motions)
    call DeallocateForces (mech%forces)
    call DeallocateLoads (mech%supLoads)

    call DeallocateJoints (mech%higherPairs)
    call DeallocateJoints (mech%joints)

    call DeallocateDampers (mech%dampers)
    call DeallocateDampers (mech%baseDampers)
    call DeallocateSprings (mech%axialSprings)
    call DeallocateSprings (mech%baseSprings)
    call DeallocateSprings (mech%springYields)
    call DeallocateSprings (mech%springFailures)
    call DeallocateMasses (mech%masses)

    call DeallocateFrictionPrms (mech%frictionSets)

    call DeallocateFunctions (mech%engines)
    call DeallocateFunctions (mech%functions)

    call DeallocateUdElms (mech%elms)
    call DeallocateSupEls (mech%sups)
    call DeallocateTriads (mech%triads)

    call DeallocateId (mech%id)

  end subroutine DeallocateMechanism


  !!============================================================================
  !> @brief Initializes the mechanism with data from the solver input file.
  !>
  !> @param[in] infp File unit number for the solver input file
  !> @param[out] mech Mechanism components of the model
  !> @param[out] tIncEngineID Base ID of the function returning time step size
  !> @param[out] weight Weighting factors for the DOFs of different kind
  !> @param[out] err Error flag
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 27 May 2002

  subroutine ReadMechanism (infp,mech,tIncEngineId,weight,err)

    use IdTypeModule          , only : ldesc_p, initId, getId, ReportInputError
    use inputUtilities        , only : iuSetPosAtFirstEntry
    use progressModule        , only : lterm
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getbool

    integer            , intent(in)  :: infp
    type(MechanismType), intent(out) :: mech
    integer            , intent(out) :: tIncEngineId
    real(dp)           , intent(out) :: weight(3)
    integer            , intent(out) :: err

    !! Local variables
    logical :: allSec, allRest, allSys, allCG, allHD, allEnergy, allAlgor

    !! Define the MECHANISM namelist
    integer            :: id, extId(10), &
         &                timeIncrEngine, modalDmpFunction, saveVar(3)
    character(ldesc_p) :: extDescr
    real(dp)           :: weightTranslation, weightRotation, weightGeneralized

    namelist /MECHANISM/ id, extId, extDescr, &
         &               weightTranslation, weightRotation, weightGeneralized, &
         &               timeIncrEngine, modalDmpFunction, saveVar

    !! --- Logic section ---

    !! Default values
    id=1; extId=0; extDescr=''; timeIncrEngine=0; modalDmpFunction=0; saveVar=0

    !! Read Mechanism data from the solver input file
    if (iuSetPosAtFirstEntry(infp,'&MECHANISM')) then
       weightTranslation=0.0_dp; weightRotation=0.0_dp; weightGeneralized=0.0_dp
       read(infp,nml=MECHANISM,iostat=err)
       if (err == 0) then
          weight(1) = weightTranslation
          weight(2) = weightRotation
          weight(3) = weightGeneralized
          call initId (mech%id,id,extId,extDescr,err)
          write(lterm,*) 'Read &MECHANISM',trim(getId(mech%id))
          write(lterm,'(4X,A,1P3E13.5)') 'weigth      =',weight
          if (timeIncrEngine > 0) then
             write(lterm,*) '   timeIncr(id) =', timeIncrEngine
          end if
          if (modalDmpFunction > 0) then
             write(lterm,*) '   modalDmp(id) =', modalDmpFunction
          end if
       else
          call ReportInputError ('MECHANISM',msg='Syntax error')
       end if
    else
       extId(1) = 1
       weight = 1.0_dp
       call initId (mech%id,id,extId,extDescr,err)
    end if

    tIncEngineId = timeIncrEngine
    mDmpFuncId   = modalDmpFunction

    !! Get output options from the command-line (or option files)

    call ffa_cmdlinearg_getbool ('allSecondaryVars',allSec)
    call ffa_cmdlinearg_getbool ('allRestartVars',allRest)
    call ffa_cmdlinearg_getbool ('allSystemVars',allSys)
    call ffa_cmdlinearg_getbool ('allCGVars',allCG)
    call ffa_cmdlinearg_getbool ('allHDVars',allHD)
    call ffa_cmdlinearg_getbool ('allEnergyVars',allEnergy)
    call ffa_cmdlinearg_getbool ('allAlgorVars',allAlgor)

    !! Determine which secondary system variables will be saved. Settings in
    !! the solver input file are overruled by command-line arguments, if any.
    do id = 1, size(saveVar)
       mech%saveVar(id) = saveVar(id) > 0 .or. allSec .or. allSys
    end do
    if (allCG)     mech%saveVar(1) = .true. ! Centre of gravity
    if (allHD)     mech%saveVar(5) = .true. ! Buoyancy centre and volume
    if (allEnergy) mech%saveVar(2) = .true. ! Energy quantities
    if (allAlgor)  mech%saveVar(3) = .true. ! Algorithm parameters
    if (allRest)   mech%saveVar(3) = .true. ! Restart variables
    if (allRest)   mech%saveVar(6) = .true. ! Restart variables (inertia forces)
    mech%saveVar(4) = allSec .or. allSys .or. allAlgor ! Convergence norms

  end subroutine ReadMechanism

end module MechanismTypeModule
