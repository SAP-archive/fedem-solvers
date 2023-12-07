!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file initiateDamperTypeModule.f90
!> @brief Initialization of damper objects from the solver input file.

!!==============================================================================
!> @brief Initialization of damper objects from the solver input file.

module InitiateDamperTypeModule

  implicit none

contains

  !!============================================================================
  !> @brief Initializes damper objects with data from the solver input file.
  !>
  !> @param[in] infp File unit number for the solver input file
  !> @param[in] triads Array of all triads in the model
  !> @param[in] engines Array of all engines in the model
  !> @param[in] functions Array of all functions in the model
  !> @param[in] springs Array of all spring elements in the model
  !> @param baseDampers Array of all base damper objects in the model
  !> @param dampers Array of all damper elements in the model
  !> @param[out] err Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Karl Erik Thoresen
  !> @date Nov 1998
  !>
  !> @author Bjorn Haugen
  !> @date Apr 2000
  !>
  !> @author Knut Morten Okstad
  !> @date 1 Jul 2002

  subroutine InitiateDampers (infp,triads,engines,functions, &
       &                      springs,baseDampers,dampers,err)

    use IdTypeModule          , only : ldesc_p, initId, getId, ReportInputError
    use TriadTypeModule       , only : TriadType, GetPtrToId, allocateNodeForce
    use SpringTypeModule      , only : SpringType, GetPtrToId, dp
    use SpringTypeModule      , only : SpringBaseType, NullifySpring
    use DamperTypeModule      , only : DamperType, DamperBaseType, GetPtrToId
    use DamperTypeModule      , only : DeallocateDampers, NullifyDamper
    use FunctionTypeModule    , only : EngineType, FunctionType
    use inputUtilities        , only : iuGetNumberOfEntries
    use inputUtilities        , only : iuSetPosAtNextEntry, iuCharToInt
    use progressModule        , only : lterm
    use reportErrorModule     , only : AllocationError, reportError
    use reportErrorModule     , only : note_p, warning_p, debugFileOnly_p
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getbool

    integer                , intent(in)  :: infp
    type(TriadType)        , intent(in)  :: triads(:)
    type(EngineType)       , intent(in)  :: engines(:)
    type(FunctionType)     , intent(in)  :: functions(:)
    type(SpringType)       , intent(in)  :: springs(:)
    type(DamperBaseType)   , pointer     :: baseDampers(:)
    type(DamperType)       , pointer     :: dampers(:)
    integer                , intent(out) :: err

    !! Local variables
    integer  :: i, idIn, nDampers, nDmpBase, stat, sprTriad(2)
    real(dp) :: fDir(3)
    logical  :: allSec, allDamp, allCoeff, allLength
    logical  :: allVel, allForce, allEnergy

    type(SpringBaseType), save, target :: dummySpring

    !! Define the DAMPER and DAMPER_BASE namelists
    character(ldesc_p) :: extDescr
    integer  :: id, extId(10), dof, triad1Id, triad2Id
    integer  :: coeffScaleEngineId, coeffFuncId, forceFuncId, saveVar(5)
    real(dp) :: d0, d1
    logical  :: isDefDamper
    namelist /DAMPER/ id, extId, extDescr, dof, triad1Id, triad2Id, &
         &            coeffScaleEngineId, coeffFuncId, forceFuncId, &
         &            d0, d1, saveVar, isDefDamper
    namelist /DAMPER_BASE/ id, extId, extDescr, &
         &                 coeffScaleEngineId, coeffFuncId, forceFuncId, &
         &                 d0, d1, saveVar, isDefDamper

    !! --- Logic section ---

    nDmpBase = iuGetNumberOfEntries(infp,'&DAMPER_BASE',err)
    if (err /= 0) return
    write(lterm,*) 'Number of &DAMPER_BASE =',nDmpBase

    nDampers = iuGetNumberOfEntries(infp,'&DAMPER',err)
    if (err /= 0) return
    write(lterm,*) 'Number of &DAMPER =',nDampers

    call DeallocateDampers (dampers)
    call DeallocateDampers (baseDampers)
    allocate(dampers(nDampers),baseDampers(nDampers+nDmpBase),STAT=stat)
    if (stat /= 0) then
       err = AllocationError('initiateDampers 10')
       return
    end if

    call ffa_cmdlinearg_getbool ('allSecondaryVars',allSec)
    call ffa_cmdlinearg_getbool ('allDamperVars',allDamp)
    call ffa_cmdlinearg_getbool ('allDampCoeff',allCoeff)
    call ffa_cmdlinearg_getbool ('allLengthVars',allLength)
    call ffa_cmdlinearg_getbool ('allVelVars',allVel)
    call ffa_cmdlinearg_getbool ('allForceVars',allForce)
    call ffa_cmdlinearg_getbool ('allEnergyVars',allEnergy)

    call NullifySpring (dummySpring)

    do idIn = 1, nDampers

       call NullifyDamper (dampers(idIn))
       call NullifyDamper (baseDampers(idIn))
       if (.not. iuSetPosAtNextEntry(infp,'&DAMPER')) then
          err = err - 1
          call ReportInputError ('DAMPER',idIn)
          cycle
       end if

       !! Default values
       id=0; extId=0; extDescr=''; dof=0; triad1Id=0; triad2Id=0
       coeffScaleEngineId=0; coeffFuncId=0; forceFuncId=0
       d0=0.0_dp; d1=0.0_dp; saveVar=0; isDefDamper=.false.

       read(infp,nml=DAMPER,iostat=stat)
       if (stat /= 0) then
          err = err - 1
          call ReportInputError ('DAMPER',idIn)
          cycle
       end if

       dampers(idIn)%dmp => baseDampers(idIn)
       call initId (dampers(idIn)%dmp%id,id,extId,extDescr,stat)
       dampers(idIn)%dmp%dmp0 = d0
       dampers(idIn)%dmp%dmp1 = d1

       !! Connect the specified engines and functions to this damper
       stat = err
       call ConnectFunctions (coeffFuncId,forceFuncId,coeffScaleEngineId, &
            &                 dampers(idIn)%dmp)
       if (err < stat) then
          call ReportInputError ('DAMPER',idIn,dampers(idIn)%dmp%id)
          cycle
       end if

       if (triad1Id > 0 .and. triad2Id > 0) then

          !! Axial damper
          dampers(idIn)%triad1 => GetPtrToId(triads,triad1Id)
          dampers(idIn)%triad2 => GetPtrToId(triads,triad2Id)
          if (.not. associated(dampers(idIn)%triad1) .or. &
              .not. associated(dampers(idIn)%triad2)) then
             err = err - 1
             call ReportInputError ('DAMPER',idIn,dampers(idIn)%dmp%id)
             cycle
          else if (dampers(idIn)%triad1%nDOFs < dampers(idIn)%triad2%nDOFs) then
             !! Swap triads such that the one with most DOFs always is triad1
             dampers(idIn)%triad1 => dampers(idIn)%triad2
             dampers(idIn)%triad2 => GetPtrToId(triads,triad1Id)
          end if

          allocate(dampers(idIn)%dmp%length, &
               &   dampers(idIn)%dmp%velocity, &
               &   dampers(idIn)%forceDir(3), STAT=stat)
          if (stat /= 0) then
             err = AllocationError('initiateDampers 20')
             return
          end if

          fDir = dampers(idIn)%triad2%ur(:,4) - dampers(idIn)%triad1%ur(:,4)
          dampers(idIn)%dmp%length   = sqrt(dot_product(fDir,fDir))
          dampers(idIn)%dmp%velocity = 0.0_dp
          dampers(idIn)%forceDir     = 0.0_dp
          dampers(idIn)%nDofs        = dampers(idIn)%triad1%nDofs &
               &                     + dampers(idIn)%triad2%nDofs

          !! Set samNodNum temporarily to flag that the triad should receive
          !! a proper node number later, also when the triad is grounded
          dampers(idIn)%triad1%samNodNum = 1
          dampers(idIn)%triad2%samNodNum = 1

       else if (triad1Id > 0) then

          !! Damper to ground
          dampers(idIn)%triad1 => GetPtrToId(triads,triad1Id)
          if (.not. associated(dampers(idIn)%triad1)) then
             err = err - 1
             call ReportInputError ('DAMPER',idIn,dampers(idIn)%dmp%id)
             cycle
          end if

          allocate(dampers(idIn)%dmp%length,STAT=stat)
          if (stat /= 0) then
             err = AllocationError('initiateDampers 21')
             return
          end if

          dampers(idIn)%dmp%length   = 0.0_dp
          dampers(idIn)%dmp%velocity =>dampers(idIn)%triad1%urd(dof)
          dampers(idIn)%nDofs        = dampers(idIn)%triad1%nDofs
          dampers(idIn)%dmp%dof      = dof

       else if (triad2Id > 0 .or. triad1Id < 0) then
          err = err - 1
          call ReportInputError ('DAMPER',idIn,dampers(idIn)%dmp%id)
          cycle
       end if

       !! If the damper triads also are connected to a superelement,
       !! allocate separate arrays for storage of the total nodal force
       stat = err
       if (associated(dampers(idIn)%triad1)) then
          call allocateNodeForce (dampers(idIn)%triad1,err)
          if (err < stat) exit
       end if
       if (associated(dampers(idIn)%triad2)) then
          call allocateNodeForce (dampers(idIn)%triad2,err)
          if (err < stat) exit
       end if

       !! Consistency check: Deactivate all completely grounded dampers
       if (dampers(idIn)%nDofs > 0) then
          dampers(idIn)%dmp%isActive = .true.
       else if (triad1Id > 0) then
          call reportError (warning_p, &
               'Axial Damper'//trim(getId(dampers(idIn)%dmp%id))// &
               ' is only connected to grounded Triads (ignored).')
       end if

       if (isDefDamper) then

          !! This is a deformational damper.
          !! Connect to axial spring using the same triads if axial damper.
          if (triad1Id > 0 .and. triad2Id > 0) then
             do i = 1, size(springs)

                if (size(springs(i)%triads) == 2) then
                   sprTriad(1) = springs(i)%triads(1)%p%id%baseId
                   sprTriad(2) = springs(i)%triads(2)%p%id%baseId
                   if ( all(sprTriad == (/triad1Id,triad2Id/)) .or. &
                        all(sprTriad == (/triad2Id,triad1Id/)) ) then
                      dampers(idIn)%dmp%spr => springs(i)%spr(1)%p
                      exit
                   end if
                end if

             end do
          else if (triad1Id == 0) then

             !! Temporarily assign the spring pointer to a dummy object.
             !! The actual associated joint spring will be connected later
             !! by the ReadJoints::InitJointDof subroutine.
             dampers(idIn)%dmp%spr => dummySpring

          end if

          if (.not. associated(dampers(idIn)%dmp%spr)) then
             call reportError (note_p, &
                  'Switching off deformational damping for Axial Damper'// &
                  trim(getId(dampers(idIn)%dmp%id))// &
                  ' since no associated spring')
          end if

       end if

       !! Determine which secondary damper variables will be saved. Settings in
       !! the solver input file are overruled by command-line arguments, if any.
       call setSaveVar (baseDampers(idIn),saveVar)

    end do
    rewind(infp)
    do idIn = nDampers+1, nDampers+nDmpBase

       call NullifyDamper (baseDampers(idIn))
       if (.not. iuSetPosAtNextEntry(infp,'&DAMPER_BASE')) then
          err = err - 1
          call ReportInputError ('DAMPER_BASE',idIn)
          cycle
       end if

       !! Default values
       id=0; extId=0; extDescr=''
       coeffScaleEngineId=0; coeffFuncId=0; forceFuncId=0
       d0=0.0_dp; d1=0.0_dp; saveVar=0; isDefDamper=.false.

       read(infp,nml=DAMPER_BASE,iostat=stat)
       if (stat /= 0) then
          err = err - 1
          call ReportInputError ('DAMPER_BASE',idIn)
          cycle
       end if

       call initId (baseDampers(idIn)%id,id,extId,extDescr,stat)
       baseDampers(idIn)%dmp0 = d0
       baseDampers(idIn)%dmp1 = d1

       !! Connect the specified engines and functions to this damper
       stat = err
       call ConnectFunctions (coeffFuncId,forceFuncId,coeffScaleEngineId, &
            &                 baseDampers(idIn))
       if (err < stat) then
          call ReportInputError ('DAMPER_BASE',idIn,baseDampers(idIn)%id)
       end if

       if (isDefDamper) then
          !! Temporarily assign the spring pointer to a dummy object.
          !! The actual associated joint spring will be connected later
          !! by the ReadContactElements subroutine.
          baseDampers(idIn)%spr => dummySpring
       end if

       !! Determine which secondary damper variables will be saved. Settings in
       !! the solver input file are overruled by command-line arguments, if any.
       call setSaveVar (baseDampers(idIn),saveVar)

    end do

    if (err < 0) call reportError (debugFileOnly_p,'InitiateDampers')

  contains

    !> @brief Connects engines and functions to the given damper object.
    subroutine ConnectFunctions (coeffId,forceId,scaleId,dmp)
      use FunctionTypeModule, only : GetPtrToId

      integer             , intent(in)    :: coeffId, forceId, scaleId
      type(DamperBaseType), intent(inout) :: dmp

      if (coeffId > 0) then
         !! Connect the coefficient-velocity function
         dmp%coeffFunction => GetPtrToId(functions,coeffId)
         if (.not. associated(dmp%coeffFunction)) err = err - 1
      else if (forceId > 0) then
         !! Connect the force-velocity function
         dmp%forceFunction => GetPtrToId(functions,forceId)
         if (.not. associated(dmp%forceFunction)) err = err - 1
      end if

      if (scaleId > 0) then
         !! Connect the coefficient scale engine
         dmp%coeffScaleEngine => GetPtrToId(engines,scaleId)
         if (.not. associated(dmp%coeffScaleEngine)) err = err - 1
      end if

    end subroutine ConnectFunctions

    !> @brief Determines which secondary damper variables will be saved.
    subroutine setSaveVar (dmp,var)
      type(DamperBaseType), intent(inout) :: dmp
      integer             , intent(in)    :: var(:)
      do i = 1, size(saveVar)
         dmp%saveVar(i) = var(i) > 0 .or. allSec .or. allDamp
      end do
      if (allCoeff)  dmp%saveVar(1) = .true. ! Damper coefficient
      if (allLength) dmp%saveVar(2) = .true. ! Damper length
      if (allVel)    dmp%saveVar(3) = .true. ! Damper velocity
      if (allForce)  dmp%saveVar(4) = .true. ! Damper force
      if (allEnergy) dmp%saveVar(5) = .true. ! Damper energy
    end subroutine setSaveVar

  end subroutine InitiateDampers

end module InitiateDamperTypeModule
