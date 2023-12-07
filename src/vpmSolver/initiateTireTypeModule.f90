!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module InitiateTireTypeModule

  implicit none

contains

  subroutine ReadTires (infp,joints,roads,tires,err)

    !!==========================================================================
    !! Initiates the tire type with data from the solver input file.
    !!
    !! Programmer : Bjorn Haugen                         date/rev : Aug 2002/1.0
    !!==========================================================================

    use TireTypeModule            , only : TireType, interfaces_p, tireTypes_p
    use TireTypeModule            , only : DeallocateTires
    use TireTypeModule            , only : NullifyTire, NullifySTIapi, STI_p
    use RoadTypeModule            , only : RoadType, GetPtrToId, dp, lfnam_p
    use MasterSlaveJointTypeModule, only : MasterSlaveJointType, GetPtrToId
    use MasterSlaveJointTypeModule, only : REVOLUTE_p
    use TriadTypeModule           , only : allocateNodeForce
    use IdTypeModule              , only : ldesc_p, initId, ReportInputError
    use inputUtilities            , only : iuGetNumberOfEntries
    use inputUtilities            , only : iuSetPosAtNextEntry, iuCharToInt
    use progressModule            , only : lterm
    use reportErrorModule         , only : AllocationError
    use reportErrorModule         , only : reportError, debugFileOnly_p
    use FFaFilePathInterface      , only : FFa_checkPath
    use FFaCmdLineArgInterface    , only : ffa_cmdlinearg_getbool
    use FFaCmdLineArgInterface    , only : ffa_cmdlinearg_getint

    integer                   , intent(in)  :: infp
    type(MasterSlaveJointType), intent(in)  :: joints(:)
    type(RoadType)            , intent(in)  :: roads(:)
    type(TireType)            , pointer     :: tires(:)
    integer                   , intent(out) :: err

    !! Local Variables
    integer :: i, idIn, nTires, stat, overrideTireType
    logical :: allSec, allTire, allLength, allDef, allVel, allForce, allEnergy
    type(MasterSlaveJointType), pointer, save :: pJoint

    !! Define the TIRE namelist
    character(ldesc_p) :: extDescr
    integer            :: id, extId(10), roadId, jointId, tireChar, saveVar(9)
    integer            :: WCYalongZ
    character(len=30)  :: api, type
    character(lfnam_p) :: tireDataFile
    real(dp)           :: Zoffset, radialStiff, radialDamp
    namelist /TIRE/ id, extId, extDescr, api, type, tireDataFile, &
         &          roadId, jointId, tireChar, Zoffset, &
         &          radialStiff, radialDamp, saveVar, WCYalongZ

    !! --- Logic section ---

    nTires = iuGetNumberOfEntries(infp,'&TIRE',err)
    if (err /= 0) return
    write(lterm,*) 'Number of &TIRE =',nTires

    call DeallocateTires (tires)
    allocate(tires(nTires),STAT=stat)
    if (stat /= 0) then
       err = AllocationError('ReadTires 10')
       return
    end if

    call ffa_cmdlinearg_getbool ('allSecondaryVars',allSec)
    call ffa_cmdlinearg_getbool ('allTireVars',allTire)
    call ffa_cmdlinearg_getbool ('allLengthVars',allLength)
    call ffa_cmdlinearg_getbool ('allDefVars',allDef)
    call ffa_cmdlinearg_getbool ('allVelVars',allVel)
    call ffa_cmdlinearg_getbool ('allForceVars',allForce)
    call ffa_cmdlinearg_getbool ('allEnergyVars',allEnergy)
    call ffa_cmdlinearg_getint  ('overrideTireType',overrideTireType)

    do idIn = 1, nTires

       call NullifyTire (tires(idIn))
       if (.not. iuSetPosAtNextEntry(infp,'&TIRE')) then
          err = err - 1
          call ReportInputError ('TIRE',idIn)
          cycle
       end if

       !! Default values
       id=0; extId=0; extDescr=''
       tireDataFile=''; api=''; type=''
       roadId=0; jointId=0; Zoffset=0.0_dp
       radialStiff=0.0_dp; radialDamp=0.0_dp; saveVar=0
       WCYalongZ=1

       read(infp,nml=TIRE,iostat=stat)
       if (stat /= 0) then
          err = err -1
          call ReportInputError ('TIRE',idIn)
          cycle
       end if

       call initId (tires(idIn)%id,id,extId,extDescr,stat)
       tires(idIn)%idIn = idIn
       tires(idIn)%api  = iuCharToInt(api,interfaces_p,stat)
       if (overrideTireType > 0) then
          tires(idIn)%type = overrideTireType
       else
          tires(idIn)%type = iuCharToInt(type,tireTypes_p,stat)
       end if
       if (stat < 0) then
          err = err - 1
          call ReportInputError ('TIRE',idIn,tires(idIn)%id)
          cycle
       end if

       call FFa_checkPath (tireDataFile)
       tires(idIn)%tireDataFileName      = trim(tireDataFile)
       tires(idIn)%nCharTireDataFileName = len_trim(tireDataFile)
       tires(idIn)%Zoffset   = Zoffset
       tires(idIn)%WCYalongZ = WCYalongZ == 1
       tires(idIn)%hasEccVec = abs(ZoffSet) > tiny(1.0_dp)
       tires(idIn)%rStiff(1) = radialStiff
       tires(idIn)%rDamp     = radialDamp

       if (tires(idIn)%api == STI_p) then
          allocate(tires(idIn)%sti,STAT=stat)
          if (stat /= 0) then
             err = AllocationError('ReadTires 20')
             return
          end if
          call nullifySTIapi (tireChar,tires(idIn)%sti)
       end if

       !! Connect the specified road to this tire
       tires(idIn)%road => GetPtrToId(roads,roadId)
       if (.not. associated(tires(idIn)%road)) then
          err = err - 1
          call ReportInputError ('TIRE',idIn,tires(idIn)%id)
          cycle
       end if

       !! Connect the specified revolute joint to this tire
       pJoint => GetPtrToId(joints,jointId)
       if (.not. associated(pJoint)) then
          err = err - 1
          call ReportInputError ('TIRE',idIn,tires(idIn)%id)
          cycle
       else if (pJoint%type /= REVOLUTE_p) then
          err = err - 1
          call ReportInputError ('TIRE',idIn,tires(idIn)%id, &
               msg='The joint specified is not a revolute joint')
          cycle
       end if

       tires(idIn)%WCtriad  => pJoint%JMTriads(1)%triad
       tires(idIn)%RimTriad => pJoint%Striad
       tires(idIn)%pJvar    => pJoint%jointDofs(1)%jVar

       !! If the rim triad is also connected to a superelement,
       !! allocate a separate array for storage of the total nodal force
       stat = err
       call allocateNodeForce (tires(idIn)%RimTriad,err)
       if (err < stat) exit

       !! Determine which secondary tire variables will be saved. Settings in
       !! the solver input file are overruled by command-line arguments, if any.
       do i = 1, size(saveVar)
          tires(idIn)%saveVar(i) = saveVar(i) > 0 .or. allSec .or. allTire
       end do
       if (allLength) tires(idIn)%saveVar(1) = .true. ! Tire angles
       if (allLength) tires(idIn)%saveVar(2) = .true. ! Tire rolling radius
       if (allForce)  tires(idIn)%saveVar(3) = .true. ! Contact force
       if (allDef)    tires(idIn)%saveVar(5) = .true. ! Tire deflection
       if (allVel)    tires(idIn)%saveVar(6) = .true. ! Deflection velocity
       if (allForce)  tires(idIn)%saveVar(8) = .true. ! Wheel carrier force
       if (allEnergy) tires(idIn)%saveVar(9) = .true. ! Tire energy

    end do
    if (err < 0) call reportError (debugFileOnly_p,'ReadTires')

  end subroutine ReadTires


  subroutine InitiateTires (sys,tires,roads,ierr)

    !!==========================================================================
    !! Initiates all tires.
    !!
    !! Programmer : Knut Morten Okstad                    date/rev: Oct 2002/1.0
    !!==========================================================================

    use SystemTypeModule      , only : SystemType, dp
    use TireTypeModule        , only : TireType, JD_TIRE_p
    use TireRoutinesModule    , only : rtm_gTires, scaleTo
    use TireRoutinesModule    , only : initiateTire, updateTire
    use RoadRoutinesModule    , only : rtm_gRoads, RoadType, scaleToSI
    use profilerModule        , only : startTimer, stopTimer, ini_p, upd_p
    use reportErrorModule     , only : reportError, debugFileOnly_p
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getdouble

    type(SystemType), intent(in)           :: sys
    type(TireType)  , intent(inout),target :: tires(:)
    type(RoadType)  , intent(in), target   :: roads(:)
    integer         , intent(out)          :: ierr

    !! Local variables
    integer :: i, stat, nTires

    integer, parameter :: okCTI_p = 3 ! Steady-state Ftire calculation

    !! --- Logic section ---

    ierr = 0
    nTires = size(tires)
    rtm_gRoads => roads
    rtm_gTires => tires

    !! Initiate scaling factors to SI-units for the tire API
    call ffa_cmdlinearg_getdouble ('scaleToM',scaleTo%m)
    call ffa_cmdlinearg_getdouble ('scaleToS',scaleTo%s)
    call ffa_cmdlinearg_getdouble ('scaleToKG',scaleTo%kg)
    scaleTo%N = scaleTo%kg * scaleTo%m/(scaleTo%s*scaleTo%s)
    scaleToSI = scaleTo%m ! Used by the road call-back

    do i = 1, nTires

       call InitiateTire (tires(i),sys%tStart,stat)
       if (stat < 0) then
          ierr = stat
          exit
       else if (stat > 0) then
          ierr = ierr - 1
          cycle
       end if

       call stopTimer (ini_p)  ! We are currently measuring initialization time,
       call startTimer (upd_p) ! but tire calculations should be in the update.

       !! Compute the initial tire response (to avoid 0,0 points in curve plots)
       !!TODO,bh: should one have selected integration scheme here as well?
       call updateTire (1,okCTI_p,sys,tires(i),ierr)

       call stopTimer (upd_p)  ! Stop the update timer
       call startTimer (ini_p) ! Restart the initialization timer again

       !! If this is John Deere tires, ignore the scaling to SI

       if (tires(i)%type == JD_TIRE_p) then
          scaleTo%m  = 1.0_dp
          scaleTo%s  = 1.0_dp
          scaleTo%kg = 1.0_dp
          scaleTo%N  = 1.0_dp
          scaleToSI  = 1.0_dp
       end if

    end do

    if (ierr < 0) call reportError (debugFileOnly_p,'InitiateTires')

  end subroutine InitiateTires

end module InitiateTireTypeModule
