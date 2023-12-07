!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file motionTypeModule.f90
!> @brief Prescribed motion data container.

!!==============================================================================
!> @brief Module with data types representing prescribed motion objects.
!>
!> @details The module also contains subroutines for accessing or manipulating
!> the prescribed motion data, and for reading them from the solver input file.

module MotionTypeModule

  use KindModule                , only : dp
  use IdTypeModule              , only : IdType
  use TriadTypeModule           , only : TriadType
  use FunctionTypeModule        , only : EngineType
  use MasterSlaveJointTypeModule, only : MasterSlaveJointType

  implicit none

  integer, parameter :: prescribedDeflection_p   = 0 !< Deflection type value
  integer, parameter :: prescribedVelocity_p     = 1 !< Velocity type value
  integer, parameter :: prescribedAcceleration_p = 2 !< Acceleration type value

  !> Prescribed motion type names
  character(len=12), parameter :: motionType_p(0:2) = (/ 'deflection  ', &
       &                                                 'velocity    ', &
       &                                                 'acceleration' /)

  !> @brief Data type representing a prescribed motion object.
  type MotionType

     type(IdType) :: id   !< General identification data
     integer      :: type !< Prescribed motion type flag

     !> Points to the joint which uses this prescribed motion
     type(MasterSlaveJointType), pointer :: joint
     !> Points to the triad which uses this prescribed motion
     type(TriadType)           , pointer :: triad
     integer :: dof !< Local DOF that is prescribed in the joint or triad
     integer :: ipd !< Index into the prescribed displacement values array

     type(EngineType), pointer :: engine !< Function giving current motion value
     real(dp) :: d0 !< Time-independent prescribed motion value
     real(dp) :: d1 !< Scaling factor for the prescribed motion function

     integer           :: sDof  !< System DOF that is prescribed
     real(dp), pointer :: cc    !< Constraint coefficient (points into TTCC)
     real(dp)          :: D     !< Current motion value
     real(dp), pointer :: F     !< Current reaction force (points into RFk)
     real(dp)          :: Fprev !< Previous reaction force value
     real(dp)          :: eInp  !< Accumulated input energy from this motion

     logical :: saveVar(3) !< Flags indicating which variables should be saved

  end type MotionType

  !> @brief Standard routine for writing an object to file.
  interface WriteObject
     module procedure WriteMotionType
  end interface


contains

  !!============================================================================
  !> @brief Initializes prescribed motions with data from the solver input file.
  !>
  !> @param[in] infp File unit number for the solver input file
  !> @param[in] dispFile Name of file to read prescribed displacements from
  !> @param[in] triads Array of all triads in the model
  !> @param[in] joints Array of all joints in the model
  !> @param[in] engines All general functions in the model
  !> @param[out] motions Array of all prescribed motions in the model
  !> @param[out] err Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 11 Oct 2004

  subroutine InitiateMotions (infp,dispFile,triads,joints,engines,motions,err)

    use IdTypeModule              , only : ldesc_p, initId, ReportInputError
    use TriadTypeModule           , only : GetPtrToId
    use FunctionTypeModule        , only : GetPtrToId
    use MasterSlaveJointTypeModule, only : GetPtrToId, getJointVar
    use PrescribedMotionModule    , only : openMotionFile
    use inputUtilities            , only : iuGetNumberOfEntries
    use inputUtilities            , only : iuSetPosAtNextEntry, iuCharToInt
    use progressModule            , only : lterm
    use reportErrorModule         , only : AllocationError
    use reportErrorModule         , only : reportError, debugFileOnly_p
    use FFaCmdLineArgInterface    , only : ffa_cmdlinearg_getbool

    integer                   , intent(in)  :: infp
    character(len=*)          , intent(in)  :: dispFile
    type(TriadType)           , intent(in)  :: triads(:)
    type(MasterSlaveJointType), intent(in)  :: joints(:)
    type(EngineType)          , intent(in)  :: engines(:)
    type(MotionType)          , pointer     :: motions(:), tmp(:)
    integer                   , intent(out) :: err

    !! Local variables
    integer :: i, idIn, nMotions, nCmp, stat
    logical :: allSec, allMotion(3), allForce, allEnergy
    integer, pointer :: pdNodes(:) => null()

    !! Define the MOTION namelist
    integer            :: id, extId(10), triadId, jointId, nodeId
    integer            :: lDof, motionEngineId
    character(ldesc_p) :: extDescr
    character(len=20)  :: type
    real(dp)           :: d0, d1
    integer            :: saveVar(3)
    namelist /MOTION/ id, extId, extDescr, triadId, jointId, nodeId, &
         &            lDof, type, motionEngineId, d0, d1, saveVar

    !! --- Logic section ---

    nMotions = iuGetNumberOfEntries(infp,'&MOTION',err)
    if (err /= 0 .or. nMotions < 1) goto 900
    write(lterm,*) 'Number of &MOTION =',nMotions

    if (dispFile /= '') then
       !! Read prescribed motions from file
       call openMotionFile (dispFile,pdNodes,nCmp,err)
       if (err < 0) goto 900
    end if

    allocate(tmp(nMotions),STAT=err)
    if (err == 0) then
       call DeallocateMotions (motions)
       motions => tmp
    else
       err = AllocationError('InitiateMotions 10')
       if (associated(pdNodes)) deallocate(pdNodes)
       return
    end if

    call ffa_cmdlinearg_getbool ('allSecondaryVars',allSec)
    call ffa_cmdlinearg_getbool ('allDefVars',allMotion(1))
    call ffa_cmdlinearg_getbool ('allVelVars',allMotion(2))
    call ffa_cmdlinearg_getbool ('allAccVars',allMotion(3))
    call ffa_cmdlinearg_getbool ('allForceVars',allForce)
    call ffa_cmdlinearg_getbool ('allEnergyVars',allEnergy)

    do idIn = 1, nMotions

       call NullifyMotion (motions(idIn))
       if (.not. iuSetPosAtNextEntry(infp,'&MOTION')) then
          err = err - 1
          call ReportInputError ('MOTION',idIn)
          cycle
       end if

       !! Default values
       id=0; extId=0; extDescr=''
       triadId=0; jointId=0; nodeId=0; lDof=0
       type=''; motionEngineId=0; d0=0.0_dp; d1=0.0_dp; saveVar=0

       read(infp,nml=MOTION,iostat=stat)
       if (stat /= 0) then
          err = err - 1
          call ReportInputError ('MOTION',idIn)
          cycle
       end if

       call initId (motions(idIn)%id,id,extId,extDescr,stat)

       motions(idIn)%type = iuCharToInt(type,motionType_p,startIndex=0)
       if (motions(idIn)%type < 0) then
          err = err - 1
          call ReportInputError ('MOTION',idIn,motions(idIn)%id, &
               &                 'Invalid motion type')
          cycle
       end if

       if (triadId > 0) then

          !! Prescribed motion in a Triad DOF
          motions(idIn)%triad => GetPtrToId(triads,triadId)
          if (.not. associated(motions(idIn)%triad)) then
             err = err - 1
             call ReportInputError ('MOTION',idIn,motions(idIn)%id, &
                  &                 'Invalid triad ID')
             cycle
          else if (lDof < 1 .or. lDof > motions(idIn)%triad%nDOFs) then
             err = err - 1
             call ReportInputError ('MOTION',idIn,motions(idIn)%id, &
                  &                 'Invalid local DOF in referred triad')
             cycle
          else
             motions(idIn)%dof = lDof
          end if

          select case (motions(idIn)%type)
          case (prescribedVelocity_p)
             motions(idIn)%D = motions(idIn)%triad%urd(lDof)
          case (prescribedAcceleration_p)
             motions(idIn)%D = motions(idIn)%triad%urdd(lDof)
          end select

          if (nodeId > 0 .and. associated(pdNodes) .and. lDof <= nCmp) then
             do i = 1, size(pdNodes)
                if (pdNodes(i) == nodeId) then
                   motions(idIn)%ipd = nCmp*(i-1) + lDof
                   exit
                end if
             end do
          end if

       else if (jointId > 0) then

          !! Prescribed motion in a Joint DOF
          motions(idIn)%joint => GetPtrToId(joints,jointId)
          if (.not. associated(motions(idIn)%joint)) then
             err = err - 1
             call ReportInputError ('MOTION',idIn,motions(idIn)%id, &
                  &                 'Invalid joint ID')
             cycle
          end if

          do i = 1, motions(idIn)%joint%nJointDofs
             if (motions(idIn)%joint%jointDofs(i)%lDof == lDof) then
                motions(idIn)%dof = i
                motions(idIn)%joint%jointDofs(i)%loadIdx = -idIn
                exit
             end if
          end do
          if (motions(idIn)%dof == 0) then
             err = err - 1
             call ReportInputError ('MOTION',idIn,motions(idIn)%id, &
                  &                 'Invalid local DOF in referred joint')
             cycle
          end if

          !! Initial motion value should equal the initial joint variable value.
          !! This is the value associated with the stress free modelling config.
          motions(idIn)%D = getJointVar(motions(idIn)%joint, &
               &                        motions(idIn)%dof, &
               &                        motions(idIn)%type+1)

       else
          err = err - 1
          call ReportInputError ('MOTION',idIn,motions(idIn)%id)
          cycle
       end if

       !! Set the motion engine pointer
       if (motionEngineId > 0 .and. motions(idIn)%ipd == 0) then
          motions(idIn)%engine => GetPtrToId(engines,motionEngineId)
          if (.not. associated(motions(idIn)%engine)) then
             err = err - 1
             call ReportInputError ('MOTION',idIn,motions(idIn)%id, &
                  &                 'Invalid engine ID')
             cycle
          end if
       end if

       motions(idIn)%d0 = d0
       motions(idIn)%d1 = d1

       !! Determine which secondary triad variables will be saved. Settings in
       !! the solver input file are overruled by command-line arguments, if any.
       do i = 1, size(saveVar)
          motions(idIn)%saveVar(i) = saveVar(i) > 0 .or. allSec
       end do
       i = motions(idIn)%type+1
       if (allMotion(i)) motions(idIn)%saveVar(1) = .true. ! Motion value
       if (allForce)     motions(idIn)%saveVar(2) = .true. ! Force value
       if (allEnergy)    motions(idIn)%saveVar(3) = .true. ! Energy quantities

    end do

    if (associated(pdNodes)) deallocate(pdNodes)
900 if (err /= 0) call reportError (debugFileOnly_p,'InitiateMotions')

  end subroutine InitiateMotions


  !!============================================================================
  !> @brief Initializes pointers to reaction forces for the motion type objects.
  !>
  !> @param motions Array of all prescribed motions in the model
  !> @param[in] mpreac Matrix of Pointers to REACtion forces
  !> @param[in] RF Array of reaction force values
  !> @param[out] err Error flag
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 5 Oct 2005

  subroutine InitiateMotions2 (motions,mpreac,RF,err)

    use reportErrorModule, only : internalError

    type(MotionType), intent(inout) :: motions(:)
    integer         , intent(in)    :: mpreac(:)
    real(dp), target, intent(in)    :: RF(:)
    integer         , intent(out)   :: err

    !! Local variables
    integer :: idIn, idof, ip

    !! --- Logic section ---

    do idIn = 1, size(motions)
       idof = motions(idIn)%sdof
       if (idof > 0 .and. idof <= size(mpreac)) then
          ip = mpreac(idof)
          if (ip > 0 .and. ip <= size(RF)) then
             motions(idIn)%F => RF(ip)
          else
             err = internalError('InitiateMotions2: Invalid mpreac array')
             return
          end if
       else
          err = internalError('InitiateMotions2: Invalid global dof number')
          return
       end if
    end do

    err = 0

  end subroutine InitiateMotions2


  !!============================================================================
  !> @brief Standard routine for writing an object to file.
  !>
  !> @param[in] motion The motiontypemodule::motiontype object to write
  !> @param[in] io File unit number to write to
  !> @param[in] complexity If present, the value indicates the amount of print
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 11 Oct 2004

  subroutine WriteMotionType (motion,io,complexity)

    use IdTypeModule, only : writeId

    type(MotionType) , intent(in) :: motion
    integer          , intent(in) :: io
    integer, optional, intent(in) :: complexity

    !! --- Logic section ---

    write(io,'(A)') 'Motion','{'
    call writeId (motion%id,io)

    write(io,*) 'type        = prescribed ', motionType_p(motion%type)
    if (associated(motion%triad)) then
       call writeRefId ('triad%id    =',motion%triad%id)
    else if (associated(motion%joint)) then
       call writeRefId ('joint%id    =',motion%joint%id)
    end if
    write(io,*) 'dof         =', motion%dof
    write(io,*) 'ipd         =', motion%ipd
    write(io,*) 'sDof        =', motion%sDof

    if (associated(motion%engine)) then
       call writeRefId ('engine%id   =',motion%engine%id)
    end if
    write(io,*) 'd0, d1      =', motion%d0, motion%d1

    write(io,*) 'saveVar     =', motion%saveVar

    if (present(complexity)) then
       if (complexity >= 2) then
          if (associated(motion%cc)) then
             write(io,*) 'cc          =', motion%cc
          end if
          write(io,*) 'd           =', motion%D
          if (associated(motion%F)) then
             write(io,*) 'F           =', motion%F
          end if
          write(io,*) 'eInp        =', motion%eInp
       end if
    end if

    write(io,'(A)') '}'

  contains

    !> @brief Writes a reference to another object.
    subroutine writeRefId (text,id)
      character(len=*), intent(in) :: text
      type(IdType)    , intent(in) :: id
      write(io,*) text,id%baseId,id%userId
    end subroutine writeRefId

  end subroutine WriteMotionType


  !!============================================================================
  !> @brief Initializes the MotionType object.
  !>
  !> @param motion The motiontypemodule::motiontype object to initialize
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 11 Oct 2004

  subroutine NullifyMotion (motion)

    use IdTypeModule, only : nullifyId

    type(MotionType), intent(out) :: motion

    !! --- Logic section ---

    call nullifyId (motion%id)
    motion%type = -1

    nullify(motion%triad)
    nullify(motion%joint)
    motion%dof = 0
    motion%ipd = 0

    nullify(motion%engine)
    motion%d0 = 0.0_dp
    motion%d1 = 0.0_dp

    motion%sDof = 0
    nullify(motion%cc)

    motion%D = 0.0_dp
    nullify(motion%F)
    motion%Fprev = 0.0_dp
    motion%eInp = 0.0_dp
    motion%saveVar = .false.

  end subroutine NullifyMotion


  !!============================================================================
  !> @brief Deallocates all prescribed moption objects.
  !>
  !> @param motions Array of all motiontypemodule::motiontype objects
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Jan 2017

  subroutine DeallocateMotions (motions)

    use IdTypeModule, only : deallocateId

    type(MotionType), pointer :: motions(:)

    !! Local variables
    integer :: i

    !! --- Logic section ---

    do i = 1, size(motions)
       call deallocateId (motions(i)%id)
    end do

    deallocate(motions)
    nullify(motions)

  end subroutine DeallocateMotions

end module MotionTypeModule
