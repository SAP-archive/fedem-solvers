!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module SupElLoadTypeModule

  use SupElTypeModule   , only : SupElType, IdType
  use FunctionTypeModule, only : EngineType, dp

  implicit none

  private

  type SupElLoadType

     type(IdType) :: id  ! General identification data

     integer :: loadCase ! Index into the supElType%S vector

     type(EngineType), pointer :: engine ! Engine giving the current amplitude
     real(dp)                  :: f0, f1 ! F = f0 + f1*engineVal
     real(dp)                  :: delay  ! Engine argument delay

     type(SupElType) , pointer :: sup    ! the superelement this load acts on

     real(dp) :: F, Fprev   ! Current and previous load amplitude value
     real(dp) :: eInp       ! Accumulated input energy from this load

     logical  :: saveVar(3) ! Flags indicating which variables should be saved

  end type SupElLoadType


  interface WriteObject
     module procedure WriteLoadType
  end interface

  public :: SupElLoadType, InitiateLoads, WriteObject, DeallocateLoads


contains

  subroutine InitiateLoads (infp,sups,engines,loads,err)

    !!==========================================================================
    !! Initiates the superelement load type with data from the solver input file
    !!
    !! Programmer : Knut Morten Okstad                 date/rev : 9 Apr 2008/1.0
    !!==========================================================================

    use inputUtilities        , only : iuGetNumberOfEntries, iuSetPosAtNextEntry
    use progressModule        , only : lterm
    use SupElTypeModule       , only : GetPtrToId
    use FunctionTypeModule    , only : GetPtrToId
    use IdTypeModule          , only : ldesc_p, initId, ReportInputError
    use reportErrorModule     , only : AllocationError
    use reportErrorModule     , only : reportError, debugFileOnly_p
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getbool

    integer            , intent(in)  :: infp
    type(SupElType)    , intent(in)  :: sups(:)
    type(EngineType)   , intent(in)  :: engines(:)
    type(SupElLoadType), pointer     :: loads(:)
    integer            , intent(out) :: err

    !! Local variables
    integer :: i, idIn, nLoads, stat
    logical :: allSec, allForce, allLoad, allEnergy

    !! Define the SUPEL_LOAD namelist
    integer            :: id, extId(10), supElId, loadCase, loadEngineId
    character(ldesc_p) :: extDescr
    real(dp)           :: f0, f1, delay
    integer            :: saveVar(3)
    namelist /SUPEL_LOAD/ id, extId, extDescr, supElId, loadCase, &
         &                loadEngineId, f0, f1, delay, saveVar

    !! --- Logic section ---

    nLoads = iuGetNumberOfEntries(infp,'&SUPEL_LOAD',err)
    if (err /= 0) goto 900
    write(lterm,*) 'Number of &SUPEL_LOAD =',nLoads

    call DeallocateLoads (loads)
    allocate(loads(nLoads),STAT=stat)
    if (stat /= 0) then
       stat = AllocationError('InitiateLoads')
       return
    end if

    call ffa_cmdlinearg_getbool ('allSecondaryVars',allSec)
    call ffa_cmdlinearg_getbool ('allLoadVars',allLoad)
    call ffa_cmdlinearg_getbool ('allForceVars',allForce)
    call ffa_cmdlinearg_getbool ('allEnergyVars',allEnergy)

    do idIn = 1, nLoads

       call NullifyLoad (loads(idIn))
       if (.not. iuSetPosAtNextEntry(infp,'&SUPEL_LOAD')) then
          err = err - 1
          call ReportInputError ('SUPEL_LOAD',idIn)
          cycle
       end if

       !! Default values
       id=0; extId=0; extDescr=''
       loadCase=0; supElId=0; loadEngineId=0
       f0=0.0_dp; f1=0.0_dp; delay=0.0_dp; saveVar=0

       read(infp,nml=SUPEL_LOAD,iostat=stat)
       if (stat /= 0) then
          err = err - 1
          call ReportInputError ('SUPEL_LOAD',idIn)
          cycle
       end if

       call initId (loads(idIn)%id,id,extId,extDescr,stat)

       !! Set the pointer to the superelement the load is acting on
       loads(idIn)%sup => GetPtrToId(sups,supElId)
       if (.not. associated(loads(idIn)%sup)) then
          err = err - 1
          call ReportInputError ('SUPEL_LOAD',idIn,loads(idIn)%id, &
               &                 'Invalid superelement ID')
          cycle
       end if

       if (loadCase > 0 .and. loadCase <= size(loads(idIn)%sup%S,2)) then
          loads(idIn)%loadCase = loadCase
       else
          err = err - 1
          call ReportInputError ('SUPEL_LOAD',idIn,loads(idIn)%id, &
               &                 'Invalid load case number')
       end if

       !! Set the load engine pointer
       if (loadEngineId > 0) then
          loads(idIn)%engine => GetPtrToId(engines,loadEngineId)
          if (.not. associated(loads(idIn)%engine)) then
             err = err - 1
             call ReportInputError ('SUPEL_LOAD',idIn,loads(idIn)%id, &
                  &                 'Invalid engine ID')
          end if
       end if

       loads(idIn)%f0 = f0
       loads(idIn)%f1 = f1
       loads(idIn)%delay = delay

       !! Determine which secondary load variables will be saved. Settings in
       !! the solver input file are overruled by command-line arguments, if any.
       do i = 1, size(saveVar)
          loads(idIn)%saveVar(i) = saveVar(i) > 0 .or. allSec .or. allLoad
       end do
       if (allForce)  loads(idIn)%saveVar(1) = .true. ! Resultant force
       if (allForce)  loads(idIn)%saveVar(2) = .true. ! Resultant moment
       if (allEnergy) loads(idIn)%saveVar(3) = .true. ! Energy quantities

    end do

900 if (err < 0) call reportError (debugFileOnly_p,'InitiateLoads')

  end subroutine InitiateLoads


  subroutine WriteLoadType (load,io,complexity)

    !!==========================================================================
    !! Standard routine for writing an object to io.
    !!
    !! Programmer : Knut Morten Okstad                 date/rev : 9 Apr 2008/1.0
    !!==========================================================================

    use IdTypeModule, only : writeId

    type(SupElLoadType), intent(in) :: load
    integer            , intent(in) :: io
    integer,   optional, intent(in) :: complexity

    !! --- Logic section ---

    write(io,'(A)') 'SuperelementLoad','{'
    call writeId (load%id,io)

    write(io,*) 'loadCase    =', load%loadCase
    write(io,*) 'sup%id      =', load%sup%id%baseId, load%sup%id%userId

    if (associated(load%engine)) then
       write(io,*) 'engine%id   =', load%engine%id%baseId, load%engine%id%userId
    end if
    write(io,*) 'f0, f1      =', load%f0, load%f1
    write(io,*) 'delay       =', load%delay

    write(io,*) 'saveVar     =', load%saveVar

    if (present(complexity)) then
       if (complexity >= 2) then
          write(io,*) 'F           =', load%F
          write(io,*) 'eInp        =', load%eInp
       end if
    end if

    write(io,'(A)') '}'

  end subroutine WriteLoadType


  subroutine NullifyLoad (load)

    !!==========================================================================
    !! Initialize the SupElLoadType object.
    !!
    !! Programmer : Knut Morten Okstad                 date/rev : 9 Apr 2008/1.0
    !!==========================================================================

    use IdTypeModule, only : nullifyId

    type(SupElLoadType), intent(out) :: load

    !! --- Logic section ---

    call nullifyId (load%id)

    load%loadCase = 0

    nullify(load%engine)
    load%f0 = 0.0_dp
    load%f1 = 0.0_dp
    load%delay = 0.0_dp

    nullify(load%sup)

    load%F = 0.0_dp
    load%Fprev = 0.0_dp
    load%eInp = 0.0_dp
    load%saveVar = .false.

  end subroutine NullifyLoad


  subroutine DeallocateLoads (loads)

    !!==========================================================================
    !! Deallocate all SupElLoadType objects.
    !!
    !! Programmer : Knut Morten Okstad                date/rev : 23 Jan 2017/1.0
    !!==========================================================================

    use IdTypeModule, only : deallocateId

    type(SupElLoadType), pointer :: loads(:)

    !! Local variables
    integer :: i

    !! --- Logic section ---

    do i = 1, size(loads)
       call deallocateId (loads(i)%id)
    end do
    deallocate(loads)
    nullify(loads)

  end subroutine DeallocateLoads

end module SupElLoadTypeModule
