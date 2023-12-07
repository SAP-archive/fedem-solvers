!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file saveModule.f90
!>
!> @brief Save of time history variables to result database files.

!!==============================================================================
!> @brief Module with subroutines for saving of time history results.
!>
!> @details This module contains a set of results database file containers
!> for the FEDEM Dynamics Solver, and a collection of subroutines for writing
!> time history results data for the various object types of a FEDEM mechanism
!> to the results database.

module saveModule

  use RDBModule, only : RDBType

  implicit none

  private

  !> @brief Variable- and item group indices for response data files.
  type HeaderId
     integer :: idPos,idPosMat
     integer :: idVelVec,idAngVelVec
     integer :: idAccVec,idAngAccVec
     integer :: idForceVec,idMomentVec
     integer :: idLength,idVel,idAcc
     integer :: idAngle,idAngVel,idAngAcc
     integer :: idClen,idCG,idCGwAM,idNon
     integer :: idEstr,idEkin,idEpot,idEdmp,idEinp,idEext,idEtot
     integer :: idForce,idStiff,idSecStiff,idCoeff,idLength0,idDefl
     integer :: idYieldDefl,idElasticDefl
     integer :: idMoment,idAngStiff,idAngSecStiff,idAngCoeff,idAngle0,idAngDefl
     integer :: idYieldAngDefl,idElasticAngDefl,idAccContDis
     integer, pointer :: idGenDOF(:), idMap(:,:)
     integer, pointer :: idSpringMap(:,:), idDamperMap(:,:), idFrictionMap(:,:)
     logical :: firstCall(5)
  end type HeaderId

  type(RDBType), save :: res1DB !< File for primary response variables
  type(RDBType), save :: res2DB !< File for secondary response variables
  type(RDBType), save :: modeDB !< File for modal data
  type(RDBType), save :: ctrlDB !< File for control system data
  type(RDBType), save :: freqDB !< File for frequency domain response variables

#ifdef FT_DEBUG
  integer, parameter :: nVarCTI_p = 100 !< Number of CTI variables in debug mode
#else
  integer, parameter :: nVarCTI_p = 49  !< Number of CTI variables
#endif
  integer, parameter :: max_UDE_p = 100 !< Max. user-defined element variables

  public :: res1DB, res2DB, modeDB, ctrlDB, freqDB
  public :: initiateSaveModule, writeSolverHeaders
  public :: writeSolverDB1, writeSolverDB2, writeSolverDB5
  public :: flushResultFiles, closeResultFiles


contains

  !!============================================================================
  !> @brief Writes an item group reference to the file header.
  !>
  !> @param[in] iFile File unit number to write to
  !> @param[in] itemGroup Item group ID
  !>
  !> @callergraph

  subroutine writeItgRef (iFile,itemGroup)

    use RDBModule, only : idatd

    integer, intent(in) :: iFile, itemGroup

    !! --- Logic section ---

    select case (itemGroup)
    case ( 0: 9); write(iFile,601,advance='NO') itemGroup
    case (10:99); write(iFile,602,advance='NO') itemGroup
    case default; write(iFile,603,advance='NO') itemGroup
    end select
    if (iFile == idatd) write(iFile,"(']')",advance='NO')
601 format('[',i1)
602 format('[',i2)
603 format('[',i3)

  end subroutine writeItgRef


  !!============================================================================
  !> @brief Initializes the solver result database file.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 4 Feb 2008

  subroutine initiateSaveModule ()

    use RDBModule, only : nullifyRDB

    !! --- Logic section ---

    call nullifyRDB (res1DB)
    call nullifyRDB (res2DB)
    call nullifyRDB (modeDB)
    call nullifyRDB (ctrlDB)
    call nullifyRDB (freqDB)

  end subroutine initiateSaveModule


  !!============================================================================
  !> @brief Administers writing of results database headers.
  !>
  !> @param[in] cResponseFiles Array of file names
  !> @param[in] mech Container of all mechanism objects
  !> @param[in] control Control system data
  !> @param[in] modes Eigenmode data
  !> @param[in] bufRat Relative buffer sizes for each result file
  !> @param[in] neq Number of equations (system vector dimension)
  !> @param[in] writeAsDouble Flags for writing double vs. single precision
  !> @param[out] ierr Error flag
  !> @param[in] modelfile Model file name from which current model is generated
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date September 2000

  subroutine writeSolverHeaders (cResponseFiles,mech,control,modes,bufRat, &
       &                         neq,writeAsDouble,ierr,modelFile)

    use kindModule         , only : dp, i8
    use MechanismTypeModule, only : MechanismType
    use ControlTypeModule  , only : ControlType
    use ControlTypeModule  , only : hasControlElements, writeCtrlSysHeader
    use ModesTypeModule    , only : ModesType, writeModesHeader
    use RDBModule          , only : openHeaderFiles, openRDBfile, closeRDBfile
    use RDBModule          , only : writeTimeStepHeader
    use RDBModule          , only : writeVarDef, writeItgDef
    use RDBModule          , only : ivard, iitem, idatd, nbi_p, nbs_p, nbd_p
    use profilerModule     , only : sav_p, startTimer, stopTimer
    use reportErrorModule  , only : allocationError
    use reportErrorModule  , only : reportError, debugFileOnly_p

    character(len=*)   , intent(in)  :: cResponseFiles(:)
    type(MechanismType), intent(in)  :: mech
    type(ControlType)  , intent(in)  :: control
    type(ModesType)    , intent(in)  :: modes
    real(dp)           , intent(in)  :: bufRat(:)
    integer            , intent(in)  :: neq
    logical            , intent(in)  :: writeAsDouble(:)
    integer            , intent(out) :: ierr
    character(len=*), optional, intent(in) :: modelFile

    !! Local variables
    type(HeaderId) :: header
    integer        :: i, i2, maxGenDOFs, nSupGenDOF, nMap, nBits

    !! --- Logic section ---

    call startTimer (sav_p)

    ! --- Initialize response header variables

    maxGenDOFs = 0
    nSupGenDOF = 0
    do i = 1, size(mech%sups)
       if (associated(mech%sups(i)%genDOFs) .and. mech%sups(i)%saveVar(2)) then
          nSupGenDOF = nSupGenDOF + 1
          maxGenDOFs = max(maxGenDOFs,mech%sups(i)%genDOFs%nDOFs)
       end if
    end do
    nMap = max(min(63,6*size(mech%joints)+3*size(mech%cElems)), &
         &     min(14,2*size(mech%triads)),nSupGenDOF)
    allocate(header%idGenDOF(maxGenDOFs), header%idMap(2,nMap), &
         &   header%idSpringMap(2,min(255,size(mech%baseSprings))), &
         &   header%idDamperMap(2,min(255,size(mech%baseDampers))), &
         &   header%idFrictionMap(2,min(15,size(mech%frictions))), stat=ierr)
    if (ierr /= 0) then
       ierr = AllocationError('writeSolverHeaders')
       return
    end if

    call initHeader ()

    ierr = size(cResponseFiles) - 1
    if (ierr < 0) goto 500

    if (ierr == 0) then
       i2 = 1 ! Write primary and secondary variables to a single file
    else
       i2 = 2 ! Check if we are writing a separate primary variable file
       do i = 1, size(mech%triads)
          if (mech%triads(i)%savePos) goto 10
       end do
       do i = 1, size(mech%sups)
          if (mech%sups(i)%savePos) goto 10
       end do
    end if
    goto 100


    ! --- Open the temporary header files

 10 continue
    call openHeaderFiles ('fedem_solver',ierr, &
         &                'response data base file',modelFile)
    if (ierr < 0) goto 500

    if (writeAsDouble(1)) then
       nBits = 64
    else
       nBits = 32
    end if

    ! --- Write response data definitions for the primary variables

    call writeTimeStepHeader (res1DB)

    do i = 1, size(mech%triads)
       call writeTriadHeader1 (res1DB,header,mech%triads(i),nBits)
    end do

    header%idMap = 0
    do i = 1, size(mech%sups)
       call writeSupElHeader1 (res1DB,header,mech%sups(i),nBits,ierr)
       if (ierr < 0) goto 500
    end do

    ! --- Set buffer size for the primary variable file

    res1DB%nBstep  = int(res1DB%nBytes)
    res1DB%bufSize = int(res1DB%nBytes*bufRat(1))
    res1DB%nBytes  = 0_i8

    ! --- Open the response database file and copy the file header to it

    call openRDBfile (res1DB,ierr,cResponseFiles(1),'#FEDEM response data')
    if (ierr < 0) goto 500


    ! --- Open the temporary header files

100 continue
    header%idPos    = 0
    header%idPosMat = 0
    call openHeaderFiles ('fedem_solver',ierr, &
         &                'response data base file',modelFile)
    if (ierr < 0) goto 500

    if (writeAsDouble(i2)) then
       nBits = 64
    else
       nBits = 32
    end if

    ! --- Write response data definitions for the secondary variables

    call writeTimeStepHeader (res2DB)

    if (associated(mech%env%waveFunc) .and. mech%saveVar(5)) then
       i = 0
       res2DB%nBytes = res2DB%nBytes + nbs_p
       call writeVarDef (res2DB,i,idatd,'Wave elevation at origin','LENGTH')
       write(idatd,"()")
    end if

    header%idMap = 0
    do i = 1, size(mech%triads)
       call writeTriadHeader2 (res2DB,header,mech%triads(i),nBits,i2==1,ierr)
       if (ierr < 0) goto 500
    end do

    header%idMap = 0
    header%idGenDOF = 0
    do i = 1, size(mech%sups)
       call writeSupElHeader2 (res2DB,header,mech%sups(i),nBits,i2==1,ierr)
       if (ierr < 0) goto 500
    end do
    do i = 1, size(mech%elms)
       call writeUserdefElHeader (res2DB,header,mech%elms(i),nBits)
    end do

    header%idSpringMap = 0
    do i = 1, size(mech%axialSprings)
       call writeSpringElementHeader (res2DB,header,mech%axialSprings(i), &
            &                         nBits,ierr)
       if (ierr < 0) goto 500
    end do

    header%idDamperMap = 0
    do i = 1, size(mech%dampers)
       if (mech%dampers(i)%dmp%dof > 0)  cycle ! Only axial dampers written here
       if (mech%dampers(i)%samElNum < 1) cycle ! Ignore unused joint dampers
       call writeAxialDamperHeader (res2DB,header,mech%dampers(i),nBits,ierr)
       if (ierr < 0) goto 500
    end do

    do i = 1, size(mech%forces)
       call writeForceHeader (res2DB,header,mech%forces(i),nBits)
    end do

    header%idMap = 0
    header%idFrictionMap = 0
    do i = 1, size(mech%joints)
       call writeJointHeader (res2DB,header,mech%joints(i), &
            &                 mech%forces,mech%motions,nBits,ierr)
       if (ierr < 0) goto 500
    end do

    do i = 1, size(mech%cElems)
       call writeContactHeader (res2DB,header,mech%cElems(i),nBits,ierr)
       if (ierr < 0) goto 500
    end do

    header%idMap = 0
    do i = 1, size(mech%tires)
       call writeTireHeader (res2DB,header,mech%tires(i),nBits,ierr)
       if (ierr < 0) goto 500
    end do

    do i = 1, size(mech%engines)
       call writeEngineHeader (res2DB,header,mech%engines(i))
    end do

    call writeControlHeader (res2DB,header,control)

    if (associated(mech%turbine)) then
       call writeWindTurbineHeader (res2DB,mech%turbine)
    end if

    call writeMechanismHeader (res2DB,header,mech,neq)

    deallocate(header%idGenDOF)
    deallocate(header%idSpringMap,header%idDamperMap,header%idFrictionMap)
    nullify(header%idGenDOF)
    nullify(header%idSpringMap)
    nullify(header%idDamperMap)
    nullify(header%idFrictionMap)

    ! --- Set buffer size for the secondary variable file

    res2DB%nBstep  = int(res2DB%nBytes)
    res2DB%bufSize = int(res2DB%nBytes*bufRat(i2))
    res2DB%nBytes  = 0_i8

    ! --- Open the response database file and copy the file header to it

    if (res2DB%nBstep > nbi_p + nbd_p) then
       call openRDBfile (res2DB,ierr,cResponseFiles(i2),'#FEDEM response data')
       if (ierr < 0) goto 500
    else
       !! No secondary variable output
       close(ivard)
       close(iitem)
       close(idatd)
    end if
    if (size(cResponseFiles) < 3) goto 500
    if (modes%nModes < 1) goto 200


    ! --- Open the temporary header files

    call openHeaderFiles ('fedem_solver',ierr, &
         &                'modes data base file',modelFile)
    if (ierr < 0) goto 500

    if (writeAsDouble(1)) then
       nBits = 64
    else
       nBits = 32
    end if

    ! --- Write modal data definitions

    call writeTimeStepHeader (modeDB)
    call writeModesHeader (modes,mech%id,mech%triads,mech%sups,modeDB,nBits)

    ! --- Set buffer size for the modes file

    modeDB%nBstep  = int(modeDB%nBytes)
    modeDB%bufSize = int(modeDB%nBytes*bufRat(3))
    modeDB%nBytes  = 0_i8

    ! --- Open the modes database file and copy the file header to it

    call openRDBfile (modeDB,ierr,cResponseFiles(3),'#FEDEM modal data')
200 if (ierr < 0 .or. size(cResponseFiles) < 4) goto 500
    if (.not.(hasControlElements(control) .and. control%saveVar)) goto 300


    ! --- Open the temporary header files

    call openHeaderFiles ('fedem_solver',ierr, &
         &                'control system data base file',modelFile)
    if (ierr < 0) goto 500

    ! --- Write control system data definitions

    call writeTimeStepHeader (ctrlDB)
    call writeCtrlSysHeader (control,mech%id,ctrlDB)

    ! --- Set buffer size for the control system file

    ctrlDB%nBstep  = int(ctrlDB%nBytes)
    ctrlDB%bufSize = int(ctrlDB%nBytes*bufRat(4))
    ctrlDB%nBytes  = 0_i8

    ! --- Open the control system database file and copy the file header to it

    call openRDBfile (ctrlDB,ierr,cResponseFiles(4),'#FEDEM control data')
300 if (ierr < 0 .or. size(cResponseFiles) < 5) goto 500
    if (len_trim(cResponseFiles(5)) < 1) goto 500 ! No frequency domain output


    ! --- Open the temporary header files

    call openHeaderFiles ('fedem_solver',ierr, &
         &                'frequency response data base file',modelFile)
    if (ierr < 0) goto 500

    ! --- Write response data definitions for the frequency response variables

    nBits = 32
    call initHeader ()

    call writeTimeStepHeader (freqDB)

    do i = 1, size(mech%triads)
       call writeTriadHeader1 (freqDB,header,mech%triads(i),nBits)
       call writeTriadHeader2 (freqDB,header,mech%triads(i),nBits,ierr=ierr)
       if (ierr < 0) goto 500
    end do

    header%idMap = 0
    do i = 1, size(mech%sups)
       call writeSupElHeader1 (freqDB,header,mech%sups(i),nBits,ierr)
       if (ierr < 0) goto 500
       call writeSupElHeader2 (freqDB,header,mech%sups(i),nBits,ierr=ierr)
       if (ierr < 0) goto 500
    end do

    ! --- Set buffer size for the frequency response file

    freqDB%nBstep  = int(freqDB%nBytes)
    freqDB%bufSize = int(freqDB%nBytes*bufRat(5))
    freqDB%nBytes  = 0_i8

    ! --- Open the response database file and copy the file header to it

    call openRDBfile (freqDB,ierr,cResponseFiles(5),'#FEDEM response data')


500 if (ierr < 0) call reportError (debugFileOnly_p,'writeSolverHeaders')
    deallocate(header%idMap)
    call stopTimer (sav_p)

  contains

    !> @brief Initializes all response header variables.
    subroutine initHeader
      header%idPos       = 0
      header%idPosMat    = 0
      header%idVelVec    = 0
      header%idAngVelVec = 0
      header%idAccVec    = 0
      header%idAngAccVec = 0
      header%idForceVec  = 0
      header%idMomentVec = 0
      header%idLength    = 0
      header%idVel       = 0
      header%idAcc       = 0
      header%idAngle     = 0
      header%idAngVel    = 0
      header%idAngAcc    = 0
      header%idCLen      = 0
      header%idCG        = 0
      header%idCGwAM     = 0
      header%idNon       = 0
      header%idEstr      = 0
      header%idEkin      = 0
      header%idEpot      = 0
      header%idEdmp      = 0
      header%idEinp      = 0
      header%idEext      = 0
      header%idEtot      = 0
      header%idForce     = 0
      header%idStiff     = 0
      header%idSecStiff  = 0
      header%idLength0   = 0
      header%idDefl      = 0
      header%idYieldDefl = 0
      header%idElasticDefl = 0
      header%idCoeff     = 0
      header%idMoment    = 0
      header%idAngStiff  = 0
      header%idAngSecStiff = 0
      header%idAngle0    = 0
      header%idAngDefl   = 0
      header%idYieldAngDefl = 0
      header%idElasticAngDefl = 0
      header%idAngCoeff  = 0
      header%idAccContDis = 0
      if (associated(header%idGenDOF))      header%idGenDOF    = 0
      if (associated(header%idMap))         header%idMap       = 0
      if (associated(header%idSpringMap))   header%idSpringMap = 0
      if (associated(header%idDamperMap))   header%idDamperMap = 0
      if (associated(header%idFrictionMap)) header%idFrictionMap = 0
      header%firstCall = .true.
    end subroutine initHeader

  end subroutine writeSolverHeaders


  !!============================================================================
  !> @brief Flushes the results database file(s) to disk.
  !>
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date January 2003

  subroutine flushResultFiles (ierr)

    use RDBModule        , only : flushRDBfile
    use profilerModule   , only : sav_p, startTimer, stopTimer
    use reportErrorModule, only : reportError, debugFileOnly_p

    integer, intent(out) :: ierr

    !! --- Logic section ---

    ierr = 0

    call startTimer (sav_p)
    call flushRDBfile (res1DB,ierr)
    call flushRDBfile (res2DB,ierr)
    call flushRDBfile (modeDB,ierr)
    call flushRDBfile (ctrlDB,ierr)
    call flushRDBfile (freqDB,ierr)
    call stopTimer (sav_p)

    if (ierr < 0) call reportError (debugFileOnly_p,'flushResultFiles')

  end subroutine flushResultFiles


  !!============================================================================
  !> @brief Closes the results database file(s).
  !>
  !> @param[out] haveData Flags telling whether each file have any data or not
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date September 2000

  subroutine closeResultFiles (haveData,ierr)

    use RDBModule        , only : closeRDBfile
    use profilerModule   , only : sav_p, startTimer, stopTimer
    use reportErrorModule, only : reportError, debugFileOnly_p

    logical, intent(out) :: haveData(5)
    integer, intent(out) :: ierr

    !! --- Logic section ---

    ierr = 0

    call startTimer (sav_p)
    call closeRDBfile (res1DB,ierr,haveData=haveData(1))
    call closeRDBfile (res2DB,ierr,haveData=haveData(2))
    call closeRDBfile (modeDB,ierr,haveData=haveData(3))
    call closeRDBfile (ctrlDB,ierr,haveData=haveData(4))
    call closeRDBfile (freqDB,ierr,haveData=haveData(5))
    call stopTimer (sav_p)

    if (ierr < 0) call reportError (debugFileOnly_p,'closeResultFiles')

  end subroutine closeResultFiles


  !!============================================================================
  !> @brief Administers writing of primary response results to database file.
  !>
  !> @param[in] sys System level data
  !> @param[in] triads Array of all triad objects
  !> @param[in] sups Array of all superelement objects
  !> @param[in] writeAsDouble Flag for writing double vs. single precision
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date September 2000

  subroutine writeSolverDB1 (sys,triads,sups,writeAsDouble,ierr)

    use SystemTypeModule , only : SystemType
    use TriadTypeModule  , only : TriadType
    use SupElTypeModule  , only : SupElType
    use RDBModule        , only : writeTimeStepDB, checkStepSize
    use reportErrorModule, only : reportError, debugFileOnly_p

    type(SystemType), intent(in)  :: sys
    type(TriadType) , intent(in)  :: triads(:)
    type(SupElType) , intent(in)  :: sups(:)
    logical         , intent(in)  :: writeAsDouble
    integer         , intent(out) :: ierr

    !! Local variables
    integer :: i, nBytesPerStep = 0

    !! --- Logic section ---

    ierr = 0
    if (res1DB%fileNumb < 0) return

    if (nBytesPerStep == 0) nBytesPerStep = -int(res1DB%nBytes)

    call writeTimeStepDB (res1DB,sys%nStep,sys%time,ierr)
    if (ierr < 0) goto 500

    do i = 1, size(triads)
       call writeTriadDB1 (res1DB,triads(i),writeAsDouble,ierr)
       if (ierr < 0) goto 500
    end do

    do i = 1, size(sups)
       call writeSupelDB1 (res1DB,sups(i),writeAsDouble,ierr)
       if (ierr < 0) goto 500
    end do

    if (nBytesPerStep < 1) then
       !! Check consistency between header and the actual data
       nBytesPerStep = nBytesPerStep + int(res1DB%nBytes)
       call checkStepSize (res1DB,nBytesPerStep,ierr)
       if (ierr < 0) goto 500
    end if

    return

500 call reportError (debugFileOnly_p,'writeSolverDB1')

  end subroutine writeSolverDB1


  !!============================================================================
  !> @brief Administers writing of secondary response results to database file.
  !>
  !> @param[in] sys System level data
  !> @param[in] mech Container of all mechanism objects
  !> @param[in] control Control system data
  !> @param[in] scaleTo Scaling factors for tire model data
  !> @param[in] eta Current wave elevation at global origin
  !> @param[in] ws Current wind speed at user-defined points
  !> @param[in] writeAsDouble Flag for writing double vs. single precision
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date September 2000

  subroutine writeSolverDB2 (sys,mech,control,scaleTo,eta,ws,writeAsDouble,ierr)

    use SystemTypeModule   , only : SystemType
    use MechanismTypeModule, only : MechanismType
    use ControlTypeModule  , only : ControlType
    use TireTypeModule     , only : TireType, ScaleToType
    use KindModule         , only : dp, epsDiv0_p
    use RDBModule          , only : writeTimeStepDB, writeRDB, checkStepSize
    use reportErrorModule  , only : reportError, debugFileOnly_p

    type(SystemType)   , intent(in)  :: sys
    type(MechanismType), intent(in)  :: mech
    type(ControlType)  , intent(in)  :: control
    type(ScaleToType)  , intent(in)  :: scaleTo
    real(dp)           , intent(in)  :: eta, ws(:,:)
    logical            , intent(in)  :: writeAsDouble
    integer            , intent(out) :: ierr

    !! Local variables
    integer  :: i, j, nBytesPerStep = 0
    real(dp) :: sysMass(3), sysCG(3), sysVb, sysCB(3)

    !! --- Logic section ---

    ierr = 0
    if (res2DB%fileNumb < 0) return

    if (nBytesPerStep == 0) nBytesPerStep = -int(res2DB%nBytes)

    call writeTimeStepDB (res2DB,sys%nStep,sys%time,ierr)
    if (ierr < 0) goto 500

    if (associated(mech%env%waveFunc) .and. mech%saveVar(5)) then
       call writeRDB (res2DB,eta,ierr)
       if (ierr < 0) goto 500
    end if

    do i = 1, size(mech%triads)
       call writeTriadDB2 (res2DB,mech%triads(i),writeAsDouble,ierr)
       if (ierr < 0) goto 500
       if (mech%triads(i)%saveVar(8)) then
          call writeTriadForcesDB (res2DB,mech%triads(i),ierr)
          if (ierr < 0) goto 500
       end if
    end do

    sysMass = 0.0_dp
    sysCG = 0.0_dp
    sysVb = 0.0_dp
    sysCB = 0.0_dp
    do i = 1, size(mech%sups)
       call writeSupelDB2 (res2DB,mech%sups(i),mech%masses, &
            &              sysMass,sysCG,sysVb,sysCB,writeAsDouble,ierr)
       if (ierr < 0) goto 500
    end do

    do i = 1, size(mech%elms)
       call writeUserdefElDB (res2DB,mech%elms(i),writeAsDouble,ierr)
       if (ierr < 0) goto 500
    end do

    do i = 1, size(mech%axialSprings)
       do j = 1, size(mech%axialSprings(i)%spr)
          if (associated(mech%axialSprings(i)%spr(j)%p)) then
             call writeSpringDB (res2DB,mech%axialSprings(i)%spr(j)%p, &
                  &              writeAsDouble,ierr)
             if (ierr < 0) goto 500
          end if
       end do
    end do

    do i = 1, size(mech%dampers)
       if (mech%dampers(i)%dmp%dof > 0)  cycle ! Only axial dampers written here
       if (mech%dampers(i)%samElNum < 1) cycle ! Ignore unused joint dampers
       call writeDamperDB (res2DB,mech%dampers(i)%dmp,writeAsDouble,ierr)
       if (ierr < 0) goto 500
    end do

    do i = 1, size(mech%forces)
       call writeForceDB (res2DB,mech%forces(i),writeAsDouble,ierr)
       if (ierr < 0) goto 500
    end do

    do i = 1, size(mech%joints)
       call writeJointDB (res2DB,mech%joints(i),mech%forces,mech%motions, &
            &             writeAsDouble,ierr)
       if (ierr < 0) goto 500
    end do

    do i = 1, size(mech%cElems)
       call writeContactDB (res2DB,mech%cElems(i),writeAsDouble,ierr)
       if (ierr < 0) goto 500
    end do

    do i = 1, size(mech%tires)
       call writeTireDB (res2DB,mech%tires(i),scaleTo,sysMass,sysCG, &
            &            writeAsDouble,ierr)
       if (ierr < 0) goto 500
    end do

    do i = 1, size(mech%engines)
       call writeEngineDB (res2DB,mech%engines(i),ierr)
       if (ierr < 0) goto 500
    end do

    call writeControlDB (res2DB,control,ierr)
    if (ierr < 0) goto 500

    if (associated(mech%turbine)) then
       call writeRDB (res2DB,ws,ierr)
       if (ierr < 0) goto 500
    end if

    do i = 1, 3
       if (sysMass(i) > epsDiv0_p) sysCG(i) = sysCG(i) / sysMass(i)
    end do
    if (sysVb > epsDiv0_p) sysCB = sysCB / sysVb
    call writeMechanismDB (res2DB,mech,sys,sysCG,sysCB,sysVb,ierr)
    if (ierr < 0) goto 500

    if (nBytesPerStep < 1) then
       !! Check consistency between header and the actual data
       nBytesPerStep = nBytesPerStep + int(res2DB%nBytes)
       call checkStepSize (res2DB,nBytesPerStep,ierr)
       if (ierr < 0) goto 500
    end if

    return

500 call reportError (debugFileOnly_p,'writeSolverDB2')

  end subroutine writeSolverDB2


  !!============================================================================
  !> @brief Administers writing of frequency response results to database file.
  !>
  !> @param[in] sys System level data
  !> @param[in] triads Array of all triad objects
  !> @param[in] sups Array of all superelement objects
  !> @param[in] masses All point masses in the model
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date June 2019

  subroutine writeSolverDB5 (sys,triads,sups,masses,ierr)

    use SystemTypeModule , only : SystemType, dp
    use TriadTypeModule  , only : TriadType
    use SupElTypeModule  , only : SupElType
    use MassTypeModule   , only : MassType
    use RDBModule        , only : writeTimeStepDB, checkStepSize
    use reportErrorModule, only : reportError, debugFileOnly_p

    type(SystemType), intent(in)  :: sys
    type(TriadType) , intent(in)  :: triads(:)
    type(SupElType) , intent(in)  :: sups(:)
    type(MassType)  , intent(in)  :: masses(:)
    integer         , intent(out) :: ierr

    !! Local variables
    integer  :: i, nBytesPerStep = 0
    real(dp) :: sysMass(3), sysCG(3), sysVb, sysCB(3)

    !! --- Logic section ---

    ierr = 0
    if (freqDB%fileNumb < 0) return

    if (nBytesPerStep == 0) nBytesPerStep = -int(freqDB%nBytes)

    call writeTimeStepDB (freqDB,sys%nStep,sys%time,ierr)
    if (ierr < 0) goto 500

    do i = 1, size(triads)
       call writeTriadDB1 (freqDB,triads(i),.false.,ierr)
       if (ierr < 0) goto 500
       call writeTriadDB2 (freqDB,triads(i),.false.,ierr)
       if (ierr < 0) goto 500
       if (triads(i)%saveVar(8)) then
          call writeTriadForcesDB (freqDB,triads(i),ierr)
          if (ierr < 0) goto 500
       end if
    end do

    sysMass = 0.0_dp
    sysCG = 0.0_dp
    sysVb = 0.0_dp
    sysCB = 0.0_dp
    do i = 1, size(sups)
       call writeSupelDB1 (freqDB,sups(i),.false.,ierr)
       if (ierr < 0) goto 500
       call writeSupelDB2 (freqDB,sups(i),masses, &
            &              sysMass,sysCG,sysVb,sysCB,.false.,ierr)
       if (ierr < 0) goto 500
    end do

    if (nBytesPerStep < 1) then
       !! Check consistency between header and the actual data
       nBytesPerStep = nBytesPerStep + int(freqDB%nBytes)
       call checkStepSize (freqDB,nBytesPerStep,ierr)
       if (ierr < 0) goto 500
    end if

    return

500 call reportError (debugFileOnly_p,'writeSolverDB5')

  end subroutine writeSolverDB5


  !!============================================================================
  !> @brief Helper function inserting an integer pair into a map.
  !>
  !> @param[in] iVal The first (key value) integer of the pair to insert
  !> @param iMap The integer map array
  !> @param nId Current size of the map
  !> @param[out] id The second value of the inserted integer pair
  !>
  !> @details The non-zero integer value @a iVal is mapped into an integer
  !> value @a id in the range [1,nId], where @a nId is the number of distinct
  !> values. Returns .true. unless the value is already in the mapping.
  !> If the provided array is too small, the variable @a id is set to -1.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 30 Apr 2002

  function insertIdMap (iVal,iMap,nId,id)

    use reportErrorModule, only : internalError

    integer, intent(in)    :: iVal
    integer, intent(inout) :: iMap(:,:), nId
    integer, intent(out)   :: id

    logical :: insertIdMap
    integer :: i

    do i = 1, size(iMap,2)
       if (iMap(1,i) == iVal) then
          !! iVal is already in the mapping
          insertIdMap = .false.
          id = iMap(2,i)
          return
       else if (iMap(1,i) == 0) then
          !! Append iVal to the mapping
          insertIdMap = .true.
          nId = nId + 1
          id  = nId
          iMap(1,i) = iVal
          iMap(2,i) = id
          return
       end if
    end do

    insertIdMap = .false.
    id = internalError('insertIdMap: Array iMap is too small')

  end function insertIdMap


  !!============================================================================
  !> @brief Writes an item group reference for a joint or triad DOF.
  !>
  !> @param[in] itemGroup Item group ID
  !> @param[in] jDof Local DOF index
  !> @param[in] cIG Item group name tag
  !>
  !> @details The beginning of the item group definition for a DOF configuration
  !> is written to the item group definition section of the results file header,
  !> and a reference to it is written to the associated data section.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 30 Apr 2002

  subroutine writeDofIG (itemGroup,jDof,cIG)

    use RDBModule, only : iitem, idatd

    integer         , intent(in)           :: itemGroup
    integer         , intent(in), optional :: jDof
    character(len=*), intent(in), optional :: cIG

    if (present(jDof) .and. present(cIG)) then
       !! This variable configuration is new so write the item group definition
       call writeItgRef (iitem,itemGroup)
       if (jDof <= 3) then
          write(iitem,600,advance='NO') 'T'//char(ichar('x')+jDof-1),cIG
       else if (jDof <= 6) then
          write(iitem,600,advance='NO') 'R'//char(ichar('x')+jDof-4),cIG
       else
          write(iitem,600,advance='NO') 'Slider',cIG
       end if
    end if

    !! Write reference to the (new or existing) item group in the data section
    call writeItgRef (idatd,itemGroup)

600 format(';"',a,1x,a,' variables";')

  end subroutine writeDofIG


  !!============================================================================
  !> @brief Writes primary results database file header for a superelement.
  !>
  !> @param respDB File for primary response variables
  !> @param header Variable- and item group index container
  !> @param[in] supel The superelement to write file header for
  !> @param[in] nBits Size of a real variables (32 or 64 bit)
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph

  subroutine writeSupElHeader1 (respDB,header,supel,nBits,ierr)

    use SupElTypeModule  , only : SupElType, IsBeam, GetSupElId
    use IdTypeModule     , only : writeIdHeader
    use RDBModule        , only : ivard, idatd, writeVarDef
    use reportErrorModule, only : debugFileOnly_p, reportError

    type(RDBType)  , intent(inout) :: respDB
    type(HeaderId) , intent(inout) :: header
    type(SupElType), intent(in)    :: supel
    integer        , intent(in)    :: nBits
    integer        , intent(out)   :: ierr

    integer :: idDis, nDOFs

    ierr = 0
    if (.not. supel%savePos) return ! Output is switched off

    if (IsBeam(supel)) then
       call writeIdHeader ('Beam',supel%id,idatd)
    else
       call writeIdHeader ('Part',supel%id,idatd)
    end if

    respDB%nBytes = respDB%nBytes + 12*nBits/8
    call writeVarDef (respDB,header%idPosMat,idatd, &
         &            'Position matrix','NONE','TMAT34', &
         &            '"x","y","z"','"i1","i2","i3","position"',nBits)

    if (associated(supel%genDOFs)) then

       nDOFs = supel%genDOFs%nDOFs
       respDB%nBytes = respDB%nBytes + nDOFs*nBits/8
       if (insertIdMap(nDOFs,header%idMap,respDB%nVar,idDis)) then
          select case (idDis)
          case  (0: 9); write(ivard,601,advance='NO') idDis
          case (10:99); write(ivard,602,advance='NO') idDis
          case default; write(ivard,603,advance='NO') idDis
          end select
          select case (nDOFs)
          case  (0: 9); write(ivard,611) nBits,nDOFs
          case (10:99); write(ivard,612) nBits,nDOFs
          case default; write(ivard,613) nBits,nDOFs
          end select
       else if (idDis < 0) then
          ierr = idDis
          call reportError (debugFileOnly_p,'writeSupElHeader1: '// &
               &            getSupElId(supel))
          return
       end if

       select case (idDis)
       case  (0: 9); write(idatd,601,advance='NO') idDis
       case (10:99); write(idatd,602,advance='NO') idDis
       case default; write(idatd,603,advance='NO') idDis
       end select

       write(idatd,"('>}')")
    else
       write(idatd,"('}')")
    end if

601 format('<',i1)
602 format('<',i2)
603 format('<',i3)
611 format(';"Generalized displacement";LENGTH;FLOAT;',i2,';VECTOR;(',i1,')>')
612 format(';"Generalized displacement";LENGTH;FLOAT;',i2,';VECTOR;(',i2,')>')
613 format(';"Generalized displacement";LENGTH;FLOAT;',i2,';VECTOR;(',i3,')>')

  end subroutine writeSupElHeader1


  !!============================================================================
  !> @brief Writes secondary results database file header for a superelement.
  !>
  !> @param respDB File for secondary response variables
  !> @param header Variable- and item group index container
  !> @param[in] supel The superelement to write file header for
  !> @param[in] nBits Size of a real variables (32 or 64 bit)
  !> @param[in] includePosMat If .true., include the position matrix also
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph

  subroutine writeSupElHeader2 (respDB,header,supel,nBits,includePosMat,ierr)

    use SupElTypeModule  , only : SupElType, IsBeam, GetSupElId
    use IdTypeModule     , only : writeIdHeader
    use RDBModule        , only : iitem, idatd, nbs_p, writeVarDef, writeItgDef
    use reportErrorModule, only : debugFileOnly_p, reportError

    type(RDBType)   , intent(inout) :: respDB
    type(HeaderId)  , intent(inout) :: header
    type(SupElType) , intent(in)    :: supel
    integer         , intent(in)    :: nBits
    logical,optional, intent(in)    :: includePosMat
    integer,optional, intent(out)   :: ierr

    logical :: newGroup, lBeam
    integer :: i, idGen, iFile, nDOFs
    integer, save :: idHD(-1:9), idSF(4)
    integer, save :: idAmp = 0, idCG = 0, idSen = 0

    if (present(ierr)) ierr = 0

    if (.not. any(supel%saveVar)) then
       return ! No secondary results for this superelement
    else if (IsBeam(supel)) then
       call writeIdHeader ('Beam',supel%id,idatd)
       lBeam = .true.
    else
       call writeIdHeader ('Part',supel%id,idatd)
       lBeam = .false.
    end if

    if (header%firstCall(1)) then
       idHD  = 0
       idSF  = 0
       idAmp = 0
       idCG  = 0
       idSen = 0
    end if

    if (present(includePosMat)) then
       if (includePosMat) then
          !! For creation of response_pos.frs only
          respDB%nBytes = respDB%nBytes + 12*nbs_p
          call writeVarDef (respDB,header%idPosMat,idatd, &
               &            'Position matrix','NONE','TMAT34', &
               &            '"x","y","z"','"i1","i2","i3","position"')
       end if
    end if

    if (supel%saveVar(2) .or. supel%saveVar(4)) then
       if (associated(supel%scoord)) then

          respDB%nBytes = respDB%nBytes + nbs_p
          call writeVarDef (respDB,header%idCLen,idatd,'Curve length','LENGTH')

       end if
    end if
    if (supel%saveVar(1)) then

       respDB%nBytes = respDB%nBytes + 3*nbs_p
       call writeVarDef (respDB,header%idCG,idatd, &
            &            'Center of Gravity', &
            &            'LENGTH','VEC3','"x_cg","y_cg","z_cg"')

       if (supel%addedMass) then
          respDB%nBytes = respDB%nBytes + 3*nbs_p
          call writeVarDef (respDB,header%idCGwAM,idatd, &
               &            'Center of Gravity including additional mass', &
               &            'LENGTH','VEC3','"x_cg","y_cg","z_cg"')
       end if

       if ( supel%triads(1)%p%id%baseId < 0 .and. &
            supel%triads(1)%p%nDOFs > 0 ) then

          !! Results at the center of gravity dummy triad
          respDB%nBytes = respDB%nBytes + 24*nBits/8
          call writeItGDef (respDB,idCG,iFile,'Center of Gravity Triad')
          call writeVarDef (respDB,header%idPosMat,iFile, &
               &            'Position matrix','NONE','TMAT34', &
               &            '"x","y","z"','"i1","i2","i3","position"',nBits)
          call writeVarDef (respDB,header%idVelVec,iFile, &
               &            'Velocity','LENGTH/TIME','VEC3', &
               &            '"v_x","v_y","v_z"',nBits=nBits)
          call writeVarDef (respDB,header%idAngVelVec,iFile, &
               &            'Angular velocity','ANGLE/TIME','ROT3', &
               &            '"omega_x","omega_y","omega_z"',nBits=nBits)
          call writeVarDef (respDB,header%idAccVec,iFile, &
               &            'Acceleration','LENGTH/TIME^2','VEC3', &
               &            '"a_x","a_y","a_z"',nBits=nBits)
          call writeVarDef (respDB,header%idAngAccVec,iFile, &
               &            'Angular acceleration','ANGLE/TIME^2','ROT3', &
               &            '"alpha_x","alpha_y","alpha_z"',nBits=nBits)
          if (iFile > 0) write(iFile,"(']')")

       end if

    end if
    if (supel%saveVar(2) .and. lBeam) then

       respDB%nBytes = respDB%nBytes + 12*nbs_p
       call writeVarDef (respDB,idSF(1),idatd,'Sectional force, end 1', &
            &            'FORCE','VEC3','"F_x","F_y","F_z"')
       call writeVarDef (respDB,idSF(2),idatd,'Sectional moment, end 1', &
            &            'FORCE*LENGTH','VEC3','"M_x","M_y","M_z"')
       call writeVarDef (respDB,idSF(3),idatd,'Sectional force, end 2', &
            &            'FORCE','VEC3','"F_x","F_y","F_z"')
       call writeVarDef (respDB,idSF(4),idatd,'Sectional moment, end 2', &
            &            'FORCE*LENGTH','VEC3','"M_x","M_y","M_z"')

    end if
    if (supel%saveVar(4) .and. associated(supel%hydyn)) then

       if (lBeam) then
          call writeItGDef (respDB,idHD(-1),iFile,'Hydrodynamic quantities')
       else
          call writeItGDef (respDB,idHD(0),iFile,'Hydrodynamic quantities')
       end if
       call writeVarDef (respDB,idHD(1),iFile,'Center of Buoyancy', &
            &            'LENGTH','VEC3','"x_cb","y_cb","z_cb"')
       call writeVarDef (respDB,idHD(2),iFile,'Buoyancy force', &
            &            'FORCE','VEC3','"Fb_x","Fb_y","Fb_z"')
       if (lBeam) then
          respDB%nBytes = respDB%nBytes + 15*nbs_p
          call writeVarDef (respDB,idHD(9),iFile,'Added mass force', &
               &            'FORCE','VEC3','"Fm_x","Fm_y","Fm_z"')
          call writeVarDef (respDB,idHD(4),iFile,'Drag force', &
               &            'FORCE','VEC3','"Fd_x","Fd_y","Fd_z"')
          call writeVarDef (respDB,idHD(7),iFile,'Slam force', &
               &            'FORCE','VEC3','"Fs_x","Fs_y","Fs_z"')
       else
          respDB%nBytes = respDB%nBytes + 24*nbs_p
          call writeVarDef (respDB,idHD(3),iFile,'Buoyancy moment', &
               &            'FORCE*LENGTH','ROT3','"Mb_x","Mb_y","Mb_z"')
          call writeVarDef (respDB,idHD(4),iFile,'Drag force', &
               &            'FORCE','VEC3','"Fd_x","Fd_y","Fd_z"')
          call writeVarDef (respDB,idHD(5),iFile,'Drag moment', &
               &            'FORCE*LENGTH','ROT3','"Md_x","Md_y","Md_z"')
          call writeVarDef (respDB,idHD(6),iFile,'Slam attack point', &
               &            'LENGTH','VEC3','"x_sl","y_sl","z_sl"')
          call writeVarDef (respDB,idHD(7),iFile,'Slam force', &
               &            'FORCE','VEC3','"Fs_x","Fs_y","Fs_z"')
          call writeVarDef (respDB,idHD(8),iFile,'Slam moment', &
               &            'FORCE*LENGTH','ROT3','"Ms_x","Ms_y","Ms_z"')
       end if
       if (iFile > 0) write(iFile,"(']')")

    end if
    if (supel%saveVar(2) .and. associated(supel%genDOFs)) then

       do i = 1, size(header%idGenDOF)
          if (header%idGenDOF(i) == 0) then
             respDB%nIG = respDB%nIG + 1
             header%idGenDOF(i) = respDB%nIG
             call writeItgRef (iitem,header%idGenDOF(i))
             write(iitem,600,advance='NO') i
             call writeVarDef (respDB,idAmp,iitem, &
                  &            'Value','LENGTH',nBits=nBits)
             call writeVarDef (respDB,header%idVel,iitem, &
                  &            'Velocity','LENGTH/TIME',nBits=nBits)
             call writeVarDef (respDB,header%idAcc,iitem, &
                  &            'Acceleration','LENGTH/TIME^2',nBits=nBits)
             if (supel%saveVar(3)) then
                call writeVarDef (respDB,header%idEkin,iitem, &
                     &            'Kinetic energy','ENERGY')
                call writeVarDef (respDB,header%idEstr,iitem, &
                     &            'Strain energy','ENERGY')
                call writeVarDef (respDB,header%idEdmp,iitem, &
                     &            'Energy loss','ENERGY')
             end if
             write(iitem,"(']')")
          end if
       end do

       nDOFs = supel%genDOFs%nDOFs
       newGroup = insertIdMap(nDOFs,header%idMap,respDB%nIG,idGen)
       respDB%nBytes = respDB%nBytes + 3*nDOFs*nBits/8
       if (supel%saveVar(3)) respDB%nBytes = respDB%nBytes + 3*nDOFs*nbs_p

       if (newGroup) then
          call writeItgRef (iitem,idGen)
          write(iitem,"(a)",advance='NO') ';"Component modes";'
          do i = 1, nDOFs
             if (nDOFs > 10 .and. mod(i,10) == 1) then
                write(iitem,'(/2X)',advance='NO')
             end if
             call writeItgRef (iitem,header%idGenDOF(i))
             write(iitem,"(']')",advance='NO')
          end do
          if (nDOFs > 10) then
             write(iitem,"(/']')")
          else
             write(iitem,"(']')")
          end if
       else if (idGen < 0 .and. present(ierr)) then
          ierr = idGen
          call reportError (debugFileOnly_p,'writeSupElHeader2: '// &
               &            getSupElId(supel))
          return
       end if
       call writeItgRef (idatd,idGen)

    end if
    if (supel%saveVar(3)) then

       respDB%nBytes = respDB%nBytes + 4*nbs_p
       call writeItGDef (respDB,idSen        ,iFile,'Structural energy')
       call writeVarDef (respDB,header%idEstr,iFile,'Strain energy','ENERGY')
       call writeVarDef (respDB,header%idEkin,iFile,'Kinetic energy','ENERGY')
       call writeVarDef (respDB,header%idEpot,iFile,'Potential energy','ENERGY')
       call writeVarDef (respDB,header%idEdmp,iFile,'Energy loss','ENERGY')
       if (iFile > 0) write(iFile,"(']')")

    end if
    write(idatd,"('}')")

    header%firstCall(1) = .false.

600 format(';"Mode',i3,'";')

  end subroutine writeSupElHeader2


  !!============================================================================
  !> @brief Writes results database file header for a user-defined element.
  !>
  !> @param respDB File for secondary response variables
  !> @param header Variable- and item group index container
  !> @param[in] elm The user-defined element to write file header for
  !> @param[in] nBits Size of a real variables (32 or 64 bit)
  !>
  !> @callgraph @callergraph

  subroutine writeUserDefElHeader (respDB,header,elm,nBits)

    use UserdefElTypeModule, only : UserdefElType
    use IdTypeModule       , only : writeIdHeader
    use RDBModule          , only : idatd, nbs_p, writeVarDef
    use FiUserElmInterface , only : Fi_UDE4

    type(RDBType)      , intent(inout) :: respDB
    type(HeaderId)     , intent(inout) :: header
    type(UserdefElType), intent(in)    :: elm
    integer            , intent(in)    :: nBits

    integer           :: i, nvar
    integer, save     :: idUDR(max_UDE_p)
    character(len=64) :: name

    call writeIdHeader ('User-defined element',elm%id,idatd)

    if (header%firstCall(5)) idUDR = 0

    respDB%nBytes = respDB%nBytes + 12*nBits/8
    call writeVarDef (respDB,header%idPosMat,idatd, &
         &            'Position matrix','NONE','TMAT34', &
         &            '"x","y","z"','"i1","i2","i3","position"',nBits)

    do i = 1, max_UDE_p

       call Fi_UDE4 (elm%id%baseId,elm%type,i, &
            &        elm%iwork(1),elm%rwork(1),name,nvar)
       if (nvar < i) exit

       respDB%nBytes = respDB%nBytes + nbs_p
       call writeVarDef (respDB,idUDR(i),idatd,trim(name))

    end do

    write(idatd,"('}')")

    header%firstCall(5) = .false.

  end subroutine writeUserDefElHeader


  !!============================================================================
  !> @brief Writes primary results database file header for a triad.
  !>
  !> @param respDB File for primary response variables
  !> @param header Variable- and item group index container
  !> @param[in] triad The triad to write file header for
  !> @param[in] nBits Size of a real variables (32 or 64 bit)
  !>
  !> @callgraph @callergraph

  subroutine writeTriadHeader1 (respDB,header,triad,nBits)

    use TriadTypeModule, only : TriadType
    use IdTypeModule   , only : writeIdHeader
    use RDBModule      , only : idatd, writeVarDef

    type(RDBType)  , intent(inout) :: respDB
    type(HeaderId) , intent(inout) :: header
    type(TriadType), intent(in)    :: triad
    integer        , intent(in)    :: nBits

    if (triad%nDOFs < 1)     return ! Earth-linked triad, no results at all
    if (triad%id%baseId < 0) return ! Dummy CG triad for generic parts
    if (.not. triad%savePos) return ! Output is switched off for this triad

    call writeIdHeader ('Triad',triad%id,idatd)

    if (triad%nDOFs == 6) then

       respDB%nBytes = respDB%nBytes + 12*nBits/8
       call writeVarDef (respDB,header%idPosMat,idatd, &
            &            'Position matrix','NONE','TMAT34', &
            &            '"x","y","z"','"i1","i2","i3","position"', &
            &            nBits=nBits)

    else

       respDB%nBytes = respDB%nBytes + 3*nBits/8
       call writeVarDef (respDB,header%idPos,idatd, &
            &            'Position','LENGTH','VEC3','"r_x","r_y","r_z"', &
            &            nBits=nBits)

    end if

    write(idatd,"('}')")

  end subroutine writeTriadHeader1


  !!============================================================================
  !> @brief Writes secondary results database file header for a triad.
  !>
  !> @param respDB File for primary response variables
  !> @param header Variable- and item group index container
  !> @param[in] triad The triad to write file header for
  !> @param[in] nBits Size of a real variables (32 or 64 bit)
  !> @param[in] includePosMat If .true., include the position matrix also
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph

  subroutine writeTriadHeader2 (respDB,header,triad,nBits,includePosMat,ierr)

    use TriadTypeModule    , only : TriadType
    use IdTypeModule       , only : writeIdHeader, getId
    use HydroDynamicsModule, only : getCalculatedFluidMotion
    use RDBModule          , only : iitem, idatd, nbs_p
    use RDBModule          , only : writeVarDef, writeItgDef
    use reportErrorModule  , only : reportError, debugFileOnly_p

    type(RDBType)   , intent(inout) :: respDB
    type(HeaderId)  , intent(inout) :: header
    type(TriadType) , intent(in)    :: triad
    integer         , intent(in)    :: nBits
    logical,optional, intent(in)    :: includePosMat
    integer,optional, intent(out)   :: ierr

    integer, parameter :: digits_p(5) = (/ 1,2,4,8,16 /)

    logical       :: newGroup
    integer       :: iFile, iPattern, idLocal, nBv
    integer, save :: idVar(7), idF(5), idM(5), idHD(3), idAD(0:3)
    integer, save :: id3D = 0, id6D = 0

    if (present(ierr)) ierr = 0
    if (triad%nDOFs < 1)          return ! Earth-linked triad, no results at all
    if (triad%id%baseId < 0)      return ! Dummy CG triad for generic parts
    if (.not. any(triad%saveVar)) return ! No secondary results for this triad

    call writeIdHeader ('Triad',triad%id,idatd)

    if (header%firstCall(2)) then
       idVar = 0
       idF   = 0
       idM   = 0
       idHD  = 0
       idAD  = 0
       id3D  = 0
       id6D  = 0
    end if

    if (present(includePosMat)) then
       if (includePosMat) then
          !! For creation of response_pos.frs only
          respDB%nBytes = respDB%nBytes + 12*nbs_p
          call writeVarDef (respDB,header%idPosMat,idatd, &
               &            'Position matrix','NONE','TMAT34', &
               &            '"x","y","z"','"i1","i2","i3","position"')
       end if
    end if

    nBv = 3*nBits/8
    if (triad%saveVar(7) .and. associated(triad%ur_def)) then

       respDB%nBytes = respDB%nBytes + nBv
       call writeVarDef (respDB,idVar(4),idatd, &
            &            'Deformational displacement','LENGTH','VEC3', &
            &            '"d_x","d_y","d_z"',nBits=nBits)

       if (triad%nDOFs == 6 .and. size(triad%ur_def) >= 6) then
          respDB%nBytes = respDB%nBytes + nBv
          call writeVarDef (respDB,idVar(5),idatd, &
               &            'Deformational rotation','ANGLE','ROT3', &
               &            '"theta_x","theta_y","theta_z"',nBits=nBits)
       end if

       if (size(triad%ur_def) >= 12) then
          respDB%nBytes = respDB%nBytes + 2*nBv
          call writeVarDef (respDB,idVar(6),idatd, &
               &            'Relative velocity','LENGTH','VEC3', &
               &            '"v_x","v_y","v_z"',nBits=nBits)
          call writeVarDef (respDB,idVar(7),idatd, &
               &            'Relative angular velocity','LENGTH','VEC3', &
               &            '"omega_x","omega_y","omega_z"',nBits=nBits)
       end if

    end if
    if (triad%saveVar(1)) then

       respDB%nBytes = respDB%nBytes + nBv
       call writeVarDef (respDB,header%idVelVec,idatd, &
            &            'Velocity','LENGTH/TIME','VEC3', &
            &            '"v_x","v_y","v_z"',nBits=nBits)

       if (triad%nDOFs == 6) then
          respDB%nBytes = respDB%nBytes + nBv
          call writeVarDef (respDB,header%idAngVelVec,idatd, &
               &            'Angular velocity','ANGLE/TIME','ROT3', &
               &            '"omega_x","omega_y","omega_z"',nBits=nBits)
       end if

    end if
    if (triad%saveVar(2)) then

       respDB%nBytes = respDB%nBytes + nBv
       call writeVarDef (respDB,header%idAccVec,idatd, &
            &            'Acceleration','LENGTH/TIME^2','VEC3', &
            &            '"a_x","a_y","a_z"',nBits=nBits)

       if (triad%nDOFs == 6) then
          respDB%nBytes = respDB%nBytes + nBv
          call writeVarDef (respDB,header%idAngAccVec,idatd, &
               &            'Angular acceleration','ANGLE/TIME^2','ROT3', &
               &            '"alpha_x","alpha_y","alpha_z"',nBits=nBits)
       end if

    end if
    if (triad%saveVar(3)) then
       if (associated(triad%scoord)) then
          iPattern = 32
          respDB%nBytes = respDB%nBytes + 6*nbs_p
       else
          iPattern = 0
       end if
       if (associated(triad%supForce) .and. triad%isFEnode) then
          iPattern = iPattern + 64
          respDB%nBytes = respDB%nBytes + nBv
          if (triad%nDOFs == 6) then
             iPattern = iPattern + 128
             respDB%nBytes = respDB%nBytes + nBv
          end if
       end if
       if (iPattern > 0) then
          newGroup = insertIdMap(iPattern,header%idMap,respDB%nIG,idLocal)
          if (newGroup) then
             call writeItgRef (iitem,idLocal)
             write(iitem,600,advance='NO') 'Beam quantities'

             if (associated(triad%scoord)) then

                call writeVarDef (respDB,header%idCLen,iitem, &
                     &            'Curve length','LENGTH')
                call writeVarDef (respDB,idVar(1),iitem, &
                     &            'Curvature','NONE/LENGTH', &
                     &            'VEC3','"kappa_x","kappa_y","kappa_z"')
                call writeVarDef (respDB,idVar(2),iitem, &
                     &            'My(kappa)','FORCE/LENGTH')
                call writeVarDef (respDB,idVar(3),iitem, &
                     &            'Mz(kappa)','FORCE/LENGTH')

             end if
             if (associated(triad%supForce) .and. triad%isFEnode) then

                call writeVarDef (respDB,header%idForceVec,iitem, &
                     &            'Force','FORCE','VEC3', &
                     &            '"F_x","F_y","F_z"',nBits=nBits)

                if (triad%nDOFs == 6) then
                   call writeVarDef (respDB,header%idMomentVec,iitem, &
                        &            'Moment','FORCE*LENGTH','ROT3', &
                        &            '"M_x","M_y","M_z"',nBits=nBits)
                end if

             end if
             write(iitem,"(']')")

          else if (idLocal < 0) then
             ierr = idLocal
             call reportError (debugFileOnly_p,'writeTriadHeader2: Triad'// &
                  &            getId(triad%id))
             return
          end if
          call writeItgRef (idatd,idLocal)

       end if
       if (associated(triad%supForce) .and. .not. triad%isFEnode) then

          respDB%nBytes = respDB%nBytes + nBv
          call writeVarDef (respDB,header%idForceVec,idatd, &
               &            'Force','FORCE','VEC3', &
               &            '"F_x","F_y","F_z"',nBits=nBits)

          if (triad%nDOFs == 6) then
             respDB%nBytes = respDB%nBytes + nBv
             call writeVarDef (respDB,header%idMomentVec,idatd, &
                  &            'Moment','FORCE*LENGTH','ROT3', &
                  &            '"M_x","M_y","M_z"',nBits=nBits)
          end if

       end if
       if (associated(triad%nodeForce)) then

          respDB%nBytes = respDB%nBytes + nBv
          call writeVarDef (respDB,idF(1),idatd,'Total force','FORCE', &
               &            'VEC3','"F_x","F_y","F_z"',nBits=nBits)

          if (triad%nDOFs == 6) then
             respDB%nBytes = respDB%nBytes + nBv
             call writeVarDef (respDB,idM(1),idatd,'Total moment', &
                  &            'FORCE*LENGTH','ROT3', &
                  &            '"M_x","M_y","M_z"',nBits=nBits)
          end if

       end if
    end if
    if (any(triad%saveVar(4:6))) then

       if (triad%saveVar(4)) then
          respDB%nBytes = respDB%nBytes + nBv
          if (triad%nDOFs == 6) respDB%nBytes = respDB%nBytes + nBv
       end if
       if (triad%saveVar(5)) then
          respDB%nBytes = respDB%nBytes + nBv
          if (triad%nDOFs == 6) respDB%nBytes = respDB%nBytes + nBv
       end if
       if (triad%saveVar(6) .and. associated(triad%supForce)) then
          if (.not. triad%isFEnode) then ! not relevant for section forces
             respDB%nBytes = respDB%nBytes + nBv
             if (triad%nDOFs == 6) respDB%nBytes = respDB%nBytes + nBv
          end if
       end if
       if (triad%saveVar(6) .and. associated(triad%nodeForce)) then
          respDB%nBytes = respDB%nBytes + nBv
          if (triad%nDOFs == 6) respDB%nBytes = respDB%nBytes + nBv
       end if

       iPattern = sum(digits_p(1:3),triad%saveVar(4:6))
       if (triad%isFEnode .or. .not.associated(triad%supForce)) then
          iPattern = iPattern - digits_p(3)
       end if
       if (triad%saveVar(6) .and. associated(triad%nodeForce)) then
          iPattern = iPattern + digits_p(4)
       end if
       if (triad%nDOFs == 6) then
          iPattern = iPattern + digits_p(5)
       end if
       newGroup = insertIdMap(iPattern,header%idMap,respDB%nIG,idLocal)

       if (newGroup) then
          call writeItgRef (iitem,idLocal)
          write(iitem,600,advance='NO') 'Local coordinates'

          if (triad%saveVar(4)) then

             call writeVarDef (respDB,header%idVelVec,iitem, &
                  &            'Velocity','LENGTH/TIME', &
                  &            'VEC3','"v_x","v_y","v_z"', &
                  &            nBits=nBits)

             if (triad%nDOFs == 6) then
                call writeVarDef (respDB,header%idAngVelVec,iitem, &
                     &            'Angular velocity','ANGLE/TIME', &
                     &            'ROT3','"v_x","v_y","v_z"', &
                     &            nBits=nBits)
             end if

          end if
          if (triad%saveVar(5)) then

             call writeVarDef (respDB,header%idAccVec,iitem, &
                  &            'Acceleration','LENGTH/TIME^2', &
                  &            'VEC3','"a_x","a_y","a_z"', &
                  &            nBits=nBits)

             if (triad%nDOFs == 6) then
                call writeVarDef (respDB,header%idAngAccVec,iitem, &
                     &            'Angular acceleration','ANGLE/TIME^2', &
                     &            'ROT3','"alpha_x","alpha_y","alpha_z"', &
                     &            nBits=nBits)
             end if

          end if
          if (triad%saveVar(6) .and. associated(triad%supForce)) then
             if (.not. triad%isFEnode) then ! not relevant for section forces

                call writeVarDef (respDB,header%idForceVec,iitem, &
                     &            'Force','FORCE', &
                     &            'VEC3','"F_x","F_y","F_z"', &
                     &            nBits=nBits)

                if (triad%nDOFs == 6) then
                   call writeVarDef (respDB,header%idMomentVec,iitem, &
                        &            'Moment','FORCE*LENGTH', &
                        &            'ROT3','"M_x","M_y","M_z"', &
                        &            nBits=nBits)
                end if

             end if
          end if
          if (triad%saveVar(6) .and. associated(triad%nodeForce)) then

             call writeVarDef (respDB,idF(1),iitem, &
                  &            'Total force','FORCE', &
                  &            'VEC3','"F_x","F_y","F_z"', &
                  &            nBits=nBits)

             if (triad%nDOFs == 6) then
                call writeVarDef (respDB,idM(1),iitem, &
                     &            'Total moment','FORCE*LENGTH', &
                     &            'ROT3','"M_x","M_y","M_z"', &
                     &            nBits=nBits)
             end if

          end if
          write(iitem,"(']')")

       else if (idLocal < 0) then
          ierr = idLocal
          call reportError (debugFileOnly_p,'writeTriadHeader2: Triad'// &
               &            getId(triad%id))
          return
       end if
       call writeItgRef (idatd,idLocal)

    end if
    if (triad%saveVar(9)) then
       if (getCalculatedFluidMotion(triad)) then

          respDB%nBytes = respDB%nBytes + 7*nbs_p
          call writeVarDef (respDB,idHD(1),idatd,'Wave elevation','LENGTH')
          call writeVarDef (respDB,idHD(2),idatd, &
               &            'Fluid velocity','LENGTH/TIME', &
               &            'VEC3','"v_x","v_y","v_z"')
          call writeVarDef (respDB,idHD(3),idatd, &
               &            'Fluid acceleration','LENGTH/TIME^2', &
               &            'VEC3','"a_x","a_y","a_z"')

       end if
    end if
    if (triad%saveVar(6) .and. associated(triad%aeroForce)) then

       respDB%nBytes = respDB%nBytes + 3*nbs_p
       call writeItGDef (respDB,idAD(0),iFile,'Aerodynamic forces')
       call writeVarDef (respDB,idAD(1),iFile,'Normal force','FORCE')
       call writeVarDef (respDB,idAD(2),iFile,'Tangential force','FORCE')
       call writeVarDef (respDB,idAD(3),iFile,'Pitching moment','FORCE*LENGTH')
       if (iFile > 0) write(iFile,"(']')")

    end if
    if (triad%saveVar(8)) then

       if (triad%nDOFs == 6) then
          respDB%nBytes = respDB%nBytes + 24*nbs_p
          call writeItGDef (respDB,id6D,iFile,'Force contributions')
       else
          respDB%nBytes = respDB%nBytes + 12*nbs_p
          call writeItGDef (respDB,id3D,iFile,'Force contributions')
       end if

       call writeVarDef (respDB,idF(2),iFile, &
            &            'F_s','FORCE','VEC3','"F_x","F_y","F_z"')
       if (triad%nDOFs == 6) then
          call writeVarDef (respDB,idM(2),iFile, &
               &            'M_s','FORCE*LENGTH','ROT3','"M_x","M_y","M_z"')
       end if

       call writeVarDef (respDB,idF(3),iFile,'F_d', &
            &            'FORCE','VEC3','"F_x","F_y","F_z"')
       if (triad%nDOFs == 6) then
          call writeVarDef (respDB,idM(3),iFile, &
               &            'M_d','FORCE*LENGTH','ROT3','"M_x","M_y","M_z"')
       end if

       call writeVarDef (respDB,idF(4),iFile, &
            &            'F_i','FORCE','VEC3','"F_x","F_y","F_z"')
       if (triad%nDOFs == 6) then
          call writeVarDef (respDB,idM(4),iFile, &
               &            'M_i','FORCE*LENGTH','ROT3','"M_x","M_y","M_z"')
       end if

       call writeVarDef (respDB,idF(5),iFile, &
            &            'F_ex','FORCE','VEC3','"F_x","F_y","F_z"')
       if (triad%nDOFs == 6) then
          call writeVarDef (respDB,idM(5),iFile, &
               &            'M_ex','FORCE*LENGTH','ROT3','"M_x","M_y","M_z"')
       end if

       if (iFile > 0) write(iFile,"(']')")

    end if

    write(idatd,"('}')")

    header%firstCall(2) = .false.

600 format(';"',a,'";')

  end subroutine writeTriadHeader2


  !!============================================================================
  !> @brief Writes results database file header for an axial spring.
  !>
  !> @param respDB File for secondary response variables
  !> @param header Variable- and item group index container
  !> @param[in] spring The axial spring to write file header for
  !> @param[in] nBits Size of a real variables (32 or 64 bit)
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph

  subroutine writeSpringElementHeader (respDB,header,spring,nBits,ierr)

    use SpringTypeModule, only : SpringType
    use IdTypeModule    , only : writeIdHeader
    use RDBModule       , only : idatd

    type(RDBType)   , intent(inout) :: respDB
    type(HeaderId)  , intent(inout) :: header
    type(SpringType), intent(in)    :: spring
    integer         , intent(in)    :: nBits
    integer         , intent(out)   :: ierr

    logical :: noResults, isAxialSpring
    integer :: i, nSpr

    !! Check if results are requested for this spring element
    ierr = 0
    noResults = .true.
    isAxialSpring = .false.
    nSpr = size(spring%spr)
    do i = 1, nSpr
       if (associated(spring%spr(i)%p)) then
          if (any(spring%spr(i)%p%saveVar)) then
             noResults = .false.
             isAxialSpring = i*nSpr == 1 .and. spring%spr(i)%p%dof < 1
          end if
       end if
    end do
    if (noResults) return

    if (isAxialSpring) then

       !! Axial spring, goes directly into the data section
       call writeIdHeader ('Axial spring',spring%id,idatd)
       call writeSpringHeader (respDB,header,spring%spr(1)%p,0,nBits,ierr)

    else

       !! Write the actual spring element to the data section
       call writeIdHeader ('Free joint',spring%id,idatd)
       do i = 1, nSpr
          if (associated(spring%spr(i)%p)) then
             call writeSpringHeader (respDB,header,spring%spr(i)%p,i,nBits,ierr)
             if (ierr < 0) exit
          end if
       end do

    end if
    write(idatd,"('}')")

  end subroutine writeSpringElementHeader


  !!============================================================================
  !> @brief Writes results database file header for a spring.
  !>
  !> @param respDB File for secondary response variables
  !> @param header Variable- and item group index container
  !> @param[in] spring The spring to write file header for
  !> @param[in] sDof Local DOF index
  !> @param[in] nBits Size of a real variables (32 or 64 bit)
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph

  subroutine writeSpringHeader (respDB,header,spring,sDof,nBits,ierr)

    use SpringTypeModule , only : SpringBaseType
    use IdTypeModule     , only : getId
    use RDBModule        , only : iitem, idatd, nbs_p, writeVarDef
    use reportErrorModule, only : reportError, debugFileOnly_p

    type(RDBType)       , intent(inout) :: respDB
    type(HeaderId)      , intent(inout) :: header
    type(SpringBaseType), intent(in)    :: spring
    integer             , intent(in)    :: sDof, nBits
    integer             , intent(out)   :: ierr

    integer, parameter :: digits_p(10) = &
         (/ 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096 /)

    integer :: iFile, iPattern, idSpring
    logical :: saveVar(10), haveYield

    ierr = 0
    if (.not. any(spring%saveVar)) return ! No results requested for this spring

    haveYield = spring%unLoadType == 1 .or. spring%unLoadType == 2 .or. &
         &      associated(spring%yield)

    saveVar(1) = spring%saveVar(1)
    saveVar(2) = spring%saveVar(2) .and. associated(spring%length0Engine)
    saveVar(3) = spring%saveVar(2)
    saveVar(4) = spring%saveVar(3)
    saveVar(5) = spring%saveVar(3) .and. haveYield
    saveVar(6) = spring%saveVar(4)
    saveVar(7) = spring%saveVar(5)
    saveVar(8) = spring%saveVar(5) .and. associated(spring%length0Engine)
    saveVar(9) = spring%saveVar(5) .and. haveYield
    saveVar(10) = spring%saveVar(1) .and. spring%unLoadType == 3

    respDB%nBytes = respDB%nBytes + count(saveVar(1:6))*nBits/8
    if (saveVar(5)) respDB%nBytes = respDB%nBytes + nBits/8
    respDB%nBytes = respDB%nBytes + count(saveVar(7:9))*nbs_p
    if (saveVar(10)) respDB%nBytes = respDB%nBytes + nBits/8

    if (sDof > 0) then

       !! This is a joint- or contact element spring, check if an item group
       !! has been defined for the current variable configuration already
       iPattern = sDof + sum(digits_p,saveVar)
       if (insertIdMap(iPattern,header%idSpringMap,respDB%nIG,idSpring)) then
          !! Write new item group definition and reference to the new group
          call writeDofIG (idSpring,sDof,'spring')
          iFile = iitem
       else if (idSpring > 0) then
          !! Write reference to the existing item group
          call writeDofIG (idSpring)
          return
       else
          ierr = idSpring
          call reportError (debugFileOnly_p,'writeSpringHeader: Spring'// &
               &            getId(spring%id))
          return
       end if

    else
       iFile = idatd ! Axial spring, no item groups
    end if

    !! Write spring variable definitions

    if (sDof <= 3 .or. sDof > 6) then

       !! Translational and slider spring variables

       if (saveVar(1)) then

          call writeVarDef (respDB,header%idStiff,iFile, &
               &            'Spring stiffness','FORCE/LENGTH',nBits=nBits)

       end if
       if (saveVar(10)) then

          call writeVarDef (respDB,header%idSecStiff,iFile, &
               &            'Secant stiffness','FORCE/LENGTH',nBits=nBits)

       end if
       if (saveVar(2)) then

          call writeVarDef (respDB,header%idLength0,iFile, &
               &            'Stress free length','LENGTH',nBits=nBits)

       end if
       if (saveVar(3)) then

          call writeVarDef (respDB,header%idLength,iFile, &
               &            'Length','LENGTH',nBits=nBits)

       end if
       if (saveVar(4)) then

          call writeVarDef (respDB,header%idDefl,iFile, &
               &            'Deflection','LENGTH',nBits=nBits)

       end if
       if (saveVar(5)) then

          call writeVarDef (respDB,header%idYieldDefl,iFile, &
               &            'Yield deflection','LENGTH',nBits=nBits)
          call writeVarDef (respDB,header%idElasticDefl,iFile, &
               &            'Elastic deflection','LENGTH',nBits=nBits)

       end if
       if (saveVar(6)) then

          call writeVarDef (respDB,header%idForce,iFile, &
               &            'Force value','FORCE',nBits=nBits)

       end if

    else

       !! Rotational spring variables

       if (saveVar(1)) then

          call writeVarDef (respDB,header%idAngStiff,iFile, &
               &            'Spring stiffness','FORCE*LENGTH/ANGLE',nBits=nBits)

       end if
       if (saveVar(10)) then

          call writeVarDef (respDB,header%idAngSecStiff,iFile, &
               &            'Secant stiffness','FORCE*LENGTH/ANGLE',nBits=nBits)

       end if
       if (saveVar(2)) then

          call writeVarDef (respDB,header%idAngle0,iFile, &
               &            'Stress free angle','ANGLE',nBits=nBits)

       end if
       if (saveVar(3)) then

          call writeVarDef (respDB,header%idAngle,iFile, &
               &            'Angle','ANGLE',nBits=nBits)

       end if
       if (saveVar(4)) then

          call writeVarDef (respDB,header%idAngDefl,iFile, &
               &            'Angular deflection','ANGLE',nBits=nBits)

       end if
       if (saveVar(5)) then

          call writeVarDef (respDB,header%idYieldAngDefl,iFile, &
               &            'Angular yield deflection','ANGLE',nBits=nBits)
          call writeVarDef (respDB,header%idElasticAngDefl,iFile, &
               &            'Angular elastic deflection','ANGLE',nBits=nBits)

       end if
       if (saveVar(6)) then

          call writeVarDef (respDB,header%idMoment,iFile, &
               &            'Moment value','FORCE*LENGTH',nBits=nBits)

       end if

    end if

    if (saveVar(7)) then

       call writeVarDef (respDB,header%idEstr,iFile, &
            &            'Strain energy','ENERGY')

    end if
    if (saveVar(8)) then

       call writeVarDef (respDB,header%idEinp,iFile, &
            &            'Energy input','ENERGY')

    end if
    if (saveVar(9)) then

       call writeVarDef (respDB,header%idEdmp,iFile, &
            &            'Energy loss','ENERGY')

    end if

    if (sDof > 0) write(iFile,"(']')")

  end subroutine writeSpringHeader


  !!============================================================================
  !> @brief Writes results database file header for an axial damper.
  !>
  !> @param respDB File for secondary response variables
  !> @param header Variable- and item group index container
  !> @param[in] damper The axial damper to write file header for
  !> @param[in] nBits Size of a real variables (32 or 64 bit)
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph

  subroutine writeAxialDamperHeader (respDB,header,damper,nBits,ierr)

    use DamperTypeModule, only : DamperType
    use IdTypeModule    , only : writeIdHeader
    use RDBModule       , only : idatd

    type(RDBType)   , intent(inout) :: respDB
    type(HeaderId)  , intent(inout) :: header
    type(DamperType), intent(in)    :: damper
    integer         , intent(in)    :: nBits
    integer         , intent(out)   :: ierr

    if (any(damper%dmp%saveVar)) then

       !! Axial damper, goes directly into the data section
       call writeIdHeader ('Axial damper',damper%dmp%id,idatd)
       call writeDamperHeader (respDB,header,damper%dmp,0,nBits,ierr)
       write(idatd,"('}')")

    end if

  end subroutine writeAxialDamperHeader


  !!============================================================================
  !> @brief Writes results database file header for a damper.
  !>
  !> @param respDB File for secondary response variables
  !> @param header Variable- and item group index container
  !> @param[in] damper The damper to write file header for
  !> @param[in] dDof Local DOF index
  !> @param[in] nBits Size of a real variables (32 or 64 bit)
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph

  subroutine writeDamperHeader (respDB,header,damper,dDof,nBits,ierr)

    use DamperTypeModule , only : DamperBaseType
    use IdTypeModule     , only : getId
    use RDBModule        , only : iitem, idatd, nbs_p, writeVarDef
    use reportErrorModule, only : reportError, debugFileOnly_p

    type(RDBType)       , intent(inout) :: respDB
    type(HeaderId)      , intent(inout) :: header
    type(DamperBaseType), intent(in)    :: damper
    integer             , intent(in)    :: dDof, nBits
    integer             , intent(out)   :: ierr

    integer, parameter :: digits_p(5) = (/ 8,16,32,64,128 /)

    integer :: iFile, iPattern, idDamper

    ierr = 0
    if (.not. any(damper%saveVar)) return ! No results requested for this damper

    respDB%nBytes = respDB%nBytes + count(damper%saveVar(1:4))*nBits/8
    if (damper%saveVar(5)) respDB%nBytes = respDB%nBytes + nbs_p

    if (dDof > 0) then

       !! This is a joint- or contact element damper, check if an item group
       !! has been defined for the current variable configuration already
       iPattern = dDof + sum(digits_p,damper%saveVar)
       if (insertIdMap(iPattern,header%idDamperMap,respDB%nIG,idDamper)) then
          !! Write new item group definition and reference to the new group
          call writeDofIG (idDamper,dDof,'damper')
          iFile = iitem
       else if (idDamper > 0) then
          !! Write reference to the existing item group
          call writeDofIG (idDamper)
          return
       else
          ierr = idDamper
          call reportError (debugFileOnly_p,'writeDamperHeader: Damper'// &
               &            getId(damper%id))
          return
       end if

    else
       iFile = idatd ! Axial damper, no item groups
    end if

    ! Write damper variable definitions

    if (dDof <= 3 .or. dDof > 6) then

       ! Translational and slider damper variables

       if (damper%saveVar(1)) then

          call writeVarDef (respDB,header%idCoeff,iFile, &
               &            'Damper coefficient','FORCE*TIME/LENGTH', &
               &            nBits=nBits)

       end if
       if (damper%saveVar(2)) then

          call writeVarDef (respDB,header%idLength,iFile, &
               &            'Length','LENGTH',nBits=nBits)

       end if
       if (damper%saveVar(3)) then

          call writeVarDef (respDB,header%idVel,iFile, &
               &            'Velocity','LENGTH/TIME',nBits=nBits)

       end if
       if (damper%saveVar(4)) then

          call writeVarDef (respDB,header%idForce,iFile, &
               &            'Force value','FORCE',nBits=nBits)

       end if

    else

       ! Rotational damper variables

       if (damper%saveVar(1)) then

          call writeVarDef (respDB,header%idAngCoeff,iFile, &
               &            'Damper coefficient','FORCE*LENGTH*TIME/ANGLE', &
               &            nBits=nBits)

       end if
       if (damper%saveVar(2)) then

          call writeVarDef (respDB,header%idAngle,iFile, &
               &            'Angle','ANGLE',nBits=nBits)

       end if
       if (damper%saveVar(3)) then

          call writeVarDef (respDB,header%idAngVel,iFile, &
               &            'Angular velocity','ANGLE/TIME',nBits=nBits)

       end if
       if (damper%saveVar(4)) then

          call writeVarDef (respDB,header%idMoment,iFile, &
               &            'Moment value','FORCE*LENGTH',nBits=nBits)

       end if

    end if

    if (damper%saveVar(5)) then

       call writeVarDef (respDB,header%idEdmp,iFile, &
            &            'Energy loss','ENERGY')

    end if

    if (dDof > 0) write(iFile,"(']')")

  end subroutine writeDamperHeader


  !!============================================================================
  !> @brief Writes results database file header for a friction.
  !>
  !> @param respDB File for secondary response variables
  !> @param header Variable- and item group index container
  !> @param[in] friction The friction to write file header for
  !> @param[in] jointDof Local DOF index
  !> @param[in] nBits Size of a real variables (32 or 64 bit)
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph

  subroutine writeFrictionHeader (respDB,header,friction,jointDof,nBits,ierr)

    use FrictionTypeModule, only : FrictionType
    use IdTypeModule      , only : getId
    use RDBModule         , only : iitem, nbs_p, nbi_p, writeVarDef
    use reportErrorModule , only : reportError, debugFileOnly_p

    type(RDBType)     , intent(inout) :: respDB
    type(HeaderId)    , intent(inout) :: header
    type(FrictionType), intent(in)    :: friction
    integer           , intent(in)    :: jointDof, nBits
    integer           , intent(out)   :: ierr

    integer, parameter :: digits_p(3) = (/ 8,16,32 /)

    integer       :: iPattern, idFriction
    integer, save :: idFricVar(9)

    ierr = 0
    if (.not. any(friction%saveVar)) return ! No results requested for friction

    if (header%firstCall(3)) idFricVar = 0
    if (friction%saveVar(1)) respDB%nBytes = respDB%nBytes + nBits/8 + 2*nbs_p
    if (friction%saveVar(2)) respDB%nBytes = respDB%nBytes + nbs_p
    if (friction%saveVar(3)) respDB%nBytes = respDB%nBytes + 3*nbs_p + nbi_p

    !! Check if an item group has been defined for the current friction variable
    !! configuration already
    iPattern = jointDof + sum(digits_p,friction%saveVar)
    if (insertIdMap(iPattern,header%idFrictionMap,respDB%nIG,idFriction)) then
       !! Write new item group definition and reference to the new group
       call writeDofIG (idFriction,jointDof,'friction')
    else if (idFriction > 0) then
       !! Write reference to the existing item group
       call writeDofIG (idFriction)
       return
    else
       ierr = idFriction
       call reportError (debugFileOnly_p,'writeSpringHeader: Friction'// &
            &            getId(friction%param%id))
       return
    end if

    ! Write friction variable definitions

    if (jointDof <= 3 .or. jointDof > 6) then

       ! Translational and slider friction variables

       if (friction%saveVar(1)) then

          call writeVarDef (respDB,header%idForce,iitem, &
               &            'Force value','FORCE',nBits=nBits)
          call writeVarDef (respDB,idFricVar(1),iitem, &
               &            'Equivalent normal load','FORCE')
          call writeVarDef (respDB,idFricVar(2),iitem,'Friction coefficient')

       end if
       if (friction%saveVar(3)) then

          call writeVarDef (respDB,idFricVar(3),iitem,'Fmax','FORCE')
          call writeVarDef (respDB,idFricVar(4),iitem,'pos' ,'LENGTH')
          call writeVarDef (respDB,idFricVar(5),iitem,'vel' ,'LENGTH/TIME')

       end if

    else

       ! Rotational joint friction variables

       if (friction%saveVar(1)) then

          call writeVarDef (respDB,header%idMoment,iitem, &
               &            'Moment value','FORCE*LENGTH',nBits=nBits)
          call writeVarDef (respDB,idFricVar(1),iitem, &
               &            'Equivalent normal load','FORCE')
          call writeVarDef (respDB,idFricVar(2),iitem,'Friction coefficient')

       end if
       if (friction%saveVar(3)) then

          call writeVarDef (respDB,idFricVar(6),iitem,'Mmax','FORCE*LENGTH')
          call writeVarDef (respDB,idFricVar(7),iitem,'pos' ,'ANGLE')
          call writeVarDef (respDB,idFricVar(8),iitem,'vel' ,'ANGLE/TIME')

       end if

    end if

    if (friction%saveVar(2)) then

       call writeVarDef (respDB,header%idEdmp,iitem, &
            &            'Energy loss','ENERGY')

    end if
    if (friction%saveVar(3)) then

       call writeVarDef (respDB,idFricVar(9),iitem,'Stick/Slip','NONE','NUMBER')

    end if

    write(iitem,"(']')")

    header%firstCall(3) = .false.

  end subroutine writeFrictionHeader


  !!============================================================================
  !> @brief Writes results database file header for an external force.
  !>
  !> @param respDB File for secondary response variables
  !> @param header Variable- and item group index container
  !> @param[in] force The external force to write file header for
  !> @param[in] nBits Size of a real variables (32 or 64 bit)
  !>
  !> @callgraph @callergraph

  subroutine writeForceHeader (respDB,header,force,nBits)

    use ForceTypeModule, only : ForceType, forceVec_p, momentVec_p
    use IdTypeModule   , only : writeIdHeader
    use RDBModule      , only : idatd, nbs_p, writeVarDef

    type(RDBType)  , intent(inout) :: respDB
    type(HeaderId) , intent(inout) :: header
    type(ForceType), intent(in)    :: force
    integer        , intent(in)    :: nBits

    if (.not. any(force%saveVar)) return ! No results requested for this force

    if (force%dof == forceVec_p) then

       call writeIdHeader ('Force',force%id,idatd)
       if (force%saveVar(1)) then

          respDB%nBytes = respDB%nBytes + 3*nBits/8
          call writeVarDef (respDB,header%idForceVec,idatd, &
               &            'Force','FORCE','VEC3', &
               &            '"F_x","F_y","F_z"',nBits=nBits)

       end if
       if (force%saveVar(2)) then

          respDB%nBytes = respDB%nBytes + nBits/8
          call writeVarDef (respDB,header%idForce,idatd, &
               &            'Force value','FORCE',nBits=nBits)

       end if

    else if (force%dof == momentVec_p) then

       call writeIdHeader ('Torque',force%id,idatd)
       if (force%saveVar(1)) then

          respDB%nBytes = respDB%nBytes + 3*nBits/8
          call writeVarDef (respDB,header%idMomentVec,idatd, &
               &            'Moment','FORCE*LENGTH','ROT3', &
               &            '"M_x","M_y","M_z"',nBits=nBits)

       end if
       if (force%saveVar(2)) then

          respDB%nBytes = respDB%nBytes + nBits/8
          call writeVarDef (respDB,header%idMoment,idatd, &
               &            'Moment value','FORCE*LENGTH',nBits=nBits)

       end if

    else
       return ! Joint and triad DOF forces are not saved
    end if

    if (force%saveVar(3)) then

       respDB%nBytes = respDB%nBytes + nbs_p
       call writeVarDef (respDB,header%idEinp,idatd,'Energy input','ENERGY')

    end if

    write(idatd,"('}')")

  end subroutine writeForceHeader


  !!============================================================================
  !> @brief Writes results database file header for a joint.
  !>
  !> @param respDB File for secondary response variables
  !> @param header Variable- and item group index container
  !> @param[in] joint The joint to write file header for
  !> @param[in] forces All external forces in the model
  !> @param[in] motions All prescribed motions in the model
  !> @param[in] nBits Size of a real variables (32 or 64 bit)
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph

  subroutine writeJointHeader (respDB,header,joint,forces,motions,nBits,ierr)

    use MasterSlaveJointTypeModule, only : MasterSlaveJointType, JointDofType
    use MasterSlaveJointTypeModule, only : jointTypename_p, getJointId
    use ForceTypeModule           , only : ForceType
    use MotionTypeModule          , only : MotionType
    use IdTypeModule              , only : writeIdHeader
    use RDBModule                 , only : iitem, idatd, nbs_p, writeVarDef
    use reportErrorModule         , only : reportError, debugFileOnly_p

    type(RDBType)             , intent(inout) :: respDB
    type(HeaderId)            , intent(inout) :: header
    type(MasterSlaveJointType), intent(in)    :: joint
    type(ForceType)           , intent(in)    :: forces(:)
    type(MotionType)          , intent(in)    :: motions(:)
    integer                   , intent(in)    :: nBits
    integer                   , intent(out)   :: ierr

    logical :: newLine, noResults
    integer :: iDof, jDof, nJDofs

    ierr = 0
    nJDofs = size(joint%jointDofs)

    ! Check if results are requested for this joint
    noResults = .true.
    do iDof = 1, nJDofs
       if (any(joint%jointDofs(iDof)%saveVar)) noResults = .false.
       if (associated(joint%jointDofs(iDof)%spring)) then
          if (any(joint%jointDofs(iDof)%spring%saveVar)) noResults = .false.
       end if
       if (associated(joint%jointDofs(iDof)%damper)) then
          if (any(joint%jointDofs(iDof)%damper%saveVar)) noResults = .false.
       end if
       if (associated(joint%jointDofs(iDof)%friction)) then
          if (any(joint%jointDofs(iDof)%friction%saveVar)) noResults = .false.
       end if
    end do
    if (noResults) return

    ! Write the actual joint to the data section
    if (joint%type > 0 .and. joint%type <= size(jointTypeName_p)) then
       call writeIdHeader (jointTypeName_p(joint%type),joint%id,idatd,.true.)
    else
       call writeIdHeader ('Joint',joint%id,idatd,.true.)
    end if
    write(idatd,"(2X)",advance='NO')

    ! Write the joint variables
    do jDof = 1, nJDofs
       iDof = joint%dofOutputOrder(jDof)
       call writeJointVarHeader (joint%jointDofs(iDof)%lDof, &
            &                    joint%jointDofs(iDof))
       if (ierr < 0) goto 900
    end do

    ! Write the joint springs
    newLine = .true.
    do jDof = 1, nJDofs
       iDof = joint%dofOutputOrder(jDof)
       if (associated(joint%jointDofs(iDof)%spring)) then
          if (newLine) write(idatd,"(/2X)",advance='NO')
          call writeSpringHeader (respDB,header, &
               &                  joint%jointDofs(iDof)%spring, &
               &                  joint%jointDofs(iDof)%lDof,nBits,ierr)
          if (ierr < 0) goto 900
          newLine = .false.
       end if
    end do

    ! Write the joint dampers
    newLine = .true.
    do jDof = 1, nJDofs
       iDof = joint%dofOutputOrder(jDof)
       if (associated(joint%jointDofs(iDof)%damper)) then
          if (newLine) write(idatd,"(/2X)",advance='NO')
          call writeDamperHeader (respDB,header, &
               &                  joint%jointDofs(iDof)%damper, &
               &                  joint%jointDofs(iDof)%lDof,nBits,ierr)
          if (ierr < 0) goto 900
          newLine = .false.
       end if
    end do

    ! Write the joint friction
    newLine = .true.
    do jDof = 1, nJDofs
       iDof = joint%dofOutputOrder(jDof)
       if (associated(joint%jointDofs(iDof)%friction)) then
          if (newLine) write(idatd,"(/2X)",advance='NO')
          call writeFrictionHeader (respDB,header, &
               &                    joint%jointDofs(iDof)%friction, &
               &                    joint%jointDofs(iDof)%lDof,nBits,ierr)
          if (ierr < 0) goto 900
          newLine = .false.
       end if
    end do

    write(idatd,"(/'}')")
    return

900 call reportError (debugFileOnly_p,'writeJointHeader')

  contains

    !> @brief Writes an item group definition for a joint variable.
    subroutine writeJointVarHeader (jDof,jointDof)
      integer           , intent(in) :: jDof
      type(JointDofType), intent(in) :: jointDof

      integer, parameter :: digits_p(5) = (/ 8,16,32,64,128 /)

      integer :: iPattern, idJoint
      logical :: saveVar(5)

      idJoint = jointDof%loadIdx
      saveVar(1:4) = jointDof%saveVar
      if (idJoint > 0 .and. idJoint <= size(forces)) then
         saveVar(4:5) = forces(idJoint)%saveVar(2:3)
      else if (idJoint < 0 .and. -idJoint <= size(motions)) then
         saveVar(4:5) = motions(-idJoint)%saveVar(2:3)
      else if (associated(jointDof%F)) then
         saveVar(5) = .false.
      else
         saveVar(4:5) = .false.
      end if
      if (.not. any(saveVar)) return ! No results for this jointDof

      respDB%nBytes = respDB%nBytes + count(saveVar(1:4))*nBits/8
      if (saveVar(5)) respDB%nBytes = respDB%nBytes + nbs_p

      !! Check if an item group has been defined
      !! for the current joint variable configuration already
      iPattern = jDof + sum(digits_p,saveVar)
      if (insertIdMap(iPattern,header%idMap,respDB%nIG,idJoint)) then
         !! Write new item group definition and reference to the new group
         call writeDofIG (idJoint,jDof,'joint')
      else if (idJoint > 0) then
         !! Write reference to the existing item group
         call writeDofIG (idJoint)
         return
      else
         ierr = idJoint
         call reportError (debugFileOnly_p,'writeJointVarHeader: '// &
              &            getJointId(joint))
         return
      end if

      if (jDof <= 3 .or. jDof > 6) then

         if (saveVar(1)) then

            call writeVarDef (respDB,header%idLength,iitem, &
                 &            'Length','LENGTH',nBits=nBits)

         end if
         if (saveVar(2)) then

            call writeVarDef (respDB,header%idVel,iitem, &
                 &            'Velocity','LENGTH/TIME',nBits=nBits)

         end if
         if (saveVar(3)) then

            call writeVarDef (respDB,header%idAcc,iitem, &
                 &            'Acceleration','LENGTH/TIME^2',nBits=nBits)

         end if
         if (saveVar(4)) then

            call writeVarDef (respDB,header%idForce,iitem, &
                 &            'Force value','FORCE',nBits=nBits)

         end if

      else

         if (saveVar(1)) then

            call writeVarDef (respDB,header%idAngle,iitem, &
                 &            'Angle','ANGLE',nBits=nBits)

         end if
         if (saveVar(2)) then

            call writeVarDef (respDB,header%idAngVel,iitem, &
                 &            'Angular velocity','ANGLE/TIME',nBits=nBits)

         end if
         if (saveVar(3)) then

            call writeVarDef (respDB,header%idAngAcc,iitem, &
                 &            'Angular acceleration','ANGLE/TIME^2',nBits=nBits)

         end if
         if (saveVar(4)) then

            call writeVarDef (respDB,header%idMoment,iitem, &
                 &            'Moment value','FORCE*LENGTH',nBits=nBits)

         end if

      end if

      if (saveVar(5)) then

         call writeVarDef (respDB,header%idEinp,iitem,'Energy input','ENERGY')

      end if

      write(iitem,"(']')")

    end subroutine writeJointVarHeader

  end subroutine writeJointHeader


  !!============================================================================
  !> @brief Writes results database file header for a contact element.
  !>
  !> @param respDB File for secondary response variables
  !> @param header Variable- and item group index container
  !> @param[in] cElem The contact element to write file header for
  !> @param[in] nBits Size of a real variables (32 or 64 bit)
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph

  subroutine writeContactHeader (respDB,header,cElem,nBits,ierr)

    use ContactElementTypeModule, only : ContactElementType
    use IdTypeModule            , only : writeIdHeader, getId
    use RDBModule               , only : iitem, idatd, nbs_p, writeVarDef
    use reportErrorModule       , only : reportError, debugFileOnly_p

    type(RDBType)           , intent(inout) :: respDB
    type(HeaderId)          , intent(inout) :: header
    type(ContactElementType), intent(in)    :: cElem
    integer                 , intent(in)    :: nBits
    integer                 , intent(out)   :: ierr

    integer, parameter :: digits_p(2) = (/ 8,16 /)

    logical :: newLine, noResults
    integer :: iDof, jDof, iPattern, idJoint

    ! Check if results are requested for this contact element
    ierr = 0
    noResults = .not.any(cElem%saveVar)
    if (noResults) then
       do iDof = 1, size(cElem%springs)
          if (associated(cElem%springs(iDof)%p)) then
             if (any(cElem%springs(iDof)%p%saveVar)) noResults = .false.
          end if
       end do
       do iDof = 1, size(cElem%dampers)
          if (associated(cElem%dampers(iDof)%p)) then
             if (any(cElem%dampers(iDof)%p%saveVar)) noResults = .false.
          end if
       end do
    end if
    if (noResults .and. associated(cElem%friction)) then
       if (any(cElem%friction%saveVar)) noResults = .false.
    end if
    if (noResults) return

    ! Write the actual cElem to the data section
    call writeIdHeader ('Cam joint',cElem%id,idatd,.true.)
    write(idatd,"(2X)",advance='NO')

    if (any(cElem%saveVar)) then

       ! Write contact variables
       if (cElem%saveVar(1)) respDB%nBytes = respDB%nBytes + 3*nBits/8 + nbs_p
       if (cElem%saveVar(2)) respDB%nBytes = respDB%nBytes + 3*nBits/8

       do iDof = 1, 3
          jDof = iDof
          if (iDof == 3) jDof = 7 ! Slider dof

          !! Check if an item group has been defined
          !! for the current joint variable configuration already
          iPattern = jDof + sum(digits_p,cElem%saveVar)
          if (insertIdMap(iPattern,header%idMap,respDB%nIG,idJoint)) then

             !! Write new item group definition and reference to the new group
             call writeDofIG (idJoint,jDof,'joint')

             if (cElem%saveVar(1)) then

                call writeVarDef (respDB,header%idLength,iitem, &
                     &            'Length','LENGTH',nBits=nBits)

             end if
             if (cElem%saveVar(2)) then

                call writeVarDef (respDB,header%idVel,iitem, &
                     &            'Velocity','LENGTH/TIME',nBits=nBits)

             end if

             write(iitem,"(']')")

          else if (idJoint > 0) then

             !! Write reference to the existing item group
             call writeDofIG (idJoint)

          else
             ierr = idJoint
             call reportError (debugFileOnly_p,'writeContactHeader: '// &
                  &            'Contact Element'//getId(cElem%id))
             return
          end if

       end do
       if (cElem%saveVar(1)) then

          call writeVarDef (respDB,header%idAccContDis,idatd, &
               &            'Accumulated contact distance','LENGTH')

       end if
    end if

    ! Write the contact springs
    newLine = .true.
    do iDof = 1, size(cElem%springs)
       if (associated(cElem%springs(iDof)%p)) then
          jDof = iDof
          if (iDof == 3) jDof = 7 ! Slider dof
          if (newLine) write(idatd,"(/2X)",advance='NO')
          call writeSpringHeader (respDB,header,cElem%springs(iDof)%p, &
               &                  jDof,nBits,ierr)
          if (ierr < 0) goto 900
          newLine = .false.
       end if
    end do

    ! Write the contact dampers
    newLine = .true.
    do iDof = 1, size(cElem%dampers)
       if (associated(cElem%dampers(iDof)%p)) then
          jDof = iDof
          if (iDof == 3) jDof = 7 ! Slider dof
          if (newLine) write(idatd,"(/2X)",advance='NO')
          call writeDamperHeader (respDB,header,cElem%dampers(iDof)%p, &
               &                  jDof,nBits,ierr)
          if (ierr < 0) goto 900
          newLine = .false.
       end if
    end do

    ! Write the contact friction
    if (associated(cElem%friction)) then
       jDof = 7
       write(idatd,"(/2X)",advance='NO')
       call writeFrictionHeader (respDB,header,cElem%friction,jDof,nBits,ierr)
       if (ierr < 0) goto 900
    end if

    write(idatd,"(/'}')")
    return

900 call reportError (debugFileOnly_p,'writeContactHeader')

  end subroutine writeContactHeader


  !!============================================================================
  !> @brief Writes results database file header for a tire.
  !>
  !> @param respDB File for secondary response variables
  !> @param header Variable- and item group index container
  !> @param[in] tire The tire to write file header for
  !> @param[in] nBits Size of a real variables (32 or 64 bit)
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph

  subroutine writeTireHeader (respDB,header,tire,nBits,ierr)

    use TireTypeModule   , only : TireType, STI_p, MF_p
    use IdTypeModule     , only : writeIdHeader, getId
    use RDBModule        , only : idatd, nbs_p, writeVarDef, writeItgDef
#ifdef FT_DEBUG
    use TireTypeModule   , only : CTI_p
    use RDBModule        , only : nbi_p
#endif
    use reportErrorModule, only : reportError, debugFileOnly_p

    type(RDBType) , intent(inout) :: respDB
    type(HeaderId), intent(inout) :: header
    type(TireType), intent(in)    :: tire
    integer       , intent(in)    :: nBits
    integer       , intent(out)   :: ierr

    integer, parameter :: digits_p(6) = (/ 1,2,4,8,16,32 /)

    integer       :: iFile, iPattern, idTire
    integer, save :: idU(4), idTireWC(2), idVarInf(0:35)
    integer, save :: idEnergy = 0
#ifdef FT_DEBUG
    integer, save :: idIsSliding = 0, idIsContact = 0
    integer, save :: idVarCTI(10:nVarCTI_p) = (/ (0, idTire = 10, nVarCTI_p) /)
    character(len=8) :: chVar
#endif

    ierr = 0
    if (.not. any(tire%saveVar)) return ! No results requested for this tire

    call writeIdHeader ('Tire',tire%id,idatd,.false.)

    if (header%firstCall(4)) then
       idU      = 0
       idTireWC = 0
       idVarInf = 0
       idEnergy = 0
    end if

    if (tire%saveVar(1)) then

       respDB%nBytes = respDB%nBytes + 3*nbs_p
       call writeVarDef (respDB,idVarInf(1),idatd, &
            &            'Slip angle','ANGLE')
       call writeVarDef (respDB,idVarInf(2),idatd, &
            &            'Longitudinal slip','LENGTH')
       call writeVarDef (respDB,idVarInf(3),idatd, &
            &            'Camber angle','ANGLE')

    end if
    if (tire%saveVar(2)) then

       respDB%nBytes = respDB%nBytes + nbs_p
       call writeVarDef (respDB,idVarInf(12),idatd, &
            &            'Effective rolling radius','LENGTH')

    end if
    if (any(tire%saveVar(3:4))) then

       !! Write contact point item group
       iPattern = sum(digits_p(1:2),tire%saveVar(3:4))
       if (tire%api == STI_p) iPattern = iPattern + digits_p(5)
       if (insertIdMap(iPattern,header%idMap,respDB%nIG,idTire)) then
          idTire = -idTire
       else if (idTire < 0) then
          goto 900
       end if
       call writeItGDef (respDB,idTire,iFile,'Contact point')

       if (tire%saveVar(3)) then
          respDB%nBytes = respDB%nBytes + 6*nBits/8
          call writeVarDef (respDB,header%idForceVec,iFile, &
               &            'Force','FORCE','VEC3', &
               &            '"F_x","F_y","F_z"',nBits=nBits)
          call writeVarDef (respDB,header%idMomentVec,iFile, &
               &            'Moment','FORCE*LENGTH','ROT3', &
               &            '"M_x","M_y","M_z"',nBits=nBits)
       end if
       if (tire%saveVar(4)) then
          respDB%nBytes = respDB%nBytes + 6*nbs_p
          call writeVarDef (respDB,header%idPos,iFile, &
               &            'Position','LENGTH','VEC3','"r_x","r_y","r_z"')
          call writeVarDef (respDB,idVarInf(16),iFile,'Road normal', &
               &            'LENGTH','VEC3','"n_x","n_y","n_z"')
       end if

#ifdef FT_DEBUG
       respDB%nBytes = respDB%nBytes + 2*nbi_p
       call writeVarDef (respDB,idIsContact,iFile,'IsInContact','NONE','NUMBER')
       call writeVarDef (respDB,idIsSliding,iFile,'IsSliding','NONE','NUMBER')
#endif

       if (iFile > 0) write(iFile,"(']')")

    end if
    if (any(tire%saveVar(5:6))) then

       !! Write tire deflection item group
       iPattern = sum(digits_p(3:4),tire%saveVar(5:6))
       if (tire%api == STI_p) iPattern = iPattern + digits_p(5)
       if (insertIdMap(iPattern,header%idMap,respDB%nIG,idTire)) then
          idTire = -idTire
       else if (idTire < 0) then
          goto 900
       end if
       call writeItGDef (respDB,idTire,iFile,'Tire deflection')

       if (tire%saveVar(5)) then
          if (tire%api == STI_p) then
             respDB%nBytes = respDB%nBytes + 2*nbs_p
             call writeVarDef (respDB,idU(1),iFile, &
                  &            'Longitudinal','LENGTH')
             call writeVarDef (respDB,idU(2),iFile, &
                  &            'Lateral','LENGTH')
          end if
          respDB%nBytes = respDB%nBytes + nbs_p
          call writeVarDef (respDB,idVarInf(10),iFile, &
               &            'Radial','LENGTH')
       end if
       if (tire%saveVar(6) .and. tire%api == STI_p) then
          respDB%nBytes = respDB%nBytes + 3*nbs_p
          call writeVarDef (respDB,idU(3),iFile, &
               &            'Longitudinal velocity','LENGTH/TIME')
          call writeVarDef (respDB,idU(4),iFile, &
               &            'Lateral velocity','LENGTH/TIME')
          call writeVarDef (respDB,idVarInf(11),iFile, &
               &            'Closing velocity','LENGTH/TIME')
       end if

       if (iFile > 0) write(iFile,"(']')")

    end if
    if (tire%saveVar(7) .and. tire%type == MF_p) then
       respDB%nBytes = respDB%nBytes + 16*nbs_p

       !! Write tire characteristics item group
       call writeItGDef (respDB,idVarInf(0),iFile,'Tire characteristics')
       call writeVarDef (respDB,idVarInf(19),iFile, &
            &            'Longitudinal relaxation length','LENGTH')
       call writeVarDef (respDB,idVarInf(20),iFile, &
            &            'Lateral relaxation length','LENGTH')
       call writeVarDef (respDB,idVarInf(21),iFile, &
            &            'Friction coefficient X')
       call writeVarDef (respDB,idVarInf(22),iFile, &
            &            'Friction coefficient Y')
       call writeVarDef (respDB,idVarInf(24),iFile, &
            &            'Longitudinal slip velocity','LENGTH/TIME')
       call writeVarDef (respDB,idVarInf(25),iFile, &
            &            'Lateral slip velocity','LENGTH/TIME')
       call writeVarDef (respDB,idVarInf(26),iFile, &
            &            'Wheel forward velocity','LENGTH/TIME')
       call writeVarDef (respDB,idVarInf(27),iFile, &
            &            'Longitudinal deformation velocity','LENGTH/TIME')
       call writeVarDef (respDB,idVarInf(28),iFile, &
            &            'Lateral deformation velocity','LENGTH/TIME')
       call writeVarDef (respDB,idVarInf(29),iFile, &
            &            'Longitudinal slip','LENGTH')
       call writeVarDef (respDB,idVarInf(30),iFile, &
            &            'Slip angle','ANGLE')
       call writeVarDef (respDB,idVarInf(31),iFile, &
            &            'Pneumatic trail','LENGTH')
       call writeVarDef (respDB,idVarInf(32),iFile, &
            &            'Residual moment Mz','FORCE*LENGTH')
       call writeVarDef (respDB,idVarInf(33),iFile, &
            &            'Moment arm of Fx','LENGTH')
       call writeVarDef (respDB,idVarInf(34),iFile, &
            &            'Gyroscopic moment Mz','FORCE*LENGTH')
       call writeVarDef (respDB,idVarInf(35),iFile, &
            &            'Braking induced plysteer force','FORCE')
       if (iFile > 0) write(iFile,"(']')")

    end if
#ifdef FT_DEBUG
    if (tire%saveVar(7) .and. tire%api == CTI_p) then

       !! Write item group for whole varinf array (for easy debugging)
       respDB%nBytes = respDB%nBytes + (nVarCTI_p-9)*nbs_p
       call writeItGDef (respDB,idVarInf(0),iFile,'CTI variables')
       do iPattern = 10, nVarCTI_p
          write(chVar,"('Var #',I3)") iPattern
          call writeVarDef (respDB,idVarCTI(iPattern),iFile,chVar)
          if (iFile > 0 .and. mod(iPattern,10) == 0) write(iFile,"()")
       end do
       if (iFile > 0) write(iFile,"(']')")

    end if
#endif
    if (tire%saveVar(8)) then
       respDB%nBytes = respDB%nBytes + 12*nBits/8

       !! Write wheel carrier item group
       call writeItGDef (respDB,idTireWC(1),iFile,'Wheel carrier')
       call writeVarDef (respDB,header%idForceVec,iFile, &
            &            'Force','FORCE','VEC3', &
            &            '"F_x","F_y","F_z"',nBits=nBits)
       call writeVarDef (respDB,header%idMomentVec,iFile, &
            &            'Moment','FORCE*LENGTH','ROT3', &
            &            '"M_x","M_y","M_z"',nBits=nBits)
       if (iFile > 0) write(iFile,"(']')")

       !! Write wheel carrier item group, local coordinate system
       call writeItGDef (respDB,idTireWC(2),iFile,'Wheel carrier, local system')
       call writeVarDef (respDB,header%idForceVec,iFile, &
            &            'Force','FORCE','VEC3', &
            &            '"F_x","F_y","F_z"',nBits=nBits)
       call writeVarDef (respDB,header%idMomentVec,iFile, &
            &            'Moment','FORCE*LENGTH','ROT3', &
            &            '"M_x","M_y","M_z"',nBits=nBits)
       if (iFile > 0) write(iFile,"(']')")

    end if
    if (tire%saveVar(9)) then
       respDB%nBytes = respDB%nBytes + 3*nbs_p

       !! Write tire energy item group
       call writeItGDef (respDB,idEnergy,iFile,'Tire energy')
       call writeVarDef (respDB,header%idEpot,iFile,'Potential energy','ENERGY')
       call writeVarDef (respDB,header%idEkin,iFile,'Kinetic energy','ENERGY')
       call writeVarDef (respDB,header%idEinp,iFile,'Energy input','ENERGY')
       if (iFile > 0) write(iFile,"(']')")

    end if

    write(idatd,"('}')")

    header%firstCall(4) = .false.

    return

900 ierr = idTire
    call reportError (debugFileOnly_p,'writeTireHeader: Tire'//getId(tire%id))

  end subroutine writeTireHeader


  !!============================================================================
  !> @brief Writes results database file header for a general function (engine).
  !>
  !> @param respDB File for secondary response variables
  !> @param header Variable- and item group index container
  !> @param[in] engine The general function to write file header for
  !>
  !> @callgraph @callergraph

  subroutine writeEngineHeader (respDB,header,engine)

    use FunctionTypeModule, only : EngineType
    use IdTypeModule      , only : writeIdHeader
    use RDBModule         , only : idatd, nbs_p, writeVarDef

    type(RDBType)   , intent(inout) :: respDB
    type(HeaderId)  , intent(inout) :: header
    type(EngineType), intent(in)    :: engine

    if (engine%saveVar /= 2) return ! No results requested for this engine

    respDB%nBytes = respDB%nBytes + nbs_p
    call writeIdHeader ('Engine',engine%id,idatd)
    call writeVarDef (respDB,header%idNon,idatd,'Value')
    write(idatd,"('}')")

  end subroutine writeEngineHeader


  !!============================================================================
  !> @brief Writes results database file header for the control system.
  !>
  !> @param respDB File for secondary response variables
  !> @param header Variable- and item group index container
  !> @param[in] ctrl Control system data
  !>
  !> @callgraph @callergraph

  subroutine writeControlHeader (respDB,header,ctrl)

    use ControlTypeModule, only : ControlType
    use IdTypeModule     , only : writeIdHeader
    use RDBModule        , only : idatd, nbs_p, writeVarDef

    type(RDBType)    , intent(inout) :: respDB
    type(HeaderId)   , intent(inout) :: header
    type(ControlType), intent(in)    :: ctrl

    integer :: i

    if (.not. ctrl%saveVar) return ! No control system results requested

    do i = 1, size(ctrl%vregId)
       respDB%nBytes = respDB%nBytes + nbs_p
       call writeIdHeader ('CtrlLine',ctrl%vregId(i),idatd)
       call writeVarDef (respDB,header%idNon,idatd,'Value')
       write(idatd,"('}')")
    end do

  end subroutine writeControlHeader


  !!============================================================================
  !> @brief Writes results database file header for the mechanism.
  !>
  !> @param respDB File for secondary response variables
  !> @param header Variable- and item group index container
  !> @param[in] mech Mechanism components of the model
  !> @param[in] neq Number of equations (system vector dimension)
  !>
  !> @callgraph @callergraph

  subroutine writeMechanismHeader (respDB,header,mech,neq)

    use MechanismTypeModule, only : MechanismType
    use NormTypeModule     , only : nNormTypes_p
    use IdTypeModule       , only : writeIdHeader
    use RDBModule          , only : idatd, nbi_p, nbs_p, nbd_p
    use RDBModule          , only : writeVarDef, writeItgDef

    type(RDBType)      , intent(inout) :: respDB
    type(HeaderId)     , intent(inout) :: header
    type(MechanismType), intent(in)    :: mech
    integer            , intent(in)    :: neq

    integer :: iFile, idSystemEn, idAlgor(0:5), idConv(0:19), idBuoy(2)

    if (.not. any(mech%saveVar)) return ! No mechanism results requested

    !! Note: Assuming this subroutine is invoked only once per file.
    !!       Therefore, the local id-variables do not need to be saved.
    idSystemEn = 0
    idAlgor    = 0
    idConv     = 0
    idBuoy     = 0

    call writeIdHeader ('Mechanism',mech%id,idatd)

    if (mech%saveVar(1)) then

       respDB%nBytes = respDB%nBytes + 3*nbs_p
       call writeVarDef (respDB,header%idCG,idatd, &
            &            'Center of Gravity', &
            &            'LENGTH','VEC3','"x_cg","y_cg","z_cg"')

    end if

    if (mech%saveVar(5)) then

       respDB%nBytes = respDB%nBytes + 4*nbs_p
       call writeVarDef (respDB,idBuoy(1),idatd, &
            &            'Center of Buoyancy', &
            &            'LENGTH','VEC3','"x_cb","y_cb","z_cb"')
       call writeVarDef (respDB,idBuoy(2),idatd,'Buoyancy volume','LENGTH^3')

    end if

    if (mech%saveVar(2)) then

       respDB%nBytes = respDB%nBytes + 5*nbs_p
       call writeItGDef (respDB,idSystemEn   ,iFile,'System energy')
       call writeVarDef (respDB,header%idEstr,iFile,'Strain energy','ENERGY')
       call writeVarDef (respDB,header%idEkin,iFile,'Kinetic energy','ENERGY')
       call writeVarDef (respDB,header%idEpot,iFile,'Potential energy','ENERGY')
       call writeVarDef (respDB,header%idEdmp,iFile,'Energy loss','ENERGY')
       if (mech%hasEinp) then
          respDB%nBytes = respDB%nBytes + nbs_p
          call writeVarDef (respDB,header%idEinp,iFile, &
               &            'Energy input','ENERGY')
       end if
       if (mech%hasEext) then
          respDB%nBytes = respDB%nBytes + nbs_p
          call writeVarDef (respDB,header%idEext,iFile, &
               &           'External energy','ENERGY')
       end if
       call writeVarDef (respDB,header%idEtot,iFile,'Energy check-sum','ENERGY')
       if (iFile > 0) write(iFile,"(']')")

    end if

    if (mech%saveVar(3)) then

       respDB%nBytes = respDB%nBytes + 2*nbi_p + 3*nbs_p
       call writeItGDef (respDB,idAlgor(0),iFile,'Algorithm parameters')
       call writeVarDef (respDB,idAlgor(1),iFile,'Number of iterations', &
            &                                    'NONE','NUMBER')
       call writeVarDef (respDB,idAlgor(2),iFile,'Number of system matrix updates', &
            &                                    'NONE','NUMBER')
       call writeVarDef (respDB,idAlgor(3),iFile,'Time step size','TIME')
       call writeVarDef (respDB,idAlgor(4),iFile,'Wall time','TIME')
       call writeVarDef (respDB,idAlgor(5),iFile,'CPU time','TIME')
       if (iFile > 0) write(iFile,"(']')")

    end if

    if (mech%saveVar(4)) then
       respDB%nBytes = respDB%nBytes + 2*(4*nNormTypes_p+2)*nbs_p

       !! Write the Convergence Item-group with all sub-items
       call writeItGDef (respDB,idConv(0),iFile,'Convergence')
       if (iFile > 0) then

          !! Displacement norms
          call writeItGDefSub ('Displacement correction',1,newLine=.true.)
          call writeItGDefSub ('Scaled vector norm',2,newLine=.true.)
          call writeVarDef (respDB,header%idNon,iFile,'Value')
          call writeVarDef (respDB,idConv(1),iFile,'Tolerance')
          write(iFile,"(']')")
          call writeItGDefSub ('Max translation',2)
          call writeVarDef (respDB,idConv(2),iFile,'Value','LENGTH')
          call writeVarDef (respDB,idConv(3),iFile,'Tolerance','LENGTH')
          write(iFile,"(']')")
          call writeItGDefSub ('Max rotation',2)
          call writeVarDef (respDB,idConv(4),iFile,'Value','ANGLE')
          call writeVarDef (respDB,idConv(5),iFile,'Tolerance','ANGLE')
          write(iFile,"(']')")
          call writeItGDefSub ('Max genDOF',2)
          call writeVarDef (respDB,header%idNon,iFile,'Value')
          call writeVarDef (respDB,idConv(1),iFile,'Tolerance')
          write(iFile,"(']'/'  ]')")

          !! Velocity norms
          call writeItGDefSub ('Velocity correction',1)
          call writeItGDefSub ('Scaled vector norm',2,newLine=.true.)
          call writeVarDef (respDB,header%idNon,iFile,'Value')
          call writeVarDef (respDB,idConv(1),iFile,'Tolerance')
          write(iFile,"(']')")
          call writeItGDefSub ('Max translation',2)
          call writeVarDef (respDB,idConv(6),iFile,'Value','LENGTH/TIME')
          call writeVarDef (respDB,idConv(7),iFile,'Tolerance','LENGTH/TIME')
          write(iFile,"(']')")
          call writeItGDefSub ('Max rotation',2)
          call writeVarDef (respDB,idConv(8),iFile,'Value','ANGLE/TIME')
          call writeVarDef (respDB,idConv(9),iFile,'Tolerance','ANGLE/TIME')
          write(iFile,"(']')")
          call writeItGDefSub ('Max genDOF',2)
          call writeVarDef (respDB,header%idNon,iFile,'Value')
          call writeVarDef (respDB,idConv(1),iFile,'Tolerance')
          write(iFile,"(']'/'  ]')")

          !! Acceleration norms
          call writeItGDefSub ('Acceleration correction',1)
          call writeItGDefSub ('Scaled vector norm',2,newLine=.true.)
          call writeVarDef (respDB,header%idNon,iFile,'Value')
          call writeVarDef (respDB,idConv(1),iFile,'Tolerance')
          write(iFile,"(']')")
          call writeItGDefSub ('Max translation',2)
          call writeVarDef (respDB,idConv(10),iFile,'Value','LENGTH/TIME^2')
          call writeVarDef (respDB,idConv(11),iFile,'Tolerance','LENGTH/TIME^2')
          write(iFile,"(']')")
          call writeItGDefSub ('Max rotation',2)
          call writeVarDef (respDB,idConv(12),iFile,'Value','ANGLE/TIME^2')
          call writeVarDef (respDB,idConv(13),iFile,'Tolerance','ANGLE/TIME^2')
          write(iFile,"(']')")
          call writeItGDefSub ('Max genDOF',2)
          call writeVarDef (respDB,header%idNon,iFile,'Value')
          call writeVarDef (respDB,idConv(1),iFile,'Tolerance')
          write(iFile,"(']'/'  ]')")

          !! Residual norms
          call writeItGDefSub ('Residual Force',1)
          call writeItGDefSub ('Scaled vector norm',2,newLine=.true.)
          call writeVarDef (respDB,header%idNon,iFile,'Value')
          call writeVarDef (respDB,idConv(1),iFile,'Tolerance')
          write(iFile,"(']')")
          call writeItGDefSub ('Max force',2)
          call writeVarDef (respDB,idConv(14),iFile,'Value','FORCE')
          call writeVarDef (respDB,idConv(15),iFile,'Tolerance','FORCE')
          write(iFile,"(']')")
          call writeItGDefSub ('Max torque',2)
          call writeVarDef (respDB,idConv(16),iFile,'Value','FORCE*LENGTH')
          call writeVarDef (respDB,idConv(17),iFile,'Tolerance','FORCE*LENGTH')
          write(iFile,"(']')")
          call writeItGDefSub ('Max generalized force',2)
          call writeVarDef (respDB,header%idNon,iFile,'Value')
          call writeVarDef (respDB,idConv(1),iFile,'Tolerance')
          write(iFile,"(']'/'  ]')")

          !! Energy norms
          call writeItGDefSub ('Iteration Energy',1)
          call writeItGDefSub ('Average DOF energy',2,newLine=.true.)
          call writeVarDef (respDB,idConv(18),iFile,'Value','ENERGY')
          call writeVarDef (respDB,idConv(19),iFile,'Tolerance','ENERGY')
          write(iFile,"(']')")
          call writeItGDefSub ('Max DOF energy',2)
          call writeVarDef (respDB,idConv(18),iFile,'Value','ENERGY')
          call writeVarDef (respDB,idConv(19),iFile,'Tolerance','ENERGY')
          write(iFile,"(']'/'  ]'/']')")

       end if
    end if
    if (mech%saveVar(6)) then
       respDB%nBytes = respDB%nBytes + neq*nbd_p

       write(idatd,601) neq
601    format(/'  <;"Inertia forces";LENGTH;FLOAT;64;VECTOR;(',i6,')>')

    end if

    write(idatd,"('}')")

  contains

    !> @brief Writes an indented item group definition.
    subroutine writeItGDefSub (cId,indent,newLine)
      character(len=*), intent(in) :: cId
      integer         , intent(in) :: indent
      logical,optional, intent(in) :: newLine
      character(len=4), parameter  :: blank = '    '
      if (present(newLine)) write(iFile,"()")
      write(iFile,600,advance='NO') blank(1:2*indent),cId
600   format(a,'[;"',a,'";')
    end subroutine writeItGDefSub

  end subroutine writeMechanismHeader


  !!============================================================================
  !> @brief Writes results database file header for a wind turbine.
  !>
  !> @param respDB File for secondary response variables
  !> @param[in] turbine Wind turbine object to write file header for
  !>
  !> @callgraph @callergraph

  subroutine writeWindTurbineHeader (respDB,turbine)

    use WindTurbineTypeModule, only : TurbineConfig
    use IdTypeModule         , only : writeIdHeader
    use RDBModule            , only : idatd, nbs_p
    use RDBModule            , only : writeVarDef, writeItgDef

    type(RDBType)      , intent(inout) :: respDB
    type(TurbineConfig), intent(in)    :: turbine

    integer :: i, iFile, idWP, idWS(2)

    call writeIdHeader ('Turbine',turbine%id,idatd)

    respDB%nBytes = respDB%nBytes &
         &        + 6*nbs_p*(1 + size(turbine%blade) + size(turbine%windPt,2))

    idWP = 0
    idWS = 0
    call writeItGDef (respDB,idWP,iFile,'Hub')
    call writeVarDef (respDB,idWS(1),iFile,'Wind speed',&
         &            'LENGTH/TIME','VEC3','"v_x","v_y","v_z"')
    call writeVarDef (respDB,idWS(2),iFile,'Undisturbed wind speed', &
         &            'LENGTH/TIME','VEC3','"v_x","v_y","v_z"')
    write(iFile,"(']')")

    do i = 1, size(turbine%blade)
       idWP = 0
       call writeItGDef (respDB,idWP,iFile,'Blade tip '//CHAR(ICHAR('0')+i))
       call writeVarDef (respDB,idWS(1),iFile,'WS1')
       call writeVarDef (respDB,idWS(2),iFile,'WS2')
       write(iFile,"(']')")
    end do

    do i = 1, size(turbine%windPt,2)
       idWP = 0
       call writeItGDef (respDB,idWP,iFile, &
            &            'User-defined point '//CHAR(ICHAR('0')+i))
       call writeVarDef (respDB,idWS(1),iFile,'WS1')
       call writeVarDef (respDB,idWS(2),iFile,'WS2')
       write(iFile,"(']')")
    end do

    write(idatd,"('}')")

  end subroutine writeWindTurbineHeader

  !!============================================================================
  !> @brief Writes primary variables for a triad to results database.
  !>
  !> @param respDB File for primary response variables
  !> @param[in] triad The triad to write result variables for
  !> @param[in] writeAsDouble If .true., write variables in double precision
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph

  subroutine writeTriadDB1 (respDB,triad,writeAsDouble,ierr)

    use TriadTypeModule  , only : TriadType
    use RDBModule        , only : writeRDB
    use reportErrorModule, only : reportError, debugFileOnly_p

    type(RDBType)  , intent(inout) :: respDB
    type(TriadType), intent(in)    :: triad
    logical        , intent(in)    :: writeAsDouble
    integer        , intent(out)   :: ierr

    ierr = 0
    if (triad%id%baseId < 0) return ! Dummy CG triad for generic parts
    if (.not. triad%savePos) return ! Output is switched off for this triad

    if (triad%nDOFs == 6) then
       call writeRDB (respDB,triad%ur,ierr,writeAsDouble)
    else if (triad%nDOFs == 3) then
       call writeRDB (respDB,triad%ur(:,4),ierr,writeAsDouble)
    end if

    if (ierr < 0) call reportError (debugFileOnly_p,'writeTriadDB1')

  end subroutine writeTriadDB1


  !!============================================================================
  !> @brief Writes secondary variables for a triad to results database.
  !>
  !> @param respDB File for secondary response variables
  !> @param[in] triad The triad to write result variables for
  !> @param[in] writeAsDouble If .true., write variables in double precision
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph

  subroutine writeTriadDB2 (respDB,triad,writeAsDouble,ierr)

    use TriadTypeModule    , only : TriadType, transVGlobToTriad, dp
    use TriadTypeModule    , only : transVSysToGlob, transVSysToTriad
    use HydroDynamicsModule, only : getCalculatedFluidMotion
    use RDBModule          , only : writeRDB
    use reportErrorModule  , only : reportError, debugFileOnly_p

    type(RDBType)  , intent(inout) :: respDB
    type(TriadType), intent(in)    :: triad
    logical        , intent(in)    :: writeAsDouble
    integer        , intent(out)   :: ierr

    real(dp) :: z, fvel(3), facc(3)

    ierr = 0
    if (triad%nDOFs < 1)     return ! Earth-linked triad, no results
    if (triad%id%baseId < 0) return ! Dummy CG triad for generic parts

    if (triad%saveVar(7) .and. associated(triad%ur_def)) then

       !! Global deformations (relative to superelement rigid motion)
       call writeRDB (respDB,triad%ur_def,ierr,writeAsDouble)

    end if
    if (triad%saveVar(1)) then

       !! Global velocities
       call writeRDB (respDB,triad%urd,ierr,writeAsDouble)

    end if
    if (triad%saveVar(2)) then

       !! Global accelerations
       call writeRDB (respDB,triad%urdd,ierr,writeAsDouble)

    end if
    if (triad%saveVar(3)) then
       if (associated(triad%scoord)) then

          !! Curvilinear coordinate along 1D structures
          call writeRDB (respDB,triad%scoord,ierr)
          !! Curvature along 1D structures
          call writeRDB (respDB,triad%kappa,ierr)
          !! Associated bending moments
          call writeRDB (respDB,triad%Moment,ierr)

       end if
       if (associated(triad%supForce)) then

          if (triad%isFEnode) then
             !! Local section force
             call writeRDB (respDB,triad%supForce,ierr,writeAsDouble)
          else
             !! Global superelement force
             call writeRDB (respDB,transVSysToGlob(triad,triad%supForce),ierr, &
                  &         writeAsDouble)
          end if

       end if
       if (associated(triad%nodeForce)) then

          !! Total system forces
          call writeRDB (respDB,triad%nodeForce,ierr,writeAsDouble)

       end if
    end if
    if (triad%saveVar(4)) then

       !! Local velocities
       call writeRDB (respDB,transVGlobToTriad(triad,triad%urd),ierr, &
            &         writeAsDouble)

    end if
    if (triad%saveVar(5)) then

       !! Local accelerations
       call writeRDB (respDB,transVGlobToTriad(triad,triad%urdd),ierr, &
            &         writeAsDouble)

    end if
    if (triad%saveVar(6) .and. associated(triad%supForce)) then
       if (.not. triad%isFEnode) then ! not relevant for section forces

          !! Local superelement force
          call writeRDB (respDB,transVSysToTriad(triad,triad%supForce),ierr, &
               &         writeAsDouble)

       end if
    end if
    if (triad%saveVar(6) .and. associated(triad%nodeForce)) then

       !! Local total system forces
       call writeRDB (respDB,transVSysToTriad(triad,triad%nodeForce),ierr, &
            &         writeAsDouble)

    end if
    if (triad%saveVar(9)) then
       if (getCalculatedFluidMotion(triad,z,fvel,facc)) then

          !! Fluid particle motion
          call writeRDB (respDB,z,ierr)
          call writeRDB (respDB,fvel,ierr)
          call writeRDB (respDB,facc,ierr)

       end if
    end if
    if (triad%saveVar(6) .and. associated(triad%aeroForce)) then

       !! Aerodynamic forces
       call writeRDB (respDB,triad%aeroForce,ierr)

    end if
    if (ierr < 0) call reportError (debugFileOnly_p,'writeTriadDB2')

  end subroutine writeTriadDB2


  !!============================================================================
  !> @brief Writes some force variables for a triad to results database.
  !>
  !> @param respDB File for secondary response variables
  !> @param[in] triad The triad to write result variables for
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph

  subroutine writeTriadForcesDB (respDB,triad,ierr)

    use TriadTypeModule  , only : TriadType
    use RDBModule        , only : writeRDB
    use reportErrorModule, only : reportError, debugFileOnly_p

    type(RDBType)   , intent(inout) :: respDB
    type(TriadType) , intent(in)    :: triad
    integer         , intent(out)   :: ierr

    ierr = 0

    if (associated(triad%force)) then
       call writeRDB (respDB,triad%force,ierr)
    end if

    if (ierr < 0) call reportError (debugFileOnly_p,'writeTriadForcesDB')

  end subroutine writeTriadForcesDB


  !!============================================================================
  !> @brief Writes primary variables for a superelement to results database.
  !>
  !> @param respDB File for primary response variables
  !> @param[in] supel The superelement to write result variables for
  !> @param[in] writeAsDouble If .true., write variables in double precision
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph

  subroutine writeSupelDB1 (respDB,supel,writeAsDouble,ierr)

    use SupElTypeModule  , only : SupElType
    use RDBModule        , only : writeRDB
    use reportErrorModule, only : reportError, debugFileOnly_p

    type(RDBType)  , intent(inout) :: respDB
    type(SupElType), intent(in)    :: supel
    logical        , intent(in)    :: writeAsDouble
    integer        , intent(out)   :: ierr

    ierr = 0
    if (.not. supel%savePos) return ! Output is switched off

    call writeRDB (respDB,supel%supTR,ierr,writeAsDouble)

    if (associated(supel%genDOFs)) then
       call writeRDB (respDB,supel%genDOFs%ur,ierr,writeAsDouble)
    end if

    if (ierr < 0) call reportError (debugFileOnly_p,'writeSupelDB1')

  end subroutine writeSupelDB1


  !!============================================================================
  !> @brief Writes secondary variables for a superelement to results database.
  !>
  !> @param respDB File for secondary response variables
  !> @param[in] supel The superelement to write result variables for
  !> @param[in] masses All point masses in the model
  !> @param sysMass Total mass for the mechanism
  !> @param sysCG Centre of gravity for the mechanism
  !> @param sysVb Buoyancy volume for the mechanism
  !> @param sysCB Centre of buoyancy for the mechanism
  !> @param[in] writeAsDouble If .true., write variables in double precision
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph

  subroutine writeSupelDB2 (respDB,supel,masses,sysMass,sysCG,sysVb,sysCB, &
       &                    writeAsDouble,ierr)

    use SupElTypeModule    , only : SupElType, IsBeam
    use MassTypeModule     , only : MassType
    use KindModule         , only : dp, epsDiv0_p
    use FiniteElementModule, only : GetSectionalForces
    use manipMatrixModule  , only : matmul34
    use RDBModule          , only : writeRDB
    use reportErrorModule  , only : reportError, debugFileOnly_p

    type(RDBType)  , intent(inout) :: respDB
    type(SupElType), intent(in)    :: supel
    type(MassType) , intent(in)    :: masses(:)
    real(dp)       , intent(inout) :: sysMass(3), sysCG(3), sysVb, sysCB(3)
    logical        , intent(in)    :: writeAsDouble
    integer        , intent(out)   :: ierr

    real(dp) :: mass(3), vec(3)
    integer  :: i, j, id

    ierr = 0

    if (supel%saveVar(2) .or. supel%saveVar(4)) then
       if (associated(supel%scoord)) then
          !! Curvilinear coordinate along 1D structures
          call writeRDB (respDB,supel%scoord,ierr)
       end if
    end if

    !! Compute center of gravity in global coordinates
    vec = matmul34(supel%suptr,supel%posMassCenter)
    if (supel%saveVar(1)) then
       call writeRDB (respDB,vec,ierr)
    end if

    mass = supel%mass
    if (supel%addedMass) then
       !! Account for additional masses on triads
       vec = vec*mass
       do i = 1, size(masses)
          id = masses(i)%triad%id%baseId
          do j = 1, supel%nExtNods
             if (supel%triads(j)%p%id%baseId == id) then
                mass = mass + masses(i)%mass
                vec = vec + supel%triads(j)%p%ur(:,4)*masses(i)%mass
             end if
          end do
       end do
       do i = 1, 3
          if (mass(i) > epsDiv0_p) vec(i) = vec(i) / mass(i)
       end do
       if (supel%saveVar(1)) then
          call writeRDB (respDB,vec,ierr)
       end if
    end if

    !! Accumulate the system mass and CG
    sysMass = sysMass + mass
    sysCG   = sysCG + vec*mass

    if ( supel%triads(1)%p%id%baseId < 0 .and. &
         supel%triads(1)%p%nDOFs > 0 .and. supel%saveVar(1)) then

       !! Results at the center of gravity dummy triad
       call writeRDB (respDB,supel%triads(1)%p%ur,ierr,writeAsDouble)
       call writeRDB (respDB,supel%triads(1)%p%urd,ierr,writeAsDouble)
       call writeRDB (respDB,supel%triads(1)%p%urdd,ierr,writeAsDouble)

    end if
    if (supel%saveVar(2) .and. IsBeam(supel)) then

       !! Sectional forces at the beam end points
       call writeRDB (respDB,GetSectionalForces(supel,1),ierr)
       call writeRDB (respDB,GetSectionalForces(supel,2),ierr)

    end if
    if (supel%saveVar(4) .and. associated(supel%hydyn)) then

       !! Hydrodynamic quantities
       if (supel%hydyn%Vb(1) > 0.0_dp) then
          vec   = matmul34(supel%suptr,supel%hydyn%C0b(:,1))
          sysVb = sysVb + supel%hydyn%Vb(1)
          sysCB = sysCB + supel%hydyn%Vb(1)*vec
          call writeRDB (respDB,vec,ierr)
       else ! no bouyancy volume, C0b should then be identically zero
          mass = 0.0_dp
          call writeRDB (respDB,mass,ierr)
       end if
       call writeRDB (respDB,matmul(supel%suptr(:,1:3),supel%hydyn%B(1:3)),ierr)
       if (IsBeam(supel)) then
          call writeRDB (respDB,matmul(supel%suptr(:,1:3),supel%hydyn%M(1:3)), &
               &         ierr)
          call writeRDB (respDB,matmul(supel%suptr(:,1:3),supel%hydyn%D(1:3)), &
               &         ierr)
          call writeRDB (respDB,matmul(supel%suptr(:,1:3),supel%hydyn%S(1:3)), &
               &         ierr)
       else
          call writeRDB (respDB,matmul(supel%suptr(:,1:3),supel%hydyn%B(4:6)), &
               &         ierr)
          call writeRDB (respDB,matmul(supel%suptr(:,1:3),supel%hydyn%D(1:3)), &
               &         ierr)
          call writeRDB (respDB,matmul(supel%suptr(:,1:3),supel%hydyn%D(4:6)), &
               &         ierr)
          if (any(abs(supel%hydyn%S(1:3)) > 1.0e-15_dp)) then
             call writeRDB (respDB,matmul34(supel%suptr,supel%hydyn%C0s),ierr)
          else ! no slam force, C0s should then be identically zero
             mass = 0.0_dp
             call writeRDB (respDB,mass,ierr)
          end if
          call writeRDB (respDB,matmul(supel%suptr(:,1:3),supel%hydyn%S(1:3)), &
               &         ierr)
          call writeRDB (respDB,matmul(supel%suptr(:,1:3),supel%hydyn%S(4:6)), &
               &         ierr)
       end if

    end if
    if (associated(supel%genDOFs) .and. supel%saveVar(2)) then

       do i = 1, supel%genDOFs%nDOFs
          vec(1) = supel%genDOFs%ur(i)
          vec(2) = supel%genDOFs%urd(i)
          vec(3) = supel%genDOFs%urdd(i)
          call writeRDB (respDB,vec,ierr,writeAsDouble)
          if (supel%saveVar(3)) then
             call writeRDB (respDB,supel%genDOFs%energy(:,i),ierr)
          end if
       end do

    end if
    if (supel%saveVar(3)) then

       call writeRDB (respDB,supel%Estr,ierr)
       call writeRDB (respDB,supel%Ekin,ierr)
       call writeRDB (respDB,supel%Epot,ierr)
       call writeRDB (respDB,supel%Edmp,ierr)

    end if

    if (ierr < 0) call reportError (debugFileOnly_p,'writeSupelDB2')

  end subroutine writeSupelDB2


  !!============================================================================
  !> @brief Writes results for a user-defined element to results database file.
  !>
  !> @param respDB File for secondary response variables
  !> @param[in] elm The user-defined element to write result variables for
  !> @param[in] writeAsDouble If .true., write variables in double precision
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph

  subroutine writeUserdefElDB (respDB,elm,writeAsDouble,ierr)

    use UserdefElTypeModule, only : UserdefElType, dp
    use RDBModule          , only : writeRDB
    use FiUserElmInterface , only : Fi_UDE5
    use reportErrorModule  , only : reportError, debugFileOnly_p

    type(RDBType)      , intent(inout) :: respDB
    type(UserdefElType), intent(in)    :: elm
    logical            , intent(in)    :: writeAsDouble
    integer            , intent(out)   :: ierr

    integer  :: i, nvar
    real(dp) :: value

    ierr = 0

    call writeRDB (respDB,elm%Tlg,ierr,writeAsDouble)

    do i = 1, max_UDE_p

       call Fi_UDE5 (elm%id%baseId,elm%type,i, &
            &        elm%iwork(1),elm%rwork(1),value,nvar)
       if (nvar < i .or. ierr < 0) exit

       call writeRDB (respDB,value,ierr)

    end do

    if (ierr < 0) call reportError (debugFileOnly_p,'writeUserdefElDB')

  end subroutine writeUserdefElDB


  !!============================================================================
  !> @brief Writes results for a spring to results database file.
  !>
  !> @param respDB File for secondary response variables
  !> @param[in] spring The spring to write result variables for
  !> @param[in] writeAsDouble If .true., write variables in double precision
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph

  subroutine writeSpringDB (respDB,spring,writeAsDouble,ierr)

    use SpringTypeModule , only : SpringBaseType, dp
    use RDBModule        , only : writeRDB
    use reportErrorModule, only : reportError, debugFileOnly_p

    type(RDBType)       , intent(inout) :: respDB
    type(SpringBaseType), intent(in)    :: spring
    integer             , intent(out)   :: ierr
    logical             , intent(in)    :: writeAsDouble

    logical  :: haveYield
    real(dp) :: elasticDeflection

    ierr = 0

    haveYield = spring%unLoadType == 1 .or. spring%unLoadType == 2 .or. &
         &      associated(spring%yield)

    if (spring%saveVar(1)) then

       call writeRDB (respDB,spring%stiffness,ierr,writeAsDouble)

    end if
    if (spring%saveVar(1) .and. spring%unLoadType == 3) then

       call writeRDB (respDB,spring%s0,ierr,writeAsDouble)

    end if
    if (spring%saveVar(2) .and. associated(spring%length0Engine)) then

       call writeRDB (respDB,spring%length0,ierr,writeAsDouble)

    end if
    if (spring%saveVar(2)) then

       call writeRDB (respDB,spring%length,ierr,writeAsDouble)

    end if
    if (spring%saveVar(3)) then

       call writeRDB (respDB,spring%deflection,ierr,writeAsDouble)

    end if
    if (spring%saveVar(3) .and. haveYield) then

       call writeRDB (respDB,spring%yieldDeflection,ierr,writeAsDouble)

       if (spring%hasFailed) then
          elasticDeflection = 0.0_dp
       else
          elasticDeflection = spring%deflection - spring%yieldDeflection
       end if
       call writeRDB (respDB,elasticDeflection,ierr,writeAsDouble)

    end if
    if (spring%saveVar(4)) then

       call writeRDB (respDB,spring%force,ierr,writeAsDouble)

    end if
    if (spring%saveVar(5)) then

       call writeRDB (respDB,spring%Estr,ierr)

    end if
    if (spring%saveVar(5) .and. associated(spring%length0Engine)) then

       call writeRDB (respDB,spring%Einp,ierr)

    end if
    if (spring%saveVar(5) .and. haveYield) then

       call writeRDB (respDB,spring%Eyield,ierr)

    end if

    if (ierr < 0) call reportError (debugFileOnly_p,'writeSpringDB')

  end subroutine writeSpringDB


  !!============================================================================
  !> @brief Writes results for a damper to results database file.
  !>
  !> @param respDB File for secondary response variables
  !> @param[in] damper The damper to write result variables for
  !> @param[in] writeAsDouble If .true., write variables in double precision
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph

  subroutine writeDamperDB (respDB,damper,writeAsDouble,ierr)

    use DamperTypeModule , only : DamperBaseType
    use RDBModule        , only : writeRDB
    use reportErrorModule, only : reportError, debugFileOnly_p

    type(RDBType)       , intent(inout) :: respDB
    type(DamperBaseType), intent(in)    :: damper
    logical             , intent(in)    :: writeAsDouble
    integer             , intent(out)   :: ierr

    ierr = 0

    if (damper%saveVar(1)) then

       call writeRDB (respDB,damper%coeff,ierr,writeAsDouble)

    end if
    if (damper%saveVar(2)) then

       call writeRDB (respDB,damper%length,ierr,writeAsDouble)

    end if
    if (damper%saveVar(3)) then

       call writeRDB (respDB,damper%velocity,ierr,writeAsDouble)

    end if
    if (damper%saveVar(4)) then

       call writeRDB (respDB,damper%force,ierr,writeAsDouble)

    end if
    if (damper%saveVar(5)) then

       call writeRDB (respDB,damper%Edmp,ierr)

    end if

    if (ierr < 0) call reportError (debugFileOnly_p,'writeDamperDB')

  end subroutine writeDamperDB


  !!============================================================================
  !> @brief Writes results for a friction to results database file.
  !>
  !> @param respDB File for secondary response variables
  !> @param[in] friction The friction to write result variables for
  !> @param[in] writeAsDouble If .true., write variables in double precision
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph

  subroutine writeFrictionDB (respDB,friction,writeAsDouble,ierr)

    use FrictionTypeModule, only : FrictionType, ROT_FRICTION_p
    use KindModule        , only : dp, epsDiv0_p
    use RDBModule         , only : writeRDB
    use reportErrorModule , only : reportError, debugFileOnly_p

    type(RDBType)     , intent(inout) :: respDB
    type(FrictionType), intent(in)    :: friction
    logical           , intent(in)    :: writeAsDouble
    integer           , intent(out)   :: ierr

    real(dp) :: R

    ierr = 0

    if (friction%saveVar(1)) then

       call writeRDB (respDB,friction%force,ierr,writeAsDouble)
       R = friction%param%typeDepParams(1) ! Bearing radius
       if (friction%param%type == ROT_FRICTION_p .and. R > epsDiv0_p) then
          call writeRDB (respDB,friction%Fequ/R,ierr)
       else
          call writeRDB (respDB,friction%Fequ,ierr)
       end if
       if (friction%Fequ > epsDiv0_p) then
          call writeRDB (respDB,friction%force/friction%Fequ,ierr)
       else
          call writeRDB (respDB,0.0_dp,ierr)
       end if

    end if
    if (friction%saveVar(3)) then

       call writeRDB (respDB,friction%Fmax,ierr)
       call writeRDB (respDB,friction%pos,ierr)
       call writeRDB (respDB,friction%vel,ierr)

    end if
    if (friction%saveVar(2)) then

       call writeRDB (respDB,friction%Edmp,ierr)

    end if
    if (friction%saveVar(3)) then

       call writeRDB (respDB,friction%lSlip(1)*10+friction%lSlip(2),ierr)

    end if

    if (ierr < 0) call reportError (debugFileOnly_p,'writeFrictionDB')

  end subroutine writeFrictionDB


  !!============================================================================
  !> @brief Writes results for an external force to results database file.
  !>
  !> @param respDB File for secondary response variables
  !> @param[in] force The external force to write result variables for
  !> @param[in] writeAsDouble If .true., write variables in double precision
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph

  subroutine writeForceDB (respDB,force,writeAsDouble,ierr)

    use ForceTypeModule  , only : ForceType, dp
    use TriadTypeModule  , only : transSysToGlob
    use RDBModule        , only : writeRDB
    use reportErrorModule, only : reportError, debugFileOnly_p

    type(RDBType)  , intent(inout) :: respDB
    type(ForceType), intent(in)    :: force
    logical        , intent(in)    :: writeAsDouble
    integer        , intent(out)   :: ierr

    real(dp) :: forceVec(3)

    ierr = 0
    if (force%dof > 0) return ! Joint and triad DOF forces are not saved

    if (force%saveVar(1)) then

       forceVec = force%F * force%forceDir
       if (associated(force%triad)) call transSysToGlob (force%triad,forceVec)
       call writeRDB (respDB,forceVec,ierr,writeAsDouble)

    end if
    if (force%saveVar(2)) then

       call writeRDB (respDB,force%F,ierr,writeAsDouble)

    end if
    if (force%saveVar(3)) then

       call writeRDB (respDB,force%Einp,ierr)

    end if

    if (ierr < 0) call reportError (debugFileOnly_p,'writeForceDB')

  end subroutine writeForceDB


  !!============================================================================
  !> @brief Writes results for a joint to results database file.
  !>
  !> @param respDB File for secondary response variables
  !> @param[in] joint The joint to write result variables for
  !> @param[in] forces All external forces in the model
  !> @param[in] motions All prescribed motions in the model
  !> @param[in] writeAsDouble If .true., write variables in double precision
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph

  subroutine writeJointDB (respDB,joint,forces,motions,writeAsDouble,ierr)

    use MasterSlaveJointTypeModule, only : MasterSlaveJointType, JointDofType
    use ForceTypeModule           , only : ForceType
    use MotionTypeModule          , only : MotionType
    use RDBModule                 , only : writeRDB
    use reportErrorModule         , only : reportError, debugFileOnly_p

    type(RDBType)             , intent(inout) :: respDB
    type(MasterSlaveJointType), intent(in)    :: joint
    type(ForceType)           , intent(in)    :: forces(:)
    type(MotionType)          , intent(in)    :: motions(:)
    logical                   , intent(in)    :: writeAsDouble
    integer                   , intent(out)   :: ierr

    integer :: iDof, jDof, nJDofs

    ierr = 0
    nJDofs = size(joint%jointDofs)
    if (nJDofs < 1) return

    ! Write the joint variables
    do jDof = 1, nJDofs
       iDof = joint%dofOutputOrder(jDof)
       call writeJointVarDB (respDB,joint%jointDofs(iDof),ierr)
       if (ierr < 0) goto 500
    end do

    ! Write the joint springs
    do jDof = 1, nJDofs
       iDof = joint%dofOutputOrder(jDof)
       if (associated(joint%jointDofs(iDof)%spring)) then
          call writeSpringDB (respDB,joint%jointDofs(iDof)%spring, &
               &              writeAsDouble,ierr)
          if (ierr < 0) goto 500
       end if
    end do

    ! Write the joint dampers
    do jDof = 1, nJDofs
       iDof = joint%dofOutputOrder(jDof)
       if (associated(joint%jointDofs(iDof)%damper)) then
          call writeDamperDB (respDB,joint%jointDofs(iDof)%damper, &
               &              writeAsDouble,ierr)
          if (ierr < 0) goto 500
       end if
    end do

    ! Write the joint friction
    do jDof = 1, nJDofs
       iDof = joint%dofOutputOrder(jDof)
       if (associated(joint%jointDofs(iDof)%friction)) then
          call writeFrictionDB (respDB,joint%jointDofs(iDof)%friction, &
               &                writeAsDouble,ierr)
          if (ierr < 0) goto 500
       end if
    end do

    return

500 call reportError (debugFileOnly_p,'writeJointDB')

  contains

    !> @brief Writes a joint variable to the results database file.
    subroutine writeJointVarDB (respDB,jointDof,ierr)
      type(RDBType)     , intent(inout) :: respDB
      type(JointDofType), intent(in)    :: jointDof
      integer           , intent(out)   :: ierr

      integer :: i

      ierr = 0

      if (all(jointDof%saveVar)) then
         call writeRDB (respDB,jointDof%jvar,ierr,writeAsDouble)
      else
         do i = 1, size(jointDof%jvar)
            if (jointDof%saveVar(i)) then
               call writeRDB (respDB,jointDof%jvar(i),ierr,writeAsDouble)
            end if
         end do
      end if

      i = jointDof%loadIdx
      if (i > 0 .and. i <= size(forces)) then
         if (forces(i)%saveVar(2)) then

            call writeRDB (respDB,forces(i)%F,ierr,writeAsDouble)

         end if
         if (forces(i)%saveVar(3)) then

            call writeRDB (respDB,forces(i)%Einp,ierr)

         end if
      else if (i < 0 .and. -i <= size(motions)) then
         i = -i
         if (motions(i)%saveVar(2)) then

            call writeRDB (respDB,motions(i)%F,ierr,writeAsDouble)

         end if
         if (motions(i)%saveVar(3)) then

            call writeRDB (respDB,motions(i)%Einp,ierr)

         end if
      else if (associated(jointDof%F)) then
         if (jointDof%saveVar(4)) then

            call writeRDB (respDB,jointDof%F,ierr,writeAsDouble)

         end if
      end if

    if (ierr < 0) call reportError (debugFileOnly_p,'writeJointVarDB')

    end subroutine writeJointVarDB

  end subroutine writeJointDB


  !!============================================================================
  !> @brief Writes results for a contact element to results database file.
  !>
  !> @param respDB File for secondary response variables
  !> @param[in] cElem The contact element to write result variables for
  !> @param[in] writeAsDouble If .true., write variables in double precision
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph

  subroutine writeContactDB (respDB,cElem,writeAsDouble,ierr)

    use ContactElementTypeModule, only : ContactElementType
    use RDBModule               , only : writeRDB
    use reportErrorModule       , only : reportError, debugFileOnly_p

    type(RDBType)           , intent(inout) :: respDB
    type(ContactElementType), intent(in)    :: cElem
    logical                 , intent(in)    :: writeAsDouble
    integer                 , intent(out)   :: ierr

    integer :: j, iDof

    ierr = 0

    ! Write the contact variables
    do iDof = 1, 3
       do j = 1, 2
          if (cElem%saveVar(j)) then
             call writeRDB (respDB,cElem%cVar(iDof,j),ierr,writeAsDouble)
             if (ierr < 0) goto 500
          end if
       end do
    end do

    ! Write accumulated contact distance
    if (cElem%saveVar(1)) then
       call writeRDB (respDB,cElem%accContactDist,ierr)
       if (ierr < 0) goto 500
    end if

    ! Write the contact springs
    do iDof = 1, size(cElem%springs)
       if (associated(cElem%springs(iDof)%p)) then
          call writeSpringDB (respDB,cElem%springs(iDof)%p,writeAsDouble,ierr)
          if (ierr < 0) goto 500
       end if
    end do

    ! Write the contact dampers
    do iDof = 1, size(cElem%dampers)
       if (associated(cElem%dampers(iDof)%p)) then
          call writeDamperDB (respDB,cElem%dampers(iDof)%p,writeAsDouble,ierr)
          if (ierr < 0) goto 500
       end if
    end do

    ! Write the contact friction
    if (associated(cElem%friction)) then
       call writeFrictionDB (respDB,cElem%friction,writeAsDouble,ierr)
       if (ierr < 0) goto 500
    end if

    return

500 call reportError (debugFileOnly_p,'writeContactDB')

  end subroutine writeContactDB


  !!============================================================================
  !> @brief Writes results for a tire to results database file.
  !>
  !> @param respDB File for secondary response variables
  !> @param[in] tire The tire to write result variables for
  !> @param[in] scaleTo Scaling factors for tire model variables
  !> @param sysMass Total mass for the mechanism
  !> @param sysCG Centre of gravity for the mechanism
  !> @param[in] writeAsDouble If .true., write variables in double precision
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph

  subroutine writeTireDB (respDB,tire,scaleTo,sysMass,sysCG,writeAsDouble,ierr)

    use TireTypeModule    , only : TireType, ScaleToType, STI_p, CTI_p, MF_p
    use KindModule        , only : dp, epsDiv0_p
    use RDBModule         , only : writeRDB
    use scratchArrayModule, only : getRealScratchArray
    use reportErrorModule , only : internalError, reportError, debugFileOnly_p

    type(RDBType)    , intent(inout) :: respDB
    type(TireType)   , intent(in)    :: tire
    type(ScaleToType), intent(in)    :: scaleTo
    real(dp)         , intent(inout) :: sysMass(3), sysCG(3)
    logical          , intent(in)    :: writeAsDouble
    integer          , intent(out)   :: ierr

    real(dp), pointer :: varinf(:)
    real(dp)          :: tirePos(3)

    ! Pointers to additional output array from CTIO (from ftp2.H by M.Gipser)
    ! ERRC     1 error code (0:last integration step is ok)
    ! IKA      2 longitudinal slip [%]
    ! IAL      3 side slip angle [deg]
    ! IGA      4 tire camber angle [deg]
    ! IFFC     5 (3) forces in foot-print, in cont. frame [N]
    ! ITFC     8 (3) torques in foot-print, in cont. frame [Nm]
    ! IFF0    11 (3) forces in foot-print, in inert. f. [N]
    ! ITF0    14 (3) torques in foot-print, in inert. f. [Nm]
    ! IFRC    17 (3) forces on rim, in cont. frame [N]
    ! ITRC    20 (3) torques on rim, in cont. frame [Nm]
    ! IFRI    23 (3) forces on rim, in ISO-axis syst. [N]
    ! ITRI    26 (3) torques on rim, in ISO-axis syst. [Nm]
    ! IFFW    29 (3) forces in foot-print, in wheel f. [N]
    ! ITFW    32 (3) torques in foot-print, in wheel f. [Nm]
    ! IFRW    35 (3) forces on rim, in wheel frame [N]
    ! ITRW    38 (3) torques on rim, in wheel frame [Nm]
    ! IRR     41 mean distance rim center - road [mm]
    ! IHFP    42 mean absolute height foot-print [m]
    ! ILFP    43 length of foot-print [mm]
    ! BEXO    44 actual relative belt extension [%]
    ! PNMO    45 maximum ground pressure in foot-print [bar]
    ! PNMEO   46 mean ground pressure in foot-print [bar]
    ! ICU     47 mean longit. long-waved curv. of road prof. [1/m]
    ! ITD     48 global tire deflection [mm]
    ! ITDD    49 global tire deflection velocity [m/s]
    ! IDRM    50 maximum radial belt element displacement [mm]
    ! IRCM    51 (3) approx. foot-print center [m]
    ! TGX     54 gyroscopic overturning moment [Nm]
    ! TGZ     55 gyroscopic aligning moment [Nm]
    ! IVRX    56 longitudinal rim center velocity [m/s]
    ! IVRY    57 lateral rim center velocity [m/s]
    ! IVRZ    58 vertical rim center velocity [m/s]
    ! IVSCPX  59 longitudinal slip velocity at contact point [m/s]
    ! IVSCPY  60 lateral slip velocity at contact point [m/s]
    ! IWSCPZ  61 contact patch bore velocity [rad/s]
    ! IVOL    62 actual interior volume [m^3]
    ! IPR     63 actual interior pressure [bar]
    ! IFP     64 foot-print area cg, expr. in contact frame [mm]
    ! IXFP    65 x-shift of belt above foot-print [mm]
    ! IYFP    66 y-shift of belt above foot-print [mm]
    ! ITFP    67 torsion of belt above foot-print [deg]
    ! IBFP    68 bending of belt above foot-print [1/m]
    ! ISFP    69 long. geometrical shift of foot-print [mm]
    ! FFRC    70 flag whether rim contacts / penetrates road [-]
    ! IA0W    71 (9) transf. m. to wheel frame
    ! IA0C    80 (9) transf. m. to road cont. frame (X-axis system)
    ! IA0I    89 (9) transf. m. to ISO 8855 axis system
    ! ILON    98 longitudinal rim displacement [mm]
    ! ILAT    99 lateral rim displacement [mm]
    ! IFTSD  100 force/moment standard deviations [N],[Nm]

    ierr = 0
    if (tire%api == STI_P) then
       varinf => tire%sti%varinf
    else if (tire%api == CTI_p) then
       varinf => getRealScratchArray(nVarCTI_p,ierr)
#ifdef FT_HAS_CTI
       call CTIO (tire%idIn,nVarCTI_p,varinf)
#endif
    else
       ierr = internalError('writeTireDB: Invalid tire API')
       return
    end if

    if (tire%saveVar(1)) then

       ! Slip angle, longitudinal slip and camber angle
       if (tire%api == STI_p) then
          call writeRDB (respDB,varinf(1:3),ierr)
       else if (tire%api == CTI_p) then
          call writeRDB (respDB,varinf(7:9),ierr)
       end if

    end if
    if (tire%saveVar(2)) then

       ! Effective rolling radius
       if (tire%api == STI_p) then
          call writeRDB (respDB,varinf(12)/scaleTo%m,ierr)
       else if (tire%api == CTI_p) then
          call writeRDB (respDB,varinf(49)/scaleTo%m,ierr)
       end if

    end if
    if (tire%saveVar(3)) then

       ! Contact forces
       if (tire%api == STI_p) then
          call writeRDB (respDB,varinf(4:6)/scaleTo%N,ierr,writeAsDouble)
          call writeRDB (respDB,varinf(7:9)/(scaleTo%N*scaleTo%m),ierr, &
               &         writeAsDouble)
       else if (tire%api == CTI_p) then
          call writeRDB (respDB,varinf(1:3)/scaleTo%N,ierr,writeAsDouble)
          call writeRDB (respDB,varinf(4:6)/(scaleTo%N*scaleTo%m),ierr, &
               &         writeAsDouble)
       end if

    end if
    if (tire%saveVar(4)) then

       ! Contact point position and road normal
       if (tire%api == STI_p) then
          call writeRDB (respDB,varinf(13:15)/scaleTo%m,ierr)
          call writeRDB (respDB,varinf(16:18),ierr)
       else if (tire%api == CTI_p) then
          call writeRDB (respDB,varinf(26:28)/scaleTo%m,ierr)
          call writeRDB (respDB,varinf(31),ierr)
          call writeRDB (respDB,varinf(34),ierr)
          call writeRDB (respDB,varinf(37),ierr)
       end if

    end if
#ifdef FT_DEBUG
    if (any(tire%saveVar(3:4))) then

       call writeRDB (respDB,tire%isInContact,ierr)
       call writeRDB (respDB,tire%isSliding,ierr)

    end if
#endif
    if (tire%saveVar(5)) then

       ! Tire deflections (state variables)
       if (tire%api == STI_p) then
          call writeRDB (respDB,tire%sti%deqvar(1:2)/scaleTo%m,ierr)
          call writeRDB (respDB,varinf(10)/scaleTo%m,ierr)
       else if (tire%api == CTI_p) then
          call writeRDB (respDB,varinf(44)/scaleTo%m,ierr)
       end if

    end if
    if (tire%saveVar(6) .and. tire%api == STI_p) then

       ! Deflection velocities
       call writeRDB (respDB,tire%sti%deqder(1:2)*scaleTo%s/scaleTo%m,ierr)
       call writeRDB (respDB,varinf(11)*scaleTo%s/scaleTo%m,ierr)

    end if
    if (tire%saveVar(7) .and. tire%type == MF_p) then

       ! Tire characteristics
       call writeRDB (respDB,varinf(19:22),ierr)
       call writeRDB (respDB,varinf(24:35),ierr)

    end if
#ifdef FT_DEBUG
    if (tire%saveVar(7) .and. tire%api == CTI_p) then

       ! Dump out whole varinf "as is"
       call writeRDB (respDB,varinf(10:nVarCTI_p),ierr)

    end if
#endif
    if (tire%saveVar(8)) then

       ! Wheel Carrier force and torque
       call writeRDB (respDB,tire%forceInG,ierr,writeAsDouble)
       call writeRDB (respDB,tire%forceInWC,ierr,writeAsDouble)

    end if
    if (tire%saveVar(9)) then

       ! Tire energy
       call writeRDB (respDB,tire%ePot,ierr)
       call writeRDB (respDB,tire%eKin,ierr)
       call writeRDB (respDB,tire%eInp,ierr)

    end if

    if (tire%mass > epsDiv0_p) then
       ! Accumulate system mass and center of gravity
       sysMass = sysMass + tire%mass
       tirePos = tire%RimTriad%ur(:,4) + tire%ZoffSet*tire%RimTriad%ur(:,3)
       sysCG   = sysCG + tire%mass*tirePos
    end if

    if (ierr < 0) call reportError (debugFileOnly_p,'writeTireDB')

  end subroutine writeTireDB


  !!============================================================================
  !> @brief Writes results for a general function to results database file.
  !>
  !> @param respDB File for secondary response variables
  !> @param[in] engine The general function to write result variables for
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph

  subroutine writeEngineDB (respDB,engine,ierr)

    use FunctionTypeModule, only : EngineType
    use RDBModule         , only : writeRDB
    use reportErrorModule , only : reportError, debugFileOnly_p

    type(RDBType)    , intent(inout) :: respDB
    type(EngineType) , intent(in)    :: engine
    integer          , intent(out)   :: ierr

    !! --- Logic section ---

    ierr = 0

    if (engine%saveVar == 2) then

       call writeRDB (respDB,engine%value,ierr)

    end if

    if (ierr < 0) call reportError (debugFileOnly_p,'writeEngineDB')

  end subroutine writeEngineDB


  !!============================================================================
  !> @brief Writes results for the control system to results database file.
  !>
  !> @param respDB File for secondary response variables
  !> @param[in] ctrl Control system data
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph

  subroutine writeControlDB (respDB,ctrl,ierr)

    use ControlTypeModule, only : ControlType
    use RDBModule        , only : writeRDB
    use reportErrorModule, only : reportError, debugFileOnly_p

    type(RDBType)    , intent(inout) :: respDB
    type(ControlType), intent(in)    :: ctrl
    integer          , intent(out)   :: ierr

    integer :: i

    ierr = 0

    if (ctrl%saveVar) then

       call writeRDB (respDB,ctrl%vreg,ierr)

       !! Bugfix #323: Write control variable for all lines referring to it
       do i = 1, size(ctrl%ivar)
          call writeRDB (respDB,ctrl%vreg(ctrl%ivar(i)),ierr)
       end do

    end if

    if (ierr < 0) call reportError (debugFileOnly_p,'writeControlDB')

  end subroutine writeControlDB


  !!============================================================================
  !> @brief Writes results for the mechanism to results database file.
  !>
  !> @param respDB File for secondary response variables
  !> @param[in] mech Mechanism components of the model
  !> @param[in] sys System level model data
  !> @param[in] posCG Centre of gravity for the mechanism
  !> @param[in] posCB Centre of buoyancy for the mechanism
  !> @param[in] sysVb Buoyancy volume for the mechanism
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph

  subroutine writeMechanismDB (respDB,mech,sys,posCG,posCB,sysVB,ierr)

    use MechanismTypeModule   , only : MechanismType
    use SystemTypeModule      , only : SystemType
    use TimerModule           , only : getCurrentTime
    use NormTypeModule        , only : nNormTypes_p
    use KindModule            , only : sp, dp
    use RDBModule             , only : writeRDB
    use reportErrorModule     , only : reportError, debugFileOnly_p
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_intValue

    type(RDBType)      , intent(inout) :: respDB
    type(MechanismType), intent(in)    :: mech
    type(SystemType)   , intent(in)    :: sys
    real(dp)           , intent(in)    :: posCG(3), posCB(3), sysVB
    integer            , intent(out)   :: ierr

    integer        :: i
    real(sp)       :: deltaTime(2)
    real(sp), save :: totalTime(2)
    logical , save :: lFirst = .true.

    ierr = 0

    if (mech%saveVar(1)) then

       call writeRDB (respDB,posCG,ierr)

    end if
    if (mech%saveVar(5)) then

       call writeRDB (respDB,posCB,ierr)
       call writeRDB (respDB,sysVB,ierr)

    end if
    if (mech%saveVar(2)) then

       call writeRDB (respDB,mech%Estr,ierr)
       call writeRDB (respDB,mech%Ekin,ierr)
       call writeRDB (respDB,mech%Epot,ierr)
       call writeRDB (respDB,mech%Edmp,ierr)
       if (mech%hasEinp) call writeRDB (respDB,mech%Einp,ierr)
       if (mech%hasEext) call writeRDB (respDB,mech%Eext,ierr)
       call writeRDB (respDB,mech%Etot,ierr)

    end if
    if (mech%saveVar(3)) then

       call writeRDB (respDB,sys%nIterThisStep,ierr)
       call writeRDB (respDB,sys%nUpdaThisStep,ierr)
       call writeRDB (respDB,sys%timeStep,ierr)

       if (lFirst) then
          totalTime = 0.0
          lFirst = .false.
       end if

       deltaTime = getCurrentTime(1) - totalTime
       totalTime = totalTime + deltaTime
       call writeRDB (respDB,deltaTime,ierr)

    end if
    if (mech%saveVar(4)) then

       do i = 1, nNormTypes_p
          call writeRDB (respDB,sys%convergenceSet%disNorms(i)%value,ierr)
          call writeRDB (respDB,sys%convergenceSet%disNorms(i)%tolerance,ierr)
       end do
       do i = 1, nNormTypes_p
          call writeRDB (respDB,sys%convergenceSet%velNorms(i)%value,ierr)
          call writeRDB (respDB,sys%convergenceSet%velNorms(i)%tolerance,ierr)
       end do
       do i = 1, nNormTypes_p
          call writeRDB (respDB,sys%convergenceSet%accNorms(i)%value,ierr)
          call writeRDB (respDB,sys%convergenceSet%accNorms(i)%tolerance,ierr)
       end do
       do i = 1, nNormTypes_p
          call writeRDB (respDB,sys%convergenceSet%resNorms(i)%value,ierr)
          call writeRDB (respDB,sys%convergenceSet%resNorms(i)%tolerance,ierr)
       end do
       call writeRDB (respDB,sys%convergenceSet%energyNorms(1)%value,ierr)
       call writeRDB (respDB,sys%convergenceSet%energyNorms(1)%tolerance,ierr)
       call writeRDB (respDB,sys%convergenceSet%energyNorms(2)%value,ierr)
       call writeRDB (respDB,sys%convergenceSet%energyNorms(2)%tolerance,ierr)

    end if
    if (mech%saveVar(6)) then

       if (mod(ffa_cmdlinearg_intValue('NewmarkFlag'),10) > 0) then
          i = size(sys%FIactual)
          !! Store FIactual = Qk - FDk - FSk
          call DCOPY (i,sys%Qk(1),1,sys%FIactual(1),1)
          call DAXPY (i,-1.0_dp,sys%FDk(1),1,sys%FIactual(1),1)
          call DAXPY (i,-1.0_dp,sys%FSk(1),1,sys%FIactual(1),1)
          call writeRDB (respDB,sys%FIactual,ierr,.true.)
       else ! Store FIactual = FIk
          call writeRDB (respDB,sys%FIk,ierr,.true.)
       end if

    end if
    if (ierr < 0) call reportError (debugFileOnly_p,'writeMechanismDB')

  end subroutine writeMechanismDB

end module saveModule
