!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module ExtCtrlSysRoutinesModule

  use PointerKindModule, only : ptr

  implicit none

  ! Define variables residing in workspace
  ! --------------------------------------
  integer(ptr), save, private :: ep, timeSpan, x, simin


contains

  subroutine InitExtCtrlSys (tstart,dt,extCtrlSys,ierr)

    !!==========================================================================
    !! Initialise the external control systems before time integration.
    !!
    !! Programmer : Knut Morten Okstad                 date/rev : Jan 2003 / 3.0
    !!==========================================================================

    use ExtCtrlSysTypeModule, only : ExtCtrlSysType, MATLAB_p, dp
    use profilerModule      , only : startTimer, stopTimer, ctrl_p
    use reportErrorModule   , only : internalError

    real(dp)            , intent(in)    :: tstart, dt
    type(ExtCtrlSysType), intent(inout) :: extCtrlSys(:)
    integer             , intent(out)   :: ierr

    !! Local variables
    integer :: elm

    !! --- Logic section ---

    call startTimer (ctrl_p)

    do elm = 1, size(extCtrlSys)
       select case (extCtrlSys(elm)%sysType)

       case (MATLAB_p)
          call InitMatlabCtrlSys (tstart,dt,extCtrlSys(elm),ierr)

       case default
          ierr = internalError('Invalid external control system type')
          return

       end select
    end do

    call stopTimer (ctrl_p)

  end subroutine InitExtCtrlSys


  subroutine IterateExtCtrlSys (mode,tstart,tstop,extCtrlSys,ierr)

    !!==========================================================================
    !! Iterate the external control systems during time integration.
    !!
    !! Programmer : Knut Morten Okstad                 date/rev : Jan 2003 / 3.0
    !!==========================================================================

    use ExtCtrlSysTypeModule, only : ExtCtrlSysType, MATLAB_p, dp
    use profilerModule      , only : startTimer, stopTimer, ctrl_p
    use reportErrorModule   , only : internalError

    integer             , intent(in)    :: mode
    real(dp)            , intent(in)    :: tstart, tstop
    type(ExtCtrlSysType), intent(inout) :: extCtrlSys(:)
    integer             , intent(out)   :: ierr

    !! Local variables
    integer :: elm

    !! --- Logic section ---

    call startTimer (ctrl_p)

    do elm = 1, size(extCtrlSys)
       select case (extCtrlSys(elm)%sysType)

       case (MATLAB_p)
          call IterateMatlabCtrlSys (mode,tstart,tstop,extCtrlSys(elm),ierr)

       case default
          ierr = internalError('Invalid external control system type')
          return

       end select
    end do

    call stopTimer (ctrl_p)

  end subroutine IterateExtCtrlSys


  subroutine InitMatlabCtrlSys (tstart,dt,ctrlSys,ierr)

    !!==========================================================================
    !! Initialise an external control system defined in simulink before time
    !! integration, i.e., start Matlab session, create workspace variables and
    !! establish the command to be used in the time integration.
    !!
    !! Programmer : Haavar Johan Johnsen                date/rev: Dec 2001 / 1.0
    !! Revised    : Torleif Iversen                               Jun 2002 / 2.0
    !! Revised    : Knut Morten Okstad                            Jan 2003 / 3.0
    !! Revised    : Torleif Iversen                               May 2003 / 3.1
    !!==========================================================================

    use ExtCtrlSysTypeModule   , only : ExtCtrlSysType
    use kindModule             , only : dp, lfnam_p
    use reportErrorModule      , only : getErrorFile, reportError, error_p
    use reportErrorModule      , only : internalError, allocationError
    use EngineRoutinesModule   , only : EngineValue
    use FiMatlabPluginInterface, only : m_p_init, m_p_engOpen, m_p_engEvalString
    use FiMatlabPluginInterface, only : m_p_mxCreateDoubleMatrix
    use FiMatlabPluginInterface, only : m_p_mxCopyReal8ToPtr
    use FiMatlabPluginInterface, only : m_p_mxCopyPtrToReal8, m_p_mxGetData
    use FiMatlabPluginInterface, only : m_p_engPutVariable, m_p_engGetVariable
    use FiMatlabPluginInterface, only : m_p_mxGetM, m_p_mxGetN

    real(dp)            , intent(in)    :: tstart, dt
    type(ExtCtrlSysType), intent(inout) :: ctrlSys
    integer             , intent(out)   :: ierr

    !! Local variables
    real(dp) :: wsInput(2), wsTime(2) !! 2*1 array to be sent to Matlab
    real(dp) :: fixedStep
    integer  :: i, strPos, dim
    integer(ptr) :: wsPtr
    character(len=lfnam_p) :: path, model
#if defined(win32) || defined(win64)
    character(len=2), parameter :: pathsep = '\\'
#else
    character(len=1), parameter :: pathsep = '/'
#endif

    !! --- Logic section ---

    ! Check that the Matlab libraries are loaded
    call m_p_init(ierr)
    if (ierr /= 0) then
       call reportError (error_p,'Matlab libraries have not been loaded.', &
            &            'Ensure that Matlab is installed on the computer', &
            &            'before running simulations with external ctrl blocks')
       return
    end if

    ! Start the Matlab engine, making ep a Matlab workspace
    ep = m_p_engOpen('matlab')
    if (ep == 0) then
       ierr = internalError('Failed to start start Matlab engine')
       return
    end if

    ! Clear the workspace in case the Matlab session is already running
    if (m_p_engEvalString(ep,'clear') /= 0) then
       call ReportError (error_p,'Clearing Matlab workspace failed')
       ierr = -1
       return
    end if

    ! Create start and stop time in Matlab workspace
    wsTime(1) = tstart
    wsTime(2) = tstart + 0.01 ! Magnitude of increment is without significance
    timeSpan  = m_p_mxCreateDoubleMatrix(1,2,0)
    call m_p_mxCopyReal8ToPtr(wsTime(1),m_p_mxGetData(timeSpan),2)
    if (m_p_engPutVariable(ep,'timeSpan',timeSpan) /= 0) then
       call ReportError (error_p,'Creating timeSpan in Matlab failed')
       ierr = -1
       return
    end if

    ! Create a 1*2 array in Matlab workspace for input to the Simulink model
    simin = m_p_mxCreateDoubleMatrix(1,2,0)

    ! Read Simulink inputs from Fedem engine
    do i = 1, size(ctrlSys%engineIn)
       wsInput(1) = tstart
       wsInput(2) = EngineValue(ctrlSys%engineIn(i)%p,ierr)
       if (ierr < 0) return
       ! Copy the name in "match" to Matlab workspace and transfer value
       call m_p_mxCopyReal8ToPtr(wsInput(1),m_p_mxGetData(simin),2)
       if (m_p_engPutVariable(ep,trim(ctrlSys%match(i)),simin) /= 0) then
          call ReportError (error_p,'Failed to copy name '// &
               &            trim(ctrlSys%match(i))//' to Matlab')
          ierr = -1
          return
       end if
    end do

    ! Break up the model path into path and modelname
    strPos = index(ctrlSys%simCmd,'.mdl')
    if (strPos <= 0) then
       call ReportError (error_p,'File type is not Simulink: '//ctrlSys%simCmd)
       ierr = -1
       return
    end if
    model  = ctrlSys%simCmd(1:strPos-1)
    strPos = scan(model,pathsep,BACK=.TRUE.)
    path   = ctrlSys%simCmd(1:strPos-1)
    model  = model(strPos+1:len(model))
    if (index(trim(path),' ') > 0) then
       path = '''' // trim(adjustl(path)) // ''''
    end if
    if (3+len_trim(path) > len(ctrlSys%path)) then
       ierr = internalError('Simulink model path is too long')
       return
    end if
    ctrlSys%path = 'cd '//path(1:len_trim(path))

#ifdef FT_DEBUG
    i = getErrorFile()
    write(i,*) 'InitExtCtrlSys: path  = ',trim(path)
    write(i,*) 'InitExtCtrlSys: model = ',trim(model),'.mdl'
#endif

    ! Change to the model directory
    if (m_p_engEvalString(ep,trim(ctrlSys%path)) /= 0) then
       call ReportError (error_p,'Failure changing directory: '//trim(path))
       ierr = -1
    end if

    ! Get integration algorithm and time step settings
    ctrlSys%simCmd = 'alg = simget(''' // trim(model) // ''',''Solver'')'
    if (m_p_engEvalString(ep,trim(ctrlSys%simCmd)) /= 0) then
       call ReportError (error_p,'Getting algorithm of ext ctrl-system failed')
       ierr = -1
       return
    end if

    wsPtr = m_p_engGetVariable(ep,'alg')
    if (m_p_mxGetN(wsPtr) > 8) then
       call ReportError (error_p,'A continuous state solver must be used in Simulink')
       ierr = -1
       return
    end if

    if (m_p_mxGetN(wsPtr) > 4) then
       ctrlSys%timeDim = -1
    else
       ctrlSys%simCmd = 'timeStep = simget(''' // trim(model) // ''',''FixedStep'')'
       if (m_p_engEvalString(ep,trim(ctrlSys%simCmd)) /= 0) then
          call ReportError (error_p,'Getting time step of ext ctrl-system failed')
          ierr = -1
          return
       end if
       wsPtr = m_p_engGetVariable(ep,'timeStep')
       call m_p_mxCopyPtrToReal8(m_p_mxGetData(wsPtr), fixedStep, 1)
       ctrlSys%timeDim = int(dt/fixedStep + 0.000001) + 1
    end if

    ! Run the trim or sim command to initialise the model

    ctrlSys%simCmd = 'x = trim(''' // trim(model) // ''')'
    if (m_p_engEvalString(ep,trim(ctrlSys%simCmd)) /= 0) then
       call ReportError (error_p,'Initialisation of ext ctrl-system failed')
       ierr = -1
       return
    end if
    wsPtr = m_p_engGetVariable(ep,'x')

    if (wsPtr == 0) then
#ifdef FT_DEBUG
       write(i,*) 'InitExtCtrlSys: Internal state vector dimension is zero'
#endif
       ctrlSys%simCmd = '[tt,x] = sim(''' // trim(model) // ''')'
       if (m_p_engEvalString(ep,trim(ctrlSys%simCmd)) /= 0) then
          call ReportError (error_p,'Initialisation of ext ctrl-system failed')
          ierr = -1
          return
       end if
       wsPtr = m_p_engGetVariable(ep,'x')
       if (wsPtr <= 0) then
          call ReportError (error_p,'Failed to create internal state vector')
          ierr = -1
          return
       end if
    end if

    ! Create and populate the initial internal state of the model
    x = m_p_mxCreateDoubleMatrix(m_p_mxGetM(wsPtr),m_p_mxGetN(wsPtr),0)
    dim = m_p_mxGetM(wsPtr)*m_p_mxGetN(wsPtr)

    ! Build the sim command to be used in the time integration
    ctrlSys%simCmd = '[tt,x] = sim(''' // trim(model) // ''',' // &
         &           'timeSpan,' // &
         &           'simset(''InitialState'',x,' // &
         &           '''FinalStateName'',''x''))'

    ! Run a time integration step to get the initial output values
    if (m_p_engEvalString(ep,trim(ctrlSys%simCmd)) /= 0) then
       call ReportError (error_p,'Initialisation of ext ctrl-system failed')
       ierr = -1
       return
    end if

    ! Save only the last values of the internal state
    if (dim > 0) then
       if (m_p_engEvalString(ep,'x = x(size(x,1),:)') /= 0) then
          call ReportError (error_p,'Getting internal state of ext ctrl-system failed')
          ierr = -1
          return
       end if
    end if

    ! Create and populate the initial internal state of the model
    wsPtr = m_p_engGetVariable(ep, 'x')
    if (wsPtr <= 0) then
       call ReportError (error_p,'No initial state vector in workspace')
       ierr = -1
       return
    end if
    x = m_p_mxCreateDoubleMatrix(m_p_mxGetM(wsPtr),m_p_mxGetN(wsPtr),0)

    ! Allocate space for the internal Simulink variables
    dim = m_p_mxGetM(wsPtr)*m_p_mxGetN(wsPtr)
    allocate(ctrlSys%state(dim),STAT=ierr)
    if (ierr /= 0) then
       ierr = AllocationError('Array of internal Simulink states')
       return
    end if

    ! Save internal Simulink variables in the state array
    if (dim > 0) then
       call m_p_mxCopyPtrToReal8(m_p_mxGetData(wsPtr),ctrlSys%state(1),dim)
    end if

    ! Transfer the Simulink model outputs to the Fedem sensors
    do i = 1, size(ctrlSys%sensorOut)
       call transfer2FedemSensor (ctrlSys%sensorOut(i)%p,ierr)
    end do

  end subroutine InitMatlabCtrlSys


  subroutine IterateMatlabCtrlSys (mode,tstart,tstop,ctrlSys,ierr)

    !!==========================================================================
    !! Returning the output from the first run.
    !! Running for all but the first iteration.
    !!
    !! Programmer : Haavar Johan Johnsen                date/rev: Dec 2001 / 1.0
    !! Revised    : Torleif Iversen                               Jun 2002 / 2.0
    !! Revised    : Knut Morten Okstad                            Jan 2003 / 3.0
    !! Revised    : Torleif Iversen                               May 2003 / 3.1
    !!==========================================================================

    use ExtCtrlSysTypeModule   , only : ExtCtrlSysType, dp
    use reportErrorModule      , only : internalError, getErrorFile
    use reportErrorModule      , only : reportError, error_p, warningFileOnly_p
    use EngineRoutinesModule   , only : EngineValue
    use FiMatlabPluginInterface, only : m_p_engEvalString
    use FiMatlabPluginInterface, only : m_p_mxCopyReal8ToPtr
    use FiMatlabPluginInterface, only : m_p_mxCopyPtrToReal8, m_p_mxGetData
    use FiMatlabPluginInterface, only : m_p_engPutVariable, m_p_engGetVariable
    use FiMatlabPluginInterface, only : m_p_mxGetM, m_p_mxGetN

    integer             , intent(in)    :: mode
    real(dp)            , intent(in)    :: tstart, tstop
    type(ExtCtrlSysType), intent(inout) :: ctrlSys
    integer             , intent(out)   :: ierr

    !! Local variables
    integer      :: i, dimXws, dimState, nIntVars
    integer(ptr) :: wsPtr, ttPtr
    real(dp)     :: wsInput(2), wsTime(2) !! 2*1 array to be sent to Matlab

    !! --- Logic section ---

    ierr = 0

    ! Copy start and stop time to Matlab workspace
    wsTime(1) = tstart
    wsTime(2) = tstop
    call m_p_mxCopyReal8ToPtr(wsTime(1),m_p_mxGetData(timeSpan),2)
    if (m_p_engPutVariable(ep,'timeSpan',timeSpan) /= 0) then
       call ReportError (error_p,'Creating timeSpan in Matlab failed')
       ierr = -1
       return
    end if

    ! Read Simulink inputs from Fedem engine
    do i = 1, size(ctrlSys%engineIn)
       wsInput(1) = tstop
       wsInput(2) = EngineValue(ctrlSys%engineIn(i)%p,ierr)
       if (ierr < 0) return
       ! Copy the name in "match" to Matlab workspace and transfer value
       call m_p_mxCopyReal8ToPtr(wsInput(1),m_p_mxGetData(simin),2)
       if (m_p_engPutVariable(ep,trim(ctrlSys%match(i)),simin) /= 0) then
          call ReportError (error_p,'Transfer of '//trim(ctrlSys%match(i))// &
               &                    ' value failed')
          ierr = -1
          return
       end if
    end do

    ! Change to the model directory
    if (m_p_engEvalString(ep,trim(ctrlSys%path)) /= 0) then
       call ReportError (error_p,trim(ctrlSys%path)//' failed')
       ierr = -1
       return
    end if

    ! At first iteration get internal state from vector x in Workspace
    ! If the size of x has changed, re-allocate the saved internal state,
    ! unless the change in dimension is caused by a truncation of x

    dimState = size(ctrlSys%state)
    if (dimState > 0) then
       wsPtr = m_p_engGetVariable(ep, 'x')
       nIntVars = m_p_mxGetN(wsPtr)
       dimXws = nIntVars*m_p_mxGetM(wsPtr)

       if (dimXws /= dimState) then
          call ReportError (error_p,'Incorrect state dimensions in workspace')
          ierr = -1
          i = getErrorFile()
          write(i,"(10X,'t =',1PE12.5,' mode =',I3,' dimension',I8,' ->',I8)") &
                   &     tstop, mode, dimState, dimXws
          return
       end if

       if (mode == 4) then

          ! check whether the state dimension has been changed
          if (ctrlSys%timeDim > 0) then
             ttPtr = m_p_engGetVariable(ep, 'tt')
             if (m_p_mxGetM(ttPtr) < ctrlSys%timeDim) then
                call ReportError (warningFileOnly_p,'Simulink has cut off internal state', &
                   'The results may not be quite reliable', &
                   'Consider shorter Simulink time steps or a variable-step algorithm')
             end if
          end if

          ! save internal state in workspace in ctrlSys%state
          call m_p_mxCopyPtrToReal8(m_p_mxGetData(wsPtr),ctrlSys%state(1),dimXws)

       else if (mode == 5) then

          ! At succeeding iterations, load internal state into vector x in Workspace
          call m_p_mxCopyReal8ToPtr(ctrlSys%state(1), &
               &                    m_p_mxGetData(wsPtr),size(ctrlSys%state))

          if (m_p_engPutVariable(ep,'x',wsPtr) /= 0) then
             call ReportError (error_p,'Copy internal state to Matlab failed')
             ierr = -1
             return
          end if

       else
          ierr = internalError('IterateMatlabCtrlSys: Illegal mode value')
          return
       end if

    end if

    ! Run the Simulink model to update the state
    if (m_p_engEvalString(ep,trim(ctrlSys%simCmd)) /= 0) then
       call ReportError (error_p,'Failure executing simulink command '// &
            &            ctrlSys%simCmd)
       ierr = -1
       return
    end if

    ! Save only the last values of the internal state
    if (dimState > 0) then
       if (m_p_engEvalString(ep,'x = x(size(x,1),:)') /= 0) then
          call ReportError (error_p,'Getting internal state of ext ctrl-system failed')
          ierr = -1
          return
       end if
    end if

    ! Transfer the Simulink model outputs to the Fedem sensors
    do i = 1, size(ctrlSys%sensorOut)
       call transfer2FedemSensor (ctrlSys%sensorOut(i)%p,ierr)
    end do

  end subroutine IterateMatlabCtrlSys


  !!============================================================================
  !> @brief Transfers Simulink model output to a Fedem sensor.

  subroutine transfer2FedemSensor (sensor,ierr)

    use SensorTypeModule       , only : SensorType
    use reportErrorModule      , only : reportError, error_p
    use FiMatlabPluginInterface, only : m_p_engGetVariable, m_p_mxCopyPtrToReal8
    use FiMatlabPluginInterface, only : m_p_mxGetM, m_p_mxGetData

    type(SensorType), intent(inout) :: sensor
    integer         , intent(inout) :: ierr

    !! Local variables
    integer(ptr) :: wsPtr

    !! --- Logic section ---

    wsPtr = m_p_engGetVariable(ep,trim(sensor%id%descr))

    if (wsPtr <= 0_ptr) then
       call ReportError (error_p,'Found no sensor '//sensor%id%descr, &
            &            addString='transfer2FedemSensor')
       ierr = ierr - 1
    else if (m_p_mxGetM(wsPtr) > 0) then
       call m_p_mxCopyPtrToReal8 (m_p_mxGetData(wsPtr),sensor%value,1)
    else
       call ReportError (error_p,'No output values obtained for sensor '// &
            &            sensor%id%descr, &
            &            'Check that Output options is set to Refine output', &
            &            addString='transfer2FedemSensor')
       ierr = ierr - 1
    end if

  end subroutine transfer2FedemSensor

end module ExtCtrlSysRoutinesModule
