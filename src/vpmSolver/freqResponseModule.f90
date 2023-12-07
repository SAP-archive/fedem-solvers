!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file freqResponseModule.f90
!>
!> @brief Frequency domain implementation for the FEDEM Dynamics Solver.
!>
!> @details Subroutines of this file are related to frequency domain analysis.
!>
!> @author Guenter Glanzer, SAP SE
!>
!> @date Feb 2019

!!==============================================================================
!> @brief Module with subroutines for frequency domain analysis.

module FreqResponseModule

  use kindModule   , only : dp, pi_p
  use pyplot_module, only : pyplot

  implicit none

  private

  !> @brief Data type for frequency-domain load/motions.
  type dloadMotion
     integer :: motion   !< 0=load, 1=displacement, 2=velocity, 3=acceleration
     logical :: id_fft   !< If .true., input is given in frequency domain (FFT)
     integer :: triad_id !< Triad id where load/motion acts
     integer :: dof      !< Triad's local DOF where load/motion acts (1 to 6)
     real(dp)   , allocatable :: x(:)    !< Abscissa values
     real(dp)   , allocatable :: y(:)    !< Ordinate values
     complex(dp), allocatable :: data(:) !< Complex-valued data collector
  end type dloadMotion

  !> All frequency-domain loads and motions
  type(dloadMotion), pointer, save :: plm(:) => null()

  integer , save :: nrModes        !< number of modes used in the calculation
  integer , save :: windowSize     !< describes the window size in samples
  integer , save :: sweep_dof      !< sweep analysis (sweep location for acting force)
  real(dp), save :: sweep_range(2) !< sweep range (start and end frequency)
  real(dp), save :: fs             !< sampling frequency

  integer, save :: output_triad_dof(20) !< postprocessing, contains triad id and local dof number
  integer, save :: out_dof(10)          !< internal equation number (based on output_triad_dof)

  type(pyplot), pointer, save :: plt => null() !< python plot handler

  integer, save :: ccw = 0         !< total number of windows
  integer, save :: seg_start = 0   !< segstart (total position)

  public :: solveFreqDomain, deallocateFreq


contains

  !!============================================================================
  !> @brief Interface to the frequency-domain solver module in Fedem.
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph

  subroutine solveFreqDomain (ierr)

    use SolverModule          , only : sam, sys, mech
    use SolverModule          , only : closeAll
    use SolverRoutinesModule  , only : restoreLastStep
    use ProgressModule        , only : writeProgress, lterm
    use ReportErrorModule     , only : allocationError, getErrorFile
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_isTrue
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getdouble

    integer, intent(out) :: ierr

    !! Local variables
    real(dp)       :: inc, nextTime
    real(dp), save :: lastFDA = -1.0_dp
    logical , save :: inp_read = .false.      !< flag for reading input file
    logical , save :: run_eigenvalue = .true. !< flag for eigenvalue solution
    logical , save :: run_sweep = .true.      !< flag for sweep analysis


    !! --- Logic section ---

    ierr = 0
    if (ffa_cmdlinearg_isTrue('frequency_domain')) then
       call ffa_cmdlinearg_getdouble ('eiginc',inc)
       if (inc < 0.1_dp*sys%minInc) then
          call ffa_cmdlinearg_getdouble ('timeInc',inc)
       end if
    else
       inc = -1.0_dp
    end if
    if (inc < 0.0_dp) return ! No frequency domain analysis in this execution

#ifdef FT_HAS_FFTPACK
    if (lastFDA < 0.0_dp) lastFDA = -inc
#else
    write(*,*) '*** solveFreqDomain: Built without FFTPACK, frequency domain analysis unavailable.'
    ierr = -999
    call closeAll (ierr,'solveFreqDomain')
    return
#endif

    !! Check if we have reached the time for the next frequency domain analysis
    nextTime = sys%time + 0.1_dp*sys%minInc
    if (nextTime < lastFDA+inc .or. nextTime > sys%tEnd) return ! Not yet

    lastFDA = sys%time

    ! reading parameters for analysis (only once)
    if (.not. inp_read) then
       if (ffa_cmdlinearg_isTrue('pyplot')) then
          allocate(plt,STAT=ierr)
          if (ierr /= 0) ierr = allocationError('solveFreqDomain')
       end if
       call read_input (mech)
       if (nrModes > sam%neq) nrModes = sam%neq
       inp_read = .true.
    end if

    ! perform an eigenvalue analysis (only once)
    if (run_eigenvalue .and. nrModes > 0 .and. ierr == 0) then
       call writeProgress (' --> EIGENVALUE ANALYSIS')
       call freq_eigenvalues (sam%neq,lterm,getErrorFile(),ierr)
       run_eigenvalue = .false.
    end if

    ! perform a sweep analysis (only once)
    if (run_sweep .and. sweep_dof > 0 .and. ierr == 0) then
       call writeProgress (' --> SWEEP ANALYSIS')
       call freq_sweep (sam%neq,ierr)
       run_sweep = .false.
    end if

    ! perform the frequency domain analysis
    if (ierr == 0) then
       call writeProgress (' --> FREQUENCY-DOMAIN ANALYSIS')
       call freq_seg_analysis (sam,sys,mech,inc,lterm,getErrorFile(),ierr)
       call restoreLastStep (sam,sys,mech)
    end if

    if (ierr /= 0) then
       call deAllocateFreq ()
       call closeAll (ierr,'solveFreqDomain')
    end if

  end subroutine solveFreqDomain


  !!============================================================================
  !> @brief Deallocates the dynamic structures of the frequency domain analysis.
  !>
  !> @callergraph

  subroutine deallocateFreq ()

    !! Local variables
    integer :: i


    !! --- Logic section ---

    if (associated(plm)) then
       do i = 1, size(plm)
          deallocate(plm(i)%x, plm(i)%y, plm(i)%data)
       end do
       deallocate(plm)
       nullify(plm)
    end if

    if (associated(plt)) then
       deallocate(plt)
       nullify(plt)
    end if

  end subroutine deallocateFreq


  !!============================================================================
  !> @brief Reads input specifications for the frequency domain analysis.
  !> @param[in] mech Mechanism components of the model
  !>
  !> @callgraph @callergraph

  subroutine read_input (mech)

    use MechanismTypeModule   , only : MechanismType
    use SolverModule          , only : objectEquations
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getint
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getints
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getdouble
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getdoubles

    type(MechanismType), intent(in) :: mech

    !! Local variables
    integer :: i, j, k     !< loop variables
    integer :: triad       !< triad base id
    integer :: ndof        !< number of degrees of freedom (triad related)
    integer :: meqn(6)     !< global equation numbers at the triad
    integer :: inputDof(2) !< integer pair (triad_id/dof_id for sweep force def)


    !! --- Logic section ---

    call ffa_cmdlinearg_getint ('nrModes',nrModes)
    call ffa_cmdlinearg_getint ('windowSize',windowSize)
    call ffa_cmdlinearg_getdouble ('sample_freq',fs)
    call ffa_cmdlinearg_getdoubles ('sweep_range',sweep_range,2)
    call ffa_cmdlinearg_getints ('sweep_input',inputDof,2)
    if (associated(plt)) then ! needed only with -pyplot specified
       call ffa_cmdlinearg_getints ('freq_output',output_triad_dof,20)
    else
       output_triad_dof = 0
    end if

    ! initialisation
    out_dof   = 0
    sweep_dof = 0

    ! search corresponding internal equation number:
    k = 1
    do i = 1, size(mech%triads)
      ndof = mech%triads(i)%nDOFs
      triad = mech%triads(i)%id%baseId
      do j = 1, size(output_triad_dof), 2
        if(triad == output_triad_dof(j)) then
          call objectEquations(triad, meqn, ndof)
          out_dof(k) = meqn(output_triad_dof(j+1))
          k = k+1
        end if
      end do

      if (triad == inputDof(1)) then
        call objectEquations(triad, meqn, ndof)
        sweep_dof = meqn(inputDof(2))
        if (sweep_range(1) >= sweep_range(2)) then
          ! reset values to default
          sweep_range(1) =  0.0_dp
          sweep_range(2) = 50.0_dp
        end if
      end if
    end do

    if (k == 1 .and. associated(plt)) then
      ! output triad(s) not existing, deactivate plotting (always)
      deallocate(plt)
      nullify(plt)
    end if

  end subroutine read_input


  !!============================================================================
  !> @brief Calculates the eigenvalues and eigenvectors (optional).
  !> @param[in]  ndim Dimension of the global equation system
  !> @param[in]  IOC  File unit number for console output
  !> @param[in]  IOF  File unit number for file output
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph

  subroutine freq_eigenvalues (ndim, IOC, IOF, ierr)

    use SolverModule     , only : getSystemMatrix, resFilFM
    use DenseMatrixModule, only : solveEigenvalues
    use ProgressModule   , only : writeProgress
    use ReportErrorModule, only : allocationError
    use ReportErrorModule, only : reportError, debugFileOnly_p

    integer, intent(in)  :: ndim, IOC, IOF
    integer, intent(out) :: ierr

    !! Local variables
    integer               :: i              !< loop variable
    real(dp), allocatable :: M(:,:)         !< global mass matrix
    real(dp), allocatable :: K(:,:)         !< global stiffness matrix
    real(dp), allocatable :: eigenValues(:) !< eigenvalues of the system
    real(dp), allocatable :: eigenVec(:,:)  !< eigenvectors of the system


    !! --- Logic section ---

    allocate(M(ndim, ndim), K(ndim, ndim), &
         &   eigenValues(nrModes), eigenVec(0,0), STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('freq_eigenvalues')
       return
    end if

    call getSystemMatrix(M, 2, ierr) ! 2 returns the mass matrix
    if (ierr /= 0) goto 900

    call getSystemMatrix(K, 1, ierr) ! 1 returns the stiffness matrix
    if (ierr /= 0) goto 900

    call solveEigenvalues (K, M, eigenValues, eigenVec, ndim, nrModes, 0, ierr)
    if (ierr /= 0) goto 900

    eigenvalues = sqrt(eigenValues)/(2.0_dp*pi_p)
    write(IOC,100) (i,eigenValues(i),i=1,size(eigenValues))
    write(IOF,100) (i,eigenValues(i),i=1,size(eigenValues))
100 format('Natural frequencies (angular, no damping):' / (I8,F33.6))
    if (resFilFM == 1) resFilFM = 2 ! Print new convergence heading for the TDA
    goto 999

900 call reportError (debugFileOnly_p,'freq_eigenvalues')
999 deallocate(M,K,eigenValues,eigenVec)

  end subroutine freq_eigenvalues


  !!============================================================================
  !> @brief Runs a sweep analysis for detecting the eigenfrequencies (optional).
  !> @param[in]  ndim Dimension of the global equation system
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph

  subroutine freq_sweep (ndim, ierr)

    use SolverModule     , only : getSystemMatrix
    use ProgressModule   , only : writeProgress
    use ProfilerModule   , only : startTimer, stopTimer, fra_p
    use ReportErrorModule, only : allocationError
    use ReportErrorModule, only : reportError, debugFileOnly_p

    integer, intent(in)  :: ndim
    integer, intent(out) :: ierr

    !! Local variables
    integer                  :: nc          !< loop variable
    integer                  :: inc         !< total number of increments

    real(dp)   , allocatable :: omega(:)    !< circular frequency

    real(dp)   , allocatable :: M(:,:)      !< global mass matrix
    real(dp)   , allocatable :: K(:,:)      !< global stiffness matrix
    real(dp)   , allocatable :: C(:,:)      !< global damping matrix

    real(dp)                 :: w, dw       !< frequency, frequency increment
    complex(dp), allocatable :: SysMat(:,:) !< global system matrix
    complex(dp), allocatable :: eqVec(:)    !< equation vector
    complex(dp), allocatable :: disp(:,:)   !< result matrix (result vector for all increments)


    !! --- Logic section ---

    allocate(M(ndim, ndim), K(ndim, ndim), C(ndim, ndim), STAT=ierr)
    if (ierr /= 0) goto 999

    call getSystemMatrix(M, 2, ierr) ! 2 returns the global mass matrix
    if (ierr /= 0) goto 998

    call getSystemMatrix(K, 1, ierr) ! 1 returns the global stiffness matrix
    if (ierr /= 0) goto 998

    call getSystemMatrix(C, 3, ierr) ! 3 returns the global damping matrix
    if (ierr /= 0) goto 998

    allocate(SysMat(ndim,ndim), eqVec(ndim), STAT=ierr)
    if (ierr /= 0) goto 999

    w  = 2.0_dp*pi_p*sweep_range(1) ! start frequency
    dw = 0.1_dp ! set frequency increment to 0.1 [Hz]

    inc = int((2.0_dp*pi_p*sweep_range(2)-w)/dw)
    if (associated(plt)) then
       allocate(disp(ndim,inc), omega(inc), STAT=ierr)
       if (ierr /= 0) goto 999
       omega = 0.0_dp
    end if

    call startTimer (fra_p)

    do nc = 1, inc

      eqVec = (0.0_dp, 0.0_dp)
      eqVec(sweep_dof) = cmplx(1,0,dp) ! unit load at sweep_dof
      SysMat = cmplx(-w*w*M + K, w*C, dp)
      call complexEqSolv(ndim, SysMat, eqVec, ierr)
      if (ierr /= 0) exit

      if (associated(plt)) then
         disp(:,nc) = eqVec
         omega(nc) = w/(2.0_dp*pi_p)
      end if
      w = w + dw

      if (.false.) call writeToFile("sweep.txt", omega(nc), eqVec(sweep_dof))

      call writeProgress (nc,inc)
    end do

    if (associated(plt) .and. ierr == 0) then
      call createPyPlotCurves(omega, disp, "sweep", '$\ sweep- $', "sweep", "frequency[w]", type=3, savefigure=.true.)
    end if

    call stopTimer (fra_p)

    if (associated(plt)) then
       deallocate(disp,omega)
    end if
    deallocate(SysMat,eqVec)
998 deallocate(M,K,C)

    if (ierr /= 0) call reportError (debugFileOnly_p,'freq_sweep')
    return

999 ierr = allocationError('freq_sweep')

  contains

    !> @brief Writes data records to the file @a filename.
    subroutine writeToFile (filename,abscissa,ordinate)
      character(len=*), intent(in) :: filename
      real(dp)        , intent(in) :: abscissa
      complex(dp)     , intent(in) :: ordinate
      if (nc == 1) then ! Create new file, replacing existing
         open(31,FILE=filename,FORM='FORMATTED',STATUS='REPLACE')
      else ! Append to existing file
         open(31,FILE=filename,FORM='FORMATTED',STATUS='OLD',ACCESS='APPEND')
      end if
      write(31,600) abscissa,dreal(ordinate),dimag(ordinate)
      close(31)
600   format(F12.6,',',F26.14,',',F26.14)
    end subroutine writeToFile

  end subroutine freq_sweep


  !!============================================================================
  !> @brief Load and motion definition for the frequency domain analysis.
  !> @param[in] forces All force objects in the model
  !> @param[inout] npt Number of sample points
  !> @param[inout] inc Time increment [s]
  !> @param[in]  start Start time for window
  !> @param[in]  IOC   File unit number for console output
  !> @param[out] ierr  Error flag
  !>
  !> @callgraph @callergraph

  subroutine freq_input (forces, npt, inc, start, IOC, ierr)

    use ForceTypeModule    , only : ForceType
    use FunctionTypeModule , only : FunctionValue
    use FileUtilitiesModule, only : findUnitNumber
    use ReportErrorModule  , only : reportError, allocationError
    use ReportErrorModule  , only : note_p, error_p, debugFileOnly_p

    type(ForceType), intent(in)    :: forces(:)
    integer        , intent(inout) :: npt
    real(dp)       , intent(inout) :: inc
    real(dp)       , intent(in)    :: start
    integer        , intent(in)    :: IOC
    integer        , intent(out)   :: ierr

    !! Local variables
    integer             :: i, j, cc, col !< loop variables
    integer             :: e_npt         !< extended data set
    integer             :: io            !< error indicator
    integer             :: iffp          !< file unit number for Fourier file
    character(len=1024) :: line          !< data line
    real(dp)            :: num_array(50) !< help array


    !! --- Logic section ---

    !! Count number of relevant cases
    cc = 0
    e_npt = npt
    do i = 1, size(forces)
       if (forces(i)%loadType > 0) then
          cc = cc + 1
          !! In case of FFT input don't use windowing, reset windowSize
          if (forces(i)%fftFile /= '' .and. windowSize /= 0) then
             windowSize = 0
             call reportError (note_p,'Windowing is deactivated '// &
                  &            'for frequency-dependent load/motion input')
          end if
       end if
    end do

    ! recompute windowing parameters (inc,npt, windowSize)
    if(windowSize > 0) then
      inc = inc*real(npt,dp)
      npt = 2**nint(log(real(npt,dp))/log(2.0_dp))
      windowSize = 2**nint(log(real(windowSize,dp))/log(2.0_dp))
      if(windowSize >= (npt/2)) windowSize = npt/2
      inc = inc/real(npt,dp)
      e_npt = 2*npt
    end if

    !! Initialising data structure
    allocate(plm(cc), STAT=ierr)
    if (ierr /= 0) goto 999

    cc = 0
    num_array = 0.0_dp
    do i = 1, size(forces)
       if (forces(i)%loadType > 0) then

          cc = cc + 1
          allocate(plm(cc)%x(e_npt), plm(cc)%y(e_npt), &
               &   plm(cc)%data(e_npt), STAT=ierr)
          if (ierr /= 0) goto 999

          ! motion: 1 displ, 2 velocity, 3 acceleration
          plm(cc)%motion   = forces(i)%loadType - 1
          plm(cc)%id_fft   = forces(i)%fftFile /= ''
          plm(cc)%triad_id = forces(i)%triad%id%baseId
          plm(cc)%dof      = forces(i)%dof

          if (plm(cc)%id_fft) then

             !! Reading FFT from file
             iffp = findUnitNumber(20)
             open(iffp,FILE=trim(forces(i)%fftFile),STATUS='OLD',IOSTAT=ierr)
             if (ierr /= 0) then
                call reportError (error_p,'Failed to open '//forces(i)%fftFile)
                goto 998
             end if
             write(IOC,*) 'Reading from file '//trim(forces(i)%fftFile)
             ReadLoop: do j = 1, e_npt
                line = '#comment'
                ! read over lines with '#' at the beginning (comment lines)
                do while (index(line,'#') > 0)
                   read(iffp,'(A)',IOSTAT=ierr) line
                   if (ierr /= 0) goto 998
                end do
                ! file contains for each line 3 columns (freq,real,imag)
                ! if not - select the 2 last columns (real,imag)
                ! ccw starts with 0 and is incremented by 1 in the windows loop
                read(line,*,iostat=io) num_array(1:2)
                if (io < 0) cycle ReadLoop
                do col = 3, min(3*(ccw+1),size(num_array))
                   read(line,*,iostat=io) num_array(1:col)
                   if (io == -1) then
                      plm(cc)%x(j) = num_array(col-2)
                      plm(cc)%y(j) = num_array(col-1)
                      cycle ReadLoop
                   else if (io == 0 .and. col == 3*(ccw+1)) then
                      plm(cc)%x(j) = num_array(col-1)
                      plm(cc)%y(j) = num_array(col)
                      cycle ReadLoop
                   end if
                end do
                call reportError (error_p,'Invalid FFT input file '// &
                     &            trim(forces(i)%fftFile), &
                     &            'Last line read: '//trim(line))
                goto 998
             end do ReadLoop
             write(IOC,*) 'Last record -> x: ', plm(cc)%x(e_npt), &
                  &                      'y: ', plm(cc)%y(e_npt)
             close(iffp)
             plm(cc)%data = cmplx(plm(cc)%x,plm(cc)%y,dp)

          else if (associated(forces(i)%engine)) then

             plm(cc)%x(1) = start
             do j = 1, e_npt
                plm(cc)%y(j) = FunctionValue(forces(i)%engine%func,plm(cc)%x(j),ierr)
                if (ierr < 0) then
                   call reportError (debugFileOnly_p,'freq_input')
                   return
                else if (j < e_npt) then
                   plm(cc)%x(j+1) = plm(cc)%x(j) + inc
                end if
             end do

             ! assign values to global data structure data
             plm(cc)%data = cmplx(plm(cc)%y,0,dp)

          end if
       end if
    end do

    return

998 call reportError (debugFileOnly_p,'freq_input')
    return

999 ierr = allocationError('freq_input')

  end subroutine freq_input


  !!============================================================================
  !> @brief Updates mechanism based on displacement, velocity and acceleration.
  !> @param[in]    sam  Data for managing system matrix assembly
  !> @param[inout] sys  System level model data
  !> @param[inout] mech Mechanism components of the model
  !> @param[in]    dis  Displacement increment
  !> @param[in]    vel  Velocity increment
  !> @param[in]    acc  Acceleration increment
  !> @param[out]   ierr Error flag
  !>
  !> @details All position variables are updated with respect to the previous
  !> configuration to not intervene with the time domain response. The current
  !> configuration will then be restored back to the previous one after the
  !> results are saved.
  !>
  !> @callgraph @callergraph

  subroutine freqIncAndUpdate (sam,sys,mech,dis,vel,acc,ierr)

    use SamModule          , only : SamType
    use SystemTypeModule   , only : SystemType
    use MechanismTypeModule, only : MechanismType

    use TriadTypeModule    , only : SetTriadsVelAcc, IncTriadsPos
    use TriadTypeModule    , only : TransSysVecToGlobal, clearTriadForces
    use SupElRoutinesModule, only : SetSupElsVelAcc, IncSupElsGenDofs
    use SupElRoutinesModule, only : updateSupElsStatic
    use MasterSlaveJointTypeModule    , only : SetJointsVelAcc
    use MasterSlaveJointRoutinesModule, only : IncJointsVar
    use SolExtensionModule , only : csExpand

    type(SamType)      , intent(in)    :: sam
    type(SystemType)   , intent(inout) :: sys
    type(MechanismType), intent(inout) :: mech
    real(dp)           , intent(in)    :: dis(:), vel(:), acc(:)
    integer            , intent(out)   :: ierr


    !! --- Logic section ---

    !! Expand the equation-ordered displacements to DOF-order
    call csExpand (sam,dis,sys%rinc)
    call TransSysVecToGlobal (mech%triads,sys%rinc)

    !! Increment all position variables
    call IncTriadsPos     (mech%triads,sys%rinc,.true.)
    call IncSupElsGenDofs (mech%sups,sys%rinc,.true.)
    call IncJointsVar     (mech%joints,sys%rinc,mech%motions,.true.)

    !! Expand the equation-ordered velocities/accelerations to DOF-order
    call csExpand (sam,vel,sys%urd)
    call csExpand (sam,acc,sys%urdd)
    call TransSysVecToGlobal (mech%triads,sys%urd)
    call TransSysVecToGlobal (mech%triads,sys%urdd)

    !! Transfer the updated velocities and accelerations down to object level
    call SetTriadsVelAcc (mech%triads,sys%urd,sys%urdd)
    call SetSupElsVelAcc (mech%sups,sam,sys%urd,sys%urdd)
    call SetJointsVelAcc (mech%joints,sys%urd,sys%urdd)

    !! Update mechanism objects
    call clearTriadForces   (mech%triads)
    call updateSupElsStatic (mech%sups,mech%supLoads,mech%env, &
         &                   sys%time,sys%nIterThisStep,1,ierr)

  end subroutine freqIncAndUpdate


  !!============================================================================
  !> @brief Handles the data segments including modal or direct solution.
  !> @param[inout] data Segment data (on input; loads and constraints)
  !> @param[inout] sUd  Segment displacement data
  !> @param[inout] sUv  Segment velocity data
  !> @param[inout] sUa  Segment acceleration data
  !> @param[in]    mDmp Frequency-dependent modal damping function
  !> @param[in]    del  Step size for time series
  !> @param[in]    npt  Number of sampling points
  !> @param[in]    sdim System dimension
  !> @param[out]   ierr Error indicator
  !>
  !> @callgraph @callergraph

  subroutine freq_segment (data, sUd, sUv, sUa, mDmp, del, npt, sdim, ierr)

    use SolverModule      , only : objectEquations, getSystemMatrix
    use DenseMatrixModule , only : solveEigenvalues
    use ProfilerModule    , only : startTimer, stopTimer, fra_p
    use FunctionTypeModule, only : FunctionType, FunctionValue
    use ReportErrorModule , only : allocationError
    use ReportErrorModule , only : reportError, debugFileOnly_p

    complex(dp), intent(inout) :: data(:,:), sUd(:,:), sUv(:,:), sUa(:,:)

    type(FunctionType), pointer, intent(in) :: mDmp

    real(dp), intent(in)  :: del
    integer , intent(in)  :: npt, sdim
    integer , intent(out) :: ierr

    !! Local variables
    integer                             :: cc                     !< counter variable
    integer                             :: i, j, k                !< loop variables
    integer                             :: ndim                   !< dimension of equation system
    integer                             :: meqn(6)                !< array contains triads global dof numbers
    integer                             :: ndof                   !< number of dgree of freedoms
    integer                             :: pos                    !< position within the array
    integer                             :: nrMotion               !< nr of motion cases
    integer             , allocatable   :: dof_pattern(:)         !< dof scheme for equation system
    integer             , allocatable   :: pattern(:)             !< dof scheme for RHS

    real(dp)            , allocatable   :: Mn(:,:)                !< mass matrix
    real(dp)            , allocatable   :: Cn(:,:)                !< damping matrix
    real(dp)            , allocatable   :: Kn(:,:)                !< stiffness matrix

    real(dp)            , allocatable   :: M1(:,:)                !< storing mass matrix columns (constraint handling)
    real(dp)            , allocatable   :: C1(:,:)                !< storing damping matrix columns (constraint handling)
    real(dp)            , allocatable   :: K1(:,:)                !< storing stiffness matrix columns (constraint handling)

    real(dp)            , allocatable   :: tA(:,:)                !< temporary matrix
    real(dp)            , allocatable   :: mM(:,:)                !< modal mass matrix (eVec^T*M*eVec)
    real(dp)            , allocatable   :: mK(:,:)                !< modal stiffness matrix (eVec^T*K*eVec)
    real(dp)            , allocatable   :: mC(:,:)                !< modal damping matrix (eVec^T*C*eVec)

    complex(dp)         , allocatable   :: mF(:)                  !< modal force vector
    complex(dp)         , allocatable   :: mU(:)                  !< modal displacement vector
    complex(dp)         , allocatable   :: mSys(:)                !< modal equation system

    complex(dp)         , allocatable   :: Sys(:,:)               !< system matrix
    complex(dp)         , allocatable   :: Rhs(:)                 !< right hand side vector

    real(dp)            , allocatable   :: eVal(:)                !< eigenvalue
    real(dp)            , allocatable   :: eVec(:,:)              !< mode shapes

    real(dp)                            :: w                      !< circular frequency
    real(dp)                            :: zeta                   !< modal damping ratio


    !! --- Logic section ---

    ! fft for all FRA input loads
    do i = 1, size(data,1)
       if (plm(i)%id_fft) cycle
       call fft_complex_1D(npt, data(i,:), del, "forward", .false., ierr)
       if (ierr /= 0) goto 998
    end do

    nrMotion = 0
    do i = 1, size(plm)
       if (plm(i)%motion > 0) nrMotion = nrMotion + 1
    end do

    ndim = sdim

    allocate(Mn(ndim,ndim), Cn(ndim,ndim), Kn(ndim,ndim), Rhs(sdim), STAT=ierr)
    if (ierr /= 0) goto 999

    call stopTimer (fra_p) ! because measuring assembly time below

    call getSystemMatrix(Kn, 1, ierr) ! 1 returns the stiffness matrix
    if (ierr < 0) goto 998

    call getSystemMatrix(Mn, 2, ierr) ! 2 returns the mass matrix
    if (ierr < 0) goto 998

    call getSystemMatrix(Cn, 3, ierr) ! 3 returns the damping matrix
    if (ierr < 0) goto 998

    call startTimer (fra_p) ! resume timing frequency response analysis

    if(nrMotion > 0) then
      allocate(K1(ndim,nrMotion), STAT=ierr)
      if (ierr /= 0) goto 999
      allocate(M1(ndim,nrMotion), STAT=ierr)
      if (ierr /= 0) goto 999
      allocate(C1(ndim,nrMotion), STAT=ierr)
      if (ierr /= 0) goto 999

      allocate(dof_pattern(ndim), pattern(ndim), STAT=ierr)
      if (ierr /= 0) goto 999

      pattern     = (/ (i, i = 1, ndim) /)
      dof_pattern = (/ (i, i = 1, ndim) /)

      ! decrease size of the global equation system
      ! set 0 for exchanged column and decrease numbers behind by 1
      k = 0
      do i = 1, size(plm)
        if(plm(i)%motion > 0) then
          k = k+1
          call objectEquations(plm(i)%triad_id, meqn, ndof)
          pos = meqn(plm(i)%dof)

          ndim = ndim - 1
          call remove_RowCol_Matrix(Kn, dof_pattern(pos), K1(:,k))
          if (ierr /= 0) goto 999
          call remove_RowCol_Matrix(Mn, dof_pattern(pos), M1(:,k))
          if (ierr /= 0) goto 999
          call remove_RowCol_Matrix(Cn, dof_pattern(pos), C1(:,k))
          if (ierr /= 0) goto 999

          dof_pattern(pos) = 0
          do j = pos+1, size(dof_pattern)
            if(dof_pattern(j) > 0) dof_pattern(j) = dof_pattern(j) - 1
          end do
        end if
      end do
    else
       ! Allocate zero-size arrays.
       ! Just to avoid compiler warning on maybe-uninitialized below.
       allocate(K1(0,0),M1(0,0),C1(0,0),pattern(0),dof_pattern(0))
    end if

    ! in case of model decomposition compute eigenvalues and
    ! eigenvectors for the current segment
    if(nrModes > 0) then
      allocate(mK(ndim,ndim), mM(ndim,ndim), &
           &   eVal(ndim), eVec(ndim, nrModes), STAT=ierr)
      if (ierr /= 0) goto 999

      mK = Kn
      mM = Mn
      call solveEigenvalues (mK, mM, eVal, eVec, ndim, nrModes, nrModes, ierr)
      deallocate(mK,mM)
      if (ierr /= 0) goto 998

      allocate(tA(nrModes, ndim), mK(nrModes, nrModes), mM(nrModes, nrModes), &
           &   mC(nrModes, nrModes), STAT=ierr)
      if (ierr /= 0) goto 999

      !! create generalized mass and stiffness matrix
      tA = matmul(transpose(eVec), Mn)
      mM = matmul(tA, eVec)
      tA = matmul(transpose(eVec), Kn)
      mK = matmul(tA, eVec)

      if(associated(mDmp)) then
        ! calculate damping matrix based on damping for each mode
        ! zeta_k percentage of critical damping for each mode (user input)
        ! zeta_k = f(omega) values for eigenValues are extra/interpolated - straight line interpolation
        ! mC(k,k) = 2.0_dp*zeta(k)*mM(k,k)*sqrt(eVal(k))  (sqrt because eigenvalue solver delivers w^2)
        do k = 1, nrModes
          zeta = FunctionValue(mDmp, dsqrt(eVal(k)), ierr)
          if(ierr < 0) goto 998
          mC(k,k) = 2.0_dp*zeta*mM(k,k)*dsqrt(eVal(k))
        end do
      else
        ! use rayleigh damping
        tA = matmul(transpose(eVec), Cn)
        mC = matmul(tA, eVec)
      end if

      allocate(mSys(nrModes), mF(nrModes), mU(nrModes), STAT=ierr)
      if (ierr /= 0) goto 999
    else
      allocate(Sys(ndim, ndim), STAT=ierr)
      if (ierr /= 0) goto 999
    endif

    ! loop over samples
    do cc = 1, npt/2+1

      ndim = sdim  ! reset dimension

      ! circular frequency
      w = 2.0_dp*pi_p*real(cc-1,dp)/real(npt,dp)/del

      ! building rhs (loads, constraints)
      Rhs = (0.0_dp, 0.0_dp)
      j = 0
      do i = 1, size(plm)
         if (plm(i)%motion == 0) then
            call objectEquations(plm(i)%triad_id, meqn, ndof)
            pos = meqn(plm(i)%dof)
            Rhs(pos) = Rhs(pos) + data(i,cc)
         else if (plm(i)%motion > 0) then
            j = j + 1
            ! velocity
            if (w > 0.0_dp .and. plm(i)%motion == 2) data(i,cc) = data(i,cc)/cmplx(   0,w,dp)
            ! acceleration
            if (w > 0.0_dp .and. plm(i)%motion == 3) data(i,cc) = data(i,cc)/cmplx(-w*w,0,dp)
            do k = 1, ndim
               Rhs(k) = Rhs(k) - cmplx(-w*w*M1(k,j) + K1(k,j), w*C1(k,j), dp)*data(i,cc)
            end do
         end if
      end do

      ! decrease size of the rhs system (prescribed)
      if(nrMotion > 0) then
        pattern = (/ (i, i = 1, ndim) /)
        do i = 1, size(dof_pattern)
          if(dof_pattern(i) == 0) then
            ndim = ndim - 1
            call remove_Row_Vec(Rhs,pattern(i))
            if (ierr /= 0) goto 999
            pattern(i) = 0
            do j = i+1, size(pattern)
              if(pattern(j) > 0) pattern(j)=pattern(j)-1
            enddo
          endif
        enddo
      end if

      ! solution vector for modal and direct solution
      if (nrModes > 0) then
        mF  = matmul(transpose(eVec), Rhs)
        do k = 1, nrModes
          mSys(k) = cmplx(-w*w + eVal(k), w*mC(k,k), dp)
        end do
        ! modal displacement vector
        mU = mF/mSys
        ! transform modal displacement vector to displacement vector
        Rhs = matmul(eVec, mU)
      else
        ! direct solution
        Sys = cmplx(-w*w*Mn + Kn, w*Cn, dp)
        call complexEqSolv(ndim, Sys, Rhs, ierr)
        if (ierr /= 0) goto 998
      end if

      ! rebuilding solution vector (in case of prescribed motions)
      if(nrMotion > 0) then
        do i = 1, size(dof_pattern)
          if (dof_pattern(i) == 0) then
            do k = 1, size(plm)
              if(plm(k)%motion > 0) then
                call objectEquations(plm(k)%triad_id, meqn, ndof)
                if (i == meqn(plm(k)%dof)) then
                   ndim = ndim + 1
                   call add_Row_Vec(data(k,cc),Rhs,i)
                   if (ierr /= 0) goto 999
                end if
              end if
            end do
          end if
        end do
      end if

      ! solution matrices
      sUd(:,cc) = Rhs                    ! solution matrix for displacements
      sUv(:,cc) = Rhs * cmplx(   0,w,dp) ! solution matrix for velocities
      sUa(:,cc) = Rhs * cmplx(-w*w,0,dp) ! solution matrix for accelerations

    end do

    if (nrModes > 0) then
       deallocate(mSys,mF,mU)
       deallocate(tA,mK,mM,mC)
       deallocate(eVal,eVec)
    else
       deallocate(Sys)
    end if

    deallocate(dof_pattern,pattern)
    deallocate(K1,M1,C1)
    deallocate(Kn,Mn,Cn,Rhs)

    ! iFFT transform results to the time domain
    do i = 1, ndim
      do j = npt/2+2, npt
        pos = npt+2-j
        sUd(i,j) = conjg(sUd(i, pos))
        sUv(i,j) = conjg(sUv(i, pos))
        sUa(i,j) = conjg(sUa(i, pos))
      end do

      call fft_complex_1D(npt, sUd(i,:), del, "backward", ierr=ierr)
      if (ierr /= 0) goto 998
      call fft_complex_1D(npt, sUv(i,:), del, "backward", ierr=ierr)
      if (ierr /= 0) goto 998
      call fft_complex_1D(npt, sUa(i,:), del, "backward", ierr=ierr)
      if (ierr /= 0) goto 998
    end do

    return

998 call reportError (debugFileOnly_p,'freq_segment')
    return

999 ierr = allocationError('freq_segment')

  contains

    !> @brief Removes a column and corresponding row in a quadratic matrix.
    subroutine remove_RowCol_Matrix (A,i,sCol)
      integer              , intent(in)    :: i
      real(dp), allocatable, intent(inout) :: A(:,:)
      real(dp)             , intent(inout) :: sCol(:)
      real(dp), allocatable                :: B(:,:)

      allocate(B(ndim,ndim),STAT=ierr)
      if (ierr /= 0) return

      B(1:i-1,1:i-1) = A(1:i-1,1:i-1)
      B(1:i-1,i:)    = A(1:i-1,i+1:)
      B(i:,1:i-1)    = A(i+1:,1:i-1)
      B(i:,i:)       = A(i+1:,i+1:)

      sCol = A(:,i)

      deallocate(A)
      allocate(A(ndim,ndim),STAT=ierr)
      if (ierr /= 0) return

      call move_alloc(B,A) ! move B to A, B is automatically deallocated

    end subroutine remove_RowCol_Matrix

    !> @brief Extends the 1D array @a A by inserting a value at position @a i.
    subroutine add_Row_Vec (value,A,i)
      integer                 , intent(in)    :: i
      complex(dp)             , intent(in)    :: value
      complex(dp), allocatable, intent(inout) :: A(:)
      complex(dp), allocatable                :: B(:)

      allocate(B(ndim),STAT=ierr)
      if (ierr /= 0) return

      B(1:i-1) = A(1:i-1)
      B(i)     = value
      B(i+1:)  = A(i:)

      deallocate(A)
      allocate(A(ndim),stat=ierr)
      if (ierr /= 0) return

      call move_alloc(B,A) ! move B to A, B is automatically deallocated

    end subroutine add_Row_Vec

    !> @brief Reduces the 1D array @a A by deleting the @a i'th element.
    subroutine remove_Row_Vec (A,i)
      integer                 , intent(in)    :: i
      complex(dp), allocatable, intent(inout) :: A(:)
      complex(dp), allocatable                :: B(:)

      allocate(B(ndim),stat=ierr)
      if (ierr /= 0) return

      B(1:i-1) = A(1:i-1)
      B(i:)    = A(i+1:)

      deallocate(A)
      allocate(A(ndim),stat=ierr)
      if (ierr /= 0) return

      call move_alloc(B,A) ! move B to A, B is automatically deallocated

    end subroutine remove_Row_Vec

  end subroutine freq_segment


  !!============================================================================
  !> @brief Main driver for the segmented frequency response analysis.
  !> @param[in]    sam  Data for managing system matrix assembly
  !> @param[inout] sys  System level model data
  !> @param[inout] mech Mechanism components of the model
  !> @param[in]    t    time interval
  !> @param[in]    IOC  File unit number for console output
  !> @param[in]    IOF  File unit number for file output
  !> @param[out]   ierr Error flag
  !>
  !> @details Data segmenting and windowing.
  !>
  !> @callgraph @callergraph

  subroutine freq_seg_analysis (sam, sys, mech, t, IOC, IOF, ierr)

    use SamModule          , only : SamType
    use SystemTypeModule   , only : SystemType
    use MechanismTypeModule, only : MechanismType, mDmpFuncId
    use FunctionTypeModule , only : FunctionType, GetPtrToId
    use SaveModule         , only : writeSolverDB5
    use SolverModule       , only : resFilFM
    use ProfilerModule     , only : startTimer, stopTimer, fra_p, upd_p, sav_p
    use ReportErrorModule  , only : allocationError
    use ReportErrorModule  , only : reportError, debugFileOnly_p

    type(SamType)      , intent(in)    :: sam
    type(SystemType)   , intent(inout) :: sys
    type(MechanismType), intent(inout) :: mech
    real(dp)           , intent(in)    :: t
    integer            , intent(in)    :: IOC, IOF
    integer            , intent(out)   :: ierr

    !! Local variables
    integer  :: n, m                 !< number of sample points (total and per window)
    integer  :: dseg_start, dseg_end !< segment start and end
    integer  :: i, k, iw, jw         !< loop variables
    integer  :: cw                   !< window/segment counter (local)
    integer  :: sdim                 !< size of the system
    integer  :: skip, nw, pos, w_end !< windowing variables
    real(dp) :: start, inc           !< start time and sample period

    complex(dp), allocatable       :: Ud(:,:)    !< displacement matrix
    complex(dp), allocatable       :: Uv(:,:)    !< velocity matrix
    complex(dp), allocatable       :: Ua(:,:)    !< acceleration matrix
    complex(dp), allocatable       :: data(:,:)  !< data used in segmenting
    complex(dp), allocatable, save :: sUd(:,:)   !< displacements (segment/over segment border)
    complex(dp), allocatable, save :: sUv(:,:)   !< velocities (segment/over segment border)
    complex(dp), allocatable, save :: sUa(:,:)   !< accelerations (segment/over segment border)
    real(dp)   , allocatable       :: hwin(:)    !< window function (hamming or rectangular)
    real(dp)   , allocatable       :: x(:)       !< abscissa values
    real(dp)   , allocatable       :: ax(:)      !< abscissa values
    real(dp)   , allocatable       :: time(:)    !< time array

    type(FunctionType), pointer :: modalDmp !< frequency-dependent modal damping function


    !! --- Logic section ---

    call startTimer (fra_p)

    write(IOF,600) sys%time
600 format(/' Frequency response analysis at time [s]:', 1PE12.5)
    if (resFilFM == 1) resFilFM = 2 ! Print new convergence heading for the TDA

    inc     = 1.0_dp/fs
    start   = sys%time
    n       = 2*int(0.5_dp*t*fs) ! t*fs should be an even number

    ! initialisation of load vectors/prescribed motions
    call freq_input(mech%forces, n, inc, start, IOC, ierr)
    if (ierr < 0) goto 998

    write(IOC,*) "starttime: ", start, "endtime: ", start+t, "inc: ", inc
    write(IOC,*) "time_window: ", t
    write(IOC,*) "number of increments: ", n

    if (nrModes > 0 .and. mDmpFuncId > 0) then
       ! find the global modal damping ratio function
       modalDmp => GetPtrToId(mech%functions,mDmpFuncId)
       if (associated(modalDmp)) then
          write(IOC,*) "Modal Damping activated"
       else
          mDmpFuncId = 0
       end if
    else
       nullify(modalDmp)
    end if

    sdim = sam%neq
    allocate(Ud(sdim, n), Uv(sdim, n), Ua(sdim, n), STAT=ierr)
    if (ierr /= 0) goto 999
    Ud = (0.0_dp, 0.0_dp)
    Uv = (0.0_dp, 0.0_dp)
    Ua = (0.0_dp, 0.0_dp)

    ! windowSize vs. samples
    if (windowSize < 2) then
       m = n
       skip = m
    else
       m = windowSize
       skip = m/2
    end if
    w_end = n-skip

    ! calculate the number of windows/slices
    nw = 1 + (n-m)/skip

    allocate(hwin(m+1), data(size(plm),m), x(n+1), STAT=ierr)
    if (ierr /= 0) goto 999

    if(ccw == 0) then
      write(IOC,*) "ALLOCATE memory for data in intersection zone"
      allocate(sUd(sdim, m), sUv(sdim, m), sUa(sdim, m), STAT=ierr)
      if (ierr /= 0) goto 999
      sUd = (0.0_dp, 0.0_dp)
      sUv = (0.0_dp, 0.0_dp)
      sUa = (0.0_dp, 0.0_dp)
    end if

    x = (/ (real(i,dp)*inc,i=0,n) /)

    if(windowSize /= 0) then
      hwin = 0.5_dp*(1.0_dp-cos((2.0_dp*pi_p/(real(m,dp)*inc))*x(1:m+1))) ! Hamming window function
      if (.false.) then
        allocate(ax(m+1),STAT=k)
      else
        k = 1
      end if
      if (k == 0) then
        ax = (/ (i,i=0,m) /)
        call createPyPlot(ax, hwin, 'hann','$\ signal (x)$','n_fft', 'points', istat=k)
        deallocate(ax)
      end if
      write(IOF,*) "Hann window used"
    else
      hwin = 1.0_dp
      write(IOF,*) "Rectangular window used"
    end if

    ! loop over all segments/windows
    cw = 0
    dseg_end = m
    do i = 0, w_end, skip
      dseg_start = i+1
      dseg_end   = i+m

      cw  = cw+1
      ccw = ccw+1

      ! add values from the border segment into the first half of the new segment
      if (cw == 1 .and. nw > 1) then
         iw = dseg_start
         jw = m/2+1
         do k = 1, m/2
            call ZAXPY (sdim,cmplx(1,0,dp),sUd(1,jw),1,Ud(1,iw),1)
            call ZAXPY (sdim,cmplx(1,0,dp),sUv(1,jw),1,Uv(1,iw),1)
            call ZAXPY (sdim,cmplx(1,0,dp),sUa(1,jw),1,Ua(1,iw),1)
            iw = iw + 1
            jw = jw + 1
         end do
      end if

      write(IOC,610) ccw, seg_start+dseg_start, seg_start+dseg_end
      write(IOF,610) ccw, seg_start+dseg_start, seg_start+dseg_end
610   format(5X,'seg[',I4,' ]  start - end:', I10, ' --', I10)

      ! extract segment data
      do k = 1, size(plm)
         data(k,:) = plm(k)%data(dseg_start:dseg_end)*hwin(1:m)
      end do

      call freq_segment (data, sUd, sUv, sUa, modalDmp, inc, m, sdim, ierr)
      if (ierr /= 0) goto 998

      ! manage to add values from the first half of the border segment into the
      ! second half of the last segment
      if (cw > nw .and. nw > 1) then
        pos = m/2+1
        dseg_start = dseg_start-skip
        dseg_end   = dseg_end-skip
      else
        pos = 1
      end if

      iw = dseg_start
      jw = 1
      do k = pos, m
         call ZAXPY (sdim,cmplx(1,0,dp),sUd(1,jw),1,Ud(1,iw),1)
         call ZAXPY (sdim,cmplx(1,0,dp),sUv(1,jw),1,Uv(1,iw),1)
         call ZAXPY (sdim,cmplx(1,0,dp),sUa(1,jw),1,Ua(1,iw),1)
         iw = iw + 1
         jw = jw + 1
      end do

    end do

    seg_start = seg_start + dseg_end
    call stopTimer (fra_p)

    ! output plots (using python facility)
    if (associated(plt)) then
       allocate(time(n),STAT=k)
    else
       k = 1
    end if
    if (k == 0) then
       do i = 1, n
          time(i) = start + real(i-1,dp)*inc
       end do
       call createPyPlotCurves(time, Ud, "displacement",'$\ dis- $', "displacement_seg", 'time', &
            &                  type=1, savefigure=.true.)
       if (.false.) then
          call createPyPlotCurves(time, Uv, "velocity",'$\ vel- $', "velocity_seg", 'time', &
               &                  type=1, savefigure=.true.)
          call createPyPlotCurves(time, Ua, "acceleration",'$\ acc- $', "acceleration_seg", 'time', &
               &                  type=1, savefigure=.true.)
       end if
       deallocate(time)
    end if

    ! Update and save the equivalent time-domain response to results database
    do i = 1, n

       call startTimer (upd_p)
       sys%time = start + real(i-1,dp)*inc
       call freqIncAndUpdate (sam,sys,mech, &
            &                 dreal(Ud(:,i)), &
            &                 dreal(Uv(:,i)), &
            &                 dreal(Ua(:,i)), ierr)
       call stopTimer (upd_p)
       if (ierr /= 0) goto 998

       call startTimer (sav_p)
       call writeSolverDB5 (sys,mech%triads,mech%sups,mech%masses,ierr)
       call stopTimer (sav_p)
       if (ierr < 0) goto 998

    end do
    sys%time = start
    ierr = 0

    if (start+t+0.1_dp*sys%minInc > sys%tEnd .and. allocated(sUd)) then
       deallocate(sUd,sUv,sUa)
    end if
    deallocate(hwin,data,x)
    deallocate(Ud,Uv,Ua)

    return

998 call reportError (debugFileOnly_p,'freq_seg_analysis')
    return

999 ierr = allocationError('freq_seg_analysis')

  end subroutine freq_seg_analysis


  !!============================================================================
  !> @brief Solving a complex system of linear equations.
  !> @param[in]    n    matrix dimension
  !> @param[inout] A    contains complex matrix A
  !> @param[inout] B    contains rhs and the solution vector
  !> @param[out]   ierr error flag
  !>
  !> @details Number of right hand sides is set to 1.
  !>
  !> @callgraph @callergraph

  subroutine complexEqSolv(n, A, B, ierr)

    use ScratchArrayModule, only : getIntegerScratchArray
    use ReportErrorModule , only : reportError, debugFileOnly_p, error_p

    integer    , intent(in)    :: n
    complex(dp), intent(inout) :: A(:,:), B(:)
    integer    , intent(out)   :: ierr

    !! Local variables
    integer, parameter :: nrhs = 1
    integer, pointer   :: ipiv(:)
    character(len=32)  :: cdiag


    !! --- Logic section ---

    ipiv => getIntegerScratchArray(n,ierr)
    if (ierr < 0) then
       call reportError (debugFileOnly_p,'complexEqSolv')
       return
    end if

    !! solve the equations A*X = B.
    call zgesv(n, nrhs, A(1,1), n, ipiv(1), B(1), n, ierr)

    !! check for the exact singularity.
    if (ierr > 0) then
       write(cdiag,"('U(',I6,' ,',I6,' ) is zero,')") ierr,ierr
       call reportError (error_p,'The diagonal element of the triangular'// &
            &            ' factor of A, '//cdiag,'so that A is singular;'// &
            &            ' the solution could not be computed.', &
            &            addString='complexEqSolv')
       ierr = -ierr
    else if (ierr < 0) then
       call reportError (error_p,'LAPACK::ZGESV failed',IERR=ierr, &
            &            addString='complexEqSolv')
    end if

  end subroutine complexEqSolv


  !!============================================================================
  !> @brief Fourier transformation for a single periodic sequence within a complex array.
  !> @param[in] n         number of points
  !> @param[inout] data   input data (complex)
  !> @param[in] del       data spacing
  !> @param[in] transform forward or backward fourier transformation
  !> @param[in] viewing   optional parameter
  !> @param[out] ierr     error flag
  !>
  !> @callgraph @callergraph

  subroutine fft_complex_1D(n, data, del, transform, viewing, ierr)

    use ReportErrorModule, only : allocationError, internalError
    use ReportErrorModule, only : reportError, error_p

    integer         , intent(in)    :: n
    real(dp)        , intent(in)    :: del
    complex(dp)     , intent(inout) :: data(:)
    character(len=*), intent(in)    :: transform
    logical,optional, intent(in)    :: viewing
    integer         , intent(out)   :: ierr

    !! Local variables
    logical  :: view_image
    integer  :: inc, lenc, lensav, lenwrk, k
    real(dp) :: t

    real(dp)   , allocatable, dimension(:) :: wsave, work, freq
    complex(dp), allocatable, dimension(:) :: c


    !! --- Logic section ---

    view_image = .true.
    if(present(viewing)) view_image = viewing

    ! allocate and initialize the wsave array
    lenwrk = 2*n
    lensav = 2*n + int(log(real(N,dp))/log(2.0_dp)) + 4
    allocate(wsave(lensav), work(lenwrk), c(n), freq(n), STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('fft_complex_1D')
       return
    end if

    c = data  ! copy data into c

    ! initialisation
#ifdef FT_HAS_FFTPACK
    call cfft1i(n, wsave(1), lensav, ierr)
    if (ierr /= 0) goto 999
#endif

    ! compute the (forward) FFT coefficients.
    inc  = 1
    lenc = n

    if(transform == "forward") then
#ifdef FT_HAS_FFTPACK
       call cfft1f(n, inc, c, lenc, wsave(1), lensav, work(1), lenwrk, ierr)
#endif

       do k = 1, n/2+1
          freq(k) = real(k-1,dp)/(n*del)
       end do

      ! write coefficients to file
      if(.false.) then
      open(1, file = "fftpack_fourier.txt")
      do k = 1, n/2+1
        write(1, '(f12.6,A, 4X, f22.11, A, 4X, f22.11)') freq(k), ',', dreal(c(k)), ',', dimag(c(k))
      end do
      close(1)
      end if

      ! output pyplot
      if (view_image .and. associated(plt)) then
        call createPyPlot(freq(1:n/2+1), abs(c(1:n/2+1)), 'fft load','$\ fft (x)$', 'fftpack', 'frequency', istat=ierr)
      end if

    else if(transform == "backward") then

      ! compute backward(inverse) FFT
#ifdef FT_HAS_FFTPACK
      call cfft1b(n, inc, c, lenc, wsave(1), lensav, work(1), lenwrk, ierr)
#endif

      if(.false.) then
      open(20, file = "fftpack_invfourier.txt")
      t = 0.0_dp
      do k = 1, n
        write(20, '(f12.6,A, 4X, f16.11, A, 4X, f16.11)') t, ',', dreal(c(k)), ',', dimag(c(k))
        t = t + del
      end do
      close(20)
      end if

    else
      ierr = internalError('fft_complex_1D: Invalid FFT-type '//transform)
    end if

    data = c

#ifdef FT_HAS_FFTPACK
999 continue
#endif
    if (ierr > 0) then
       call reportError(error_p,'FFT failure',addString='fft_complex_1D')
    end if
    deallocate(wsave,work,c,freq)

  end subroutine fft_complex_1D

  !!============================================================================
  !> @brief Creates output file for curve plotting in python.
  !> @param[in] abscissa   abscissa values
  !> @param[in] ordinate   ordinate values
  !> @param[in] atitle     graph plot title
  !> @param[in] alabel     describes the plot
  !> @param[in] afilename  filename for writing the plot values
  !> @param[in] axlabel    describes the x axis
  !> @param[in] savefigure flag for saving the plot as *.png
  !> @param[out] istat     status flag
  !>
  !> @details Optionally also creates *.png graphics.
  !>
  !> @callgraph @callergraph

  subroutine createPyPlot(abscissa, ordinate, atitle, alabel, afilename, axlabel, savefigure, istat)

    real(dp)        , intent(in)  :: abscissa(:), ordinate(:)
    character(len=*), intent(in)  :: atitle, alabel, afilename, axlabel
    logical,optional, intent(in)  :: savefigure
    integer         , intent(out) :: istat


    !! --- Logic section ---

#ifdef FT_HAS_PYPLOT
    call plt%initialize (grid=.true., xlabel=axlabel, figsize = [14,10], &
         &               title=atitle, legend=.true., axis_equal=.false., &
         &               tight_layout=.true., istat=istat)
    if (istat /= 0) return

    call plt%add_plot (abscissa, ordinate, label=alabel, linestyle='b-o', &
         &             markersize=5, linewidth=2, istat=istat)

    if (present(savefigure)) then
       if (savefigure) call plt%savefig(afilename//'.png', pyfile=afilename//'.py', istat=istat)
    end if

    call plt%showfig(pyfile=afilename//'.py', istat=istat)
#else
    print *,' ** createPyPlot: Not included', size(abscissa), size(ordinate), &
         trim(atitle), trim(alabel), trim(afilename), trim(axlabel), &
         present(savefigure)
    istat = 0
#endif

  end subroutine createPyPlot


  !!============================================================================
  !> @brief Creates output file for curve plotting in python.
  !> @param[in] abscissa   abscissa values
  !> @param[in] ordinate   ordinate values
  !> @param[in] atitle     graph plot title
  !> @param[in] alabel     describes the plot
  !> @param[in] afilename  filename for writing the plot values
  !> @param[in] axlabel    describes the x axis
  !> @param[in] type       real(1), imaginary(2) and absolute plot(3)
  !> @param[in] savefigure flag for saving the plot as *.png
  !>
  !> @details Optionally also creates *.png graphics.
  !>
  !> @callgraph @callergraph

  subroutine createPyPlotCurves(abscissa, ordinate, atitle, alabel, afilename, axlabel, type, savefigure)

    real(dp)        , intent(in) :: abscissa(:)
    complex(dp)     , intent(in) :: ordinate(:,:)
    character(len=*), intent(in) :: atitle, alabel, afilename, axlabel
    integer,optional, intent(in) :: type
    logical,optional, intent(in) :: savefigure

#ifdef FT_HAS_PYPLOT
    !! Local variables
    integer                     :: i, k, z, id_dof, istat, plot_type
    character(len=5)            :: s1
    character(len=4), parameter :: mdof(6) = (/ "(Dx)","(Dy)","(Dz)","(Rx)","(Ry)","(Rz)" /)
    character(len=1), parameter :: mcol(7) = (/ "b","g","r","c","m","y","b" /)
    character(len=5), parameter :: c1 = 'triad'


    !! --- Logic section ---

    if (present(type)) then
       plot_type = type
    else
       plot_type = 0
    end if

    call plt%initialize (grid=.true., xlabel=axlabel, figsize = [14,10], &
         &               title=atitle, legend=.true., axis_equal=.false., &
         &               tight_layout=.true., istat=istat)
    if (istat /= 0) return

    z = 1
    do i = 1, size(out_dof)
      k = out_dof(i)
      if (k > 0) then
        write(s1,"(I5)") output_triad_dof(2*i-1)
        id_dof = output_triad_dof(2*i)

        if (plot_type == 1) then
          call plt%add_plot(abscissa, dreal(ordinate(k,:)), &
               &            label=alabel//c1//s1//mdof(id_dof), linestyle=mcol(z), &
               &            markersize=5, linewidth=2, istat=istat)
        else if (plot_type == 2) then
          call plt%add_plot(abscissa, dimag(ordinate(k,:)), &
               &            label=alabel//c1//s1//mdof(id_dof), linestyle=mcol(z), &
               &            markersize=5, linewidth=2, istat=istat)
        else
          call plt%add_plot(abscissa, abs(ordinate(k,:)), &
               &            label=alabel//c1//s1//mdof(id_dof), linestyle=mcol(z), &
               &            markersize=5, linewidth=2, istat=istat)
        end if
        z = z+1
        if(z == size(mcol)) z = 1
      end if
    end do

    if (present(savefigure)) then
       if (savefigure) call plt%savefig(afilename//'.png', pyfile=afilename//'.py', istat=istat)
    end if

    call plt%showfig(pyfile=afilename//'.py', istat=istat)
#else
    print *,' ** createPyPlot: Not included', size(abscissa), size(ordinate), &
         trim(atitle), trim(alabel), trim(afilename), trim(axlabel), &
         present(type), present(savefigure)
#endif

  end subroutine createPyPlotCurves

end module FreqResponseModule
