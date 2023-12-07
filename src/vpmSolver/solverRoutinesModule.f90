!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file solverRoutinesModule.f90
!>
!> @brief Generic solver utilities.

!!==============================================================================
!> @brief Module with some generic solver utilities.
!>
!> @details This module contains utility subroutines for the Dynamics Solver
!> that are not related to any specific solution algorithm or mechanism part.
!> Many of the subroutines and functions of this module are not documented.
!> You have to configure doxygen with the option ENABLED_SECTIONS = FULL_DOC
!> to extract detailed documentation of those subroutines and functions.

module SolverRoutinesModule

  use KindModule   , only : dp
  use sprKindModule, only : ik
  use SaveVTFModule, only : openVTF

  implicit none

  integer(ik), allocatable, save :: meqErr(:) !< All singular equations found

  real(dp), private, save :: lasTime(4) = -1.0e99_dp !< Previous save times

  private :: nodeEntity


contains

  !> @cond FULL_DOC
  !!============================================================================
  !> @brief Returns a description of the entity a given node is associated with.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 18 Apr 2002

  function nodeEntity (node,ldof,sups,triads,joints,writePos) result(entity)

    use KindModule                , only : pi_p
    use SupElTypeModule           , only : SupElType
    use TriadTypeModule           , only : TriadType
    use MasterSlaveJointTypeModule, only : MasterSlaveJointType, jointTypeName_p
    use IdTypeModule              , only : getId
    use rotationModule            , only : FFa_glbEulerZYX

    integer                   , intent(in) :: node, ldof
    type(SupElType)           , intent(in) :: sups(:)
    type(TriadType)           , intent(in) :: triads(:)
    type(MasterSlaveJointType), intent(in) :: joints(:)
    logical, optional         , intent(in) :: writePos
    character(len=128)                     :: entity

    !! Local variables
    integer           :: i
    real(dp)          :: vec(3)
    character(len=72) :: cVal

    character(len=6), parameter :: dofSymbol_p(7) = (/ 'Tx    ', &
         &                                             'Ty    ', &
         &                                             'Tz    ', &
         &                                             'Rx    ', &
         &                                             'Ry    ', &
         &                                             'Rz    ', &
         &                                             'Slider' /)

    !! --- Logic section ---

    do i = 1, size(triads)
       if (triads(i)%samNodNum == node) then
          if (lDof < 1 .or. lDof > size(dofSymbol_p)) exit
          entity = trim(dofSymbol_p(lDof))//' in Triad'//getId(triads(i)%id)
          if (present(writePos)) then
             if (writePos) then
                call FFa_glbEulerZYX (triads(i)%ur(:,1:3),vec)
                write(cVal,600) triads(i)%ur(:,4),vec*180.0_dp/pi_p
                entity(44:) = cVal
             end if
          end if
          return
       end if
    end do

    do i = 1, size(joints)
       if (joints(i)%samNodNum == node) then
          if (lDof < 1 .or. lDof > size(joints(i)%jointDofs)) exit
          entity = trim(dofSymbol_p(joints(i)%jointDofs(lDof)%lDof))//' in '// &
               &   trim(jointTypeName_p(joints(i)%type))//getId(joints(i)%id)
          if (present(writePos)) then
             if (writePos) then
                write(cVal,601) joints(i)%jointDofs(lDof)%jVar(1)
                entity(56:) = cVal
             end if
          end if
          return
       end if
    end do

    do i = 1, size(sups)
       if (associated(sups(i)%genDOFs)) then
          if (sups(i)%genDOFs%samNodNum == node) then
             if (lDof < 1 .or. lDof > sups(i)%genDOFs%nDOFs) exit
             write(cVal,"('Component mode',I3)") lDof
             entity = trim(cVal)//' in Part'//getId(sups(i)%id)
             return
          end if
       end if
    end do

    entity = '(invalid node/dof number)'

600 format(' (pos =',1P,3E12.4,'  Eul =',0P,3F7.1,')')
601 format(' (var =',1P,E12.4,')')

  end function nodeEntity


  !!============================================================================
  !> @brief Prints an error message on singular equation system.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Mar 2001

  subroutine reportEquationError (sam,sups,triads,joints,ierr,addMsg)

    use SamModule                 , only : SamType, nodeNumFromEqNum
    use SupElTypeModule           , only : SupElType
    use TriadTypeModule           , only : TriadType
    use MasterSlaveJointTypeModule, only : MasterSlaveJointType
    use ReportErrorModule         , only : reportError, error_p, empty_p

    type(SamType)             , intent(in) :: sam
    type(SupElType)           , intent(in) :: sups(:)
    type(TriadType)           , intent(in) :: triads(:)
    type(MasterSlaveJointType), intent(in) :: joints(:)
    integer                   , intent(in) :: ierr
    character(len=*),optional , intent(in) :: addMsg

    !! Local variables
    integer            :: i, ldof, node
    character(len=128) :: errMsg

    !! --- Logic section ---

    if (ierr > -3 .or. ierr < -5) return

    if (size(meqErr) > 1 .and. ierr == -3) then
       if (meqErr(1) <= 0_ik) return

       !! All singular equations were found in one go
       call reportError (error_p,'Singular or ill-conditioned system matrix'// &
            &                    ' detected at the following DOF(s):')
       do i = 1, size(meqErr)
          ldof = int(meqErr(i))
          if (ldof <= 0) exit
          node = nodeNumFromEqNum(sam,ldof)
          if (node <= 0) exit
          errMsg = nodeEntity(node,ldof,sups,triads,joints)
          write(errMsg(41:),102) meqErr(i), sam%minex(node), ldof
          call reportError (empty_p,errMsg)
       end do
       if (present(addMsg)) call reportError (empty_p,addMsg)
       return

    end if

    ldof = int(meqErr(1))
    node = nodeNumFromEqNum(sam,ldof)

    if (node <= 0) then
       return
    else if (ierr == -3) then
       write(errMsg,103) meqErr(1),sam%minex(node),ldof
    else if (ierr == -4) then
       write(errMsg,104) meqErr(1),sam%minex(node),ldof
    else
       write(errMsg,105) meqErr(1),sam%minex(node),ldof
    end if

    call reportError (error_p,trim(errMsg)//' '// &
         &            nodeEntity(node,ldof,sups,triads,joints),addMsg)

102 format('  ( Equation', I6, ', Node', I6, ', dof', I3, ' )')
103 format('Ill-conditioned Equation', I6, ', Node', I6, ', dof', I3, ',')
104 format('Zero pivot in Equation', I6, ', Node', I6, ', dof', I3, ',')
105 format('System matrix is not positive definite, Equation', I6, &
         & ', Node', I6, ', dof ', I6, ',')

  end subroutine reportEquationError


  !!============================================================================
  !> @brief Prints a warning message on negative pivots in the equation system.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 7 Dec 2005

  subroutine reportNegativePivots (sam,sups,triads,joints,numneg,time,addMsg)

    use SamModule                 , only : SamType, nodeNumFromEqNum
    use SupElTypeModule           , only : SupElType
    use TriadTypeModule           , only : TriadType, dp
    use MasterSlaveJointTypeModule, only : MasterSlaveJointType
    use ReportErrorModule         , only : reportError, warning_p, empty_p
    use ReportErrorModule         , only : getErrorFile

    type(SamType)             , intent(in) :: sam
    type(SupElType)           , intent(in) :: sups(:)
    type(TriadType)           , intent(in) :: triads(:)
    type(MasterSlaveJointType), intent(in) :: joints(:)
    integer         , optional, intent(in) :: numneg
    real(dp)        , optional, intent(in) :: time
    character(len=*), optional, intent(in) :: addMsg

    !! Local variables
    integer            :: i, ldof, node, lpu
    integer, parameter :: maxPrint_p = 20
    character(len=128) :: errMsg

    !! --- Logic section ---

    if (meqErr(1) <= 0_ik) return

    if (present(numneg) .and. present(time)) then
       write(errMsg,100) numneg,time
    else if (present(numneg)) then
       write(errMsg,101) numneg
    else
       errMsg = 'Negative pivot(s) in the system matrix were detected'
    end if
    call reportError (warning_p,adjustl(errMsg))
    lpu = getErrorFile()
    do i = 1, min(maxPrint_p,size(meqErr))
       ldof = int(meqErr(i))
       if (ldof <= 0) exit
       node = nodeNumFromEqNum(sam,ldof)
       if (node <= 0) exit
       errMsg = nodeEntity(node,ldof,sups,triads,joints)
       write(lpu,102) errMsg(1:40), meqErr(i), sam%minex(node), ldof
    end do
    if (numneg > maxPrint_p) write(lpu,103)
    if (present(addMsg)) call reportError (empty_p,addMsg)

100 format(I8,' negative pivot(s) in the system matrix were detected,', &
         & ' at time =',1PE12.5)
101 format(I8,' negative pivot(s) in the system matrix were detected')
102 format(10X,A,'  ( Equation', I6, ', Node', I6, ', dof', I3, ' )')
103 format(10X,'...')

  end subroutine reportNegativePivots
  !> @endcond

  !!============================================================================
  !> @brief Calculates all iteration norms defined in the system.
  !>
  !> @param[in]  sam Data for managing system matrix assembly
  !> @param      sys System level model data
  !> @param[in]  iter Iteration counter
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !> @date 12 Apr 2003
  !>
  !> @author Knut Morten Okstad
  !> @date 29 Jan 2004

  subroutine CalculateIterationNorms (sam,sys,iter,ierr)

    use SamModule         , only : SamType, dp
    use SystemTypeModule  , only : SystemType
    use NormTypeModule    , only : iVecNorm_p, iInfTra_p, iInfRot_p, iInfGen_p
    use NormTypeModule    , only : checkDivergence
    use NormRoutinesModule, only : InfiniteNorms, EnergyNorms, ScaledNorm
    use ReportErrorModule , only : debugFileOnly_p, reportError

    type(SamType)   , intent(in)    :: sam
    type(SystemType), intent(inout) :: sys
    integer         , intent(in)    :: iter
    integer         , intent(out)   :: ierr

    !! Local variables
    logical  :: findWD
    real(dp) :: dN, vN, aN, rN

    !! --- Logic section ---

    findWD = iter > sys%convergenceSet%startMonitor ! Find the worst DOFs?

    !! Calculate the scaled vector norms

    dN = ScaledNorm(sam,sys%sinc(:,1),sys%wDisp,2,ierr=ierr)
    if (ierr < 0) goto 900
    vN = ScaledNorm(sam,sys%sinc(:,2),sys%wDisp,2,ierr=ierr)
    if (ierr < 0) goto 900
    aN = ScaledNorm(sam,sys%sinc(:,3),sys%wDisp,2,ierr=ierr)
    if (ierr < 0) goto 900
    rN = ScaledNorm(sam,sys%residual,sys%wForce,2,.true.,ierr)
    if (ierr < 0) goto 900

    if (iter > 1 .and. .not.findWD) then
       !! Check for possibly starting divergence and find the worst DOFs if so
       call checkDivergence (sys%convergenceSet%disNorms(iVecNorm_p),dN,findWD)
       call checkDivergence (sys%convergenceSet%velNorms(iVecNorm_p),vN,findWD)
       call checkDivergence (sys%convergenceSet%accNorms(iVecNorm_p),aN,findWD)
       call checkDivergence (sys%convergenceSet%resNorms(iVecNorm_p),rN,findWD)
    end if

    sys%convergenceSet%disNorms(iVecNorm_p)%value = dN
    sys%convergenceSet%velNorms(iVecNorm_p)%value = vN
    sys%convergenceSet%accNorms(iVecNorm_p)%value = aN
    sys%convergenceSet%resNorms(iVecNorm_p)%value = rN
    sys%convergenceSet%doingWell = .not.findWD

    !! Calculate the infinite norms and find the worst DOFs, if requested

    call InfiniteNorms (sam, sys%sinc(:,1), findWD, &  ! Displacement correction
         &              sys%convergenceSet%disNorms(iInfTra_p)%value, &
         &              sys%convergenceSet%disNorms(iInfRot_p)%value, &
         &              sys%convergenceSet%disNorms(iInfGen_p)%value, &
         &              sys%convergenceSet%worstDOFs(:,1), ierr=ierr)
    if (ierr < 0) goto 900

    call InfiniteNorms (sam, sys%sinc(:,2), .false., & ! Velocity correction
         &              sys%convergenceSet%velNorms(iInfTra_p)%value, &
         &              sys%convergenceSet%velNorms(iInfRot_p)%value, &
         &              sys%convergenceSet%velNorms(iInfGen_p)%value, &
         &              ierr=ierr)
    if (ierr < 0) goto 900

    call InfiniteNorms (sam, sys%sinc(:,3), .false., & ! Acceleration correction
         &              sys%convergenceSet%accNorms(iInfTra_p)%value, &
         &              sys%convergenceSet%accNorms(iInfRot_p)%value, &
         &              sys%convergenceSet%accNorms(iInfGen_p)%value, &
         &              ierr=ierr)
    if (ierr < 0) goto 900

    call InfiniteNorms (sam, sys%residual, findWD, &   ! Force residual
         &              sys%convergenceSet%resNorms(iInfTra_p)%value, &
         &              sys%convergenceSet%resNorms(iInfRot_p)%value, &
         &              sys%convergenceSet%resNorms(iInfGen_p)%value, &
         &              sys%convergenceSet%worstDOFs(:,2), .true., ierr=ierr)
    if (ierr < 0) goto 900

    !! Calculate the energy norms and find the worst energy DOFs, if requested

    call EnergyNorms (sam, sys%sinc(:,1), sys%residual, findWD, &
         &            sys%convergenceSet%energyNorms(1)%value, &
         &            sys%convergenceSet%energyNorms(2)%value, &
         &            sys%convergenceSet%worstDOFs(:,3), &
         &            sys%convergenceSet%worstEnerg, ierr)
    if (ierr == 0) return

900 call reportError (debugFileOnly_p,'CalculateIterationNorms')

  end subroutine CalculateIterationNorms


  !!============================================================================
  !> @brief Calculates all iteration tolerances defined in the system.
  !>
  !> @param[in] sam Data for managing system matrix assembly
  !> @param     sys System level model data
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 12 Apr 2003

  subroutine CalculateIterationTolerances (sam,sys)

    use SamModule         , only : SamType
    use SystemTypeModule  , only : SystemType
    use NormTypeModule    , only : nNormTypes_p, iVecNorm_p
    use NormRoutinesModule, only : CalcTolerance

    type(SamType)   , intent(in)    :: sam
    type(SystemType), intent(inout) :: sys

    !! Local variables
    integer :: i

    !! --- Logic section ---

    !! Set the (possibly) velocity proportional tolerances
    !!TODO,bh: Should most of these be velocity proportional ????

    do i = 1, nNormTypes_p
       call CalcTolerance (sys%convergenceSet%disNorms(i))
       call CalcTolerance (sys%convergenceSet%velNorms(i))
       call CalcTolerance (sys%convergenceSet%accNorms(i))
       call CalcTolerance (sys%convergenceSet%resNorms(i))
    end do

    !! Recalculate for velocity scaled norm since this can be velocity
    !! proportional (for historic reasons)

    call CalcTolerance (sys%convergenceSet%velNorms(iVecNorm_p), &
         &              sam,sys%urd,sys%wDisp)

    !! Energy tolerances

    call CalcTolerance (sys%convergenceSet%energyNorms(1))
    call CalcTolerance (sys%convergenceSet%energyNorms(2))

  end subroutine CalculateIterationTolerances


  !!============================================================================
  !> @brief Accelerates the solution increment in case of line search.
  !>
  !> @param[in]  iter Iteration counter
  !> @param[in]  res Current force residual
  !> @param      del Current solution increment
  !> @param[out] ierr Error flag
  !>
  !> @details The solution increment is scaled if the previous and current
  !> residual vectors are "almost" tangential, as are previous and current
  !> solution vectors (kind-of line search).
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 29 Mar 2002

  subroutine IterationAccelerator (iter,res,del,ierr)

    use reportErrorModule, only : allocationError
#ifdef bjorn
    use fileUtilitiesModule, only : getDBGfile
#endif

    integer , intent(in)    :: iter
    real(dp), intent(in)    :: res(:)
    real(dp), intent(inout) :: del(:)
    integer , intent(out)   :: ierr

    !! Local variables
    logical, save :: firstCall = .true.
    integer  :: n
    real(dp), allocatable, save :: resPrev(:), resPar(:), resOrt(:)
    real(dp), allocatable, save :: delPrev(:), delPar(:), delOrt(:)
    real(dp) :: fResPar, fResOrt, fDelPar, fDelOrt, factor
    real(dp) :: lRes, lDel, lResPrev, lDelPrev
    real(dp), parameter :: eps = 1.0e-20

#ifdef bjorn
    integer :: dbg2
    dbg2 = getDBGfile(2)
#endif

    !! --- Logic section ---

    if (firstCall) then
       firstCall = .false.
       n = size(res)
       allocate(resPrev(n),resPar(n), resOrt(n), &
            &   delPrev(n),delPar(n), delOrt(n), STAT=ierr)
       if (ierr /= 0) then
          ierr = allocationError('IterationAccelerator')
          return
       end if
    end if

    ierr = 0
    if (iter < 1)         goto 100

    lRes     = sqrt(dot_product(res,res))
    if ( lRes < eps )     goto 100

    lResPrev = sqrt(dot_product(resPrev,resPrev))
    if ( lResPrev < eps ) goto 100

    lDel     = sqrt(dot_product(del,del))
    if ( lDel < eps )     goto 100

    lDelPrev = sqrt(dot_product(delPrev,delPrev))
    if ( lDelPrev < eps ) goto 100

    resPar   = (dot_product(resPrev,res)/(lResPrev**2))*resPrev
    resOrt   = res - resPar
    fResPar  = sqrt(dot_product(resPar,resPar))/lResPrev
    fResOrt  = sqrt(dot_product(resOrt,resOrt))/lResPrev

    delPar   = (dot_product(delPrev,del)/(lDelPrev**2))*delPrev
    delOrt   = del - delPar
    fDelPar  = sqrt(dot_product(delPar,delPar))/lDelPrev
    fDelOrt  = sqrt(dot_product(delOrt,delOrt))/lDelPrev

#ifdef bjorn
    write(dbg2,"('  -- Residual: lPar/lPrev, lOrt/lPrev :',1p 2e11.2)") fResPar,fResOrt
    write(dbg2,"('  -- Increm  : lPar/lPrev, lOrt/lPrev :',1p 2e11.2)") fDelPar,fDelOrt
#endif

    if ( lRes > 0.1_dp*lResPrev ) then  !Slow convergence
#ifdef bjorn
       write(dbg2,"('    -> Cond. 1 OK: lRes > 0.1*lResPrev ')")
#endif
       if ( fResPar > 10.0_dp*fResOrt ) then !Residual vector going in same dir as prev
#ifdef bjorn
          write(dbg2,"('    --> Cond. 2 OK: fResPar > 10.0*fResOrt ')")
#endif
          if ( fDelPar > 10.0_dp*fDelOrt ) then !Increment vector going in same dir as prev
#ifdef bjorn
             write(dbg2,"('    ---> Cond. 3 OK: fDelPar > 10.0*fDelOrt ')")
#endif
             factor = lResPrev/(lResPrev-lRes)
             if ( abs(factor) > 5.0_dp) factor = sign(5.0_dp,factor)

             del = factor * delPar + delOrt
#ifdef bjorn
             write(dbg2,"('    ----> ===== Acceleration with factor :',1pe12.5)") factor
#endif
          end if
       end if
    end if

100 continue
    resPrev = res
    delPrev = del

  end subroutine IterationAccelerator


  !!============================================================================
  !> @brief Saves response variables of current time step to file.
  !>
  !> @param[in]  sys System level model data
  !> @param      mech Mechanism components of the model
  !> @param[in]  ctrl Control system data
  !> @param[out] ierr Error flag
  !> @param[in]  checkTimeStep If present and .false., save the secondary
  !> variables for each time step (ignoring the -saveinc2 option).
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 4 Jun 2002

  subroutine saveStep (sys,mech,ctrl,ierr,checkTimeStep)

    use SystemTypeModule         , only : SystemType, dp
    use MechanismTypeModule      , only : MechanismType
    use ControlTypeModule        , only : ControlType
    use ControlTypeModule        , only : hasControlElements, writeCtrlSysDB
    use TireRoutinesModule       , only : scaleTo
    use EngineRoutinesModule     , only : UpdateEnginesForSave
    use HydrodynamicsModule      , only : getWaveElevation
    use WindTurbineRoutinesModule, only : getWindSpeed, getTipWindSpeed
    use WindTurbineRoutinesModule, only : getHubWindSpeed, getBladeDeflections
    use FiniteElementModule      , only : updateBeamSectionForces
    use FiniteElementModule      , only : writeBeamJointForces
    use SaveModule               , only : writeSolverDB1, writeSolverDB2, ctrlDB
    use SaveNastranModule        , only : writeSolverDB3
    use SaveVTFModule            , only : writeSolverVTF
    use ProfilerModule           , only : sav_p, startTimer, stopTimer
    use ReportErrorModule        , only : debugFileOnly_p, reportError
    use FFaCmdLineArgInterface   , only : ffa_cmdlinearg_getbool
    use FFaCmdLineArgInterface   , only : ffa_cmdlinearg_getdouble

    type(SystemType)   , intent(in)    :: sys
    type(MechanismType), intent(inout) :: mech
    type(ControlType)  , intent(in)    :: ctrl
    integer            , intent(out)   :: ierr
    logical, optional  , intent(in)    :: checkTimeStep

    !! Local variables
    integer, parameter :: maxWP_p = 20
    integer  :: i, nBlade, nWP
    logical  :: lDouble, dumpJoints
    real(dp) :: saveStart, saveInc, saveInc2, epsT, wave, wind(6,maxWP_p)

    !! --- Logic section ---

    epsT = 0.1_dp*sys%minInc
    call ffa_cmdlinearg_getdouble ('savestart',saveStart)
    if (sys%time+epsT < saveStart) return ! No save before saveStart

    call startTimer (sav_p)

    !! Always save primary response variables
    lasTime(1) = sys%time
    call ffa_cmdlinearg_getbool ('double1',lDouble)
    call writeSolverDB1 (sys,mech%triads,mech%sups,lDouble,ierr)
    if (ierr < 0) goto 900

    !! Write response data to the VTF-file
    call writeSolverVTF (sys,mech%sups,ierr)
    if (ierr < 0) goto 900

    saveInc2 = 0.0_dp
    if (present(checkTimeStep)) then
       !! Check if the secondary variables should be saved at current time step
       if (checkTimeStep) then
          call ffa_cmdlinearg_getdouble ('saveinc2',saveInc2)
          if (.not. saveNext(lasTime(2),saveInc2)) goto 100
       end if
    end if

    call updateBeamSectionForces (mech%sups,ierr)
    if (ierr < 0) goto 900

    call updateEnginesForSave (mech%engines,2,ierr)
    if (ierr < 0) goto 900

    if (associated(mech%env%waveFunc)) then
       wind = 0.0_dp
       wave = getWaveElevation(mech%env,wind(1:3,1),sys%time,ierr)
       if (ierr < 0) goto 900
    end if

    if (associated(mech%turbine)) then
       nBlade = size(mech%turbine%blade)
       nWP = 1 + nBlade + size(mech%turbine%windPt,2)
       call getBladeDeflections (mech%turbine)
       call getHubWindSpeed (sys%time,mech%turbine,wind(1:3,1),wind(4:6,1),ierr)
       if (ierr < 0) goto 900
       do i = 2, min(1+nBlade,maxWP_p)
          call getTipWindSpeed (i-1,sys%time,mech%turbine, &
               &                wind(1:3,i),wind(4:6,i),ierr)
          if (ierr < 0) goto 900
       end do
       nWP = min(1+nBlade+size(mech%turbine%windPt,2),maxWP_p)
       do i = 2+nBlade, nWP
          call getWindSpeed (mech%turbine%windPt(:,i-1-nBlade),sys%time, &
               &             wind(1:3,i),wind(4:6,i),ierr)
          if (ierr < 0) goto 900
       end do
    else
       nWP = 1
    end if

    !! Save secondary variables at current time step
    call ffa_cmdlinearg_getbool ('double2',lDouble)
    call writeSolverDB2 (sys,mech,ctrl,scaleTo,wave,wind(:,1:nWP),lDouble,ierr)
    if (ierr < 0) goto 900

    !! Write beam joint forces for external fatigue analysis
    call ffa_cmdlinearg_getbool ('saveBeamJoints',dumpJoints)
    if (dumpJoints) call writeBeamJointForces (sys%time,ierr)

    !! Check if variables for external stress recovery
    !! should be saved at current time step
100 call ffa_cmdlinearg_getdouble ('saveinc3',saveInc)
    if (saveNext(lasTime(3),saveInc)) then
       !! Save results for external stress recovery, etc.
       call writeSolverDB3 (sys,mech,ierr)
       if (ierr < 0) goto 900
    end if

    if (hasControlElements(ctrl) .and. ctrl%saveVar) then
       !! Check if control system should be saved at current time step
       call ffa_cmdlinearg_getdouble ('saveinc4',saveInc)
       if (abs(saveInc) < 1.0e-12_dp) saveInc = saveInc2 ! Default use saveInc2
       if (saveNext(lasTime(4),saveInc)) then
          !! Save control system data
          call writeCtrlSysDB (ctrlDB,ctrl,sys%nStep,sys%time,ierr)
          if (ierr < 0) goto 900
       end if
    end if

    goto 990

900 call reportError (debugFileOnly_p,'saveStep')
990 call stopTimer (sav_p)

  contains

    !> @brief Convenience function checking if next step should be saved or not.
    logical function saveNext (lastTime,saveSkip)
      real(dp), intent(inout) :: lastTime
      real(dp), intent(in)    :: saveSkip
      saveNext = .false.
      if (saveSkip < 0.0_dp) then
         return
      else if (lastTime > -1.0e99_dp) then
         if (sys%time+epsT < lastTime+saveSkip) return
      else
         if (sys%time+epsT < saveStart) return
      end if
      saveNext = .true.
      lastTime = sys%time
    end function saveNext

  end subroutine saveStep


  !!============================================================================
  !> @brief Terminates a time step by doing some final calculations.
  !>
  !> @param[in]  sam Data for managing system matrix assembly
  !> @param      sys System level model data
  !> @param      mech Mechanism components of the model
  !> @param[in]  ctrl Control system data
  !> @param[in]  linearStatic If .true., we are doing a linear static analysis
  !> @param[in]  finalStep If .true., we have finished the final time step
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine performs some final calculations on the converged
  !> solution state, updates the previous state variables, etc.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date Sep 2000

  subroutine terminateStep (sam,sys,mech,ctrl,linearStatic,finalStep,ierr)

    use KindModule                , only : i8
    use SamModule                 , only : SamType
    use SystemTypeModule          , only : SystemType
    use MechanismTypeModule       , only : MechanismType
    use ControlTypeModule         , only : ControlType
    use IdTypeModule              , only : StrId
    use SpringTypeModule          , only : checkUnLoading, checkFailure
    use SpringTypeModule          , only : updateAtConvergence
    use SupElTypeModule           , only : updateAtConvergence
    use HydroDynamicsModule       , only : updateAtConvergence
    use MasterSlaveJointTypeModule, only : updateAtConvergence
    use ContactElementTypeModule  , only : updateAtConvergence
    use FrictionTypeModule        , only : updateAtConvergence
    use TimeStepModule            , only : pushAccelStack
    use RotationModule            , only : orthonorm3
    use ReportErrorModule         , only : debugFileOnly_p, reportError
    use FFaCmdLineArgInterface    , only : ffa_cmdlinearg_getint

    type(SamType)      , intent(in)    :: sam
    type(SystemType)   , intent(inout) :: sys
    type(MechanismType), intent(inout) :: mech
    type(ControlType)  , intent(in)    :: ctrl
    logical            , intent(in)    :: linearStatic, finalStep
    integer            , intent(out)   :: ierr

    !! Local variables
    integer :: i, iprint

    !! --- Logic section ---

    ierr = 0
    do i = 1, size(mech%triads)
       if (mech%triads(i)%nDOFs == 6) then
          call orthonorm3 (mech%triads(i)%ur(:,1:3))
       end if
    end do

    do i = 1, size(mech%cElems)
       call updateAtConvergence (mech%cElems(i))
    end do

    !! If failure criterion has been exceeded, set hasFailed for next timestep
    do i = 1, size(mech%baseSprings)
       call checkUnloading (mech%baseSprings(i),ierr)
       call checkFailure (mech%baseSprings(i))
    end do
    if (ierr < 0) goto 900

    do i = 1, size(mech%joints)
       if (associated(mech%joints(i)%springEl)) then
          call checkFailure (mech%joints(i)%springEl)
       end if
    end do

    if (sys%varInc > 0) then
       call pushAccelStack (sam,sys%urdd,ierr)
       if (ierr < 0) goto 900
    end if

    call ffa_cmdlinearg_getint ('debug',iprint)
    if (finalStep .and. iprint > 0) iprint = -int(sys%nStep)
    if (sys%nStep > 0_i8 .and. iprint == -int(sys%nStep)) then
       !! Debug output, to compare against in restart, etc.
       call dumpStep (sys,mech,ctrl,sam, &
            &         fName='datastructure.s'//adjustl(StrId(int(sys%nStep))))
    end if

    !! Now update all previous state variables for the next time step.
    !! We do these tasks here such that the dumpStep call above correctly
    !! prints the "previous" values.

    !! Skip this in the final step and when doing linear statics
    if (.not. (linearStatic .or. finalStep)) then
       call updatePreviousState (sam,sys,mech,ierr)
       if (ierr < 0) goto 900
    end if

    return

900 call reportError (debugFileOnly_p,'terminateStep')

  end subroutine terminateStep


  !!============================================================================
  !> @brief Updates the previous state variables in the mechanism objects.
  !>
  !> @param[in] sam Data for managing system matrix assembly
  !> @param sys System level model data
  !> @param mech Mechanism components of the model
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 11 Aug 2023

  subroutine updatePreviousState (sam,sys,mech,ierr)

    use SamModule                 , only : SamType
    use SystemTypeModule          , only : SystemType
    use MechanismTypeModule       , only : MechanismType
    use SpringTypeModule          , only : updateAtConvergence
    use SupElTypeModule           , only : updateAtConvergence
    use HydroDynamicsModule       , only : updateAtConvergence
    use MasterSlaveJointTypeModule, only : updateAtConvergence
    use ContactElementTypeModule  , only : updateAtConvergence
    use FrictionTypeModule        , only : updateAtConvergence
    use ReportErrorModule         , only : debugFileOnly_p, reportError

    type(SamType)      , intent(in)    :: sam
    type(SystemType)   , intent(inout) :: sys
    type(MechanismType), intent(inout) :: mech
    integer            , intent(out)   :: ierr

    !! Local variables
    integer :: i

    !! --- Logic section ---

    ierr = 0

    call DCOPY (sam%ndof,sys%urd(1) ,1,sys%urdp(1) ,1)
    call DCOPY (sam%ndof,sys%urdd(1),1,sys%urddp(1),1)

    sys%nIterPrevStep = real(sys%nIterThisStep,dp)

    do i = 1, size(mech%cElems)
       mech%cElems(i)%wasActive = mech%cElems(i)%isActive
       if (mech%cElems(i)%isActive) then
          mech%cElems(i)%cVarPrev = mech%cElems(i)%cVar(1:3,1)
       end if
    end do

    do i = 1, size(mech%triads)
       if (mech%triads(i)%nDOFs > 0) then
          mech%triads(i)%urPrev = mech%triads(i)%ur
       end if
    end do

    do i = 1, size(mech%sups)
       call updateAtConvergence (mech%sups(i))
       if (associated(mech%sups(i)%hydyn)) then
          call updateAtConvergence (mech%sups(i)%hydyn,mech%sups(i)%supTr,ierr)
          if (ierr < 0) goto 900
       end if
    end do

    do i = 1, size(mech%elms)
       if (associated(mech%elms(i)%hydyn)) then
          call updateAtConvergence (mech%elms(i)%hydyn,mech%elms(i)%Tlg,ierr)
          if (ierr < 0) goto 900
       end if
    end do

    do i = 1, size(mech%joints)
       call updateAtConvergence (mech%joints(i))
    end do

    do i = 1, size(mech%baseSprings)
       call updateAtConvergence (mech%baseSprings(i))
    end do

    do i = 1, size(mech%frictions)
       if (associated(mech%frictions(i)%p)) then
          call updateAtConvergence (mech%frictions(i)%p)
       end if
    end do

    return

900 call reportError (debugFileOnly_p,'updatePreviousState')

  end subroutine updatePreviousState


  !!============================================================================
  !> @brief Restores all state variables from the last converged time step.
  !>
  !> @param[in] sam Data for managing system matrix assembly
  !> @param sys System level model data
  !> @param mech Mechanism components of the model
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date Nov 2008

  subroutine restoreLastStep (sam,sys,mech)

    use SamModule                 , only : SamType
    use SystemTypeModule          , only : SystemType
    use MechanismTypeModule       , only : MechanismType
    use ControlTypeModule         , only : ControlType
    use SupElTypeModule           , only : restoreFromLastStep
    use MasterSlaveJointTypeModule, only : restoreFromLastStep
    use SpringTypeModule          , only : restoreFromLastStep
    use FrictionTypeModule        , only : restoreFromLastStep

    type(SamType)      , intent(in)    :: sam
    type(SystemType)   , intent(inout) :: sys
    type(MechanismType), intent(inout) :: mech

    !! Local variables
    integer :: i

    !! --- Logic section ---

    call DCOPY (sam%ndof,sys%urdp(1) ,1,sys%urd(1) ,1)
    call DCOPY (sam%ndof,sys%urddp(1),1,sys%urdd(1),1)

    do i = 1, size(mech%cElems)
       mech%cElems(i)%cVar(1:3,1) = mech%cElems(i)%cVarPrev
       mech%cElems(i)%isActive = mech%cElems(i)%wasActive
    end do

    do i = 1, size(mech%triads)
       if (mech%triads(i)%nDOFs > 0) then
          mech%triads(i)%ur = mech%triads(i)%urPrev
       end if
    end do

    do i = 1, size(mech%sups)
       call restoreFromLastStep (mech%sups(i))
    end do

    do i = 1, size(mech%joints)
       call restoreFromLastStep (mech%joints(i))
    end do

    do i = 1, size(mech%baseSprings)
       call restoreFromLastStep (mech%baseSprings(i))
    end do

    do i = 1, size(mech%frictions)
       if (associated(mech%frictions(i)%p)) then
          call restoreFromLastStep (mech%frictions(i)%p)
       end if
    end do

  end subroutine restoreLastStep


  !!============================================================================
  !> @brief Clears velocity and acceleration variables in the triads and joints.
  !>
  !> @param     triads All triads in the model
  !> @param     joints All joints in the model
  !> @param[in] notify If .true., give message if a non-zero variable is cleared
  !>
  !> @details This subroutine will (re)set the velocity- and accelerations
  !> variables in all triads and joints to zero. It mainly is used to nullify
  !> potentially non-zero initial conditions in the case of quasi-static load
  !> incrementation, in which inertia and damping terms should be zero, such
  !> that a possible dynamic restart will be performed from a resting state.
  !> Notice that the velocity/accelerations of the generalized DOFs in super-
  !> elements do not need to be reset, since those variables are pointers into
  !> respective system vectors which are assumed to be initialized separately.
  !>
  !> @callergraph
  !>
  !> @author Guenter Glanzer
  !>
  !> @date 24 Aug 2022

  subroutine clearVelAcc (triads,joints,notify)

    use TriadTypeModule           , only : TriadType
    use MasterSlaveJointTypeModule, only : MasterSlaveJointType

    type(TriadType)           , intent(inout) :: triads(:)
    type(MasterSlaveJointType), intent(inout) :: joints(:)
    logical                   , intent(in)    :: notify

    !! Local variables
    integer :: i, j
    logical :: nulled(2)

    !! --- Logic section ---

    do i = 1, size(triads)
       if (associated(triads(i)%urd)) then
          if (notify .and. any(nonZero(triads(i)%urd))) then
             call notifyClear ('Triad',triads(i)%id,'velocity')
          end if
          triads(i)%urd = 0.0_dp
       end if
       if (associated(triads(i)%urdd)) then
          if (notify .and. any(nonZero(triads(i)%urdd))) then
             call notifyClear ('Triad',triads(i)%id,'acceleration')
          end if
          triads(i)%urdd = 0.0_dp
       end if
    end do

    do i = 1, size(joints)
       nulled = .false.
       do j = 1, joints(i)%nJointDOFs
          if (notify) then
             nulled = nulled .or. nonZero(joints(i)%jointDOFs(j)%jVar(2:3))
          end if
          joints(i)%jointDOFs(j)%jVar(2:3) = 0.0_dp
       end do
       if (nulled(1)) then
          call notifyClear ('Joint',joints(i)%id,'velocity')
       end if
       if (nulled(2)) then
          call notifyClear ('Joint',joints(i)%id,'acceleration')
       end if
    end do

  contains

    !> @brief Checks if a value is nonzero or not with some tolerance.
    elemental logical function nonZero (value)
      real(dp), intent(in) :: value
      nonZero = abs(value) > 1.0e-15_dp
    end function nonZero

    !> @brief Notifies on clearing velocity or acceleration variables.
    subroutine notifyClear (objType,objId,var)
      use IdTypeModule     , only : IdType, getId
      use reportErrorModule, only : reportError, noteFileOnly_p
      character(len=*), intent(in) :: objType, var
      type(IdType)    , intent(in) :: objId
      call reportError (noteFileOnly_p,'Clearing non-zero '// var //' in '// &
           &            objType // getId(objId))
    end subroutine notifyClear

  end subroutine clearVelAcc


  !!============================================================================
  !> @brief Closes all solver data base files before program termination.
  !>
  !> @param      mech Mechanism components of the model
  !> @param[out] anyRes Indicates whether the result files contain data or not
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 30 Jan 2006

  subroutine closeSolverDB (mech,anyRes,ierr)

    use MechanismTypeModule, only : MechanismType
    use SaveVTFModule      , only : closeVTF
    use SaveModule         , only : closeResultFiles
    use TireRoutinesModule , only : closeTireModules
    use ReportErrorModule  , only : debugFileOnly_p, reportError

    type(MechanismType), intent(inout) :: mech
    logical            , intent(out)   :: anyRes(:)
    integer            , intent(out)   :: ierr

    !! --- Logic section ---

    lasTime = -1.0e99_dp ! In case solver is started again by the same process

    call closeVTF (size(mech%sups),ierr)
    if (ierr < 0) goto 900

    call closeResultFiles (anyRes,ierr)
    call closeTireModules (mech%tires)
    if (ierr >= 0) return

900 call reportError (debugFileOnly_p,'closeSolverDB')

  end subroutine closeSolverDB


  !!============================================================================
  !> @brief Prints a summary of the mass distribution to the res-file.
  !>
  !> @param[in]  sups All superelements in the model
  !> @param      masses All point masses in the model
  !> @param[in]  tires All tires in the model
  !> @param[in]  elms All user-defined elements in the model
  !> @param[in]  gravity Gravitation vector
  !> @param[out] totMass Total mass of the model
  !> @param[in]  lpu File unit number for res-file output
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 17 Apr 2002

  subroutine printMassDistribution (sups,masses,tires,elms,gravity,totMass,lpu)

    use KindModule                  , only : epsDiv0_p
    use SupElTypeModule             , only : SupElType
    use MassTypeModule              , only : MassType
    use TireTypeModule              , only : TireType
    use UserDefElTypeModule         , only : UserdefElType
    use FiUserElmInterface          , only : Fi_UDE6
    use FiniteElementModule         , only : massGrowth, massIntFluid
    use MassRoutinesModule          , only : updateMasses
    use IdTypeModule                , only : getId, getAssId
    use ManipMatrixModule           , only : matmul34
    use ScratchArrayModule          , only : getRealScratchArray
    use ReportErrorModule           , only : reportError, noteFileOnly_p
    use FFaTensorTransformsInterface, only : traTensor, traInertia

    type(SupElType)    , intent(in)    :: sups(:)
    type(MassType)     , intent(inout) :: masses(:)
    type(TireType)     , intent(in)    :: tires(:)
    type(UserdefElType), intent(in)    :: elms(:)
    real(dp)           , intent(in)    :: gravity(:)
    real(dp)           , intent(out)   :: totMass
    integer            , intent(in)    :: lpu

    !! Local variables
    integer, parameter :: maxAss_p = 100
    character(len=32)  :: assId(maxAss_p)
    real(dp)           :: assMass(4,maxAss_p)
    real(dp), pointer  :: Xnod(:)
    logical            :: massFlag(size(masses))

    integer  :: i, j, k, trId, nAssemb, nenod
    real(dp) :: m, totAddM, totUDEM, linkMass, totCG(3), linkCG(3)
    real(dp) :: evec(3), inertia(6), totInertia(6), gdir(3)

    !! --- Logic section ---

    m = gravity(1)*gravity(1) + gravity(2)*gravity(2) + gravity(3)*gravity(3)
    if (m < epsDiv0_p) then
       gdir = (/ 0.0_dp, 0.0_dp, 1.0_dp /)
    else
       gdir = gravity/sqrt(m)
    end if

    write(lpu,600)

    if (size(masses) > 0) then
       !! Flags to ensure the additional masses are counted only once
       !! In case the triads are connected to more than one superelement
       massFlag = .false.
       !! Update additional masses, in case of time-dependency and/or scaling
       call updateMasses (masses,gravity,.true.)
    end if

    nAssemb = 0
    totMass = 0.0_dp
    totCG = 0.0_dp
    do i = 1, size(sups)
       linkMass = sups(i)%mass
       linkCG   = matmul34(sups(i)%suptr,sups(i)%posMassCenter)
       totMass  = totMass + linkMass
       totCG    = totCG   + linkMass*linkCG
       if (sups(i)%addedMass) then
          !! Sum up the additional mass for this superelement
          linkCG = linkCG*linkMass
          do k = 1, size(masses)
             if (.not. massFlag(k)) then
                trId = masses(k)%triad%id%baseId
                do j = 1, sups(i)%nExtNods
                   if (sups(i)%triads(j)%p%id%baseId == trId) then
                      m = getMass(masses(k),gDir)
                      linkMass = linkMass + m
                      linkCG   = linkCG   + m*sups(i)%triads(j)%p%ur(:,4)
                      massFlag(k) = .true.
                      exit
                   end if
                end do
             end if
          end do
          if (linkMass > epsDiv0_p) linkCG = linkCG / linkMass
          write(lpu,601) getId(sups(i)%id),linkCG,linkMass,sups(i)%mass
       else
          write(lpu,601) getId(sups(i)%id),linkCG,linkMass
       end if
       call getAssMass (getAssId(sups(i)%id),linkMass,linkCG)
    end do

    totUDEM = 0.0_dp
    do i = 1, size(elms)
       nenod = size(elms(i)%triads)
       Xnod => getRealScratchArray(3*nenod,j)
       if (j < 0) cycle
       linkCG = 0.0_dp
       do j = 1, nenod
          evec = elms(i)%triads(j)%p%ur(:,4)
          Xnod(3*j-2:3*j) = evec
          linkCG = linkCG + evec/real(nenod,dp)
       end do
       call Fi_UDE6 (elms(i)%id%baseId,elms(i)%type,nenod,Xnod(1), &
            &        elms(i)%iwork(1),elms(i)%rwork(1),linkMass,j)
       if (linkMass > 0.0_dp .and. j >= 0) then
          totUDEM  = totUDEM + linkMass
          totMass  = totMass + linkMass
          totCG    = totCG   + linkMass*linkCG
          write(lpu,601) getId(elms(i)%id),linkCG,linkMass
          call getAssMass (getAssId(elms(i)%id),linkMass,linkCG)
       end if
    end do

    totAddM = 0.0_dp
    do k = 1, size(masses)
       m = getMass(masses(k),gDir)
       totAddM = totAddM + m
       totCG   = totCG   + m*masses(k)%triad%ur(:,4)
    end do
    do k = 1, size(tires)
       totAddM = totAddM + tires(k)%mass
       linkCG  = tires(k)%RimTriad%ur(:,4) &
            &  + tires(k)%RimTriad%ur(:,3)*tires(k)%ZoffSet
       totCG   = totCG + linkCG*tires(k)%mass
       call getAssMass (getAssId(tires(k)%id),tires(k)%mass,linkCG)
    end do

    if (totMass+totAddM > epsDiv0_p) totCG = totCG / (totMass+totAddM)

    do i = nAssemb, 1, -1
       k = index(assId(i),',')
       if (k > 0) then
          do j = 1, i-1
             if (assId(j) == assId(i)(k+1:)) then
                assMass(:,j) = assMass(:,j) + assMass(:,i)
                exit
             end if
          end do
       end if
    end do
    do j = 1, nAssemb
       if (assMass(4,j) > epsDiv0_p) then
          assMass(1:3,j) = assMass(1:3,j) / assMass(4,j)
          write(lpu,601) 'Sub-assembly '//assId(j),assMass(:,j)
       end if
    end do

    !! Compute the total inertia of the mechanism

    totInertia = 0.0_dp
    do i = 1, size(sups)
       evec = matmul34(sups(i)%suptr,sups(i)%posMassCenter) - totCG
       call mat2tensor (sups(i)%inertia(4:6,4:6),inertia)
       call traTensor  (3,inertia,sups(i)%suptr)
       call traInertia (3,inertia,evec,sups(i)%mass)
       totInertia = totInertia + inertia
    end do

    do k = 1, size(masses)
       !! Account for possibly direction-dependent mass
       evec = matmul(masses(k)%triad%ur(:,4)-totCG, &
            &        masses(k)%triad%ur(:,1:3)) * sqrt(masses(k)%mass)
       call mat2tensor (masses(k)%II,inertia)
       call traInertia (3,inertia,evec,1.0_dp)
       call traTensor  (3,inertia,masses(k)%triad%ur)
       totInertia = totInertia + inertia
    end do

    do k = 1, size(tires)
       totAddM = totAddM + tires(k)%mass
       evec = tires(k)%RimTriad%ur(:,4) &
          & + tires(k)%RimTriad%ur(:,3)*tires(k)%ZoffSet - totCG
       inertia(1) = tires(k)%Iyy
       inertia(2) = tires(k)%Iyy
       inertia(3) = tires(k)%Izz
       inertia(4:6) = 0.0_dp
       call traTensor  (3,inertia,tires(k)%RimTriad%ur)
       call traInertia (3,inertia,evec,tires(k)%mass)
       totInertia = totInertia + inertia
    end do

    totMass = totMass - massIntFluid - massGrowth
    totAddM = totAddM + massIntFluid + massGrowth
    if (totAddM > 0.0_dp) then
       write(lpu,602) totCG,totMass+totAddM,totMass
       if (massIntFluid > 0.0_dp) write(lpu,603) massIntFluid
       if (massGrowth   > 0.0_dp) write(lpu,604) massGrowth
       totMass = totMass + totAddM
    else
       write(lpu,602) totCG,totMass
    end if
    if (totUDEM > 0.0_dp) then
       write(lpu,605) totUDEM
       call reportError (noteFileOnly_p,'The global inertia reported below '// &
            &            'is excluding masses from the user-defined elements')
    end if
    write(lpu,610) totInertia

600 format( 5X,'Part/Element',20X,'Centre of Gravity', &
         & 23X,'Mass',9X,'excl. add. mass')
601 format( 5X,A30,1P,3E13.5,1X,2E13.5)
602 format( 5X,'Whole mechanism',15X,1P,3E13.5,1X,2E13.5)
603 format( 5X,'Mass of internal fluid',48X,1PE13.5)
604 format( 5X,'Mass of marine growth ',48X,1PE13.5)
605 format( 5X,'Mass of all user-defined elements',37X,1PE13.5/)
610 format(/5X,'Global inertia about Centre of Gravity' &
         / 10X,'Ixx =',1PE13.5, &
         / 10X,'Iyy =',1PE13.5, &
         / 10X,'Izz =',1PE13.5, &
         / 10X,'Ixy =',1PE13.5, &
         / 10X,'Ixz =',1PE13.5, &
         / 10X,'Iyz =',1PE13.5 //)

  contains

    !> @brief Extracts equivalent mass in a given direction for a mass element.
    function getMass (mass,gdir)
      use TriadTypeModule, only : transVGlobToTriad
      type(MassType), intent(in) :: mass
      real(dp)      , intent(in) :: gdir(3)
      real(dp) :: tdir(3), getMass
      !! Direction vector in local triad system tdir = {t1,t2,t3}
      tdir = transVGlobToTriad(mass%triad,gdir)
      !! Mass = {t1*mx,t2*my,t3*mz}*{t1,t2,t3}
      getMass = dot_product(mass%mass*tdir,tdir) ! assuming tdir has unit length
    end function getMass

    !> @brief Converts a 3x3 matrix into a second-order tensor.
    !> @details The matrix is assumed symmetric.
    !> The tensor storage scheme is: T11, T22, T33, T12, T13, T23.
    !> @sa FFaTensor3.
    subroutine mat2tensor (M,T)
      real(dp), intent(in)  :: M(3,3)
      real(dp), intent(out) :: T(6)
      T(1) = M(1,1)
      T(2) = M(2,2)
      T(3) = M(3,3)
      T(4) = M(1,2)
      T(5) = M(1,3)
      T(6) = M(2,3)
    end subroutine mat2tensor

    !> @brief Sums up masses and CG for sub-assemblies.
    subroutine getAssMass (assemb,mass,CG)
      character(len=*), intent(in) :: assemb
      real(dp)        , intent(in) :: mass, CG(3)
      if (assemb == ' ') return
      do j = 1, nAssemb
         if (assemb == assId(j)) then
            assMass(4,j)   = assMass(4,j) + mass
            assMass(1:3,j) = assMass(1:3,j) + mass*CG
            return
         end if
      end do
      if (nAssemb < size(assId)) then
         nAssemb = nAssemb + 1
         assId(nAssemb) = assemb
         assMass(4,nAssemb) = mass
         assMass(1:3,nAssemb) = mass*CG
      end if
    end subroutine getAssMass

  end subroutine printMassDistribution


  !!============================================================================
  !> @brief Prints basic simulation variables to the res-file.
  !>
  !> @param[in] time Current simulation time
  !> @param[in] triads All triads in the model
  !> @param[in] sups All superelements in the model
  !> @param[in] iprint Print switch; the higher value the more print is produced
  !> @param[in] lpu File unit number for res-file output
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 18 Jun 2001

  subroutine printResults (time,triads,sups,iprint,lpu)

    use TriadTypeModule       , only : TriadType
    use SupElTypeModule       , only : SupElType, dp
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getints

    real(dp)       , intent(in) :: time
    type(TriadType), intent(in) :: triads(:)
    type(SupElType), intent(in) :: sups(:)
    integer        , intent(in) :: iprint, lpu

    !! Local variables
    integer :: i, j, k, triadId(10)
    logical :: lFirst

    !! --- Logic section ---

    call ffa_cmdlinearg_getints ('printTriad',triadId,10)
    if (iprint > 0 .or. any(triadId > 0)) then

       write(lpu,600) time
       do i = 1, size(triads)
          if (.not. printThis(triads(i)%id)) cycle

          if ( .not. associated(triads(i)%nodeForce) .and. &
               .not. associated(triads(i)%supForce) ) then
             if (triads(i)%nDOFs == 6) then
                write(lpu,612) triads(i)%id%baseId, &
                     &       ((triads(i)%ur(k,j),j=1,4), &
                     &         triads(i)%urd(k),triads(i)%urd(3+k), &
                     &         triads(i)%urdd(k),triads(i)%urdd(3+k), k=1,3)
             else if (triads(i)%nDOFs == 3) then
                write(lpu,613) triads(i)%id%baseId, &
                     &       ((triads(i)%ur(k,j),j=1,4), &
                     &         triads(i)%urd(k),triads(i)%urdd(k), k=1,3)
             end if
          else if ( associated(triads(i)%supForce) .and. &
               &    .not. associated(triads(i)%nodeForce) ) then
             if (triads(i)%nDOFs == 6) then
                write(lpu,610) triads(i)%id%baseId, &
                     &       ((triads(i)%ur(k,j),j=1,4), &
                     &         triads(i)%urd(k),triads(i)%urd(3+k), &
                     &         triads(i)%urdd(k),triads(i)%urdd(3+k), &
                     &         triads(i)%supForce(k), &
                     &         triads(i)%supForce(3+k), k=1,3)
             else if (triads(i)%nDOFs == 3) then
                write(lpu,611) triads(i)%id%baseId, &
                     &       ((triads(i)%ur(k,j),j=1,4), &
                     &         triads(i)%urd(k),triads(i)%urdd(k), &
                     &         triads(i)%supForce(k), k=1,3)
             end if
          else if ( .not. associated(triads(i)%supForce) .or. &
               &    associated(triads(i)%nodeForce,triads(i)%supForce) ) then
             if (triads(i)%nDOFs == 6) then
                write(lpu,610) triads(i)%id%baseId, &
                     &       ((triads(i)%ur(k,j),j=1,4), &
                     &         triads(i)%urd(k),triads(i)%urd(3+k), &
                     &         triads(i)%urdd(k),triads(i)%urdd(3+k), &
                     &         triads(i)%nodeForce(k), &
                     &         triads(i)%nodeForce(3+k), k=1,3)
             else if (triads(i)%nDOFs == 3) then
                write(lpu,611) triads(i)%id%baseId, &
                     &       ((triads(i)%ur(k,j),j=1,4), &
                     &         triads(i)%urd(k),triads(i)%urdd(k), &
                     &         triads(i)%nodeForce(k), k=1,3)
             end if
          else ! we have both nodeForce and supForce !=> nodeForce
             if (triads(i)%nDOFs == 6) then
                write(lpu,614) triads(i)%id%baseId, &
                     &       ((triads(i)%ur(k,j),j=1,4), &
                     &         triads(i)%urd(k),triads(i)%urd(3+k), &
                     &         triads(i)%urdd(k),triads(i)%urdd(3+k), &
                     &         triads(i)%supForce(k), &
                     &         triads(i)%supForce(3+k), &
                     &         triads(i)%nodeForce(k), &
                     &         triads(i)%nodeForce(3+k), k=1,3)
             else if (triads(i)%nDOFs == 3) then
                write(lpu,615) triads(i)%id%baseId, &
                     &       ((triads(i)%ur(k,j),j=1,4), &
                     &         triads(i)%urd(k),triads(i)%urdd(k), &
                     &         triads(i)%supForce(k), &
                     &         triads(i)%nodeForce(k), k=1,3)
             end if

          end if
       end do

    end if
    if (iprint > 0) then

       lFirst = .true.
       do i = 1, size(sups)
          if (associated(sups(i)%genDOFs)) then
             if (lFirst) write(lpu,620)
             write(lpu,621) sups(i)%id%baseId, &
                  &        (sups(i)%genDOFs%ur(k), &
                  &         sups(i)%genDOFs%urd(k), &
                  &         sups(i)%genDOFs%urdd(k), k=1,sups(i)%genDOFs%nDOFs)
             lFirst = .false.
          end if
       end do
       if (.not.lFirst) write(lpu,622)

    end if

600 format(/140('*') /' Triad    T =', 1PE13.5, &
         '    ----------POSITION-----------   --------VELOCITY--------', &
         '   ------ACCELERATION------   ---------FORCE----------')
610 format(/I6,1X,1P,4E13.5,3(1X,2E13.5)     / (7X,1P,4E13.5,3(1X,2E13.5)))
611 format(/I6,1X,1P,4E13.5,3(1X, E13.5,13X) / (7X,1P,4E13.5,3(1X, E13.5,13X)))
612 format(/I6,1X,1P,4E13.5,2(1X,2E13.5)     / (7X,1P,4E13.5,2(1X,2E13.5)))
613 format(/I6,1X,1P,4E13.5,2(1X, E13.5,13X) / (7X,1P,4E13.5,2(1X, E13.5,13X)))
614 format(/I6,1X,1P,4E13.5,4(1X,2E13.5)     / (7X,1P,4E13.5,4(1X,2E13.5)))
615 format(/I6,1X,1P,4E13.5,4(1X, E13.5,13X) / (7X,1P,4E13.5,4(1X, E13.5,13X)))
620 format(/140('*') /' SupEl', &
         '     DISPLACEMENT      VELOCITY      ACCELERATION')
621 format(/I6,1X,1P,3(3X,E13.5) / (7X,1P,3(3X,E13.5)))
622 format(/140('*'))

  contains

    !> @brief Checks if a given triad should be printed or not.
    function printThis (tId)
      use IdTypeModule, only : IdType
      type(IdType), intent(in) :: tId
      logical :: printThis
      printThis = .true.
      if (iprint > 0) then
         return ! Debug mode, print out for all triads
      else if (.not. associated(tId%assId)) then ! Top-level triads only
         do j = 1, size(triadId)
            if (triadId(j) == tId%userId) then
               return ! Print this triad
            else if (triadId(j) <= 0) then
               exit
            end if
         end do
      end if
      printThis = .false.
    end function printThis

  end subroutine printResults


  !!============================================================================
  !> @brief Prints convergence information to the res-file.
  !>
  !> @param[in] lpu File unit number for res-file output
  !> @param[in] fileFormat File format option,
  !>            &lt; 0 : old format (R4.0-compatible),
  !>               = 0 : new format
  !>            &gt; 0 : new format and write heading first
  !> @param[in] iStep Time increment counter
  !> @param[in] iter Iteration counter
  !> @param[in] tanUpdate If .true., a new tangent is computed in this iteration
  !> @param[in] time Current simulation time
  !> @param[in] timeStep Time increment size
  !> @param[in] res Force residual
  !> @param[in] del Solution increment
  !> @param[in] convData Iteration norm tolerances and values
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[in] mech Mechanism components of the model
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 29 Jan 2004

  subroutine printConvergence (lpu,fileFormat,iStep,iter,tanUpdate, &
       &                       time,timeStep,del,res,convData,sam,mech)

    use KindModule         , only : i8
    use NormTypeModule     , only : TestSetType, TestItemType, iVecNorm_p
    use SamModule          , only : SamType, nodeNumFromDofNum
    use MechanismTypeModule, only : MechanismType

    integer            , intent(in) :: lpu, fileFormat, iter
    integer(i8)        , intent(in) :: iStep
    logical            , intent(in) :: tanUpdate
    real(dp)           , intent(in) :: time, timeStep, del(:), res(:)
    type(TestSetType)  , intent(in) :: convData
    type(SamType)      , intent(in) :: sam
    type(MechanismType), intent(in) :: mech

    !! Local variables
    integer :: i, ieq, lDof, node, nWorst

    !! --- Logic section ---

    nWorst = size(convData%worstDOFs,1)

    if (fileFormat < 0) then
       if (iStep > 9999999999_i8) return ! To avoid format overflow

       !! Old res-file format (retained mostly for historic reasons
       !! and to enable diff'ing with old res-files of previous releases)
       write(lpu,600,advance='NO') iStep,iter,time,timeStep
       if (convData%velNorms(iVecNorm_p)%relTol > 0.0_dp) then
          write(lpu,601,advance='NO') convData%velNorms(iVecNorm_p)%tolerance, &
               &                      convData%velNorms(iVecNorm_p)%value, &
               &                      convData%resNorms(iVecNorm_p)%value
       else
          write(lpu,602,advance='NO') convData%velNorms(iVecNorm_p)%value, &
               &                      convData%resNorms(iVecNorm_p)%value
       end if

    else

       !! New res-file format
       if (fileFormat > 0 .and. iter == 0) then
          write(lpu,'(/A)',advance='NO') '  TimeStep   Time / h    Iter '
          call writeNormHeadings (convData%disNorms)
          call writeNormHeadings (convData%velNorms)
          call writeNormHeadings (convData%accNorms)
          call writeNormHeadings (convData%resNorms)
          call writeNormHeadings (convData%energyNorms)
          write(lpu,'(/)')
       end if
       if (iter == 0 .or. (.not.convData%doingWell .and. nWorst > 0)) then
          if (iStep < 10000000000_i8) then
             write(lpu,603,advance='NO') iStep,time,iter
          else
             write(lpu,604,advance='NO') time,iter
          end if
       else if (iter == 1) then
          write(lpu,604,advance='NO') timeStep,iter
       else
          write(lpu,605,advance='NO') iter
       end if
       call writeNormValues (convData%disNorms)
       call writeNormValues (convData%velNorms)
       call writeNormValues (convData%accNorms)
       call writeNormValues (convData%resNorms)
       call writeNormValues (convData%energyNorms)

    end if

    if (tanUpdate) then
       write(lpu,"('')")
    else
       write(lpu,'(A)') ' (no tangent update)'
    end if

    if (convData%doingWell) return ! We are converging nicely
    if (nWorst < 1) return ! Suppress all convergence problem output

    !! Report convergence problems encountered

    write(lpu,610)
    if (convData%isActive(1)) then
       !! Write out the worst displacement DOFs
       write(lpu,620) 'displacement'
       do i = 1, nWorst
          lDof = convData%worstDOFs(i,1)
          ieq  = sam%meqn(lDof)
          node = nodeNumFromDofNum(sam,lDof)
          write(lpu,630) convData%worstDOFs(i,1), del(ieq), &
           trim(nodeEntity(node,lDof,mech%sups,mech%triads,mech%joints,.true.))
       end do
    end if

    if (convData%isActive(2)) then
       !! Write out the worst residual DOFs
       write(lpu,620) 'force residual'
       do i = 1, nWorst
          lDof = convData%worstDOFs(i,2)
          ieq  = sam%meqn(lDof)
          node = nodeNumFromDofNum(sam,lDof)
          write(lpu,630) convData%worstDOFs(i,2), res(ieq), &
           trim(nodeEntity(node,lDof,mech%sups,mech%triads,mech%joints,.true.))
       end do
    end if

    if (convData%isActive(3)) then
       !! Write out the worst energy DOFs
       write(lpu,620) 'energy'
       do i = 1, nWorst
          lDof = convData%worstDOFs(i,3)
          node = nodeNumFromDofNum(sam,lDof)
          write(lpu,630) convData%worstDOFs(i,3), convData%worstEnerg(i), &
           trim(nodeEntity(node,lDof,mech%sups,mech%triads,mech%joints,.true.))
       end do
    end if
    write(lpu,"('')")

600 format('     NSTP =', I10, '  IT =', I3, '  T =', 1P,E12.5, '  H =', E12.5 )
601 format('  ERRM =', 1PE12.5, '  DELN =', E12.5, '  RESN =', E12.5 )
602 format('  DELN =', 1PE12.5, '  RESN =', E12.5 )
603 format(I10,1PE13.5,I6,1X)
604 format(10X,1PE13.5,I6,1X)
605 format(        23X,I6,1X)
610 format('            Slow convergence or possible divergence is detected.')
620 format('              The worst ', A, ' DOF(s) in this iteration:')
630 format(I20,1PE13.5,2X,A)

  contains

    !> @brief Writes the heading of the iteration norms.
    subroutine writeNormHeadings (norms)
      type(TestItemType), intent(in) :: norms(:)
      do i = 1, size(norms)
         if (norms(i)%code > 0) then
            write(lpu,'(2X,A11)',advance='NO') norms(i)%title
         end if
      end do
    end subroutine writeNormHeadings

    !> @brief Writes the iteration norms values.
    subroutine writeNormValues (norms)
      type(TestItemType), intent(in) :: norms(:)
      character(len=12) :: cnorm
      do i = 1, size(norms)
         if (norms(i)%code > 0) then
            write(cnorm,"(1PE12.5)") norms(i)%value
            if (norms(i)%value < norms(i)%tolerance) then
               cnorm(1:1) = '*' ! Converged norm
            else if (norms(i)%nIterIncr > 0) then
               cnorm(1:1) = '#' ! Possibly diverging norm
            end if
            write(lpu,"(1X,A12)",advance='NO') cnorm
         end if
      end do
    end subroutine writeNormValues

  end subroutine printConvergence

  !> @cond FULL_DOC
  !!============================================================================
  !> @brief Debug print of the equation system iteration process.
  !>
  !> @param[in] lpu File unit number for res-file output
  !> @param[in] iter Iteration counter
  !> @param[in] time Physical time
  !> @param[in] res Force residual
  !> @param[in] del Solution increment
  !> @param[in] FIk Inertia forces
  !> @param[in] FDk Damping forces
  !> @param[in] FSk Stiffness forces
  !> @param[in] Qk External forces
  !> @param[in] FIprev Inertia forces from previous time step
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[in] mech Mechanism components of the model
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 29 Mar 2002

  subroutine printSysSolveProcess (lpu,iter,time,res,del, &
       &                           FIk,FDk,FSk,Qk,FIprev, &
       &                           sam,mech)

    use KindModule            , only : epsDiv0_p
    use SamModule             , only : SamType, nodeNumFromEqNum
    use MechanismTypeModule   , only : MechanismType
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getints
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_isTrue

    integer , intent(in) :: lpu, iter
    real(dp), intent(in) :: time, res(:), del(:)
    real(dp), intent(in) :: FIk(:), FDk(:), FSk(:), Qk(:), FIprev(:)

    type(SamType)      , intent(in) :: sam
    type(MechanismType), intent(in) :: mech

    !! Local variables
    integer  :: ieq, iDof, lDof, node, printSys(2)
    real(dp) :: maxInc, maxRes, relMax(2)

    !! --- Logic section ---

    call ffa_cmdlinearg_getints ('printSys',printSys,2)
    if (printSys(2) > size(del)) printSys(2) = size(del)
    if (printSys(2) < 1 .or. printSys(2) < printSys(1)) return

    maxRes = 0.0_dp
    maxInc = 0.0_dp
    do ieq = printSys(1), printSys(2)
       if (abs(res(ieq)) > abs(maxRes)) maxRes = res(ieq)
       if (abs(del(ieq)) > abs(maxInc)) maxInc = del(ieq)
    end do

    if (iter < 1) write(lpu,600) time
    if (ffa_cmdlinearg_isTrue('printSysDofs')) then
       write(lpu,601) iter, '  idof'
    else
       write(lpu,601) iter, '   ieq'
    end if
    do iDof = printSys(1), printSys(2)
       if (ffa_cmdlinearg_isTrue('printSysDofs')) then
          ieq = sam%meqn(iDof)
          if (ieq <= 0) cycle
       else
          ieq = iDof
       end if
       relMax = 0.0_dp
       if (maxRes > epsDiv0_p) relMax(1) = res(ieq)/maxRes
       if (maxInc > epsDiv0_p) relMax(2) = del(ieq)/maxInc
       lDof = ieq
       node = nodeNumFromEqNum(sam,lDof)
       write(lpu,610) iDof, res(ieq), del(ieq), relMax, &
            &         FIk(ieq), FDk(ieq), FSk(ieq), Qk(ieq), FIprev(ieq), &
            lDof, trim(nodeEntity(node,lDof,mech%sups,mech%triads,mech%joints))
    end do

    write(lpu,611) sqrt(dot_product(res,res)), &
         &         sqrt(dot_product(del,del)), &
         &         sqrt(dot_product(FIk,FIk)), &
         &         sqrt(dot_product(FDk,FDk)), &
         &         sqrt(dot_product(FSk,FSk)), &
         &         sqrt(dot_product(Qk,Qk)), &
         &         sqrt(dot_product(FIprev,FIprev))

600 format(/'=======================================' &
         & /'===== New time-step at Time =',1p,e10.3)
601 format(/'== residual and increment, it =',i3 &
         & / a6, '    residual     dispInc      %res      %del   ', &
         &  '   FI         FD         FS          Q        FIprev')
610 format(i6,1p,2e12.3,2p,2(f9.2,'%'),1p,5e11.2,'  <-dof',i2,1x,a,' -')
611 format('L2-Norm',1p,e11.3,e12.3,20X,5e11.2)

  end subroutine printSysSolveProcess


  !!============================================================================
  !> @brief Debug print of all current solver data to file.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 28 Oct 2008

  subroutine dumpStep (sys,mech,ctrl,sam,lpu,fName)

    use KindModule         , only : i8
    use SystemTypeModule   , only : SystemType, writeObject
    use MechanismTypeModule, only : MechanismType, writeObject
    use ControlTypeModule  , only : ControlType, writeObject
    use SamModule          , only : SamType, writeObject
    use FileUtilitiesModule, only : findUnitNumber

    type(SystemType)          , intent(in) :: sys
    type(MechanismType)       , intent(in) :: mech
    type(ControlType)         , intent(in) :: ctrl
    type(SamType)             , intent(in) :: sam
    integer         , optional, intent(in) :: lpu
    character(len=*), optional, intent(in) :: fName

    !! Local variables
    integer :: idata

    !! --- Logic section ---

    if (present(lpu)) then
       idata = lpu
    else if (present(fName)) then
       idata = findUnitNumber(10)
       open(idata,FILE=trim(fName))
    else
       return
    end if

    if (sys%nStep < 100000_i8) then
       write(idata,"('* Dump of solver data after step',I6,' *'/)") sys%nStep
    else
       write(idata,"('* Dump of solver data at time',1PE12.5,' *'/)") sys%time
    end if
    call writeObject (sys,idata,2)
    write (idata,'(/)')
    call writeObject (mech,idata,2)
    write (idata,'(/)')
    call writeObject (ctrl,idata,2)
    write (idata,'(/)')
    call writeObject (sam,idata,2)
    close(idata)

  end subroutine dumpStep
  !> @endcond

end module SolverRoutinesModule
