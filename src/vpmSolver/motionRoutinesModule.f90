!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file motionRoutinesModule.f90
!> @brief Subroutines for prescribed motion calculations.

!!==============================================================================
!> @brief Module with subroutines for prescribed motion calculations.
!>
!> @details This module contains a set of subroutines for updating the state-
!> dependent prescribed motion objects.

module MotionRoutinesModule

  implicit none

  private

  !> @brief Updates the prescribed motions.
  interface updatePrescribedMotions
     module procedure updatePrescribedMotionsA
     module procedure updatePrescribedMotionsB
  end interface

  public :: updatePrescribedMotions, initPrescribedVelAcc


contains

  !!============================================================================
  !> @brief Updates the prescribed motions.
  !>
  !> @param[in] t Current simulation time
  !> @param[in] h Time increment size
  !> @param[in] beta Newmark time integration parameter
  !> @param[in] gamma Newmark time integration parameter
  !> @param[in] sysVel System velocity vector
  !> @param[in] sysAcc System acceleration vector
  !> @param motions Array of all prescribed motions in the model
  !> @param[in] linearStatic If .true., multi-load case linear static analysis
  !> @param[in] lpu File unit number for debug output of motions from file
  !> @param[out] ierr Error flag
  !>
  !> @details A prescribed motion can be updated either by evaluating the
  !> general function connected to it, or by reading new values from the binary
  !> file provided. The latter approach is typically used in sub-model analysis.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date Oct 2004

  subroutine updatePrescribedMotionsA (t,h,beta,gamma,sysVel,sysAcc, &
       &                               motions,linearStatic,ierr,lpu)

    use MotionTypeModule      , only : MotionType, dp
    use MotionTypeModule      , only : prescribedDeflection_p
    use MotionTypeModule      , only : prescribedVelocity_p
    use MotionTypeModule      , only : prescribedAcceleration_p
    use IdTypeModule          , only : getId, StrId
    use TriadTypeModule       , only : getSysValue
    use EngineRoutinesModule  , only : Evaluate
    use PrescribedMotionModule, only : readMotionFile
    use PrescribedMotionModule, only : pDisFile, nextTime, pDispl
    use reportErrorModule     , only : internalError, reportError
    use reportErrorModule     , only : error_p, debugFileOnly_p
#ifdef FT_DEBUG
    use dbgUnitsModule        , only : dbgSolve
#endif

    real(dp)        , intent(in)    :: t, h, beta, gamma, sysVel(:), sysAcc(:)
    type(MotionType), intent(inout) :: motions(:)
    logical         , intent(in)    :: linearStatic
    integer,optional, intent(in)    :: lpu
    integer         , intent(out)   :: ierr

    !! Local variables
    integer  :: i, j, k, lerr, last
    logical  :: newData
    real(dp) :: delta, dPrev, vel, acc

    !! --- Logic section ---

    ierr = 1
    newData = .false.
    do while (ierr > 0 .and. t + 1.0e-12_dp > nextTime)
       !! Read prescribed displacements for next time step
       call readMotionFile (ierr,lpu)
       if (ierr > 0) then
          newData = .true.
       else if (ierr < 0) then
          call reportError (debugFileOnly_p,'updatePrescribedMotions')
          return
       end if
    end do

    ierr = 0
    lerr = 0
    last = 0
    do i = 1, size(motions)

       Dprev = motions(i)%D
       if (motions(i)%ipd == 0) then
          !! Get the function value for this motion in current state
          motions(i)%D = Evaluate(motions(i)%engine, &
               &                  motions(i)%d1,motions(i)%d0,ierr)
       else if (newData) then
          !! Update motion value from binary file
          motions(i)%D = pDispl(motions(i)%ipd)
       end if
       if (ierr < lerr) then
          lerr = ierr
          call reportError (error_p,'Failed to evaluate Function for '// &
               &                    'prescribed motion'//getId(motions(i)%id))
          cycle
       else if (newData .and. present(lpu)) then
          if (lpu > 0 .and. associated(motions(i)%triad)) then
             if (motions(i)%triad%id%baseId /= last) then
                if (last > 0) write(lpu,"()")
                write(lpu,600,advance='NO') motions(i)%triad%id%userId
                last = motions(i)%triad%id%baseId
             end if
             write(lpu,601,advance='NO') motions(i)%dof,motions(i)%D
600          format(5X,'Triad Id',I6,' :')
601          format(' u',I1,' = ',1PE12.5)
          end if
       end if

       !! Get velocity and acceleration in system directions, if needed
       if (motions(i)%type > 0 .and. associated(motions(i)%triad)) then
          j = motions(i)%triad%sysDOF
          k = j + motions(i)%triad%nDOFs-1
          vel = getSysValue(motions(i)%triad,motions(i)%dof,sysVel(j:k))
          acc = getSysValue(motions(i)%triad,motions(i)%dof,sysAcc(j:k))
       else
          vel = sysVel(motions(i)%sDof)
          acc = sysAcc(motions(i)%sDof)
       end if

       delta = motions(i)%D - Dprev
       select case (motions(i)%type)

       case (prescribedDeflection_p)
          if (linearStatic) then
             !! Apply total displacement state
             motions(i)%cc = motions(i)%D
          else
             !! Apply displacement increment when doing nonlinear analysis
             motions(i)%cc = delta
          end if

       case (prescribedVelocity_p)
          motions(i)%cc = (vel + acc*h*(0.5_dp-beta/gamma) + delta*beta/gamma)*h

       case (prescribedAcceleration_p)
          motions(i)%cc = (vel + acc*h*0.5_dp + delta*h*beta)*h

       case default
          ierr = internalError('updatePrescribedMotions: Invalid motion type')
          return

       end select

#ifdef FT_DEBUG
       if (dbgSolve <= 0) cycle
       if (i == 1) write(dbgSolve,"(/'   --- in updatePrescribedMotions')")
       write(dbgSolve,"(7X,'Motion',A,', DOF =',I6,' :')") &
            trim(getId(motions(i)%id)), motions(i)%sDof
       write(dbgSolve,"(7X,'D, Dprev, cc = ',1P3E13.5)") &
            motions(i)%D, Dprev, motions(i)%cc
       if (motions(i)%type > 0) then
          write(dbgSolve,"(7X,'vel, acc, dt = ',1P2E13.5,0P3F8.4)") &
               vel, acc, h, beta, gamma
       end if
#endif

    end do
    if (present(lpu) .and. last > 0) then
       write(lpu,"()")
    end if

    if (pDisFile < 0 .and. allocated(pDispl)) then
       deallocate(pDispl)
    end if

    if (ierr < 0) call reportError (debugFileOnly_p,'updatePrescribedMotions')

  end subroutine updatePrescribedMotionsA


  !!============================================================================
  !> @brief Updates the prescribed motions.
  !>
  !> @param motions Array of all prescribed motions in the model
  !> @param[out] ierr Error flag
  !>
  !> @details This version only (re)evaluates the prescribed motion functions,
  !> and is typically used after the initial configuration has been defined
  !> in order to ensure that all functions depending on other state variables
  !> are up to date before starting the time integration.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date Sep 2010

  subroutine updatePrescribedMotionsB (motions,ierr)

    use MotionTypeModule    , only : MotionType
    use IdTypeModule        , only : getId
    use EngineRoutinesModule, only : Evaluate
    use reportErrorModule   , only : reportError, error_p, debugFileOnly_p

    type(MotionType), intent(inout) :: motions(:)
    integer         , intent(out)   :: ierr

    !! Local variables
    integer :: i, lerr

    !! --- Logic section ---

    ierr = 0
    lerr = 0
    do i = 1, size(motions)

       !! Get the function value for this motion in current state
       motions(i)%D = Evaluate(motions(i)%engine, &
            &                  motions(i)%d1,motions(i)%d0,ierr)
       if (ierr < lerr) then
          lerr = ierr
          call reportError (error_p,'Failed to evaluate Function for '// &
               &                    'prescribed motion'//getId(motions(i)%id))
       end if

    end do

    if (ierr < 0) call reportError (debugFileOnly_p,'updatePrescribedMotions')

  end subroutine updatePrescribedMotionsB


  !!==========================================================================
  !> @brief Inserts prescribed motion values into the system vectors.
  !>
  !> @param[in] motions Array of all prescribed motions in the model
  !> @param sysVel System velocity vector
  !> @param sysAcc System acceleration vector
  !>
  !> @details Only for prescribed velocities and accelerations.
  !> Does nothing for prescribed displacements.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date Oct 2004

  subroutine initPrescribedVelAcc (motions,sysVel,sysAcc)

    use MotionTypeModule, only : MotionType, dp
    use MotionTypeModule, only : prescribedVelocity_p, prescribedAcceleration_p

    type(MotionType), intent(in)    :: motions(:)
    real(dp)        , intent(inout) :: sysVel(:), sysAcc(:)

    !! Local variables
    integer :: i, idof

    !! --- Logic section ---

    do i = 1, size(motions)
       select case (motions(i)%type)

       case (prescribedVelocity_p)
          if (associated(motions(i)%triad)) then
             call initPrescribed (motions(i)%triad,motions(i)%D, &
                  &               motions(i)%sDof,motions(i)%dof,sysVel)
          else
             sysVel(motions(i)%sDof) = motions(i)%D
          end if

       case (prescribedAcceleration_p)
          if (associated(motions(i)%triad)) then
             call initPrescribed (motions(i)%triad,motions(i)%D, &
                  &               motions(i)%sDof,motions(i)%dof,sysAcc)
          else
             sysAcc(motions(i)%sDof) = motions(i)%D
          end if

       end select
    end do

  contains

    !> @brief Inserts motion values for triad DOFs into a system vector.
    subroutine initPrescribed (triad,D,sysDof,locDof,sysVec)
      use TriadTypeModule, only : TriadType, hasLocalDirections
      use TriadTypeModule, only : transGlobToSys, transSysToGlob
      type(TriadType), intent(in)    :: triad
      real(dp)       , intent(in)    :: D
      integer        , intent(in)    :: sysDof, locDof
      real(dp)       , intent(inout) :: sysVec(:)
      if (hasLocalDirections(triad)) then
         idof = sysDof - locDof + 1
         if (locDof > 3) idof = idof + 3
         call transGlobToSys (triad,sysVec(idof:idof+2))
         sysVec(sysDof) = D
         call transSysToGlob (triad,sysVec(idof:idof+2))
      else
         sysVec(sysDof) = D
      end if
    end subroutine initPrescribed

  end subroutine initPrescribedVelAcc

end module MotionRoutinesModule
