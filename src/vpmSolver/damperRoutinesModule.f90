!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file damperRoutinesModule.f90
!>
!> @brief Subroutines for damper calculations.

!!==============================================================================
!> @brief Module with subroutines for damper calculations.
!>
!> @details This module contains a set of subroutines for performing various
!> computation tasks on the dampertypemodule::dampertype objects in the model
!> (axial- and joint dampers) during the dynamic or quasi-static simulation.

module DamperRoutinesModule

  implicit none

  private :: calcDamperForces


contains

  !!============================================================================
  !> @brief Updates all axial- and joint dampers in the model.
  !>
  !> @param dampers All damper elements in the model
  !> @param[in] timeStep Current time increment size
  !> @param[in] restart If .true., we are initializing after restart
  !> @param ierr Error flag
  !> @param[in] updateLV If .true., update axial damper lengths and velocities
  !> @param[in] updateVar If .true., update all other damper variables
  !>
  !> @details The updated length and velocity (and direction) for axial dampers
  !> are calculated first. Then all the remaining damper variables are updated.
  !> This is to ensure that sensors measuring damper length/velocity are up to
  !> date before they are used by engines. Note that for joint dampers,
  !> the damper length (and velocity unless deformational damper) coincide
  !> with the corresponding joint variables of the owner joint,
  !> and need therefore not to be updated here.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !> @date 11 Jun 2002
  !>
  !> @author Bjorn Haugen
  !> @date 14 Jan 2004

  subroutine UpdateDampers (dampers,timeStep,restart,ierr,updateLV,updateVar)

    use DamperTypeModule    , only : DamperType, dp
    use DamperTypeModule    , only : updateDamperVelocity
    use EngineRoutinesModule, only : updateDamperBase
    use TriadTypeModule     , only : updateNodeForce
    use reportErrorModule   , only : reportError, debugFileOnly_p

    type(DamperType), intent(inout) :: dampers(:)
    real(dp)        , intent(in)    :: timeStep
    logical         , intent(in)    :: restart
    integer         , intent(inout) :: ierr
    logical,optional, intent(in)    :: updateLV, updateVar

    !! Local variables
    logical  :: doUpdate
    integer  :: i, lerr
    real(dp) :: eV(6)

    !! --- Logic section ---

    if (present(updateLV)) then
       doUpdate = updateLV
    else
       doUpdate = .true.
    end if
    if (doUpdate) then

       !! Update all axial damper lengths first in case some are used by engines
       do i = 1, size(dampers)
          if (dampers(i)%dmp%isActive .and. dampers(i)%dmp%dof == 0) then

             call updateLenVel (dampers(i)%dmp, &
                  &             dampers(i)%triad1,dampers(i)%triad2, &
                  &             dampers(i)%forceDir)

          end if

          !! Update deformational velocity, if used
          call updateDamperVelocity (dampers(i)%dmp,timeStep)
       end do

    end if

    if (present(updateVar)) then
       doUpdate = updateVar
    else
       doUpdate = .true.
    end if
    if (doUpdate) then
       lerr = ierr

       !! Now update all damper variables
       do i = 1, size(dampers)
          call updateDamperBase (dampers(i)%dmp,ierr)
          if (restart .or. dampers(i)%dmp%dof > 0 .or. ierr < lerr) cycle

          call calcDamperForces (dampers(i),eV)
          call updateNodeForce (dampers(i)%triad1,eV(1:3))
          call updateNodeForce (dampers(i)%triad2,eV(4:6))
       end do

       if (ierr < lerr) call reportError (debugFileOnly_p,'updateDampers')
    end if

  contains

    !> @brief Updates the length, velocity and direction for an axial damper.
    subroutine UpdateLenVel (dmp,triad1,triad2,fDir)

      use kindModule       , only : epsDiv0_p
      use DamperTypeModule , only : DamperBaseType, TriadType
      use IdTypeModule     , only : getId
      use reportErrorModule, only : warning_p

      type(DamperBaseType), intent(inout) :: dmp
      type(TriadType)     , intent(in)    :: triad1, triad2
      real(dp)            , intent(out)   :: fDir(:)

      real(dp) :: velocity(3)

      fDir = triad2%ur(:,4) - triad1%ur(:,4)

      dmp%length = sqrt(dot_product(fDir,fDir))
      if (dmp%length < epsDiv0_p) then
         call reportError (warning_p,'Axial damper'//getId(dmp%id), &
              'The two points defining the damper direction coincide.')
         fDir = 0.0_dp
      else
         fDir = fDir / dmp%length
      end if

      if (associated(dmp%spr)) return ! Skip deformational dampers

      if (triad1%nDOFs >= 3 .and. triad2%nDOFs >= 3) then
         velocity = triad2%urd(1:3) - triad1%urd(1:3)
      else if (triad1%nDOFs >= 3) then
         velocity = -triad1%urd(1:3)
      else if (triad2%nDOFs >= 3) then
         velocity = triad2%urd(1:3)
      else
         return
      end if
      dmp%velocity = dot_product(fDir,velocity)

    end subroutine UpdateLenVel

  end subroutine UpdateDampers


  !!============================================================================
  !> @brief Calculates system force vector contributions from a damper element.
  !>
  !> @param[in] damper The damper element to calculate damping forces for
  !> @param[out] eV Element damping force vector
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 18 Feb 2008

  subroutine calcDamperForces (damper,eV)

    use DamperTypeModule, only : DamperType, dp
    use TriadTypeModule , only : transGlobToSys

    type(DamperType), intent(in)  :: damper
    real(dp)        , intent(out) :: eV(:)

    !! Local variables
    integer :: dofs1

    !! --- Logic section ---

    eV = 0.0_dp
    if (.not. damper%dmp%isActive) return

    if (damper%dmp%dof == 0) then

       !! Axial damper
       if (size(eV) == 6) then
          dofs1 = 3
       else
          dofs1 = damper%triad1%nDOFs
       end if
       eV(1:3) = -damper%dmp%Force * damper%forceDir
       eV(dofs1+1:dofs1+3) = -eV(1:3)

       !! Transform the element vector to system directions
       call transGlobToSys (damper%triad1,eV(1:3))
       call transGlobToSys (damper%triad2,eV(dofs1+1:dofs1+3))

    else

       !! Joint damper
       dofs1     = damper%dmp%dof
       eV(dofs1) = damper%dmp%Force

    end if

  end subroutine calcDamperForces


  !!============================================================================
  !> @brief Assembles system force vector contributions from a damper element.
  !>
  !> @param FD System damping force vector
  !> @param RF System reaction forces associated with constrained DOFs
  !> @param[in] damper The damper element to calculate damping forces for
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Karl Erik Thoresen
  !>
  !> @date Feb 1999

  subroutine addInDamperForces (FD,RF,damper,sam,ierr)

    use SamModule         , only : SamType, dp
    use DamperTypeModule  , only : DamperType
    use AsmExtensionModule, only : csAddEV
#ifdef FT_DEBUG
    use IdTypeModule      , only : getId
    use manipMatrixModule , only : writeObject
    use dbgUnitsModule    , only : dbgSolve
#endif
    use reportErrorModule , only : reportError, debugFileOnly_p

    real(dp)        , intent(inout) :: FD(:), RF(:)
    type(DamperType), intent(in)    :: damper
    type(SamType)   , intent(in)    :: sam
    integer         , intent(inout) :: ierr

    !! Local variables
    integer  :: err, nDOFs
    real(dp) :: eV(12)

    !! --- Logic section ---

    if (.not. damper%dmp%isActive) return

    call calcDamperForces (damper,eV)

    nDOFs = damper%nDOFs
#ifdef FT_DEBUG
    if (dbgSolve > 0) then
       call writeObject (eV(1:nDOFs),dbgSolve,'Force vector'// &
            &            ' for Damper'//trim(getId(damper%dmp%id)))
    end if
#endif

    call csAddEV (sam, damper%samElnum, nDOFs, eV(1:nDOFs), FD, RF, err)
    if (err == 0) return

    ierr = ierr - 1
    call reportError (debugFileOnly_p,'addInDamperForces')

  end subroutine addInDamperForces


  !!============================================================================
  !> @brief Assembles system damping matrix contributions from a damper element.
  !>
  !> @param[in] includeStressStiff If .true., include geometric stiffness
  !> @param[in] scaleC Damping matrix scaling factor
  !> @param[in] scaleK Stiffness matrix scaling factor
  !> @param Nmat System newton matrix
  !> @param[in] damper The damper element to calculate damping matrix for
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[out] ierr Error flag
  !> @param Rhs System right-hand-side vector associated with the Newton matrix
  !>
  !> @details This subroutine calculates the damping matrix for a damper element
  !> and adds it into the system Newton matrix, multiplied by a factor, scaleC.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Karl Erik Thoresen
  !>
  !> @date Feb 1999

  subroutine addInDamperMat (includeStressStiff,scaleC,scaleK,Nmat, &
       &                     damper,sam,ierr,Rhs)

    use SamModule          , only : SamType, dp
    use SysMatrixTypeModule, only : SysMatrixType
    use DamperTypeModule   , only : DamperType
    use TriadTypeModule    , only : transGlobToSys
    use AsmExtensionModule , only : csAddElM
    use manipMatrixModule  , only : dyadic_product
#ifdef FT_DEBUG
    use IdTypeModule       , only : getId
    use manipMatrixModule  , only : writeObject
    use dbgUnitsModule     , only : dbgSolve
#endif
    use reportErrorModule  , only : reportError, debugFileOnly_p

    logical            , intent(in)    :: includeStressStiff
    real(dp)           , intent(in)    :: scaleC, scaleK
    type(SysMatrixType), intent(inout) :: Nmat
    type(DamperType)   , intent(in)    :: damper
    type(SamType)      , intent(in)    :: sam
    integer            , intent(out)   :: ierr
    real(dp), optional , intent(inout) :: Rhs(:)

    !! Local variables
    integer  :: j, dofs1, nDofs
    real(dp) :: Elmat(12,12), KgNom

    !! --- Logic section ---

    ierr = 0
    if (.not. damper%dmp%isActive) return

    elMat = 0.0_dp
    nDofs = damper%nDofs

    if (damper%dmp%dof == 0) then ! Axial Damper

       elMat(1:3,1:3) = dyadic_product(damper%forceDir) &
            &         * damper%dmp%coeff * scaleC

       if (includeStressStiff) then
          !! Add geometric stiffness
          KgNom = scaleK * damper%dmp%force / damper%dmp%length
          do j = 1, 3
             elMat(j,j) = elMat(j,j) + KgNom
          end do
          elMat(1:3,1:3) = elMat(1:3,1:3) &
               &         - KgNom * dyadic_product(damper%forceDir)
       end if

       dofs1 = damper%triad1%nDofs
       elMat(dofs1+1 : dofs1+3 , dofs1+1 : dofs1+3) =  elMat(1:3,1:3)
       elMat(      1 :       3 , dofs1+1 : dofs1+3) = -elMat(1:3,1:3)
       elMat(dofs1+1 : dofs1+3 ,       1 :       3) = -elMat(1:3,1:3)

       !! Transform the element matrix to system directions
       call transGlobToSys (damper%triad1, elMat, 1, 3)
       call transGlobToSys (damper%triad2, elMat, dofs1+1, 4)

    else ! Joint Damper

       dofs1 = damper%dmp%dof
       elMat(dofs1,dofs1) = damper%dmp%coeff * scaleC

    end if

#ifdef FT_DEBUG
    if (dbgSolve > 0) then
       call writeObject (elMat(1:nDOFs,1:nDOFs),dbgSolve,'Damping matrix'// &
            &            ' for Damper'//trim(getId(damper%dmp%id)))
    end if
#endif

    call csAddElM (sam, damper%samElnum, nDofs, &
         &         elMat(1:nDofs,1:nDofs), Nmat, ierr, Rhs)
    if (ierr < 0) call reportError (debugFileOnly_p,'addInDamperMat')

  end subroutine addInDamperMat

end module DamperRoutinesModule
