!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file massRoutinesModule.f90
!>
!> @brief Subroutines for discrete mass calculations.

!!==============================================================================
!> @brief Module with subroutines for calculation of discrete mass elements.
!>
!> @details This module contains a set of subroutines for performing various
!> computation tasks on the masstypemodule::masstype objects in the model
!> during the dynamic or quasi-static simulation.

module MassRoutinesModule

  implicit none

  private

  public :: updateMasses
  public :: addInGravitationForces, addInInertiaForces, addInMassMat


contains

  !!============================================================================
  !> @brief Updates all mass elements.
  !>
  !> @param masses All mass elements of the model
  !> @param[in] gravity Gravitation vector
  !> @param[in] restart If .true., skip force calculation
  !> @param ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 30 May 2002

  subroutine updateMasses (masses,gravity,restart,ierr)

    use MassTypeModule   , only : MassType, dp
    use reportErrorModule, only : reportError, debugFileOnly_p

    type(MassType)   , intent(inout) :: masses(:)
    real(dp)         , intent(in)    :: gravity(3)
    logical          , intent(in)    :: restart
    integer, optional, intent(inout) :: ierr

    !! Local variables
    integer :: i, lerr

    !! --- Logic section ---

    if (present(ierr)) then
       lerr = ierr
    else
       lerr = 0
    end if

    do i = 1, size(masses)
       if (masses(i)%triad%nDOFs > 0) then
          call updateMass (masses(i),gravity,restart,ierr)
       end if
    end do

    if (present(ierr)) then
       if (ierr < lerr) call reportError (debugFileOnly_p,'updateMasses')
    end if

  end subroutine updateMasses


  !!============================================================================
  !> @brief Updates a mass element.
  !>
  !> @param mass The mass element to update
  !> @param[in] gravity Gravitation vector
  !> @param[in] restart If .true., skip force calculation
  !> @param ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 30 Jan 2008

  subroutine updateMass (mass,gravity,restart,ierr)

    use MassTypeModule      , only : MassType, dp
    use TriadTypeModule     , only : updateNodeForce
    use EngineRoutinesModule, only : EngineValue
    use reportErrorModule   , only : reportError, debugFileOnly_p

    type(MassType)   , intent(inout) :: mass
    real(dp)         , intent(in)    :: gravity(3)
    logical          , intent(in)    :: restart
    integer, optional, intent(inout) :: ierr

    !! Local variables
    integer  :: lerr
    real(dp) :: eV(6)

    !! --- Logic section ---

    lerr = 0

    if (associated(mass%massEngine)) then
       mass%mass = mass%m0 + mass%m1*EngineValue(mass%massEngine,lerr)
    else
       mass%mass = mass%m0
    end if

    if (associated(mass%IIEngine)) then
       mass%II = mass%II0 + mass%II1*EngineValue(mass%IIEngine,lerr)
    else
       mass%II = mass%II0
    end if

    if (lerr < 0) then
       if (present(ierr)) then
          ierr = ierr + lerr
          call reportError (debugFileOnly_p,'updateMass')
       end if
       return
    else if (restart) then
       return
    end if

    if (.not. mass%addedMass) then ! No gravity force for added mass
       call calcGravitationForces (mass, gravity, eV)
       call updateNodeForce (mass%triad, -eV)
    end if

    call calcInertiaForces (mass, eV)
    call updateNodeForce (mass%triad, eV)

  end subroutine updateMass


  !!============================================================================
  !> @brief Calculates the gravitation forces for a mass element.
  !>
  !> @param Q External force vector (gravity loads)
  !> @param RF System reaction forces associated with constrained DOFs
  !> @param[in] mass The mass element to calculate gravitation forces for
  !> @param[in] gravity Gravitation vector
  !> @param[in] sam Data for managing system matrix assembly
  !> @param ierr Error flag
  !>
  !> @details The calculated forces are added into the system force vector @a Q
  !> with contributions to the reaction forces @a RF for the constrained DOFs.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Karl Erik Thoresen
  !> @date Feb 1999
  !>
  !> @author Bjorn Haugen
  !> @date Nov 2001
  !>
  !> @author Knut Morten Okstad
  !> @date 30 Jan 2008

  subroutine addInGravitationForces (Q,RF,mass,gravity,sam,ierr)

    use SamModule         , only : SamType, dp
    use MassTypeModule    , only : MassType
    use AsmExtensionModule, only : csAddEV
    use reportErrorModule , only : reportError, debugFileOnly_p

    real(dp)      , intent(inout) :: Q(:), RF(:)
    type(MassType), intent(in)    :: mass
    real(dp)      , intent(in)    :: gravity(3)
    type(SamType) , intent(in)    :: sam
    integer       , intent(inout) :: ierr

    !! Local variables
    integer  :: lerr, ndofs
    real(dp) :: eV(6)

    !! --- Logic section ---

    ndofs = mass%triad%nDOFs
    if (ndofs < 1 .or. mass%addedMass) return ! No gravity force for added mass

    call calcGravitationForces (mass, gravity, eV)
    call csAddEV (sam, -mass%samElnum, ndofs, eV(1:ndofs), Q, RF, lerr)
    if (lerr == 0) return

    ierr = ierr - 1
    call reportError (debugFileOnly_p,'addInGravitationForces')

  end subroutine addInGravitationForces


  !!============================================================================
  !> @brief Calculates the inertia forces for a mass element.
  !>
  !> @param FI Inertial force vector
  !> @param RF System reaction forces associated with constrained DOFs
  !> @param[in] mass The mass element to calculate inertia forces for
  !> @param[in] sam Data for managing system matrix assembly
  !> @param ierr Error flag
  !>
  !> @details The calculated forces are added into the system force vector @a Q
  !> with contributions to the reaction forces @a RF for the constrained DOFs.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Karl Erik Thoresen
  !> @date Feb 1999
  !>
  !> @author Bjorn Haugen
  !> @date Nov 2001
  !>
  !> @author Knut Morten Okstad
  !> @date 30 Jan 2008

  subroutine addInInertiaForces (FI,RF,mass,sam,ierr)

    use SamModule         , only : SamType, dp
    use MassTypeModule    , only : MassType
    use AsmExtensionModule, only : csAddEV
    use reportErrorModule , only : reportError, debugFileOnly_p

    real(dp)      , intent(inout) :: FI(:), RF(:)
    type(MassType), intent(in)    :: mass
    type(SamType) , intent(in)    :: sam
    integer       , intent(inout) :: ierr

    !! Local variables
    integer  :: lerr, ndofs
    real(dp) :: eV(6)

    !! --- Logic section ---

    ndofs = mass%triad%nDOFs
    if (ndofs < 1) return

    call calcInertiaForces (mass, eV)
    call csAddEV (sam, mass%samElnum, ndofs, eV(1:ndofs), FI, RF, lerr)
    if (lerr == 0) return

    ierr = ierr - 1
    call reportError (debugFileOnly_p,'addInInertiaForces')

  end subroutine addInInertiaForces


  !!============================================================================
  !> @brief Calculates the gravitation forces for a mass element.
  !>
  !> @param[in] mass The mass element to calculate gravitation forces for
  !> @param[in] gravity Gravitation vector
  !> @param[out] eV Element force vector
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 30 Jan 2008

  subroutine calcGravitationForces (mass,gravity,eV)

    use MassTypeModule , only : MassType, dp
    use TriadTypeModule, only : transVGlobToTriad, transTriadToSys

    type(MassType), intent(in)  :: mass
    real(dp)      , intent(in)  :: gravity(3)
    real(dp)      , intent(out) :: eV(6)

    !! --- Logic section ---

    eV(1:3) = mass%mass * transVGlobToTriad(mass%triad,gravity)
    eV(4:6) = 0.0_dp

    call transTriadToSys (mass%triad,eV(1:3))

  end subroutine calcGravitationForces


  !!============================================================================
  !> @brief Calculates the inertia forces for a mass element.
  !>
  !> @param[in] mass The mass element to calculate inertia forces for
  !> @param[out] eV Element force vector
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 30 Jan 2008

  subroutine calcInertiaForces (mass,eV)

    use MassTypeModule   , only : MassType, dp
    use TriadTypeModule  , only : transVGlobToTriad, transVTriadToGlob
    use TriadTypeModule  , only : transTriadToGlob, transGlobToSys
    use manipMatrixModule, only : cross_product

    type(MassType), intent(in)  :: mass
    real(dp)      , intent(out) :: eV(6)

    !! Local variables
    real(dp) :: vel(3), acc(3), II(3,3)

    !! --- Logic section ---

    !! Translational masses: F = m*a
    acc     = transVGlobToTriad(mass%triad,mass%triad%urdd(1:3))
    eV(1:3) = transVTriadToGlob(mass%triad,mass%mass*acc)

    !! Rotational masses: M = II*dw + wx(II*w)  (w == omega, dw == domega/dt)
    if (mass%triad%nDOFs == 6) then
       II  = mass%II
       vel = mass%triad%urd(4:6)
       acc = mass%triad%urdd(4:6)
       call transTriadToGlob (mass%triad,II,1,3)
       eV(4:6) = matmul(II,acc) + cross_product(vel,matmul(II,vel))
    else
       eV(4:6) = 0.0_dp
    end if

    call transGlobToSys (mass%triad,eV(1:mass%triad%nDOFs))

  end subroutine calcInertiaForces


  !!============================================================================
  !> @brief Calculates the mass matrix for a mass element.
  !>
  !> @param[in] scaleI Mass matrix scaling factor
  !> @param Nmat System Newton matrix
  !> @param[in] mass The mass element to calculate mass matrix for
  !> @param[in] sam Data for managing system matrix assembly
  !> @param ierr Error flag
  !> @param Rhs System right-hand-side vector associated with the Newton matrix
  !>
  !> @details This subroutine calculates the mass matrix for a mass element
  !> and adds it into the system Newton matrix,
  !> multiplied by the given scale factor @a scaleI.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Karl Erik Thoresen
  !>
  !> @date Feb 1999

  subroutine addInMassMat (scaleI,Nmat,mass,sam,ierr,Rhs)

    use SamModule          , only : SamType, dp
    use SysMatrixTypeModule, only : SysMatrixType
    use MassTypeModule     , only : MassType
    use TriadTypeModule    , only : transTriadToSys
    use AsmExtensionModule , only : csAddElm
    use reportErrorModule  , only : reportError, debugFileOnly_p

    real(dp)           , intent(in)    :: scaleI
    type(SysMatrixType), intent(inout) :: Nmat
    type(MassType)     , intent(in)    :: mass
    type(SamType)      , intent(in)    :: sam
    integer            , intent(out)   :: ierr
    real(dp), optional , intent(inout) :: Rhs(:)

    !! Local variables
    integer  :: ndofs
    real(dp) :: elMat(6,6)

    !! --- Logic section ---

    ierr  = 0
    ndofs = mass%triad%nDOFs
    if (ndofs < 1) return

    elMat = 0.0_dp
    elMat(1,1) = mass%mass(1) * scaleI
    elMat(2,2) = mass%mass(2) * scaleI
    elMat(3,3) = mass%mass(3) * scaleI
    if (ndofs > 3) elMat(4:6,4:6) = mass%II * scaleI

    call transTriadToSys (mass%triad,elMat,1,ndofs)
    call csAddElM (sam,mass%samElNum,ndofs,elMat(1:ndofs,1:ndofs),Nmat,ierr,Rhs)
    if (ierr < 0) call reportError (debugFileOnly_p,'addInMassMat')

  end subroutine addInMassMat

end module MassRoutinesModule
