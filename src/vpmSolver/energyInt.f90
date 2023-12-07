!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file energyInt.f90
!>
!> @brief Subroutines for structural energy integration.

!!==============================================================================
!> @brief Module with subroutines for structural energy integration.
!>
!> @details This module contains a set of subroutines for calculating the
!> structural energy contribution for each object of the meachanism.
!>
!> @author Bjorn Haugen
!>
!> @date 1 Oct 1997

module EnergyIntegrationModule

  implicit none

  private

  public :: energyIntegration, energyInitialization


contains

  !!============================================================================
  !> @brief Prints out a warning if negative energy loss is detected.
  !>
  !> @param[in] objType Type name of the object to print message for
  !> @param[in] objId Id of the object to print message for
  !> @param[in] time Current simulation time
  !> @param[in] Eloss Energy loss value
  !> @param[in] F Force value causing the energy loss
  !> @param[in] du Displacement increment causing the energy loss
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 11 Jul 2014

  subroutine checkLoss (objType,objId,time,Eloss,F,du)

    use IdTypeModule     , only : IdType, getId
    use kindModule       , only : dp
    use reportErrorModule, only : getErrorFile, reportError
    use reportErrorModule, only : note_p, warning_p, warningFileOnly_p

    character(len=*), intent(in) :: objType
    type(IdType)    , intent(in) :: objId
    real(dp)        , intent(in) :: time, Eloss, F, du

    !! Local variables
    integer, parameter :: maxWarnings_p = 30
    integer, save      :: nWarnings = 0
    integer            :: msgUnit
    character(len=24)  :: ctime

    !! --- Logic section ---

    if (F*du < -1.0e-6_dp*max(1.0e-6_dp,Eloss)) then
       write(ctime,"(' at T =',1PE12.5)") time
       nWarnings = nWarnings + 1
       if (nWarnings > maxWarnings_p) then
          msgUnit = warningFileOnly_p
       else
          msgUnit = warning_p
       end if
       call reportError (msgUnit,'Negative energy loss detected in '// &
            &            objType//trim(getId(objId))//ctime)
       write(getErrorFile(),600) F,du,F*du
600    format(10X,'F =',1PE13.5,' du =',E13.5,' dE =',E13.5)
       if (nWarnings == maxWarnings_p) call reportError (note_p, &
           'No further "Negative energy loss" messages will be given.', &
           'See fedem_solver.res file for the complete history.')
    end if

  end subroutine checkLoss


  !!============================================================================
  !> @brief Integrates the system energies for current time step.
  !>
  !> @param mech Mechanism components of the model
  !> @param[in] time Current simulation time
  !> @param[in] dTime Time increment size
  !> @param[in] linearStatic If .true., we are doing a linear static analysis
  !> @param[in] doStiffDampE If .true., use mean deformational displacement
  !> velocity when calculating energy loss due to stiffness proportional damping
  !> @param[out] ierr Error flag
  !>
  !> @details If @a dTime is zero, only the potential, kinetical and strain
  !> energies are computed. The remaining quantities are then identically zero.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 2 Jun 2004

  subroutine energyIntegration (mech,time,dTime,linearStatic,doStiffDampE,ierr)

    use MechanismTypeModule, only : MechanismType, dp
    use reportErrorModule  , only : reportError, error_p

    type(MechanismType), intent(inout) :: mech
    real(dp)           , intent(in)    :: time, dTime
    logical            , intent(in)    :: linearStatic, doStiffDampE
    integer            , intent(out)   :: ierr

    !! Local variables
    integer             :: i, j, nDofs(2)
    real(dp), parameter :: rigidFactor_p = 0.1_dp
    real(dp)            :: Edens(2)
    character(len=128)  :: errMsg

    !! --- Logic section ---

    ierr      = 0
    mech%Epot = 0.0_dp
    mech%Ekin = 0.0_dp
    mech%Estr = 0.0_dp
    mech%Edmp = 0.0_dp
    mech%Einp = 0.0_dp
    mech%Eext = 0.0_dp
    Edens     = 0.0_dp
    nDofs     = 0

    !! === ENERGY CONTRIBUTION FROM SUPERELEMENTS

    do i = 1, size(mech%sups)
       call SupElEnergies (mech%sups(i), &
            &              mech%gravity,mech%env%seaLevel,dTime,doStiffDampE, &
            &              mech%Epot,mech%Ekin,mech%Estr,mech%Edmp,ierr)
       if (mech%sups(i)%rigidFlag == 2) then
          Edens(1) = Edens(1) + mech%sups(i)%Estr
          nDofs(1) = nDofs(1) + mech%sups(i)%nTotDofs
       else
          Edens(2) = Edens(2) + mech%sups(i)%Estr
          nDofs(2) = nDofs(2) + mech%sups(i)%nTotDofs
       end if
    end do

    !! === ENERGY CONTRIBUTION FROM SUPERELEMENT LOADS

    do i = 1, size(mech%supLoads)
       call SupElLoadEnergies (mech%supLoads(i),linearStatic,dTime,mech%Einp)
    end do

    !! === ENERGY CONTRIBUTION FROM EXTERNAL FORCES

    do i = 1, size(mech%forces)
       call ForceEnergies (mech%forces(i),linearStatic,dTime,mech%Einp)
    end do

    !! === ENERGY CONTRIBUTION FROM PRESCRIBED MOTIONS

    do i = 1, size(mech%motions)
       call MotionEnergies (mech%motions(i),linearStatic,dTime,mech%Einp)
    end do

    !! === POTENTIAL AND KINETIC ENERGY FROM LUMPED MASSES

    do i = 1, size(mech%masses)
       call MassEnergies (mech%masses(i),mech%gravity,dTime,mech%Epot,mech%Ekin)
    end do

    !! === STRAIN AND INPUT ENERGY FROM SPRINGS

    do i = 1, size(mech%baseSprings)
       call SpringEnergies (mech%baseSprings(i),linearStatic,time,dTime, &
            &               mech%Estr,mech%Einp,mech%Edmp,ierr)
       if (mech%baseSprings(i)%isActive .and. mech%baseSprings(i)%dof >= 0) then
          Edens(2) = Edens(2) + mech%baseSprings(i)%Estr
          nDofs(2) = nDofs(2) + 1
       end if
    end do

    !! === ENERGY LOSS FROM DISCRETE DAMPERS

    do i = 1, size(mech%baseDampers)
       if (mech%baseDampers(i)%isActive) then
          call DamperEnergies (mech%baseDampers(i),time,dTime,mech%Edmp)
       end if
    end do

    !! === ENERGY LOSS FROM FRICTION

    do i = 1, size(mech%frictions)
       call FrictionEnergies (mech%frictions(i)%p,time,dTime,mech%Edmp)
    end do

    !! === ENERGY CONTRIBUTION FROM TIRES

    do i = 1, size(mech%tires)
       call TireEnergies (mech%tires(i),mech%gravity,dTime, &
            &             mech%Epot,mech%Ekin,mech%Eext)
    end do

    !! === ENERGY CONTRIBUTION FROM AERODYNAMIC FORCES

    if (dTime > 0.0_dp .and. associated(mech%turbine)) then
       do i = 1, size(mech%turbine%node,1)
          do j = 1, size(mech%turbine%node,2)
             call AeroEnergy (mech%turbine%node(i,j),mech%Eext)
          end do
       end do
    end if

    !! Check that the strain energy from the rigid parts (if any) is small
    if (nDofs(1) > 0) then
       if (nDofs(2) > 0) then ! compare against the flexible parts
          if (Edens(1)/nDofs(1) > rigidFactor_p * Edens(2)/nDofs(2)) then
             call rigidWarning (Edens(1)/nDofs(1),Edens(2)/nDofs(2),nDofs(2))
          end if
       else ! no flexible parts, check against the total pot+kin energy instead
          if (Edens(1) > rigidFactor_p * abs(mech%Epot+mech%Ekin)) then
             call rigidWarning (Edens(1),abs(mech%Epot+mech%Ekin),nDofs(2))
          end if
       end if
    end if

    if (ierr < 0) then
       write(errMsg,"('Energy integration failed at t =',1PE12.5)") time
       call reportError (error_p,errMsg,addString='energyIntegration')
    end if

  contains

    !> @brief Prints a warning on too much strain energy for semi-rigid parts.
    subroutine rigidWarning (rigEn,flexEn,nDofF)
      use reportErrorModule, only : warning_p, warningFileOnly_p
      use reportErrorModule, only : empty_p, noteFileOnly_p
      real(dp), intent(in) :: rigEn, flexEn
      integer , intent(in) :: nDofF
      integer , save       :: nWarning = 0, msgUnit = warning_p
      nWarning = nWarning + 1
      if (nWarning > 100) return
      if (nWarning == 2) msgUnit = warningFileOnly_p
      if (nDofF > 0) then
         write(errMsg,600) rigEn,flexEn
      else
         write(errMsg,610) rigEn,flexEn
      end if
      call reportError (msgUnit, &
           'The rigid parts of the mechanism have significant strain energy.', &
           addString=errMsg)
600   format('Strain energy density in rigid parts:',1PE12.5,2X, &
           & 'In flexible parts:',E12.5)
610   format('Strain energy in rigid parts:',1PE12.5,2X, &
           & 'Total potential+kinetic energy:', E12.5)
      if (nWarning == 1) then
         call reportError (empty_p, &
              'This may indicate a topology error in the mechanism assembly.', &
              'Please verify that the mechanism is correctly modelled,', &
              'or manually increase the pseudo-stiffness of the rigid part(s).')
      else if (nWarning == 100) then
         call reportError (noteFileOnly_p, &
              '100 insidences of the above warning have been encountered.', &
              'Any further occurances of this warning will be ommitted.')
      end if
    end subroutine rigidWarning

  end subroutine energyIntegration


  !!============================================================================
  !> @brief Initializes the energy integration variables.
  !>
  !> @param mech Mechanism components of the model
  !> @param[in] time Current simulation time
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 2 Jun 2004

  subroutine energyInitialization (mech,time,ierr)

    use MechanismTypeModule, only : MechanismType, dp

    type(MechanismType), intent(inout) :: mech
    real(dp)           , intent(in)    :: time
    integer            , intent(out)   :: ierr

    !! Local variables
    integer :: i

    !! --- Logic section ---

    !! There is no need to update the mechanism objects here, as that already
    !! has been done as the final step of the static equilibrium iterations.
    !! Since the Engine update function now works recursively, depth-first
    !! (it has been like that since the Newmark cleanup prior to R3.0),
    !! each mechanism object needs only a single update call to be consistent
    !! with the current configuration (kmo).

    !! Compute the initial energies with zero time step, i.e.,
    !! only potential, kinetical and strain energies are computed here.
    call energyIntegration (mech,time,0.0_dp,.false.,.false.,ierr)
    if (ierr < 0) return

    !! Store the initial potential energies for links, masses and tires.
    !! The potential energies computed later are then relative to these.
    do i = 1, size(mech%sups)
       mech%sups(i)%ePot0 = mech%sups(i)%ePot
       mech%sups(i)%ePot  = 0.0_dp
    end do

    do i = 1, size(mech%masses)
       mech%masses(i)%ePot0 = mech%masses(i)%ePot
       mech%masses(i)%ePot  = 0.0_dp
    end do

    do i = 1, size(mech%tires)
       mech%tires(i)%ePot0 = mech%tires(i)%ePot
       mech%tires(i)%ePot  = 0.0_dp
    end do

    !! Total energy of the initial configuration
    mech%Etot = mech%Ekin + mech%Estr
    mech%Epot = 0.0_dp

  end subroutine energyInitialization


  !!============================================================================
  !> @brief Calculates the system energy contributions from a superelement.
  !>
  !> @param supel Superelement to calculate energies for
  !> @param[in] gravVec Gravitation vector
  !> @param[in] seaLevel Current sea level
  !> @param[in] dTime Time increment size
  !> @param[in] doStiffDampE If .true., use mean deformational displacement
  !> velocity when calculating energy loss due to stiffness proportional damping
  !> @param Epot Accumulated potential energy for the mechanism
  !> @param Ekin Accumulated kinetic energy for the mechanism
  !> @param Estr Accumulated elastic strain energy for the mechanism
  !> @param Edmp Energy loss due to damping for the mechanism
  !> @param ierr Error flag
  !>
  !> @details The following contributions are included:
  !>  - Potential energy
  !>  - Kinetic energy
  !>  - Elastic strain energy
  !>  - Damping loss
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 13 Oct 1999

  subroutine SupElEnergies (supel,gravVec,seaLevel,dTime,doStiffDampE, &
       &                    Epot,Ekin,Estr,Edmp,ierr)

    use SupElTypeModule           , only : SupElType, dp
    use TriadTypeModule           , only : GetTriadIncrement
    use MassMatrixCorrectionModule, only : mmcMassMatrixWarning
    use scratchArrayModule        , only : getRealScratchArray
    use manipMatrixModule         , only : matmul34
    use rotationModule            , only : deltaRot, vec_to_mat
    use reportErrorModule         , only : reportError, debugFileOnly_p

    type(SupElType), intent(inout) :: supel
    real(dp)       , intent(in)    :: gravVec(3), seaLevel, dTime
    logical        , intent(in)    :: doStiffDampE
    real(dp)       , intent(inout) :: Epot, Ekin, Estr, Edmp
    integer        , intent(inout) :: ierr

    !! Local variables
    integer           :: i, j, n
    logical           :: haveStructDamping(2)
    real(dp), pointer :: work(:), v1(:), v2(:), v3(:)
    real(dp)          :: dTimeInv, vec(3), dRot(3,3), supTmid(3,3)

    !! --- Logic section ---

    if (associated(supel%genDOFs)) then
       j = supel%genDOFs%nDOFs
    else
       j = 0
    end if
    work => getRealScratchArray(supel%nTotDofs*2+j,i)
    if (i < 0) goto 915

    v1 => work(1:supel%nTotDofs)
    v2 => work(1+supel%nTotDofs:supel%nTotDofs*2)
    v3 => work(1+supel%nTotDofs*2:)

    if (dTime > 0.0_dp) then

       !! Build mean velocity vector between previous and current time step
       dTimeInv = 1.0_dp/dTime

       do i = 1, supel%nExtNods
          j = supel%triads(i)%firstDOF
          n = supel%triads(i)%p%nDOFs
          v1(j:j+n-1) = GetTriadIncrement(supel%triads(i)%p)*dTimeInv
       end do

       if (associated(supel%genDOFs)) then
          j = supel%genDOFs%firstDOF
          n = supel%genDOFs%nDOFs
          v1(j:j+n-1) = (supel%genDOFs%ur-supel%genDOFs%urPrev)*dTimeInv
       end if

       !! Compute mid-point transformation matrix for the element
       vec = deltaRot(supel%supTrPrev(:,1:3),supel%supTr(:,1:3))
       call vec_to_mat (0.5_dp*vec,dRot)
       supTmid = matmul(dRot,supel%supTrPrev(:,1:3))

       !! Transform total velocities for the superelement dofs
       !! from global to local coordinate system of the superelement
       do i = 1, supel%nExtNods
          j = supel%triads(i)%firstDOF
          v1(j:j+2) = matmul(v1(j:j+2),supTmid)
          if (supel%triads(i)%p%nDOFs == 6) then
             v1(j+3:j+5) = matmul(v1(j+3:j+5),supTmid)
          end if
       end do

    else
       dTimeInv = 0.0_dp
    end if

    !! Compute strain energy, using current deformational
    !! displacements stored in supel%finit
    supel%eStr = 0.5_dp*dot_product(matmul(supel%KmMat,supel%finit),supel%finit)

    if (associated(supel%Mmat)) then
       !! Calculate kinetic energy at current time
       if (associated(supel%Mamat)) then
          v2 = matmul(supel%Mmat+supel%Mamat,supel%uld)
       else
          v2 = matmul(supel%Mmat,supel%uld)
       end if
       supel%eKin = 0.5_dp*dot_product(v2,supel%uld)
    end if

    !! Calculate damping energy loss from previous time step to current
    haveStructDamping = .false.
    if (dTime > 0.0_dp) then

       !! Check if we actually have any structural damping
       if (supel%mDmpFactor > 0.0_dp) then
          haveStructDamping(1) = .true.
       else if (associated(supel%genDOFs)) then
          haveStructDamping(1) = any(supel%genDOFs%alpha1 > 0.0_dp)
       end if
       if (supel%kDmpFactor > 0.0_dp) then
          haveStructDamping(2) = .true.
       else if (associated(supel%genDOFs)) then
          haveStructDamping(2) = any(supel%genDOFs%alpha2 > 0.0_dp)
       end if

       if (associated(supel%Cdmat)) then
          !! Energy loss due to fluid drag forces is computed using the mean
          !! of the velocity from previous to current time step (stored in v1)
          !! Note: This is maybe not 100% correct since the damping matrix Cdmat
          !! also depends on the relative structure velocity, and its current
          !! value is computed at current time step. Ideally, we should use the
          !! mean value from this and previous increment for that one too.
          supel%eDmp = supel%eDmp + dTime*dot_product(matmul(supel%Cdmat,v1),v1)
       end if

       if (associated(supel%genDOFs)) then
          n = supel%genDOFs%firstDOF-1
       else
          n = supel%nTotDofs
       end if

       if (haveStructDamping(1) .and. associated(supel%Mmat)) then
          !! Energy loss due to mass proportional damping is computed using the
          !! mean of the velocity from previous to current time step (v1)
          v1(1:n) = v1(1:n) * sqrt(supel%mDmpFactor)
          if (associated(supel%genDOFs)) then
             v1(n+1:) = v1(n+1:) * sqrt(supel%genDOFs%alpha1)
          end if
          supel%eDmp = supel%eDmp + dTime*dot_product(matmul(supel%Mmat,v1),v1)
       end if

       if (haveStructDamping(2) .and. doStiffDampE) then
          !! Energy loss due to stiffness proportional damping is computed
          !! through the mean deformational displacement velocity between
          !! previous and current time step
          v2 = (supel%finit - supel%finitPrev)*dTimeInv
          v2(1:n) = v2(1:n) * sqrt(supel%kDmpFactor)
          if (associated(supel%genDOFs)) then
             v2(n+1:) = v2(n+1:) * sqrt(supel%genDOFs%alpha2)
          end if
          supel%eDmp = supel%eDmp + dTime*dot_product(matmul(supel%KmMat,v2),v2)
       end if

    end if

    if (associated(supel%genDOFs)) then
       !! Calculate energy contributions from each component mode
       j  = supel%genDOFs%firstDOF
       i  = j + supel%genDOFs%nDOFs-1
       if (associated(supel%Mmat)) then
          v3 = matmul(supel%Mmat(j:i,:),supel%uld)
          supel%genDOFs%energy(1,:) = v3 * supel%uld(j:i)
       end if
       v3 = matmul(supel%KmMat(j:i,:),supel%finit)
       supel%genDOFs%energy(2,:) = v3 * supel%finit(j:i)
       if (haveStructDamping(1) .and. associated(supel%Mmat)) then
          v3 = dTime*matmul(supel%Mmat(j:i,:),v1)
          supel%genDOFs%energy(3,:) = supel%genDOFs%energy(3,:) + v3 * v1(j:i)
       end if
       if (haveStructDamping(2) .and. doStiffDampE) then
          v3 = dTime*matmul(supel%KmMat(j:i,:),v2)
          supel%genDOFs%energy(3,:) = supel%genDOFs%energy(3,:) + v3 * v2(j:i)
       end if
    end if

    !! Find the global mass-center position
    vec = matmul34(supel%supTr,supel%posMassCenter)

    !! Potential energy as product of gravitational vector and
    !! mass center position vector
    supel%ePot = -dot_product(gravVec,vec)*supel%mass
    if (dTime > 0.0_dp) supel%ePot = supel%ePot - supel%ePot0

    !! Add potential energy due to buoyancy, if any
    if (associated(supel%hydyn)) then
       !! Global sea level position
       dRot(:,1) = matmul(supel%supTr(:,1:3),supel%hydyn%wn*seaLevel)
       !! Global buoyancy center
       dRot(:,2) = matmul34(supel%supTr,supel%hydyn%C0b(:,1))
       !! Global buoyancy force
       vec = matmul(supel%supTr(:,1:3),supel%hydyn%B(1:3))
       !! Potential energy as product of buoyancy force and relative position
       supel%ePot = supel%ePot + dot_product(dRot(:,1)-dRot(:,2),vec)
    end if

    !! Update mechanism energies
    Epot = Epot + supel%Epot
    Ekin = Ekin + supel%Ekin
    Estr = Estr + supel%Estr
    Edmp = Edmp + supel%Edmp

    !! The previous superelement state values are updated for the next time step
    !! in subroutine updateAtConvergence invoked from terminateStep

    if (.not. supel%saveVar(3)) return
    if (.not. associated(supel%Mmat)) return

    !! Print warning for large moment error for high speed rotation.
    !! Note: This will destroy the contents of v1 and v2 since the same
    !!       scratch array is used locally within mmcMassMatrixWarning.
    call mmcMassMatrixWarning (supel,supel%uld,i)
    if (i == 0) return

915 ierr = ierr - 1
    call reportError (debugFileOnly_p,'SupElEnergies')

  end subroutine SupElEnergies


  !!============================================================================
  !> @brief Calculates the incremental work done by a superelement load.
  !>
  !> @param supLoad Superelement load to calculate energy for
  !> @param[in] linearStatic If .true., we are doing a linear static analysis
  !> @param[in] dTime Time increment size
  !> @param Einp Accumulated energy done by external forces
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 17 Apr 2008

  subroutine SupElLoadEnergies (supLoad,linearStatic,dTime,Einp)

    use SupElLoadTypeModule, only : SupElLoadType
    use TriadTypeModule    , only : GetTriadIncrement, dp

    type(SupElLoadType), intent(inout) :: supLoad
    logical            , intent(in)    :: linearStatic
    real(dp)           , intent(in)    :: dTime
    real(dp)           , intent(inout) :: Einp

    !! Local variables
    integer  :: i, j, n, lc
    real(dp) :: u(6), dEinp

    !! --- Logic section ---

    dEinp = 0.0_dp
    if (linearStatic .or. dTime > 0.0_dp) then

       !! Compute the incremental input energy with unit load amplitude

       u = 0.0_dp
       lc = supLoad%loadCase
       do i = 1, supLoad%sup%nExtNods
          j = supLoad%sup%triads(i)%firstDOF
          n = supLoad%sup%triads(i)%p%nDOFs
          if (n > 0) then
             u(1:n) = GetTriadIncrement(supLoad%sup%triads(i)%p)
             dEinp = dEinp + dot_product(supLoad%sup%S(j:j+n-1,lc),u(1:n))
          end if
       end do

       if (associated(supLoad%sup%genDOFs)) then
          j = supLoad%sup%genDOFs%firstDOF
          do i = 1, supLoad%sup%genDOFs%nDOFs
             u(1) = supLoad%sup%genDOFs%ur(i) - supLoad%sup%genDOFs%urPrev(i)
             dEinp = dEinp + supLoad%sup%S(j+i-1,lc)*u(1)
          end do
       end if

    end if
    if (linearStatic) then
       !! Calculate the input energy using current amplitude
       supLoad%eInp = supLoad%F * dEinp
    else if (dTime > 0.0_dp) then
       !! Accumulate incremental input energy, using the mean amplitude
       supLoad%eInp = supLoad%eInp + 0.5_dp*(supLoad%F+supLoad%Fprev) * dEinp
    else
       supLoad%eInp = 0.0_dp
    end if

    !! Update mechanism energies
    Einp = Einp + supLoad%eInp

    !! Update the previous load value for the next time step
    supLoad%Fprev = supLoad%F

  end subroutine SupElLoadEnergies


  !!============================================================================
  !> @brief Calculates the incremental work done by an external force.
  !>
  !> @param force External force to calculate energy for
  !> @param[in] linearStatic If .true., we are doing a linear static analysis
  !> @param[in] dTime Time increment size
  !> @param Einp Accumulated energy done by external forces
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 14 Oct 1997

  subroutine ForceEnergies (force,linearStatic,dTime,Einp)

    use ForceTypeModule, only : ForceType, forceVec_p, momentVec_p, dp
    use rotationModule , only : deltaRot

    type(ForceType), intent(inout) :: force
    logical        , intent(in)    :: linearStatic
    real(dp)       , intent(in)    :: dTime
    real(dp)       , intent(inout) :: Einp

    !! Local variables
    integer  :: ldof
    real(dp) :: fMean, fVecMean(3), uVec(3), dEinp

    !! --- Logic section ---

    if (force%loadType > 0) return ! Skip frequency-domain loads

    dEinp = 0.0_dp
    if (linearStatic .or. dTime > 0.0_dp) then

       !! Compute displacement vector from previous to current

       if (associated(force%triad)) then

          select case (force%dof)

          case (forceVec_p, 1:3)  ! vector force or X, Y or Z force
             uVec = force%triad%ur(:,4) -  force%triad%urPrev(:,4)
             ldof = force%dof

          case (momentVec_p, 4:6) ! vector moment or X, Y or Z moment
             uVec = deltaRot(force%triad%urPrev(:,1:3),force%triad%ur(:,1:3))
             ldof = force%dof - 3

          case default
             return

          end select

       else if (associated(force%joint)) then

          uVec(1) = force%joint%jointDofs(force%dof)%jVar(1) &
               &  - force%joint%jointDofs(force%dof)%jVarPrev
          ldof = 1

       else
          return
       end if

       !! Compute incremental input energy

       select case (force%dof)

       case (forceVec_p, momentVec_p)
          if (linearStatic) then
             dEinp = force%F*dot_product(force%forceDir,uVec)
          else
             fVecMean = force%F*force%forceDir + force%Fprev*force%forceDirPrev
             dEinp    = 0.5_dp*dot_product(fVecMean,uVec)
          end if

       case (1:6) ! works for joint forces as well
          if (linearStatic) then
             dEinp = force%F*uVec(ldof)
          else
             fMean = force%F + force%Fprev
             dEinp = 0.5_dp*fMean*uVec(ldof)
          end if

       end select

    end if
    if (linearStatic) then
       force%eInp = dEinp
    else if (dTime > 0.0_dp) then
       force%eInp = force%eInp + dEinp
    else
       force%eInp = 0.0_dp
    end if

    !! Update mechanism energies
    Einp = Einp + force%eInp

    !! Update the previous force values for the next time step
    force%Fprev        = force%F
    force%forceDirPrev = force%forceDir

  end subroutine ForceEnergies


  !!============================================================================
  !> @brief Calculates the incremental work done by a prescribed motion.
  !>
  !> @param motion Prescribed motion to calculate energy for
  !> @param[in] linearStatic If .true., we are doing a linear static analysis
  !> @param[in] dTime Time increment size
  !> @param Einp Accumulated energy done by external forces and motions
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 4 Oct 2005

  subroutine MotionEnergies (motion,linearStatic,dTime,Einp)

    use MotionTypeModule, only : MotionType, dp

    type(MotionType), intent(inout) :: motion
    logical         , intent(in)    :: linearStatic
    real(dp)        , intent(in)    :: dTime
    real(dp)        , intent(inout) :: Einp

    !! --- Logic section ---

    if (linearStatic) then
       !! Calculate the input energy
       motion%eInp = motion%F * motion%cc
    else if (dTime > 0.0_dp) then
       !! Accumulate incremental input energy
       motion%eInp = motion%eInp + 0.5_dp*(motion%F + motion%Fprev) * motion%cc
    else
       motion%eInp = 0.0_dp
    end if

    !! Update mechanism energies
    Einp = Einp + motion%eInp

    !! Update the previous force value for the next time step
    motion%Fprev = motion%F

  end subroutine MotionEnergies


  !!============================================================================
  !> @brief Calculates energy contributions from a tire.
  !>
  !> @param tire Tire to calculate energies for
  !> @param[in] gravVec Gravitation vector
  !> @param[in] dTime Time increment size
  !> @param Epot Accumulated potential energy for the mechanism
  !> @param Ekin Accumulated kinetic energy for the mechanism
  !> @param Einp Accumulated input energy from external forces and motions
  !>
  !> @details The kinetic and potential energy from the tire mass is calculated,
  !> as well as the incremental work done by the tire force from previous to
  !> current time step.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 31 Oct 2002

  subroutine TireEnergies (tire,gravVec,dTime,Epot,Ekin,Einp)

    use TireTypeModule , only : TireType, dp
    use TriadTypeModule, only : transVGlobToTriad, GetTriadIncrement

    type(TireType), intent(inout) :: tire
    real(dp)      , intent(in)    :: gravVec(3), dTime
    real(dp)      , intent(inout) :: Epot, Ekin, Einp

    !! Local variables
    real(dp) :: veloc(6), II(3), fMean(6), dU(6)

    !! --- Logic section ---

    veloc = transVGlobToTriad(tire%RimTriad,tire%RimTriad%urd)

    !! Kinetic energy from translations
    tire%Ekin = 0.5_dp*tire%mass*dot_product(veloc(1:3),veloc(1:3))

    !! Kinetic energy from rotations
    II(1) = tire%Iyy
    II(2) = tire%Iyy
    II(3) = tire%Izz
    tire%Ekin = tire%Ekin + 0.5_dp*dot_product(veloc(4:6)*II,veloc(4:6))

    !! Potential energy from tire masses
    dU(1:3) = tire%RimTriad%ur(:,4) + tire%ZoffSet*tire%RimTriad%ur(:,3)
    tire%Epot = -tire%mass*dot_product(gravVec,dU(1:3))

    if (dTime > 0.0_dp) then

       tire%Epot = tire%Epot - tire%Epot0

       !! Mean force and moment between previous and current time step
       fMean = 0.5_dp*(tire%forceInG + tire%forceInGPrev)

       !! Compute displacement vector from previous to current configuration
       dU = GetTriadIncrement(tire%RimTriad)

       !! Incremental energy from force and moment
       tire%eInp = tire%eInp + dot_product(fMean,dU)

    end if

    !! Update mechanism energies
    Epot = Epot + tire%ePot
    Ekin = Ekin + tire%eKin
    Einp = Einp + tire%eInp

    !! Update the previous tire values for the next time step
    tire%forceInGPrev = tire%forceInG

  end subroutine TireEnergies


  !!============================================================================
  !> @brief Computes strain energy, input energy and energy loss for a spring.
  !>
  !> @param spr Spring to calculate energies for
  !> @param[in] linearStatic If .true., we are doing a linear static analysis
  !> @param[in] time Current simulation time
  !> @param[in] dTime Time increment size
  !> @param Estr Accumulated elastic strain energy for the mechanism
  !> @param Einp Accumulated input energy from external forces and motions
  !> @param Edmp Energy loss due to damping for the mechanism
  !> @param ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 1 Oct 1997

  subroutine SpringEnergies (spr,linearStatic,time,dTime,Estr,Einp,Edmp,ierr)

    use SpringTypeModule    , only : SpringBaseType, springEnergy, dp
    use EngineRoutinesModule, only : EngineValue
    use IdTypeModule        , only : getId
    use reportErrorModule   , only : reportError, error_p

    type(SpringBaseType), intent(inout) :: spr
    logical             , intent(in)    :: linearStatic
    real(dp)            , intent(in)    :: time, dTime
    real(dp)            , intent(inout) :: Estr, Einp, Edmp
    integer             , intent(inout) :: ierr

    !! Local variables
    integer  :: lerr
    real(dp) :: forceMean, deltaExt, deltaInt, deltaYield, scale

    !! --- Logic section ---

    if (.not. spr%isActive) return

    lerr = ierr

    !! Elongation of spring externally (L - Lprev)
    deltaExt = spr%length - spr%lengthPrev

    !! Elongation of spring internally (L0 - L0prev)
    deltaInt = spr%deflection - spr%deflectionPrev - deltaExt

    if (linearStatic) then
       !! Calculate input energy
       spr%Einp = spr%force*deltaInt
       spr%Eyield = 0.0_dp
    else if (dTime > 0.0_dp) then

       !! Mean force between previous and current
       forceMean = 0.5_dp*(spr%force + spr%forcePrev)

       !! Integrate input energy
       spr%Einp = spr%Einp + forceMean*deltaInt

       !! Integrate energy loss due to yielding
       deltaYield = spr%yieldDeflection - spr%yieldDeflectionPrev
       spr%Eyield = spr%Eyield + forceMean*deltaYield
       call checkLoss ('Spring',spr%id,time,Edmp,forceMean,deltaYield)

    else
       spr%Einp = 0.0_dp
       spr%Eyield = 0.0_dp
    end if

    !! Compute the elastic strain energy

    if (spr%deflection >= 0.0_dp) then
       if (associated(spr%stiffScaleEnginePos)) then
          scale = EngineValue(spr%stiffScaleEnginePos,ierr)
       else
          scale = 1.0_dp
       end if
    else if (associated(spr%stiffScaleEngineNeg)) then
       scale = EngineValue(spr%stiffScaleEngineNeg,ierr)
    else
       scale = 1.0_dp
    end if

    if (abs(scale) <= 1.0e-16_dp .or. spr%hasFailed) then
       spr%Estr = 0.0_dp
    else
       deltaInt = spr%deflection - spr%yieldDeflection
       spr%Estr = scale * springEnergy(spr,deltaInt,ierr)
    end if

    if (ierr < lerr) then
       call reportError (error_p,'Failed to compute strain energy for Spring'//&
            &            getId(spr%id),addString='SpringEnergies')
    end if

    !! Update mechanism energies
    Estr = Estr + spr%eStr
    Einp = Einp + spr%eInp
    Edmp = Edmp + spr%eYield

    !! The previous spring state values are updated for the next time step
    !! in subroutine updateAtConvergence invoked from terminateStep

  end subroutine SpringEnergies


  !!============================================================================
  !> @brief Calculates the incremental work done by a discrete damper.
  !>
  !> @param dmp Damper to calculate energies for
  !> @param[in] time Current simulation time
  !> @param[in] dTime Time increment size
  !> @param Edmp Energy loss due to damping for the mechanism
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 7 Dec 1999

  subroutine DamperEnergies (dmp,time,dTime,Edmp)

    use DamperTypeModule, only : DamperBaseType, dp

    type(DamperBaseType), intent(inout) :: dmp
    real(dp)            , intent(in)    :: time, dTime
    real(dp)            , intent(inout) :: Edmp

    !! Local variables
    real(dp) :: forceMean, delta

    !! --- Logic section ----

    if (.not. dmp%isActive) return

    if (dTime > 0.0_dp) then

       forceMean = 0.5_dp*(dmp%force + dmp%forcePrev)
       delta     = dmp%length - dmp%lengthPrev
       dmp%eDmp  = dmp%eDmp   + forceMean*delta
       call checkLoss ('Damper',dmp%id,time,Edmp,forceMean,delta)

       !! Update mechanism energies
       Edmp = Edmp + dmp%eDmp

    end if

    !! Update the previous damper values for the next time step
    dmp%lengthPrev = dmp%length
    dmp%forcePrev  = dmp%force

  end subroutine DamperEnergies


  !!============================================================================
  !> @brief Calculates the kinetic and potential energy from a lumped mass.
  !>
  !> @param mass Point mass to calculate energies for
  !> @param[in] gravVec Gravitation vector
  !> @param[in] dTime Time increment size
  !> @param Epot Accumulated potential energy for the mechanism
  !> @param Ekin Accumulated kinetic energy for the mechanism
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 14 Oct 1997

  subroutine MassEnergies (mass,gravVec,dTime,Epot,Ekin)

    use MassTypeModule, only : MassType, dp

    type(MassType), intent(inout) :: mass
    real(dp)      , intent(in)    :: gravVec(3), dTime
    real(dp)      , intent(inout) :: Epot, Ekin

    !! Local variables
    real(dp) :: vel(3), Tmat(3,3)

    !! --- Logic section ---

    if (mass%triad%nDOFs < 1) return

    !! Kinetic energy from translations
    Tmat = transpose(mass%triad%ur(:,1:3))
    vel  = matmul(Tmat,mass%triad%urd(1:3))
    mass%Ekin = 0.5_dp*dot_product(mass%mass*vel,vel)

    if (mass%triad%nDOFs == 6) then
       !! Kinetic energy from rotations
       vel = matmul(Tmat,mass%triad%urd(4:6))
       mass%Ekin = mass%Ekin + 0.5_dp*dot_product(matmul(mass%II,vel),vel)
    end if

    if (.not. mass%addedMass) then
       !! Potential energy from lumped masses
       vel = matmul(Tmat,mass%triad%ur(:,4)) ! position in local directions
       mass%Epot = -dot_product(mass%mass*matmul(Tmat,gravVec),vel)
       if (dTime > 0.0_dp) mass%Epot = mass%Epot - mass%Epot0
    end if

    !! Update mechanism energies
    Epot = Epot + mass%ePot
    Ekin = Ekin + mass%eKin

  end subroutine MassEnergies


  !!============================================================================
  !> @brief Calculate the incremental work done by a joint friction.
  !>
  !> @param friction Friction to calculate energies for
  !> @param[in] time Current simulation time
  !> @param[in] dTime Time increment size
  !> @param Edmp Energy loss due to damping for the mechanism
  !>
  !> @callgraph @callergraph
  !>
  !> @author Karl Erik Thoresen
  !> @date 1 Apr 1998
  !>
  !> @author Bjorn Haugen
  !> @date 7 Dec 1999

  subroutine FrictionEnergies (friction,time,dTime,Edmp)

    use FrictionTypeModule, only : FrictionType, dp

    type(FrictionType), intent(inout) :: friction
    real(dp)          , intent(in)    :: time, dTime
    real(dp)          , intent(inout) :: Edmp

    !! Local variables
    real(dp) :: fricMean, delta

    !! --- Logic section ---

    if (.not. friction%isActive) return

    if (dTime > 0.0_dp) then

       if (associated(friction%spr)) then
          friction%eDmp = friction%spr%eYield
       else
          fricMean      = 0.5_dp*(friction%force + friction%forcePrev)
          delta         = friction%pos - friction%posPrev
          friction%eDmp = friction%eDmp + fricMean*delta
          call checkLoss ('Friction',friction%param%id,time,Edmp,fricMean,delta)

          !! Update mechanism energies
          Edmp = Edmp + friction%eDmp
       end if

    end if

    !! Update the previous friction values for the next time step
    friction%forcePrev = friction%force
    friction%posPrev   = friction%pos

  end subroutine FrictionEnergies


  !!============================================================================
  !> @brief Calculates the incremental work done by the aerodynamic force on a
  !> wind turbine blade node.
  !>
  !> @param aNode Blade node to calculate energies for
  !> @param Einp Accumulated input energy from external forces and motions
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 22 Jan 2010

  subroutine AeroEnergy (aNode,Einp)

    use TriadTypeModule      , only : TriadType, dp
    use WindTurbineTypeModule, only : BladeNode
    use TriadTypeModule      , only : transTriadToGlob, GetTriadIncrement

    type(BladeNode), intent(inout) :: aNode
    real(dp)       , intent(inout) :: Einp

    !! Local variables
    real(dp) :: fa(6), dU(6)

    !! --- Logic section ---

    !! Aerodynamic forces in global coordinates
    fa(1:2) = aNode%aeroForce(1:2)
    fa(3:5) = 0.0_dp
    fa(6)   = aNode%aeroForce(3)
    call transTriadToGlob (aNode%triad,fa)

    !! Compute displacement vector from previous to current configuration
    dU = GetTriadIncrement(aNode%triad)

    !! Incremental energy from force and moment
    Einp = Einp + dot_product(0.5_dp*(fa+aNode%aeroForcePrev),dU)

    !! Update the previous aerodynamic force for the next time step
    aNode%aeroForcePrev = fa

  end subroutine AeroEnergy

end module EnergyIntegrationModule
