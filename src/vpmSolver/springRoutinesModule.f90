!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file springRoutinesModule.f90
!>
!> @brief Subroutines for spring calculations.

!!==============================================================================
!> @brief Module with subroutines for spring calculations.
!>
!> @details This module contains a set of subroutines for performing various
!> computation tasks on the springtypemodule::springtype objects in the model
!> (axial- and joint spring) during the dynamic or quasi-static simulation.

module SpringRoutinesModule

  use kindModule, only : dp

  implicit none

  private

  !> @brief Assembles system force vector contributions from a spring element.
  interface addInSpringForces
     module procedure addInStiffnessForces
     module procedure addInDampingForces
  end interface

  public :: UpdateSprings, UpdateSpringYields
  public :: addInSpringForces, addInSpringStiffMat

  integer , parameter :: maxDofs_p = 60         !< Max number of DOFs pr element
  real(dp), parameter :: epsZero_p = 1.0e-16_dp !< Equal to zero tolerance


contains

  !!============================================================================
  !> @brief Update all springs in the model.
  !>
  !> @param springs Array of all spring elements in the model
  !> @param joints Array of all joints in the model
  !> @param[in] restart If .true., axial spring forces are not updated
  !> @param ierr Error flag
  !> @param[in] updateLength If .true. (default), all spring length variables
  !> are updated first, before updating the other variables
  !> @param[in] updateVar If .true. (default), all spring variables are updated
  !>
  !> @details The length (and direction) of all axial springs are updated first.
  !> Then, all the remaining spring variables are updated.
  !> This is to ensure that sensors measuring spring length are up to date
  !> before they are used by engines. Note that for joint springs, the spring
  !> length coincides with the corresponding joint variable of the owner joint
  !> and thus needs not to be updated here.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 11 June 2002

  subroutine UpdateSprings (springs,joints,restart,ierr,updateLength,updateVar)

    use SpringTypeModule          , only : SpringType
    use TriadTypeModule           , only : updateNodeForce
    use EngineRoutinesModule      , only : updateSpringBase
    use MasterSlaveJointTypeModule, only : MasterSlaveJointType
    use IdTypeModule              , only : getId
#ifdef FT_DEBUG
    use SpringTypeModule          , only : writeObject
    use dbgUnitsModule            , only : dbgSolve
#endif
    use reportErrorModule         , only : reportError, error_p, debugFileOnly_p

    type(SpringType)          , intent(inout) :: springs(:)
    type(MasterSlaveJointType), intent(inout) :: joints(:)
    logical                   , intent(in)    :: restart
    integer                   , intent(inout) :: ierr
    logical, optional         , intent(in)    :: updateLength, updateVar

    !! Local variables
    logical  :: doUpdate, firstRot
    integer  :: i, j, lerr, n
    real(dp) :: theta(3), eV(60)

    !! --- Logic section ---

    lerr = ierr

    if (present(updateLength)) then
       doUpdate = updateLength
    else
       doUpdate = .true.
    end if
    if (doUpdate) then

       !! Update all spring lengths first in case some are used by engines
       do i = 1, size(springs)
          if (.not. associated(springs(i)%triads)) cycle

          firstRot = .true.
          do j = 1, size(springs(i)%spr)
             if (.not. associated(springs(i)%spr(j)%p)) cycle
             if (.not. springs(i)%spr(j)%p%isActive) cycle

             if (springs(i)%spr(j)%p%dof == 0) then

                !! Axial spring
                call updateAxialSpringLen (springs(i)%spr(j)%p, &
                     &                     springs(i)%triads, &
                     &                     springs(i)%forceDir)

             else if (size(springs(i)%triads) == 2) then

                !! Global spring
                call updateGlobalSpringLen (springs(i)%spr(j)%p, &
                     &                      springs(i)%triads(1)%p, &
                     &                      springs(i)%triads(2)%p)

             else

                !! Spring to ground
                call updateGroundedSpringLen (springs(i)%spr(j)%p, &
                     &                        springs(i)%triads(1)%p)

             end if

          end do
       end do
    end if

    if (present(updateVar)) then
       doUpdate = updateVar
    else
       doUpdate = .true.
    end if
    if (doUpdate) then

       !! Now update all spring variables
       do i = 1, size(springs)
          do j = 1, size(springs(i)%spr)
             if (associated(springs(i)%spr(j)%p)) then
                call updateSpringBase (springs(i)%spr(j)%p,ierr)
             end if
          end do
          if (associated(springs(i)%couplStiff)) then
             call coupleSpringForces (springs(i)%spr,springs(i)%couplStiff)
          end if
          if (.not.restart .and. associated(springs(i)%triads)) then
             call calcSpringForces (springs(i)%spr,springs(i)%triads, &
                  &                 springs(i)%forceDir,eV)
             n = 1
             do j = 1, size(springs(i)%triads)
                call updateNodeForce (springs(i)%triads(j)%p, &
                     &                eV(n:n+springs(i)%triads(j)%p%nDOFs-1))
                n = n + springs(i)%triads(j)%p%nDOFs
             end do
          end if
#ifdef FT_DEBUG
          if (dbgSolve > 0) then
             do j = 1, size(springs(i)%spr)
                if (associated(springs(i)%spr(j)%p)) then
                   call writeObject (springs(i)%spr(j)%p,dbgSolve,2)
                end if
             end do
          end if
#endif
       end do
       do i = 1, size(joints)
          if (associated(joints(i)%springEl)) then
             do j = 1, size(joints(i)%springEl%spr)
                if (associated(joints(i)%springEl%spr(j)%p)) then
                   call updateSpringBase (joints(i)%springEl%spr(j)%p,ierr)
                end if
             end do
             if (associated(joints(i)%springEl%lConn)) then
                call updateSpringInterConn (joints(i)%springEl,.true.,ierr)
             end if
          end if
          if (associated(joints(i)%sprFric)) then
             do j = 1, size(joints(i)%sprFric%spr)
                if (associated(joints(i)%sprFric%spr(j)%p)) then
                   call updateSpringBase (joints(i)%sprFric%spr(j)%p,ierr)
                end if
             end do
          end if
       end do

    end if

    if (ierr < lerr) call reportError (debugFileOnly_p,'updateSprings')

  contains

    !> @brief Updates the length and direction of an axial spring.
    subroutine updateAxialSpringLen (spr,triads,sDir)

      use kindModule      , only : epsDiv0_p
      use TriadTypeModule , only : TriadPtrType
      use SpringTypeModule, only : SpringBaseType

      type(SpringBaseType), intent(inout) :: spr
      type(TriadPtrType)  , intent(in)    :: triads(:)
      real(dp)            , intent(out)   :: sDir(:,:)

      spr%length = 0.0_dp
      do n = 1, size(triads)-1

         sDir(1:3,n) = triads(n+1)%p%ur(:,4) - triads(n)%p%ur(:,4)
         sDir(4,n)   = dot_product(sDir(1:3,n),sDir(1:3,n))
         if (sDir(4,n) < epsDiv0_p) then
            call reportError (error_p,'Axial spring'//getId(spr%id), &
                 'The two points defining the spring direction coincide.', &
                 'Triad'//trim(getId(triads(n)%p%id))// &
                 'and Triad'//trim(getId(triads(n+1)%p%id)))
            sDir(:,n) = 0.0_dp
            ierr = ierr - 1
         else
            sDir(4,n)   = sqrt(sDir(4,n))
            sDir(1:3,n) = sDir(1:3,n)/sDir(4,n)
            spr%length  = spr%length + sDir(4,n)
         end if

      end do

#ifdef FT_DEBUG
      if (dbgSolve > 0) then
         write(dbgSolve,600) trim(getId(spr%id)),(sDir(:,n),n=1,size(sDir,2))
         if (size(sDir,2) > 1) write(dbgSolve,610) spr%length
600      format('Force direction and length for Spring',A,':' &
              /('   Dir =',1P3E13.5,'   Length =',E13.5))
610      format(44X,'Total Length =',1PE13.5)
      end if
#endif
    end subroutine updateAxialSpringLen

    !> @brief Updates the length of a global spring.
    subroutine updateGlobalSpringLen (spr,triad1,triad2)
      use TriadTypeModule , only : TriadType
      use SpringTypeModule, only : SpringBaseType
      type(SpringBaseType), intent(inout) :: spr
      type(TriadType)     , intent(in)    :: triad1, triad2
      if (spr%dof < 4) then
         spr%length = triad2%ur(spr%dof,4) - triad1%ur(spr%dof,4)
      else
         if (firstRot) then
            call FFa_eulerZYX (triad1%ur(:,1:3),triad2%ur(:,1:3),theta)
            firstRot = .false.
         end if
         spr%length = theta(spr%dof-3)
      end if
    end subroutine updateGlobalSpringLen

    !> @brief Updates the length of a grounded spring.
    subroutine updateGroundedSpringLen (spr,triad1)
      use TriadTypeModule , only : TriadType
      use SpringTypeModule, only : SpringBaseType
      type(SpringBaseType), intent(inout) :: spr
      type(TriadType)     , intent(in)    :: triad1
      if (spr%dof < 4) then
         spr%length = triad1%ur(spr%dof,4)
      else
         if (firstRot) then
            call FFa_glbEulerZYX (triad1%ur(:,1:3),theta)
            firstRot = .false.
         end if
         spr%length = theta(spr%dof-3)
      end if
    end subroutine updateGroundedSpringLen

  end subroutine UpdateSprings


  !!============================================================================
  !> @brief Corrects spring force and stiffness according to interconnectivity.
  !>
  !> @param springEl The spring element to update stiffness and force for
  !> @param[in] includeStressStiff If .true., account for geometric stiffness
  !> @param ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !> @date 27 May 2004
  !>
  !> @date Knut Morten Okstad
  !> @date 27 Jun 2014

  subroutine UpdateSpringInterConn (springEl,includeStressStiff,ierr)

    use kindModule          , only : epsDiv0_p
    use SpringTypeModule    , only : SpringType, springForce, springStiff
    use EngineRoutinesModule, only : EngineValue
    use manipMatrixModule   , only : dyadic_product
#ifdef FT_DEBUG
    use IdTypeModule        , only : getId
    use manipMatrixModule   , only : writeObject
    use dbgUnitsModule      , only : dbgSolve
    use reportErrorModule   , only : note_p
#endif
    use reportErrorModule   , only : reportError, error_p

    type(SpringType), intent(inout) :: springEl
    logical         , intent(in)    :: includeStressStiff
    integer         , intent(inout) :: ierr

    !! Local variables
    integer  :: i, j, lerr, nDOFs
    real(dp) :: cplDef(2), dirCos(2,6), cplForce(2), cplStiff(2), sgnDef, intPol

    !! --- Logic section ---

    nDOFs = size(springEl%spr)

    cplDef   = 0.0_dp
    dirCos   = 0.0_dp
    cplForce = 0.0_dp
    cplStiff = 0.0_dp

    do i = 1, 2

       !! Establish the interconnected (radial) deflections
       do j = 1, nDOFs
          if (springEl%lConn(i,j)) then
             cplDef(i) = cplDef(i) + springEl%spr(j)%p%deflection**2.0_dp
          end if
       end do
       cplDef(i) = sqrt(cplDef(i))

       if (cplDef(i) > epsDiv0_p) then
          !! Establish directional cosines
          do j = 1, nDOFs
             if (springEl%lConn(i,j)) then
                dirCos(i,j) = springEl%spr(j)%p%deflection / cplDef(i)
             end if
          end do
       else if (any(springEl%lConn(i,:))) then
          cplDef(i) = -1.0_dp ! No interconnection when close to zero
       end if

    end do
#ifdef FT_DEBUG
    if (any(cplDef < 0.0_dp)) then
       call reportError (note_p,'Disabled spring interconnections for Joint'// &
            &            trim(getId(springEl%id))//' due to zero deflections.')
    end if
    if (dbgSolve > 0 .and. any(springEl%lConn(1,:))) then
       write(dbgSolve,100) trim(getId(springEl%id)),'Tra',cplDef(1),dirCos(1,:)
    end if
    if (dbgSolve > 0 .and. any(springEl%lConn(2,:))) then
       write(dbgSolve,100) trim(getId(springEl%id)),'Rot',cplDef(2),dirCos(2,:)
    end if
100 format(/'Spring',A / 5X,'cpl',A3,'   = ',1PE12.5 / 5X,'dirCos   =',1P6E13.5)
#endif

    lerr = 0
    do i = 1, 2
       if (cplDef(i) < 0.0_dp) cycle ! No interconnection

       !! Establish the interpolated force and stiffness
       do j = 1, nDOFs
          if (abs(dirCos(i,j)) <= epsZero_p) cycle

          !! Only elastic springs are implemented so far.
          !! We would need to account for the current yield deflection
          !! in the below calculation of the coupled (elastic) deflection.
          if (associated(springEl%spr(j)%p%yield)) lerr = lerr + 1
          if (springEl%spr(j)%p%unLoadType > 0) lerr = lerr + 1

          !! Possible scaling of the stiffness and force
          if (springEl%spr(j)%p%deflection >= 0.0_dp) then
             if (associated(springEl%spr(j)%p%stiffScaleEnginePos)) then
                intPol = EngineValue(springEl%spr(j)%p%stiffScaleEnginePos,ierr)
             else
                intPol = 1.0_dp
             end if
          else
             if (associated(springEl%spr(j)%p%stiffScaleEngineNeg)) then
                intPol = EngineValue(springEl%spr(j)%p%stiffScaleEngineNeg,ierr)
             else
                intPol = 1.0_dp
             end if
          end if

          intPol = intPol*dirCos(i,j)*dirCos(i,j)
          sgnDef = sign(cplDef(i),springEl%spr(j)%p%deflection)
          cplForce(i) = cplForce(i) &
               &      + intPol*springForce(springEl%spr(j)%p,sgnDef,ierr) &
               &      * sign(1.0_dp,springEl%spr(j)%p%deflection)
          cplStiff(i) = cplStiff(i) &
               &      + intPol*springStiff(springEl%spr(j)%p,sgnDef,ierr)
       end do

    end do
    if (lerr > 0) then
       ierr = ierr - lerr
       call reportError (error_p,'Interconnection of non-elastic '// &
            &            'springs not yet implemented, sorry...', &
            &            addString='UpdateSpringInterConn')
       return
    end if

#ifdef FT_DEBUG
    if (dbgSolve > 0) then
       write(dbgSolve,"(5X,'cplForce =',1P2E13.5)") cplForce
       write(dbgSolve,"(5X,'cplStiff =',1P2E13.5)") cplStiff
    end if
#endif

    !! Add the stiffness and force contributions
    springEl%stiffMat = 0.0_dp

    do i = 1, 2
       if (cplDef(i) > 0.0_dp) then
          springEl%stiffMat = springEl%stiffMat + &
               &              dyadic_product(dirCos(i,:))*cplStiff(i)
          if (i == 1 .and. includeStressStiff) then
             !! Add geometric stiffness terms for translational springs
             call addGeoStiffness (springEl%stiffMat(1:3,1:3),dirCos(1,1:3), &
                  &                springEl%lConn(1,:),cplForce(1)/cplDef(1))
          end if
          do j = 1, nDOFs
             if (springEl%lConn(i,j) .and. associated(springEl%spr(j)%p)) then
                springEl%spr(j)%p%force = dirCos(i,j)*cplForce(i)
                springEl%spr(j)%p%stiffness = springEl%stiffMat(j,j)
             end if
          end do
       else
          do j = 3*i-2, min(3*i,nDOFs)
             if (associated(springEl%spr(j)%p)) then
                springEl%stiffMat(j,j) = springEl%spr(j)%p%stiffness
             end if
          end do
       end if
    end do

#ifdef FT_DEBUG
    if (dbgSolve > 0) then
       call writeObject (springEl%stiffMat,dbgSolve, &
            &            'stiffMat for Spring'//trim(getId(springEl%id)))
    end if
#endif

  contains

    !> @brief Adds geometric stiffness terms for a translational spring
    subroutine addGeoStiffness (Kt,forceDir,activeSpr,KgNom)
      real(dp), intent(inout) :: Kt(3,3)
      logical , intent(in)    :: activeSpr(3)
      real(dp), intent(in)    :: forceDir(3), KgNom
      real(dp) :: dir(3)
      if (count(activeSpr) == 3) then
         do j = 1, 3
            Kt(j,j) = Kt(j,j) + KgNom
         end do
         Kt = Kt - KgNom*dyadic_product(forceDir)
      else if (count(activeSpr) == 2) then
         if (.not. activeSpr(1)) then
            dir(1) =  0.0_dp
            dir(2) = -forceDir(3)
            dir(3) =  forceDir(2)
         else if (.not. activeSpr(2)) then
            dir(1) = -forceDir(3)
            dir(2) =  0.0_dp
            dir(3) =  forceDir(1)
         else
            dir(1) = -forceDir(2)
            dir(2) =  forceDir(1)
            dir(3) =  0.0_dp
         end if
         Kt = Kt + KgNom*dyadic_product(dir)
      end if
    end subroutine addGeoStiffness

  end subroutine UpdateSpringInterConn


  !!============================================================================
  !> @brief Updates the spring yield variables.
  !>
  !> @param yields Array of all spring yield objects in the model
  !> @param ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 11 Jul 2005

  subroutine UpdateSpringYields (yields,ierr)

    use SpringTypeModule    , only : SpringYieldType
    use IdTypeModule        , only : getId
    use EngineRoutinesModule, only : Evaluate
    use reportErrorModule   , only : reportError, error_p

    type(SpringYieldType), intent(inout) :: yields(:)
    integer              , intent(inout) :: ierr

    !! Local variables
    integer :: i, lerr

    !! --- Logic section ---

    lerr = ierr
    do i = 1, size(yields)

       call updateYieldLimit (yields(i)%pos)
       call updateYieldLimit (yields(i)%neg)
       if (ierr < lerr) then
          lerr = ierr
          call reportError (error_p,'Failed to evaluate the yield limit '// &
               &            'engine(s) for Spring Yield'//getId(yields(i)%id), &
               &            addString='UpdateSpringYields')
       end if

    end do

  contains

    !> @brief Updates the yield limit by evaluating the yield limit engine.
    subroutine updateYieldLimit (yld)
      use SpringTypeModule, only : YieldLimitType
      type(YieldLimitType), intent(inout) :: yld
      if (yld%active) then
         yld%yieldLimit = Evaluate(yld%engine,yld%ysf(2),yld%ysf(1),ierr)
      end if
    end subroutine updateYieldLimit

  end subroutine UpdateSpringYields


  !!============================================================================
  !> @brief Calculates system force vector contributions from a spring element.
  !>
  !> @param[in] spr Array of all base springs in current (axial) spring element
  !> @param[in] triads Array of triads connected to the axial spring
  !> @param[in] sDir Direction vector for axial spring
  !> @param[out] eV Element force vector
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 31 Oct 2014

  subroutine calcSpringForces (spr,triads,sDir,eV)

    use SpringTypeModule, only : SpringPtrType
    use TriadTypeModule , only : TriadPtrType, transGlobToSys

    type(SpringPtrType), intent(in)  :: spr(:)
    type(TriadPtrType) , intent(in)  :: triads(:)
    real(dp)           , intent(in)  :: sDir(:,:)
    real(dp)           , intent(out) :: eV(:)

    !! Local variables
    integer  :: i, j, ip1, ip2
    real(dp) :: F

    !! --- Logic section ---

    eV = 0.0_dp

    do j = 1, size(spr)
       if (.not. associated(spr(j)%p)) cycle
       if (.not. spr(j)%p%isActive)    cycle

       F = spr(j)%p%Force
       if (spr(j)%p%dof <= 0) then

          !! Axial spring
          ip1 = 1
          do i = 1, size(triads)-1
             ip2 = ip1 + triads(i)%p%nDOFs
             if (triads(i)%p%nDOFs >= 3) then
                eV(ip1:ip1+2) = eV(ip1:ip1+2) - F*sDir(1:3,i)
             end if
             if (triads(i+1)%p%nDOFs >= 3) then
                eV(ip2:ip2+2) = eV(ip2:ip2+2) + F*sDir(1:3,i)
             end if
             ip1 = ip2
          end do

       else if (size(triads) == 2) then

          !! Global spring
          ip1     =  spr(j)%p%dof
          ip2     =  ip1 + triads(1)%p%nDOFs
          eV(ip1) = -F
          eV(ip2) =  F

       else

          !! Grounded global spring
          ip1     =  spr(j)%p%dof
          eV(ip1) =  F

       end if

    end do

    !! Transform the element vector to system directions
    ip1 = 1
    do i = 1, size(triads)
       if (triads(i)%p%nDOFs > 0) then
          ip2 = ip1 + triads(i)%p%nDOFs
          call transGlobToSys (triads(i)%p, eV(ip1:ip2-1))
          ip1 = ip2
       end if
    end do

  end subroutine calcSpringForces


  !!============================================================================
  !> @brief Calculates the stiffness matrix contributions from a spring element.
  !>
  !> @param[in] spr Array of all base springs in current (axial) spring element
  !> @param[in] triads Array of triads connected to the axial spring
  !> @param[in] sDir Direction vector for axial spring
  !> @param[in] cplStiff Explicit coupling stiffness coefficients
  !> @param[in] scaleK Stiffness scaling factor
  !> @param[out] eM Element stiffness matrix
  !> @param[in] includeStressStiff If .true., account for geometric stiffness
  !> @param haveStiffness If .true., a non-zero spring stiffness was detected
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 31 Oct 2014

  subroutine calcSpringStiffness (spr,triads,sDir,cplStiff,scaleK,eM, &
       &                          includeStressStiff,haveStiffness)

    use SpringTypeModule , only : SpringPtrType
    use TriadTypeModule  , only : TriadPtrType, transGlobToSys
    use manipMatrixModule, only : dyadic_product

    type(SpringPtrType), intent(in)    :: spr(:)
    type(TriadPtrType) , intent(in)    :: triads(:)
    real(dp)           , intent(in)    :: sDir(:,:), scaleK
    real(dp), pointer  , intent(in)    :: cplStiff(:)
    real(dp)           , intent(out)   :: eM(:,:)
    logical            , intent(in)    :: includeStressStiff
    logical            , intent(inout) :: haveStiffness

    !! Local variables
    integer  :: i, j, ip1, ip2, jp1, jp2
    real(dp) :: EAoL, KgNom, KmNom, KtMat(3,3)

    !! --- Logic section ---

    eM = 0.0_dp
    KgNom = 0.0_dp

    do j = 1, size(spr)
       if (.not. associated(spr(j)%p)) cycle
       if (.not. spr(j)%p%isActive)    cycle

       haveStiffness = .true.
       KmNom = spr(j)%p%stiffness * scaleK

       if (spr(j)%p%dof <= 0) then

          !! Axial spring
          ip1  = 1
          EAoL = KmNom
          do i = 1, size(triads)-1
             KmNom = EAoL * spr(j)%p%length/sDir(4,i)
             if (includeStressStiff) then
                !! Add geometric stiffness
                KgNom = scaleK * spr(j)%p%force/sDir(4,i)
             end if
             KtMat = dyadic_product(sDir(1:3,i)) * (KmNom-KgNom)
             KtMat(1,1) = KtMat(1,1) + KgNom
             KtMat(2,2) = KtMat(2,2) + KgNom
             KtMat(3,3) = KtMat(3,3) + KgNom
             ip2 = ip1 + triads(i)%p%nDOFs
             jp1 = ip1 + 2
             jp2 = ip2 + 2
             if (triads(i)%p%nDOFs >= 3) then
                eM(ip1:jp1,ip1:jp1) = eM(ip1:jp1,ip1:jp1) + KtMat
             end if
             if (triads(i+1)%p%nDOFs >= 3) then
                eM(ip2:jp2,ip2:jp2) = eM(ip2:jp2,ip2:jp2) + KtMat
             end if
             if (triads(i)%p%nDOFs >= 3 .and. triads(i+1)%p%nDOFs >= 3) then
                eM(ip1:jp1,ip2:jp2) = eM(ip1:jp1,ip2:jp2) - KtMat
                eM(ip2:jp2,ip1:jp1) = eM(ip2:jp2,ip1:jp1) - KtMat
             end if
             ip1 = ip2
          end do

       else if (size(triads) == 2) then

          !! Global spring
          ip1 = spr(j)%p%dof
          ip2 = ip1 + triads(1)%p%nDOFs
          eM(ip1,ip1) =  KmNom
          eM(ip2,ip2) =  KmNom
          eM(ip1,ip2) = -KmNom
          eM(ip2,ip1) = -KmNom

       else

          !! Grounded global spring
          ip1 = spr(j)%p%dof
          eM(ip1,ip1) = KmNom

       end if

    end do
    if (associated(cplStiff)) then

       !! Include explicit coupling stiffness for global springs
       if (size(triads) > 1) then
          call coupleSpringStiffness (spr,cplStiff,eM,scaleK,triads(1)%p%nDOFs)
       else
          call coupleSpringStiffness (spr,cplStiff,eM,scaleK)
       end if

    end if
    if (haveStiffness) then

       !! Transform the element matrix to system directions
       ip1 = 1
       do i = 1, size(triads)
          if (triads(i)%p%nDOFs > 0) then
             call transGlobToSys (triads(i)%p, eM, ip1, triads(i)%p%nDOFs)
             ip1 = ip1 + triads(i)%p%nDOFs
          end if
       end do

    end if

  end subroutine calcSpringStiffness


  !!==========================================================================
  !> @brief Updates spring element forces due to explicit coupling.
  !>
  !> @param spr Array of all base springs in current (axial) spring element
  !> @param[in] cplStiff Explicit coupling stiffness coefficients
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 2 Jul 2010

  subroutine coupleSpringForces (spr,cplStiff)

    use SpringTypeModule, only : SpringPtrType

    type(SpringPtrType), intent(inout) :: spr(:)
    real(dp)           , intent(in)    :: cplStiff(:)

    !! Local variables
    integer :: i, j, k

    !! --- Logic section ---

    k = 0
    do i = 1, 6
       do j = i+1, 6
          k = k + 1
          if (k > size(cplStiff)) exit
          if (j > size(spr) .or. abs(cplStiff(k)) <= epsZero_p) cycle
          if (.not. associated(spr(i)%p)) cycle
          if (.not. associated(spr(j)%p)) cycle

          if (spr(i)%p%isActive) then
             spr(i)%p%Force = spr(i)%p%Force + cplStiff(k)*spr(j)%p%deflection
          end if

          if (spr(j)%p%isActive) then
             spr(j)%p%Force = spr(j)%p%Force + cplStiff(k)*spr(i)%p%deflection
          end if

       end do
       if (k > size(cplStiff)) exit
    end do

  end subroutine coupleSpringForces


  !!============================================================================
  !> @brief Updates a spring element stiffness matrix due to explicit coupling.
  !>
  !> @param spr Array of all base springs in current (axial) spring element
  !> @param[in] cplStiff Explicit coupling stiffness coefficients
  !> @param elMat Element stiffness matrix
  !> @param[in] scaleK Stiffness scaling factor
  !> @param[in] d1 Number of DOFs in the first triad connected to spring element
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 2 Jul 2010

  subroutine coupleSpringStiffness (spr,cplStiff,elMat,scaleK,d1)

    use SpringTypeModule, only : SpringPtrType

    type(SpringPtrType), intent(in)    :: spr(:)
    real(dp)           , intent(in)    :: cplStiff(:), scaleK
    real(dp)           , intent(inout) :: elMat(:,:)
    integer , optional , intent(in)    :: d1

    !! Local variables
    integer :: i, j, k, idof, jdof

    !! --- Logic section ---

    k = 0
    do i = 1, 6
       do j = i+1, 6
          k = k + 1
          if (k > size(cplStiff)) exit
          if (j > size(spr) .or. abs(cplStiff(k)) <= epsZero_p) cycle
          if (.not. associated(spr(i)%p)) cycle
          if (.not. associated(spr(j)%p)) cycle
          if (.not. spr(i)%p%isActive)    cycle
          if (.not. spr(j)%p%isActive)    cycle

          idof = spr(i)%p%dof
          jdof = spr(j)%p%dof
          elMat(idof,jdof) = elMat(idof,jdof) + cplStiff(k)*scaleK
          elMat(jdof,idof) = elMat(jdof,jdof) + cplStiff(k)*scaleK
          if (.not. present(d1)) cycle

          elMat(d1+idof,jdof) = elMat(d1+idof,jdof) - cplStiff(k)*scaleK
          elMat(jdof,d1+idof) = elMat(jdof,d1+idof) - cplStiff(k)*scaleK
          elMat(idof,d1+jdof) = elMat(idof,d1+jdof) - cplStiff(k)*scaleK
          elMat(d1+jdof,idof) = elMat(d1+jdof,idof) - cplStiff(k)*scaleK
          elMat(d1+idof,d1+jdof) = elMat(d1+idof,d1+jdof) + cplStiff(k)*scaleK
          elMat(d1+jdof,d1+idof) = elMat(d1+jdof,d1+idof) + cplStiff(k)*scaleK

       end do
       if (k > size(cplStiff)) exit
    end do

  end subroutine coupleSpringStiffness


  !!============================================================================
  !> @brief Assembles system force vector contributions from a spring element.
  !>
  !> @param FS System stiffness force vector
  !> @param RF System reaction forces associated with constrained DOFs
  !> @param[in] spring The spring element to calculate stiffness forces for
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine calculates the stiffness force contributions to
  !> the system force vector from the spring element, and adds them into
  !> the system vectors @a FS and @a RF.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Karl Erik Thoresen
  !>
  !> @date Feb 1999

  subroutine addInStiffnessForces (FS,RF,spring,sam,ierr)

    use SamModule         , only : SamType
    use SpringTypeModule  , only : SpringType
    use AsmExtensionModule, only : csAddEV
#ifdef FT_DEBUG
    use IdTypeModule      , only : getId, StrId
    use manipMatrixModule , only : writeObject
    use dbgUnitsModule    , only : dbgSolve
#endif
    use reportErrorModule , only : reportError, debugFileOnly_p

    real(dp)        , intent(inout) :: FS(:), RF(:)
    type(SpringType), intent(in)    :: spring
    type(SamType)   , intent(in)    :: sam
    integer         , intent(inout) :: ierr

    !! Local variables
    integer  :: j, jdof, nDOFs, err
    real(dp) :: eV(maxDofs_p)

    !! --- Logic section ---

    nDOFs = spring%nDOFs

    if (associated(spring%triads)) then

       !! Axial or global spring
       call calcSpringForces (spring%spr,spring%triads,spring%forceDir,eV)

    else

       !! Joint spring
       eV = 0.0_dp
       do j = 1, size(spring%spr)
          if (.not. associated(spring%spr(j)%p)) cycle
          if (.not. spring%spr(j)%p%isActive) cycle

          jdof = spring%spr(j)%p%dof
          if (jdof > 0 .and. jdof <= maxDofs_p) then
             eV(jdof) = spring%spr(j)%p%Force
          end if
       end do

    end if

#ifdef FT_DEBUG
    if (dbgSolve > 0) then
       call writeObject (eV(1:nDOFs),dbgSolve, &
            &            'Force vector for Spring'//trim(getId(spring%id))// &
            &            ' (element'//trim(StrId(spring%samElnum))//')')
    end if
#endif

    call csAddEV (sam, spring%samElnum, nDOFs, eV(1:nDOFs), FS, RF, err)
    if (err == 0) return

    ierr = ierr  - 1
    call reportError (debugFileOnly_p,'addInStiffnessForces')

  end subroutine addInStiffnessForces


  !!============================================================================
  !> @brief Assembles system force vector contributions from a spring element.
  !>
  !> @param FD System damping force vector
  !> @param RF System reaction forces associated with constrained DOFs
  !> @param[in] spring The spring element to calculate damping forces for
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[in] alpha2 Stiffness-proportional damping coefficient
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine calculates the damping force contributions to
  !> the system force vector from the spring element assuming Rayleigh damping,
  !> and adds them into the system vectors @a FD and @a RF.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 31 Oct 2014

 subroutine addInDampingForces (FD,RF,spring,sam,alpha2,ierr)

    use SamModule         , only : SamType
    use SpringTypeModule  , only : SpringType
    use TriadTypeModule   , only : transGlobToSys
    use AsmExtensionModule, only : csAddEV
#ifdef FT_DEBUG
    use IdTypeModule      , only : getId, StrId
    use manipMatrixModule , only : writeObject
    use dbgUnitsModule    , only : dbgSolve
#endif
    use reportErrorModule , only : reportError, debugFileOnly_p

    real(dp)        , intent(inout) :: FD(:), RF(:)
    type(SpringType), intent(in)    :: spring
    type(SamType)   , intent(in)    :: sam
    real(dp)        , intent(in)    :: alpha2
    integer         , intent(inout) :: ierr

    !! Local variables
    logical  :: haveStiffness
    integer  :: i, jdof, nDOFs, err
    real(dp) :: eV(maxDofs_p), eC(maxDofs_p,maxDofs_p)

    !! --- Logic section ---

    if (abs(alpha2) <= epsZero_p .or. .not.associated(spring%triads)) return

    !! Calculate stiffness-proportional damping matrix for axial spring
    call calcSpringStiffness (spring%spr,spring%triads,spring%forceDir, &
         &                    spring%couplStiff,alpha2,eC,.false.,haveStiffness)
    if (.not. haveStiffness) return

    !! Extract velocity vector for the spring element
    jdof = 1
    do i = 1, size(spring%triads)
       if (spring%triads(i)%p%nDOFs >= 3) then
          eV(jdof:jdof+2) = spring%triads(i)%p%urd(1:3)
          call transGlobToSys (spring%triads(i)%p, eV(jdof:jdof+2))
       end if
       jdof = jdof + spring%triads(i)%p%nDOFs
    end do

    nDOFs = spring%nDOFs
#ifdef FT_DEBUG
    if (dbgSolve > 0) then
       call writeObject (matmul(eC(1:nDOFs,1:nDOFs),eV(1:nDOFs)),dbgSolve, &
            &            'Damping forces for Spring'//trim(getId(spring%id))// &
            &            ' (element'//trim(StrId(spring%samElnum))//')')
    end if
#endif

    !! Add in damping forces, Fd = C*v
    call csAddEV (sam, spring%samElnum, nDOFs, &
         &        matmul(eC(1:nDOFs,1:nDOFs),eV(1:nDOFs)), FD, RF, err)
    if (err == 0) return

    ierr = ierr  - 1
    call reportError (debugFileOnly_p,'addInDampingForces')

  end subroutine addInDampingForces


  !!============================================================================
  !> @brief Assembles system stiffness contributions from a spring element.
  !>
  !> @param[in] includeStressStiff If .true., include geometric stiffness
  !> @param[in] scaleK Stiffness matrix scaling factor
  !> @param Nmat System newton matrix
  !> @param[in] spring The spring element to calculate stiffness matrix for
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[out] ierr Error flag
  !> @param Rhs System right-hand-side vector associated with the Newton matrix
  !>
  !> @details This subroutine calculates the stiffness matrix for a spring
  !> element and adds it into the system Newton matrix @a Nmat, multiplied by
  !> a scaling factor, @a scaleK.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Karl Erik Thoresen
  !>
  !> @date Feb 1999

  subroutine addInSpringStiffMat (includeStressStiff,scaleK,Nmat, &
       &                          spring,sam,ierr,Rhs)

    use SamModule          , only : SamType
    use SysMatrixTypeModule, only : SysMatrixType
    use SpringTypeModule   , only : SpringType
    use AsmExtensionModule , only : csAddElM
#ifdef FT_DEBUG
    use IdTypeModule       , only : getId, StrId
    use manipMatrixModule  , only : writeObject
    use dbgUnitsModule     , only : dbgSolve
#endif
    use reportErrorModule  , only : reportError, debugFileOnly_p

    logical            , intent(in)    :: includeStressStiff
    real(dp)           , intent(in)    :: scaleK
    type(SysMatrixType), intent(inout) :: Nmat
    type(SpringType)   , intent(in)    :: spring
    type(SamType)      , intent(in)    :: sam
    integer            , intent(out)   :: ierr
    real(dp), optional , intent(inout) :: Rhs(:)

    !! Local variables
    logical  :: haveStiffness
    integer  :: i, j, iDof, jDof, nDOFs
    real(dp) :: elMat(maxDofs_p,maxDofs_p)

    !! --- Logic section ---

    ierr = 0
    nDOFs = spring%nDOFs
    haveStiffness = .false.

    if (associated(spring%lConn)) then

       !! Interconnected joint spring
       elMat = 0.0_dp
       do i = 1, nDOFs
          do j = 1, nDOFs
             if (any(spring%lConn(:,i)) .and. any(spring%lConn(:,j))) then
                haveStiffness = .true.
                iDof = spring%spr(i)%p%dof
                jDof = spring%spr(j)%p%dof
                elMat(iDof,jDof) = spring%stiffMat(i,j) * scaleK
             end if
          end do
       end do

    else if (associated(spring%triads)) then

       !! Axial or global spring
       call calcSpringStiffness (spring%spr,spring%triads,spring%forceDir, &
            &                    spring%couplStiff,scaleK,elMat, &
            &                    includeStressStiff,haveStiffness)

    else

       !! Joint spring without interconnection
       elMat = 0.0_dp
       do j = 1, size(spring%spr)
          if (.not. associated(spring%spr(j)%p)) cycle
          if (.not. spring%spr(j)%p%isActive) cycle

          haveStiffness = .true.
          iDof = spring%spr(j)%p%dof
          elMat(iDof,iDof) = elMat(iDof,iDof) + &
               &             spring%spr(j)%p%stiffness * scaleK
       end do

    end if

    if (.not. haveStiffness) return

#ifdef FT_DEBUG
    if (dbgSolve > 0) then
       call writeObject (elMat(1:nDOFs,1:nDOFs),dbgSolve,'Stiffness matrix'// &
            &            ' for Spring'//trim(getId(spring%id))// &
            &            ' (element'//trim(StrId(spring%samElnum))//')')
    end if
#endif

    call csAddElM (sam, spring%samElnum, nDOFs, &
         &         elMat(1:nDOFs,1:nDOFs), Nmat, ierr, Rhs)
    if (ierr < 0) call reportError (debugFileOnly_p,'addInSpringStiffMat')

  end subroutine addInSpringStiffMat

end module SpringRoutinesModule
