!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module FrictionRoutinesModule

  implicit none

  private

  public :: updateFrictions, addInFrictionForces, findForceInFriction

  interface findForceInFriction
     module procedure findJointForceInFriction
     module procedure findContactElmForceInFriction
  end interface


contains

  subroutine updateFrictions (joints, cElems, dt, iter, useRealVel, ierr)

    !!==========================================================================
    !! Updates the friction forces in all joints and contact elements.
    !!
    !! Programmer : Knut Morten Okstad               date/rev: 30 May 2008 / 1.0
    !!==========================================================================

    use kindModule                , only : dp
    use MasterSlaveJointTypeModule, only : MasterSlaveJointType
    use ContactElementTypeModule  , only : ContactElementType
    use reportErrorModule         , only : reportError, debugFileOnly_p
    use FFaCmdLineArgInterface    , only : ffa_cmdlinearg_intValue

    type(MasterSlaveJointType), intent(inout) :: joints(:)
    type(ContactElementType)  , intent(inout) :: cElems(:)
    real(dp)                  , intent(in)    :: dt
    integer                   , intent(in)    :: iter
    logical                   , intent(in)    :: useRealVel
    integer                   , intent(inout) :: ierr

    !! Local variables
    integer  :: i, j, lerr
    real(dp) :: deltaT

    !! --- Logic section ---

    if (iter < ffa_cmdlinearg_intValue('nFricSkip')) then
       return
    else if (.not. useRealVel) then
       deltaT = dt ! Compute velocity based on position increment
    else if (iter < ffa_cmdlinearg_intValue('nFricVelUpdateSkip')) then
       deltaT = dt ! Compute velocity based on position increment
    else
       deltaT = 0.0_dp ! Use actual velocity
    end if

    lerr = ierr
    do i = 1, size(joints)
       if (associated(joints(i)%multiDofFriction)) then
          !! Handle multi-dof friction in separate routine. In this case the
          !! joint DOF friction is governed by this joint-wide master friction.
          call updateMultiDofFriction (ierr, dt, deltaT, &
               &                       joints(i)%multiDofFriction, joints(i))
       else
          !! Potentially one friction for each joint DOF
          do j = 1, joints(i)%nJointDOFs
             if (associated(joints(i)%jointDofs(j)%friction)) then
                call updateFriction (ierr, dt, deltaT, &
                     &               joints(i)%jointDofs(j)%friction, &
                     &               jDof=j, joint=joints(i))
             end if
          end do
       end if
    end do

    do i = 1, size(cElems)
       if (associated(cElems(i)%friction)) then
          call updateFriction (ierr, dt, deltaT, &
               &               cElems(i)%friction, cElem=cElems(i))
       end if
    end do

    if (ierr < lerr) call reportError (debugFileOnly_p,'updateFrictions')

  end subroutine updateFrictions


  subroutine updateFriction (ierr, dt, deltaT, friction, jDof, joint, cElem)

    !!==========================================================================
    !! Updates friction forces in the given joint or contact element.
    !!
    !! Programmer : Hans Petter Hildre                         date: 20 Dec 1990
    !! Revised    : Karl E. Thoresen / Bjorn Haugen (to F90)   date: around 2000
    !! Revised    : Knut Morten Okstad (completely rewritten)  date: 28 Aug 2002
    !!==========================================================================

    use MasterSlaveJointTypeModule, only : MasterSlaveJointType
    use MasterSlaveJointTypeModule, only : transVSlaveToJoint
    use ContactElementTypeModule  , only : ContactElementType
    use SpringTypeModule    , only : SpringBaseType, setYieldLimit
    use DamperTypeModule    , only : DamperBaseType
    use FrictionTypeModule  , only : FrictionType
    use FrictionTypeModule  , only : fricForceTol, UpdateFrictionAtStart
    use FrictionTypeModule  , only : ROT_FRICTION_p, TRANS_FRICTION_p
    use FrictionTypeModule  , only : CAM_FRICTION_p, GENERIC_ENGINE_p
    use FrictionTypeModule  , only : BALL_JNT_FRICTION2_p
    use TriadTypeModule     , only : transSysToGlob
    use EngineRoutinesModule, only : UpdateSpringBase, EngineValue
    use kindModule          , only : dp, epsDiv0_p
    use reportErrorModule   , only : internalError
#ifdef FT_DEBUG
    use IdTypeModule        , only : StrId
    use fileUtilitiesModule , only : getDBGfile
    use dbgUnitsModule      , only : dbgFric
#endif

    integer                   , intent(inout)        :: ierr
    real(dp)                  , intent(in)           :: dt, deltaT
    type(FrictionType)        , intent(inout)        :: friction
    integer                   , intent(in), optional :: jDof
    type(MasterSlaveJointType), intent(in), optional :: joint
    type(ContactElementType)  , intent(in), optional :: cElem

    !! Local variables
    type(SpringBaseType), pointer :: spring
    type(DamperBaseType), pointer :: damper

    integer  :: i, j, k, fricType
    real(dp) :: a, R, Y, RF(6), RFX1, RFY1, RFX2, RFY2, RFA, RFR

    !! --- Logic section ---

    if (.not. friction%isActive) then
       friction%force       = 0.0_dp
       friction%forcePrevIt = 0.0_dp
       return
    end if

    i  = 1
    j  = 2
    k  = 3
    RF = 0.0_dp

    if (present(joint) .and. present(jDof)) then

#ifdef FT_DEBUG
       dbgFric = getDBGfile(9,'friction.dbg')
#endif
       !! Find in which direction the friction force acts (i)      i
       !! and the two orthogonal directions (j) and (k)           / \
       j = joint%jointDofs(jDof)%lDof                 !          k - j
       if (j == 7) then ! slider dof (=Tz)
          i = 3
       else
          i = 1 + mod(j-1,3)
       end if
       j = 1 + mod(i,3)
       k = 1 + mod(j,3)

       !! Slave/joint forces. Note that the nodeForce values are referred to
       !! the system directions of the triad. We must transform them to the
       !! local joint coordinate system before calculating equivalent force.
       RF(1:joint%STriad%nDOFs) = joint%STriad%nodeForce
       call transSysToGlob (joint%STriad,RF(1:joint%STriad%nDOFs))
       RF = transVSlaveToJoint(joint,RF(1:joint%STriad%nDOFs))

       friction%pos = joint%jointDofs(jDof)%jvar(1)
       friction%vel = joint%jointDofs(jDof)%jvar(2)

    else if (present(cElem)) then

#ifdef FT_DEBUG
       !! Use a separate debug file for each contact element
       dbgFric = getDBGfile(20+cElem%id%userid,'friction_'// &
            &               trim(adjustl(StrId(cElem%id%userid)))//'.dbg')
#endif

       friction%posPrev = cElem%cVarPrev(3) ! Bugfix #392
       friction%pos = cElem%cVar(3,1)
       friction%vel = cElem%cVar(3,2)

    end if

    if (deltaT > epsDiv0_p) then
       friction%vel = (friction%pos - friction%pos0) / deltaT
    end if

    !! Calculation of equivalent joint force
    if (associated(friction%eng)) then
       fricType = GENERIC_ENGINE_p
    else
       fricType = friction%param%type
    end if
    select case (fricType)

    case (ROT_FRICTION_p)

       R = friction%param%typeDepParams(1) ! Bearing radius
       a = friction%param%typeDepParams(2) ! Bearing width
       Y = friction%param%typeDepParams(3) ! Bearing constant

       if (a < epsDiv0_p) then
          RFR  = sqrt(RF(j)*RF(j) + RF(k)*RF(k))
       else
          RFX1 = 0.5_dp*RF(j) - RF(3+k)/a
          RFY1 = 0.5_dp*RF(k) - RF(3+j)/a
          RFX2 = 0.5_dp*RF(j) + RF(3+k)/a
          RFY2 = 0.5_dp*RF(k) + RF(3+j)/a
          RFR  = sqrt(RFX1*RFX1 + RFY1*RFY1) + sqrt(RFX2*RFX2 + RFY2*RFY2)
       end if

       if (R > 0.0_dp) then
          if (Y > 0.0_dp) then
             RFA = abs(RF(i))
             friction%Fequ = (RFR + RFA*Y)*R ! Equivalent torque
          else
             friction%Fequ = RFR*R ! Equivalent torque
          end if
       end if

    case (TRANS_FRICTION_p)

       R = friction%param%typeDepParams(1) ! Distance to locking device
       Y = friction%param%typeDepParams(3) ! Bearing constant

       if (R < epsDiv0_p) then
          RFA = 0.0_dp
       else
          RF(j) = RF(j) - RF(3+i)/R
          RFA = RF(3+i)/R
       end if
       RFR = sqrt(RF(j)*RF(j) + RF(k)*RF(k))

       if (Y > 0.0_dp) then
          friction%Fequ = abs(RFR + RFA*Y) ! Equivalent force
       else
          friction%Fequ = RFR
       end if

    case (CAM_FRICTION_p)

       do j = 1, 2
          if (present(joint)) then
             spring => joint%jointDofs(j)%spring
             damper => joint%jointDofs(j)%damper
          else
             spring => cElem%springs(j)%p
             damper => cElem%dampers(j)%p
          end if
          if (associated(spring)) then
             RF(j) = spring%force
             !! Add damper contribution only if the spring force is non-zero
             if (abs(RF(j)) > 1.0e-16_dp .and. associated(damper)) then
                RF(j) = RF(j) + damper%force
             end if
          end if
       end do

       friction%Fequ = sqrt(RF(1)*RF(1) + RF(2)*RF(2)) ! Equivalent force

    case (BALL_JNT_FRICTION2_p)

       friction%Fequ = sqrt(RF(1)*RF(1) + RF(2)*RF(2) + RF(3)*RF(3))

    case (GENERIC_ENGINE_p)

       friction%Fequ = EngineValue(friction%eng,ierr)
       if (friction%param%type == ROT_FRICTION_p) then
          R = friction%param%typeDepParams(1) ! Bearing radius
          friction%Fequ = friction%Fequ * R   ! Equivalent torque
       end if

    case default
       ierr = internalError('updateFriction: Invalid friction type')
       return
    end select

    call FindMaxFrictionForce (friction%Fmax, friction%Fequ, friction%Fext, &
         &                     friction%vel, friction%lSlip(1), &
         &                     friction%param%PrestressLoad, &
         &                     friction%param%CoulombCoeff, &
         &                     friction%param%StribeckMagn, &
         &                     friction%param%StribeckSpeed)

    if (associated(friction%spr)) then ! Spring-based friction model

       !bh: probably better to do this only at the end of the step?
       call setYieldLimit (friction%spr,abs(friction%Fmax))

       !! Update the friction spring force
       call updateSpringBase (friction%spr,ierr)
       friction%force = friction%spr%force

    else if (friction%lInit == 2) then

       call UpdateFrictionAtStart (friction)

    else

       call FindFrictionForce (friction%force, friction%stiff, &
            &                  friction%posPrevIt, friction%velPrevIt, &
            &                  friction%forcePrevIt, friction%df, &
            &                  friction%lSlip(2), friction%lInit, dt, &
            &                  friction%Fmax, friction%Fext, fricForceTol, &
            &                  friction%pos0, friction%pos, friction%vel)

    end if

  end subroutine updateFriction


  subroutine updateMultiDofFriction (ierr, dt, deltaT, friction, joint)

    !!==========================================================================
    !! Updates multidof friction forces in the given joint.
    !!
    !! Programmer : Bjorn Haugen                               date: 25 Apr 2003
    !!==========================================================================

    use kindModule                , only : dp, epsDiv0_p
    use FrictionTypeModule        , only : FrictionType, BALL_JNT_FRICTION_p
    use FrictionTypeModule        , only : fricForceTol, UpdateFrictionAtStart
    use MasterSlaveJointTypeModule, only : MasterSlaveJointType
    use reportErrorModule         , only : internalError

    integer                   , intent(inout) :: ierr
    real(dp)                  , intent(in)    :: dt, deltaT
    type(FrictionType)        , intent(inout) :: friction
    type(MasterSlaveJointType), intent(inout) :: joint

    !! Local variables
    integer  :: i, iFree, iMat, lRot
    real(dp) :: pos, RFX, RFY, RFZ
    real(dp) :: omegaVec(3), momentVec(3), Rmat(3,3), rotAxis(3)

    !! --- Logic section ---

    if (.not. friction%isActive) then
       friction%force       = 0.0_dp
       friction%forcePrevIt = 0.0_dp
       do iFree = 1, joint%nJointDOFs
          if (associated(joint%jointDofs(iFree)%friction)) then
             joint%jointDofs(iFree)%friction%isActive    = .false.
             joint%jointDofs(iFree)%friction%force       = 0.0_dp
             joint%jointDofs(iFree)%friction%forcePrevIt = 0.0_dp
          end if
       end do
       return
    end if

    !! Slave/joint forces, directions do not matter here
    RFX = joint%STriad%nodeForce(1)
    RFY = joint%STriad%nodeForce(2)
    RFZ = joint%STriad%nodeForce(3)

    !! Calculation of equivalent joint force
    select case (friction%param%type)

    case (BALL_JNT_FRICTION_p)

       pos = 0.0_dp
       do iFree = 1, joint%nJointDOFs
          pos = pos + joint%jointDofs(iFree)%jvar(1)**2.0_dp
       end do
       friction%pos = sqrt(pos)
       !bh friction%pos = abs(joint%jointDofs(1)%jvar(1)) !TODO,bh: testing

       !! Calculate friction force based on relative velocity and joint force
       if (joint%JMTriads(1)%triad%nDOFs >= 6) then
          omegaVec = joint%STriad%urd(4:6) - joint%JMTriads(1)%triad%urd(4:6)
       else
          omegaVec = joint%STriad%urd(4:6) ! Master is triad grounded
       end if

       friction%vel = sqrt(dot_product(omegaVec,omegaVec))
       !bh friction%vel = abs(joint%jointDofs(1)%jvar(2)) !TODO,bh: testing

       friction%Fequ = sqrt(RFX*RFX + RFY*RFY + RFZ*RFZ) ! Equivalent force

    case default
       ierr = ierr + internalError('updateFriction: Invalid friction type')
       return
    end select

    if (deltaT > epsDiv0_p) then
       friction%vel = (friction%pos - friction%pos0) / deltaT
    end if

    call FindMaxFrictionForce (friction%Fmax, friction%Fequ, friction%Fext, &
         &                     friction%vel, friction%lSlip(1), &
         &                     friction%param%PrestressLoad, &
         &                     friction%param%CoulombCoeff, &
         &                     friction%param%StribeckMagn, &
         &                     friction%param%StribeckSpeed)

    !!TODO,bh: Handle friction spring here

    if (friction%lInit == 2) then
       call UpdateFrictionAtStart (friction)
       do iFree = 1, joint%nJointDOFs
          joint%jointDofs(iFree)%friction%pos = friction%pos
          joint%jointDofs(iFree)%friction%vel = friction%vel
          call UpdateFrictionAtStart (joint%jointDofs(iFree)%friction)
       end do
       return
    else
       call FindFrictionForce (friction%force, friction%stiff, &
            &                  friction%posPrevIt, friction%velPrevIt, &
            &                  friction%forcePrevIt, friction%df, &
            &                  friction%lSlip(2), friction%lInit, dt, &
            &                  friction%Fmax, friction%Fext, fricForceTol, &
            &                  friction%pos0, friction%pos, friction%vel)
    end if

    !! Set the friction force on all the joint DOFs

    if (friction%vel > epsDiv0_p) then
       !! We have motion and apply friciton force in the direction
       !! opposing the rotation
       omegaVec = omegaVec/friction%vel
    else
       !! No motion, apply stick friction force opposing the triad forces
       omegaVec = joint%STriad%nodeForce(4:6)
       RFZ = sqrt(dot_product(omegaVec,omegaVec))
       if (RFZ > epsDiv0_p) then
          omegaVec = omegaVec/RFZ
       else
          omegaVec = 0.0_dp
       end if
    end if

    !! Resultant global moment vector to be applied at the joint rotation DOFs
    momentVec = abs(friction%force)*omegaVec

    do iFree = 1, joint%nJointDOFs

       !! Establish global direction of each rotational joint DOF
       !! Set the position of slave triad

       iMat = joint%jointDofs(iFree)%iMat
       lRot = joint%jointDofs(iFree)%lDof - 3
       if (lRot < 1) cycle

       Rmat = joint%STPos0InJ(:,1:3)
       do i = 1, joint%nMats
          Rmat = matmul(joint%PosFromJVars(:,1:3,i),Rmat)
       end do

       rotAxis = matmul(joint%JPosInG(:,1:3),Rmat(lRot,:))

       joint%jointDofs(iFree)%friction%force = dot_product(rotAxis,momentVec)

       !! Set position and velocity of jdof friction for energy calculations
       joint%jointDofs(iFree)%friction%pos  = joint%jointDofs(iFree)%jvar(1)
       joint%jointDofs(iFree)%friction%vel  = joint%jointDofs(iFree)%jvar(2)
       joint%jointDofs(iFree)%friction%Fequ = friction%Fequ
       joint%jointDofs(iFree)%friction%Fmax = friction%Fmax

    end do

  end subroutine updateMultiDofFriction


  subroutine FindMaxFrictionForce (Fmax, Fequ, Fext, Vel, StickSlip, &
       &                           F0, Coulomb, Stribeck, Vc)

    !!==========================================================================
    !! Calculates the maximum friction force.
    !!
    !! Programmer: Karl Erik Thoresen                          date: 11 Mar 1998
    !! Revised:    Knut Morten Okstad (F90 version)            date: 09 Jan 2001
    !! Revised:    Knut Morten Okstad (completely rewritten)   date: 22 Aug 2002
    !!==========================================================================

    use kindModule            , only : dp, epsDiv0_p
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getint

    real(dp), intent(out)   :: Fmax
    integer , intent(inout) :: StickSlip
    real(dp), intent(in)    :: Fequ, Fext, Vel, F0, Coulomb, Stribeck, Vc

    !! Local variables
    integer  :: fricForm
    real(dp) :: VelSgn, TolMin, TolMax, Rvel, Strb

    !! --- Logic section ---

    VelSgn = sign(1.0_dp,Vel)
    TolMin = 0.1_dp*Vc
    TolMax = 2.5_dp*Vc

    call ffa_cmdlinearg_getint('fricForm',fricForm)

    !! Stick-slip or sliding and hysteresis effect

    if (Vc < epsDiv0_p) then
       StickSlip = 0 ! Stick-slip
    else if (abs(Vel) < TolMin .and. fricForm < 10) then
       StickSlip = 0 ! Stick-slip
    else if (abs(Vel) > TolMax) then
       StickSlip = int(VelSgn) ! Sliding
    end if

    if (Vel*StickSlip < 0.0_dp .and. fricForm < 10) StickSlip = 0 ! Hysteresis

    !! Maximum friction force

    if (StickSlip == 0 .and. abs(Vel) < 10.0_dp*Vc) then
       !! Stick-slip
       Rvel = Vel/Vc
       Strb = 1.0_dp + Stribeck*exp(-Rvel*Rvel)
       if (mod(fricForm,10) > 1) then
          !! Let the applied force determine the direction of the friction force
          !! when the velocity is smaller than the stick threshold
          Fmax = (F0 + Coulomb*Fequ*Strb) * sign(1.0_dp,Fext)
       else
          Fmax = (F0 + Coulomb*Fequ*Strb) * VelSgn
       end if
    else
       !! Sliding
       Fmax = (F0 + Coulomb*Fequ) * VelSgn
    end if

  end subroutine FindMaxFrictionForce


  subroutine FindFrictionForce (Ffric, Kfric, Xn_1, Vn_1, Fold, dF, &
       &                        lStick, lInit, dT, Fmax, Fext, Ftol, X0, Xn, Vn)

    !!==========================================================================
    !! Calculates the actual friction force.
    !!
    !! Programmer : Karl Erik Thoresen                         date: 11 Mar 1998
    !! Revised    : Knut Morten Okstad (F90 version)           date: 09 Jan 2001
    !! Revised    : Knut Morten Okstad (completely rewritten)  date: 22 Aug 2002
    !!==========================================================================

    use kindModule            , only : dp, epsDiv0_p
#ifdef FT_DEBUG
    use dbgUnitsModule        , only : dbgFric, dbgTime
#endif
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_intValue

    real(dp), intent(out)   :: Ffric
    integer , intent(inout) :: lStick, lInit
    real(dp), intent(inout) :: Kfric, Xn_1, Vn_1, Fold, dF
    real(dp), intent(in)    :: dT, Fmax, Fext, Ftol, X0, Xn, Vn

    !! Local variables
    real(dp) :: epsFr

    !! --- Logic section ---

    !! Calculation of friction force needed to prevent movement (stick)

    if (lInit > 0) then
#ifdef FT_DEBUG
       write(dbgFric,"(/'=== Initializing friction at time =',1pe10.3)") dbgTime
#endif

       !! Initialising at the beginning of each time step
       lInit = 0
       if (lStick > 0) lStick = 2 ! previous time step was sliding

       if (mod(ffa_cmdlinearg_intValue('fricForm'),2) == 1) then
          epsFr = 0.0001_dp/min(max(dT,0.0001_dp),1.0_dp)
       else
          epsFr = 0.001_dp
       end if
       Ffric = Fold + epsFr*Fold
#ifdef FT_DEBUG
       write(dbgFric,"('In 1: Ffric = Fold*',F6.4)") 1.0_dp+epsFr
#endif

    else if (abs(Vn-Vn_1) > epsDiv0_p) then
!!$bh: test    else if (abs(Xn-Xn_1) > epsDiv0_p) then

       !! During iterations, search for friction force that gives zero velocity.
       !! This is the condition for "stick" force.

       if (abs(dF) > epsDiv0_p) then
          Kfric = abs(dF/(Vn-Vn_1))
!!$bh: test          Kfric = dF/(Xn-Xn_1)
!!$bh: test          if (Kfric < 1.0e-5_dp) Kfric = 1.0e-5_dp
          Ffric = Fold + Kfric * Vn
!!$bh: test          Ffric = Fold + sign(Kfric*(Vn-Vn_1),Vn)
!!$bh: test          Ffric = Fold - Kfric*(Xn-X0)
#ifdef FT_DEBUG
          write(dbgFric,"('In 3: Kfric=',1pe11.3)") Kfric
#endif
       else if (abs(Fold) < epsDiv0_p .and. abs(Vn) > epsDiv0_p) then
          Ffric = 0.1_dp*Fmax
#ifdef FT_DEBUG
          write(dbgFric,"('In 4: Kfric=',1pe11.3)") Kfric
#endif
       else
          Ffric = Fold + Kfric * Vn
!!$bh: test          Ffric = Fold + sign(Kfric*(Vn-Vn_1),Vn)
!!$bh: test          Ffric = Fold - Kfric*(Xn-X0)
#ifdef FT_DEBUG
          write(dbgFric,"('In 5: Kfric=',1pe11.3)") Kfric
#endif
       end if

    else ! No position change in this iteration

       Ffric = Fold
#ifdef FT_DEBUG
       write(dbgFric,"('In 6')")
#endif
    end if

#ifdef FT_DEBUG
    write(dbgFric,"('Fmax,lStick=',1p,e11.3,i3)") Fmax,lStick
    write(dbgFric,"('Ffric,Fold =',1p2e11.3, &
         &          '  Vn,dVn =',2e11.3,'  Xn,dXn,dX =',3e11.3)") &
         &           Ffric, Fold, Vn, Vn-Vn_1, Xn, Xn-Xn_1, Xn-X0
#endif

    if (lStick == 2) then

       !! Don't stick if the velocity has not been zero
       if ( (Xn-X0)*Fold > 0.0_dp .or. &
            ffa_cmdlinearg_intValue('fricForm') >= 10) then
          Ffric = Fmax
       else
          lStick = 0 ! stick detected
       end if

    end if

    !! Stick-slip or sliding movement
    if (abs(Ffric) >= abs(Fmax)) then
       Ffric = sign(Fmax,Ffric)
       if (lStick == 0) lStick = 1 ! slip detected
    else if (abs(Fext) < Ftol .and. lStick == 0) then
       Ffric = 0.0_dp
    end if

    Xn_1 = Xn
    Vn_1 = Vn
    dF   = Ffric - Fold
    Fold = Ffric

  end subroutine FindFrictionForce


  subroutine addInFrictionForces (joint, sam, SysForce, SysReac, ierr)

    !!==========================================================================
    !! Adds friction forces for the given joint into the given system vectors.
    !!
    !! Programmer : Bjorn Haugen                     date/rev:    Nov 2001 / 1.0
    !! Revised    : Knut Morten Okstad               date/rev: 25 Jun 2002 / 2.0
    !!==========================================================================

    use MasterSlaveJointTypeModule, only : MasterSlaveJointType
    use SamModule                 , only : SamType
    use kindModule                , only : dp
    use AsmExtensionModule        , only : csAddNV
    use reportErrorModule         , only : reportError, debugFileOnly_p
#ifdef FT_DEBUG
    use fileUtilitiesModule       , only : getDBGfile
    use dbgUnitsModule            , only : dbgFric
#endif

    type(MasterSlaveJointType), intent(in)    :: joint
    type(SamType)             , intent(in)    :: sam
    real(dp)                  , intent(inout) :: SysForce(:), SysReac(:)
    integer                   , intent(inout) :: ierr

    !! Local variables
    integer  :: i, err
    logical  :: hasFriction
    real(dp) :: forceVec(joint%nJointDOFs) ! Dynamic stack allocation

    !! --- Logic section ---

#ifdef FT_DEBUG
    dbgFric = getDBGfile(9,'friction.dbg')
#endif

    forceVec = 0.0_dp
    hasFriction = .false.
    do i = 1, joint%nJointDOFs
       if (associated(joint%jointDofs(i)%friction)) then
          if (joint%jointDofs(i)%friction%isActive) then
             !! Friction routines calculate the friction force to be positive
             !! for positive velocity, but since it is resisting motion it
             !! should have a negative value here
             if ( associated(joint%jointDofs(i)%friction%spr) ) then
                !! Adding in will be handled by spring element at the joint
#ifdef FT_DEBUG
                write(dbgFric,"(' --spr-- Fric spring deflection',i3,1p2e12.3)") &
                     i, joint%jointDofs(i)%friction%spr%deflection, &
                     joint%jointDofs(i)%friction%spr%yieldDeflection
                write(dbgFric,"(' --spr-- Fric spring force',8x,1pe12.3)") joint%jointDofs(i)%friction%spr%force
#endif
             else
                forceVec(i) = -joint%jointDofs(i)%friction%force
                hasFriction = .true.
             end if
          end if
       end if
    end do

    if (hasFriction) then
       call csAddNV (sam,joint%samNodNum,forceVec,SysForce,SysReac,err)
       if (err == 0) return

       ierr = ierr - 1
       call reportError (debugFileOnly_p,'addInFrictionForces')
    end if

  end subroutine addInFrictionForces


  subroutine findJointForceInFriction (joint, sam, RF, Q, Fs, Fd, Fi)

    !!==========================================================================
    !! Computes the total applied force in a joint dof with friction.
    !! This is equal to the sum of all force contributions to that dof,
    !! except for the friction force itself.
    !!
    !! Programmer : Knut Morten Okstad               date/rev: 17 Sep 2007 / 1.0
    !!==========================================================================

    use MasterSlaveJointTypeModule, only : MasterSlaveJointType
    use SamModule                 , only : SamType
    use kindModule                , only : dp

    type(MasterSlaveJointType), intent(inout) :: joint
    type(SamType)             , intent(in)    :: sam
    real(dp)                  , intent(in)    :: RF(:), Q(:), Fs(:)
    real(dp), optional        , intent(in)    :: Fd(:), Fi(:)

    !! Local variables
    integer  :: i, ieq

    !! --- Logic section ---

    do i = 1, joint%nJointDOFs
       if (associated(joint%jointDofs(i)%friction)) then
          ieq = sam%meqn(joint%jointDofs(i)%sysDOF)
          if (ieq > 0) then
             if (present(Fd) .and. present(Fi)) then
                joint%jointDofs(i)%friction%Fext = Q(ieq) - Fs(ieq) &
                     &                           - Fd(ieq) - Fi(ieq)
             else
                joint%jointDofs(i)%friction%Fext = Q(ieq) - Fs(ieq)
             end if
          else if (associated(sam%mpreac)) then
             if (joint%jointDofs(i)%sysDOF <= size(sam%mpreac)) then
                ieq = sam%mpreac(joint%jointDofs(i)%sysDOF)
                if (ieq > 0 .and. ieq <= size(RF)) then
                   joint%jointDofs(i)%friction%Fext = RF(ieq)
                end if
             end if
          end if
       end if
    end do

  end subroutine findJointForceInFriction


  subroutine findContactElmForceInFriction (cElem, sam, RF, Q, Fs, Fd, Fi)

    !!==========================================================================
    !! Computes the total applied force in a contact element dof with friction.
    !! This is equal to the sum of all force contributions to that dof,
    !! except for the friction force itself.
    !!
    !! Programmer : Knut Morten Okstad               date/rev: 17 Sep 2007 / 1.0
    !!==========================================================================

    use ContactElementTypeModule, only : ContactElementType
    use SamModule               , only : SamType
    use kindModule              , only : dp

    type(ContactElementType), intent(inout) :: cElem
    type(SamType)           , intent(in)    :: sam
    real(dp)                , intent(in)    :: RF(:), Q(:), Fs(:)
    real(dp), optional      , intent(in)    :: Fd(:), Fi(:)

    !! TODO...
    print *,' ** findContactElmForceInFriction not yet implemented! ', &
         cElem%id%userId,sam%nel,size(RF),size(Q),size(Fs), &
         present(Fd),present(Fi)

  end subroutine findContactElmForceInFriction

end module FrictionRoutinesModule
