!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file masterSlaveJointRoutinesModule.f90
!>
!> @brief Subroutines for joint constraint calculation calculations.

!!==============================================================================
!> @brief Module with subroutines for joint constraint calculation calculations.

module MasterSlaveJointRoutinesModule

  use kindModule, only : dp

  implicit none

  real(dp), parameter :: epsZero_p = 1.0e-15_dp !< Zero tolerance

  private

  public :: moveJointInPosition, updateDependencies
  public :: incJointsVar, updateJoints, initiateJointVars


contains

  !!==========================================================================
  !> @brief Increments all joint DOFs in the model.
  !>
  !> @author Bjorn Haugen
  !> @date 16 Nov 2001
  !>
  !> @author Knut Morten Okstad
  !> @date 21 Jun 2002

  subroutine IncJointsVar (joints,inc,motions,useTotalInc)

    use MasterSlaveJointTypeModule, only : MasterSlaveJointType
    use MotionTypeModule          , only : MotionType, prescribedDeflection_p
    use RotationModule            , only : vec_to_mat

    type(MasterSlaveJointType),target, intent(inout) :: joints(:)
    real(dp)                         , intent(in)    :: inc(:)
    type(MotionType)                 , intent(in)    :: motions(:)
    logical, optional                , intent(in)    :: useTotalInc

    !! Local variables
    logical           :: fromPrevious
    integer           :: i, j, lDof, iMat, ip
    real(dp), pointer :: jVar

    !! --- Logic section ---

    if (present(useTotalInc)) then
       fromPrevious = useTotalInc
    else
       fromPrevious = .false.
    end if

    do i = 1, size(joints)

       do j = 1, joints(i)%nJointDOFs
          jVar => joints(i)%jointDofs(j)%jVar(1)
          lDof =  joints(i)%jointDofs(j)%lDof
          iMat =  joints(i)%jointDofs(j)%iMat
          ip   = -joints(i)%jointDofs(j)%loadIdx
          if (ip < 1 .or. ip > size(motions)) then
             ip = 0
          else if (motions(ip)%type /= prescribedDeflection_p) then
             ip = 0
          end if
          if (ip > 0) then
             jVar = motions(ip)%D
          else if (fromPrevious) then
             jVar = joints(i)%jointDofs(j)%jVarPrev &
                  + inc(joints(i)%jointDofs(j)%sysDof)
          else
             jVar = jVar + inc(joints(i)%jointDofs(j)%sysDof)
          end if

          if (lDof <= 3) then
             joints(i)%PosFromJVars(lDof,4,iMat) = jVar
          else if (lDof <= 6) then
             joints(i)%thetaVec(lDof-3,iMat)     = jVar
          end if
       end do

       do j = 1, joints(i)%nMats
          call vec_to_mat (joints(i)%thetaVec(:,j), &
               &           joints(i)%PosFromJVars(:,1:3,j))
       end do

    end do

  end subroutine IncJointsVar


  !!============================================================================
  !> @brief Updates joint triad positions and associated constraint variables.
  !>
  !> @author Karl Erik Thoresen
  !> @date Nov 1999
  !>
  !> @author Bjorn Haugen
  !> @date Sep 2000
  !>
  !> @author Knut Morten Okstad
  !> @date May 2008

  subroutine updateJoints (joints,motions,ierr)

    use MasterSlaveJointTypeModule, only : MasterSlaveJointType, CAM_p
    use MotionTypeModule          , only : MotionType
#ifdef FT_DEBUG
    use dbgUnitsModule            , only : dbgJoint, dbgCurve, dbgTime, dbgIter
#endif
    use MasterSlaveJointTypeModule, only : getJointId
    use ReportErrorModule         , only : reportError, error_p

    type(MasterSlaveJointType), intent(inout) :: joints(:)
    type(MotionType)          , intent(in)    :: motions(:)
    integer                   , intent(inout) :: ierr

    !! Local variables
    integer :: i, err, lerr

    !! --- Logic section ---

#ifdef FT_DEBUG
    if (dbgJoint > 0) write(dbgJoint,600) dbgTime,dbgIter
    if (dbgCurve > 0) write(dbgCurve,600) dbgTime,dbgIter
600 format(/'   --- in updateJoints at t =',1PE12.5,I4)
#endif

    lerr = ierr
    do i = 1, size(joints)

       call moveJointInPosition (joints(i),err)
       if (err < 0) ierr = ierr + err

       call updateDependencies (joints(i),motions)

       if (joints(i)%type == CAM_p .and. joints(i)%nJointDOFs == 6) then
          call SetSpringStatus (joints(i)%jointDofs,joints(i)%slider)
       end if

       if (ierr < lerr) then
          lerr = ierr
          call reportError (error_p,'Failed to update '// &
               &            getJointId(joints(i)),addString='updateJoints')
       end if

    end do

  end subroutine updateJoints


  !!============================================================================
  !> @brief Sets the status of contact springs past ends and outside the active
  !> radius region for Cam joints.
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 12 Apr 2001

  subroutine SetSpringStatus (jDOFs,slider)

    use MasterSlaveJointTypeModule, only : JointDofType, SliderType

    type(JointDofType), intent(inout) :: jDOFs(:)
    type(SliderType)  , intent(in)    :: slider

    !! Local variables
    logical  :: isActive
    integer  :: i
    real(dp) :: slideBegin, slideEnd, jVar

    !! --- Logic section ---

    !! Turn of springs if x-distance larger than contact thickness
    jVar = jDOFs(1)%jVar(1)
    if (slider%thickness > 0.0_dp .and. abs(jVar) > slider%thickness) then
       isActive = .false.

    else if (.not. slider%master%isLooping) then

       !! Turn off cam past ends of the non-looping cam
       slideBegin = slider%master%CPoints(1)%slideVar
       slideEnd   = slider%master%CPoints(size(slider%master%CPoints))%slideVar
       jVar       = slider%slideDof%jVar(1)
       isActive   = jVar >= slideBegin .and. jVar <= slideEnd

    else
       isActive = .true.
    end if

    !! Actually set the spring-, damper- and friction status for each DOF
    do i = 1, size(jDOFs)
       if (associated(jDOFs(i)%spring))   jDOFs(i)%spring%isActive   = isActive
       if (associated(jDOFs(i)%damper))   jDOFs(i)%damper%isActive   = isActive
       if (associated(jDOFs(i)%friction)) jDOFs(i)%friction%isActive = isActive
    end do

  end subroutine SetSpringStatus


  !!============================================================================
  !> @brief Updates dependent triad position and the dependent variables.
  !>
  !> @details To counteract any nonlinearities that are trying to move the
  !> triad away from its constrained position, we have to move it to the
  !> correct position.
  !>
  !> @author Karl Erik Thoresen / Bjorn Haugen
  !> @date Des 1999
  !>
  !> @author Bjorn Haugen
  !> @date Oct 2000

  subroutine moveJointInPosition (joint,ierr)

    use MasterSlaveJointTypeModule, only : MasterSlaveJointType, AXIAL_p, CAM_p
    use MasterSlaveJointTypeModule, only : GetJointId
    use CurveTypeModule           , only : GetPosOnCurve
    use ManipMatrixModule         , only : matmul34, trans1V
    use ReportErrorModule         , only : reportError, debugFileOnly_p
    use ReportErrorModule         , only : getErrorFile

    type(MasterSlaveJointType), intent(inout) :: joint
    integer                   , intent(out)   :: ierr

    !! Local Variables
    integer  :: i
    real(dp) :: Xaxis(3), P(3,4)

    !! --- Logic section ---

    ierr = 0

    if (joint%type == AXIAL_p) then

       !! The X-axis should extend from the independent joint triad
       !! to the dependent triad
       Xaxis = joint%STriad%ur(:,4) - joint%JMTriads(1)%triad%ur(:,4)
       joint%JPosInG(:,1:3) = trans1V(Xaxis,getErrorfile(),ierr)
       joint%JPosInG(:,4)   = joint%JMTriads(1)%triad%ur(:,4)
       if (ierr < 0) goto 900

    else if (size(joint%JMTriads) == 1) then

       !! Set the position of the point-to-point joint
       joint%JPosInG = matmul34(joint%JMTriads(1)%triad%ur, &
            &                   joint%JMTriads(1)%JPosInT)

       call SetLinearDependencies (joint%JMTriads(1),joint%STriad)

    else ! Set the position for a point-to-path joint

       if (joint%type == CAM_p) then
          call GetPosOnCurve (GetJointId(joint), joint%version, &
               &              joint%slider%master%CPoints, &
               &              joint%slider%slideDof%jVar(1), &
               &              joint%slider%master%loopLength, &
               &              joint%slider%master%isExtended, &
               &              joint%slider%master%isLooping, &
               &              joint%JPosInG, joint%slider%tangent, &
               &              joint%slider%coeff, joint%slider%rotGrad, ierr)
       else
          call GetPosOnCurve (GetJointId(joint), joint%version, &
               &              joint%slider%master%CPoints, &
               &              joint%slider%slideDof%jVar(1), &
               &              joint%slider%master%loopLength, &
               &              joint%slider%master%isExtended, &
               &              joint%slider%master%isLooping, &
               &              joint%JPosInG, joint%slider%tangent, &
               &              joint%slider%coeff, ierr=ierr)
       end if
       if (ierr < 0) goto 900

       do i = 1, size(joint%JMTriads)
          if (abs(joint%slider%coeff(i)) > epsZero_p) then
             call SetLinearDependencies (joint%JMTriads(i),joint%STriad)
          end if
       end do

    end if

    !! Set the position of dependent joint triad
    P = joint%STPos0InJ
    do i = 1, joint%nMats
       P = matmul34(joint%PosFromJVars(:,:,i),P)
    end do

    if (size(joint%slaveDofs) <= 3) then
       joint%STriad%ur(:,4) = matmul34(joint%JPosInG,P(:,4))
    else
       joint%STriad%ur = matmul34(joint%JPosInG,P)
    end if

    return

900 call reportError (debugFileOnly_p,'moveJointInPosition')

  contains

    !> @brief Updates the constraint equation coeffecicients of the joint.
    !> @details The linear dependency matrices in Global to Global system are:
    !> @verbatim
    !>    | Us | = | I : Spin(e) || Um | = | TT TR | | Um |
    !>    | Rs |   | 0 :   I     || Rm |   | RT RR | | Rm |
    !> @endverbatim
    !> Only the @b TR matrix needs to be established explicitly.
    !> @b TT and @b RR are always unity whereas @b RT always is zero.
    subroutine SetLinearDependencies (jointTriad,depTriad)

      use MasterSlaveJointTypeModule, only : JMTriadType, TriadType
      use MasterSlaveJointTypeModule, only : PRISMATIC_p, CYLINDRIC_p, FREE_p
      use RotationModule            , only : spin

      type(JMTriadType), intent(inout) :: jointTriad
      type(TriadType)  , intent(in)    :: depTriad

      !! Local variables
      real(dp) :: eVec(3)

      !! --- Logic section ---

      if (joint%type == PRISMATIC_p .or. joint%type == CYLINDRIC_p) then

         jointTriad%TR = 0.0_dp ! No rotation to translation coupling

      else if (joint%type == FREE_p .and. joint%version < 0) then

         jointTriad%TR = 0.0_dp ! No rotation to translation coupling

      else

         !! Dependent joint triad position relative to current independent one
         !TODO,bh: eVec should be normal to secant when point-to path joint
         !TODO,bh: eVec should be to contact point, not dependent triad ??????
         eVec = depTriad%ur(:,4) - jointTriad%triad%ur(:,4)

         !! TR = -spin(e) where e is position of dependent joint triad
         !! relative to independent triad, measured in global system
         call spin (-eVec,jointTriad%TR) ! Rotation to translation coupling

      end if

    end subroutine SetLinearDependencies

  end subroutine moveJointInPosition


  !!============================================================================
  !> @brief Updatea the linear dependencies in a joint.
  !>
  !> @details The linear dependencies are stores in the @a TTCC table in SAM.
  !> Prescribed motions in joint triad DOFs and/or directly in the joint DOFs
  !> are accounted for by multiplying the current prescribed value with the
  !> associated terms of the joint constraint equations, and adding the results
  !> to the constant term of the constraint equations (stored in @a TCC(1)
  !> for each dependent DOF).
  !>
  !> @author Karl Erik Thoresen
  !> @date 27 Sep 1998
  !>
  !> @author Knut Morten Okstad
  !> @date 12 Oct 2004

  subroutine updateDependencies (joint,motions)

    use MasterSlaveJointTypeModule, only : MasterSlaveJointType, AXIAL_p
    use MotionTypeModule          , only : MotionType
    use ManipMatrixModule         , only : unify, cross_product, matmul34
    use RotationModule            , only : dOmega_dTheta
#ifdef FT_DEBUG
    use MasterSlaveJointTypeModule, only : getJointId
    use TriadTypeModule           , only : writeObject
    use dbgUnitsModule            , only : dbgJoint
#endif

    type(MasterSlaveJointType), intent(inout) :: joint
    type(MotionType)          , intent(in)    :: motions(:)

    !! Local variables
    logical  :: isDependent
    integer  :: nSdofs, nMdofs, mDof, lm, ls, mt, ip, jp, im, iMat
    real(dp) :: Tcpl(6,6), RMS(3,3), transposeRSS(3,3), Pacc(3,4), Jmat(3,3)
    real(dp) :: dOmegaG(3), dUG(3), dUS(6), coeff
    real(dp) :: CSysNG(3,3,joint%nMats) ! Dynamic scratch allocation
    real(dp) :: PNG(3,4,joint%nMats)    ! Dynamic scratch allocation
    real(dp), pointer :: tcc(:), ctcc(:)

    !! --- Logic section ---

    nSdofs = size(joint%slaveDofs) ! Number of dependent DOFs (3 or 6)

    !! Initialize the coefficients of the constraint equations to zero
    do ls = 1, nSdofs
       joint%slaveDofs(ls)%tcc = 0.0_dp
    end do

    !! System directions for dependent triad
    if (associated(joint%STriad%SysDirInG)) then
       !! Calculate the transpose since that is used mostly
       transposeRSS = transpose(joint%STriad%SysDirInG)
    else
       call unify (transposeRSS) ! The dependent system directions are global
    end if
    Tcpl = 0.0_dp

    !! === Connection to the independent joint triad DOFs ===

    do mt = 1, size(joint%JMTriads)

       nMdofs = joint%JMTriads(mt)%nDOFs
       if (nMdofs < 1) cycle

       if (associated(joint%slider)) then
          coeff = joint%slider%coeff(mt)
          if (abs(coeff) <= epsZero_p) cycle
       else
          coeff = 1.0_dp
       end if

       Tcpl(1:3,1:3) = transposeRSS*coeff ! TT is unity
       if (joint%type /= AXIAL_p) then ! No rotational coupling for AXIAL joints
          Tcpl(1:3,4:6) = matmul(transposeRSS,joint%JMTriads(mt)%TR)*coeff
          Tcpl(4:6,4:6) = transposeRSS*coeff ! RR is unity
       end if

       if (associated(joint%JMTriads(mt)%triad%SysDirInG)) then
          !! Joint triad system directions
          RMS = joint%JMTriads(mt)%triad%SysDirInG
          Tcpl(:,1:3) = matmul(Tcpl(:,1:3),RMS)
          Tcpl(:,4:6) = matmul(Tcpl(:,4:6),RMS)
       end if

       !! Check if this triad also is dependent in another joint (chaining)
       if (associated(joint%chain)) then
          isDependent = associated(joint%chain%STriad,joint%JMTriads(mt)%triad)
       else
          isDependent = .false.
       end if

       !! Insert the contributions from current triad
       !! into the constraint equations associated with each dependent DOF
       ip = joint%JMTriads(mt)%ipMceq
       do lm = 1, joint%JMTriads(mt)%nDOFs
          if (isDependent) then
             !! This triad is a dependent joint triad in joint%chain,
             !! use the constraint equation of that joint (assuming it has
             !! already been processed) to find the resulting coefficients
             !! for this joint
             ctcc => joint%chain%slaveDofs(lm)%tcc
             jp = ip + size(ctcc)-2
             do ls = 1, nSdofs
                tcc => joint%slaveDofs(ls)%tcc
                tcc(1)     = tcc(1)     + Tcpl(ls,lm)*ctcc(1)
                tcc(ip:jp) = tcc(ip:jp) + Tcpl(ls,lm)*ctcc(2:)
             end do
          else
             do ls = 1, nSdofs
                joint%slaveDofs(ls)%tcc(ip) = Tcpl(ls,lm)
             end do
             ip = ip + 1
             im = joint%JMTriads(mt)%triad%eqNumber(lm)
             if (im < 0) then ! Prescribed joint triad DOF
                !! Update the constant coefficient of the constraint equations
                !! with contributions from the current prescribed motion
                coeff = motions(-im)%cc
                do ls = 1, nSdofs
                   joint%slaveDofs(ls)%tcc(1) = joint%slaveDofs(ls)%tcc(1) &
                        &                     + Tcpl(ls,lm)*coeff
                end do
             end if
          end if
       end do

    end do


    !! === Connection to the hierarchical joint variables ===

    !! Establish accumulative positions in Global System
    Pacc = joint%JPosInG
    do iMat = joint%nMats, 1, -1
       CSysNG(:,:,iMat) = Pacc(:,1:3) ! Coordinate system of dofs in iMat
       Pacc = matmul34(Pacc,joint%PosFromJVars(:,:,iMat))
       PNG(:,:,iMat) = Pacc           ! Position and orientation of body iMat
    end do

    do mDof = 1, size(joint%jointDofs)
       iMat = joint%jointDofs(mDof)%iMat
       lm   = joint%jointDofs(mDof)%ldof
       ip   = joint%jointDofs(mDof)%ipMceq
       im   = joint%jointDofs(mDof)%loadIdx

       if (lm <= 3) then ! Translational independent DOF

          dUS(1:3) = matmul(transposeRSS,CSysNG(:,lm,iMat)) ! Line of RSS'*RMS
          dUS(4:6) = 0.0_dp ! No coupling to the rotational dependent DOFs

       else if (lm <= 6) then ! Rotational independent DOF

          if (joint%numRot(iMat) > 1) then
             !! Establish dOmega/dTheta
             call dOmega_dTheta (joint%thetaVec(:,iMat),Jmat)
             dOmegaG = matmul(CSysNG(:,:,iMat),Jmat(:,lm-3))
          else
             !! When only one rotational DOF the matrix dOmega/dTheta simplifies
             dOmegaG = CSysNG(:,lm-3,iMat)
          end if

          dUG = joint%STriad%ur(:,4) - PNG(:,4,iMat) ! Eccentricity to dependent
          !                                          ! triad in global system
          dUG = cross_product(dOmegaG,dUG) ! Displacement gradient of dependent
          !                                ! triad in the global system, and
          dUS(1:3) = matmul(transposeRSS,dUG) ! in the dependent triad system
          dUS(4:6) = matmul(transposeRSS,dOmegaG) ! Rotational dependent DOFs

       else if (lm == 7) then ! Slider DOF

          dUS(1:3) = matmul(transposeRSS,joint%slider%tangent)
          dUS(4:6) = matmul(transposeRSS,joint%slider%rotGrad)

       else
          cycle
       end if

       !! Insert the contributions from current joint DOF
       !! into the constraint equations associated with each dependent DOF
       coeff = joint%jointDofs(mDof)%coeff
       do ls = 1, nSdofs
          joint%slaveDofs(ls)%tcc(ip) = dUS(ls)*coeff
       end do
       if (im < 0) then ! Prescribed joint DOF
          !! Update the constant coefficient of the constraint equations
          !! with contributions from the current prescribed motion
          coeff = joint%jointDofs(mDof)%coeff * motions(-im)%cc
          do ls = 1, nSdofs
             joint%slaveDofs(ls)%tcc(1) = joint%slaveDofs(ls)%tcc(1) &
                  &                     + dUS(ls)*coeff
          end do
       end if

    end do

#ifdef FT_DEBUG
    if (dbgJoint == 0) return
    write(dbgJoint,600) trim(getJointId(joint))
    do ls = 1, nSdofs
       ip = 1
       jp = min(13,size(joint%slaveDofs(ls)%mceq))
       write(dbgJoint,610) ls,joint%slaveDofs(ls)%mceq(ip:jp)
       write(dbgJoint,620)    joint%slaveDofs(ls)%tcc(ip:jp)
       do while (jp < size(joint%slaveDofs(ls)%mceq))
          ip = jp + 1
          jp = min(ip+11,size(joint%slaveDofs(ls)%mceq))
          write(dbgJoint,611) joint%slaveDofs(ls)%mceq(ip:jp)
          write(dbgJoint,621) joint%slaveDofs(ls)%tcc(ip:jp)
       end do
    end do
    write(dbgJoint,"(/'Dependent Triad ')",advance='NO')
    call writeObject (joint%STriad,dbgJoint)
600 format(/'Constraint equations associated with ',A)
610 format(5X,'Dependent Triad',I2,':',I8,12I12)
611 format(31X,12I12)
620 format(23X,1P13E12.4)
621 format(35X,1P12E12.4)
#endif

  end subroutine updateDependencies


  !!============================================================================
  !> @brief Initializes the joint variables (velocities and accelerations).
  !>
  !> @details The initialization is based on the current state of the
  !> joint triads, or the dependent DOF vel/acc
  !> if the joint variable vel/acc are given instead.
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 24 Oct 2002

  subroutine initiateJointVars (joint,ierr)

    use MasterSlaveJointTypeModule, only : MasterSlaveJointType
    use MasterSlaveJointTypeModule, only : hasZeroVelAcc, getJointId, getId
    use ScratchArrayModule        , only : getRealScratchMatrix
    use ReportErrorModule         , only : internalError, getErrorFile
    use ReportErrorModule         , only : reportError, error_p, warning_p
    use ReportErrorModule         , only : debugFileOnly_p, warningFileOnly_p

    type(MasterSlaveJointType), intent(inout) :: joint
    integer                   , intent(out)   :: ierr

    !! Local variables
    logical             :: hasInitVelAcc
    integer             :: i, j, k, nJ, nS, lpu
    real(dp)            :: RHS(6,2), Work(12)
    real(dp), pointer   :: Tcc(:,:)
    real(dp), parameter :: zero_p = 10.0_dp*tiny(1.0_dp)
    character(len=64)   :: errMsg

    !! --- Logic section ---

    nS = size(joint%slaveDofs)
    nJ = joint%nJointDOFs

    if (nS < 1) then

       !! The dependent triad is grounded (should normally not happen)
       ierr = internalError('initiateJointVars: Grounded dependent triad')
       return

    else if (hasZeroVelAcc(joint)) then

       !! In many cases we start from zero velocities/accelerations in all DOFs
       !! and we can then skip further computations to save initialization time
       ierr = 0
       return

    else if (associated(joint%chain)) then

       ierr = -1
       call reportError (error_p,'Initial conditions in chained joints '// &
            &            'are currently not supported.', &
            &            'Detected for '//getJointId(joint), &
            &            addString='initiateJointVars')
       return

    end if

    !! Get the matrix of constraint coefficients for this joint
    Tcc => getRealScratchMatrix(6,size(joint%slaveDofs(1)%tcc),ierr)
    if (ierr < 0) goto 100

    do i = 1, nS
       Tcc(i,:) = joint%slaveDofs(i)%tcc
    end do

    !! Compute contributions to the dependent DOF velocities/accelerations
    !! from initial velocities/accelerations at the joint triad DOFs
    RHS(:,1) = Tcc(:,1)
    RHS(:,2) = Tcc(:,1)
    do i = 1, size(joint%JMTriads)
       j = joint%JMTriads(i)%ipMceq
       do k = 1, joint%JMTriads(i)%nDOFs
          if (joint%JMTriads(i)%triad%eqNumber(k) /= 0) then
             RHS(:,1) = RHS(:,1) + Tcc(:,j)*joint%JMTriads(i)%triad%urd(k)
             RHS(:,2) = RHS(:,2) + Tcc(:,j)*joint%JMTriads(i)%triad%urdd(k)
          end if
          j = j + 1
       end do
    end do

    !! Check if initial joint velocities and/or accelerations are given.
    !! They will override initial conditions specified for the dependent triad.
    hasInitVelAcc = .false.
    do i = 1, nJ
       if ( abs(joint%jointDofs(i)%jVar(2)) > zero_p .or. &
            abs(joint%jointDofs(i)%jVar(3)) > zero_p ) hasInitVelAcc = .true.
    end do

    if (hasInitVelAcc .or. nJ < 1) then

       !! Compute dependent DOF values of initial velocities and accelerations
       !! based on the joint triad DOF values, and the joint variable values
       do i = 1, nJ
          j = joint%jointDofs(i)%ipMceq
          if ( joint%slaveDofs(1)%mceq(j) > 0 .or. & ! free DOF
               joint%jointDofs(i)%loadIdx < 0 ) then ! prescribed DOF
             RHS(:,1) = RHS(:,1) + Tcc(:,j)*joint%jointDofs(i)%jVar(2)
             RHS(:,2) = RHS(:,2) + Tcc(:,j)*joint%jointDofs(i)%jVar(3)
          end if
       end do

       lpu = getErrorFile()
       write(lpu,600) 'dependent Triad', trim(getId(joint%STriad%id))
       write(lpu,610) 'Velocity     :', RHS(1:nS,1)
       write(lpu,610) 'Acceleration :', RHS(1:nS,2)

       if (dot_product(joint%STriad%urd,joint%STriad%urd) > zero_p) then
          call warningOverride (joint%STriad%id,'initial velocities')
       end if
       if (dot_product(joint%STriad%urdd,joint%STriad%urdd) > zero_p) then
          call warningOverride (joint%STriad%id,'initial accelerations')
       end if

       joint%STriad%urd(1:nS)  = RHS(1:nS,1)
       joint%STriad%urdd(1:nS) = RHS(1:nS,2)

    else

       !! There are no initial velocities or accelerations at the joint DOFs,
       !! but some non-zero initial values are given at the joint triad DOFs.
       !! Try to derive the equivalent initial joint variables by
       !! inverting the constraint equations of the joint.

       !! Check for suppressed joint DOFs and squeeze out the respective columns
       !! in the Tcc array such that it only represents the free joint DOFs
       nJ = 0
       do i = 1, joint%nJointDOFs
          j = joint%jointDofs(i)%ipMceq
          if (joint%slaveDofs(1)%mceq(j) > 0) then
             !! This joint DOF is free
             nJ = nJ + 1
             Tcc(:,nJ) = Tcc(:,j)
          end if
       end do
       if (nJ > 0) then

          !! Establish the right-hand-sides of the system of equations
          !! in which the joint DOFs are unknowns instead of the dependent DOFs
          RHS(1:nS,1) = joint%STriad%urd(1:nS)  - RHS(1:nS,1)
          RHS(1:nS,2) = joint%STriad%urdd(1:nS) - RHS(1:nS,2)

          !! Skip solution if only zero RHS (the solution is then zero as well)
          if ( dot_product(RHS(1:nS,1),RHS(1:nS,1)) <= zero_p .and. &
               dot_product(RHS(1:nS,2),RHS(1:nS,2)) <= zero_p ) nJ = 0

       end if
       if (nJ < 1) return ! Skip solution if only prescribed DOFs

       !! Solve the possibly over-determined system using least squares
       lpu = getErrorFile()
       call DGELS ('N',nS,nJ,2,Tcc(1,1),6,RHS,6,Work(1),size(Work),ierr)
       if (ierr /= 0) then
          write(errMsg,"('Error return from LAPACK::DGELS, INFO =',I6)") ierr
          call reportError (error_p,'Failed to solve for initial velocities'// &
               ' and accelerations for '//getJointId(joint),addString=errMsg)
          if (ierr == -10) write(lpu,"(10X,'nS,nJ,nW =',3I6)") nS,nJ,size(Work)
          goto 100
       end if

       !! Assign initial vel/acc conditions to the free joint variables
       nJ = 0
       do i = 1, joint%nJointDOFs
          j = joint%jointDofs(i)%ipMceq
          if (joint%slaveDofs(1)%mceq(j) > 0) then
             nJ = nJ + 1
             joint%jointDofs(i)%jVar(2:3) = RHS(nJ,:)
          end if
       end do

       write(lpu,600) trim(getJointId(joint))
       write(lpu,610) 'based on the inverse constraint equation '// &
            &         'and the triad DOFs.'
       write(lpu,610) 'Velocity     :', RHS(1:nJ,1)
       write(lpu,610) 'Acceleration :', RHS(1:nJ,2)

    end if
    write(lpu,"('')")

    return

100 call reportError (debugFileOnly_p,'initiateJointVars')
600 format(/5X,'>>> Computed initial conditions for ',2A )
610 format( 9X,A,1P10E13.5 )

  contains

    !> @brief Convenience routine for printing a warning message.
    subroutine warningOverride (triadId,var)
      use MasterSlaveJointTypeModule, only : IdType, getId, getJointId, RIGID_p
      type(idType)    , intent(in) :: triadId
      character(len=*), intent(in) :: var
      if (joint%type == RIGID_p) then
         call reportError (warningFileOnly_p,'The non-zero '//trim(var)// &
              ' specified for dependent Triad'//getId(triadId), &
              'of '//trim(getJointId(joint))// &
              ' are overridden by the independent joint triad DOFs')
      else
         call reportError (warning_p,'The non-zero '//trim(var)// &
              ' specified for dependent Triad'//getId(triadId), &
              'of '//trim(getJointId(joint))// &
              ' are overridden by the joint- and independent triad DOFs')
      end if
    end subroutine warningOverride

  end subroutine initiateJointVars

end module MasterSlaveJointRoutinesModule
