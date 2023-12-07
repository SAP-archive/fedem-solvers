!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module MasterSlaveJointTypeModule

  use IdTypeModule        , only : IdType, getId
  use TriadTypeModule     , only : TriadType, dp
  use SpringTypeModule    , only : SpringBaseType, SpringType
  use DamperTypeModule    , only : DamperBaseType
  use FrictionTypeModule  , only : FrictionType
  use ContactSurfaceModule, only : GliderCurveType

  implicit none

  integer, parameter :: FOLLOWER_AXIS_p   = 1, &
       &                ORTHOGONAL_AXIS_p = 2, &
                        ROT_AXIS_p        = 3

  character(len=15), parameter :: rotParamTypes_p(3) = (/ 'FOLLOWER_AXIS  ', &
       &                                                  'ORTHOGONAL_AXIS', &
       &                                                  'ROT_AXIS       ' /)

  integer, parameter :: REVOLUTE_p  = 1, &
       &                BALL_p      = 2, &
       &                RIGID_p     = 3, &
       &                FREE_p      = 4, &
       &                PRISMATIC_p = 5, &
       &                CYLINDRIC_p = 6, &
       &                CAM_p       = 7, &
       &                AXIAL_p     = 8

  character(len=15), parameter :: jointTypeName_p(8) = (/ 'Revolute joint ', &
       &                                                  'Ball joint     ', &
       &                                                  'Rigid joint    ', &
       &                                                  'Free joint     ', &
       &                                                  'Prismatic joint', &
       &                                                  'Cylindric joint', &
       &                                                  'Cam joint      ', &
       &                                                  'Axial joint    ' /)


  type JMTriadType
     type(TriadType), pointer :: triad        ! This is the master triad
     real(dp)       , pointer :: JPosInT(:,:) ! Relative joint position
     integer                  :: nDOFs        ! <= 3 for AXIAL_p, else <= 6
     integer                  :: ipMceq       ! Index in mceq for master
     real(dp)                 :: TR(3,3)      ! Translation-Rotation coupling
  end type JMTriadType


  type SlaveDofType
     real(dp), pointer :: tcc(:)  ! Table of Constraint equation Coefficients
     integer , pointer :: mceq(:) ! mceq(1) is global DOF of current slave DOF
  end type SlaveDofType


  type JointDofType
     integer  :: iMat       !! Which position matrix this DOF is associated with
     integer  :: lDOF       !! Which DOF (1-6) within the position matrix
     integer  :: sysDOF     !! System DOF number in SAM
     integer  :: ipMceq     !! Position of this master DOF in mceq
     integer  :: loadIdx    !! Index to load/prescribed motion, if any
     logical  :: fixed      !! Is this dof fixed? Needed for save/restart
     real(dp) :: coeff      !! = 1.0, unless master is from Higher Pairs
     real(dp) :: jVar(3)    !! Pos., vel., and accel. of the joint variable
     real(dp) :: jVarPrev   !! Joint variable value at previous time step
     real(dp), pointer :: F !! Pointer to reaction force or applied force
     logical  :: saveVar(4) !! Flags indicating which variables should be saved

     type(SpringBaseType), pointer :: spring
     type(DamperBaseType), pointer :: damper
     type(FrictionType)  , pointer :: friction
  end type JointDofType


  type SliderType
     type(GliderCurveType), pointer :: master   !! Glider curve definition
     type(JointDofType)   , pointer :: slideDOF !! Slide DOF along master curve
     real(dp)             , pointer :: coeff(:) !! Scaling of each master triad
     real(dp) :: tangent(3) !! Unit tangent vector in slide direction
     real(dp) :: rotGrad(3) !! Rate of rotation w.r.t. slide variable
     real(dp) :: thickness  !! Thickness of cam domain in local x-direction
  end type SliderType


  type MasterSlaveJointType

     type(IdType) :: id       !! General identification data

     integer  :: type         !! The joint type
     integer  :: version      !! Joint type version
     integer  :: rotParam     !! Rotational parametrization type
     integer  :: samNodNum    !! Node number for SAM reference (madof)
     integer  :: nJointDOFs   !! Number of degrees of freedom in the joint

     real(dp) :: JPosInG(3,4) !! Position of joint relative to the global system

     integer, pointer :: BC(:) ! 0 is fixed, 1 is free, 2 is fixed during
     !                         ! initial static equilibrium iterations and
     !                         ! eigenvalue calculations, and free otherwise

     type(TriadType), pointer :: STriad         !! The slave triad
     real(dp)                 :: STPos0InJ(3,4) !! Position of slave triad
     !                                          !! in joint system when all
     !                                          !! joint variables are zero

     type(JMTriadType) , pointer :: JMTriads(:)

     type(SlaveDofType), pointer :: slaveDOFs(:)

     type(JointDofType), pointer :: jointDOFs(:)

     type(SliderType)  , pointer :: slider !! Data for multi-master joints only

     type(MasterSlaveJointType), pointer :: chain !! If one of the master triads
     !                                            !! also is a slave this points
     !                                            !! to the other joint where
     !                                            !! that triad is slave

     type(SpringType), pointer :: springEl !< For interconnected springs
     type(SpringType), pointer :: sprFric  !< Friction spring element

     !!TODO,bh: This should be removed or recoded (unused code/work in progress)
     type(FrictionType), pointer :: multiDofFriction !! For ball joints, etc.

     integer :: nMats                         !! Number of intermediate matrices
     real(dp), pointer :: PosFromJVars(:,:,:) !! dim(3,4,nMats)
     real(dp), pointer :: thetaVec(:,:)       !! dim(3,nMats)
     integer , pointer :: numRot(:)           !! dim(nMats)

     integer , pointer :: dofOutputOrder(:)   !! DOF-order on the frs-file

  end type MasterSlaveJointType


  type HigherPairType
     type(IdType)                        :: id
     type(MasterSlaveJointType), pointer :: slaveJoint   , masterJoint
     integer                             :: slaveJointDOF, masterJointDOF
     real(dp)                            :: coeff ! slave = coeff * master
  end type HigherPairType


  interface GetPtrToId
     module procedure GetPtrToIdJoint
  end interface

  interface GetPtrToOwner
     module procedure GetPtrToJointWithSpring
     module procedure GetPtrToJointWithDamper
  end interface

  interface WriteObject
     module procedure WriteJoint
  end interface

  interface DeallocateJoints
     module procedure DeallocateMasterSlaveJoints
     module procedure DeallocateHigherPairs
  end interface

  interface updateAtConvergence
     module procedure updatePreviousJointValues
  end interface

  interface restoreFromLastStep
     module procedure restorePreviousJointValues
  end interface

  private :: GetPtrToIdJoint, GetPtrToJointWithSpring, GetPtrToJointWithDamper
  private :: WriteJoint, updatePreviousJointValues, restorePreviousJointValues


contains

  function GetPtrToIdJoint (array,id,index,jointType) result(ptr)

    !!==========================================================================
    !! Return pointer to (first) element with specified id or NULL if not found.
    !! Optionally, the userID for a specified jointType is searched for
    !! (default is to search for baseID and any joint type).
    !!
    !! Programmer : Bjorn Haugen / Knut Morten Okstad
    !! date/rev   : 17 Dec 2009 / 1.0
    !!==========================================================================

    type(MasterSlaveJointType), pointer               :: ptr
    type(MasterSlaveJointType), intent(in) , target   :: array(:)
    integer                   , intent(in)            :: id
    integer                   , intent(out), optional :: index
    integer                   , intent(in) , optional :: jointType

    !! Local variables
    integer :: i

    !! --- Logic section ---

    do i = 1, size(array)
      if (present(jointType)) then
         if (array(i)%id%userId /= id .or. array(i)%type /= jointType) cycle
      else
         if (array(i)%id%baseId /= id) cycle
      end if
      ptr => array(i)
      if (present(index)) index = i
      return
    end do

    nullify(ptr)
    if (present(jointType)) then
       write(*,*) '*** GetPtrToIdJoint returned nullified, userId =',id
    else
       write(*,*) '*** GetPtrToIdJoint returned nullified, baseId =',id
    end if

  end function GetPtrToIdJoint


  function GetPtrToJointWithSpring (array,spring) result(ptr)

    !!==========================================================================
    !! Return pointer to joint connected to given spring, or NULL if not found.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 1 Jul 2002 / 1.0
    !!==========================================================================

    type(MasterSlaveJointType), pointer            :: ptr
    type(MasterSlaveJointType), intent(in), target :: array(:)
    type(SpringBaseType)      , intent(in), target :: spring

    !! Local variables
    integer :: i, j

    !! --- Logic section ---

    do i = 1, size(array)
       do j = 1, size(array(i)%jointDOFs)
          if (associated(array(i)%jointDOFs(j)%spring,spring)) then
             ptr => array(i)
             return
          end if
       end do
    end do

    nullify(ptr)
    write(*,*) '*** GetPtrToJointWithSpring returned nullified, springId ='// &
         &     getId(spring%id)

  end function GetPtrToJointWithSpring


  function GetPtrToJointWithDamper (array,damper) result(ptr)

    !!==========================================================================
    !! Return pointer to joint connected to given damper, or NULL if not found.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 1 Jul 2002 / 1.0
    !!==========================================================================

    type(MasterSlaveJointType), pointer            :: ptr
    type(MasterSlaveJointType), intent(in), target :: array(:)
    type(DamperBaseType)      , intent(in), target :: damper

    !! Local variables
    integer :: i, j

    !! --- Logic section ---

    do i = 1, size(array)
       do j = 1, size(array(i)%jointDOFs)
          if (associated(array(i)%jointDOFs(j)%damper,damper)) then
             ptr => array(i)
             return
          end if
       end do
    end do

    nullify(ptr)
    write(*,*) '*** GetPtrToJointWithDamper returned nullified, damperId ='// &
         &     getId(damper%id)

  end function GetPtrToJointWithDamper


  subroutine WriteJoint (joint,io,complexity)

    !!==========================================================================
    !! Standard routine for writing an object to io.
    !!
    !! Programmer : Karl Erik Thoresen
    !! date/rev   : 27 Sep 1998 / 1.0
    !!==========================================================================

    use IdTypeModule, only : writeId

    type(MasterSlaveJointType), intent(in) :: joint
    integer                   , intent(in) :: io
    integer, optional         , intent(in) :: complexity

    !! Local variables
    integer :: i, j

    !! --- Logic section ---

    write(io,'(A)') 'Joint','{'
    call writeId (joint%id,io)

    if (joint%type > 0 .and. joint%type <= size(jointTypeName_p)) then
       write(io,*) 'type        = ',jointTypeName_p(joint%type)
    else
       write(io,*) 'type        =', joint%type
    end if
    write(io,*) 'version     =', joint%version
    if (joint%rotParam > 0 .and. joint%rotParam <= size(rotParamTypes_p)) then
       write(io,*) 'rotParam    = ',rotParamTypes_p(joint%rotParam)
    else
       write(io,*) 'rotParam    =', joint%rotParam
    end if
    write(io,*) 'samNodNum   =', joint%samNodNum
    write(io,*) 'nJointDOFs  =', joint%nJointDOFs
    if (associated(joint%BC)) then
       write(io,*) 'BC          =', joint%BC
    end if

    if (associated(joint%chain)) then
       write(io,*) 'chain joint =', joint%chain%id%baseId
    end if

    if (associated(joint%STriad)) then
       write(io,*) 'slave(id)   =', joint%STriad%id%baseId
    end if
    if (associated(joint%slaveDofs)) then
       write(io,*) 'nSlaveDOFs  =', size(joint%slaveDofs)
    end if
    if (associated(joint%JMTriads)) then
       write(io,*) 'masters(id) =',(joint%JMTriads(i)%triad%id%baseId, &
            &                       i=1,size(joint%JMTriads))
       write(io,*) 'nMasterDOFs =',(joint%JMTriads(i)%nDOFs, &
            &                       i=1,size(joint%JMTriads))
    end if

    if (present(complexity)) then
       if (complexity >= 1 .and. associated(joint%JMTriads)) then
          write(io,601) ' STPos0InJ   =',(joint%STPos0InJ(i,:),i=1,3)
          do j = 1, size(joint%JMTriads)
             write(io,602) j,(joint%JMTriads(j)%JPosInT(i,:),i=1,3)
          end do
       end if
       if (complexity >= 2) then
          write(io,601) ' JPosInG     =',(joint%JPosInG(i,:),i=1,3)
          do i = 1, size(joint%jointDOFs)
             write(io,603) i,joint%jointDOFs(i)%lDOF, &
                  &          joint%jointDOFs(i)%sysDOF, &
                  &          joint%jointDOFs(i)%jVar
          end do
       end if
    end if

    write(io,'(A)') '}'

601 format(A14,1P,4E13.5/(14X,4E13.5))
602 format(' JPosInT',I3,'  =',1P,4E13.5/(14X,4E13.5))
603 format(' Joint var',I2,' =',I3,I10,1P,3E13.5)

  end subroutine WriteJoint


  subroutine NullifyJoint (joint)

    !!==========================================================================
    !! Initialize the MasterSlaveJointType object.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 24 Jun 2002 / 1.0
    !!==========================================================================

    use IdTypeModule, only : nullifyId

    type(MasterSlaveJointType), intent(out) :: joint

    !! --- Logic section ---

    call nullifyId (joint%id)

    joint%type = 0
    joint%version = 0
    joint%rotParam = 0
    joint%samNodNum = 0
    joint%nJointDOFs = 0
    joint%JPosInG = 0.0_dp

    nullify(joint%STriad)
    joint%STPos0InJ = 0.0_dp

    nullify(joint%BC)
    nullify(joint%JMTriads)
    nullify(joint%slaveDOFs)
    nullify(joint%jointDOFs)
    nullify(joint%slider)
    nullify(joint%chain)
    nullify(joint%springEl)
    nullify(joint%sprFric)
    nullify(joint%multiDofFriction)

    joint%nMats = 0
    nullify(joint%PosFromJVars)
    nullify(joint%thetaVec)
    nullify(joint%numRot)

    nullify(joint%dofOutputOrder)

  end subroutine NullifyJoint


  subroutine NullifyJointDof (jointDof)

    !!==========================================================================
    !! Initialize the JointDofType object.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 4 Oct 2005 / 1.0
    !!==========================================================================

    type(JointDofType), intent(out) :: jointDof

    !! --- Logic section ---

    jointDof%iMat     = 0
    jointDof%lDof     = 0
    jointDof%sysDof   = 0
    jointDof%ipMceq   = 0
    jointDof%loadIdx  = 0
    jointDof%fixed    = .false.
    jointDof%coeff    = 1.0_dp
    jointDof%jVar     = 0.0_dp
    jointDof%jVarPrev = 0.0_dp
    jointDof%saveVar  = .false.

    nullify(jointDof%F)
    nullify(jointDof%spring)
    nullify(jointDof%damper)
    nullify(jointDof%friction)

  end subroutine NullifyJointDof


  subroutine DeallocateJoint (joint)

    !!==========================================================================
    !! Deallocate the MasterSlaveJointType object.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 23 Jan 2017 / 1.0
    !!==========================================================================

    use IdTypeModule    , only : deallocateId
    use SpringTypeModule, only : deallocateSpring

    type(MasterSlaveJointType), intent(inout) :: joint

    !! Local variables
    integer :: i

    !! --- Logic section ---

    call deallocateId (joint%id)

    if (associated(joint%BC)) deallocate(joint%BC)
    if (associated(joint%JMTriads)) then
       if ( size(joint%JMTriads) == 1 .and. &
            associated(joint%JMTriads(1)%triad%scoord) ) then
          nullify(joint%STriad%scoord) ! Points to the master triad's scoord
       end if
       if (.not. associated(joint%slider)) then
          do i = 1, size(joint%JMTriads)
             deallocate(joint%JMTriads(i)%JPosInT)
          end do
       end if
       deallocate(joint%JMTriads)
    end if
    if (associated(joint%slaveDOFs)) deallocate(joint%slaveDOFs)
    if (associated(joint%jointDOFs)) then
       do i = 1, size(joint%jointDOFs)
          if (associated(joint%jointDOFs(i)%friction)) then
             deallocate(joint%jointDOFs(i)%friction)
          end if
       end do
       deallocate(joint%jointDOFs)
    end if
    if (associated(joint%slider)) then
       if (associated(joint%slider%coeff)) deallocate(joint%slider%coeff)
       deallocate(joint%slider)
    end if
    if (associated(joint%springEl)) then
       nullify(joint%springEl%id)
       call deallocateSpring (joint%springEl)
       deallocate(joint%springEl)
    end if
    if (associated(joint%sprFric)) then
       nullify(joint%sprFric%id)
       call deallocateSpring (joint%sprFric)
       deallocate(joint%sprFric)
    end if
    if (associated(joint%multiDofFriction)) deallocate(joint%multiDofFriction)
    if (associated(joint%PosFromJVars))     deallocate(joint%PosFromJVars)
    if (associated(joint%thetaVec))         deallocate(joint%thetaVec)
    if (associated(joint%numRot))           deallocate(joint%numRot)
    if (associated(joint%dofOutputOrder))   deallocate(joint%dofOutputOrder)

    call nullifyJoint (joint)

  end subroutine DeallocateJoint


  subroutine deallocateMasterSlaveJoints (joints)
    type(MasterSlaveJointType), pointer :: joints(:)
    integer :: i
    do i = 1, size(joints)
       call DeallocateJoint (joints(i))
    end do
    deallocate(joints)
    nullify(joints)
  end subroutine deallocateMasterSlaveJoints


  subroutine deallocateHigherPairs (higherPairs)
    use IdTypeModule, only : deallocateId
    type(HigherPairType), pointer :: higherPairs(:)
    integer :: i
    do i = 1, size(higherPairs)
       call DeallocateId (higherPairs(i)%id)
    end do
    deallocate(higherPairs)
    nullify(higherPairs)
  end subroutine deallocateHigherPairs


  subroutine SetJointsVelAcc (joints,velGlobal,accGlobal)

    !!==========================================================================
    !! Set velocity and acceleration for all joint DOFs from the system vectors.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 20 Jun 2002 / 1.0
    !!==========================================================================

    type(MasterSlaveJointType), intent(inout) :: joints(:)
    real(dp)                  , intent(in)    :: velGlobal(:), accGlobal(:)

    !! Local variables
    integer :: i, j, sysDOF

    !! --- Logic section ---

    do i = 1, size(joints)
       do j = 1, joints(i)%nJointDOFs
          sysDOF = joints(i)%jointDOFs(j)%sysDOF
          joints(i)%jointDOFs(j)%jVar(2) = velGlobal(sysDOF)
          joints(i)%jointDOFs(j)%jVar(3) = accGlobal(sysDOF)
       end do
    end do

  end subroutine SetJointsVelAcc


  subroutine GetJointsVelAcc (joints,velGlobal,accGlobal)

    !!==========================================================================
    !! Fill system velocity- and acceleration vectors with joint DOF values.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 18 Oct 2002 / 1.0
    !!==========================================================================

    type(MasterSlaveJointType), intent(in)  :: joints(:)
    real(dp)                  , intent(out) :: velGlobal(:), accGlobal(:)

    !! Local variables
    integer :: i, j, sysDOF

    !! --- Logic section ---

    do i = 1, size(joints)
       do j = 1, joints(i)%nJointDOFs
          sysDOF = joints(i)%jointDOFs(j)%sysDOF
          velGlobal(sysDOF) = joints(i)%jointDOFs(j)%jVar(2)
          accGlobal(sysDOF) = joints(i)%jointDOFs(j)%jVar(3)
       end do
    end do

  end subroutine GetJointsVelAcc


  function GetJointVar (joint,dofInd,type)

    !!==========================================================================
    !! Return the current value of a joint variable.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 8 Apr 2005 / 1.0
    !!==========================================================================

    type(MasterSlaveJointType), intent(in) :: joint
    integer                   , intent(in) :: dofInd, type
    real(dp)                               :: GetJointVar

    !! --- Logic section ---

    if (dofInd > 0 .and. dofInd <= size(joint%jointDOFs)) then
       if (type > 0 .and. type <= size(joint%jointDOFs(dofInd)%jVar)) then
          GetJointVar = joint%jointDOFs(dofInd)%jVar(type)
       else
          GetJointVar = 0.0_dp
       end if
    else
       GetJointVar = 0.0_dp
    end if

  end function GetJointVar


  function GetJointId (joint)

    !!==========================================================================
    !! Return the full id (type name, user id and description) of a joint.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 27 Oct 2005 / 1.0
    !!==========================================================================

    use IdTypeModule, only : lId_p, getId

    type(MasterSlaveJointType), intent(in) :: joint
    character(len=15+lId_p)                :: GetJointId

    !! --- Logic section ---

    GetJointId = trim(jointTypeName_p(joint%type)) // getId(joint%id)

  end function GetJointId


  recursive function GetNumberOfMasterDOFs (joint) result (nDOFs)

    !!==========================================================================
    !! Return the total number of independent master DOFs of a joint.
    !! The function works only after the joint chaining has been resolved.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 1 Nov 2005 / 1.0
    !!==========================================================================

    type(MasterSlaveJointType), intent(in) :: joint
    integer                                :: nDOFs

    !! Local variables
    integer :: mt

    !! --- Logic section ---

    nDOFs = joint%nJointDOFs
    do mt = 1, size(joint%JMTriads)
       if (.not. associated(joint%chain)) then
          nDOFs = nDOFs + joint%JMTriads(mt)%nDOFs
       else if (associated(joint%chain%STriad,joint%JMTriads(mt)%triad)) then
          nDOFs = nDOFs + GetNumberOfMasterDOFs(joint%chain)
       else
          nDOFs = nDOFs + joint%JMTriads(mt)%nDOFs
       end if
    end do

  end function GetNumberOfMasterDOFs


  recursive function HasZeroVelAcc (joint) result(zeroVA)

    !!==========================================================================
    !! Checks if the given joint has zero velocities and accelerations
    !! in all of its master- and slave DOFs.
    !! The function works only after the joint chaining has been resolved.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 7 Nov 2005 / 1.0
    !!==========================================================================

    type(MasterSlaveJointType), intent(in) :: joint
    logical                                :: zeroVA

    !! Local variables
    integer :: i, mt
    real(dp), parameter :: eps_p = 1.0e-16_dp

    !! --- Logic section ---

    zeroVA = .false.

    !! Check the joint DOFs
    do i = 1, joint%nJointDOFs
       if (abs(joint%jointDofs(i)%jVar(2)) > eps_p) return
       if (abs(joint%jointDofs(i)%jVar(3)) > eps_p) return
    end do

    !! Check the slave triad DOFs
    do i = 1, size(joint%slaveDofs)
       if (abs(joint%STriad%urd(i))  > eps_p) return
       if (abs(joint%STriad%urdd(i)) > eps_p) return
    end do

    !! Check the master triad DOFs
    do mt = 1, size(joint%JMTriads)
       if (associated(joint%chain)) then
          if (associated(joint%chain%STriad,joint%JMTriads(mt)%triad)) then
             if (HasZeroVelAcc(joint%chain)) then
                cycle
             else
                return
             end if
          end if
       end if
       do i = 1, joint%JMTriads(mt)%nDOFs
          if (abs(joint%JMTriads(mt)%triad%urd(i))  > eps_p) return
          if (abs(joint%JMTriads(mt)%triad%urdd(i)) > eps_p) return
       end do
    end do

    zeroVA = .true.

  end function HasZeroVelAcc


  function transVSlaveToJoint (joint,u) result(v)

    !!==========================================================================
    !! Transform vector u from global directions in slave triad to joint
    !! directions, accounting for possible eccentricity between the joint
    !! position and the slave triad position.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 23 May 2008/1.0
    !!==========================================================================

    use RotationModule, only : EccExpand

    type(MasterSlaveJointType), intent(in) :: joint
    real(dp)                  , intent(in) :: u(:)

    !! Local variables
    real(dp) :: v(6)

    !! --- Logic section ---

    if (size(u) < 3) return

    !! Transform the vector to the joint position accounting for eccentricity
    call EccExpand (joint%JPosInG(:,4)-joint%STriad%ur(:,4),u(1:3),v)
    if (size(u) >= 6) v(4:6) = v(4:6) + u(4:6)

    !! Transform to local joint directions
    v(1:3) = matmul(v(1:3),joint%JPosInG(:,1:3))
    v(4:6) = matmul(v(4:6),joint%JPosInG(:,1:3))

  end function transVSlaveToJoint


  subroutine updatePreviousJointValues (joint)

    !!==========================================================================
    !! Updates the state variables pertaining to the previous time step.
    !! This subroutine is invoked once after convergence as been achieved.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 28 Oct 2008 / 1.0
    !!==========================================================================

    use FrictionTypeModule, only : updateAtConvergence

    type(MasterSlaveJointType), intent(inout) :: joint

    !! Local variables
    integer :: j

    !! --- Logic section ---

    do j = 1, size(joint%jointDofs)
       joint%jointDofs(j)%jVarPrev = joint%jointDofs(j)%jVar(1)
    end do

    if (associated(joint%multiDofFriction)) then
       call updateAtConvergence (joint%multiDofFriction)
    end if

  end subroutine updatePreviousJointValues


  !!============================================================================
  !> @brief Restores the state variables from the last converged time step.
  !>
  !> @param joint The joint to restore the state variables for
  !>
  !> @details This subroutine is invoked when doing iteration cut-back.
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 2 Nov 2008

  subroutine restorePreviousJointValues (joint)

    use FrictionTypeModule, only : restoreFromLastStep
    use RotationModule    , only : vec_to_mat

    type(MasterSlaveJointType), intent(inout) :: joint

    !! Local variables
    integer  :: jDof, lDof, iMat
    real(dp) :: jVar

    !! --- Logic section ---

    do jDof = 1, joint%nJointDOFs
       jVar = joint%jointDofs(jDof)%jVarPrev
       lDof = joint%jointDofs(jDof)%lDof
       iMat = joint%jointDofs(jDof)%iMat
       joint%jointDofs(jDof)%jVar(1) = jVar
       if (lDof <= 3) then
          joint%PosFromJVars(lDof,4,iMat) = jVar
       else if (lDof <= 6) then
          joint%thetaVec(lDof-3,iMat) = jVar
       end if
    end do

    do iMat = 1, joint%nMats
       call vec_to_mat (joint%thetaVec(:,iMat),joint%PosFromJVars(:,1:3,iMat))
    end do

    if (associated(joint%multiDofFriction)) then
       call restoreFromLastStep (joint%multiDofFriction)
    end if

  end subroutine restorePreviousJointValues

end module MasterSlaveJointTypeModule
