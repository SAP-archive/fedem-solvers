!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

subroutine response_pos (ierr)

  !!============================================================================
  !! This program generates the response_pos.frs file used by fedem_main.
  !!
  !! Programmer : Knut Morten Okstad
  !! date/rev   : May 2001/1.0
  !!============================================================================

  use TireTypeModule            , only : STI_p, MF_p, dp
  use MasterSlaveJointTypeModule, only : FREE_p, BALL_p, REVOLUTE_p
  use MasterSlaveJointTypeModule, only : PRISMATIC_p, CYLINDRIC_p, CAM_p
  use MechanismTypeModule       , only : MechanismType
  use ControlTypeModule         , only : ControlType
  use ModesTypeModule           , only : ModesType
  use MechanismTypeModule       , only : nullifyMechanism, deallocateMechanism
  use ControlTypeModule         , only : nullifyCtrl, deallocateCtrl
  use ModesTypeModule           , only : nullifyModes, deallocateModes
  use TriadTypeModule           , only : nullifyTriad, allocateTriadForces
  use SupElTypeModule           , only : nullifySupEl
  use ForceTypeModule           , only : nullifyForce
  use SpringTypeModule          , only : nullifySpring
  use DamperTypeModule          , only : nullifyDamper
  use MasterSlaveJointTypeModule, only : nullifyJoint, nullifyJointDof
  use TireTypeModule            , only : nullifyTire
  use IdTypeModule              , only : nullifyId
  use ContactElementTypeModule  , only : nullifyContactElement
  use FunctionTypeModule        , only : nullifyEngine
  use FrictionTypeModule        , only : nullifyFriction
  use RDBModule                 , only : RDBType, nullifyRDB, closeRDBfile
  use saveModule                , only : writeSolverHeaders, res2DB
  use StrainRosetteModule       , only : StrainElementType, nullifyRosette
  use saveStrainGageModule      , only : writeStrainGageHeader
  use reportErrorModule         , only : setErrorFile
  use timerModule               , only : initTime, showTime

  implicit none

  !! Local variables
  integer                    :: i, j, ierr
  type(MechanismType),target :: mech
  type(ControlType)          :: ctrl
  type(ModesType)            :: modes
  type(StrainElementType)    :: rosettes(1)
  type(RDBType)              :: rdb

  !! --- Logic section ---

  call initTime()
  call setErrorFile(6)

  !! Establish a dummy mechanism containing one instance of each object type,
  !! such that its frs-header will cover all possibilities of solver results.

  call nullifyMechanism(mech)
  call nullifyModes(modes)
  call nullifyCtrl(ctrl)
  mech%saveVar(1:5) = .true.
  mech%hasEinp = .true.
  mech%hasEext = .true.

  !! To avoid memory leakage due to the 0-length allocation in nullifyMechanism
  deallocate(mech%triads,mech%sups,mech%joints,mech%cElems, &
       &     mech%axialSprings,mech%baseSprings,mech%dampers,mech%baseDampers, &
       &     mech%frictions,mech%forces,mech%tires,mech%engines)

  allocate(mech%triads(2),mech%sups(1),mech%joints(6),mech%cElems(1), &
       &   mech%axialSprings(1),mech%baseSprings(7),mech%dampers(7), &
       &   mech%baseDampers(7),mech%frictions(3),mech%forces(3),mech%tires(1), &
       &   mech%engines(1),ctrl%vreg(1),ctrl%vregId(1), STAT=ierr)
  if (ierr /= 0) stop 'Allocation error'

  do i = 1, size(mech%triads)
     call nullifyTriad(mech%triads(i))
     mech%triads(i)%nDOFS = 6
     mech%triads(i)%saveVar(1:6) = .true.
  end do
  mech%triads(1)%id%baseId = -1 ! Dummy CG triad

  do i = 1, size(mech%sups)
     call nullifySupEl(mech%sups(i))
     allocate(mech%sups(i)%triads(2), STAT=ierr)
     if (ierr /= 0) stop 'Allocation error'
     mech%sups(i)%addedMass = .true.
     mech%sups(i)%saveVar = .true.
     do j = 1, size(mech%sups(i)%triads)
        mech%sups(i)%triads(j)%p => mech%triads(j)
        call allocateTriadForces(mech%sups(i)%triads(j)%p,ierr)
        if (ierr /= 0) stop 'Allocation error'
     end do
  end do

  do i = 1, size(mech%forces)
     call nullifyForce(mech%forces(i))
     mech%forces(i)%dof = -i
     mech%forces(i)%saveVar = .true.
  end do
  mech%forces(3)%dof = 1
  mech%forces(3)%joint => mech%joints(1)

  do i = 1, size(mech%axialSprings)
     call nullifySpring(mech%axialSprings(i))
     allocate(mech%axialSprings(i)%id, mech%axialSprings(i)%spr(1), STAT=ierr)
     if (ierr /= 0) stop 'Allocation error'
     call nullifyId(mech%axialSprings(i)%id)
     mech%axialSprings(i)%samElNum = i ! Axial spring
     mech%axialSprings(i)%spr(1)%p => mech%baseSprings(7)
  end do

  do i = 1, size(mech%baseSprings)
     call nullifySpring(mech%baseSprings(i))
     if (i <= 6) mech%baseSprings(i)%dof = i ! Joint spring
     mech%baseSprings(i)%length0Engine => mech%engines(1)
     mech%baseSprings(i)%saveVar = .true.
  end do

  do i = 1, size(mech%dampers)
     call nullifyDamper(mech%dampers(i))
     call nullifyDamper(mech%baseDampers(i))
     mech%dampers(i)%dmp => mech%baseDampers(i)
     if (i > 6) then ! Axial damper
        mech%dampers(i)%samElNum = i
     else            ! Joint damper
        mech%dampers(i)%dmp%dof = i
     end if
     mech%dampers(i)%dmp%saveVar = .true.
  end do

  do i = 1, size(mech%joints)
     call nullifyJoint(mech%joints(i))
  end do

  do i = 1, size(mech%cElems)
     call nullifyContactElement(mech%cElems(i))
     mech%cElems(i)%saveVar = .true.
  end do

  do i = 1, size(mech%tires)
     call nullifyTire(mech%tires(i))
     mech%tires(i)%api = STI_p
     mech%tires(i)%type = MF_p
     mech%tires(i)%saveVar = .true.
  end do

  do i = 1, size(mech%engines)
     call nullifyEngine(mech%engines(i))
     mech%engines(i)%saveVar = 2
  end do

  do i = 1, size(mech%frictions)
     allocate(mech%frictions(i)%p, STAT=ierr)
     if (ierr /= 0) stop 'Allocation error'
     call nullifyFriction(mech%frictions(i)%p)
     mech%frictions(i)%p%saveVar(1:2) = .true.
  end do

  do i = 1, size(ctrl%vregId)
     call nullifyId(ctrl%vregId(i))
     ctrl%vregId(i)%baseId = i
  end do
  ctrl%saveVar = .true.

  allocate(mech%joints(1)%jointDofs(6), &
       &   mech%joints(2)%jointDofs(3), &
       &   mech%joints(3)%jointDofs(1), &
       &   mech%joints(4)%jointDofs(3),mech%joints(4)%slider, &
       &   mech%joints(5)%jointDofs(4),mech%joints(5)%slider, &
       &   mech%joints(6)%jointDofs(6),mech%joints(6)%slider, STAT=ierr)
  if (ierr /= 0) stop 'Allocation error'

  do i = 1, size(mech%joints)
     mech%joints(i)%nJointDOFs = size(mech%joints(i)%jointDofs)
     allocate(mech%joints(i)%dofOutputOrder(mech%joints(i)%nJointDOFs))
     do j = 1, mech%joints(i)%nJointDOFs
        mech%joints(i)%dofOutputOrder(j) = j
        call nullifyJointDof(mech%joints(i)%jointDofs(j))
        if (i <= 3) then
           mech%joints(i)%jointDofs(j)%spring => mech%baseSprings(j)
           mech%joints(i)%jointDofs(j)%damper => mech%dampers(j)%dmp
        end if
        mech%joints(i)%jointDofs(j)%loadIdx = 3
        mech%joints(i)%jointDofs(j)%saveVar = .true.
        !TODO,kmo: Cam joints are now represented by the ContactElement type.
        !TODO,kmo: The old master-slave-type cam formulation is now beta feature
        !TODO,kmo: only, and the associated response possibilities are disabled.
        !TODO,kmo: Remove the following line if it should be enabled again,
        !TODO,kmo: and reinstate the out-commented code below.
        if (i == 6) mech%joints(i)%jointDofs(j)%saveVar = .false.
     end do
  end do

  mech%joints(1)%type = FREE_p
  mech%joints(1)%jointDofs(1)%lDof = 1
  mech%joints(1)%jointDofs(2)%lDof = 2
  mech%joints(1)%jointDofs(3)%lDof = 3
  mech%joints(1)%jointDofs(4)%lDof = 4
  mech%joints(1)%jointDofs(5)%lDof = 5
  mech%joints(1)%jointDofs(6)%lDof = 6

  mech%joints(2)%type = BALL_p
  mech%joints(2)%jointDofs(1)%lDof = 4
  mech%joints(2)%jointDofs(2)%lDof = 5
  mech%joints(2)%jointDofs(3)%lDof = 6

  mech%joints(3)%type = REVOLUTE_p
  mech%joints(3)%jointDofs(1)%lDof = 6
  mech%joints(3)%jointDofs(1)%friction => mech%frictions(1)%p

  mech%joints(4)%type = PRISMATIC_p
  mech%joints(4)%jointDofs(1)%lDof = 4
  mech%joints(4)%jointDofs(2)%lDof = 5
  mech%joints(4)%jointDofs(3)%lDof = 7
  mech%joints(4)%slider%slideDof => mech%joints(4)%jointDofs(3)
  mech%joints(4)%slider%slideDof%spring => mech%baseSprings(3)
  mech%joints(4)%slider%slideDof%damper => mech%dampers(3)%dmp
  mech%joints(4)%slider%slideDof%friction => mech%frictions(2)%p
  nullify(mech%joints(4)%slider%coeff)

  mech%joints(5)%type = CYLINDRIC_p
  mech%joints(5)%jointDofs(1)%lDof = 4
  mech%joints(5)%jointDofs(2)%lDof = 5
  mech%joints(5)%jointDofs(3)%lDof = 6
  mech%joints(5)%jointDofs(3)%spring => mech%baseSprings(1)
  mech%joints(5)%jointDofs(3)%damper => mech%dampers(1)%dmp
  mech%joints(5)%jointDofs(4)%lDof = 7
  mech%joints(5)%slider%slideDof => mech%joints(5)%jointDofs(4)
  mech%joints(5)%slider%slideDof%spring => mech%baseSprings(3)
  mech%joints(5)%slider%slideDof%damper => mech%dampers(3)%dmp
  nullify(mech%joints(5)%slider%coeff)

  mech%joints(6)%type = CAM_p
  !mech%joints(6)%jointDofs(1)%lDof = 1
  !TODO,kmo: See TODO note above
  !mech%joints(6)%jointDofs(1)%spring => mech%springs(1)
  !mech%joints(6)%jointDofs(1)%damper => mech%dampers(1)
  !mech%joints(6)%jointDofs(2)%lDof = 2
  !mech%joints(6)%jointDofs(2)%spring => mech%springs(2)
  !mech%joints(6)%jointDofs(2)%damper => mech%dampers(2)
  !mech%joints(6)%jointDofs(3)%lDof = 4
  !mech%joints(6)%jointDofs(4)%lDof = 5
  !mech%joints(6)%jointDofs(5)%lDof = 6
  !mech%joints(6)%jointDofs(6)%lDof = 7
  !mech%joints(6)%slider%slideDof => mech%joints(6)%jointDofs(6)
  !mech%joints(6)%slider%slideDof%friction => mech%frictions(2)%p
  nullify(mech%joints(6)%slider%coeff)

  mech%cElems(1)%springs(1)%p => mech%baseSprings(1)
  mech%cElems(1)%dampers(1)%p => mech%dampers(1)%dmp
  mech%cElems(1)%springs(2)%p => mech%baseSprings(2)
  mech%cElems(1)%dampers(2)%p => mech%dampers(2)%dmp
  mech%cElems(1)%springs(3)%p => mech%baseSprings(3)
  mech%cElems(1)%dampers(3)%p => mech%dampers(3)%dmp
  mech%cElems(1)%friction => mech%frictions(3)%p

  do i = 1, size(rosettes)
     call nullifyRosette(rosettes(i))
     allocate(rosettes(1)%data(1),STAT=ierr)
     if (ierr /= 0) stop 'Allocation error'
     allocate(rosettes(1)%data(1)%gages(3),STAT=ierr)
     if (ierr /= 0) stop 'Allocation error'
  end do

  !! Now write the frs-file header to file response_pos.frs
  call writeSolverHeaders ((/'response_pos.frs'/),mech,ctrl,modes, &
       &                   (/0.0_dp/),0,(/.false./),ierr)
  if (ierr /= 0) stop 'writeSolverHeaders failed'

  !! Now write the gage frs-file header to file gage_pos.frs
  call nullifyRDB(rdb)
  rdb%nVar = res2DB%nVar-2
  rdb%nIG  = res2DB%nIG
  call writeStrainGageHeader ('gage_pos.frs','dummy.fmm',rdb,rosettes, &
       &                      0.0_dp,ierr,.false.)
  if (ierr /= 0) stop 'writeStrainGageHeader failed'

  !! Close files and clean up memory
  call closeRDBfile (rdb,ierr,.true.)
  call closeRDBfile (res2DB,ierr,.true.)
  do i = 1, size(rosettes)
     deallocate(rosettes(1)%data(1)%gages)
     deallocate(rosettes(1)%data)
  end do
  call deallocateModes(modes)
  call deallocateCtrl(ctrl)
  call deallocateMechanism(mech)
  call showTime(6)

end subroutine response_pos
