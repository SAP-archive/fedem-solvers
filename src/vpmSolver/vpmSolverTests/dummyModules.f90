!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file dummyModules.f90
!> @brief Dummy modules used to reduce dependencies in unit tests.
!>
!> @details This file contains dummy implementations of some modules,
!> to facilitate unit testing with a minimum of dependencies. The modules here
!> only contain some entries with empty/no bodies to avoid undefined symbols.
!>
!> Dummies for the following modules and entries are included here:
!> - @ref environmenttypemodule
!>      - environmenttypemodule::environmenttype
!>      - environmenttypemodule::nullifyenvironment
!>      - environmenttypemodule::glob2sea
!> - @ref functiontypemodule
!>      - functiontypemodule::enginetype
!>      - functiontypemodule::getptrtoid
!> - @ref engineroutinesmodule
!>      - engineroutinesmodule::evaluate
!> - @ref masterslavejointtypemodule
!>      - masterslavejointtypemodule::jointdoftype
!>      - masterslavejointtypemodule::masterslavejointtype
!>      - masterslavejointtypemodule::getptrtoid
!>      - masterslavejointtypemodule::getjointvar
!>
!> @cond NO_DOCUMENTATION

!!==============================================================================
!> @brief Dummy environment data module.
!> @details This module only contains items needed by the finiteelementmodule
!> such that unit tests for that module can be made with minimum dependencies.

module EnvironmentTypeModule

  use kindModule, only : dp

  implicit none

  !> @brief Dummy environmental data type containing marine growth items only.
  type EnvironmentType
     real(dp) :: rhoG, tG, zG(2)
  end type EnvironmentType

contains

  subroutine NullifyEnvironment (env)
    type(EnvironmentType), intent(out) :: env
    env%rhoG = 0.0_dp
    env%tG = 0.0_dp
    env%zG = 0.0_dp
  end subroutine NullifyEnvironment

  function glob2Sea (env,xg)
    type(EnvironmentType), intent(in) :: env
    real(dp) :: xg(3), glob2Sea(3)
    print*,'Dummy glob2Sea, rhoG =',env%rhoG
    glob2Sea = xg
  end function glob2Sea

end module EnvironmentTypeModule

!!==============================================================================
!> @brief Dummy function data module.
!> @details This module only contains the items needed by the motion unit test.

module FunctionTypeModule

  use IdTypeModule, only : IdType

  implicit none

  type EngineType
     type(IdType) :: id
  end type EngineType

  interface GetPtrToId
     module procedure GetPtrToIdEngine
  end interface

  private :: IdType, GetPtrToIdEngine

contains

  function GetPtrToIdEngine (array,id) result(ptr)
    type(EngineType), pointer :: ptr
    type(EngineType), target  :: array(:)
    integer :: id, i
    do i = 1, size(array)
       if (array(i)%id%baseId == id) then
          ptr => array(i)
          return
       end if
    end do
    nullify(ptr)
  end function GetPtrToIdEngine

end module FunctionTypeModule

!!==============================================================================
!> @brief Dummy function evaluation module.
!> @details This module only contains the items needed by the motion unit test.

module EngineRoutinesModule

  implicit none

contains

  function Evaluate (engine,E1,E0,ierr) result(eVal)
    use kindModule        , only : dp
    use FunctionTypeModule, only : EngineType
    type(EngineType), pointer :: engine
    integer  :: ierr
    real(dp) :: E1, E0, eVal
    print*,'Dummy Evaluate',engine%id%baseId,E0,E1
    eVal = E0
    ierr = 0
  end function Evaluate

end module EngineRoutinesModule


#ifdef FT_DEBUG
module dbgUnitsModule
  integer, save :: dbgSolve = 0
end module dbgUnitsModule
#endif

!!==============================================================================
!> @brief Dummy joint data module.
!> @details This module only contains the items needed by the motion unit test.

module MasterSlaveJointTypeModule

  use IdTypeModule, only : IdType

  implicit none

  !> @brief Dummy joint DOF data type.
  type JointDofType
     integer :: lDOF, loadIdx
  end type JointDofType

  !> @brief Dummy joint data type.
  type MasterSlaveJointType
     type(IdType)                :: id
     integer                     :: nJointDOFs
     type(JointDofType), pointer :: jointDOFs(:)
  end type MasterSlaveJointType

  interface GetPtrToId
     module procedure GetPtrToIdJoint
  end interface

  private :: IdType, GetPtrToIdJoint

contains

  function GetPtrToIdJoint (array,id) result(ptr)
    type(MasterSlaveJointType), pointer :: ptr
    type(MasterSlaveJointType), target  :: array(:)
    integer :: id, i
    do i = 1, size(array)
       if (array(i)%id%baseId == id) then
          ptr => array(i)
          return
       end if
    end do
    nullify(ptr)
  end function GetPtrToIdJoint

  function GetJointVar (joint,dofInd,type)
    use kindModule, only : dp
    type(MasterSlaveJointType), intent(in) :: joint
    integer, intent(in) :: dofInd, type
    real(dp) :: GetJointVar
    print*,'Dummy GetJointVar:', joint%id%baseId,dofInd,type
    GetJointVar = 0.0_dp
  end function GetJointVar

end module MasterSlaveJointTypeModule

!> @endcond
