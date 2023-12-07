!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file dummy2Modules.f90
!> @brief Dummy modules used to reduce dependencies in unit tests.
!>
!> @details This file contains dummy implementations of some modules,
!> to facilitate unit testing with a minimum of dependencies. The modules here
!> only contain some entries with empty/no bodies to avoid undefined symbols.
!>
!> Dummies for the following modules and entries are included here:
!> - @ref environmenttypemodule
!>      - environmenttypemodule::environmenttype
!> - @ref hydrodynamicsmodule
!>      - hydrodynamicsmodule::getseastate
!> - @ref windturbineroutinesmodule
!>      - windturbineroutinesmodule::getwindspeed
!>
!> @cond NO_DOCUMENTATION

!!==============================================================================
!> @brief Dummy environment data module.

module EnvironmentTypeModule

  implicit none

  type EnvironmentType
     integer :: dummy
  end type EnvironmentType

end module EnvironmentTypeModule


!!==============================================================================
!> @brief Dummy hydrodynamics module.

module HydroDynamicsModule

  implicit none

contains

  subroutine getSeaState (env,time,Xg,sea,err)
    use kindModule           , only : dp
    use EnvironmentTypeModule, only : EnvironmentType
    type(EnvironmentType), intent(in)  :: env
    real(dp)             , intent(in)  :: time, Xg(3)
    real(dp)             , intent(out) :: sea(:)
    integer              , intent(out) :: err
    print*,'Dummy getSeaState',time,Xg,env%dummy
    sea = 0.0_dp
    err = -999
  end subroutine getSeaState

end module HydroDynamicsModule


!!==============================================================================
!> @brief Dummy wind turbine routines module.

module WindTurbineRoutinesModule

  implicit none

contains

  subroutine getWindSpeed (pos,time,dws,uws,err)
    use kindModule, only  : dp
    real(dp), intent(in)  :: pos(3), time
    real(dp), intent(out) :: dws(3), uws(3)
    integer , intent(out) :: err
    print*,'Dummy getWindSpeed',pos,time
    dws = 0.0_dp
    uws = 0.0_dp
    err = -999
  end subroutine getWindSpeed

end module WindTurbineRoutinesModule

!> @endcond
