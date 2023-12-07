!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file vpmCommon/environmentTypeModule.f90
!> @brief Namelist for reading environment data.

!!==============================================================================
!> @brief Module with a namelist for reading environment data.

module EnvironmentNamelistModule

  use KindModule, only : dp

  implicit none

  !! Define the ENVIRONMENT namelist
  real(dp) :: rhoAir      !< Mass density for the ambient air
  real(dp) :: rhoWater    !< Mass density for sea water
  real(dp) :: rhoInternal !< Mass density for internal fluids
  real(dp) :: rhoGrowth   !< Mass density for marine growth
  real(dp) :: tGrowth     !< Marine growth thickness
  real(dp) :: zGrowth(2)  !< Depth level (start and stop) for the marine growth
  real(dp) :: gravity(3)  !< Gravitation vector
  real(dp) :: seaDepth    !< Water depth
  real(dp) :: sea0        !< Global Z-coordinate of the sea surface
  real(dp) :: seaCS(4,3)  !< Coordinate system for wave and current

  integer :: waveTheory      !< Wave theory identifier
  integer :: waveFunction    !< BaseId of the wave function
  integer :: waveScaleEngine !< BaseId of the wave scale function
  integer :: hdfScaleEngine  !< BaseId of the hydrodynamic force scale function
  integer :: currFunction    !< BaseId of the current function
  integer :: currDirFunction !< BaseId of the current direction function
  integer :: currScaleEngine !< BaseId of the current scale function

  namelist /ENVIRONMENT/ rhoAir, rhoWater, rhoInternal, rhoGrowth, &
       &                 tGrowth, zGrowth, gravity, seaDepth, sea0, seaCS, &
       &                 waveTheory, waveFunction, waveScaleEngine, &
       &                 currFunction, currDirFunction, currScaleEngine, &
       &                 hdfScaleEngine

contains

  !> @cond FULL_DOC
  !> @brief Reads the ENVIRONMENT namelist from the given file unit.
  subroutine read_ENVIRONMENT (infp,stat)

    integer, intent(in)  :: infp
    integer, intent(out) :: stat

    !! Default values
    rhoAir=0.0_dp; rhoWater=0.0_dp; rhoInternal=0.0_dp; rhoGrowth=0.0_dp
    tGrowth=0.0_dp; zGrowth=0.0_dp; gravity=0.0_dp
    seaDepth=0.0_dp; sea0=0.0_dp; seaCS=0.0_dp
    waveTheory=1; waveFunction=0; waveScaleEngine=0
    currFunction=0; currDirFunction=0; currScaleEngine=0; hdfScaleEngine=0

    read(infp,nml=ENVIRONMENT,iostat=stat)

  end subroutine read_ENVIRONMENT
  !> @endcond

end module EnvironmentNamelistModule
