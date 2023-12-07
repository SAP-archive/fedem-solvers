!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file headingNamelistModule.f90
!> @brief Namelist for reading model file name and version.

!!==============================================================================
!> @brief Module with a namelist for reading model file name and version.

module HeadingNamelistModule

  use KindModule, only : sp, lfnam_p

  implicit none

  !! Define the HEADING namelist
  character(lfnam_p) :: modelFile !< Model file name
  real(sp)           :: version   !< Model file version
  namelist /HEADING/ modelFile, version

contains

  !> @cond FULL_DOC
  !> @brief Reads the HEADING namelist from the given file unit.
  subroutine read_HEADING (infp,stat)

    integer, intent(in)  :: infp
    integer, intent(out) :: stat

    !! Default values
    modelFile=''; version=0.0_sp

    read(infp,nml=HEADING,iostat=stat)

  end subroutine read_HEADING
  !> @endcond

end module HeadingNamelistModule
