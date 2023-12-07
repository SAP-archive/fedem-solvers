!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file dummy_pyplot.f90
!> @brief Dummy module used when built without python plotting support
!> @cond NO_DOCUMENTAION
module pyplot_module
  type pyplot
    integer :: dummy
  end type pyplot
end module pyplot_module
!> @endcond
