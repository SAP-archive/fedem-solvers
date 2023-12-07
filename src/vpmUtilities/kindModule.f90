!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file kindModule.f90
!> @brief Global kind- and variable-size parameters.

!!==============================================================================
!> @brief Module with kind-parameters and parameters for size of variable types.

module KindModule

  implicit none

  integer , parameter :: i4 = selected_int_kind(9)   !< 4-byte integer
#ifdef FT_HAS_INT4_ONLY
  integer , parameter :: i8 = selected_int_kind(9)   !< 4-byte integer
#else
  integer , parameter :: i8 = selected_int_kind(18)  !< 8-byte integer
#endif
  integer , parameter :: dp = kind(1.0D0)            !< 8-byte real (double)
  integer , parameter :: sp = kind(1.0)              !< 4-byte real (single)

#ifdef FT_USE_INT8
  integer , parameter :: nbi_p = 8 !< Number of bytes in default integer kind
#else
  integer , parameter :: nbi_p = 4 !< Number of bytes in default integer kind
#endif
  integer , parameter :: nbs_p = 4 !< Number of bytes in single precision real
  integer , parameter :: nbd_p = 8 !< Number of bytes in double precision real

  integer , parameter :: lfnam_p = 256            !< Max. length of file names

  integer , parameter :: maxInt_p  = huge(1)      !< The largest integer value
  real(dp), parameter :: hugeVal_p = huge(1.0_dp) !< The largest real number

  real(dp), parameter :: pi_p = 4.0_dp*atan(1.0_dp) !< The value of &pi;

  !> Tolerance to be checked for division by zero.
  !> Usage: if (abs(denominator) > abs(numerator*epsDiv0_p)) then division OK
  real(dp), parameter :: epsDiv0_p = epsilon(1.0_dp)

end module KindModule
