!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module dbgUnitsModule

  !!============================================================================
  !! This module contains all file unit numbers used for various debug print.
  !! They are gathered in a separate module to better control their deallocation
  !! on program termination (without actual termination of the running process).
  !!============================================================================

  use KindModule, only : dp

  implicit none

  integer, save :: dbgShadowPos = 0
  integer, save :: dbgJoint = 0
  integer, save :: dbgCurve = 0
  integer, save :: dbgFric = 0
  integer, save :: dbgSolve = 0
  integer, save :: dbgTire = 0
  integer, save :: dbgRoad = 0
  integer, save :: dbgRot = 0
  integer, save :: dbgUDE = 0
  integer, save :: dbgAD = 0
  integer, save :: dbgHD = 0

  integer , save :: dbgIter = 0
  real(dp), save :: dbgTime = 0.0_dp

contains

  subroutine closeDbgUnits

    !!==========================================================================
    !! Close all debug files
    !!
    !! Programmer : Knut Morten Okstad                 date/rev : 1 Feb 2017/1.0
    !!==========================================================================

    use FileUtilitiesModule, only : closeDBGfiles

    call closeDBGfiles ()

    dbgShadowPos = 0
    dbgJoint = 0
    dbgCurve = 0
    dbgFric = 0
    dbgSolve = 0
    dbgTire = 0
    dbgRoad = 0
    dbgRot = 0
    dbgUDE = 0
    dbgAD = 0
    dbgHD = 0

  end subroutine closeDBGUnits

end module dbgUnitsModule
