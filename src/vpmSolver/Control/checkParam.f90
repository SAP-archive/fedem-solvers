!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file checkParam.f90
!> @brief Checks the validity of control parameter inputs.

!!==============================================================================
!> @brief Checks the validity of control parameter inputs.
!> @details This global function is used to determine whether the number of
!> real parameters, @a nRealData, and the number of variables, @a nVar,
!> specified on the solver input file, are correct for the specified control
!> element type, and that the element type number itself, @a elmType, is valid.

function checkCtrlParams (elmType,nRealData,nVar)
  !DEC$ ATTRIBUTES DLLEXPORT :: checkCtrlParams

  implicit none

  integer, intent(inout) :: elmType
  integer, intent(in)    :: nRealData, nVar

  integer :: nR, nV
  logical :: checkPrm, checkCtrlParams

  checkPrm(nR,nV) = nRealData == nR-1 .and. nVar == nV

  checkCtrlParams = .true.

  select case (elmType)
  case (1,2,5)
     if (checkPrm(2,3)) return
  case (3,4,7,11)
     if (checkPrm(2,2)) return
  case (6)
     if (checkPrm(1,3)) return
  case (12)
     if (checkPrm(7,2)) return
  case (21)
     if (checkPrm(4,2)) return
  case (22,23,41)
     if (checkPrm(3,2)) return
  case (24)
     if (checkPrm(10,2)) return
  case (31)
     if (checkPrm(3,3)) return
  case (32,34,43)
     if (checkPrm(4,4)) return
  case (33)
     if (checkPrm(3,4)) return
  case (35)
     if (checkPrm(4,5)) return
  case (36)
     if (checkPrm(5,5)) return
  case (37)
     if (checkPrm(6,5)) return
  case (42)
     if (checkPrm(4,3)) return
  case (44)
     if (checkPrm(6,6)) return
  case (45,46)
     return ! TODO: more testing
  case default
     elmType = -abs(elmType) ! Unknown element type
  end select

  checkCtrlParams = .false.

end function checkCtrlParams
