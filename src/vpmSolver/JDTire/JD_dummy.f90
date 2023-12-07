!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file JD_dummy.f90
!> @brief Dummy module used when built without JD tire support
!> @cond NO_DOCUMENTAION
module JD_tireRoutinesModule
  implicit none
contains
  subroutine jdt_stiffness (lDyn,tire,s,scaleK,eM)
    use TireTypeModule, only : TireType, dp
    logical       , intent(in)  :: lDyn
    type(TireType), intent(in)  :: tire
    real(dp)      , intent(in)  :: s(:), scaleK
    real(dp)      , intent(out) :: eM(:,:)
    print *,' ** jdt_stiffness dummy:',lDyn,tire%id%baseId,size(s),scaleK
    eM = 0.0_dp
  end subroutine jdt_stiffness
end module JD_tireRoutinesModule
!> @endcond
