!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file curvePointTypeModule.f90
!>
!> @brief Representation of contact curve/surface points.

!!==============================================================================
!> @brief Module with a data type representing contact curve/surface points.

module CurvePointTypeModule

  use TriadTypeModule, only : TriadType, dp

  implicit none

  !> @brief Data type representing a point on a contact curve/surface.
  type CurvePointType
     type(TriadType), pointer :: triad !< Pointer to triad object of this point
     real(dp) :: CposInT(3,4) !< Curve point position relative to this triad
     real(dp) :: slideVar     !< Slide variable value at this curve point
     real(dp) :: curvature    !< Curvature of segment starting at this point
     real(dp) :: h0           !< Initial arc heigth from secant line
     real(dp) :: upVecInT1(3) !< Up-direction of curvature when first node
     real(dp) :: upVecInT2(3) !< Up-direction of curvature when second node
     real(dp) :: tan1G(3)     !< Tangent vector in global when first node
     real(dp) :: tan2G(3)     !< Tangent vector in global when second node
  end type CurvePointType

end module CurvePointTypeModule
