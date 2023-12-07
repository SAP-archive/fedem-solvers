!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file isoMatModule.f90
!> @brief Constitutive matrices for linear-elastic isotropic materials.

!!==============================================================================
!> @brief Module with subroutines for isotropic constitutive matrices.

module isoMatModule

  implicit none

contains

  !> @brief 2D constitutive matrix for linear-elastic isotropic materials.
  subroutine isoMat2D (Emod,Rnu,C)

    use KindModule, only : dp

    real(dp), intent(in)  :: Emod,Rnu
    real(dp), intent(out) :: C(:,:)

    C = 0.0_dp

    C(1,1) = Emod / (1.0_dp - Rnu*Rnu)
    C(1,2) = Rnu * C(1,1)

    C(2,1) = C(1,2)
    C(2,2) = C(1,1)

    C(3,3) = 0.5_dp * Emod / (1.0_dp + Rnu)

  end subroutine isoMat2D


  !> @brief Invers 2D constitutive matrix for linear-elastic isotropic material.
  subroutine isoMat2Dinv (Emod,Rnu,C)

    use KindModule, only : dp

    real(dp), intent(in)  :: Emod,Rnu
    real(dp), intent(out) :: C(:,:)

    C = 0.0_dp

    C(1,1) = 1.0_dp / Emod
    C(1,2) = -Rnu / Emod

    C(2,1) = C(1,2)
    C(2,2) = C(1,1)

    C(3,3) = 2.0_dp * (1.0_dp + Rnu) / Emod

  end subroutine isoMat2Dinv


  !> @brief 3D constitutive matrix for linear-elastic isotropic materials.
  subroutine isoMat3D (Emod,Rnu,C)

    use KindModule, only : dp

    real(dp), intent(in)  :: Emod,Rnu
    real(dp), intent(out) :: C(6,6)

    real(dp) :: fac

    C   = 0.0_dp
    fac = Emod / ((1.0_dp + Rnu) * (1.0_dp - Rnu - Rnu))

    C(1,1) = (1.0_dp - Rnu) * fac
    C(2,1) = Rnu * fac
    C(3,1) = C(2,1)

    C(1,2) = C(2,1)
    C(2,2) = C(1,1)
    C(3,2) = C(2,1)

    C(1,3) = C(2,1)
    C(2,3) = C(2,1)
    C(3,3) = C(1,1)

    C(4,4) = (0.5_dp - Rnu) * fac
    C(5,5) = C(4,4)
    C(6,6) = C(4,4)

  end subroutine isoMat3D


  !> @brief Invers 3D constitutive matrix for linear-elastic isotropic material.
  subroutine isoMat3Dinv (Emod,Rnu,C)

    use KindModule, only : dp

    real(dp), intent(in)  :: Emod,Rnu
    real(dp), intent(out) :: C(6,6)

    C = 0.0_dp

    C(1,1) = 1.0_dp / Emod
    C(2,1) = -Rnu / Emod
    C(3,1) = C(2,1)

    C(1,2) = C(2,1)
    C(2,2) = C(1,1)
    C(3,2) = C(2,1)

    C(1,3) = C(2,1)
    C(2,3) = C(2,1)
    C(3,3) = C(1,1)

    C(4,4) = 2.0_dp * (1.0_dp + Rnu) / Emod
    C(5,5) = C(4,4)
    C(6,6) = C(4,4)

  end subroutine isoMat3Dinv

end module isoMatModule
