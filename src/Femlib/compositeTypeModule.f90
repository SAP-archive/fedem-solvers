!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module CompositeTypeModule

  !!============================================================================
  !! This module contains a data types representing a composite materiels, and
  !! subroutines for calculating constitutive matrices for composite shells.
  !!
  !! The code is developed by Bjorn Haugen (NTNU), April 2012.
  !! Adapted to FEDEM standard by Knut Morten Okstad, April 2016.
  !!============================================================================

  use KindModule, only : dp

  implicit none

  type OrthoMatType ! Contains parameters concerning each ply
     ! Parameter names are according to the Nastran bulk data entry MAT8
     integer  :: MID  ! Material identification number
     real(dp) :: E1   ! Youngs modulus, longitudinal
     real(dp) :: E2   ! Youngs modulus, transverse
     real(dp) :: NU12 ! Poisson's ratio 12
     real(dp) :: NU21 ! Poisson's ratio 21
     real(dp) :: G12  ! Shear modulus
     real(dp) :: RHO  ! Mass density
  end type OrthoMatType

  type OrthoMatPtrType
     type(OrthoMatType), pointer :: p
  end type OrthoMatPtrType

  type CompositeType ! Layered composite element property
     ! Parameter names are according to the Nastran bulk data entry PCOMP
     integer  :: PID  ! Property identification number
     integer  :: nPly ! Number of plies
     real(dp) :: Z0   ! Distance from the reference plane to the bottom surfac
     real(dp), pointer :: T(:)     ! Thickness of each ply
     real(dp), pointer :: THETA(:) ! Fiber direction aangle of each ply
     type(OrthoMatType), pointer :: pMat(:) ! Material properties of each ply
  end type CompositeType

  private :: trans_strain2D_deg, makeQ


contains

  subroutine trans_strain2D_deg (a_deg,T)

    !***************************************************************************
    ! Creates a 3x3 transformation matrix for strains in 2D.
    !***************************************************************************

    use KindModule, only : pi_p

    real(dp), intent(in)  :: a_deg
    real(dp), intent(out) :: T(3,3)

    real(dp) :: a, c, s

    a = a_deg*pi_p/180.0_dp
    s = sin(a)
    c = cos(a)

    T(1,1) =  c*c
    T(1,2) =  s*s
    T(1,3) = -s*c*2.0_dp
    T(2,1) =  T(1,2)
    T(2,2) =  T(1,1)
    T(2,3) = -T(1,3)
    T(3,1) =  s*c
    T(3,2) = -T(3,1)
    T(3,3) =  T(1,1) - T(1,2)

  end subroutine trans_strain2D_deg


  subroutine makeQ (Q,e1,e2,nu12,g12,theta)

    !***************************************************************************
    ! Makes the Q matrix for a ply, ref faserverbundwerkstoff (3.35) and (3.51).
    ! Input: e1, e2 - Young's modulus in both directions
    !        g12, nu12 - Shear modulus and poisson's ratio
    !        theta - Fiber direction of the ply
    ! Output: Q - The Q-matrix
    !***************************************************************************

    use KindModule, only : epsDiv0_p

    real(dp), intent(out) :: Q(3,3)
    real(dp), intent(in)  :: e1, e2, nu12, g12, theta

    real(dp) :: nu21, denom, QT(3,3), T(3,3)

    nu21  = nu12*(e2/e1)
    denom = 1.0_dp - nu12*nu21

    Q(1,1) = e1 / denom
    Q(1,2) = nu21*e1 / denom
    Q(1,3) = 0.0_dp
    Q(2,1) = Q(1,2)
    Q(2,2) = e2 / denom
    Q(2,3) = 0.0_dp
    Q(3,1) = 0.0_dp
    Q(3,2) = 0.0_dp
    Q(3,3) = g12

    if (abs(theta) < epsDiv0_p) return

    call trans_strain2D_deg (Theta,T)

    QT = matmul(Q,transpose(T))
    Q  = matmul(T,QT)

  end subroutine makeQ


  subroutine Cmatrix (laminate,thetae,Cmat)

    !***************************************************************************
    ! Makes the C matrix for a material, ref. faserverbundwerkstoff (4.17-19).
    !***************************************************************************

    type(CompositeType), intent(in)  :: laminate
    real(dp)           , intent(in)  :: thetae
    real(dp)           , intent(out) :: Cmat(6,6)

    real(dp) :: C(6,6), Q(3,3), trmat(6,6), hk, hkm1
    integer  :: k

    C = 0.0_dp
    do k = 1, laminate%nPly

       ! Make the Q-matrix for the given fiber direction, theta
       call makeQ (Q,laminate%pMat(k)%E1,laminate%pMat(k)%E2, &
            &      laminate%pMat(k)%NU12,laminate%pMat(k)%G12, &
            &      laminate%THETA(k))

       ! Calculate the thickness of each ply
       if (k > 1) then
          hkm1 = laminate%Z0 + sum(laminate%T(1:k-1))
       else
          hkm1 = laminate%Z0
       end if
       hk = hkm1 + laminate%T(k)

       ! Membrane stiffness (A)
       C(1:3,1:3) = C(1:3,1:3) + Q*(hk - hkm1)
       ! Coupling terms (B)**2
       C(1:3,4:6) = C(1:3,4:6) + Q*(hk*hk - hkm1*hkm1)
       ! Bending stiffness (C)**3
       C(4:6,4:6) = C(4:6,4:6) + Q*(hk*hk*hk - hkm1*hkm1*hkm1)

    end do

    C(1:3,4:6) = C(1:3,4:6)/2.0_dp ! adjust coupling terms
    C(4:6,1:3) = transpose(C(1:3,4:6))
    C(4:6,4:6) = C(4:6,4:6)/3.0_dp ! adjust bending terms

    ! Create the transformation matrix
    call trans_strain2D_deg (thetae,trmat(1:3,1:3))
    trmat(1:3,4:6) = 0.0_dp
    trmat(4:6,1:3) = 0.0_dp
    trmat(4:6,4:6) = trmat(1:3,1:3)

    ! Transform to element orientation
    Cmat = matmul(trmat,matmul(C,transpose(trmat)))

  end subroutine Cmatrix

end module CompositeTypeModule
