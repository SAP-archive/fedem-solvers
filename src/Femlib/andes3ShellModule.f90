!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module Andes3ShellModule

  !!============================================================================
  !! This module contains subroutines for calculating the element stiffness
  !! matrix of a 3-noded shell element based on the ANDES formulation.
  !!
  !! The code is developed by Bjorn Haugen (NTNU), April 2012.
  !! Adapted to FEDEM standard by Knut Morten Okstad, April 2016.
  !!============================================================================

  implicit none

  private

  public :: Andes3shell_stiffmat


contains

  subroutine StrainDispMatrix3ben (lumpType,alphaH,x,y,zeta,Bmatrix,ierr)

    !***************************************************************************
    ! Computes the strain-displacement matrix for a 3-noded bending element.
    ! The strains are computed using the projection rule.
    !
    ! Input arguments:
    ! lumpType : Order of normal rotation along edges
    ! alphaH   : Scaling factor for the higher-order part
    ! x, y     : Delta coordinates
    ! zeta     : Area coordinates of current integration point
    !
    ! Output arguments:
    ! Bmatrix  : Strain-displacement matrix with cartesian curvatures
    !            Nodal dof ordering: 3 x-rot., 3 y-rot., 3 z-disp.
    ! ierr     : Error flag
    !***************************************************************************

    use KindModule, only : dp

    integer , intent(in)  :: lumpType
    real(dp), intent(in)  :: alphaH, x(3,3), y(3,3), zeta(3)
    real(dp), intent(out) :: Bmatrix(3,9)
    integer , intent(out) :: ierr

    real(dp) :: a, Lt(3,9), Bh(3,9)

    ! Area of the triangle
    a = 0.5_dp*(x(2,1)*y(3,1) - x(3,1)*y(2,1))

    ! Compute the transpose of the lumping matrix
    call lump3ben (lumpType,x,y,Lt)

    ! Compute the deviatoric part of the strain-displacement matrix
    call bh3ben (x,y,zeta,Bh,ierr)
    if (ierr < 0) return

    ! Total strain-displacement matrix as the sum of the lumping matrix divided
    ! by the area, and the scaled higher-order strain-displacement matrix
    Bmatrix = Lt/a + Bh*alphaH

  end subroutine StrainDispMatrix3ben


  subroutine StrainDispMatrix3mem (alphaNR,alphaH,x,y,zeta,Bmatrix,ierr)

    !***************************************************************************
    ! Computes the strain-displacement matrix for a 3-noded membrane element.
    !
    ! Input arguments:
    ! alphaNR : Scaling factor for the normal rotation part
    ! alphaH  : Scaling factor for the higher-order part
    ! x, y    : Delta coordinates
    ! zeta    : Area coordinates of current integration point
    !
    ! Output arguments:
    ! Bmatrix : strain-displacement matrix with natural membrane strains
    !           Nodal dof ordering: 3 x-disp., 3 y-disp., 3 normal rotations
    ! ierr    : Error flag
    !***************************************************************************

    use KindModule, only : dp

    real(dp), intent(in)  :: alphaNR, alphaH, x(3,3), y(3,3), zeta(3)
    real(dp), intent(out) :: Bmatrix(3,9)
    integer , intent(out) :: ierr

    real(dp) :: a, Lt(3,9), Bh(3,9)

    ! Area of the triangle
    a = 0.5_dp*(x(2,1)*y(3,1) - x(3,1)*y(2,1))

    ! Compute the transpose of the lumping matrix
    call lump3mem (x,y,alphaNR,Lt)

    ! Compute the deviatoric part of the strain-displacement matrix
    call bh3mem (x,y,zeta,Bh,ierr)
    if (ierr < 0) return

    ! Total strain-displacement matrix as the sum of the lumping matrix divided
    ! by the area, and the scaled higher-order strain-displacement matrix
    Bmatrix = Lt/a + Bh*alphaH

  end subroutine StrainDispMatrix3mem


  subroutine lump3ben (ltype,x,y,lt)

    !***************************************************************************
    ! Computes the transposed lumping matrix for a 3-noded bending element.
    !
    ! Input arguments:
    ! ltype : = 1 Linear normal rotation along edges
    !             (i.e. C0 or Mindlin-Riesner like shape functions)
    !             Otherwise: Quadratic rotation variation
    !             (i.e. Beam-like hermitean shape functions)
    ! x, y  : Delta coordinates (x(i,j) == xi - xj)
    !
    ! Output arguments:
    ! lt    : Transposed lumping matrix
    !         Nodal dof ordering: 3 x-rot., 3 y-rot., 3 normal displacements
    !***************************************************************************

    use KindModule, only : dp

    integer , intent(in)  :: ltype
    real(dp), intent(in)  :: x(3,3), y(3,3)
    real(dp), intent(out) :: lt(3,9)

    integer  :: i, j
    real(dp) :: s(3,3), c(3,3), len

    select case (ltype)

    case (1) ! Linear normal rotation variation

       ! x - rotations
       lt(1,1) = x(3,2)*0.5_dp
       lt(1,2) = x(1,3)*0.5_dp
       lt(1,3) = x(2,1)*0.5_dp
       lt(2,1) = 0.0_dp
       lt(2,2) = 0.0_dp
       lt(2,3) = 0.0_dp
       lt(3,1) = y(2,3)*0.5_dp
       lt(3,2) = y(3,1)*0.5_dp
       lt(3,3) = y(1,2)*0.5_dp

       ! y - rotations
       lt(1,4) = 0.0_dp
       lt(1,5) = 0.0_dp
       lt(1,6) = 0.0_dp
       lt(2,4) = y(3,2)*0.5_dp
       lt(2,5) = y(1,3)*0.5_dp
       lt(2,6) = y(2,1)*0.5_dp
       lt(3,4) = x(2,3)*0.5_dp
       lt(3,5) = x(3,1)*0.5_dp
       lt(3,6) = x(1,2)*0.5_dp

       ! normal displacements
       lt(1,7) = 0.0_dp
       lt(1,8) = 0.0_dp
       lt(1,9) = 0.0_dp
       lt(2,7) = 0.0_dp
       lt(2,8) = 0.0_dp
       lt(2,9) = 0.0_dp
       lt(3,7) = 0.0_dp
       lt(3,8) = 0.0_dp
       lt(3,9) = 0.0_dp

    case default ! Quadratic normal rotation variation

       do i = 1, 3
          do j = 1, 3
             if (j /= i)  then
                len = sqrt(x(i,j)*x(i,j) + y(i,j)*y(i,j))
                c(i,j) = x(i,j)/len
                s(i,j) = y(i,j)/len
             end if
          end do
          s(i,i) = 0.0_dp
          c(i,i) = 0.0_dp
       end do

       ! x - rotations
       lt(1,1) =  (s(1,2)*s(1,2)*x(1,2) + s(3,1)*s(3,1)*x(3,1))*0.5_dp
       lt(1,2) =  (s(1,2)*s(1,2)*x(1,2) + s(2,3)*s(2,3)*x(2,3))*0.5_dp
       lt(1,3) =  (s(2,3)*s(2,3)*x(2,3) + s(3,1)*s(3,1)*x(3,1))*0.5_dp
       lt(2,1) =  (c(1,2)*c(1,2)*x(1,2) + c(3,1)*c(3,1)*x(3,1))*0.5_dp
       lt(2,2) =  (c(1,2)*c(1,2)*x(1,2) + c(2,3)*c(2,3)*x(2,3))*0.5_dp
       lt(2,3) =  (c(2,3)*c(2,3)*x(2,3) + c(3,1)*c(3,1)*x(3,1))*0.5_dp
       lt(3,1) = -(c(1,2)*c(1,2)*y(1,2) + c(3,1)*c(3,1)*y(3,1))
       lt(3,2) = -(c(1,2)*c(1,2)*y(1,2) + c(2,3)*c(2,3)*y(2,3))
       lt(3,3) = -(c(2,3)*c(2,3)*y(2,3) + c(3,1)*c(3,1)*y(3,1))

       ! y - rotations
       lt(1,4) =  (s(1,2)*s(1,2)*y(1,2) + s(3,1)*s(3,1)*y(3,1))*0.5_dp
       lt(1,5) =  (s(1,2)*s(1,2)*y(1,2) + s(2,3)*s(2,3)*y(2,3))*0.5_dp
       lt(1,6) =  (s(2,3)*s(2,3)*y(2,3) + s(3,1)*s(3,1)*y(3,1))*0.5_dp
       lt(2,4) =  (c(1,2)*c(1,2)*y(1,2) + c(3,1)*c(3,1)*y(3,1))*0.5_dp
       lt(2,5) =  (c(1,2)*c(1,2)*y(1,2) + c(2,3)*c(2,3)*y(2,3))*0.5_dp
       lt(2,6) =  (c(2,3)*c(2,3)*y(2,3) + c(3,1)*c(3,1)*y(3,1))*0.5_dp
       lt(3,4) = -(s(1,2)*s(1,2)*x(1,2) + s(3,1)*s(3,1)*x(3,1))
       lt(3,5) = -(s(1,2)*s(1,2)*x(1,2) + s(2,3)*s(2,3)*x(2,3))
       lt(3,6) = -(s(2,3)*s(2,3)*x(2,3) + s(3,1)*s(3,1)*x(3,1))

       ! normal displacements
       lt(1,7) = -c(1,2)*s(1,2) + c(3,1)*s(3,1)
       lt(1,8) = -c(2,3)*s(2,3) + c(1,2)*s(1,2)
       lt(1,9) = -c(3,1)*s(3,1) + c(2,3)*s(2,3)
       lt(2,7) = -lt(1,7)
       lt(2,8) = -lt(1,8)
       lt(2,9) = -lt(1,9)
       lt(3,7) = -(s(1,2)*s(1,2)-c(1,2)*c(1,2)) + (s(3,1)*s(3,1)-c(3,1)*c(3,1))
       lt(3,8) = -(s(2,3)*s(2,3)-c(2,3)*c(2,3)) + (s(1,2)*s(1,2)-c(1,2)*c(1,2))
       lt(3,9) = -(s(3,1)*s(3,1)-c(3,1)*c(3,1)) + (s(2,3)*s(2,3)-c(2,3)*c(2,3))

    end select

  end subroutine lump3ben


  subroutine bh3ben (x,y,zeta,Bben,ierr)

    !***************************************************************************
    ! Computes strain-displacement matrix for a 3-noded bending element.
    ! The strains are computed using the projection rule.
    !
    ! Input arguments:
    ! x, y : Delta coordinates
    ! zeta : Area coordinates of current integration point
    !
    ! Output arguments:
    ! Bben : Strain-displacement matrix with cartesian curvatures
    !        Nodal dof ordering: 3 x-rot., 3 y-rot., 3 z-disp.
    ! ierr : Error flag (ierr < 0: Edge -ierr has zero or negative length)
    !***************************************************************************

    use kindModule, only : dp, epsDiv0_p

    real(dp), intent(in)  :: x(3,3), y(3,3), zeta(3)
    real(dp), intent(out) :: Bben(3,9)
    integer , intent(out) :: ierr

    real(dp) :: len2, Tmat(3,3), Bbenc(3,6), Qmat(6,9), eta(3,3), zetabar(3)
    integer  :: i, j, k

    ! Compute side lenghts, heights and side lenght projections
    ! Compute side angle cosines and deviatoric part of the natural coordinates
    ierr = 0
    eta  = 0.0_dp
    do i = 1, 3
       j = 1 + mod(i,3)
       k = 1 + mod(j,3)
       len2 = x(k,j)*x(k,j) + y(k,j)*y(k,j)
       if (len2 <= epsDiv0_p) then
          ierr = -i
          return
       end if
       eta(i,j) = (x(i,k)*x(j,k) + y(i,k)*y(j,k))/len2
       eta(i,k) = (x(i,j)*x(k,j) + y(i,j)*y(k,j))/len2
       zetabar(i) = zeta(i) - 1.0_dp/3.0_dp
    end do

    ! Compute transformation from natural to cartesian strains
    call tran_bstrain (x,y,Tmat)

    ! Form natural coordinate b-matrix from hierarcical rotations
    Bben = 0.0_dp
    Bben(1,1) = zetabar(1) + eta(3,1)*zetabar(3)
    Bben(1,2) = zetabar(2) + eta(3,2)*zetabar(3)
    Bben(2,3) = zetabar(2) + eta(1,2)*zetabar(1)
    Bben(2,4) = zetabar(3) + eta(1,3)*zetabar(1)
    Bben(3,5) = zetabar(3) + eta(2,3)*zetabar(2)
    Bben(3,6) = zetabar(1) + eta(2,1)*zetabar(2)

    ! Transform to cartesian strains
    do i = 1, 3
       Bbenc(i,1) = Tmat(i,1)*Bben(1,1)
       Bbenc(i,2) = Tmat(i,1)*Bben(1,2)
       Bbenc(i,3) = Tmat(i,2)*Bben(2,3)
       Bbenc(i,4) = Tmat(i,2)*Bben(2,4)
       Bbenc(i,5) = Tmat(i,3)*Bben(3,5)
       Bbenc(i,6) = Tmat(i,3)*Bben(3,6)
    end do

    ! Form transformation matrix to visible DOF
    qmat(1,1) = -4.0_dp*y(2,1)
    qmat(2,1) =  2.0_dp*y(2,1)
    qmat(3,1) =  0.0_dp
    qmat(4,1) =  0.0_dp
    qmat(5,1) = -2.0_dp*y(1,3)
    qmat(6,1) =  4.0_dp*y(1,3)

    qmat(1,2) = -2.0_dp*y(2,1)
    qmat(2,2) =  4.0_dp*y(2,1)
    qmat(3,2) = -4.0_dp*y(3,2)
    qmat(4,2) =  2.0_dp*y(3,2)
    qmat(5,2) =  0.0_dp
    qmat(6,2) =  0.0_dp

    qmat(1,3) =  0.0_dp
    qmat(2,3) =  0.0_dp
    qmat(3,3) = -2.0_dp*y(3,2)
    qmat(4,3) =  4.0_dp*y(3,2)
    qmat(5,3) = -4.0_dp*y(1,3)
    qmat(6,3) =  2.0_dp*y(1,3)

    qmat(1,4) =  4.0_dp*x(2,1)
    qmat(2,4) = -2.0_dp*x(2,1)
    qmat(3,4) =  0.0_dp
    qmat(4,4) =  0.0_dp
    qmat(5,4) =  2.0_dp*x(1,3)
    qmat(6,4) = -4.0_dp*x(1,3)

    qmat(1,5) =  2.0_dp*x(2,1)
    qmat(2,5) = -4.0_dp*x(2,1)
    qmat(3,5) =  4.0_dp*x(3,2)
    qmat(4,5) = -2.0_dp*x(3,2)
    qmat(5,5) =  0.0_dp
    qmat(6,5) =  0.0_dp

    qmat(1,6) =  0.0_dp
    qmat(2,6) =  0.0_dp
    qmat(3,6) =  2.0_dp*x(3,2)
    qmat(4,6) = -4.0_dp*x(3,2)
    qmat(5,6) =  4.0_dp*x(1,3)
    qmat(6,6) = -2.0_dp*x(1,3)

    qmat(1,7) = -6.0_dp
    qmat(2,7) =  6.0_dp
    qmat(3,7) =  0.0_dp
    qmat(4,7) =  0.0_dp
    qmat(5,7) =  6.0_dp
    qmat(6,7) = -6.0_dp

    qmat(1,8) =  6.0_dp
    qmat(2,8) = -6.0_dp
    qmat(3,8) = -6.0_dp
    qmat(4,8) =  6.0_dp
    qmat(5,8) =  0.0_dp
    qmat(6,8) =  0.0_dp

    qmat(1,9) =  0.0_dp
    qmat(2,9) =  0.0_dp
    qmat(3,9) =  6.0_dp
    qmat(4,9) = -6.0_dp
    qmat(5,9) = -6.0_dp
    qmat(6,9) =  6.0_dp

    ! Transform to visible DOF
    Bben = matmul(Bbenc,qmat)

  end subroutine bh3ben


  subroutine lump3mem (x,y,alpha,lt)

    !***************************************************************************
    ! Computes the transposed lumping matrix for a 3-noded membrane element.
    !
    ! Input arguments:
    ! x, y  : Delta coordinates
    ! alpha : Scaling factor for rotational contribution
    !
    ! Output arguments:
    ! lt    : Transposed lumping matrix
    !         Nodal dof ordering: 3 x-disp., 3 y-disp., 3 normal rotations
    !***************************************************************************

    use KindModule, only : dp

    real(dp), intent(in)  :: x(3,3), y(3,3), alpha
    real(dp), intent(out) :: lt(3,9)

    ! x - displacements
    lt(1,1) = y(2,3)*0.5_dp
    lt(1,2) = y(3,1)*0.5_dp
    lt(1,3) = y(1,2)*0.5_dp
    lt(2,1) = 0.0_dp
    lt(2,2) = 0.0_dp
    lt(2,3) = 0.0_dp
    lt(3,1) = x(3,2)*0.5_dp
    lt(3,2) = x(1,3)*0.5_dp
    lt(3,3) = x(2,1)*0.5_dp

    ! y - displacements
    lt(1,4) = 0.0_dp
    lt(1,5) = 0.0_dp
    lt(1,6) = 0.0_dp
    lt(2,4) = x(3,2)*0.5_dp
    lt(2,5) = x(1,3)*0.5_dp
    lt(2,6) = x(2,1)*0.5_dp
    lt(3,4) = y(2,3)*0.5_dp
    lt(3,5) = y(3,1)*0.5_dp
    lt(3,6) = y(1,2)*0.5_dp

    ! normal rotations
    lt(1,7) =  y(2,3)*(y(1,3)-y(2,1))*alpha/12.0_dp
    lt(2,7) =  x(3,2)*(x(3,1)-x(1,2))*alpha/12.0_dp
    lt(3,7) = (x(3,1)*y(1,3)-x(1,2)*y(2,1))*alpha/6.0_dp
    lt(1,8) =  y(3,1)*(y(2,1)-y(3,2))*alpha/12.0_dp
    lt(2,8) =  x(1,3)*(x(1,2)-x(2,3))*alpha/12.0_dp
    lt(3,8) = (x(1,2)*y(2,1)-x(2,3)*y(3,2))*alpha/6.0_dp
    lt(1,9) =  y(1,2)*(y(3,2)-y(1,3))*alpha/12.0_dp
    lt(2,9) =  x(2,1)*(x(2,3)-x(3,1))*alpha/12.0_dp
    lt(3,9) = (x(2,3)*y(3,2)-x(3,1)*y(1,3))*alpha/6.0_dp

  end subroutine lump3mem


  subroutine bh3mem (x,y,zeta,Bmem,ierr)

    !***************************************************************************
    ! Computes the strain-displacement matrix for a 3-noded membrane element.
    !
    ! Input arguments:
    ! x, y : Delta coordinates
    ! zeta : Area coordinates of current integration point
    !
    ! Output arguments:
    ! Bmem : Strain-displacement matrix with natural membrane strains
    !        Nodal dof ordering: 3 x-disp., 3 y-disp., 3 normal rotations
    ! ierr : Error flag (ierr < 0: Edge -ierr has zero or negative length)
    !***************************************************************************

    use KindModule, only : dp

    real(dp), intent(in)  :: x(3,3), y(3,3), zeta(3)
    real(dp), intent(out) :: Bmem(3,9)
    integer , intent(out) :: ierr

    real(dp) :: length(3), area2, ht(3,9), Tmat(3,3), Bcar(3,3)
    real(dp) :: xfac12, xfac23, xfac31
    integer  :: i

    ! Compute side lengths and strain transformation matrix
    call tran_mstrain (x,y,length,Tmat,ierr)
    if (ierr < 0) return

    ! Form natural coordinate B-matrix from hierarcical rotations
    area2  = x(2,1)*y(3,1) - x(3,1)*y(2,1)
    xfac12 = area2 / (3.0_dp*length(3)*length(3))
    xfac23 = area2 / (3.0_dp*length(1)*length(1))
    xfac31 = area2 / (3.0_dp*length(2)*length(2))
    Bmem(1,1) = xfac12*(   zeta(1) - 2*zeta(2) +   zeta(3))
    Bmem(1,2) = xfac12*( 2*zeta(1) -   zeta(2) -   zeta(3))
    Bmem(1,3) = xfac12*(   zeta(1) -   zeta(2)            )
    Bmem(2,1) = xfac23*(               zeta(2) -   zeta(3))
    Bmem(2,2) = xfac23*(   zeta(1) +   zeta(2) - 2*zeta(3))
    Bmem(2,3) = xfac23*(  -zeta(1) + 2*zeta(2) -   zeta(3))
    Bmem(3,1) = xfac31*(  -zeta(1) -   zeta(2) + 2*zeta(3))
    Bmem(3,2) = xfac31*(  -zeta(1)             +   zeta(3))
    Bmem(3,3) = xfac31*(-2*zeta(1) +   zeta(2) +   zeta(3))

    ! Transform to cartesian strains
    Bcar = matmul(Tmat,Bmem(:,1:3))

    ! Form transformation matrix to visible DOF
    do i = 1, 3
       ht(i,1)   = x(3,2)*0.5_dp/area2
       ht(i,2)   = x(1,3)*0.5_dp/area2
       ht(i,3)   = x(2,1)*0.5_dp/area2
       ht(i,4)   = y(3,2)*0.5_dp/area2
       ht(i,5)   = y(1,3)*0.5_dp/area2
       ht(i,6)   = y(2,1)*0.5_dp/area2
       ht(i,7:9) = 0.0_dp
       ht(i,6+i) = 1.0_dp
    end do

    ! Transform to global DOFs
    Bmem = matmul(Bcar,ht)

  end subroutine bh3mem


  subroutine tran_bstrain (x,y,Tmat)

    !***************************************************************************
    ! Computes transformation matrix from Natural to Cartestian curvatures.
    !
    ! Input arguments:
    ! x, y : Delta coordinates (x(i,j) = xi - xj)
    !
    ! Output arguments:
    ! Tmat : Transformation matrix from from natural to cartesian curvatures
    !***************************************************************************

    use KindModule, only : dp

    real(dp), intent(in)  :: x(3,3), y(3,3)
    real(dp), intent(out) :: Tmat(3,3)

    real(dp) :: area2

    Tmat(1,1) =   y(2,3)*y(1,3)
    Tmat(1,2) =   y(3,1)*y(2,1)
    Tmat(1,3) =   y(1,2)*y(3,2)
    Tmat(2,1) =   x(2,3)*x(1,3)
    Tmat(2,2) =   x(3,1)*x(2,1)
    Tmat(2,3) =   x(1,2)*x(3,2)
    Tmat(3,1) = -(y(2,3)*x(1,3) + x(2,3)*y(1,3))
    Tmat(3,2) = -(y(3,1)*x(2,1) + x(3,1)*y(2,1))
    Tmat(3,3) = -(y(1,2)*x(3,2) + x(1,2)*y(3,2))

    ! Area of triangle times 2.0
    area2 = x(2,1)*y(3,1) - x(3,1)*y(2,1)

    Tmat = Tmat/(area2*area2)

  end subroutine tran_bstrain


  subroutine tran_mstrain (x,y,length,Tmat,ierr)

    !***************************************************************************
    ! Computes transformation matrix from Natural to Cartesian membrane strains.
    !
    ! Input arguments:
    ! x, y   : Delta coordinates (x(i,j) = xi - xj)
    !
    ! Output arguments:
    ! length : Side lengths
    ! Tmat   : Transformation matrix from natural- to xy-strains
    ! ierr   : Error flag
    !***************************************************************************

    use kindModule       , only : dp, epsDiv0_p
    use manipMatrixModule, only : invert33
    use reportErrorModule, only : getErrorFile

    real(dp), intent(in)  :: x(3,3), y(3,3)
    real(dp), intent(out) :: length(3), tmat(3,3)
    integer , intent(out) :: ierr

    integer  :: i, j, k
    real(dp) :: len2, Tinv(3,3)

    ierr = 0
    do i = 1, 3
       j = 1 + mod(i,3)
       k = 1 + mod(j,3)
       len2 = x(k,j)*x(k,j) + y(k,j)*y(k,j)
       if (len2 <= epsDiv0_p) then
          ierr = -i
          return
       end if
       length(i) = sqrt(len2)
       Tinv(j,1) = x(k,j)*x(k,j)/len2
       Tinv(j,2) = y(k,j)*y(k,j)/len2
       Tinv(j,3) = x(k,j)*y(k,j)/len2
    end do

    Tmat = invert33(Tinv,getErrorFile(),ierr)
    ierr = ierr*10

  end subroutine tran_mstrain


  subroutine Bmatrix (lumpType,alphaH,alphaNR,dx,dy,zeta,Bmat,ierr)

    !***************************************************************************
    ! Creates the final B-matrix, based on the membrane and bending parts.
    !
    ! Input arguments:
    ! lumpType : Order of normal rotation along edges
    ! alphaH   : Scaling factor for the higher-order part
    ! alphaNR  : Scaling factor for the normal rotation part
    ! dx, dy   : Delta coordinates
    ! zeta     : Area coordinates of current integration point
    !
    ! Output arguments:
    ! Bmat     : Final B-matrix with dof ordering
    !            u1, v1, w1, x-rot1, y-rot1, z-rot1, u2, v2, w2 ... z-rot3
    ! ierr     : Error flag
    !***************************************************************************

    use KindModule, only : dp

    integer , intent(in)  :: lumpType
    real(dp), intent(in)  :: dx(3,3), dy(3,3), alphaH, alphaNR, zeta(3)
    real(dp), intent(out) :: Bmat(6,18)
    integer , intent(out) :: ierr

    integer  :: i
    real(dp) :: Bb(3,9), Bm(3,9)

    Bb = 0.0_dp ! Should bot be needed, but gcc -Wall -Wextra complains
    Bm = 0.0_dp ! Should bot be needed, but gcc -Wall -Wextra complains

    ! Create Bb
    call StrainDispMatrix3ben (lumpType,alphaH,dx,dy,zeta,Bb,ierr)
    if (ierr < 0) return

    ! Create Bm
    call StrainDispMatrix3mem (alphaNR,alphaH,dx,dy,zeta,Bm,ierr)
    if (ierr < 0) return

    ! Create the final B-matrix(6,18)
    Bmat = 0.0_dp
    do i = 1, 3
       Bmat(i  ,1 ) = Bm(i, 1)
       Bmat(i  ,2 ) = Bm(i, 4)
       Bmat(i+3,3 ) = Bb(i, 7)
       Bmat(i+3,4 ) = Bb(i, 1)
       Bmat(i+3,5 ) = Bb(i, 4)
       Bmat(i  ,6 ) = Bm(i, 7)
       Bmat(i  ,7 ) = Bm(i, 2)
       Bmat(i  ,8 ) = Bm(i, 5)
       Bmat(i+3,9 ) = Bb(i, 8)
       Bmat(i+3,10) = Bb(i, 2)
       Bmat(i+3,11) = Bb(i, 5)
       Bmat(i  ,12) = Bm(i, 8)
       Bmat(i  ,13) = Bm(i, 3)
       Bmat(i  ,14) = Bm(i, 6)
       Bmat(i+3,15) = Bb(i, 9)
       Bmat(i+3,16) = Bb(i, 3)
       Bmat(i+3,17) = Bb(i, 6)
       Bmat(i  ,18) = Bm(i, 9)
    end do

  end subroutine Bmatrix


  subroutine Andes3shell_stiffmat (xl,yl,Cmat,alpha,beta,lumpType,Kmat,ierr)

    !***************************************************************************
    ! Creates the stiffness matrix for a 3-noded shell element.
    !***************************************************************************

    use kindModule       , only : dp, epsDiv0_p
    use reportErrorModule, only : reportError, error_p, debugFileOnly_p

    real(dp), intent(in)  :: xl(3), yl(3), Cmat(6,6), alpha, beta
    real(dp), intent(out) :: Kmat(18,18)
    integer , intent(in)  :: lumpType
    integer , intent(out) :: ierr

    integer, parameter :: nmult = 3
    integer            :: i, j
    real(dp)           :: Bmat(6,18), lstzeta(3,nmult), aw, dx(3,3), dy(3,3)
    character(len=64)  :: errMsg

    ! Triangle quadrature points
    lstzeta = 0.5_dp
    lstzeta(1,1) = 0.0_dp
    lstzeta(2,2) = 0.0_dp
    lstzeta(3,3) = 0.0_dp

    ! Delta values for the nodal coordinates
    do i = 1, 3
       do j = 1, 3
          dx(i,j) = xl(i) - xl(j)
          dy(i,j) = yl(i) - yl(j)
       end do
    end do

    ! Element area scaled by integration point weight (1/3)
    aw = (dx(2,1)*dy(3,1) - dx(3,1)*dy(2,1))/6.0_dp
    if (aw < epsDiv0_p) then
       write(errMsg,"('Degenerated element, non-positive area',1PE13.5)") aw
       call reportError (error_p,errMsg,addString='Andes3shell_stiffmat')
       ierr = -9
       return
    end if

    ! Summation of the stiffness matrix (B_i^t*C*B_i)
    ierr = 0
    Kmat = 0.0_dp
    do i = 1, nmult
       call Bmatrix (lumpType,beta,alpha,dx,dy,lstzeta(:,i),Bmat,ierr)
       if (ierr < 0) exit
       Kmat = Kmat + matmul(transpose(Bmat),matmul(Cmat*aw,Bmat))
    end do

    if (ierr < -3) then
       call reportError (debugFileOnly_p,'Andes3shell_stiffmat')
    else if (ierr < 0) then
       write(errMsg,"('Degenerated element, non-positive edge',I2)") -ierr
       call reportError (error_p,errMsg,addString='Andes3shell_stiffmat')
    end if

  end subroutine Andes3shell_stiffmat

end module Andes3ShellModule
