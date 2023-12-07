!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module ipri6Module

  !!============================================================================
  !! This module provides routines to generate the element stiffness- and mass-
  !! matrix as well as consistent load vectors, and to compute element stresses
  !! and strains for a 6-noded isoparametric prism element with 18 DOFs.
  !!
  !! Coded by: Hans Magnus Haug, 21 Apr 1999
  !! Revised : Trond Arne Svidal, 8 Oct 1999 / 1.0
  !!           Converted to f90.
  !!           Fixed mass matrix (2 -> 6 Gauss points).
  !!           Return strain and stress for all nodes in the element.
  !!============================================================================

  implicit none

  private :: ipri6bmat, pdvn


contains

  subroutine ipri6stiff (x, y, z, cmat, stiffm, iel, iw, ipsw, ierr)

    !!==========================================================================
    !! This subroutine generates the stiffness matrix for
    !! the 6-noded isoparametric prism element.
    !! Integration is done by use of Gauss quadrature with
    !! one sample point in triangle plane and two in zeta-direction.
    !!
    !! Input:
    !!    x      - array of x-coordinates for all six nodes
    !!    y      - array of y-coordinates
    !!    z      - array of z-coordinates
    !!    cmat   - constitutive matrix
    !!    iel    - element number
    !!    iw     - print unit
    !!    ipsw   - print switch
    !!
    !! Output:
    !!    stiffm - element stiffness matrix for current element,
    !!             dof organized as vxi,vyi,vzi for each node
    !!    ierr   - error flag, negative value signals error
    !!
    !! Coded by : Hans Magnus Haug
    !! Revised  : Trond Arne Svidal,  8 Oct 1999 / 1.0
    !!==========================================================================

    use KindModule, only : dp

    integer , intent(in)  :: iel, iw, ipsw
    real(dp), intent(in)  :: x(6), y(6), z(6), cmat(6,6)
    real(dp), intent(out) :: stiffm(18,18)
    integer , intent(out) :: ierr

    !! Local variables
    integer  :: i, j, k
    real(dp) :: bmat(6,18), cb(6), spntxi(3), spntzeta(2), detJ

    !! --- Logic section ---

    if (ipsw >= 1) then
       write(iw,*) '===> Entering ipri6stiff with input ==='
       write(iw,9000) 'element number:', iel
       write(iw,9001) 'nodal x-coord: ', x
       write(iw,9001) 'nodal z-coord: ', y
       write(iw,9001) 'nodal y-coord: ', z
       write(iw,*) 'constitutive matrix:'
       write(iw,9002) (cmat(i,:),i=1,6)
       write(iw,*) '======================================='
    end if

    !! Initialization
    ierr   = 0
    stiffm = 0.0_dp

    !! Integration point coordinates
    spntxi      =  1.0_dp/3.0_dp
    spntzeta(1) = -1.0_dp/sqrt(3.0_dp)
    spntzeta(2) = -spntzeta(1)

    !! Start the Gauss quadrature loop, k=SUM(1,2)[1*BT*C*B*detJ]
    do k = 1, 2

       call ipri6bmat (x,y,z, spntxi,spntzeta(k), bmat, detJ, iel, iw,ipsw,ierr)
       if (ierr < 0) then
          write(iw,*) '*** Error return from ipri6stiff'
          write(iw,*) '    ierr = ', ierr
          return
       end if

       do j = 1, 18
          cb = matmul(cmat,bmat(:,j))
          do i = 1, 18
             stiffm(i,j) = stiffm(i,j) + dot_product(bmat(:,i),cb)*detJ
          end do
       end do

    end do

    if (ipsw >= 2) then
       write(iw,*) '<=== Leaving ipri6stiff with output ==='
       write(iw,*) 'Stiffness matrix:'
       write(iw,9003) (stiffm(i,:),i=1,18)
       write(iw,*) '======================================='
    end if

9000 format(1X,A,I6)
9001 format(1X,A,1P,4E12.4)
9002 format(1P,6E12.4)
9003 format(1P,18E12.4)

  end subroutine ipri6stiff


  subroutine ipri6mass (x, y, z, rho, massm, iel, iw, ipsw, ierr)

    !!==========================================================================
    !! This subroutine generates the consistent mass matrix for
    !! the 6-noded isoparametric prism element.
    !! Integration is done by use of Gauss quadrature with
    !! three sample points in triangle plane and two in zeta-direction.
    !!
    !! Input:
    !!    x      - array of x-coordinates for all six nodes
    !!    y      - array of y-coordinates
    !!    z      - array of z-coordinates
    !!    rho    - material density
    !!    iel    - element number
    !!    iw     - print unit
    !!    ipsw   - print switch
    !!
    !! Output:
    !!    massm  - element mass matrix
    !!             dof organized as vxi,vyi,vzi for each node
    !!    ierr   - error flag, negative value signals error
    !!
    !! Coded by: Hans Magnus Haug
    !! Revised : Trond Arne Svidal,  8 Oct 1999 / 1.0
    !!==========================================================================

    use KindModule, only : dp

    integer , intent(in)  :: iel, iw, ipsw
    real(dp), intent(in)  :: x(6), y(6), z(6), rho
    real(dp), intent(out) :: massm(18,18)
    integer , intent(out) :: ierr

    !! Local variables
    integer  :: i, numZ, numX, numN, numR, numC
    real(dp) :: spntxi(3,3), spntzeta(2), n(6)
    real(dp) :: dnxi1(6), dnxi2(6), dnzeta(6), inja(3,3), detJ, nMat(3,18)

    !! --- Logic section ---

    if (ipsw >= 1) then
       write(iw,*) '===> Entering ipri6mass with input ==='
       write(iw,9000) 'element number:', iel
       write(iw,9001) 'nodal x-coord: ', x
       write(iw,9001) 'nodal z-coord: ', y
       write(iw,9001) 'nodal y-coord: ', z
       write(iw,9001) ' mass density: ', rho
       write(iw,*) '======================================'
    end if

    !! Initialization
    ierr   = 0
    massm  = 0.0_dp
    nMat   = 0.0_dp

    !! Integration point coordinates
    spntxi      = 0.5_dp
    spntxi(1,1) = 0.0_dp
    spntxi(2,2) = 0.0_dp
    spntxi(3,3) = 0.0_dp
    spntzeta(1) = -1.0_dp / sqrt(3.0_dp)
    spntzeta(2) = -spntzeta(1)

    !! Start the Gauss quadrature loop, m=SUM(1,6)[1*1/3*NT*N*detJ]
    do numZ = 1, 2
       do numX = 1, 3

          !! Compute partial derivatives of shape functions
          call pdvn (spntxi(:,numX), spntzeta(numZ), dnxi1, dnxi2, dnzeta)

          !! Compute Jacobi determinant, detJ
          call JACI31 (inja, detJ, dnxi1, dnxi2, dnzeta, &
               &       x, y, z, iel, 6, iw, ipsw, ierr)
          if (ierr < 0) then
             write(iw,*) '*** Error return from ipri6mass'
             write(iw,*) '    ierr = ', ierr
             return
          end if

          !! === Calculate shape function matrix, nMat
          !!
          !!            |n(1) 0   0  n(2) 0   0   ...  n(6) 0   0  |
          !!      nMat =| 0  n(1) 0   0  n(2) 0   ...   0  n(6) 0  |
          !!            | 0   0  n(1) 0   0  n(2) ...   0   0  n(6)|

          n(1:3) = spntxi(:,numX) * (1.0_dp - spntzeta(numZ))*0.5_dp
          n(4:6) = spntxi(:,numX) * (1.0_dp + spntzeta(numZ))*0.5_dp

          do numN = 0, 5
             do numR = 1, 3
                numC = 3*numN + numR
                nMat(numR,numC) = n(numN + 1)
             end do
          end do

          !! Integrate mass matrix
          massm = massm + matmul(transpose(nMat),nMat) * rho*detJ/6.0_dp

       end do
    end do

    if (ipsw >= 2) then
       write(iw,*) '<=== Leaving ipri6mass with output ==='
       write(iw,*) 'Mass matrix:'
       write(iw,9002) (massm(i,:),i=1,18)
       write(iw,*) '======================================'
    end if

9000 format(1X,A,I6)
9001 format(1X,A,1P,4E12.3)
9002 format(1P,18E12.3)

  end subroutine ipri6mass


  subroutine ipri6pload (x, y, z, p, iface, fvec, iw, ipsw, ierr)

    !!==========================================================================
    !! This subroutine generates the consistent force vector for
    !! the 6-noded isoparametric prism element, due to a pressure load.
    !!
    !! Input:
    !!    x      - array of x-coordinates for all six nodes
    !!    y      - array of y-coordinates
    !!    z      - array of z-coordinates
    !!    p      - surface pressure intensities at three face nodes
    !!    iface  - index of the pressurized face, [1,5]
    !!    iw     - print unit
    !!    ipsw   - print switch
    !!
    !! Output:
    !!    fvec   - element force vector, organized as fxi,fyi,fzi for each node
    !!    ierr   - error flag, negative value signals error
    !!
    !! Coded by: Knut Morten Okstad
    !! Date/ver: 1 Apr 2008 / 1.0
    !!==========================================================================

    use KindModule       , only : dp
    use ManipMatrixModule, only : cross_product

    integer , intent(in)  :: iface, iw, ipsw
    real(dp), intent(in)  :: x(6), y(6), z(6), p(3,4)
    real(dp), intent(out) :: fvec(3,6)
    integer , intent(out) :: ierr

    !! Local variables
    integer  :: i, j, i1, i2, i3, i4
    real(dp) :: s12(3), s13(3), vn(3), pw(3), xi, eta, dS, N1, N2, N3, N4

    !! Four-point triangle quadrature (see Hughes 1987, Table 3.I.1, p. 173).
    real(dp), parameter :: w4_p = 1.5625_dp/3.0_dp
    real(dp), parameter :: w(4) = (/    -0.5625_dp,  w4_p,   w4_p,   w4_p  /)
    real(dp), parameter :: r(4) = (/ 1.0_dp/3.0_dp, 0.6_dp, 0.2_dp, 0.2_dp /)
    real(dp), parameter :: s(4) = (/ 1.0_dp/3.0_dp, 0.2_dp, 0.6_dp, 0.2_dp /)

    !! --- Logic section ---

    if (ipsw >= 1) then
       write(iw,*) '===> Entering ipri6pload with input ==='
       write(iw,9001) 'nodal x-coord: ', x
       write(iw,9001) 'nodal z-coord: ', y
       write(iw,9001) 'nodal y-coord: ', z
       write(iw,9002) 'face index = ', iface
       write(iw,*) 'face load intensities:'
       write(iw,9003) p
       write(iw,*) '======================================='
    end if

    ierr = 0
    fvec = 0.0_dp

    select case (iface)
    case (1); i1 = 1; i2 = 2; i3 = 5; i4 = 4
    case (2); i1 = 2; i2 = 3; i3 = 6; i4 = 5
    case (3); i1 = 1; i2 = 4; i3 = 6; i4 = 3
    case (4); i1 = 1; i2 = 3; i3 = 2
    case (5); i1 = 4; i2 = 5; i3 = 6
    case default
       ierr = -1
       write(iw,*) '*** Error return from ipri6pload'
       write(iw,*) '    invalid face index = ', iface
       return
    end select

    if (iface > 3) then ! triangular face

       !! Form side vector from node 1 to 2
       s12(1) = x(i2) - x(i1)
       s12(2) = y(i2) - y(i1)
       s12(3) = z(i2) - z(i1)

       !! Form side vector from node 1 to 3
       s13(1) = x(i3) - x(i1)
       s13(2) = y(i3) - y(i1)
       s13(3) = z(i3) - z(i1)

       !! Find the face area, dS = 0.5 * || s12 x s13 ||
       vn = cross_product(s12,s13)
       dS = 0.5_dp*sqrt(dot_product(vn,vn))

       !! Numerical integration loop
       do i = 1, 4
          N1 = r(i)
          N2 = s(i)
          N3 = 1.0_dp - N1 - N2
          pw = (p(:,1)*N1 + p(:,2)*N2 + p(:,3)*N3) * w(i)*dS
          fvec(:,i1) = fvec(:,i1) + pw*N1
          fvec(:,i2) = fvec(:,i2) + pw*N2
          fvec(:,i3) = fvec(:,i3) + pw*N3
       end do

    else ! quadrilateral face

       !! Numerical integration loop (2x2 Gauss)
       do i = -1, 1, 2
          do j = -1, 1, 2
             xi  = real(i,dp)/sqrt(3.0_dp)
             eta = real(j,dp)/sqrt(3.0_dp)

             !! Find the surface dilatation, dS = || s12 x s13 ||, where
             !! s12 and s13 are the two tangent vectors {X},xi and {X},eta
             N1 = -0.25_dp * (1.0_dp - eta)
             N2 =  0.25_dp * (1.0_dp - eta)
             N3 =  0.25_dp * (1.0_dp + eta)
             N4 = -0.25_dp * (1.0_dp + eta)
             s12(1) = N1*x(i1) + N2*x(i2) + N3*x(i3) + N4*x(i4)
             s12(2) = N1*y(i1) + N2*y(i2) + N3*y(i3) + N4*y(i4)
             s12(3) = N1*z(i1) + N2*z(i2) + N3*z(i3) + N4*z(i4)
             N1 = -0.25_dp * (1.0_dp - xi)
             N2 = -0.25_dp * (1.0_dp + xi)
             N3 =  0.25_dp * (1.0_dp + xi)
             N4 =  0.25_dp * (1.0_dp - xi)
             s13(1) = N1*x(i1) + N2*x(i2) + N3*x(i3) + N4*x(i4)
             s13(2) = N1*y(i1) + N2*y(i2) + N3*y(i3) + N4*y(i4)
             s13(3) = N1*z(i1) + N2*z(i2) + N3*z(i3) + N4*z(i4)
             vn = cross_product(s12,s13)
             dS = sqrt(dot_product(vn,vn))

             !! Evaluate the surface pressure at current integration point
             N1 = 0.25_dp * (1.0_dp - xi)*(1.0_dp - eta)
             N2 = 0.25_dp * (1.0_dp + xi)*(1.0_dp - eta)
             N3 = 0.25_dp * (1.0_dp + xi)*(1.0_dp + eta)
             N4 = 0.25_dp * (1.0_dp - xi)*(1.0_dp + eta)
             pw = (p(:,1)*N1 + p(:,2)*N2 + p(:,3)*N3 + p(:,4)*N4) * dS

             !! Integrate the nodal force vector
             fvec(:,i1) = fvec(:,i1) + pw*N1
             fvec(:,i2) = fvec(:,i2) + pw*N2
             fvec(:,i3) = fvec(:,i3) + pw*N3
             fvec(:,i4) = fvec(:,i4) + pw*N4
          end do
       end do

    end if

    if (ipsw >= 2) then
       write(iw,*) '<=== Leaving ipri6pload with output ==='
       write(iw,*) 'Force vector:'
       write(iw,9003) fvec
       write(iw,*) '======================================='
    end if

9001 format(1X,A,1P,6E12.4)
9002 format(1X,A,I3)
9003 format(1P,3E12.4)

  end subroutine ipri6pload


  subroutine ipri6bmat (x, y, z, xi, zeta, bmat, detjac, iel, iw, ipsw, ierr)

    !!==========================================================================
    !! This subroutine generates the strain-displacement matrix for
    !! the 6-noded isoparametric prism element.
    !!
    !! Input:
    !!    x      - array of x-coordinates for all six nodes
    !!    y      - array of y-coordinates
    !!    z      - array of z-coordinates
    !!    xi     - array of area-coordinates. Dimension(3)
    !!    zeta   - interpolation point over element height
    !!    iel    - element number
    !!    iw     - print unit
    !!    ipsw   - print switch
    !!
    !! Output:
    !!    detjac - Jacobi determinant
    !!    bmat   - element strain-displacement matrix
    !!             dof ordering: vxi,vyi,vzi for each node
    !!    ierr   - error flag, negative value signals error
    !!
    !! Coded by: Hans Magnus Haug
    !! Revised : Trond Arne Svidal,  8 Oct 1999 / 1.0
    !!=================================================================

    use KindModule, only : dp

    integer , intent(in)  :: iel, iw, ipsw
    real(dp), intent(in)  :: x(6), y(6), z(6), xi(3), zeta
    real(dp), intent(out) :: bmat(6,18), detjac
    integer , intent(out) :: ierr

    !! Local variables
    integer  :: i, inod, xpos, ypos, zpos
    real(dp) :: a(6), b(6), c(6), inja(3,3)

    !! --- Logic section ---

    if (ipsw >= 3) then
       write(iw,*) '===> Entering ipri6bmat with input ==='
       write(iw,9000) 'element number:', iel
       write(iw,9001) 'nodal x-coord: ', x
       write(iw,9001) 'nodal z-coord: ', y
       write(iw,9001) 'nodal y-coord: ', z
       write(iw,9001) '  area coords: ', xi
       write(iw,9001) ' height coord: ', zeta
       write(iw,*) '======================================'
    end if

    !! Initialization
    ierr = 0
    bmat = 0.0_dp

    !! Collect derivatives of interpolation functions w.r.t xi1, xi2 and zeta
    call pdvn (xi, zeta, a, b, c)

    !! Get Jacobi inverse and determinant
    call JACI31 (inja, detjac, a, b, c, x, y, z, iel, 6, iw, ipsw, ierr)
    if (ierr < 0) then
       write(iw,*) '*** Error return from ipri6bmat'
       write(iw,*) '    ierr = ', ierr
       return
    end if

    !! Build b-matrix
    do inod = 1, 6
       xpos = (inod-1) * 3 + 1
       ypos = xpos + 1
       zpos = xpos + 2

       bmat(1,xpos) = a(inod)*inja(1,1) + b(inod)*inja(1,2) + c(inod)*inja(1,3)
       bmat(2,ypos) = a(inod)*inja(2,1) + b(inod)*inja(2,2) + c(inod)*inja(2,3)
       bmat(3,zpos) = a(inod)*inja(3,1) + b(inod)*inja(3,2) + c(inod)*inja(3,3)

       bmat(4,xpos) = bmat(2,ypos)
       bmat(4,ypos) = bmat(1,xpos)
       bmat(5,ypos) = bmat(3,zpos)
       bmat(5,zpos) = bmat(2,ypos)
       bmat(6,xpos) = bmat(3,zpos)
       bmat(6,zpos) = bmat(1,xpos)
    end do

    if (ipsw >= 4) then
       write(iw,*) '<=== Leaving ipri6bmat with output ==='
       write(iw,*) 'DetJac = ', detjac
       write(iw,*) 'Strain-displacement matrix:'
       write(iw,9002) (bmat(i,:),i=1,6)
       write(iw,*) '======================================'
    end if

9000 format(1X,A,I6)
9001 format(1X,A,1P,6E12.4)
9002 format(1P,18E12.4)

  end subroutine ipri6bmat


  subroutine pdvn (xi, zeta, dnxi1, dnxi2, dnzeta)

    !!==========================================================================
    !! Computes partial derivatives of shape functions w.r.t.
    !! natural coordinates xi and zeta.
    !!
    !! Input:
    !!    xi     - areal coordinates, dimension(3)
    !!    zeta   - height coordinate
    !!
    !! Output:
    !!    dnxi1  - partial derivates w.r.t. xi(1)
    !!    dnxi2  - partial derivates w.r.t. xi(2)
    !!    dnzeta - partial derivates w.r.t. zeta
    !!
    !! Coded by: Hans Magnus Haug
    !! Revised : Trond Arne Svidal,  8 Oct 1999 / 1.0
    !!==========================================================================

    use KindModule, only : dp

    real(dp), intent(in)  :: xi(3), zeta
    real(dp), intent(out) :: dnxi1(6), dnxi2(6), dnzeta(6)

    !! --- Logic section ---

    dnxi1(1) = (1.0_dp - zeta)*0.5_dp
    dnxi1(4) = (1.0_dp + zeta)*0.5_dp
    dnxi1(3) = -dnxi1(1)
    dnxi1(6) = -dnxi1(4)
    dnxi1(2) =  0.0_dp
    dnxi1(5) =  0.0_dp

    dnxi2(2) =  dnxi1(1)
    dnxi2(5) =  dnxi1(4)
    dnxi2(3) = -dnxi1(1)
    dnxi2(6) = -dnxi1(4)
    dnxi2(1) =  0.0_dp
    dnxi2(4) =  0.0_dp

    dnzeta(1:3) = -xi*0.5_dp
    dnzeta(4:6) =  xi*0.5_dp

  end subroutine pdvn


  subroutine ipri6stress (x, y, z, cmat, disp, stress, iel, code, iw,ipsw, ierr)

    !!==========================================================================
    !! This subroutine calculates the stresses for
    !! the 6-noded isoparametric prism element.
    !!
    !! Input:
    !!    x      - array of x-coordinates for all six nodes
    !!    y      - array of y-coordinates
    !!    z      - array of z-coordinates
    !!    cmat   - constitutive matrix
    !!    disp   - element displacement vector, vxi,vyi,vzi for each node
    !!    iel    - element number
    !!    iw     - print unit
    !!    ipsw   - print switch
    !!
    !! Output:
    !!    stress - element stress vector
    !!             sig_xx, sig_yy, sig_zz, tau_xy, tau_yz, tau_zx
    !!    ierr   - error flag, negative value signals error
    !!
    !! Coded by: Hans Magnus Haug
    !! Revised : Trond Arne Svidal,  8 Oct 1999 / 1.0
    !!==========================================================================

    use KindModule, only : dp

    integer , intent(in)  :: iel, code, iw, ipsw
    real(dp), intent(in)  :: x(6), y(6), z(6), cmat(6,6), disp(18)
    real(dp), intent(out) :: stress(6,6)
    integer , intent(out) :: ierr

    !! Local variables
    integer  :: i
    real(dp) :: strain(6,6)

    !! --- Logic section ---

    if (ipsw >= 1) then
       write(iw,*) '===> Entering ipri6stress with input ==='
       write(iw,9000) 'element number:', iel
       write(iw,9001) 'nodal x-coord: ', x
       write(iw,9001) 'nodal z-coord: ', y
       write(iw,9001) 'nodal y-coord: ', z
       write(iw,*) 'constitutive matrix:'
       write(iw,9002) (cmat(i,:),i=1,6)
       write(iw,*) 'displacement vector:'
       write(iw,9003) disp
       write(iw,*) '========================================'
    end if

    call ipri6strain (x, y, z, disp, strain, iel, code, iw, ipsw, ierr)
    if (ierr < 0) then
       write(iw,*) '*** Error return from ipri6stress'
       write(iw,*) '    ierr = ', ierr
       return
    end if

    do i = 1, 6
       stress(:,i) = matmul(cmat,strain(:,i))
    end do

    if (ipsw >= 2) then
       write(iw,*) '<=== Leaving ipri6stress with output'
       write(iw,*) 'Stress vector:'
       write(iw,9004) (i,stress(:,i),i=1,6)
       write(iw,*) '========================================'
    end if

9000 format(1X,A,I6)
9001 format(1X,A,1P,6E12.4)
9002 format(1P,6E12.4)
9003 format(1P,3E12.4)
9004 format(' Node',I2,': ',1P,6E12.4)

  end subroutine ipri6stress


  subroutine ipri6strain (x, y, z, disp, strain, iel, code, iw, ipsw, ierr)

    !!==========================================================================
    !! This subroutine calculates the strains for
    !! the 6-noded isoparametric prism element.
    !!
    !! Input:
    !!    x      - array of x-coordinates for all six nodes
    !!    y      - array of y-coordinates
    !!    z      - array of z-coordinates
    !!    disp   - element displacement vector, vxi,vyi,vzi for each node
    !!    iel    - element number
    !!    code   - extrapolation option (0: direct evaluation in nodes)
    !!    iw     - print unit
    !!    ipsw   - print switch
    !!
    !! Output:
    !!    strain - element strain vector
    !!             eps_xx, eps_yy, eps_zz, gamma_xy, gamma_yz, gamma_zx
    !!    ierr   - error flag, negative value signals error
    !!
    !! Coded by: Hans Magnus Haug
    !! Revised : Trond Arne Svidal,  8 Oct 1999 / 1.0
    !! Revised : Trond Arne Svidal, 26 Oct 1999 / 2.0
    !!
    !! Comment : 1.0: Strains are calculated in the nodes directly.
    !!                Should be calculated in the Gauss points and extrapolated.
    !!           2.0: Extrapolation implemented
    !!==========================================================================

    use KindModule, only : dp

    integer , intent(in)  :: iel, code, iw, ipsw
    real(dp), intent(in)  :: x(6), y(6), z(6), disp(18)
    real(dp), intent(out) :: strain(6,6)
    integer , intent(out) :: ierr

    !! Local variables
    integer  :: i, xpnt, zpnt
    real(dp) :: bmat(6,18), auxstrain(6,6), spntxi(3,3), spntzeta(2)

    !! --- Logic section ---

    if (ipsw >= 1) then
       write(iw,*) '===> Entering ipri6strain with input ==='
       write(iw,9000) 'element number:', iel, code
       write(iw,9001) 'nodal x-coord: ', x
       write(iw,9001) 'nodal z-coord: ', y
       write(iw,9001) 'nodal y-coord: ', z
       write(iw,*) 'displacement vector:'
       write(iw,9002) disp
       write(iw,*) '========================================'
    end if

    !! Choose sampling point coordinates

    select case (code)
    case (0) ! Calculate in the nodes directly
       spntzeta(1) = -1.0_dp
       spntzeta(2) = -spntzeta(1)
       spntxi      = 0.0_dp
       do i = 1, 3
          spntxi(i,i) = 1.0_dp
       end do

    case (1) ! Extrapolate using scheme 1
       spntzeta(1) = -1.0_dp/sqrt(3.0_dp)
       spntzeta(2) = -spntzeta(1)
       spntxi      = 1.0_dp/6.0_dp
       do i = 1, 3
          spntxi(i,i) = 2.0_dp/3.0_dp
       end do

    case (2) ! Extrapolate using scheme 2
       spntzeta(1) = -1.0_dp/sqrt(3.0_dp)
       spntzeta(2) = -spntzeta(1)
       spntxi      =  0.5_dp
       do i = 1, 3
          spntxi(i,i) = 0.0_dp
       end do

    case (3) ! Same as scheme 2, but constant in zeta direction
       spntzeta = 0.0_dp
       spntxi   = 0.5_dp
       do i = 1, 3
          spntxi(i,i) = 0.0_dp
       end do

    case default
       write(iw,*) '*** Invalid extrapolation code in ipri6strain:', code
       ierr = -1
       return
    end select

    !! Calculate strains in chosen points
    do zpnt = 1, 2
       do xpnt = 1, 3
          call ipri6bmat (x, y, z, spntxi(:,xpnt), spntzeta(zpnt), Bmat, &
               &          auxstrain(1,1), iel, iw, ipsw, ierr)
          if (ierr < 0) then
             write(iw,*) '*** Error return from ipri6strain'
             write(iw,*) '    ierr = ', ierr
             return
          end if
          strain(:,xpnt+3*zpnt-3) = matmul(Bmat,disp)
       end do
       if (code == 3) exit
    end do

    !! Extrapolate, if necessary

    select case (code)
    case (1) ! Extrapolate from Gauss points to nodes using scheme 1

       !! First extrapolate to each corner of both triangles
       auxstrain(:,1) = (5.0_dp*strain(:,1) - strain(:,2) - strain(:,3))/3.0_dp
       auxstrain(:,2) = (5.0_dp*strain(:,2) - strain(:,1) - strain(:,3))/3.0_dp
       auxstrain(:,3) = (5.0_dp*strain(:,3) - strain(:,1) - strain(:,2))/3.0_dp
       auxstrain(:,4) = (5.0_dp*strain(:,4) - strain(:,5) - strain(:,6))/3.0_dp
       auxstrain(:,5) = (5.0_dp*strain(:,5) - strain(:,4) - strain(:,6))/3.0_dp
       auxstrain(:,6) = (5.0_dp*strain(:,6) - strain(:,4) - strain(:,5))/3.0_dp

       !! Then extrapolate in heigth to nodes
       call ipri6extrapolH (auxstrain, strain)

    case (2) ! Extrapolate from Gauss points to nodes using scheme 2

       !! First extrapolate to each corner of both triangles
       auxstrain(:,1) = strain(:,2) + strain(:,3) - strain(:,1)
       auxstrain(:,2) = strain(:,1) + strain(:,3) - strain(:,2)
       auxstrain(:,3) = strain(:,1) + strain(:,2) - strain(:,3)
       auxstrain(:,4) = strain(:,5) + strain(:,6) - strain(:,4)
       auxstrain(:,5) = strain(:,4) + strain(:,6) - strain(:,5)
       auxstrain(:,6) = strain(:,4) + strain(:,5) - strain(:,6)

       !! Then extrapolate to nodes
       call ipri6extrapolH (auxstrain, strain)

    case (3) ! Extrapolate from Gauss points to nodes using scheme 3

       !! Extrapolate to each corner of the center triangle
       strain(:,4) = strain(:,2) + strain(:,3) - strain(:,1)
       strain(:,5) = strain(:,1) + strain(:,3) - strain(:,2)
       strain(:,6) = strain(:,1) + strain(:,2) - strain(:,3)

       !! Assume constant strain in zeta-direction
       strain(:,1:3) = strain(:,4:6)

    end select

    if (ipsw >= 2) then
       write(iw,*) '<=== Leaving ipri6strain with output ==='
       write(iw,*) 'Strain vector:'
       write(iw,9003) (i,strain(:,i),i=1,6)
       write(iw,*) '========================================'
    end if

9000 format(1X,A,I6,I3)
9001 format(1X,A,1P,6E12.4)
9002 format(1P,3E12.4)
9003 format(' Node',i2,': ',1P,6E12.4)

  end subroutine ipri6strain


  subroutine ipri6extrapolH (auxstrain, strain)

    !!==========================================================================
    !! This subroutine extrapolates strains from the corners of the
    !! two Gauss point planes and to the nodes.
    !! The planes are at coordinate zeta = +- 1/sqrt(3)
    !!
    !! Input:
    !!    auxstrain - element strain matrix, extrapolated in plane
    !!
    !! Output:
    !!    strain - element strain matrix, extrapolated in heigth
    !!             eps_xx, eps_yy, eps_zz, gamma_xy, gamma_yz, gamma_zx
    !!
    !! Coded by: Trond Arne Svidal
    !! Revised : Trond Arne Svidal, 26 Oct 1999 / 1.0
    !!==========================================================================

    use KindModule, only : dp

    real(dp), intent(in)  :: auxstrain(6,6)
    real(dp), intent(out) :: strain(6,6)

    !! Local variables
    real(dp) :: zm1, zp1

    !! --- Logic section ---

    zm1 = 0.5_dp*(sqrt(3.0_dp) - 1.0_dp)
    zp1 = 0.5_dp*(sqrt(3.0_dp) + 1.0_dp)

    !! Extrapolate to bottom
    strain(:,1:3) = zp1*auxstrain(:,1:3) - zm1*auxstrain(:,4:6)

    !! Extrapolate to top
    strain(:,4:6) = zp1*auxstrain(:,4:6) - zm1*auxstrain(:,1:3)

  end subroutine ipri6extrapolH

end module ipri6Module
