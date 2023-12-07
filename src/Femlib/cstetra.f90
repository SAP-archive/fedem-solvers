!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module cstetraModule

  !!============================================================================
  !! This module provides routines to generate the element stiffness- and
  !! mass-matrix as well as consistent load vectors, and to compute element
  !! stresses and strains for a constant strain tetrahedron element.
  !!============================================================================

  implicit none

  private :: pdvcoor, cstetbmat, cstetvolume


contains

  subroutine cstetbmat (x,y,z,bmat,volume,iw,ipsw,ierr)

    !!==========================================================================
    !! This subroutine generates the strain-displacement matrix for
    !! the constant strain tetrahedron element.
    !!
    !! Input:
    !!    x      - array of x-coordinates for all four nodes
    !!    y      - array of y-coordinates
    !!    z      - array of z-coordinates
    !!    iw     - print unit
    !!    ipsw   - print switch
    !!
    !! Output:
    !!    bmat   - element strain-displacement matrix
    !!             dof ordering: vxi,vyi,vzi for each node
    !!    volume - element volume
    !!    ierr   - error flag, negative value signals error
    !!
    !! Coded by Bjorn Haugen
    !! ->F90 by Knut Morten Okstad
    !!==========================================================================

    use KindModule, only : dp, epsDiv0_p

    integer , intent(in)  :: iw, ipsw
    real(dp), intent(in)  :: x(4), y(4), z(4)
    real(dp), intent(out) :: bmat(6,12), volume
    integer , intent(out) :: ierr

    !! Local variables
    integer  :: i, inod, xpos, ypos, zpos
    real(dp) :: a(4), b(4), c(4), factor

    !! --- Logic section ---

    if (ipsw >= 3) then
       write(iw,*) '=== Entering cstetbmat with input ==='
       write(iw,9001) 'nodal x-coord: ', x
       write(iw,9001) 'nodal z-coord: ', y
       write(iw,9001) 'nodal y-coord: ', z
       write(iw,*) '====================================='
    end if

    ierr = 0
    bmat = 0.0_dp

    volume = cstetvolume(x,y,z)
    if (volume > epsDiv0_p) then
       factor = 1.0_dp/(6.0_dp*volume)
    else
       ierr = -1
       write(iw,*) '*** Error return from cstetbmat'
       write(iw,*) '    Negative volume = ', volume
       return
    end if

    call pdvcoor (x,y,z,a,b,c)

    do inod = 1, 4
       xpos = (inod-1)*3 + 1
       ypos = xpos + 1
       zpos = xpos + 2

       a(inod) = a(inod)*factor
       b(inod) = b(inod)*factor
       c(inod) = c(inod)*factor

       bmat(1,xpos) = a(inod)
       bmat(2,ypos) = b(inod)
       bmat(3,zpos) = c(inod)

       bmat(4,xpos) = b(inod)
       bmat(4,ypos) = a(inod)

       bmat(5,ypos) = c(inod)
       bmat(5,zpos) = b(inod)

       bmat(6,zpos) = a(inod)
       bmat(6,xpos) = c(inod)
    end do

    if (ipsw >= 4) then
       write(iw,*) '=== Leaving cstetbmat with output ==='
       write(iw,*) 'Volume = ', volume
       write(iw,*) 'Strain-displacement matrix:'
       write(iw,9002) (bmat(i,:),i=1,6)
       write(iw,*) '====================================='
    end if

9001 format(1X,A,1P,4E12.4)
9002 format(1P,12E12.4)

  end subroutine cstetbmat


  subroutine cstetstiff (x,y,z,cmat,stiffm,volume,iw,ipsw,ierr)

    !!==========================================================================
    !! This subroutine generates the stiffness matrix for
    !! the constant strain tetrahedron element.
    !!
    !! Input:
    !!    x      - array of x-coordinates for all four nodes
    !!    y      - array of y-coordinates
    !!    z      - array of z-coordinates
    !!    cmat   - constitutive matrix
    !!    iw     - print unit
    !!    ipsw   - print switch
    !!
    !! Output:
    !!    stiffm - element stiffness matrix for current element,
    !!             dof organized as vxi,vyi,vzi for each node
    !!    volume - element volume
    !!    ierr   - error flag, negative value signals error
    !!
    !! Coded by Bjorn Haugen
    !! ->F90 by Knut Morten Okstad
    !!==========================================================================

    use KindModule, only : dp

    integer , intent(in)  :: iw, ipsw
    real(dp), intent(in)  :: x(4), y(4), z(4), cmat(6,6)
    real(dp), intent(out) :: stiffm(12,12), volume
    integer , intent(out) :: ierr

    !! Local variables
    integer  :: i, j
    real(dp) :: bmat(6,12), cb(6)

    !! --- Logic section ---

    if (ipsw >= 1) then
       write(iw,*) '=== Entering cstetstiff with input ==='
       write(iw,9001) 'nodal x-coord: ', x
       write(iw,9001) 'nodal z-coord: ', y
       write(iw,9001) 'nodal y-coord: ', z
       write(iw,*) 'constitutive matrix:'
       write(iw,9002) (cmat(i,:),i=1,6)
       write(iw,*) '======================================'
    end if

    !! Get element strain-displacement matrix

    call cstetbmat (x,y,z,bmat,volume,iw,ipsw,ierr)
    if (ierr < 0) then
       write(iw,*) '*** Error return from cstetstiff'
       write(iw,*) '    ierr = ', ierr
       return
    end if

    !! Form basic stiffness  smb = bmat'*c*bmat * volume

    do j = 1, 12
       cb = matmul(cmat,bmat(:,j))
       do i = 1, 12
          stiffm(i,j) = dot_product(bmat(:,i),cb) * volume
       end do
    end do

    if (ipsw >= 2) then
       write(iw,*) '=== Leaving cstetstiff with output ==='
       write(iw,*) 'Volume = ', volume
       write(iw,*) 'Stiffness matrix:'
       write(iw,9003) (stiffm(i,:),i=1,12)
       write(iw,*) '======================================'
    end if

9001 format(1X,A,1P,4E12.4)
9002 format(1P,6E12.4)
9003 format(1P,12E12.4)

  end subroutine cstetstiff


  subroutine cstetmass (x,y,z,rho,massm,volume,iw,ipsw,ierr)

    !!==========================================================================
    !! This subroutine generates the consistent mass matrix for
    !! the constant strain tetrahedron element.
    !!
    !! Input:
    !!    x      - array of x-coordinates for all four nodes
    !!    y      - array of y-coordinates
    !!    z      - array of z-coordinates
    !!    rho    - material density
    !!    iw     - print unit
    !!    ipsw   - print switch
    !!
    !! Output:
    !!    massm  - element mass matrix
    !!             dof organized as vxi,vyi,vzi for each node
    !!    volume - volume of the element
    !!    ierr   - error flag, negative value signals error
    !!
    !! Coded by Bjorn Haugen
    !! ->F90 by Knut Morten Okstad
    !!==========================================================================

    use KindModule, only : dp

    integer , intent(in)  :: iw, ipsw
    real(dp), intent(in)  :: x(4), y(4), z(4), rho
    real(dp), intent(out) :: massm(12,12), volume
    integer , intent(out) :: ierr

    !! Local variables
    integer  :: i, j
    real(dp) :: massdiag, massoffd

    !! --- Logic section ---

    if (ipsw >= 1) then
       write(iw,*) '=== Entering cstetmass with input ==='
       write(iw,9001) 'nodal x-coord: ', x
       write(iw,9001) 'nodal z-coord: ', y
       write(iw,9001) 'nodal y-coord: ', z
       write(iw,9001) ' mass density: ', rho
       write(iw,*) '====================================='
    end if

    ierr     = 0
    massm    = 0.0_dp
    volume   = cstetvolume(x,y,z)
    massdiag = rho*volume*0.1_dp
    massoffd = massdiag*0.5_dp

    do i = 0, 9, 3
       do j = 0, 9, 3
          massm(i+1,j+1) = massoffd
          massm(i+2,j+2) = massoffd
          massm(i+3,j+3) = massoffd
       end do
    end do

    do i = 1, 12
       massm(i,i) = massdiag
    end do

    if (ipsw >= 2) then
       write(iw,*) '=== Leaving cstetmass with output ==='
       write(iw,*) 'Volume = ', volume
       write(iw,*) 'Mass matrix:'
       write(iw,9002) (massm(i,:),i=1,12)
       write(iw,*) '====================================='
    end if

9001 format(1X,A,1P,4E12.4)
9002 format(1P,12E12.4)

  end subroutine cstetmass


  subroutine cstetpload (x,y,z,p,iface,fvec,iw,ipsw,ierr)

    !!==========================================================================
    !! This subroutine generates the consistent force vector for
    !! the constant strain tetrahedron element, due to a pressure load.
    !!
    !! Input:
    !!    x      - array of x-coordinates for all four nodes
    !!    y      - array of y-coordinates
    !!    z      - array of z-coordinates
    !!    p      - surface pressure intensities at three face nodes
    !!    iface  - index of the pressurized face, [1,4]
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
    real(dp), intent(in)  :: x(4), y(4), z(4), p(3,3)
    real(dp), intent(out) :: fvec(3,4)
    integer , intent(out) :: ierr

    !! Local variables
    integer  :: i, i1, i2, i3
    real(dp) :: s12(3), s13(3), vn(3), pw(3), dS, N1, N2, N3

    !! Four-point triangle quadrature (see Hughes 1987, Table 3.I.1, p. 173).
    real(dp), parameter :: w4_p = 1.5625_dp/3.0_dp
    real(dp), parameter :: w(4) = (/    -0.5625_dp,  w4_p,   w4_p,   w4_p  /)
    real(dp), parameter :: r(4) = (/ 1.0_dp/3.0_dp, 0.6_dp, 0.2_dp, 0.2_dp /)
    real(dp), parameter :: s(4) = (/ 1.0_dp/3.0_dp, 0.2_dp, 0.6_dp, 0.2_dp /)

    !! --- Logic section ---

    if (ipsw >= 1) then
       write(iw,*) '=== Entering cstetpload with input ==='
       write(iw,9001) 'nodal x-coord: ', x
       write(iw,9001) 'nodal z-coord: ', y
       write(iw,9001) 'nodal y-coord: ', z
       write(iw,9002) 'face index = ', iface
       write(iw,*) 'face load intensities:'
       write(iw,9003) p
       write(iw,*) '======================================'
    end if

    ierr = 0
    fvec = 0.0_dp

    select case (iface)
    case (1); i1 = 1; i2 = 3; i3 = 2
    case (2); i1 = 1; i2 = 2; i3 = 4
    case (3); i1 = 2; i2 = 3; i3 = 4
    case (4); i1 = 1; i2 = 4; i3 = 3
    case default
       ierr = -1
       write(iw,*) '*** Error return from cstetpload'
       write(iw,*) '    invalid face index = ', iface
       return
    end select

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

    !! Numerical integration
    do i = 1, 4
       N1 = r(i)
       N2 = s(i)
       N3 = 1.0_dp - N1 - N2
       pw = (p(:,1)*N1 + p(:,2)*N2 + p(:,3)*N3) * w(i)*dS
       fvec(:,i1) = fvec(:,i1) + pw*N1
       fvec(:,i2) = fvec(:,i2) + pw*N2
       fvec(:,i3) = fvec(:,i3) + pw*N3
    end do

    if (ipsw >= 2) then
       write(iw,*) '=== Leaving cstetpload with output ==='
       write(iw,*) 'Force vector:'
       write(iw,9003) fvec
       write(iw,*) '======================================'
    end if

9001 format(1X,A,1P,4E12.4)
9002 format(1X,A,I3)
9003 format(1P,3E12.4)

  end subroutine cstetpload


  subroutine cstetstrain (x,y,z,disp,strain,volume,iw,ipsw,ierr)

    !!==========================================================================
    !! This subroutine calculates the strains for
    !! the constant strain tetrahedron element.
    !!
    !! Input:
    !!    x      - array of x-coordinates for all four nodes
    !!    y      - array of y-coordinates
    !!    z      - array of z-coordinates
    !!    disp   - element displacement vector, vxi,vyi,vzi for each node
    !!    iw     - print unit
    !!    ipsw   - print switch
    !!
    !! Output:
    !!    strain - element strain vector
    !!             eps_xx, eps_yy, eps_zz, gamma_xy, gamma_yz, gamma_zx
    !!    volume - volume of the element
    !!    ierr   - error flag, negative value signals error
    !!
    !! Coded by Bjorn Haugen
    !! ->F90 by Knut Morten Okstad
    !!==========================================================================

    use KindModule, only : dp

    integer , intent(in)  :: iw, ipsw
    real(dp), intent(in)  :: x(4), y(4), z(4), disp(12)
    real(dp), intent(out) :: strain(6), volume
    integer , intent(out) :: ierr

    !!  Local variables
    real(dp) :: bmat(6,12)

    !! --- Logic section ---

    if (ipsw >= 1) then
       write(iw,*) '=== Entering cstetstrain with input ==='
       write(iw,9001) 'nodal x-coord: ', x
       write(iw,9001) 'nodal z-coord: ', y
       write(iw,9001) 'nodal y-coord: ', z
       write(iw,*) 'displacement vector:'
       write(iw,9003) disp
       write(iw,*) '======================================='
    end if

    call cstetbmat (x,y,z,bmat,volume,iw,ipsw,ierr)
    if (ierr < 0) then
       write(iw,*) '*** Error return from cstetstrain'
       write(iw,*) '    ierr = ', ierr
       return
    end if

    strain = matmul(bmat,disp)

    if (ipsw >= 2) then
       write(iw,*) '=== Leaving cstetstrain with output ==='
       write(iw,*) 'Volume = ', volume
       write(iw,*) 'Strain vector:'
       write(iw,9002) strain
       write(iw,*) '======================================='
    end if

9001 format(1X,A,1P,4E12.4)
9002 format(1P,6E12.4)
9003 format(1P,3E12.4)

  end subroutine cstetstrain


  subroutine cstetstress (x,y,z,cmat,disp,stress,volume,iw,ipsw,ierr)

    !!==========================================================================
    !! This subroutine calculates the stresses for
    !! the constant strain tetrahedron element.
    !!
    !! Input:
    !!    x      - array of x-coordinates for all four nodes
    !!    y      - array of y-coordinates
    !!    z      - array of z-coordinates
    !!    cmat   - constitutive matrix
    !!    disp   - element displacement vector, vxi,vyi,vzi for each node
    !!    iw     - print unit
    !!    ipsw   - print switch
    !!
    !! Output:
    !!    stress - element stress vector
    !!             sig_xx, sig_yy, sig_zz, tau_xy, tau_yz, tau_zx
    !!    volume - volume of the element
    !!    ierr   - error flag, negative value signals error
    !!
    !! Coded by Bjorn Haugen
    !! ->F90 by Knut Morten Okstad
    !!==========================================================================

    use KindModule, only : dp

    integer , intent(in)  :: iw, ipsw
    real(dp), intent(in)  :: x(4), y(4), z(4), cmat(6,6), disp(12)
    real(dp), intent(out) :: stress(6), volume
    integer , intent(out) :: ierr

    !! Local variables
    integer  :: i
    real(dp) :: strain(6)

    !! --- Logic section ---

    if (ipsw >= 1) then
       write(iw,*) '=== Entering cstetstress with input ==='
       write(iw,9001) 'nodal x-coord: ', x
       write(iw,9001) 'nodal z-coord: ', y
       write(iw,9001) 'nodal y-coord: ', z
       write(iw,*) 'constitutive matrix:'
       write(iw,9002) (cmat(i,:),i=1,6)
       write(iw,*) 'displacement vector:'
       write(iw,9003) disp
       write(iw,*) '======================================='
    end if

    call cstetstrain (x,y,z,disp,strain,volume,iw,ipsw,ierr)
    if (ierr < 0) then
       write(iw,*) '*** Error return from cstetstress'
       write(iw,*) '    ierr = ', ierr
       return
    end if

    stress = matmul(cmat,strain)

    if (ipsw >= 2) then
       write(iw,*) '=== Leaving cstetstress with output ==='
       write(iw,*) 'Volume = ', volume
       write(iw,*) 'Stress vector:'
       write(iw,9002) stress
    end if

9001 format(1X,A,1P,4E12.4)
9002 format(1P,6E12.4)
9003 format(1P,3E12.4)

  end subroutine cstetstress


  subroutine pdvcoor (x,y,z,a,b,c)

    !!==========================================================================
    !! Gets the partial derivatives of the volume coordinates w.r.t. x, y and z.
    !!
    !! Input:
    !!    x      - array of x-coordinate for all four nodes
    !!    y      - array of y-coordinate
    !!    z      - array of z-coordinate
    !!
    !! Output:
    !!    a      - array of partial derivatives w.r.t. x multiplied with 6*vol
    !!    b      - array of partial derivatives w.r.t. y
    !!    c      - array of partial derivatives w.r.t. z
    !!
    !! Coded by Bjorn Haugen
    !! ->F90 by Knut Morten Okstad
    !!==========================================================================

    use KindModule, only : dp

    real(dp), intent(in)  :: x(4), y(4), z(4)
    real(dp), intent(out) :: a(4), b(4), c(4)

    !! --- Logic section ---

    a(1) = (y(2)-y(3))*(z(4)-z(2)) - (y(4)-y(2))*(z(2)-z(3))
    a(2) = (y(4)-y(3))*(z(1)-z(3)) - (y(1)-y(3))*(z(4)-z(3))
    a(3) = (y(4)-y(1))*(z(2)-z(4)) - (y(2)-y(4))*(z(4)-z(1))
    a(4) = (y(2)-y(1))*(z(3)-z(1)) - (y(3)-y(1))*(z(2)-z(1))

    b(1) = (z(2)-z(3))*(x(4)-x(2)) - (z(4)-z(2))*(x(2)-x(3))
    b(2) = (z(4)-z(3))*(x(1)-x(3)) - (z(1)-z(3))*(x(4)-x(3))
    b(3) = (z(4)-z(1))*(x(2)-x(4)) - (z(2)-z(4))*(x(4)-x(1))
    b(4) = (z(2)-z(1))*(x(3)-x(1)) - (z(3)-z(1))*(x(2)-x(1))

    c(1) = (x(2)-x(3))*(y(4)-y(2)) - (x(4)-x(2))*(y(2)-y(3))
    c(2) = (x(4)-x(3))*(y(1)-y(3)) - (x(1)-x(3))*(y(4)-y(3))
    c(3) = (x(4)-x(1))*(y(2)-y(4)) - (x(2)-x(4))*(y(4)-y(1))
    c(4) = (x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1))

  end subroutine pdvcoor


  function cstetvolume (x,y,z)

    !!==========================================================================
    !! Returns the volume of the tetrahedron element.
    !!
    !! Input:
    !!    x      - array of x-coordinates for all four nodes
    !!    y      - array of y-coordinates
    !!    z      - array of z-coordinates
    !!
    !! Coded by Bjorn Haugen
    !! ->F90 by Knut Morten Okstad
    !===========================================================================

    use KindModule       , only : dp
    use ManipMatrixModule, only : cross_product

    real(dp), intent(in) :: x(4), y(4), z(4)
    real(dp)             :: cstetvolume

    !! Local variables
    real(dp) :: s12(3), s13(3), s14(3)

    !! --- Logic section ---

    !! Form side vector from node 1 to 2
    s12(1) = x(2) - x(1)
    s12(2) = y(2) - y(1)
    s12(3) = z(2) - z(1)

    !! Form side vector from node 1 to 3
    s13(1) = x(3) - x(1)
    s13(2) = y(3) - y(1)
    s13(3) = z(3) - z(1)

    !! Form side vector from node 1 to 4
    s14(1) = x(4) - x(1)
    s14(2) = y(4) - y(1)
    s14(3) = z(4) - z(1)

    !! Find volume of the tetrahedron
    cstetvolume = dot_product(cross_product(s12,s13),s14) / 6.0_dp

  end function cstetvolume

end module cstetraModule
