!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module Andes4ShellModule

  !!============================================================================
  !! This module contains subroutines for calculating the element stiffness
  !! matrix of a 4-noded shell element based on the ANDES formulation.
  !!
  !! The code is developed by Bjorn Haugen (NTNU), April 2012.
  !! Adapted to FEDEM standard by Knut Morten Okstad, April 2016.
  !!============================================================================

  implicit none

  private

  public :: Andes4shell_stiffmat


contains

  subroutine StrainDispMatrix4ben (lumpType,alphaH,x,y,z,xsi,eta,Bmat,lpu,ierr)

    !***************************************************************************
    ! Computes the strain-displacement matrix for a 4-noded bending element.
    ! The strains are computed using the projection rule.
    !
    ! Input arguments:
    ! lumpType : Order of normal rotation along edges
    ! alphaH   : Scaling factor for the higher-order part
    ! x, y, z  : Delta coordinates
    ! xsi, eta : Natural coordinates of current integration point
    ! lpu      : Logical print unit
    !
    ! Output arguments:
    ! Bmat     : Strain-displacement matrix with cartesian curvatures
    !            Nodal dof ordering: 4 x-rot., 4 y-rot., 4 z-disp.
    ! ierr     : Error flag
    !***************************************************************************

    use KindModule, only : dp

    integer , intent(in)  :: lumpType
    real(dp), intent(in)  :: alphaH, x(4,4), y(4,4), z(4,4), eta, xsi
    real(dp), intent(out) :: Bmat(3,12)
    integer , intent(in)  :: lpu
    integer , intent(out) :: ierr

    real(dp) :: a, Lt(3,12), Bh(3,12)

    ! Area of the quadrilateral
    a = 0.5_dp*(x(2,1)*y(3,1) - x(3,1)*y(2,1) + x(4,3)*y(1,3) - x(1,3)*y(4,3))

    ! Compute the transpose of the lumping matrix
    call lump4ben (lumpType,x,y,Lt)

    ! Compute the deviatoric part of the strain-displacement matrix
    call bh4ben (x,y,z,xsi,eta,Bh,lpu,ierr)
    if (ierr < 0) return

    ! Total strain-displacement matrix as the sum of the lumping matrix divided
    ! by the area, and the scaled higher-order strain-displacement matrix
    Bmat = Lt/a + Bh*alphaH

  end subroutine StrainDispMatrix4ben


  subroutine StrainDispMatrix4mem (alphaNR,alphaH,x,y,z,xsi,eta,Bmat,lpu,ierr)

    !***************************************************************************
    ! Computes the strain-displacement matrix for a 4-noded membrane element.
    !
    ! Input arguments:
    ! alphaNR  : Scaling factor for the normal rotation part
    ! alphaH   : Scaling factor for the higher-order part
    ! x, y, z  : Delta coordinates
    ! xsi, eta : Natural coordinates of current integration point
    ! lpu      : Logical print unit
    !
    ! Output arguments:
    ! Bmat     : Strain-displacement matrix with natural membrane strains
    !            Nodal dof ordering: 4 x-disp., 4 y-disp., 4 normal rotations
    ! ierr     : Error flag
    !***************************************************************************

    use KindModule, only : dp

    real(dp), intent(in)  :: alphaNR, alphaH, x(4,4), y(4,4), z(4,4), xsi, eta
    real(dp), intent(out) :: Bmat(3,12)
    integer , intent(in)  :: lpu
    integer , intent(out) :: ierr

    real(dp) :: a, Lt(3,12), Bh(3,12)

    ! Area of the quadrilateral
    a = 0.5_dp*(x(2,1)*y(3,1) - x(3,1)*y(2,1) + x(4,3)*y(1,3) - x(1,3)*y(4,3))

    ! Compute the transpose of the lumping matrix
    call lump4mem (x,y,alphaNR,Lt)

    ! Compute the deviatoric part of the strain-displacement matrix
    call bh4mem (x,y,z,xsi,eta,Bh,lpu,ierr)
    if (ierr < 0) return

    ! Total strain-displacement matrix as the sum of the lumping matrix divided
    ! by the area, and the scaled higher-order strain-displacement matrix
    Bmat = Lt/a + Bh*alphaH

  end subroutine StrainDispMatrix4mem


  subroutine lump4ben (ltype,x,y,lt)

    !***************************************************************************
    ! Computes the transposed lumping matrix for a 4-noded bending element.
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
    !         Nodal dof ordering: 4 x-rot., 4 y-rot., 4 normal displacements
    !***************************************************************************

    use KindModule, only : dp

    integer , intent(in)  :: ltype
    real(dp), intent(in)  :: x(4,4), y(4,4)
    real(dp), intent(out) :: lt(3,4,3)

    integer  :: i, j, l
    real(dp) :: len2, ccij, ssij, csij, ccli, ssli, csli

    do i = 1, 4
       j = 1 + mod(i,4)
       l = 1 + mod(i+2,4)
       select case (ltype)

       case (1) ! Linear interpolation

          lt(1,i,1) =  0.0_dp
          lt(2,i,1) = -x(j,l)*0.5_dp
          lt(3,i,1) =  y(j,l)*0.5_dp
          lt(1,i,2) = -y(j,l)*0.5_dp
          lt(2,i,2) =  0.0_dp
          lt(3,i,2) =  x(j,l)*0.5_dp
          lt(1,i,3) =  0.0_dp
          lt(2,i,3) =  0.0_dp
          lt(3,i,3) =  0.0_dp

       case (2) ! Quadratic interpolation

          len2 = x(i,j)*x(i,j) + y(i,j)*y(i,j)
          ccij = x(i,j)*x(i,j)/len2
          ssij = y(i,j)*y(i,j)/len2
          csij = x(i,j)*y(i,j)/len2

          len2 = x(l,i)*x(l,i) + y(l,i)*y(l,i)
          ccli = x(l,i)*x(l,i)/len2
          ssli = y(l,i)*y(l,i)/len2
          csli = x(l,i)*y(l,i)/len2

          lt(1,i,1) = (ssij*x(i,j) + ssli*x(l,i))*0.5_dp
          lt(2,i,1) = (ccij*x(i,j) + ccli*x(l,i))*0.5_dp
          lt(3,i,1) = -ccij*y(i,j) - ccli*y(l,i)
          lt(1,i,2) = (ssij*y(i,j) + ssli*y(l,i))*0.5_dp
          lt(2,i,2) = (ccij*y(i,j) + ccli*y(l,i))*0.5_dp
          lt(3,i,2) = -ssij*x(i,j) - ssli*x(l,i)
          lt(1,i,3) = -csij+csli
          lt(2,i,3) =  csij-csli
          lt(3,i,3) = (ssli-ccli) - (ssij-ccij)

       end select

    end do

  end subroutine lump4ben


  subroutine bben_nod_quad (xij,yij,bbenn,lpu,ierr)

    !***************************************************************************
    ! Computes the nodal strain-displacement matrices for the bending strains
    ! of a quadrilateral element with respect to the twelve visible dofs.
    ! *** Element is assumed "best fit" in the xy-plane ****
    !
    ! Input:
    !   xij : delta coordinate, xij(i,j) = x(i) - x(j)
    !   yij : delta coordinate, yij(i,j) = y(i) - y(j)
    !   lpu : Logical print unit
    !
    ! Output:
    !   bbenn : nodal strain-displacement matrices
    !           The strain ordering is: kappa_xx, kappa_yy, 2*kappa_xy
    !           The nodal dof ordering is:
    !              w1, w2, w3, w4, tx1, tx2, tx3, tx4, ty1, ty2, ty3, ty4
    !***************************************************************************

    use KindModule       , only : dp
    use ManipMatrixModule, only : invert33

    real(dp), intent(in)  :: xij(4,4), yij(4,4)
    real(dp), intent(out) :: bbenn(4,3,12)
    integer , intent(in)  :: lpu
    integer , intent(out) :: ierr

    integer  :: stnod(3,4) = reshape((/4,2,3,1,3,4,2,4,1,3,1,2/),(/3,4/))
    integer  :: i, j, k, inod, jnod
    real(dp) :: lij(4,4), sij(4,4,2), nij(4,4,2), tran(4,3,3), tinv(3,3), bxy(3)

    ierr = 0

    ! Initialize nodal strain-displacement matrix
    bbenn = 0.0_dp

    ! Compute side lenghts and tangential and normal side unit vectors

    do i = 1, 4
       do k = 1, 4
          j = 1 + mod(i+k-1,4)

          lij(i,j) = sqrt(xij(i,j)*xij(i,j)+yij(i,j)*yij(i,j))
          lij(j,i) = lij(i,j)

          sij(i,j,1) =  xij(j,i)/lij(i,j)
          sij(i,j,2) =  yij(j,i)/lij(i,j)

          nij(i,j,1) =  sij(i,j,2)
          nij(i,j,2) = -sij(i,j,1)

          sij(j,i,1) = -sij(i,j,1)
          sij(j,i,2) = -sij(i,j,2)

          nij(j,i,1) = -nij(i,j,1)
          nij(j,i,2) = -nij(i,j,2)

       end do
    end do

    ! Compute nodal strain-displacement matrices
    do inod = 1, 4
       do i = 1, 3
          jnod = stnod(i,inod)

          bbenn(inod,i,inod)   = -6.0_dp / (lij(inod,jnod)*lij(inod,jnod))
          bbenn(inod,i,inod+4) = -4.0_dp * nij(inod,jnod,1)/lij(inod,jnod)
          bbenn(inod,i,inod+8) = -4.0_dp * nij(inod,jnod,2)/lij(inod,jnod)

          bbenn(inod,i,jnod)   = -bbenn(inod,i,inod)
          bbenn(inod,i,jnod+4) =  bbenn(inod,i,inod+4)*0.5_dp
          bbenn(inod,i,jnod+8) =  bbenn(inod,i,inod+8)*0.5_dp
       end do
    end do

    ! Compute nodal strain transformation matrices
    do inod = 1, 4
       do i = 1, 3
          jnod = stnod(i,inod)
          tinv(i,1) = sij(inod,jnod,1)*sij(inod,jnod,1)
          tinv(i,2) = sij(inod,jnod,2)*sij(inod,jnod,2)
          tinv(i,3) = sij(inod,jnod,1)*sij(inod,jnod,2)
       end do

       tran(inod,:,:) = invert33(tinv,lpu,ierr)
       if (ierr < 0) return

    end do

    ! Transform nodal strain-displacement matrices to Cartesian strains
    do inod = 1, 4
       do j = 1, 12
          do i = 1, 3
             bxy(i) = tran(inod,i,1)*bbenn(inod,1,j) &
                  & + tran(inod,i,2)*bbenn(inod,2,j) &
                  & + tran(inod,i,3)*bbenn(inod,3,j)
          end do
          bbenn(inod,:,j) = bxy
       end do
    end do

  end subroutine bben_nod_quad


  subroutine bh4ben (x,y,z,xsi,eta,Bben,lpu,ierr)

    !***************************************************************************
    ! Computes strain-displacement matrix for a 4-noded bending element.
    !
    ! Input arguments:
    ! x, y, z  : Delta coordinates
    ! xsi, eta : Natural coordinates of current integration point
    ! lpu      : Logical print unit
    !
    ! Output arguments:
    ! Bben : Strain-displacement matrix with cartesian curvatures
    !        Nodal dof ordering: 4 x-rot., 4 y-rot., 4 z-disp.
    ! ierr : Error flag
    !***************************************************************************

    use KindModule, only : dp

    real(dp), intent(in)  :: x(4,4), y(4,4), z(4,4), xsi, eta
    real(dp), intent(out) :: Bben(3,12)
    integer , intent(in)  :: lpu
    integer , intent(out) :: ierr

    integer, parameter :: nint = 2

    integer  :: i, j, k
    real(dp) :: Bbenn(4,3,12), Bben_mean(3,12), nv(4), xl(4,3)

    ! bbenn has dof ordering 4 z-displ, 4 x-rot, 4 y-rot
    call bben_nod_quad (x,y,Bbenn,lpu,ierr)
    if (ierr < 0) return

    ! Calculate mean
    xl(:,1) = x(:,1)
    xl(:,2) = y(:,1)
    xl(:,3) = z(:,1)
    call bmean_quad (Bbenn,xl,nint,Bben_mean)

    call shape4 (xsi,eta,nv)
    do i = 1, 3
       do j = 1, 12
          k = 1 + mod(7+j,12)
          Bben(i,k) = dot_product(nv,Bbenn(:,i,j)) - Bben_mean(i,j)
       end do
    end do

  end subroutine bh4ben


  subroutine bmean_quad (Bnods,xl,nint,Bmean)

    !***************************************************************************
    ! Computes mean strain-displacement matrix for a quadrilateral element with
    ! respect to the visible DOFs, by means of numerical integration.
    !
    ! Input arguments:
    ! Bnods   : Nodal Cartesian strain-displacement matrices
    !           Strain ordering: kappa_xx, kappa_yy, 2*kappa_xy
    !           Nodal dof ordering: w1,...,w4, tx1,...,tx4, ty1,...,ty4
    ! xl(i,j) : Coordinate j of node i.
    !           The element is assumed to be in the x-y plane, i.e., xl(:,3)=0.
    ! nint    : Number of integration points in xi- and eta-direction
    !
    ! Output arguments:
    ! Bmean   : Mean of the nodal strain-displacement matrix
    !           Strain ordering: kappa_xx, kappa_yy, 2*kappa_xy
    !           Nodal dof ordering: w1,...,w4, tx1,...,tx4, ty1,...,ty4
    !***************************************************************************

    use KindModule, only : dp

    real(dp), intent(in)  :: bnods(:,:,:), xl(4,3)
    integer , intent(in)  :: nint
    real(dp), intent(out) :: bmean(:,:)

    integer  :: i, j
    real(dp) :: intcoord(3), weight(3), n(4), area, fac

    call getGaussForm (nint,intcoord,weight)

    Bmean = 0.0_dp
    area  = 0.0_dp

    do i = 1, nint
       do j = 1, nint
          fac = jacobi_quad(intcoord(i),intcoord(j),xl)*weight(i)*weight(j)
          call shape4 (intcoord(i),intcoord(j),n)
          n = n * fac
          bmean = bmean + n(1)*bnods(1,:,:) + n(2)*bnods(2,:,:) &
               &        + n(3)*bnods(3,:,:) + n(4)*bnods(4,:,:)
          area  = area + fac
       end do
    end do

    bmean = bmean/area

  end subroutine bmean_quad


  function jacobi_quad (xsi,eta,xyz,jacobi) result(d_area)

    !***************************************************************************
    ! Computes the surface jacobian at the point given by (xi,eta).
    !
    ! Input arguments:
    ! xsi, eta : Natural coordinates
    ! xyz      : Array of nodal coordinates (xyz(i,j) is coordinate j at node i)
    !
    ! Output arguments:
    ! jacobi   : Array containing the jacobi matrix (optional)
    !            jacobi(1,:) is gradient vector in xi-direction (g_xi)
    !            jacobi(2,:) is gradient vector in eta-direction (g_eta)
    !
    ! Return value:
    ! d_area   : Differential area at current point
    !***************************************************************************

    use KindModule       , only : dp
    use ManipMatrixModule, only : cross_product

    real(dp), intent(in)            :: xsi, eta, xyz(4,3)
    real(dp), intent(out), optional :: jacobi(2,3)

    integer  :: i
    real(dp) :: d_area, nv_xsi(4), nv_eta(4), g_xsi(3), g_eta(3), n(3)

    ! Shape function values and gradients at the point (xi,eta)
    call shape4 (xsi,eta,nv_xsi=nv_xsi,nv_eta=nv_eta)

    ! Forming the jacobi matrix by gradient vectors
    do i = 1, 3
       g_xsi(i) = dot_product(nv_xsi,xyz(:,i))
       g_eta(i) = dot_product(nv_eta,xyz(:,i))
    end do
    if (present(jacobi)) then
       jacobi(1,:) = g_xsi
       jacobi(2,:) = g_eta
    end if

    ! Normal vector as cross product of g_xi and g_eta
    n = cross_product(g_xsi,g_eta)

    ! Differential area as lenght of normal vector
    d_area = sqrt(n(1)*n(1) + n(2)*n(2) + n(3)*n(3))

  end function jacobi_quad


  subroutine lump4mem (x,y,alpha,lt)

    !***************************************************************************
    ! Computes the transposed lumping matrix for a 4-noded membrane element.
    !
    ! Input arguments:
    ! x, y  : Delta coordinates
    ! alpha : Scaling factor for rotational contribution
    !
    ! Output arguments:
    ! lt    : Transposed lumping matrix
    !         Nodal dof ordering: 4 x-disp., 4 y-disp., 4 normal rotations
    !***************************************************************************

    use KindModule, only : dp

    real(dp), intent(in)  :: x(4,4), y(4,4), alpha
    real(dp), intent(out) :: lt(3,4,3)

    integer :: i, j, l

    do i = 1, 4
       j = 1 + mod(i,4)
       l = 1 + mod(i+2,4)
       lt(1,i,1) =  y(j,l)*0.5_dp
       lt(2,i,1) =  0.0_dp
       lt(3,i,1) = -x(j,l)*0.5_dp
       lt(1,i,2) =  0.0_dp
       lt(2,i,2) = -x(j,l)*0.5_dp
       lt(3,i,2) =  y(j,l)*0.5_dp
       lt(1,i,3) = (y(l,i)*y(l,i) - y(j,i)*y(j,i))*alpha/12.0_dp
       lt(2,i,3) = (x(l,i)*x(l,i) - x(j,i)*x(j,i))*alpha/12.0_dp
       lt(3,i,3) = (x(j,i)*y(j,i) - x(l,i)*y(l,i))*alpha/6.0_dp
    end do

  end subroutine lump4mem


  subroutine bmqh_nod_quad (xij,yij,bmqhn,lpu,ierr)

    !***************************************************************************
    ! Computes the nodal higher-order strain displacement matrices for the
    ! membrane strains of a quadrilateral element with respect to the
    ! six higher order dofs.
    ! *** Element is assumed "best fit" in the xy-plane ****
    !
    ! Input:
    !   xij : delta coordinate, xij(i,j) = x(i) - x(j)
    !   yij : delta coordinate, yij(i,j) = y(i) - y(j)
    !   lpu : Logical print unit
    !
    ! Output:
    !   bmqhn : higher-order strain displacement matrix
    !           The strain ordering is: e_xx, e_yy, gamma_xy
    !           The nodal dof ordering is : rd_1,rd_2,rd_3,rd_4,rtor,txi,teta
    !           rd_i = r_i - r0 where r0 is the rigid body constant strain
    !           rotation of the element (normal/drilling rotations).
    !***************************************************************************

    use KindModule       , only : dp
    use ManipMatrixModule, only : invert33

    real(dp), intent(in)  :: xij(4,4), yij(4,4)
    real(dp), intent(out) :: bmqhn(4,3,7)
    integer , intent(in)  :: lpu
    integer , intent(out) :: ierr

    real(dp), parameter :: rho1 =  0.1_dp
    real(dp), parameter :: rho2 = -rho1
    real(dp), parameter :: rho3 =  rho2
    real(dp), parameter :: rho4 =  rho1

    real(dp), parameter :: alpha = 0.1_dp
    real(dp), parameter :: beta1 = 0.5_dp + rho1
    real(dp), parameter :: beta2 = 0.0_dp

    real(dp), parameter :: rho5 =  0.0_dp
    real(dp), parameter :: rho6 =  beta1 - rho1
    real(dp), parameter :: rho7 =  0.0_dp
    real(dp), parameter :: rho8 = -rho6

    integer  :: i, j, inod
    real(dp) :: rc1, rc2, rvec(2,4)
    real(dp) :: Xi(2), Eta(2), dv13(2), dv24(2), unit13(2), unit24(2)
    real(dp) :: bxy(3), tran13(3,3), tinv13(3,3), tran24(3,3), tinv24(3,3)
    real(dp) :: len, lXi, lEta, l13, l24
    real(dp) :: facXi(4), facEta(4), facXiMean, facEtaMean
    real(dp) :: fac13, fac24, dp13xi, dp13eta, dp24xi, dp24eta

    !  Compute vectors from centroid to nodes,
    !  xi- and eta unit vectors, diagonal unit vectors etc.

    rc1 = (xij(2,1) + xij(3,1) + xij(4,1)) * 0.25_dp
    rc2 = (yij(2,1) + yij(3,1) + yij(4,1)) * 0.25_dp

    do inod = 1, 4
       rvec(1,inod) = xij(inod,1) - rc1
       rvec(2,inod) = yij(inod,1) - rc2
    end do

    Xi(1)   = (rvec(1,2) + rvec(1,3)) - (rvec(1,1) + rvec(1,4))
    Xi(2)   = (rvec(2,2) + rvec(2,3)) - (rvec(2,1) + rvec(2,4))
    len     = sqrt(Xi(1)*Xi(1) + Xi(2)*Xi(2))
    Xi      = Xi/len
    lXi     = len*0.5_dp

    Eta(1)  = (rvec(1,3) + rvec(1,4)) - (rvec(1,1) + rvec(1,2))
    Eta(2)  = (rvec(2,3) + rvec(2,4)) - (rvec(2,1) + rvec(2,2))
    len     = sqrt(Eta(1)*Eta(1) + Eta(2)*Eta(2))
    Eta     = Eta/len
    lEta    = len*0.5_dp

    dv13(1) = xij(3,1)
    dv13(2) = yij(3,1)
    l13     = sqrt(dv13(1)*dv13(1) + dv13(2)*dv13(2))
    unit13  = dv13/l13

    dv24(1) = xij(4,2)
    dv24(2) = yij(4,2)
    l24     = sqrt(dv24(1)*dv24(1) + dv24(2)*dv24(2))
    unit24  = dv24/l24

    ! Compute strain distribution factors
    do i = 1, 4
       facXi(i)  = cross2d(Xi, rvec(:,i))/lXi
       facEta(i) = cross2d(Eta,rvec(:,i))/lEta
    end do

    facXiMean  = (facXi(1)  + facXi(2)  + facXi(3)  + facXi(4) ) * 0.25_dp
    facEtaMean = (facEta(1) + facEta(2) + facEta(3) + facEta(4)) * 0.25_dp

    do i = 1, 4
       facXi(i)  = facXi(i)  - facXiMean
       facEta(i) = facEta(i) - facEtaMean
    end do

    ! Change signs so that all factors are positive
    facXi(1:2)  = -facXi(1:2)
    facEta(2:3) = -facEta(2:3)
    facXiMean   = ( facXi(1) + facXi(2)  + facXi(3)  + facXi(4) ) * 0.25_dp
    facEtaMean  = (facEta(1) + facEta(2) + facEta(3) + facEta(4)) * 0.25_dp

    ! Diagonal strains
    fac13 =  cross2d(unit13,dv24) * 0.5_dp/l13
    fac24 = -cross2d(unit24,dv13) * 0.5_dp/l24

    dp13xi  = unit13(1)* Xi(1) + unit13(2)* Xi(2)
    dp13eta = unit13(1)*Eta(1) + unit13(2)*Eta(2)
    dp24xi  = unit24(1)* Xi(1) + unit24(2)* Xi(2)
    dp24eta = unit24(1)*Eta(1) + unit24(2)*Eta(2)

    ! Compute natural coordinat strain displacement matrix at
    ! all four nodes (strains with respect to higher order dofs)

    ! Node 1
    bmqhn(1,1,1) =  rho1*facXi(1)
    bmqhn(1,1,2) =  rho2*facXi(1)
    bmqhn(1,1,3) =  rho3*facXi(1)
    bmqhn(1,1,4) =  rho4*facXi(1)
    bmqhn(1,1,5) =  alpha*lEta/lXi
    bmqhn(1,1,6) = -beta1*facXi(1)/(facXiMean*lXi)
    bmqhn(1,1,7) =  0.0_dp

    bmqhn(1,2,1) = -rho1*facEta(1)
    bmqhn(1,2,2) = -rho4*facEta(1)
    bmqhn(1,2,3) = -rho3*facEta(1)
    bmqhn(1,2,4) = -rho2*facEta(1)
    bmqhn(1,2,5) = -alpha*lXi/lEta
    bmqhn(1,2,6) =  0.0_dp
    bmqhn(1,2,7) = -beta1*facEta(1)/(facEtaMean*lEta)

    bmqhn(1,3,1) =  rho5*fac24
    bmqhn(1,3,2) =  rho6*fac24
    bmqhn(1,3,3) =  rho7*fac24
    bmqhn(1,3,4) =  rho8*fac24
    bmqhn(1,3,5) =  0.0_dp
    bmqhn(1,3,6) =  beta2*dp24xi/l24
    bmqhn(1,3,7) = -beta2*dp24eta/l24

    ! Node 2
    bmqhn(2,1,1) = -rho2*facXi(2)
    bmqhn(2,1,2) = -rho1*facXi(2)
    bmqhn(2,1,3) = -rho4*facXi(2)
    bmqhn(2,1,4) = -rho3*facXi(2)
    bmqhn(2,1,5) = -alpha*lEta/lXi
    bmqhn(2,1,6) = -beta1*facXi(2)/(facXiMean*lXi)
    bmqhn(2,1,7) =  0.0_dp

    bmqhn(2,2,1) =  rho4*facEta(2)
    bmqhn(2,2,2) =  rho1*facEta(2)
    bmqhn(2,2,3) =  rho2*facEta(2)
    bmqhn(2,2,4) =  rho3*facEta(2)
    bmqhn(2,2,5) =  alpha*lXi/lEta
    bmqhn(2,2,6) =  0.0_dp
    bmqhn(2,2,7) =  beta1*facEta(2)/(facEtaMean*lEta)

    bmqhn(2,3,1) =  rho8*fac13
    bmqhn(2,3,2) =  rho5*fac13
    bmqhn(2,3,3) =  rho6*fac13
    bmqhn(2,3,4) =  rho7*fac13
    bmqhn(2,3,5) =  0.0_dp
    bmqhn(2,3,6) = -beta2*dp13xi/l13
    bmqhn(2,3,7) =  beta2*dp13eta/l13

    ! Node 3
    bmqhn(3,1,1) =  rho3*facXi(3)
    bmqhn(3,1,2) =  rho4*facXi(3)
    bmqhn(3,1,3) =  rho1*facXi(3)
    bmqhn(3,1,4) =  rho2*facXi(3)
    bmqhn(3,1,5) =  alpha*lEta/lXi
    bmqhn(3,1,6) =  beta1*facXi(3)/(facXiMean*lXi)
    bmqhn(3,1,7) =  0.0_dp

    bmqhn(3,2,1) = -rho3*facEta(3)
    bmqhn(3,2,2) = -rho2*facEta(3)
    bmqhn(3,2,3) = -rho1*facEta(3)
    bmqhn(3,2,4) = -rho4*facEta(3)
    bmqhn(3,2,5) = -alpha*lXi/lEta
    bmqhn(3,2,6) =  0.0_dp
    bmqhn(3,2,7) =  beta1*facEta(3)/(facEtaMean*lEta)

    bmqhn(3,3,1) =  rho7*fac24
    bmqhn(3,3,2) =  rho8*fac24
    bmqhn(3,3,3) =  rho5*fac24
    bmqhn(3,3,4) =  rho6*fac24
    bmqhn(3,3,5) =  0.0_dp
    bmqhn(3,3,6) = -beta2*dp24xi/l24
    bmqhn(3,3,7) =  beta2*dp24eta/l24

    ! Node 4
    bmqhn(4,1,1) = -rho4*facXi(4)
    bmqhn(4,1,2) = -rho3*facXi(4)
    bmqhn(4,1,3) = -rho2*facXi(4)
    bmqhn(4,1,4) = -rho1*facXi(4)
    bmqhn(4,1,5) = -alpha*lEta/lXi
    bmqhn(4,1,6) =  beta1*facXi(4)/(facXiMean*lXi)
    bmqhn(4,1,7) =  0.0_dp

    bmqhn(4,2,1) =  rho2*facEta(4)
    bmqhn(4,2,2) =  rho3*facEta(4)
    bmqhn(4,2,3) =  rho4*facEta(4)
    bmqhn(4,2,4) =  rho1*facEta(4)
    bmqhn(4,2,5) =  alpha*lXi/lEta
    bmqhn(4,2,6) =  0.0_dp
    bmqhn(4,2,7) = -beta1*facEta(4)/(facEtaMean*lEta)

    bmqhn(4,3,1) =  rho6*fac13
    bmqhn(4,3,2) =  rho7*fac13
    bmqhn(4,3,3) =  rho8*fac13
    bmqhn(4,3,4) =  rho5*fac13
    bmqhn(4,3,5) =  0.0_dp
    bmqhn(4,3,6) =  beta2*dp13xi/l13
    bmqhn(4,3,7) = -beta2*dp13eta/l13

    ! Compute strain transformation matrix for node 1-3 and 2-4, and
    ! transform the nodal strain displacement matrices to x-y strains

    tinv24(1,1) = Xi(1)*Xi(1)
    tinv13(1,1) = tinv24(1,1)
    tinv24(1,2) = Xi(2)*Xi(2)
    tinv13(1,2) = tinv24(1,2)
    tinv24(1,3) = Xi(1)*Xi(2)
    tinv13(1,3) = tinv24(1,3)
    tinv24(2,1) = Eta(1)*Eta(1)
    tinv13(2,1) = tinv24(2,1)
    tinv24(2,2) = Eta(2)*Eta(2)
    tinv13(2,2) = tinv24(2,2)
    tinv24(2,3) = Eta(1)*Eta(2)
    tinv13(2,3) = tinv24(2,3)
    tinv13(3,1) = unit24(1)*unit24(1)
    tinv13(3,2) = unit24(2)*unit24(2)
    tinv13(3,3) = unit24(1)*unit24(2)
    tinv24(3,1) = unit13(1)*unit13(1)
    tinv24(3,2) = unit13(2)*unit13(2)
    tinv24(3,3) = unit13(1)*unit13(2)

    tran13 = invert33(tinv13,lpu,ierr)
    tran24 = invert33(tinv24,lpu,ierr)

    ! Tranform to Cartesian strains
    do inod = 1, 3, 2
       do j = 1, 7
          do i = 1, 3
             bxy(i) = tran13(i,1)*bmqhn(inod,1,j) &
                  & + tran13(i,2)*bmqhn(inod,2,j) &
                  & + tran13(i,3)*bmqhn(inod,3,j)
          end do
          bmqhn(inod,:,j) = bxy
       end do
    end do
    do inod = 2, 4, 2
       do j = 1, 7
          do i = 1, 3
             bxy(i) = tran24(i,1)*bmqhn(inod,1,j) &
                  & + tran24(i,2)*bmqhn(inod,2,j) &
                  & + tran24(i,3)*bmqhn(inod,3,j)
          end do
          bmqhn(inod,:,j) = bxy
       end do
    end do

  contains

    function cross2d(a,b)
      real(dp), intent(in) :: a(2), b(2)
      real(dp) :: cross2d
      cross2d = a(1)*b(2) - b(1)*a(2)
    end function cross2d

  end subroutine bmqh_nod_quad


  subroutine hmqv_quad (xij,yij,zij,hmem,lpu,ierr)

    !***************************************************************************
    ! Computes transformation from higher order membrane dofs to visible dofs.
    !
    ! Input:
    !   xij(i,j) : delta coordinate x(i) -x(j)
    !   yij(i,j) : delta coordinate y(i) -y(j)
    !   zij(i,j) : delta coordinate z(i) -z(j) (dummy so far )
    !   lpu      : Logical print unit
    !
    ! Output:
    !   hmem : transformation from dof ordering:
    !            rd_1, rd_2, rd_3, rd_4, rtor, txi, teta
    !            (rd_i = r_i-rtor-r0 where r0 is the rigid body constant strain)
    !          The visible dof ordering is
    !          vx1, vx2, vx3, vx4, vy1, vy2, vy3, vy4, rz1, rz2, rz3, rz4
    !***************************************************************************

    use KindModule       , only : dp
    use ManipMatrixModule, only : invert22

    real(dp), intent(in)  :: xij(4,4), yij(4,4), zij(4,4)
    real(dp), intent(out) :: hmem(7,12)
    integer,  intent(in)  :: lpu
    integer,  intent(out) :: ierr

    integer  :: i, inod
    real(dp) :: fac, rc1, rc2, rc3, Xi1, Xi2, Eta1, Eta2
    real(dp) :: rvec(4,3), rnat(4,3), tnat(2,2), jacobi0(2,3)

    ! Compute vectors from centroid to nodes,
    ! xi- and eta unit vectors, diagonal unit vectors etc.

    rc1 = (xij(2,1) + xij(3,1) + xij(4,1)) * 0.25_dp
    rc2 = (yij(2,1) + yij(3,1) + yij(4,1)) * 0.25_dp
    rc3 = (zij(2,1) + zij(3,1) + zij(4,1)) * 0.25_dp

    do inod = 1, 4
       rvec(inod,1) = xij(inod,1) - rc1
       rvec(inod,2) = yij(inod,1) - rc2
       rvec(inod,3) = zij(inod,1) - rc3
    end do

    Xi1 = (rvec(2,1) + rvec(3,1)) - (rvec(1,1) + rvec(4,1))
    Xi2 = (rvec(2,2) + rvec(3,2)) - (rvec(1,2) + rvec(4,2))
    fac = sqrt(Xi1*Xi1 + Xi2*Xi2)
    Xi1 = Xi1 / fac
    Xi2 = Xi2 / fac

    Eta1 = (rvec(3,1) + rvec(4,1)) - (rvec(1,1) + rvec(2,1))
    Eta2 = (rvec(3,2) + rvec(4,2)) - (rvec(1,2) + rvec(2,2))
    fac  = sqrt(Eta1*Eta1 + Eta2*Eta2)
    Eta1 = Eta1 / fac
    Eta2 = Eta2 / fac

    ! Compute transformation from visible dofs to hierarcial rotations
    ! and higher order translational dofs

    fac = 0.0625_dp / jacobi_quad(0.0_dp,0.0_dp,rvec,jacobi0)

    ! Hierarcical rotations
    do i = 1, 4
       hmem(i,1:8)  =  0.0_dp
       hmem(i,9:12) = -0.25_dp
       hmem(i,i+8)  =  0.75_dp
    end do
    hmem(5,1) = -xij(2,4)*fac
    hmem(5,2) = -xij(3,1)*fac
    hmem(5,3) = -xij(4,2)*fac
    hmem(5,4) = -xij(1,3)*fac
    hmem(5,5) = -yij(2,4)*fac
    hmem(5,6) = -yij(3,1)*fac
    hmem(5,7) = -yij(4,2)*fac
    hmem(5,8) = -yij(1,3)*fac
    hmem(5,9:12) = 0.25_dp

    ! Higher order translations
    tnat = invert22(transpose(jacobi0(:,1:2)),lpu,ierr)
    if (ierr < 0) return

    do inod = 1, 4
       do i = 1, 2
          rnat(inod,i) = tnat(i,1)*rvec(inod,1) + tnat(i,2)*rvec(inod,2)
       end do
    end do

    fac = (rnat(1,1)*rnat(1,2) + rnat(2,1)*rnat(2,2) &
         + rnat(3,1)*rnat(3,2) + rnat(4,1)*rnat(4,2)) * 0.25_dp

    hmem(6,1) = (rnat(3,1)*rnat(3,2) - fac) * Xi1
    hmem(6,2) = (rnat(4,1)*rnat(4,2) - fac) * Xi1
    hmem(6,3) = (rnat(1,1)*rnat(1,2) - fac) * Xi1
    hmem(6,4) = (rnat(2,1)*rnat(2,2) - fac) * Xi1
    hmem(6,5) = (rnat(3,1)*rnat(3,2) - fac) * Xi2
    hmem(6,6) = (rnat(4,1)*rnat(4,2) - fac) * Xi2
    hmem(6,7) = (rnat(1,1)*rnat(1,2) - fac) * Xi2
    hmem(6,8) = (rnat(2,1)*rnat(2,2) - fac) * Xi2
    hmem(6,9:12) = 0.0_dp

    hmem(7,1) = (rnat(3,1)*rnat(3,2) - fac) * Eta1
    hmem(7,2) = (rnat(4,1)*rnat(4,2) - fac) * Eta1
    hmem(7,3) = (rnat(1,1)*rnat(1,2) - fac) * Eta1
    hmem(7,4) = (rnat(2,1)*rnat(2,2) - fac) * Eta1
    hmem(7,5) = (rnat(3,1)*rnat(3,2) - fac) * Eta2
    hmem(7,6) = (rnat(4,1)*rnat(4,2) - fac) * Eta2
    hmem(7,7) = (rnat(1,1)*rnat(1,2) - fac) * Eta2
    hmem(7,8) = (rnat(2,1)*rnat(2,2) - fac) * Eta2
    hmem(7,9:12) = 0.0_dp

  end subroutine hmqv_quad


  subroutine bh4mem (x,y,z,xsi,eta,Bmem,lpu,ierr)

    !***************************************************************************
    ! Computes the strain-displacement matrix for a 4-noded membrane element.
    !
    ! Input arguments:
    ! x, y, z  : Delta coordinates
    ! xsi, eta : Natural coordinates of current integration point
    ! lpu      : Logical print unit
    !
    ! Output arguments:
    ! Bmem : Strain-displacement matrix with natural membrane strains
    !        Nodal dof ordering: 4 x-disp., 4 y-disp., 4 normal rotations
    ! ierr : Error flag (ierr < 0: Edge -ierr has zero or negative length)
    !***************************************************************************

    use KindModule, only : dp

    real(dp), intent(in)  :: x(4,4), y(4,4), z(4,4), xsi, eta
    real(dp), intent(out) :: Bmem(3,12)
    integer , intent(in)  :: lpu
    integer , intent(out) :: ierr

    integer, parameter :: nint = 2

    integer  :: i, j
    real(dp) :: Bmemn(4,3,7), Bmemq(3,7), Bmem_mean(3,7), Hmem(7,12), nv(4)
    real(dp) :: xl(4,3)

    call bmqh_nod_quad (x,y,Bmemn,lpu,ierr)

    ! Calculate mean
    xl(:,1) = x(:,1)
    xl(:,2) = y(:,1)
    xl(:,3) = z(:,1)
    call bmean_quad (bmemn,xl,nint,bmem_mean)

    call shape4 (xsi,eta,nv)
    do i = 1, 3
       do j = 1, 7
          bmemq(i,j) = dot_product(nv,Bmemn(:,i,j)) - bmem_mean(i,j)
       end do
    end do

    call hmqv_quad (x,y,z,Hmem,lpu,ierr)

    bmem = matmul(bmemq,hmem)

  end subroutine bh4mem


  subroutine Bmatrix (lumpType,alphaH,alphaNR,xij,yij,zij,xsi,eta,Bmat,lpu,ierr)

    !***************************************************************************
    ! Creates the final B-matrix, based on the membrane and bending parts.
    !
    ! Input:
    ! lumpType : Order of normal rotation along edges
    ! alphaH   : Scaling factor for the higher-order part
    ! alphaNR  : Scaling factor for the normal rotation part
    ! xij-zij  : Delta coordinates
    ! xsi,eta  : Natural coordinates of current integration point
    ! lpu      : Logical print unit
    !
    ! Output:
    ! Bmat     : Final B-matrix with dof ordering
    !            u1, v1, w1, x-rot1, y-rot1, z-rot1, u2, v2, w2 ... z-rot4
    ! ierr     : Error flag
    !***************************************************************************

    use KindModule, only : dp

    integer , intent(in)  :: lumpType, lpu
    real(dp), intent(in)  :: xij(4,4), yij(4,4), zij(4,4), alphaH, alphaNR
    real(dp), intent(in)  :: xsi, eta
    real(dp), intent(out) :: Bmat(6,24)
    integer , intent(out) :: ierr

    integer  :: i, j
    real(dp) :: Bb(3,12), Bm(3,12)

    integer, parameter :: dofm(12) =(/ 1, 7,13,19,  2, 8,14,20,  6,12,18,24 /)
    integer, parameter :: dofb(12) =(/ 4,10,16,22,  5,11,17,23,  3, 9,15,21 /)

    ! Create Bm
    call StrainDispMatrix4mem (alphaNR,alphaH,xij,yij,zij,xsi,eta,Bm,lpu,ierr)
    if (ierr < 0) return

    ! Create Bb
    call StrainDispMatrix4ben (lumptype,alphaH,xij,yij,zij,xsi,eta,Bb,lpu,ierr)
    if (ierr < 0) return

    ! Create the final B-matrix(6,24)
    Bmat = 0.0_dp
    do j = 1, 12
       do i = 1, 3
          Bmat(i,  dofm(j)) = Bm(i,j)
          Bmat(i+3,dofb(j)) = Bb(i,j)
       end do
    end do

  end subroutine Bmatrix


  subroutine getGaussForm (nmult,intPos,weight)

    !***************************************************************************
    ! Returns the correct Gauss quadrature point coordinates and weights.
    !
    ! Input:
    ! nmult  : Number of quadrature points
    !
    ! Output:
    ! intPos : Quadrature point coordinates
    ! weight : The weight of each quadrature point
    !***************************************************************************

    use kindModule, only : dp

    integer , intent(in)  :: nmult
    real(dp), intent(out) :: intPos(nmult), weight(nmult)

    select case (nMult)
    case (1)
       intPos = 0.0_dp
       weight = 2.0_dp
    case (2)
       intPos(1) = -sqrt(1.0_dp/3.0_dp)
       intPos(2) = -intPos(1)
       weight    =  1.0_dp
    case (3)
       intPos(1) = -sqrt(3.0_dp/5.0_dp)
       intPos(2) =  0.0_dp
       intPos(3) = -intPos(1)
       weight(1) =  5.0_dp/9.0_dp
       weight(2) =  8.0_dp/9.0_dp
       weight(3) =  weight(1)
    case default
       intPos = 0.0_dp
       weight = 0.0_dp
    end select

  end subroutine getGaussForm


  subroutine shape4 (xsi,eta,nv,nv_xsi,nv_eta)

    !***************************************************************************
    ! Evaluates the shape functions and/or their gradient.
    !
    ! Input arguments:
    ! xsi, eta : Natural coordinates of the point to evaluate at
    !
    ! Output arguments:
    ! nv     : Shape function values
    ! nv_xsi : Shape function derivatives w.r.t. xsi
    ! nv_eta : Shape function derivatives w.r.t. eta
    !***************************************************************************

    use kindModule, only : dp

    real(dp), intent(in)            :: xsi, eta
    real(dp), intent(out), optional :: nv(4), nv_xsi(4), nv_eta(4)

    if (present(nv)) then
       nv(1) = (1.0_dp-xsi)*(1.0_dp-eta)*0.25_dp
       nv(2) = (1.0_dp+xsi)*(1.0_dp-eta)*0.25_dp
       nv(3) = (1.0_dp+xsi)*(1.0_dp+eta)*0.25_dp
       nv(4) = (1.0_dp-xsi)*(1.0_dp+eta)*0.25_dp
    end if

    if (present(nv_xsi)) then
       nv_xsi(1) = -(1.0_dp-eta)*0.25_dp
       nv_xsi(2) =  (1.0_dp-eta)*0.25_dp
       nv_xsi(3) =  (1.0_dp+eta)*0.25_dp
       nv_xsi(4) = -(1.0_dp+eta)*0.25_dp
    end if

    if (present(nv_eta)) then
       nv_eta(1) = -(1.0_dp-xsi)*0.25_dp
       nv_eta(2) = -(1.0_dp+xsi)*0.25_dp
       nv_eta(3) =  (1.0_dp+xsi)*0.25_dp
       nv_eta(4) =  (1.0_dp-xsi)*0.25_dp
    end if

  end subroutine shape4


  subroutine Andes4shell_stiffmat (Xg,Yg,Zg,Cmat,alpha,beta,lumpType,Kmat, &
       &                           lpu,IERR)

    !***************************************************************************
    ! Creates the stiffness matrix for a 4-noded shell element.
    !***************************************************************************

    use KindModule, only : dp

    real(dp), intent(in)  :: Xg(4), Yg(4), Zg(4), Cmat(6,6), alpha, beta
    real(dp), intent(out) :: Kmat(24,24)
    integer , intent(in)  :: lumpType, lpu
    integer , intent(out) :: ierr

    integer, parameter :: nmult = 2
    integer            :: i, j
    real(dp)           :: xi(nmult), weight(nmult), wJ, Bmat(6,24)
    real(dp)           :: alpH, dx(4,4), dy(4,4), dz(4,4), xl(4,3)

    ! Gauss quadrature coordinates and weights
    call getGaussForm (nmult,xi,weight)

    ! Delta values for the nodal coordinates
    do i = 1, 4
       xl(i,1) = xg(i)
       xl(i,2) = yg(i)
       xl(i,3) = zg(i)
       do j = 1, 4
          dx(i,j) = xg(i) - xg(j)
          dy(i,j) = yg(i) - yg(j)
          dz(i,j) = zg(i) - zg(j)
       end do
    end do

    ! Summation of the stiffness matrix (B_i^t*C*B_i)
    ierr = 0
    Kmat = 0.0_dp
    alpH = sqrt(beta)
    do i = 1, nMult
        do j = 1, nMult
          wJ = jacobi_quad(xi(i),xi(j),xl)*weight(i)*weight(j)
          call Bmatrix (lumpType,alpH,alpha,dx,dy,dz,xi(i),xi(j),Bmat,lpu,ierr)
          if (ierr < 0) exit
          Kmat = Kmat + matmul(transpose(Bmat),matmul(Cmat*wJ,Bmat))
        end do
    end do

    if (ierr < 0) write(lpu,"(' *** Andes4shell_stiffmat failed')")

  end subroutine Andes4shell_stiffmat

end module Andes4ShellModule
