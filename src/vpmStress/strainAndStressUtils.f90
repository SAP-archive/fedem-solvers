!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module StrainAndStressUtilitiesModule

  implicit none

contains

  subroutine PrincipleStrains2D (epsC, eps1, eps2, gammaMax, alpha1, alphaGamma)

    !!==========================================================================
    !! Compute 2D principle strains, max shear strain, angle to principle
    !! strain and angle to max shear strain. Using Mohr's circle interpretation.
    !!
    !! NOTE: The shear strain epsC(3) is assumed to be the quantity
    !!       gamma_xy = eps_xy + eps_yx where eps_xy = eps_yx is the
    !!       tensorial shear strain component.
    !!
    !! Programmer : Bjorn Haugen
    !! date/rev   : 13 Apr 2000/1.0
    !!==========================================================================

    use KindModule, only : dp, epsDiv0_p

    real(dp),           intent(in)  :: epsC(3)
    real(dp),           intent(out) :: eps1, eps2, gammaMax
    real(dp), optional, intent(out) :: alpha1, alphaGamma

    !! Local variables
    real(dp) :: eps_12, eps_xy, origo, radius

    !! --- Logic section ---

    origo    = (epsC(1) + epsC(2))*0.5_dp
    eps_12   =  epsC(1) - epsC(2)
    eps_xy   =  epsC(3)*0.5_dp
    radius   = sqrt( eps_12*eps_12 + epsC(3)*epsC(3) )*0.5_dp

    eps1     = origo + radius
    eps2     = origo - radius
    gammaMax = radius*2.0_dp

    !! Angles to principle strain (eps1) and max shear (gammaMax)
    if (abs(eps_xy) > epsDiv0_p .or. abs(eps_12) > epsDiv0_p) then
       if (present(alpha1))     alpha1     = atan2(eps_xy,eps_12) * 0.5_dp
       if (present(alphaGamma)) alphaGamma = atan2(eps_12,eps_xy) * 0.5_dp
    else
       if (present(alpha1))     alpha1     = 0.0_dp
       if (present(alphaGamma)) alphaGamma = 0.0_dp
    end if

  end subroutine PrincipleStrains2D


  subroutine PrincipleStresses2D (sigC, sig1, sig2, tauMax, alpha1, alphaTau)

    !!==========================================================================
    !! Compute 2D principle stresses, max shear stress, angle to principle
    !! stress and angle to max shear stress. Using Mohr's circle interpretation.
    !!
    !! Programmer : Bjorn Haugen
    !! date/rev   : 13 Apr 2000/1.0
    !!==========================================================================

    use KindModule, only : dp, epsDiv0_p

    real(dp),           intent(in)  :: sigC(3)
    real(dp),           intent(out) :: sig1, sig2, tauMax
    real(dp), optional, intent(out) :: alpha1, alphaTau

    !! Local variables
    real(dp) :: sig_12, origo, radius

    !! --- Logic section ---

    origo    = (sigC(1) + sigC(2))*0.5_dp
    sig_12   =  sigC(1) - sigC(2)
    radius   = sqrt( sig_12*sig_12 + 4.0_dp*sigC(3)*sigC(3) )*0.5_dp

    sig1     = origo + radius
    sig2     = origo - radius
    tauMax   = radius

    !! Angles to principle stress (sig1) and max shear (tauMax)
    if (abs(sigC(3)) > epsDiv0_p .or. abs(sig_12) > epsDiv0_p) then
       if (present(alpha1))   alpha1   = atan2(sigC(3),sig_12) * 0.5_dp
       if (present(alphaTau)) alphaTau = atan2(sig_12,sigC(3)) * 0.5_dp
    else
       if (present(alpha1))   alpha1   = 0.0_dp
       if (present(alphaTau)) alphaTau = 0.0_dp
    end if

  end subroutine PrincipleStresses2D


  subroutine StrainDispCST (nndof, xEl, yEl, zEl, T_el, zPos, B_el)

    !!==========================================================================
    !! Compute the strain-displacement matrix for constant strain triangle.
    !!
    !! Programmer : Bjorn Haugen
    !! date/rev   : 13 Apr 2000/1.0
    !!==========================================================================

    use KindModule, only : dp, epsDiv0_p

    integer , intent(in)  :: nndof
    real(dp), intent(in)  :: xEl(:), yEl(:), zEl(:), T_el(3,3), zPos
    real(dp), intent(out) :: B_el(3,nndof,3)

    !! Local variables
    real(dp) :: xLij(3,3), yLij(3,3), vec(3), areaDobbel
    integer  :: in, jn

    !! --- Logic section ---

    B_el = 0.0_dp

    do in = 1, 3
       do jn = 1, 3
          if (in == jn) cycle
          vec(1) = xEl(in) - xEl(jn)
          vec(2) = yEl(in) - yEl(jn)
          vec(3) = zEl(in) - zEl(jn)
          xLij(in,jn) = dot_product(T_el(1,:),vec)
          yLij(in,jn) = dot_product(T_el(2,:),vec)
       end do
    end do
    areaDobbel = xLij(2,1)*yLij(3,1) - xLij(3,1)*yLij(2,1)

    !! X-displacements
    B_el(1,1,1) =  yLij(2,3)/areaDobbel
    B_el(1,1,2) =  yLij(3,1)/areaDobbel
    B_el(1,1,3) =  yLij(1,2)/areaDobbel

    B_el(3,1,1) = -xLij(2,3)/areaDobbel
    B_el(3,1,2) = -xLij(3,1)/areaDobbel
    B_el(3,1,3) = -xLij(1,2)/areaDobbel

    !! Y-displacements
    B_el(2,2,1) = -xLij(2,3)/areaDobbel
    B_el(2,2,2) = -xLij(3,1)/areaDobbel
    B_el(2,2,3) = -xLij(1,2)/areaDobbel

    B_el(3,2,1) =  yLij(2,3)/areaDobbel
    B_el(3,2,2) =  yLij(3,1)/areaDobbel
    B_el(3,2,3) =  yLij(1,2)/areaDobbel

    !! Rotational degrees of freedom
    if (nndof >= 5 .and. abs(zpos) > epsDiv0_p) then
       do in = 1, 3
          B_el(:,4,in) = -zpos*B_el(:,2,in) ! X-rotation
          B_el(:,5,in) =  zpos*B_el(:,1,in) ! Y-rotation
       end do
    end if

  end subroutine StrainDispCST


  subroutine StrainDispQuad4 (nndof, xEl, yEl, zEl, T_el, xi, eta, zPos, B_el)

    !!==========================================================================
    !! Compute the strain-displacement matrix for the 4-noded quadrilateral.
    !!
    !! Programmer : Bjorn Haugen
    !! date/rev   : 13 Apr 2000/1.0
    !!==========================================================================

    use KindModule, only : dp, epsDiv0_p

    integer , intent(in)  :: nndof
    real(dp), intent(in)  :: xEl(:), yEl(:), zEl(:), T_el(3,3), xi, eta, zPos
    real(dp), intent(out) :: B_el(3,nndof,4)

    !! Local variables
    real(dp) :: xL(4), yL(4), sfunc_x(4), sfunc_y(4), vec(3)
    integer  :: in

    !! --- Logic section ---

    B_el = 0.0_dp

    ! Local coordinates relative to node 1
    xL(1) = 0.0_dp
    yL(1) = 0.0_dp
    do in = 2, 4
       vec(1) = xEl(in) - xEl(1)
       vec(2) = yEl(in) - yEl(1)
       vec(3) = zEl(in) - zEl(1)
       xL(in) = dot_product(T_el(1,:),vec)
       yL(in) = dot_product(T_el(2,:),vec)
    end do

    call Quad4ShapeDer (xi, eta, xL, yL, sfunc_x, sfunc_y)

    !! X-displacements
    B_el(1,1,:) = sfunc_x
    B_el(3,1,:) = sfunc_y

    !! Y-displacements
    B_el(2,2,:) = sfunc_y
    B_el(3,2,:) = sfunc_x

    !! Rotational degrees of freedom
    if (nndof >= 5 .and. abs(zpos) > epsDiv0_p) then
       do in = 1, 4
          B_el(:,4,in) = -zpos*B_el(:,2,in) ! X-rotation
          B_el(:,5,in) =  zpos*B_el(:,1,in) ! Y-rotation
       end do
    end if

  end subroutine StrainDispQuad4


  subroutine Quad4Shape (xi, eta, sfunc)

    !!==========================================================================
    !! Evaluate shape functions for the 4-noded quadrilateral.
    !!
    !! Programmer : Bjorn Haugen
    !! date/rev   : 13 Apr 2000/1.0
    !!==========================================================================

    use KindModule, only : dp

    real(dp), intent(in)  :: xi, eta
    real(dp), intent(out) :: sfunc(4)

    !! --- Logic section ---

    sfunc(1) = (1.0_dp-xi)*(1.0_dp-eta)*0.25_dp
    sfunc(2) = (1.0_dp+xi)*(1.0_dp-eta)*0.25_dp
    sfunc(3) = (1.0_dp+xi)*(1.0_dp+eta)*0.25_dp
    sfunc(4) = (1.0_dp-xi)*(1.0_dp+eta)*0.25_dp

  end subroutine Quad4Shape


  subroutine Quad4ShapeDer (xi, eta, xL, yL, sfunc_x, sfunc_y)

    !!==========================================================================
    !! Evaluate shape function derivatives for the 4-noded quadrilateral.
    !!
    !! Programmer : Bjorn Haugen
    !! date/rev   : 13 Apr 2000/1.0
    !!==========================================================================

    use KindModule, only : dp

    real(dp), intent(in)  :: xi, eta, xL(4), yL(4)
    real(dp), intent(out) :: sfunc_x(4), sfunc_y(4)

    !! Local variables
    real(dp) :: sfunc_xi(4), sfunc_eta(4), jacobi(2,2), jacobi_inv(2,2), det

    !! --- Logic section ---

    !! Partial derivatives with respect to xi
    sfunc_xi(1) = -(1.0_dp-eta)*0.25_dp
    sfunc_xi(2) =  (1.0_dp-eta)*0.25_dp
    sfunc_xi(3) =  (1.0_dp+eta)*0.25_dp
    sfunc_xi(4) = -(1.0_dp+eta)*0.25_dp

    !! Partial derivatives with respect to eta
    sfunc_eta(1) = -(1.0_dp-xi)*0.25_dp
    sfunc_eta(2) = -(1.0_dp+xi)*0.25_dp
    sfunc_eta(3) =  (1.0_dp+xi)*0.25_dp
    sfunc_eta(4) =  (1.0_dp-xi)*0.25_dp

    !! Jacobi matrix  | dx/dxi   dy/dxi  |
    !!                | dx/deta  dy/deta |
    jacobi(1,1) = dot_product(sfunc_xi,xL)
    jacobi(1,2) = dot_product(sfunc_xi,yL)
    jacobi(2,1) = dot_product(sfunc_eta,xL)
    jacobi(2,2) = dot_product(sfunc_eta,yL)

    !! Invert the Jacobian, i.e.  | dxi/dx  deta/dx |
    !!                            | dxi/dy  deta/dy |
    det = jacobi(1,1)*jacobi(2,2) - jacobi(2,1)*jacobi(1,2)
    jacobi_inv(1,1) =  jacobi(2,2)/det
    jacobi_inv(2,2) =  jacobi(1,1)/det
    jacobi_inv(1,2) = -jacobi(1,2)/det
    jacobi_inv(2,1) = -jacobi(2,1)/det

    !! Partial derivatives with respect to X and Y
    sfunc_x = jacobi_inv(1,1)*sfunc_xi + jacobi_inv(1,2)*sfunc_eta
    sfunc_y = jacobi_inv(2,1)*sfunc_xi + jacobi_inv(2,2)*sfunc_eta

  end subroutine Quad4ShapeDer


  function getGlobalizedX (VZ) result(V1)

    !!==========================================================================
    !! Find the vector defined by the projection of the global X-axis onto
    !! the plane defined by the given vector VZ.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 30 Nov 2000/1.0
    !!==========================================================================

    use KindModule       , only : dp, epsDiv0_p
    use manipMatrixModule, only : cross_product

    real(dp), intent(in) :: VZ(3)
    real(dp)             :: V1(3), V2(3), VLength2
    real(dp), parameter  :: somewhatSmall_p = 0.01_dp

    !! --- Logic section ---

    if (abs(VZ(2)) > somewhatSmall_p .or. abs(VZ(3)) > somewhatSmall_p) then
       !! Define V1 by projecting global X-axis onto the plane
       V1(1) =  VZ(2)*VZ(2) + VZ(3)*VZ(3)
       V1(2) = -VZ(1)*VZ(2)
       V1(3) = -VZ(1)*VZ(3)
    else
       !! Define V2 by projecting global Y-axis onto the plane, then V1=V2xVZ
       V2(1) = -VZ(2)*VZ(1)
       V2(2) =  VZ(1)*VZ(1) + VZ(3)*VZ(3)
       V2(3) = -VZ(2)*VZ(3)
       V1    = cross_product(V2,VZ)
    end if

    VLength2 = dot_product(V1,V1)
    if (VLength2 > epsDiv0_p*epsDiv0_p) then
       V1 = V1 / sqrt(VLength2)
    else
       V1 = 0.0_dp
    end if

  end function getGlobalizedX


  subroutine getShellElementAxes (nenod, X, Y, Z, V1, V2, V3, ierr, doGlobalize)

    !!==========================================================================
    !! Compute the local- (or globalized) element axes for a thin shell element.
    !! The local-to-global transformation matrix is then T = [V1,V2,V3].
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 30 Nov 2000/1.0
    !!==========================================================================

    use KindModule       , only : dp, epsDiv0_p
    use manipMatrixModule, only : cross_product
    use reportErrorModule, only : internalError, reportError, errorFileOnly_p
    use reportErrorModule, only : getErrorFile

    integer , intent(in)  :: nenod
    real(dp), intent(in)  :: X(nenod), Y(nenod) ,Z(nenod)
    real(dp), intent(out) :: V1(3), V2(3), V3(3)
    integer , intent(out) :: ierr
    logical , optional    :: doGlobalize

    !! Local variables
    real(dp) :: VN
    logical  :: useGlobalized

    !! --- Logic section ---

    ierr = 0
    if (present(doGlobalize)) then
       useGlobalized = doGlobalize
    else
       useGlobalized = .false.
    end if

    !! Define the shell normal vector (local Z-axis)
    if (nenod == 3) then
       !! Use X12 and X13
       V1(1) = X(2) - X(1)
       V1(2) = Y(2) - Y(1)
       V1(3) = Z(2) - Z(1)
       V2(1) = X(3) - X(1)
       V2(2) = Y(3) - Y(1)
       V2(3) = Z(3) - Z(1)
    else if (nenod == 4) then
       !! Use the diagonals X13 and X24
       V1(1) = X(3) - X(1)
       V1(2) = Y(3) - Y(1)
       V1(3) = Z(3) - Z(1)
       V2(1) = X(4) - X(2)
       V2(2) = Y(4) - Y(2)
       V2(3) = Z(4) - Z(2)
    else
       ierr = internalError('getShellElementAxes: Invalid element type')
       return
    end if
    V3 = cross_product(V1,V2)
    VN = dot_product(V3,V3)
    if (VN > epsDiv0_p*epsDiv0_p) then
       V3 = V3 / sqrt(VN)
    else
       ierr = 2
       goto 900
    end if

    !! Define the local- or globalized X-axis
    if (useGlobalized) then
       V1 = getGlobalizedX(V3)
    else if (nenod == 4) then
       !! Use the vector between node 1 and 2 (equivalent to Femlib::DIRC_QUAD)
       V1(1) = X(2) - X(1)
       V1(2) = Y(2) - Y(1)
       V1(3) = Z(2) - Z(1)
       !! Project V1 onto the plane defined by the normal vector V3
       V2 = cross_product(V3,V1)
       V1 = cross_product(V2,V3)
    end if
    VN = dot_product(V1,V1)
    if (VN > epsDiv0_p*epsDiv0_p) then
       V1 = V1 / sqrt(VN)
    else
       ierr = 3
       goto 900
    end if

    !! And finally the Y-axis
    V2 = cross_product(V3,V1)

    return

900 continue
    call reportError(errorFileOnly_p,'Degenerated shell element')
    write(getErrorFile(),600) V1,V2,V3,7-2*ierr,7-2*ierr,VN,epsDiv0_p*epsDiv0_p
600 format(10X,'V1 =',1P3E13.5 / 10X,'V2 =',1P3E13.5 / 10X,'V3 =',1P3E13.5, &
         / 10X,'V',I1,'*V',I1,' = ',1PE12.5,' tol = ',E12.5)

  end subroutine getShellElementAxes


  subroutine getShellStressTrans (VX, VZ, T, ierr)

    !!==========================================================================
    !! Compute the 2D transformation matrix from the element coordinate system
    !! to the continuous stress output coordinate system for shell elements.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 30 Nov 2000/1.0
    !!==========================================================================

    use KindModule       , only : dp, epsDiv0_p
    use manipMatrixModule, only : cross_product
    use reportErrorModule, only : internalError

    real(dp), intent(in)  :: VX(3)  ! X-axis of the element coordinate system
    real(dp), intent(in)  :: VZ(3)  ! Z-axis of the element coordinate system
    real(dp), intent(out) :: T(2,2) ! in-plane transformation matrix
    integer , intent(out) :: ierr

    !! Local variables
    real(dp) :: V1(3), V2(3), CA, SA

    !! --- Logic section ---

    ierr = 0

    !! Find the local X-axis (V1) of the stress output coordinate system
    V1 = getGlobalizedX(VZ)
    if (dot_product(V1,V1) <= epsDiv0_p*epsDiv0_p) then
       ierr = internalError('getShellStressTrans: Invalid element axes')
       return
    end if

    !! Find the rotation from VX to V1 in terms of cos(angle) and sin(angle)
    V2 = cross_product(VX,V1)
    CA = dot_product(V1,VX)
    SA = sign(sqrt(dot_product(V2,V2)),dot_product(V2,VZ))

    !! Set up the in-plane rotation matrix
    T(1,1) =  CA
    T(1,2) =  SA
    T(2,1) = -SA
    T(2,2) =  CA

  end subroutine getShellStressTrans


  subroutine calcVonMises (inTensor,ncomp,nstrp,vmVec)

    !!==========================================================================
    !! Calculate von Mises stress from input stress tensor.
    !!
    !! Programmer : Tommy Stokmo Jørstad
    !! date/rev   : 14 May 2002/1.0
    !!==========================================================================

    use KindModule                  , only : dp
    use FFaTensorTransformsInterface, only : vonMises

    integer,  intent(in)  :: ncomp, nstrp
    real(dp), intent(in)  :: inTensor(ncomp,nstrp)
    real(dp), intent(out) :: vmVec(nstrp)

    !! Local variables
    integer :: i, ncalc

    !! --- Logic section ---

    select case(ncomp)
    case(1); ncalc = 1 ! Beam element
    case(3); ncalc = 2 ! Shell element
    case(6); ncalc = 3 ! Solid element
    case default; ncalc = 0
    end select

    do i = 1, nstrp
       vmVec(i) = vonMises(ncalc,inTensor(:,i))
    end do

  end subroutine calcVonMises


  subroutine calcPrincipalVals (inTensor,ncomp,nstrp,maxPVal,minPVal,maxSVal)

    !!==========================================================================
    !! Calculate max/min principal value and max shear from input tensor.
    !!
    !! Programmer : Tommy Stokmo Jørstad
    !! date/rev   : 14 May 2002/1.0
    !!==========================================================================

    use KindModule                  , only : dp
    use FFaTensorTransformsInterface, only : princval, maxshearvalue

    integer,  intent(in)  :: ncomp, nstrp
    real(dp), intent(in)  :: inTensor(ncomp,nstrp)
    real(dp), intent(out) :: maxPVal(nstrp), minPVal(nstrp), maxSVal(nstrp)

    !! Local variables
    integer  :: i, ncalc
    real(dp) :: prinValVec(3)

    !! --- Logic section ---

    select case(ncomp)
    case(1); ncalc = 1 ! Beam element
    case(3); ncalc = 2 ! Shell element
    case(6); ncalc = 3 ! Solid element
    case default; ncalc = 0
    end select

    do i = 1, nstrp
       call princval (ncalc,inTensor(:,i),prinValVec)
       maxPVal(i) = prinValVec(1)
       minPVal(i) = prinValVec(ncalc)
       call maxshearvalue (ncalc,prinValVec,maxSVal(i))
    end do

  end subroutine calcPrincipalVals

end module StrainAndStressUtilitiesModule
