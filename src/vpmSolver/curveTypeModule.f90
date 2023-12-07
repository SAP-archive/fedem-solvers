!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file curveTypeModule.f90
!>
!> @brief Subroutines for contact surface point calculations.

!!==============================================================================
!> @brief Module with subroutines for contact surface point calculations.

module CurveTypeModule

  implicit none

  private :: Cubic3D, adjustTheta


contains

  !!============================================================================
  !> @brief Computes position variables along a contact curve.
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 4 Nov 2016

  subroutine GetPosOnCurve (jointId, interpolationOrder, curvePts, &
       &                    slideVar, loopLength, isExtended, isLooping, &
       &                    PosOnCurveG, tanVec, coeff, rotGrad, ierr)

    use CurvePointTypeModule, only : CurvePointType
#ifdef FT_DEBUG
    use dbgUnitsModule   , only : dbgCurve
#endif
    use KindModule       , only : dp, epsDiv0_p
    use RotationModule   , only : deltaRot, vec_to_mat, quat_to_mat, mat_to_quat
    use ManipMatrixModule, only : matmul34, cross_product
    use ReportErrorModule, only : error_p, reportError

    character(len=*)    , intent(in)  :: jointId
    integer             , intent(in)  :: interpolationOrder
    type(CurvePointType), intent(in)  :: curvePts(:)
    real(dp)            , intent(in)  :: slideVar, loopLength
    logical             , intent(in)  :: isExtended, isLooping
    real(dp)            , intent(out) :: PosOnCurveG(3,4), tanVec(3), coeff(:)
    real(dp), optional  , intent(out) :: rotGrad(3)
    integer             , intent(out) :: ierr

    !! Local variables
    integer  :: i, i0, i1, i2, i3, nPt
    real(dp) :: Pos1G(3), Pos2G(3), R1(3,3), R2(3,3), Rmat(3,3), q(4), qPos(4)
    real(dp) :: Xn(3,4), vec1(3), vec2(3), orig(3), secVec(3), uVec(3), nVec(3)
    real(dp) :: len, lenSec, halfSoR, phiHalf, theta, Xsi, Eta, deltaSlideVar

    !! --- Logic section ---

    ierr  = 0
    coeff = 0.0_dp

    Xsi = slideVar
    if (isLooping) then
       !! Add or subtract multiple revolutions for looping joints
       do while (Xsi < 0.0_dp)
          Xsi = Xsi + loopLength
       end do
       do while (Xsi > loopLength)
          Xsi = Xsi - loopLength
       end do
    end if

    nPt = size(curvePts)
    if (Xsi > curvePts(nPt)%slideVar) then
       if (isLooping) then
          !! Slide value within last interval of looping joint
          i1 = nPt
          i2 = 1
       else if (isExtended) then
          !! Special handling of extended joint, past end of track
          i1 = nPt - 1
          i2 = nPt
       else
          ierr = -1
          call reportError (error_p,'Slider is past end of track for '//jointId)
          return
       end if
    else if (Xsi < curvePts(1)%slideVar) then
       if (isLooping) then
          !! Slide value within last interval of looping joint
          i1 = nPt
          i2 = 1
       else if (isExtended) then
          !! Special handling of extended joint, before beginning of track
          i1 = 1
          i2 = 2
       else
          ierr = -2
          call reportError (error_p, &
               &            'Slider is before beginning of track for '//jointId)
          return
       end if
    else
       !! Slide value within ordinary joint interval
       i2 = 2
       do while (Xsi > curvePts(i2)%slideVar)
          i2 = i2 + 1
       end do
       i1 = i2 - 1
    end if
    i0 = i1
    i3 = i2

    if (isLooping .and. i1 > i2) then
       deltaSlideVar = loopLength - curvePts(i1)%slideVar
    else
       deltaSlideVar = curvePts(i2)%slideVar - curvePts(i1)%slideVar
    end if

    Xsi = (Xsi - curvePts(i1)%slideVar)/deltaSlideVar
    Eta = 1.0_dp - Xsi

    Pos1G = matmul34(curvePts(i1)%triad%ur,curvePts(i1)%CPosInT(:,4))
    Pos2G = matmul34(curvePts(i2)%triad%ur,curvePts(i2)%CPosInT(:,4))

    secVec = 0.5_dp*(Pos2G - Pos1G)           ! Half of secant vector
    lenSec = sqrt(dot_product(secVec,secVec)) ! Half secant length (S/2)

    !! Curvature is stored at the first of the two curve points of the section
    halfSoR = lenSec*curvePts(i1)%curvature ! S/(2*R)
    !! Arc height divided by secant length: h/S = S/(8*R) = 0.25*halfSoR
    !! This is deduced by assuming that h is small compared to R, together with
    !! the trigonometric relations sin(phi/2) = S/(2*R) and cos(phi/2) = (R-h)/R
    if (abs(halfSoR) < 0.0004_dp) then ! Assume straight line when h/S < 0.0001

       PosOnCurveG(:,4) = Pos1G*Eta + Pos2G*Xsi

       i0 = max(i1-1,1)
       i3 = min(i2+1,nPt)
       if (interpolationOrder >= 3 .and. i3-i0 > 1) then ! Cubic interpolation

          nPt = 0
          do i = i0, i3
             nPt = nPt + 1
             Xn(:,nPt) = matmul34(curvePts(i)%triad%ur,curvePts(i)%CPosInT(:,4))
          end do
          if (i0 == i1) then
             len = 0.5_dp
          else
             len = 1.5_dp
          end if
          call Cubic3D (Xn(:,1:nPt),len,PosOnCurveG(:,4),tanVec,coeff(i0:i3))

          if (interpolationOrder == 4 .and. .not.present(rotGrad)) then
             !! Try cubic interpolation of the orientation using quaternions
             qPos = 0.0_dp
             do i = i0, i3
                Rmat = matmul(curvePts(i)%triad%ur(:,1:3), &
                     &        curvePts(i)%CPosInT(:,1:3))
                call mat_to_quat (Rmat,q)
                qPos = qPos + coeff(i)*q
             end do
             call quat_to_mat (qPos,PosOnCurveG(:,1:3))
             return
          end if

       else ! Linear interpolation

          if (lenSec > epsDiv0_p) then
             tanVec = secVec/lenSec
          else
             tanVec = 0.0_dp
          end if

          coeff(i1) = Eta
          coeff(i2) = Xsi
          nPt = 2

       end if

    else ! Circular curve segment

       uVec = matmul(curvePts(i1)%triad%ur(:,1:3),curvePts(i1)%upVecInT1)*Eta &
            + matmul(curvePts(i2)%triad%ur(:,1:3),curvePts(i2)%upVecInT2)*Xsi
       nVec = cross_product(secVec,uVec)
       len = sqrt(dot_product(nVec,nVec))
       if (len > epsDiv0_p) nVec = nVec/len      ! Unit normal to curve plane

       if (halfSoR >= 1.0_dp) then
          phiHalf = asin(1.0_dp)                 ! This is a 180 degree curve
          vec1    = 0.0_dp
       else if (halfSoR <= -1.0_dp) then
          phiHalf = asin(-1.0_dp)                ! This is a -180 degree curve
          vec1    = 0.0_dp
       else
          phiHalf = asin(halfSoR)                ! Half the angle of the curve
          vec1    = cross_product(secVec,nVec)/tan(phiHalf)
       end if

       orig = Pos1G + vec1 + secVec              ! Coordinates of the origin
       vec1 = Pos2G - orig                       ! From origin to P2

       theta = 2.0_dp*phiHalf*Eta                ! Angle from P2 to pos on curve

       vec2 = nVec*theta                         ! Rotation axis vector
       call vec_to_mat(vec2,Rmat)                ! Rotation matrix
       vec2 = matmul(Rmat,vec1)                  ! From origin to pos on curve

       PosOnCurveG(:,4) = orig + vec2

       tanVec = cross_product(vec2,nVec)
       len = sqrt(dot_product(tanVec,tanVec))
       if (len > epsDiv0_p) tanVec = tanVec/len

       coeff(i1) = Eta
       coeff(i2) = Xsi
       nPt = 2

    end if

    !! Linear interpolation of the orientation
    R1 = matmul(curvePts(i1)%triad%ur(:,1:3),curvePts(i1)%CPosInT(:,1:3))
    R2 = matmul(curvePts(i2)%triad%ur(:,1:3),curvePts(i2)%CPosInT(:,1:3))
    uVec = deltaRot(R1,R2)
    call vec_to_mat(Xsi*uVec,Rmat)
    PosOnCurveG(:,1:3) = matmul(Rmat,R1)

    if (present(rotGrad)) rotGrad = uVec/deltaSlideVar

#ifdef FT_DEBUG
    if (dbgCurve > 0) then
       if (interpolationOrder < 3 .or. nPt == 2) then
          Xn(:,1) = matmul34(curvePts(i1)%triad%ur,curvePts(i1)%CPosInT(:,4))
          Xn(:,2) = matmul34(curvePts(i2)%triad%ur,curvePts(i2)%CPosInT(:,4))
          len = Xsi
       end if
       write(dbgCurve,"(/6X,'Position on ',A)") trim(jointId)
       write(dbgCurve,"(6X,'Xn =',1P3E13.5/(10X,1P3E13.5))") Xn(:,1:nPt)
       write(dbgCurve,"(6X,'SlideVars =',4F8.3)") (curvePts(i)%slideVar,i=i0,i3)
       write(dbgCurve,"(6X,'SlidePos, Xsi =',2F8.3)") slideVar,len
       write(dbgCurve,"(6X,'PosOnCurve =',1P3E13.5)") PosOnCurveG(:,4)
       write(dbgCurve,"(6X,'X-axis =',1P3E13.5)") PosOnCurveG(:,1)
       write(dbgCurve,"(6X,'Y-axis =',1P3E13.5)") PosOnCurveG(:,2)
       write(dbgCurve,"(6X,'Z-axis =',1P3E13.5)") PosOnCurveG(:,3)
       write(dbgCurve,"(6X,'coeffs =',10F7.3/(14X,10F7.3))") coeff
    end if
#endif

  end subroutine GetPosOnCurve


  !!============================================================================
  !> @brief Checks if a given point @a Pos is on a circular curve segment.
  !>
  !> @details The curve segment is defined from two points @a CP1 and @a CP2.
  !> If the point is found to be on the segment, compute the local position and
  !> orientation (@a Xsi and @a Def) of the point w.r.t. to the curve segment.
  !> The curve tangent vectors at the curve points are also updated based on
  !> current triad configuration (for use in the subsequent corner test).
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 2001

  subroutine GetPosOnCircSection (Pos, DefPrev, CP1, CP2, &
       &                          thickness, width, radius, &
       &                          isInDomain, PosOnCurveG, Xsi, Def)

    use CurvePointTypeModule, only : CurvePointType
#ifdef FT_DEBUG
    use dbgUnitsModule   , only : dbgCurve
#endif
    use KindModule       , only : dp, pi_p
    use RotationModule   , only : FFa_eulerZYX
    use ManipMatrixModule, only : cross_product, matmul34, invert34

    real(dp)            , intent(in)    :: Pos(:,:), DefPrev(3)
    type(CurvePointType), intent(inout) :: CP1, CP2
    real(dp)            , intent(in)    :: thickness, width, radius
    logical             , intent(out)   :: isInDomain
    real(dp)            , intent(out)   :: PosOnCurveG(:,:), Xsi, Def(:)

    !! Local variables
    real(dp) :: Pos1G(3,4), Pos2G(3,4), upVec(3), curvature, vecToPos(3), &
         &      tanVec(3), len, R, &
         &      phiHalf, phi, nVec(3), secVec(3), vec1(3), vec2(3), origo(3), &
         &      OP(3), OPproj(3), sinTheta, theta, lenToPos, hOverS
    real(dp) :: unitM(3), unitSec(3), lenSec
    real(dp) :: phiHalfNew, RNew, l1, l2, d1, d2
    real(dp), parameter :: EPS = 1.0e-30_dp

    !! --- Logic section ---

    isInDomain = .false.

    Pos1G = matmul34(CP1%triad%ur,CP1%CPosInT)
    Pos2G = matmul34(CP2%triad%ur,CP2%CPosInT)
#ifdef FT_DEBUG
    if (dbgCurve > 0) then
       write(dbgCurve,"(6X,'Point 1:',1P3E13.5)") Pos1G(:,4)
       write(dbgCurve,"(6X,'Point 2:',1P3E13.5)") Pos2G(:,4)
    end if
#endif

    secVec = Pos2G(:,4) - Pos1G(:,4)
    lenSec = sqrt(dot_product(secVec,secVec))

    !! Curvature is stored at the first point of the two point curve section
    curvature = CP1%curvature
    hOverS = 0.125_dp*lenSec*curvature ! Arc height divided by secant length

    if (abs(hOverS) < 0.0001_dp) then  ! Assume straight line

#ifdef FT_DEBUG
       if (dbgCurve > 0) then
          write(dbgCurve,"('  === Testing domain USING STRAIGHT LINE')")
       end if
#endif

       tanVec = secVec/lenSec
       CP1%tan1G = tanVec
       CP2%tan2G = tanVec

       vecToPos = Pos(:,4) - Pos1G(:,4)
       lenToPos = dot_product(vecToPos,tanVec)
       Xsi      = lenToPos/lenSec

       if (Xsi >= 0.0_dp .and. Xsi <= 1.0_dp) then
          isInDomain = .true.

          !! Establish curve coordinate system and position
          vec2 = Pos1G(:,2)*(1.0_dp-Xsi) + Pos2G(:,2)*Xsi
          vec1 = cross_product(vec2,tanVec)
          PosOnCurveG(:,1) = vec1/sqrt(dot_product(vec1,vec1))
          PosOnCurveG(:,3) = tanVec
          PosOnCurveG(:,2) = cross_product(PosOnCurveG(:,3),PosOnCurveG(:,1))
          PosOnCurveG(:,4) = Pos1G(:,4) + lenToPos*tanVec
       end if

    else

#ifdef FT_DEBUG
       if (dbgCurve > 0) then
          write(dbgCurve,"('  === Testing domain USING CIRCULAR SECTION')")
       end if
#endif

       !TODO,bh: should use (1-xsi) and xsi as factors here
       upVec = matmul(CP1%triad%ur(:,1:3),CP1%upVecInT1)*0.5_dp + &
            &  matmul(CP2%triad%ur(:,1:3),CP2%upVecInT2)*0.5_dp

       R = 1.0_dp/curvature ! Curve radius

       !! Half the angle of the curve (TODO,bh: acos is better for large angles)
       if (lenSec < R+R) then
          phiHalf = asin(0.5_dp*lenSec/R)
       else
          phiHalf = asin(1.0_dp)
       end if

       if (abs(CP1%h0) > EPS) then
          phiHalfNew = atan(2.0_dp*CP1%h0/lenSec) * 2.0_dp
          RNew       = 0.5_dp*lenSec/sin(phiHalfNew)
#ifdef FT_DEBUG
          if (dbgCurve > 0) then
             write(dbgCurve,"('h0',1PE12.3)") CP1%h0
             write(dbgCurve,"('PhiHalf, orig and adjusted',1P2E12.3,1PF6.1,'%')") &
                  &    phiHalf, phiHalfNew, (PhiHalfNew-PhiHalf)/PhiHalf
             write(dbgCurve,"('Radius, fixed and adjusted',1P2E12.3,1PF6.1,'%')") &
                  &    R, Rnew, (RNew-R)/R
          end if
#endif
          phiHalf = phiHalfNew
          R       = RNew
       end if

       phi = phiHalf*2.0_dp ! Angle of the curve

       nVec = cross_product(secVec,upVec) ! Unit normal to the curve plane
       nVec = nVec/sqrt(dot_product(nVec,nVec))

       unitM = cross_product(nVec,secVec) ! Unit vector from origo to mid secant
       unitM = unitM/sqrt(dot_product(unitM,unitM))

       origo = (Pos1G(:,4) + Pos2G(:,4))*0.5_dp - unitM*cos(phiHalf)*R ! Coordinate of origo
       vec1 = Pos1G(:,4) - origo          ! Vector from origo to P1
       vec2 = Pos2G(:,4) - origo          ! Vector from origo to P2
       l1 = sqrt(dot_product(vec1,vec1))
       l2 = sqrt(dot_product(vec2,vec2))
       CP1%tan1G = cross_product(vec1,nVec)/l1
       CP2%tan2G = cross_product(vec2,nVec)/l2

#ifdef FT_DEBUG
       if (dbgCurve > 0) then
          write(dbgCurve,"('R l1 l2 :',1P3E12.3)") R, l1, l2
          write(dbgCurve,"('tan1G   :',1P3E12.3)") CP1%tan1G
          write(dbgCurve,"('tan2G   :',1P3E12.3)") CP2%tan2G
          write(dbgCurve,"('Origo   :',1P3E12.3)") origo
       end if
#endif

       OP = Pos(:,4) - origo
       OPproj = OP - nVec*dot_product(OP,nVec)

       len = sqrt(dot_product(OPproj,OPproj))
       d1 =  dot_product(OPproj-vec1,CP1%tan1G)
       d2 = -dot_product(OPproj-vec2,CP2%tan2G)

       if (len > EPS .and. d1 > 0.0_dp .and. d2 > 0.0_dp) then
          unitSec = secVec/lenSec
          sinTheta = dot_product(OPproj,unitSec)/len
          if (abs(sinTheta) >= 1.0_dp) then
             theta = sign(asin(1.0_dp),sinTheta)
          else
             theta = asin(sinTheta)
          end if

#ifdef bjorn
          if (theta > 0.0_dp) then
             phiHalf = abs(asin(dot_product(vec2,unitSec)*curvature))
          else
             phiHalf = abs(asin(dot_product(vec1,unitSec)*curvature))
          end if
#endif

#if defined(FT_DEBUG) || defined(bjorn)
          OP = Pos1G(:,4) - Pos(:,4)
          l1 = sqrt(dot_product(OP,OP))
          OP = Pos2G(:,4) - Pos(:,4)
          l2 = sqrt(dot_product(OP,OP))
          if (dbgCurve > 0) then
             write(dbgCurve,"('unitM   :',1P3E12.3)") unitM
             write(dbgCurve,"('theta, phiHalf :',1P2E12.3)") theta, phiHalf
             write(dbgCurve,"('l1, l2, lenSec :',1P3E12.3)") l1, l2, lenSec
          end if
#endif
          Xsi = (theta+phiHalf)/phi
          if (abs(theta) <= phiHalf) then
             isInDomain = .true.

             !! Define curve position
             tanVec = cross_product(OPproj,nVec) ! Z-axis in pos tangent dir
             vec2 = Pos1G(:,2)*(1.0_dp-Xsi) + Pos2G(:,2)*Xsi
             vec1 = cross_product(vec2,tanVec)
             PosOnCurveG(:,1) = vec1/sqrt(dot_product(vec1,vec1))
             PosOnCurveG(:,3) = tanVec/sqrt(dot_product(tanVec,tanVec))
             PosOnCurveG(:,2) = cross_product(PosOnCurveG(:,3),PosOnCurveG(:,1))
             PosOnCurveG(:,4) = origo + OPproj*R/len
          end if

       else
          Xsi = -99.9999_dp
       end if

    end if

    if (isInDomain) then

       !! Position of follower in the curve coordinate system
       Def(1:3) = matmul34(invert34(PosOnCurveG),Pos(:,4))
       if (size(Def) >= 6) then
          call FFa_eulerZYX (PosOnCurveG(:,1:3),Pos(:,1:3),Def(4:6))
       end if

       if (radius > 0.0_dp) then
          !! Assume circular contact, transform to polar coordinates
          R      = sqrt(Def(1)*Def(1)+Def(2)*Def(2))
          theta  = atan2(Def(2),Def(1))
          Def(1) = R
          Def(2) = theta
          call rotateAboutZ (PosOnCurveG,theta)
          !! Check that the radial/angular position is within parameter domain
          isInDomain = R < radius .and. &
               &      (abs(theta) < width*pi_p/360.0_dp .or. width <= 0.0_dp)
#ifdef FT_DEBUG
          if (dbgCurve > 0) then
             write(dbgCurve,"('Within endpoints,  Xsi = ',F7.4)") Xsi
             write(dbgCurve,"('rDef R :',1P2E12.3)") Def(1),radius
             write(dbgCurve,"('theta  :',1P2E12.3)") Def(2),width*0.5_dp
             write(dbgCurve,"('zDef   :',1P,E12.3)") Def(3)
          end if
#endif
          call adjustTheta (DefPrev(2),Def(2))
       else
          !! Assume rectangular contact
          !! Check that the normal/lateral position is within parameter domain
          isInDomain = abs(Def(1)) < thickness*0.5_dp .and. &
               &      (abs(Def(2)) < width*0.5_dp .or. width <= 0.0_dp)
#ifdef FT_DEBUG
          if (dbgCurve > 0) then
             write(dbgCurve,"('Within endpoints,  Xsi = ',F7.4)") Xsi
             write(dbgCurve,"('xDef t/2 :',1P2E12.3)") Def(1),thickness*0.5_dp
             write(dbgCurve,"('yDef w/2 :',1P2E12.3)") Def(2),width*0.5_dp
             write(dbgCurve,"('zDef     :',1P,E12.3)") Def(3)
          end if
#endif
       end if

    end if

#ifdef FT_DEBUG
    if (dbgCurve > 0) then
       write(dbgCurve,"('isInDomain :',L2,' Xsi =',1PE12.5)") isInDomain, Xsi
    end if
#endif

  end subroutine GetPosOnCircSection


  !!============================================================================
  !> @brief Checks if @a Pos is in a corner section of a curve.
  !>
  !> @details Compute position data if the point is in a corner.
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 2001

  subroutine GetPosCornerSection (Pos, DefPrev, CP1, &
       &                          thickness, width, radius, &
       &                          isInDomain, PosOnCurveG, Def)

    use CurvePointTypeModule, only : CurvePointType
#ifdef FT_DEBUG
    use dbgUnitsModule   , only : dbgCurve
#endif
    use KindModule       , only : dp, pi_p
    use RotationModule   , only : FFa_eulerZYX
    use ManipMatrixModule, only : cross_product, matmul34

    real(dp)            , intent(in)  :: Pos(:,:), DefPrev(3)
    type(CurvePointType), intent(in)  :: CP1
    real(dp)            , intent(in)  :: thickness, width, radius
    logical             , intent(out) :: isInDomain
    real(dp)            , intent(out) :: PosOnCurveG(:,:), Def(:)

    !! Local variables
    real(dp) :: Pos1G(3,4), rVec(3), rVecProj(3), zVec(3), d1, d2, R, theta
    real(dp), parameter :: EPS = 1.0e-30_dp

    !! --- Logic section ---

    isInDomain = .false.

    Pos1G = matmul34(CP1%triad%ur,CP1%CPosInT)
#ifdef FT_DEBUG
    if (dbgCurve > 0) write(dbgCurve,"(6X,'Point 1:',1P3E13.5)") Pos1G(:,4)
#endif

    rVec     = Pos(:,4) - Pos1G(:,4)
    Def(2)   = dot_product(rVec,Pos1G(:,2))
    rVecProj = rVec - Def(2)*Pos1G(:,2)
    Def(1)   = sqrt(dot_product(rVecProj,rVecProj))
    Def(3)   = 0.0_dp

    if (radius > 0.0_dp) then
       !! Assume circular contact, transform to polar coordinates
       R     = sqrt(Def(1)*Def(1)+Def(2)+Def(2))
       theta = atan2(Def(2),Def(1))
#ifdef FT_DEBUG
       if (dbgCurve > 0) then
          write(dbgCurve,"('insideRad, R, theta: ',2L6,1P2E11.3)") &
               (R<=radius .or. Def(1)<=0.0_dp), (abs(Def(2))<=radius), R, theta
       end if
#endif
       if (R > radius .and. Def(1) > 0.0_dp)                      return
       if (abs(Def(2)) > radius)                                  return
       if (abs(theta) > width*pi_p/360.0_dp .and. width > 0.0_dp) return
    else
       !! Assume rectangular contact
       R = 0.0_dp
       theta = 0.0_dp
#ifdef FT_DEBUG
       if (dbgCurve > 0) then
          write(dbgCurve,"('insideThick/Width, xDef, yDef: ',2L6,1P2E11.3)") &
               (Def(1)<=thickness*0.5_dp), (abs(Def(2))<=width*0.5_dp), Def(1:2)
       end if
#endif
       if (Def(1) > thickness*0.5_dp)                       return
       if (abs(Def(2)) > width*0.5_dp .and. width > 0.0_dp) return
    end if

    d1 = -dot_product(rVec,CP1%tan1G) ! To the left of right curve section
    d2 =  dot_product(rVec,CP1%tan2G) ! To the right of left curve section

#ifdef FT_DEBUG
    if (dbgCurve > 0) then
       write(dbgCurve,"('unitY: ',1P3E11.3)") Pos1G(:,2)
       write(dbgCurve,"('tan1G: ',1P4E11.3)") CP1%tan1G,sqrt(dot_product(CP1%tan1G,CP1%tan1G))
       write(dbgCurve,"('tan2G: ',1P4E11.3)") CP1%tan2G,sqrt(dot_product(CP1%tan2G,CP1%tan2G))
       write(dbgCurve,"('alpha: ',1P,E11.3)") dot_product(CP1%tan1G,CP1%tan2G)
       write(dbgCurve,"('inside1, inside2: ',2L6,1P2E11.3)") (d1>=0.0_dp),(d2>=0.0_dp),d1,d2
    end if
#endif
    if (d1 < 0.0_dp .or. d2 < 0.0_dp) return

    isInDomain = .true.
    PosOnCurveG(:,4) = Pos1G(:,4)             ! Position
    PosOnCurveG(:,2) = Pos1G(:,2)             ! Y-axis
    zVec             = CP1%tan1G + CP1%tan2G  ! Approximate Z-axis
    if (Def(1) > EPS) then
       PosOnCurveG(:,1) = rVecProj / Def(1)   ! X-axis
       PosOnCurveG(:,3) = cross_product(PosOnCurveG(:,1),PosOnCurveG(:,2)) ! Z-axis

       if (dot_product(zVec,PosOnCurveG(:,3)) < 0.0_dp) then
#ifdef bjorn
          if (dbgCurve > 0) write(dbgCurve,"(' switching X and Z directions')")
#endif
          Def(1) = -Def(1)
          PosOnCurveG(:,1) = -PosOnCurveG(:,1)
          PosOnCurveG(:,3) = -PosOnCurveG(:,3)
       end if
    else
       PosOnCurveG(:,3) = zVec / sqrt(dot_product(zVec,zVec))              ! Z-axis
       PosOnCurveG(:,1) = cross_product(PosOnCurveG(:,2),PosOnCurveG(:,3)) ! X-axis
    end if

    if (radius > 0.0_dp) then
       Def(1) = R
       Def(2) = theta
       call rotateAboutZ (PosOnCurveG,theta)
       call adjustTheta (DefPrev(2),Def(2))
    end if

    if (size(Def) >= 6) then
       call FFa_eulerZYX (PosOnCurveG(:,1:3),Pos(:,1:3),Def(4:6))
    end if

  end subroutine GetPosCornerSection


  !!============================================================================
  !> @brief Transforms the given position matrix by rotating it the given angle
  !> @a theta about its local Z-axis.
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date Sep 2005

  subroutine rotateAboutZ (PosMat,theta)

    use KindModule, only : dp

    real(dp), intent(inout) :: PosMat(:,:)
    real(dp), intent(in)    :: theta

    !! Local variables
    real(dp) :: cosTh, sinTh, vec1(3), vec2(3)

    !! --- Logic section ---

    cosTh = cos(theta)
    sinTh = sin(theta)

    vec1 = cosTh*PosMat(:,1) + sinTh*PosMat(:,2)
    vec2 = cosTh*PosMat(:,2) - sinTh*PosMat(:,1)

    PosMat(:,1) = vec1
    PosMat(:,2) = vec2

  end subroutine rotateAboutZ


  !!============================================================================
  !> @brief Adjust an angular quantity such that it is within the range
  !> [-2&pi;,2&pi;].
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date Oct 2005

  subroutine adjustTheta (oldTheta,newTheta)

    use KindModule, only : dp, pi_p

    real(dp), intent(in)    :: oldTheta
    real(dp), intent(inout) :: newTheta

    !! Local variables
    real(dp) :: theta

    !! --- Logic section ---

    theta = oldTheta
    do while (abs(theta) > 2.0_dp*pi_p)
       theta = theta - sign(2.0_dp,oldTheta)*pi_p
    end do
    if (newTheta*theta < 0.0_dp .and. abs(newTheta) > 0.5_dp*pi_p) then
       newTheta = oldTheta + newTheta-theta - sign(2.0_dp,newTheta)*pi_p
    else
       newTheta = oldTheta + newTheta-theta
    end if

  end subroutine adjustTheta


  !!============================================================================
  !> @brief Finds the position of a point on a cubic curve in space.
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 4 Nov 2016

  subroutine Cubic3D (PNod, S, PPt, tanVec, ShapeFunc, dirVec)

    use KindModule, only : dp

    real(dp), intent(in)           :: PNod(:,:)
    real(dp), intent(in), optional :: dirVec(3)
    real(dp), intent(inout)        :: S, PPt(3)
    real(dp), intent(out)          :: tanVec(3), ShapeFunc(:)

    !! Local variables
    real(dp) :: dVec(3), gVec(3), a(3), b(3), c(3), d(3), dS
    integer  :: i, order

    !! --- Logic section ---

    order = size(PNod,2) - 1
    if (order >= 3) then
       a = (PNod(:,4) - 3.0_dp*PNod(:,3) + 3.0_dp*PNod(:,2) - PNod(:,1))/6.0_dp
    else
       a = 0.0_dp
    end if
    if (order >= 2) then
       b = (PNod(:,3) - 2.0_dp*PNod(:,2) + PNod(:,1))*0.5_dp - 3.0_dp*a
    else
       b = 0.0_dp
    end if
    c = PNod(:,2) - PNod(:,1) - a - b
    d = PNod(:,1)

    !! Solve for position on curve with Newton iterations
    do i = 1, 10
       dVec = PPt - posVec(S)
       gVec = gradVec(S)
       if (present(dirVec)) then
          dS = dot_product(dVec,dirVec)/dot_product(gVec,dirVec)
       else
          dS = dot_product(dVec,gVec)/dot_product(gVec,gVec)
       end if
       S = S + dS
       if (abs(ds) < 1.0e-15_dp) exit
    end do

    PPt    = posVec(S)
    tanVec = gradVec(S)
    tanVec = tanVec/sqrt(dot_product(tanVec,tanVec))

    if (order <= 1) then
       ShapeFunc(1) = 1.0_dp - S
       ShapeFunc(2) = S
    else if (order == 2) then
       ShapeFunc(1) =  (S-1.0_dp)*(S-2.0_dp)/2.0_dp
       ShapeFunc(2) =          -S*(S-2.0_dp)
       ShapeFunc(3) =           S*(S-1.0_dp)/2.0_dp
    else
       ShapeFunc(1) = -(S-1.0_dp)*(S-2.0_dp)*(S-3.0_dp)/6.0_dp
       ShapeFunc(2) =           S*(S-2.0_dp)*(S-3.0_dp)/2.0_dp
       ShapeFunc(3) =          -S*(S-1.0_dp)*(S-3.0_dp)/2.0_dp
       ShapeFunc(4) =           S*(S-1.0_dp)*(S-2.0_dp)/6.0_dp
    end if

  contains

    !> @brief Evaluates a cubic curve at the parametric point @a s.
    function posVec(s) result(v)
      real(dp), intent(in) :: s
      real(dp) :: v(3)
      v = ((a*s + b)*s + c)*s + d
    end function posVec

    !> @brief Evaluates the gradient of a cubic curve at the point @a s.
    function gradVec(s) result(v)
      real(dp), intent(in) :: s
      real(dp) :: v(3)
      v = (3.0_dp*a*s + 2.0_dp*b)*s + c
    end function gradVec

  end subroutine Cubic3D

end module CurveTypeModule
