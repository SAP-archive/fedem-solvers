!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file rotationModule.f90
!> @brief Utilities for manipulation of finite rotations.

!!==============================================================================
!> @brief Module with subroutines for manipulation of finite rotations.
!> @details This module contains various utility subroutines (and functions)
!> for manipulation of finite rotations in the nonlinear solver, like convertion
!> of rotation tensors to angles, and vice versa.
!>
!> @brief Bjorn Haugen
!> @date February 1999

module RotationModule

  use KindModule, only : dp

  implicit none

  real(dp), parameter, private :: epsTh2_p = 0.0005_dp !< Angles zero-tolerance

  !> @brief Takes into account possible eccentricity.
  interface EccExpand
     module procedure EccExpandVec
     module procedure EccExpandMat
     module procedure EccExpandDiagMat
  end interface

  private :: EccExpandVec, EccExpandMat, EccExpandDiagMat


  interface

     !> @brief Computes the EulerZYX angles for a relative rotation.
     !> @param[in] a The rotation tensor to rotate from
     !> @param[in] b The rotation tensor to rotate to
     !> @param[out] angles The resulting Euler angles
     subroutine ffa_eulerzyx (a,b,angles)
       use KindModule, only  : dp
       real(dp), intent(in)  :: a(3,3), b(3,3)
       real(dp), intent(out) :: angles(3)
     end subroutine ffa_eulerzyx

     !> @brief Computes the EulerZYX angles for a global rotation.
     !> @param[in] a The rotation tensor to rotate to (from global inertia axes)
     !> @param[out] angles The resulting Euler angles
     subroutine ffa_glbeulerzyx (a,angles)
       use KindModule, only  : dp
       real(dp), intent(in)  :: a(3,3)
       real(dp), intent(out) :: angles(3)
     end subroutine ffa_glbeulerzyx

  end interface


contains

  !!============================================================================
  !> @brief Converts a rotation matrix (tensor) to EulerXYZ angles.
  !> @param[in] R The rotation matrix to convert
  !> @param[in] eulerIn Euler angle set to use as a hint
  !> @param[out] err Error flag (=1: singularity encountered, infinite
  !> solutions for Z and X rotation. = -1: Failed to find a valid solution,
  !> try orthonormalizing the coordinate system)
  !> @return The EulerXYZ rotation angles
  !>
  !> @details This function tries to find the Euler angle set that is closest
  !> to the set given as input.
  !>
  !> @note This function is no longer in use, consider remove?

  function matToEulerXYZ (R,eulerIn,err) result(euler)

    use KindModule       , only : pi_p
    use ReportErrorModule, only : reportError, warning_p, error_p

    real(dp), intent(in)  :: R(3,3), eulerIn(3)
    integer , intent(out) :: err

    !! Local Variables
    integer  :: i, bestSet, n(3)
    real(dp) :: a11, a12 , a21, a22 !! alfa
    real(dp) :: b1, b2              !! beta
    real(dp) :: g11, g12, g21, g22  !! gamma
    real(dp) :: eps, eps2, sum, bestSum
    real(dp) :: euler(3)

    !! --- Logic section ---

    err  = 0
    eps  = sqrt(sqrt(epsilon(1.0_dp)))    ! Verification limit
    eps2 = 0.001_dp*sqrt(epsilon(1.0_dp)) ! To remove singularity

    do i = 1, 3
       n(i) = int(EulerIn(i)*0.5_dp/pi_p)
       if (EulerIn(i) < 0.0_dp) n(i) = n(i) - 1
       Euler(i) = EulerIn(i) - 2.0_dp*pi_p*n(i)
    end do

    if ( abs(R(3,1)) < 1.0_dp - eps2 ) then

       !! Finding beta
       b1  = adjust2Pi(-asin(R(3,1)))
       b2  = adjust2Pi(pi_p - b1)

       !! Finding alfa
       a11 = adjust2Pi(asin(R(3,2)/cos(b1)))
       a12 = adjust2Pi(pi_p - a11)

       a21 = adjust2Pi(asin(R(3,2)/cos(b2)))
       a22 = adjust2Pi(pi_p - a21)

       !! Finding gamma
       g11 = adjust2Pi(asin(R(2,1)/cos(b1)))
       g12 = adjust2Pi(pi_p - g11)

       g21 = adjust2Pi(asin(R(2,1)/cos(b2)))
       g22 = adjust2Pi(pi_p - g21)

       !! check all sets and return the which is closest

       !! We do now have eight solution sets
       !! set 1:  a11 b1 g11
       !! set 1:  a11 b1 g12
       !! set 1:  a12 b1 g11
       !! set 1:  a12 b1 g12
       !! set 1:  a21 b2 g21
       !! set 1:  a21 b2 g22
       !! set 1:  a22 b2 g21
       !! set 1:  a22 b2 g22
       !! We then finds the best legal set and use it

       bestSet = 0
       bestSum = 100.0_dp
       if     (SetisOk(a11,b1,g11)) then
          sum  =  dist(a11,b1,g11)
          if (sum < bestSum) then; bestSum = sum; bestSet = 1
          endif
       endif
       if (SetisOk(a11,b1,g12)) then
          sum  =  dist(a11,b1,g12)
          if (sum < bestSum) then; bestSum = sum; bestSet = 2
          endif
       endif
       if (SetisOk(a12,b1,g11)) then
          sum  =  dist(a12,b1,g11)
          if (sum < bestSum) then; bestSum = sum; bestSet = 3
          endif
       endif
       if(SetisOk(a12,b1,g12)) then
          sum  =  dist(a12,b1,g12)
          if (sum < bestSum) then; bestSum = sum; bestSet = 4
          endif
       endif
       if (SetisOk(a21,b2,g21)) then
          sum  =  dist(a21,b2,g21)
          if (sum < bestSum) then; bestSum = sum; bestSet = 5
          endif
       endif
       if(SetisOk(a21,b2,g22)) then
          sum  =  dist(a21,b2,g22)
          if (sum < bestSum) then; bestSum = sum; bestSet = 6
          endif
       endif
       if (SetisOk(a22,b2,g21)) then
          sum  =  dist(a22,b2,g21)
          if (sum < bestSum) then; bestSum = sum; bestSet = 7
          endif
       endif
       if (SetisOk(a22,b2,g22)) then
          sum  =  dist(a22,b2,g22)
          if (sum < bestSum) then; bestSum = sum; bestSet = 8
          endif
       endif

       select case(BestSet)

       case(0); err = -1
          call reportError (error_p,'Failed to calculate euler angles', &
               &                    'returning previous values')
       case(1);  euler(1) = a11;  euler(2) = b1; euler(3)=g11
       case(2);  euler(1) = a11;  euler(2) = b1; euler(3)=g12
       case(3);  euler(1) = a12;  euler(2) = b1; euler(3)=g11
       case(4);  euler(1) = a12;  euler(2) = b1; euler(3)=g12
       case(5);  euler(1) = a21;  euler(2) = b2; euler(3)=g21
       case(6);  euler(1) = a21;  euler(2) = b2; euler(3)=g22
       case(7);  euler(1) = a22;  euler(2) = b2; euler(3)=g21
       case(8);  euler(1) = a22;  euler(2) = b2; euler(3)=g22

       end select

    else if ( R(3,1) < 0) then

       err = 1
       call reportError (warning_p, &
            'Singularity encountered when calculating euler angles.', &
            'Infinite number of solutions for x and z rotations.', &
            'returning previous/input z-value')

       a11 = adjust2Pi(asin(R(1,2)) + euler(3))
       a22 = adjust2Pi(pi_p - a11 + 2.0_dp*euler(3))
       if      (abs(R(1,3) - cos(a11 - euler(3))) < eps) then
          euler(1) = a11
          euler(2) = pi_p/2.0_dp
       else if (abs(R(1,3) - cos(a22 - euler(3))) < eps) then
          euler(1) = a22
          euler(2) = pi_p/2.0_dp
       else
          err = -1
          call reportError (error_p,'Failed to calculate euler angles', &
               &                    'returning previous/input values')
       endif

    else ! R(2,1) < 0

       err = 1
       call reportError (warning_p, &
            'Singularity encountered when calculating euler angles.', &
            'Infinite number of solutions for x and z rotations.', &
            'returning previous/input z-value')

       a11 = adjust2Pi(-asin(R(1,2)) - euler(3))
       a22 = adjust2Pi(pi_p - a11 - 2.0_dp*euler(3))
       if      (abs(R(1,3) + cos(a11 + euler(3))) < eps) then
          euler(1) =  a11
          euler(2) = -pi_p/2.0_dp
       else if (abs(R(1,3) + cos(a22 + euler(3))) < eps) then
          euler(1) =  a22
          euler(2) = -pi_p/2.0_dp
       else
          err = -1
          call reportError (error_p,'Failed to calculate euler angles', &
               &                    'returning previous/input values')
       endif

    endif

    !! Return the angle with the correct number of rotations
    do i = 1, 3
       if ( EulerIn(i) > 0.5_dp*(4*n(i)+3)*pi_p .and. &
            Euler(i)   < 0.5_dp*pi_p ) then
          n(i) = n(i) + 1
       else if ( EulerIn(i) < 0.5_dp*(4*n(i)+1)*pi_p .and. &
            &    Euler(i)   > 1.5_dp*pi_p) then
          n(i) = n(i) - 1
       end if
       if (n(i) /= 0) Euler(i) = Euler(i) + 2*n(i)*pi_p
    end do

  contains

    !> @brief Adjusts an angle such that it is in the interval [0,2&pi;].
    function adjust2Pi (a)
      real(dp), intent(in) :: a
      real(dp)             :: adjust2Pi
      real(dp), parameter  :: twoPi_p = 2.0_dp*pi_p
      if (a < 0.0_dp) then
         adjust2Pi = a + twoPi_p
      else if (a > twoPi_p) then
         adjust2Pi = a - twoPi_p
      else
         adjust2Pi = a
      end if
    end function adjust2Pi

    !> @brief Calculates the distance between the angles a,b,g and euler.
    function dist (a,b,g)
      real(dp), intent(in) :: a, b, g
      real(dp)             :: dist
      if      (a > 1.5_dp*pi_p .and. euler(1) < 0.5_dp*pi_p) then
         dist = euler(1) - (a-2.0*pi_p)
      else if (a < 0.5_dp*pi_p .and. euler(1) > 1.5_dp*pi_p) then
         dist = (a+2.0*pi_p) - euler(1)
      else
         dist = abs(euler(1)-a)
      end if
      if      (b > 1.5_dp*pi_p .and. euler(2) < 0.5_dp*pi_p) then
         dist = dist + euler(2)- (b-2.0*pi_p)
      else if (b < 0.5_dp*pi_p .and. euler(2) > 1.5_dp*pi_p) then
         dist = dist +(b+2.0*pi_p) - euler(2)
      else
         dist = dist + abs(euler(2)-b)
      end if
      if      (g > 1.5_dp*pi_p .and. euler(3) < 0.5_dp*pi_p) then
         dist = dist + euler(3) - (g-2.0*pi_p)
      else if (g < 0.5_dp*pi_p .and. euler(3) > 1.5_dp*pi_p) then
         dist = dist + (g+2.0*pi_p) - euler(3)
      else
         dist = dist + abs(euler(3)-g)
      end if
    end function dist

    !> @brief Checks the consistency of an Euler angle set.
    function SetIsOk (a,b,g)
      real(dp), intent(in) :: a, b, g
      logical              :: setIsOk
      SetIsOk = .false.
      if (abs(cos(b)*cos(g)                      - R(1,1)) > eps) return
      if (abs(sin(a)*sin(b)*cos(g)-cos(a)*sin(g) - R(1,2)) > eps) return
      if (abs(cos(a)*sin(b)*cos(g)+sin(a)*sin(g) - R(1,3)) > eps) return
      if (abs(sin(a)*sin(b)*sin(g)+cos(a)*cos(g) - R(2,2)) > eps) return
      if (abs(cos(a)*sin(b)*sin(g)-sin(a)*cos(g) - R(2,3)) > eps) return
      if (abs(cos(a)*cos(b)                      - R(3,3)) > eps) return
      SetIsOk = .true.
    end function SetIsOk

  end function matToEulerXYZ


  !!============================================================================
  !> @brief Converts a set of EulerXYZ angles into a rotation matrix.

  function eulerToMat (euler) result(Tr)

    real(dp)            :: Tr(3,3)
    real(dp),intent(in) :: euler(3)

    !! Local variables
    real(dp) :: cx,sx,cy,sy,cz,sz

    !! --- Logic section ---

    cx = cos(euler(1))
    sx = sin(euler(1))
    cy = cos(euler(2))
    sy = sin(euler(2))
    cz = cos(euler(3))
    sz = sin(euler(3))

    tr(1,1) =  cy*cz
    tr(2,1) =  cy*sz
    tr(3,1) = -sy
    tr(1,2) = -cx*sz + cz*sx*sy
    tr(2,2) =  cx*cz + sx*sy*sz
    tr(3,2) =  cy*sx
    tr(1,3) =  cx*cz*sy + sx*sz
    tr(2,3) = -cz*sx + cx*sy*sz
    tr(3,3) =  cx*cy

  end function eulerToMat


  !!============================================================================
  !> @brief Computes the rotation tensor from a rotation vector.
  !> @param[in]  rvec The rotation vector
  !> @param[out] rten The rotation tensor
  !> @param[in]  lpu Optional file unit number for error messages

  subroutine vec_to_mat (rvec,rten,lpu)

    real(dp), intent(in)  :: rvec(3)
    real(dp), intent(out) :: rten(3,3)
    integer , intent(in), optional :: lpu

    real(dp) :: q(4) !< quaternion representation of the rotation

    call vec_to_quat (rvec,q,lpu)
    call quat_to_mat (q,rten)

  end subroutine vec_to_mat


  !!============================================================================
  !> @brief Computes the rotation vector from a rotation tensor.
  !> @param[in]  rten The rotation tensor
  !> @param[out] rvec The rotation vector

  subroutine mat_to_vec (rten,rvec)

    real(dp), intent(in)  :: rten(3,3)
    real(dp), intent(out) :: rvec(3)

    real(dp) :: q(4) !< quaternion representation of the rotation

    call mat_to_quat (rten,q)
    call quat_to_vec (q,rvec)

  end subroutine mat_to_vec


  !!============================================================================
  !> @brief Computes the quaternion representation from a rotation vector.
  !> @param[in]  rvec The rotation vector
  !> @param[out] q The quaternion representation of the rotation
  !> @param[in]  lpu Optional file unit number for error messages

  subroutine vec_to_quat (rvec,q,lpu)

    real(dp), intent(in)  :: rvec(3)
    real(dp), intent(out) :: q(4)
    integer , intent(in), optional :: lpu

    real(dp) :: thh, sthh, cthh, f1, fac

    thh  = 0.5_dp*sqrt(sum(rvec*rvec))
    sthh = sin(thh)
    cthh = cos(thh)
    if (thh > 1.0e6_dp) then
       if (abs(1.0_dp - cthh*cthh - sthh*sthh) > 0.00001_dp) then
          if (present(lpu)) then
             write(lpu,"(A,1PE12.5,A)") ' *** vec_to_quat: theta = ', thh+thh, &
                  &                      ' too large for sin/cos evaluation'
          end if
          q(1)  = 1.0_dp
          q(2:) = 0.0_dp
          return
       end if
    end if

    if (thh < epsTh2_p) then
       f1  = thh/epsTh2_p
       fac = f1*sin(epsTh2_p)/epsTh2_p + 1.0_dp - f1
    else
       fac = sthh/thh
    end if

    q(1)   = cthh
    q(2:4) = rvec*fac*0.5_dp

    q = q / sqrt(sum(q*q))

  end subroutine vec_to_quat


  !!============================================================================
  !> @brief Computes the quaternion representation from a rotation tensor.
  !> @param[in]  rten The rotation tensor
  !> @param[out] q The rotation quaternions (Euler parameters)
  !>
  !> @details The computation is performed by means of equations from
  !> Richard A. Spurrier comments on "Singularity free Extraction of a
  !> Quaternion from a Direction-Cosine Matrix"
  !> in Journal of Rockets and Spacecraft, 1978.

  subroutine mat_to_quat (rten,q)

    real(dp), intent(in)  :: rten(3,3)
    real(dp), intent(out) :: q(4)

    integer  :: imax, i, j, k
    real(dp) :: trace

    trace = rten(1,1) +rten(2,2) +rten(3,3)

    imax = 1
    if ( rten(2,2) > rten(imax,imax) ) imax = 2
    if ( rten(3,3) > rten(imax,imax) ) imax = 3

    if ( trace > rten(imax,imax) )  then
       q(1) = sqrt(1.0_dp + trace)*0.5_dp
       q(2) = ( rten(3,2) -rten(2,3) )/(4.0_dp*q(1))
       q(3) = ( rten(1,3) -rten(3,1) )/(4.0_dp*q(1))
       q(4) = ( rten(2,1) -rten(1,2) )/(4.0_dp*q(1))
    else
       i = imax
       j = mod(imax,3)+1
       k = mod(imax+1,3)+1
       q(i+1) = sqrt(rten(i,i)*0.5_dp + (1.0_dp-trace)*0.25_dp)
       q(  1) = ( rten(k,j) -rten(j,k) )/(4.0_dp*q(i+1))
       q(j+1) = ( rten(j,i) +rten(i,j) )/(4.0_dp*q(i+1))
       q(k+1) = ( rten(k,i) +rten(i,k) )/(4.0_dp*q(i+1))
    end if

  end subroutine mat_to_quat


  !!============================================================================
  !> @brief Computes the rotation tensor from a quaternion representation.
  !> @param[in]  q The quaternion representation of the rotation
  !> @param[out] rten The rotation tensor

  subroutine quat_to_mat (q,rten)

    real(dp), intent(inout) :: q(4)
    real(dp), intent(out)   :: rten(3,3)

    q = q / sqrt(sum(q*q))

    rten(1,1) = 2.0_dp*(q(2)*q(2) + q(1)*q(1)) - 1.0_dp
    rten(2,2) = 2.0_dp*(q(3)*q(3) + q(1)*q(1)) - 1.0_dp
    rten(3,3) = 2.0_dp*(q(4)*q(4) + q(1)*q(1)) - 1.0_dp

    rten(1,2) = 2.0_dp*(q(2)*q(3) - q(4)*q(1))
    rten(1,3) = 2.0_dp*(q(2)*q(4) + q(3)*q(1))
    rten(2,3) = 2.0_dp*(q(3)*q(4) - q(2)*q(1))

    rten(2,1) = 2.0_dp*(q(3)*q(2) + q(4)*q(1))
    rten(3,1) = 2.0_dp*(q(4)*q(2) - q(3)*q(1))
    rten(3,2) = 2.0_dp*(q(4)*q(3) + q(2)*q(1))

  end subroutine quat_to_mat


  !!============================================================================
  !> @brief Computes the rotation vector from a quaternion representation.
  !> @param[in]  q The quaternion representation of the rotation
  !> @param[out] rvec The rotation vector

  subroutine quat_to_vec(q,rvec)

    real(dp), intent(inout) :: q(4)
    real(dp), intent(out)   :: rvec(3)

    real(dp) :: thh, sthh, cthh, f1, fac

    q = q / sqrt(sum(q*q))

    cthh = q(1)
    sthh = sqrt(sum(q(2:4)*q(2:4)))
    if (sthh < 0.7_dp) then
       thh = asin(sthh)
    else
       thh = acos(cthh)
    end if

    if (thh < epsTh2_p) then
       f1  = thh/epsTh2_p
       fac = f1*epsTh2_p/sin(epsTh2_p) + 1.0_dp - f1
    else if (sthh >= 1.0) then
       fac = thh
    else
       fac = thh/sthh
    end if

    rvec = q(2:4)*fac*2.0_dp

  end subroutine quat_to_vec


  !!============================================================================
  !> @brief Computes the incremental rotation vector rotating tensor T1 to T2.

  function deltaRot (T1,T2)
    real(dp), intent(in) :: T1(3,3),T2(3,3)
    real(dp)             :: deltaRot(3)
    call mat_to_vec (matmul(T2,transpose(T1)),deltaRot)
  end function deltaRot


  !!============================================================================
  !> @brief Makes a rotation tensor orthonormal through use of quaternions.

  subroutine orthonorm3 (rten)

    real(dp), intent(inout) :: rten(3,3)

    real(dp) :: q(4) !< quaternion representation of the rotation

    call mat_to_quat (rten,q)
    call quat_to_mat (q,rten)

  end subroutine orthonorm3


  !!============================================================================
  !> @brief Computes the spin (skew-symmetric) tensor from the rotation vector.

  subroutine Spin (vec,mat)

    real(dp), intent(in)  :: vec(:)
    real(dp), intent(out) :: mat(3,3)

    mat(1,1) = 0.0_dp
    mat(2,2) = 0.0_dp
    mat(3,3) = 0.0_dp

    mat(2,3) = -vec(1)
    mat(1,3) =  vec(2)
    mat(1,2) = -vec(3)

    mat(3,2) =  vec(1)
    mat(3,1) = -vec(2)
    mat(2,1) =  vec(3)

  end subroutine Spin


  !!============================================================================
  !> @brief Computes the variation of a rotation pseudo vector.
  !>
  !> @param[in] thetaVec Initial rotation vector
  !> @param[out] jacobiMat Matrix of pseudo-vector variations with respect to
  !> `w = (eps, 0, 0)`, `w = (0, eps, 0)` and `w = (0, 0, eps)`
  !> @param[in] lpu Optional file unit number for error messages
  !>
  !> @details The variation of the rotation vector `r` is computed with respect
  !> to the incremental rotation `w`, where the rotations obey the rule
  !>
  !> R(r+w) = R(w)&lowast;R(r)
  !>
  !> The variation is given as
  !>
  !> V = I - 0.5&lowast;Spin(&theta;) +
  !> &eta;&lowast;Spin(&theta;)&lowast;Spin(&theta;)
  !>
  !> where &eta; = (sin(0.5&lowast;&theta;) -
  !>   0.5&lowast;&theta;&lowast;cos(0.5&lowast;&theta;)) &frasl;
  !> (&theta;<sup>2</sup>&lowast;sin(0.5&lowast;&theta;))

  subroutine Dtheta_Domega (thetaVec,jacobiMat,lpu)

    real(dp), intent(in)  :: thetaVec(3)
    real(dp), intent(out) :: jacobiMat(3,3)
    integer , intent(in), optional :: lpu

    real(dp) :: eta, th, spn(3,3), f2, f1
    real(dp), parameter :: epsTheta_p = 0.01_dp

    th = sqrt(dot_product(thetaVec,thetaVec))

    !! Avoid 0/0 problems for small angles
    if ( th < epsTheta_p ) then
       f1  = th/epsTheta_p
       f2  = 1.0_dp - f1
       eta = etaFunc(epsTheta_p,lpu)*f1 &  ! Value at th=epsTheta_p
            &     + (1.0_dp/12.0_dp)*f2    ! Value at th=0.0
    else
       eta = etaFunc(th,lpu)
    end if

    !! Compute the spin of the vector
    call Spin( thetaVec, spn )

    !! Compute first part of var = I - Spin/2
    jacobiMat      = -0.5_dp*spn
    jacobiMat(1,1) =  1.0_dp
    jacobiMat(2,2) =  1.0_dp
    jacobiMat(3,3) =  1.0_dp

    !! Add the term  eta*Spin*Spin
    jacobiMat = jacobiMat + matmul(spn,spn)*eta

  contains

    !> @brief Evaluates the parameter &eta;.
    function etaFunc (th,lpu) result(eta)
      real(dp), intent(in) :: th
      integer , intent(in), optional :: lpu

      real(dp) :: thh, sthh, cthh, eta

      thh  = th*0.5_dp
      sthh = sin(thh)
      cthh = cos(thh)
      if ( thh > 1.0e6_dp) then
         if ( abs(1.0_dp - cthh*cthh - sthh*sthh) > 0.00001_dp ) then
            cthh = 1.0_dp; sthh = 0.0_dp
            if (present(lpu)) then
               write(lpu,"(A,1PE12.5,A)") ' *** Dtheta_Domega: theta = ', th, &
                    &                     ' too large for sin/cos evaluation'
               call flush(lpu) ! In case we get core dump on division by zero
            end if
         end if
      end if
      eta = (sthh - thh*cthh)/( th*th*sthh )
    end function etaFunc

  end subroutine Dtheta_Domega


  !!============================================================================
  !> @brief Computes the inverse of the variation of a rotation pseudo vector.
  !>
  !> @param[in] thetaVec Initial rotation vector
  !> @param[out] jacobiMat The inverse matrix of pseudo-vector variations
  !> @param[in] lpu Optional file unit number for error messages

  subroutine Domega_Dtheta (thetaVec, jacobiMat, lpu)
    use ManipMatrixModule, only : invert33, unify
    real(dp), intent(in)  :: thetaVec(3)
    real(dp), intent(out) :: jacobiMat(3,3)
    integer , intent(in), optional :: lpu
    integer :: ierr

    call Dtheta_Domega(thetaVec,jacobiMat,lpu)
    jacobiMat = invert33(jacobiMat,lpu,ierr)
    if (ierr /= 0) call unify(jacobiMat)

  end subroutine Domega_Dtheta


  !!============================================================================
  !> @brief Takes into account possible eccentricity.
  !>
  !> @details @verbatim
  !>    | vFull | = |    I    :  0 || v |
  !>    |       |   | Spin(e) :  I || 0 |
  !> @endverbatim
  !>
  !> @author Bjorn Haugen
  !> @date Feb 2002

  subroutine EccExpandVec (e,v,vFull)

    use ManipMatrixModule, only : cross_product

    real(dp), intent(in)  :: e(3), v(3)
    real(dp), intent(out) :: vFull(6)

    !! --- Logic section ---

    vFull(1:3) = v
    vFull(4:6) = cross_product(e,v)

  end subroutine EccExpandVec


  !!============================================================================
  !> @brief Takes into account possible eccentricity.
  !>
  !> @details @verbatim
  !>   | Mfull | = |    I    :  0 || M : 0 || I : Spin(e)' |
  !>   |       |   | Spin(e) :  I || 0 : 0 || 0 :   I      |
  !> @endverbatim
  !>
  !> @author Bjorn Haugen
  !> @date Feb 2002

  subroutine EccExpandMat (e,M,Mfull)

    real(dp), intent(in)  :: e(3), M(3,3)
    real(dp), intent(out) :: Mfull(6,6)

    !! Local variables
    real(dp) :: spinE(3,3)

    !! --- Logic section ---

    call spin (e,spinE)
    Mfull(1:3,1:3) = M
    Mfull(4:6,1:3) = matmul(spinE,M)
    Mfull(1:3,4:6) = transpose(Mfull(4:6,1:3))
    Mfull(4:6,4:6) = matmul(spinE,Mfull(1:3,4:6))

  end subroutine EccExpandMat


  !!============================================================================
  !> @brief Takes into account possible eccentricity, for a diagonal matrix D.
  !>
  !> @details @verbatim
  !>   | Mfull | = |    I    :  0 || D : 0 || I : Spin(e)' |
  !>   |       |   | Spin(e) :  I || 0 : 0 || 0 :   I      |
  !> @endverbatim
  !>
  !> @author Knut Morten Okstad
  !> @date 21 Jan 2004

  subroutine EccExpandDiagMat (e,D,Mfull)

    real(dp), intent(in)  :: e(3), D(3)
    real(dp), intent(out) :: Mfull(6,6)

    !! Local variables
    integer  :: i
    real(dp) :: spinE(3,3)

    !! --- Logic section ---

    Mfull = 0.0_dp
    call spin (e,spinE)
    do i = 1, 3
       Mfull(i,i) = D(i)
       Mfull(4:6,i) = spinE(1:3,i)*D(i)
       Mfull(i,4:6) = Mfull(4:6,i)
    end do
    Mfull(4:6,4:6) = matmul(spinE,Mfull(1:3,4:6))

  end subroutine EccExpandDiagMat

end module RotationModule
