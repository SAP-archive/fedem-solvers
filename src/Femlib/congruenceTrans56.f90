!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file congruenceTrans56.f90
!> @brief Special congruence transformation for thick shell element matrices.

!!==============================================================================
!> @brief Performs a congruence transformation of a symmetric matrix.
!>
!> @details This subroutine performs a congruence transformation of a symmetric
!> (stiffness) matrix with dimension (5*N)x(5*N) where N is the number of nodes.
!> Only the rotational dofs (local dof 4 and 5) are transformed, and such that
!> the transformed matrix gets the dimension (6*N)x(6*N).
!> The transformation matrix consists of 3x3 submatrices along the diagonal.
!>
!> The input array EK5 and output array EK6 can be the same physical address
!> (as long N is larger than 2) since we are looping backward over the matrices.
!> With only two nodes this will yield incorrect result as parts of the
!> input will be overwritten before it is used. Therefore, if the subroutine is
!> used for two-noded elements, EK5 and EK6 must be different addresses.
!>
!> @author Knut Morten Okstad
!>
!> @date 26 Jan 2021

subroutine congruenceTrans56 (EK5,EK6,C,N,EPS6)

  use kindModule, only : dp

  implicit none

  integer , intent(in)  :: N
  real(dp), intent(in)  :: EK5(5,N,5,N), C(N,3,3), EPS6
  real(dp), intent(out) :: EK6(6,N,6,N)

  !! Local variables
  integer  :: I, IA, J, JA, K, L
  real(dp) :: B(3,3), fictS6

  !! ---- Logic section ----

  !! Find a proper fictitious drilling DOF stiffness
  fictS6 = 0.0_dp
  if (EPS6 > 0.0_dp) then
     do I = 1, N
        fictS6 = fictS6 + EK5(4,I,4,I) + EK5(5,I,5,I)
     end do
     fictS6 = EPS6 * fictS6/real(2*N,dp)
  end if

  do I = N, 1, -1 ! Loop over nodes (columns)
     do K = 3, 0, -3 ! Loop over rotation and translation (columns)
        do J = I, 1, -1 ! Loop over nodes (rows)
           do L = 3, 0, -3 ! Loop over rotation and translation (rows)

              B(:,3) = 0.0_dp
              if (L == 3) then
                 do IA = 1, 3-K/3
                    B(:,IA) = C(J,1,:)*EK5(4,J,K+IA,I) &
                         &  + C(J,2,:)*EK5(5,J,K+IA,I)
                 end do
                 if (K == 3 .and. I == J .and. fictS6 > 0.0_dp) then
                    B(:,3) = C(J,3,:)*fictS6
                 end if
              else
                 do IA = 1, 3-K/3
                    B(:,IA) = EK5(1:3,J,K+IA,I)
                 end do
              end if

              if (K == 3) then
                 do IA = 1, 3
                    do JA = 1, 3
                       EK6(L+JA,J,K+IA,I) = B(JA,1)*C(I,1,IA) &
                            &             + B(JA,2)*C(I,2,IA) &
                            &             + B(JA,3)*C(I,3,IA)
                    end do
                 end do
              else
                 do IA = 1, 3
                    do JA = 1, 3
                       EK6(L+JA,J,K+IA,I) = B(JA,IA)
                    end do
                 end do
              end if

           end do
        end do
     end do
  end do

  !! Fill in nodal sub-matrices in lower triangle
  do I = 1, N
     do J = 1, I-1
        do K = 1, 6
           do L = 1, 6
              EK6(K,I,L,J) = EK6(L,J,K,I)
           end do
        end do
     end do
  end do

end subroutine congruenceTrans56
