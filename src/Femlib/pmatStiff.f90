!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module pmatModule

  !! The pmatStiff subroutine is now wrapped in a module, to avoid crash due to
  !! the assumed-size definition of the arguments (bugfix #477, kmo 22.06.2016).

  implicit none

  interface pmatStiff
     module procedure pmat_stiff1, pmat_stiff2
  end interface

  private :: pmat_stiff1, pmat_stiff2

contains

subroutine pmat_stiff1 (coor,pmat,lpu,ierr)

  !!****************************************************************************
  !! Computes the linear deformational projector of an element.
  !! That is, the projector that filters out the linear rigid body modes
  !! from the displacement vector.
  !!
  !! Input:
  !!    coor : coor(i,j) is coordinate component j of node i
  !!    lpu  : print unit number for error messages
  !!
  !! Output:
  !!    pmat : symmetric projector matrix that filters out the rigid body part
  !!           of a displacement vector.
  !!           dof ordering: tx, ty, tz, rx, ry, rz for each node
  !!    ierr : error flag
  !!
  !! Coded by: Bjorn Haugen
  !! To Fortran90 by: Knut Morten Okstad, 11.03.2016
  !!****************************************************************************

  use kindModule       , only : dp
  use manipMatrixModule, only : invert33

  real(dp), intent(in)  :: coor(:,:)
  real(dp), intent(out) :: pmat(:,:)
  integer , intent(in)  :: lpu
  integer , intent(out) :: ierr

  !! Local variables
  integer  :: i, j, nnod
  real(dp) :: c, coorRel(3,size(coor,2)), rmat(size(coor,2)*6,6)
  real(dp) :: rsmat(6), sub(3,3), subinv(3,3)

  !! ---- Logic section ----

  !! Compute nodal coordinates relative to centroid of element (nodal average)

  nnod = size(coor,2)
  do i = 1, 3
     c = sum(coor(i,:))/real(nnod,dp)
     coorRel(i,:) = coor(i,:) - c
  end do

  !! Initialize rigid modes matrix

  rmat = 0.0_dp
  do i = 1, nnod
     j = (i-1)*6

     !! Rigid body translations in x, y, and z direction
     rmat(j+1,1) =  1.0_dp
     rmat(j+2,2) =  1.0_dp
     rmat(j+3,3) =  1.0_dp

     !! Rigid body rotation about x-axis
     rmat(j+2,4) = -coorRel(3,i)
     rmat(j+3,4) =  coorRel(2,i)
     rmat(j+4,4) =  1.0_dp

     !! Rigid body rotation about y-axis
     rmat(j+1,5) =  coorRel(3,i)
     rmat(j+3,5) = -coorRel(1,i)
     rmat(j+5,5) =  1.0_dp

     !! Rigid body rotation about z-axis
     rmat(j+1,6) = -coorRel(2,i)
     rmat(j+2,6) =  coorRel(1,i)
     rmat(j+6,6) =  1.0_dp

  end do

  !! Normalize modal vectors

  do i = 1, 6
     c = sqrt(dot_product(rmat(:,i),rmat(:,i)))
     rmat(:,i) = rmat(:,i)/c
  end do

  !! Compute nondiagonal part of R'*R

  do i = 1, 3
     do j = 1, 3
        sub(i,j) = dot_product(rmat(:,i+3),rmat(:,j+3))
     end do
  end do
  subinv = invert33(sub,lpu,ierr)
  if (ierr < 0) return

  !! Compute RS = R*Inverse(R'*R)

  do i = 1, 6*nnod
     rsmat(1:3) = rmat(i,1:3)
     do j = 1, 3
        rsmat(3+j) = dot_product(rmat(i,4:6),subinv(1:3,j))
     end do

     !! Compute the projector matrix as P = ( I - R*(RR)R' )
     do j = 1, 6*nnod
        pmat(i,j) = -dot_product(rsmat,rmat(j,:))
     end do
     pmat(i,i) = pmat(i,i) + 1.0_dp
  end do

end subroutine pmat_stiff1


subroutine pmat_stiff2 (x,y,z,pmat,lpu,ierr)
  use kindModule, only : dp
  real(dp), intent(in)  :: x(:), y(:), z(:)
  real(dp), intent(out) :: pmat(:,:)
  integer , intent(in)  :: lpu
  integer , intent(out) :: ierr
  real(dp) :: xyz(3,size(x))
  xyz(1,:) = x
  xyz(2,:) = y
  xyz(3,:) = z
  call pmat_stiff1 (xyz,pmat,lpu,ierr)
end subroutine pmat_stiff2

end module pmatModule
