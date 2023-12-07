!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!! Test program for various rotation matrix operations.
!!
!! $ rotationTest <angleX> <angleY> <angleZ>
!!
!! Written by: Knut M. Okstad 14.02.2018

program test

  use kindModule    , only : dp
  use rotationModule, only : eulerToMat, matToEulerXYZ, ffa_glbeulerzyx

  implicit none

  integer           :: i, nargs
  real(dp)          :: angle(3), Tmat(3,3)
  character(len=16) :: carg

  nargs = iargc()
  angle = 0.0_dp
  if (nargs > 0) then
     do i = 1, nargs
        call getarg(i,carg)
        if (i <= 3) then
           read(carg,*) angle(i)
        end if
     end do
     print "('Angles:',3F10.6)",angle
  else
     call usage('<angleX> <angleY> <angleZ>')
  end if

  Tmat = eulerToMat(angle)
  print "('Tmat  :',3F10.6 /(7X,3F10.6))",(Tmat(i,:),i=1,3)

  angle = matToEulerXYZ(Tmat,(/(0.0_dp,i=1,3)/),i)
  print "(/'matToEulerXYZ  :',3F10.6)",angle

  call ffa_glbeulerzyx (Tmat,angle)
  print "('ffa_glbeulerzyx:',3F10.6)",angle

end program test
