!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

program linkTest

  !! This test program is only used to test for undefined symbols on linking.

  use CompositeTypeModule, only : dp, CompositeType, Cmatrix
  use andes3ShellModule  , only : Andes3shell_stiffmat
  use andes4ShellModule  , only : Andes4shell_stiffmat
  use cstetraModule      , only : cstetstiff, cstetmass, cstetpload, cstetstress
  use ipri6Module        , only : ipri6stiff, ipri6mass, ipri6pload, ipri6stress

  integer             :: ierr
  real(dp)            :: alpha, beta, Xl(6), Yl(6), Zl(6), dx(3,6), dy(3,6)
  real(dp)            :: CcMat(6,6), C(6,6), EK(24,24), disp(18)
  type(CompositeType) :: laminate

  call Cmatrix (laminate,beta,C)
  call Andes3shell_stiffmat (Xl,Yl,C,alpha,beta,2,EK,IERR)
  call Andes4shell_stiffmat (Xl,Yl,Zl,C,alpha,beta,2,EK,IERR)

  call cstetstiff (xl,yl,zl,CcMat,EK,alpha,6,0,ierr)
  call cstetmass (xl,yl,zl,beta,EK,alpha,6,0,ierr)
  call cstetpload (xl,yl,zl,dx,1,dy,6,0,ierr)
  call cstetstress (xl,yl,zl,CcMat,disp,dx(1,:),alpha,6,0,ierr)

  call ipri6stiff (xl,yl,zl,CcMat,EK,1,6,0,ierr)
  call ipri6mass (xl,yl,zl,beta,EK,1,6,0,ierr)
  call ipri6pload (xl,yl,zl,dx,1,dy,6,0,ierr)
  call ipri6stress (xl,yl,zl,CcMat,disp,C,1,2,6,0,ierr)

  call beam31
  call beam35

  call fqs31
  call fqs32
  call fqs33
  call fqs35

  call fts31
  call fts33
  call fts35
  call fts38
  call ftsa31
  call ftsa32

  call hexa31
  call hexa32
  call hexa33
  call hexa35

  call ihex31
  call ihex32
  call ihex33
  call ihex35

  call ipri31
  call ipri32
  call ipri33
  call ipri35

  call itet31
  call itet32
  call itet33
  call itet35

  call scqs31
  call scqs32
  call scqs33
  call scqs35

  call scts31
  call scts32
  call scts33
  call scts35

end program linkTest
