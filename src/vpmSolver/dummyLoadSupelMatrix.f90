!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

subroutine LoadSupelMatrixFromCore (supId,dtype,data,ndata,stat)
  !DEC$ ATTRIBUTES DLLEXPORT :: LoadSupelMatrixFromCore

  integer         , intent(in)  :: supId, ndata
  character(len=1), intent(in)  :: dtype
  double precision, intent(out) :: data(ndata)
  integer         , intent(out) :: stat

  !! Dummy implementation, a real model-specific instance needs to be generated
  !! using the subroutine initiateSupElTypeModule::writeSupEls2Ftn
  if (supId < 0) then
     print *,' ** LoadSupelMatrixFromCore dummy: ',dtype,supId,ndata
  end if
  data(1) = 0.0D0
  stat = 0

end subroutine LoadSupelMatrixFromCore


subroutine LoadStrainRosetteFromCore (supId,iros,data,ndata,stat)
  !DEC$ ATTRIBUTES DLLEXPORT :: LoadStrainRosetteFromCore

  integer, intent(in)  :: supId, iros, ndata
  real   , intent(out) :: data(ndata)
  integer, intent(out) :: stat

  !! Dummy implementation, a real model-specific instance needs to be generated
  !! using the subroutine stressRecoveryModule::writeRosettes2Ftn
  if (supId < 0) then
     print *,' ** LoadStrainRosetteFromCore dummy: ',supId,iros,ndata
  end if
  data(1) = 0.0
  stat = 0

end subroutine LoadStrainRosetteFromCore
