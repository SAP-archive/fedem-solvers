!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module ResStressModule

  !! Dummy module used when stress recovery routines are linked with the solver.

  implicit none


contains

  subroutine addResStress (sam,Sigma,IEL,NENOD,NSTRP,IPRINT,LPU,IERR)
    use SamModule, only : SamType, dp
    type(SamType), intent(in)    :: sam
    real(dp)     , intent(inout) :: Sigma(:)
    integer      , intent(in)    :: IEL, IPRINT, LPU
    integer      , intent(inout) :: NENOD, NSTRP, IERR
    if (iprint > 0) then
       write(lpu,*) ' ** addResStress dummy:',sam%nel,IEL,NENOD, &
            &       NSTRP,size(Sigma),IERR
    end if
  end subroutine addResStress

end module ResStressModule
