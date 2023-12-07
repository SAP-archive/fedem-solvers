!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module fatigueModule

  implicit none

contains

  subroutine fatigueInit (ierr)

    !!==========================================================================
    !! Initialize for fatigue calculations.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 30 Apr 2008/1.0
    !!==========================================================================

    use reportErrorModule  , only : reportError, error_p
    use FFpFatigueInterface, only : ffp_initFatigue

    integer, intent(out) :: ierr

    !! --- Logic section ---

    call ffp_initFatigue (ierr)
    if (ierr < 0) then
       call reportError (error_p,'Failure in fatigue initialization', &
            addString='reported by fatigueInit')
    end if

  end subroutine fatigueInit


  subroutine fatigueAddPoint (strainCoats,time)

    !!==========================================================================
    !! Accumulate time histories for the fatigue calculations.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 30 Apr 2008/1.0
    !!==========================================================================

    use kindModule         , only : dp
    use StrainCoatModule   , only : StrainCoatType
    use FFpFatigueInterface, only : ffp_addpoint

    type(StrainCoatType), intent(inout) :: strainCoats(:)
    real(dp)            , intent(in)    :: time

    !! Local variables
    integer :: i, j

    !! --- Logic section ---

    do i = 1, size(strainCoats)
       if (.not. associated(strainCoats(i)%fppProcessorHandle)) cycle

       do j = 1, size(strainCoats(i)%fppProcessorHandle)

          call ffp_addpoint (strainCoats(i)%fppProcessorHandle(j), time, &
               &             strainCoats(i)%results(j)%fatValue)

       end do
    end do

  end subroutine fatigueAddPoint


  subroutine fatigueDamage (strainCoats,gateValue,ierr)

    !!==========================================================================
    !! Perform rainflow analysis on the time histories and calculate damage.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 30 Apr 2008/1.0
    !!==========================================================================

    use kindModule         , only : dp
    use idTypeModule       , only : GetID, StrId
    use StrainCoatModule   , only : StrainCoatType
    use reportErrorModule  , only : reportError, error_p
    use FFpFatigueInterface, only : ffp_calcDamage, ffp_releaseData

    type(StrainCoatType), intent(inout) :: strainCoats(:)
    real(dp)            , intent(in)    :: gateValue
    integer             , intent(out)   :: ierr

    !! Local variables
    integer :: i, j

    !! --- Logic section ---

    ierr = 0
    do i = 1, size(strainCoats)
       if (.not. associated(strainCoats(i)%fppProcessorHandle)) cycle

       do j = 1, size(strainCoats(i)%fppProcessorHandle)

          call ffp_calcDamage (strainCoats(i)%fppProcessorHandle(j), &
               &               strainCoats(i)%results(j)%SNcurve, gateValue, &
               &               strainCoats(i)%results(j)%damage, ierr)
          if (ierr < 0) then
             call reportError (error_p,'Failure in damage calculation', &
                  'Strain coat:'//GetId(strainCoats(i)%tensorRosette%id), &
                  'ffp Handle ='//StrId(strainCoats(i)%fppProcessorHandle(j)), &
                  addString='reported by fatigueDamage')
             return
          end if

          call ffp_releaseData (strainCoats(i)%fppProcessorHandle(j))
          strainCoats(i)%fppProcessorHandle(j) = 0

       end do
    end do


  end subroutine fatigueDamage

end module fatigueModule
