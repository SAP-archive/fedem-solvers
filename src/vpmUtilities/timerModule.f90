!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file timerModule.f90
!> @brief Subroutines for time measurement of computations.

!!==============================================================================
!> @brief Module with subroutines for time measurement of computations.
!> @details This module contains various utility subroutines and functions
!> for meassurement of CPU- and Wall time of computational tasks.
!> You have to configure doxygen with the option ENABLED_SECTIONS = FULL_DOC
!> to extract the detailed documentation of the symbols in this module.
!>
!> @author Knut Morten Okstad
!> @date 10 Jan 2003

module TimerModule
  !> @cond FULL_DOC

  use KindModule, only : sp, nbi_p, nbs_p

  implicit none

  !> @brief Timing data for a module.
  type TimerType
     logical  :: isStarted
     real(sp) :: timeAtStart(2), totalTime(2)
  end type TimerType

  type(TimerType), allocatable, save, private :: timer(:) !< Module timing data

  !> Number of bytes in a #timermodule::timertype object.
  integer, parameter, private :: nbt_p = nbi_p + 4*nbs_p


contains

  !!============================================================================
  !> @brief Initializes the timer module.
  subroutine initTime (nModules)

    use AllocationModule      , only : doLogMem, logAllocMem
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getbool

    integer, intent(in), optional :: nModules
    integer                       :: i

    !! --- Logic section ---

    call CLKINI ! Initialize the wall-time clock

    if (present(nModules)) then
       allocate(timer(nModules))
    else
       allocate(timer(1))
    end if

    call ffa_cmdlinearg_getbool ('profile',doLogMem)
    if (doLogMem) call logAllocMem ('initTime',0,size(timer),nbt_p)

    do i = 1, size(timer)
       timer(i)%isStarted = .false.
       timer(i)%totalTime = 0.0_sp
    end do

    call startTimer (1)

  end subroutine initTime


  !!============================================================================
  !> @brief Starts timer for the specified module.
  subroutine startTimer (iMod)

    integer , intent(in) :: iMod
    real(sp), external   :: CLKSEC, CPUSEC

    !! --- Logic section ---

    if (iMod > 0 .and. iMod <= size(timer)) then
       if (.not. timer(iMod)%isStarted) then
          timer(iMod)%isStarted = .true.
          timer(iMod)%timeAtStart(1) = CLKSEC(0.0_sp)
          timer(iMod)%timeAtStart(2) = CPUSEC(0.0_sp)
       end if
    end if

  end subroutine startTimer


  !!============================================================================
  !> @brief Stops timer for the specified module.
  subroutine stopTimer (iMod)

    integer , intent(in) :: iMod
    real(sp), external   :: CLKSEC, CPUSEC

    !! --- Logic section ---

    if (iMod > 0 .and. iMod <= size(timer)) then
       if (timer(iMod)%isStarted) then
          timer(iMod)%isStarted = .false.
          timer(iMod)%totalTime(1) = timer(iMod)%totalTime(1) &
               &                   + CLKSEC(timer(iMod)%timeAtStart(1))
          timer(iMod)%totalTime(2) = timer(iMod)%totalTime(2) &
               &                   + CPUSEC(timer(iMod)%timeAtStart(2))
       end if
    end if

  end subroutine stopTimer


  !!============================================================================
  !> @brief Returns the total accumulated time for the specified module.
  function getTotalTime (iMod)

    integer, intent(in) :: iMod
    real(sp)            :: getTotalTime(2)

    !! --- Logic section ---

    call stopTimer (iMod)

    if (iMod > 0 .and. iMod <= size(timer)) then
       getTotalTime = timer(iMod)%totalTime
    else
       getTotalTime = 0.0_sp
    end if

  end function getTotalTime


  !!============================================================================
  !> @brief Returns time since timer start for the specified module.
  function getCurrentTime (iMod)

    integer , intent(in) :: iMod
    real(sp)             :: getCurrentTime(2)
    real(sp), external   :: CLKSEC, CPUSEC

    !! --- Logic section ---

    getCurrentTime = 0.0_sp

    if (iMod > 0 .and. iMod <= size(timer)) then
       if (timer(iMod)%isStarted) then
          getCurrentTime(1) = CLKSEC(timer(iMod)%timeAtStart(1))
          getCurrentTime(2) = CPUSEC(timer(iMod)%timeAtStart(2))
       end if
    end if

  end function getCurrentTime


  !!============================================================================
  !> @brief Prints out elapsed- and CPU-time and peak memory usage to unit lpu.
  subroutine showTime (lpu)

    use AllocationModule    , only : doLogMem, logAllocMem, writePeakMem
    use FFaProfilerInterface, only : ffa_getMemUsage

    integer, intent(in) :: lpu

    !! Local variables
    integer  :: IDAY, IHOUR, IMIN, ISEC, JSEC
    real(sp) :: RSEC, USAGE(4)

    !! --- Logic section ---

    call stopTimer (1)

    ISEC  = int(timer(1)%totalTime(1))
    RSEC  = timer(1)%totalTime(1) - real(ISEC,sp)
    IDAY  =  ISEC/86400
    IHOUR = (ISEC-(IDAY*86400))/3600
    IMIN  = (ISEC-(IDAY*86400)-(IHOUR*3600))/60
    JSEC  = (ISEC-(IDAY*86400)-(IHOUR*3600)-(IMIN*60))
    WRITE(LPU,600) IDAY,IHOUR,IMIN,JSEC,RSEC

    ISEC  = int(timer(1)%totalTime(2))
    RSEC  = timer(1)%totalTime(2) - real(ISEC,sp)
    IDAY  =  ISEC/86400
    IHOUR = (ISEC-(IDAY*86400))/3600
    IMIN  = (ISEC-(IDAY*86400)-(IHOUR*3600))/60
    JSEC  = (ISEC-(IDAY*86400)-(IHOUR*3600)-(IMIN*60))
    WRITE(LPU,610) IDAY,IHOUR,IMIN,JSEC,RSEC

    if (doLogMem) call logAllocMem ('initTime',size(timer),0,nbt_p)
    deallocate(timer)

    if (doLogMem) then
       call writePeakMem (LPU)
    else
       call ffa_getMemUsage (usage)
       if (usage(3) > 0.0_sp) write(LPU,620) usage(3)
       if (usage(4) > 0.0_sp) write(LPU,621) usage(4)
    end if
    write(LPU,630)

600 format( /4X,'-----------------------------------', &
         &  /4X,'Elapsed time :',I4,' days ',2(I2.2,':'),I2.2,F3.2 )
610 format(  4X,'CPU time     :',I4,' days ',2(I2.2,':'),I2.2,F3.2 )
620 format(  4x,'Peak work memory usage:',F9.2,' MB' )
621 format(  4x,'Peak page memory usage:',F9.2,' MB' )
630 format(  4X,'-----------------------------------' )

  end subroutine showTime

  !> @endcond
end module TimerModule
