!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file versionModule.f90
!> @brief Handling of program versioning.

!!==============================================================================
!> @brief Module with data and subroutines for program version handling.

module versionModule

  implicit none

  character(len=32), save :: progVer   !< Version ID for this release
  character(len=32), save :: buildDate !< The date this program was built


contains

  !!============================================================================
  !> @brief Returns current date and time in the format "dd Mmm yyyy hh:mm:ss".
  !>
  !> @param[out] curDate String with current data and time
  !> @param[in]  nospace If .true. use format "yyyy-Mmm-dd_hhmmss" instead
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 14 Aug 2001

  subroutine getCurrentDate (curDate,nospace)

    character(len=*), intent(out) :: curDate
    logical,optional, intent(in)  :: nospace

    !! Local variables
    integer                       :: year, month, day, hour, min, sec

    character(len=3), parameter   :: months_p(12) = &
         &                        (/ 'Jan','Feb','Mar','Apr','May','Jun', &
         &                           'Jul','Aug','Sep','Oct','Nov','Dec' /)

    !! --- Logic section ---

#if defined(linux) || defined(linux64)
    call datime (year,day,month,hour,min,sec)
#else
    call datime (year,month,day,hour,min,sec)
#endif
    if (.not. present(nospace)) then
       write(curDate,600) day,months_p(month),year,hour,min,sec
    else if (nospace) then
       write(curDate,601) year,months_p(month),day,hour,min,sec
    else
       write(curDate,600) day,months_p(month),year,hour,min,sec
    end if

600 format(i2,a4,i5,i3.2,':',i2.2,':',i2.2)
601 format(i4,'-',a3,'-',i2.2,'_',3i2.2)

  end subroutine getCurrentDate


  !!============================================================================
  !> @brief Opens the solver res-file and print header information to it.
  !>
  !> @param[in] chname Name of the solver res-file
  !> @param[in] solver Name of the solver executable
  !> @param[out] lpu   File unit number of the solver res-file
  !> @param[out] ierr  Error flag
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 14 Aug 2001

  subroutine openResFile (chname,solver,lpu,ierr)

    use fileUtilitiesModule    , only : findUnitNumber
    use reportErrorModule      , only : error_p, reportError, setErrorFile
    use computerConfigInterface, only : getComputerConfig, getUserName
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_list

    character(len=*), intent(in)  :: chname, solver
    integer         , intent(out) :: lpu, ierr

    !! Local variables
    character(len=80) :: chid
    character(len=16) :: user
    character(len=32) :: date

    !! --- Logic section ---

    lpu = findUnitNumber(10)
    open(lpu,FILE=chname,STATUS='UNKNOWN',IOSTAT=ierr)
    if (ierr /= 0) then
       call reportError (error_p,'Unable to open res-file '//chname,lpu=6)
       return
    end if

    call setErrorFile (lpu)
    call getComputerConfig (chid)
    call getUserName (user)
    call getCurrentDate (date)

    write(lpu,600)
    write(lpu,601) trim(solver)
    write(lpu,602) trim(progVer),trim(buildDate)
    write(lpu,603) trim(date)
    write(lpu,604) trim(chid)
    write(lpu,605) trim(user)
    write(lpu,600)

    !! List all specified command-line options
    call ffa_cmdlinearg_list (.true.)

600 format(/5X,55('='))
601 format(/5X,'FEDEM ',A)
602 format(/5X,'Module version: ',A,2X,A)
603 format(/5X,'Execution date: ',A)
604 format(/5X,'Executed on: ',A)
605 format(/5X,'User: ',A)

  end subroutine openResFile

end module versionModule


!!==============================================================================
!> @brief Global subroutine (callable from C++) for setting program version tag.
!>
!> @param[in] theVer Program version tag
!>
!> @author Knut Morten Okstad
!>
!> @date 14 Aug 2001

subroutine setVersion (theVer)
  use versionModule, only : progVer
  character(len=*), intent(in) :: theVer
  progVer = theVer
end subroutine setVersion

!!==============================================================================
!> @brief Global subroutine (callable from C++) for setting program build date.
!>
!> @param[in] theDate Program build date
!>
!> @author Knut Morten Okstad
!>
!> @date 14 Aug 2001

subroutine setDate (theDate)
  use versionModule, only : buildDate
  character(len=*), intent(in) :: theDate
  buildDate = theDate
end subroutine setDate
