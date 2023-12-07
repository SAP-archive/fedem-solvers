!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file wrimsg.f90
!> @brief Global subroutine (callable from C++) to output an error message.

!!==============================================================================
!> @brief Global subroutine (callable from C++) to output an error message.
!>
!> @param[in] ifile Flag telling where to write the message (see below).
!> @param[in] msg The error message to write
!> @param[in] newLine If @e .true., write newline character at the end
!>
!> @details
!> - @a ifile &lt; 0: Write to console only.
!> - @a ifile = 0: Write to console only, or to error file if no console output.
!> - @a ifile = 1: Write to error file only.
!> - @a ifile = 2: Write to both console and error file.

recursive subroutine wrimsg (ifile,msg,newLine)

  use KindModule       , only : i4
  use ReportErrorModule, only : getErrorFile

  implicit none

  integer(kind=i4), intent(in) :: ifile
  character(len=*), intent(in) :: msg
  logical(kind=i4), intent(in) :: newLine

  !! Local variables
  integer :: iunit
#if defined(irix) || defined(irix64)
  integer :: istat
#endif
  logical :: opened

  !! --- Logic section ---

  !! Fetch the file unit number
  iunit = ifile
  iunit = getErrorFile(iunit)
  if (iunit < 0 .and. ifile == 0) iunit = getErrorFile(1)

  !! Ensure that the file has been opened
  if (iunit == 0 .or. iunit == 6) then
     opened = .true. ! Using fortran unit 0 or 6 (=> console)
  else if (iunit > 0) then
     inquire(iunit,OPENED=opened)
  else
     opened = .false.
  end if

  if (opened) then
     !! Now write the message
     if (newLine .or. len(msg) < 1) then
        write(iunit,'(A)') msg
     else
        write(iunit,'(A)',advance='NO') msg
     end if
#if defined(irix) || defined(irix64)
     call flush(iunit,istat)
#else
     call flush(iunit)
#endif
  else if (ifile == 1) then
     !! The error file is not opened yet so write the message to console instead
     call wrimsg(-1,msg,newLine)
  end if

  !! Write the message to console also
  if (ifile == 2) call wrimsg(-1,msg,newLine)

end subroutine wrimsg
