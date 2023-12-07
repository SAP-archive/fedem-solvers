!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file inputUtilities.f90
!> @brief Utilities for parsing of solver input files.

!!==============================================================================
!> @brief Module with functions and subroutines for parsing solver input files.
!> @details This module contains various utility subroutines and functions
!> used for parsing solver input files based on the namelist Fortran90 feature.
!>
!> @brief Bjorn Haugen
!> @date 31 Sep 1999

module InputUtilities

  implicit none

  private :: iuToUpper


contains

  !!============================================================================
  !> @brief Copies a solver input file to a temporary scratch file.
  !>
  !> @param[in] itmp File unit number of the file to copy to
  !> @param[in] chname Name of the solver input file to copy from
  !> @param[out] ierr Error flag
  !>
  !> @details Any `\r` character at the end of each line is removed on UNIX.
  !>
  !> @author Knut Morten Okstad
  !> @date 8 Jan 2002

  subroutine iuCopyToScratch (itmp,chname,ierr)

    use KindModule         , only : lfnam_p
    use FileUtilitiesModule, only : findUnitNumber
    use reportErrorModule  , only : reportError, error_p, debugFileOnly_p

    integer         , intent(in)  :: itmp
    character(len=*), intent(in)  :: chname
    integer         , intent(out) :: ierr

    !! Local variables
    integer                   :: infp, i, n, readHeading
    character(len=32+lfnam_p) :: chline

    !! --- Logic section ---

    infp = findUnitNumber(10)
    open(infp,FILE=chname,STATUS='OLD',ACTION='READ',IOSTAT=ierr)
    if (ierr /= 0) then
       call reportError (error_p,'Unable to open input file '//chname, &
            &            addString='iuCopyToScratch')
       return
    end if

    readHeading = 0
    ReadLines : do
       read(infp,'(a)',END=100) chline
       n = len_trim(chline)
       if (n >= len(chline)) then
          ierr = ierr - 1
          call reportError (error_p,'Line too long in input file '//chname, &
               &            chline(1:80)//' ...')
       else if (n < 1) then
          cycle ! Ignore blank lines
       end if
#if !defined(win32) && !defined(win64)
       !! Remove trailing \r on UNIX systems
       if (ichar(chline(n:n)) == 13) n = n - 1
#endif
       !! Replace invalid version field in old &HEADING namelist
       if (readHeading == 0) then
          if (chline(1:8) == '&HEADING') readHeading = 1
       else if (readHeading == 1) then
          if (chline(1:11) == '  version =') then
             if (index(chline,'_') > 0) chline(12:) = ' 2.5'
             readHeading = 2
          end if
       end if
       !! If a description string contains apostrophes, add another one to avoid
       !! read failure because apostrophes also are used as string delimiters.
       !! The added apostrophe(s) will not appear when printing the string.
       if (chline(1:12) == '  extDescr =' .and. n > 15) then
          i = n-1
          do while (i > 14)
             if (chline(i:i) == "'") then
                chline(i:n+1) = "'"//chline(i:n)
                n = n + 1
             end if
             i = i - 1
          end do
       end if
       write(itmp,'(a)') chline(1:n)
    end do ReadLines

100 continue
    close(infp)
    if (ierr < 0) call reportError (debugFileOnly_p,'iuCopyToScratch')

  end subroutine iuCopyToScratch


  !!============================================================================
  !> @brief Writes an input file to a temporary scratch file.
  !>
  !> @param[in] itmp File unit number of the file to write to
  !> @param[in] chmodel Character string containing the solver model to write
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine is used the model is passed to the solver through
  !> a long text string, and not via a solver input file. That is, when the
  !> actual solver input file has been parsed by some external process.
  !>
  !> @author Knut Morten Okstad
  !> @date 12 Dec 2016

  subroutine iuWriteToScratch (itmp,chmodel,ierr)

    use KindModule       , only : lfnam_p
    use reportErrorModule, only : reportError, error_p, debugFileOnly_p

    integer         , intent(in)  :: itmp
    character(len=*), intent(in)  :: chmodel
    integer         , intent(out) :: ierr

    !! Local variables
    integer                   :: istart, i, j, k, n, nchar
    character(len=64)         :: errMsg
    character(len=33+lfnam_p) :: chrec

    !! --- Logic section ---

    j = 0
    istart = 1
    nchar = len(chmodel)
    do i = 1, nchar
       if (i <= istart .or. chmodel(i:i) /= '$') cycle

       j = j + 1
       n = i - istart
       if (n > len(chrec)) then
          ierr = ierr - 1
          write(errMsg,"('Input record',I4,' is too long:',2I4)") j,n,len(chrec)
          call reportError (error_p,errMsg,chrec(1:80)//' ...')
          n = len(chrec)
       end if
       chrec = chmodel(istart:i-1)
       !! Add an extra space in the beginning of each data line
       if (chrec(1:1) /= '&' .and. chrec(1:1) /= '/') then
          chrec = ' '//chrec(1:n)
          n = n + 1
       end if
       istart = i + 1
       !! If a description string contains apostrophes, add another one to avoid
       !! read failure because apostrophes also are used as string delimiters.
       !! The added apostrophe(s) will not appear when printing the string.
       if (chrec(1:11) == ' extDescr =' .and. n > 14) then
          k = n-1
          do while (k > 13)
             if (chrec(k:k) == "'") then
                chrec(k:) = "'"//chrec(k:n)
                n = n + 1
             end if
             k = k - 1
          end do
       end if

       write(itmp,'(a)') chrec(1:n)
    end do

    if (ierr < 0) call reportError (debugFileOnly_p,'iuWriteToScratch')

  end subroutine iuWriteToScratch


  !!============================================================================
  !> @brief Converts the text @a string to all uppercase.

  function iuToUpper (string) result(upper)

    character(len=*), intent(in) :: string

    !! Local variables
    integer                    :: i, delta
    character(len=len(string)) :: upper

    !! --- Logic section ---

    delta = ichar('A') - ichar('a')

    do i = 1, len(string)
       if (string(i:i) >= 'a' .and. string(i:i) <= 'z') then
          upper(i:i) = char(ichar(string(i:i))+delta)
       else
          upper(i:i) = string(i:i)
       end if
    end do

    !! Shift characters left if leading blanks
    if (upper(1:1) == ' ') upper = adjustl(upper)

  end function iuToUpper


  !!============================================================================
  !> @brief Counts number of entries of a text string in a solver input file.
  !>
  !> @param[in] file File unit number of the file to count entries in
  !> @param[in] string The text string to count the number of entries of
  !> @param[in] ierr Error flag

  function iuGetNumberOfEntries (file,string,ierr) result(numEntries)

    use reportErrorModule, only : reportError, error_p

    integer         , intent(in)  :: file
    character(len=*), intent(in)  :: string
    integer         , intent(out) :: ierr

    !! Local variables
    integer            :: numEntries
    character(len=256) :: chline, entryUpper

    !! --- Logic section ---

    numEntries = 0
    entryUpper = iuToUpper(string)

    rewind(file,IOSTAT=ierr)
    do while (ierr == 0)
       read(file,'(a)',IOSTAT=ierr) chline
       if (ierr == 0 .and. iuToUpper(chline) == entryUpper) then
          numEntries = numEntries + 1
       else if (ierr > 0) then
          call reportError (error_p,'Failed to search for string "'// &
               &            trim(entryUpper)//'" in input file', &
               &            addString='iuGetNumberOfEntries')
          return
       end if
    end do

    rewind(file,IOSTAT=ierr)

  end function iuGetNumberOfEntries


  !!============================================================================
  !> @brief Sets position at the first entry of a text string in a file.
  !>
  !> @param[in] file File unit number of the file to set position for
  !> @param[in] string The text string to find the file position for
  !> @param ierr Error flag
  !> @return .true. if the @a string is found, otherwise .false.

  logical function iuSetPosAtFirstEntry (file,string,ierr)

    integer         , intent(in)  :: file
    character(len=*), intent(in)  :: string
    integer,optional, intent(out) :: ierr

    !! Local variables
    integer :: stat

    !! --- Logic section ---

    rewind(file,IOSTAT=stat)
    if (stat /= 0) then
       if (present(ierr)) ierr = stat
       iuSetPosAtFirstEntry = .false.
    else
       iuSetPosAtFirstEntry = iuSetPosAtNextEntry(file,string,ierr)
    end if

  end function iuSetPosAtFirstEntry


  !!============================================================================
  !> @brief Sets position at the next entry of a text string in a file.
  !>
  !> @param[in] file File unit number of the file to set position for
  !> @param[in] string The text string to find the file position for
  !> @param ierr Error flag
  !> @return .true. if the @a string is found, otherwise .false.

  logical function iuSetPosAtNextEntry (file,string,ierr)

    integer         , intent(in)  :: file
    character(len=*), intent(in)  :: string
    integer,optional, intent(out) :: ierr

    !! Local variables
    integer            :: stat
    character(len=256) :: chline, entryUpper

    !! --- Logic section ---

    if (present(ierr)) ierr = 0
    iuSetPosAtNextEntry = .true.

    stat = 0
    entryUpper = iuToUpper(string)

    do while (stat == 0)
       read(file,'(a)',IOSTAT=stat) chline
       if (stat == 0 .and. iuToUpper(chline) == entryUpper) then
          backspace(file,IOSTAT=stat)
          if (stat == 0) return
       end if
    end do

    if (present(ierr)) ierr = stat
    iuSetPosAtNextEntry = .false.

  end function iuSetPosAtNextEntry


  !!============================================================================
  !> @brief Searches for a text string in an array of strings.
  !>
  !> @param[in] string The text string to search for
  !> @param[in] stringArray Array of text strings to search within
  !> @param ierr Error flag
  !> @param[in] startIndex Optional offset for the returned array index
  !> @return Index of the first occurence of the @a string

  integer function iuCharToInt (string,stringArray,ierr,startIndex)

    use reportErrorModule, only : reportError, errorFileOnly_p

    character(len=*), intent(in)    :: string, stringArray(:)
    integer,optional, intent(in)    :: startIndex
    integer,optional, intent(inout) :: ierr

    !! Local variables
    integer :: i, startPos

    !! --- Logic section ---

    if (present(startIndex)) then
       startPos = startIndex-1
    else
       startPos = 0
    end if

    if (string == '') then
       iuCharToInt = startPos ! A blank string is legal (always has position 0)
       return
    end if

    do i = 1, size(stringArray)
       if (iuToUpper(string) == iuToUpper(stringArray(i))) then
          iuCharToInt = startPos + i
          return
       end if
    end do

    !! If not returned before here, an error has occured
    iuCharToInt = startPos
    if (present(ierr)) ierr = ierr - 1
    call reportError (errorFileOnly_p,'Invalid string value "'//trim(string)// &
         &                            '" encountered',addString='iuCharToInt')

  end function iuCharToInt

end module InputUtilities
