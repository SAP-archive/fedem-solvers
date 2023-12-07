!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module ioModule

  implicit none

  integer, save :: IO = 0 ! File unit number for error/warning messages

contains

  subroutine initIO

    !!==========================================================================
    !! Initializes the file unit number for the solver res-file.
    !!==========================================================================
#ifdef FT_SHARED
    !! When the control module is built as a shared object library (DLL), we
    !! cannot just pass the file unit as an argument from the calling program,
    !! as it is not recognised as an open unit inside the shared library (we
    !! would then get a fort.<unit> file instead). We must check all currently
    !! open file units, for one with a matching name (should only be one .res).
    logical :: named
    character*256 :: fName
    if (IO == 0) then
       IO = 9
       named = .true.
       fName = '(none)'
       do while (named)
          if (fName(len_trim(fName)-3:) == '.res') return
          IO = IO + 1
          inquire(IO,NAMED=named,NAME=fName)
       end do
       if (.not. named) then
          !! Did not detect any open res-file, using ControlSystem.res instead
          open(IO,FILE='ControlSystem.res',STATUS='UNKNOWN',ERR=1000)
          return
1000      IO = 6 ! Opening failure, just use console output
       end if
    end if
#else
    !! Static build, just get it from reportErrorModule
    use reportErrorModule, only : getErrorFile
    if (IO == 0) IO = getErrorFile()
#endif

  end subroutine initIO

end module ioModule
