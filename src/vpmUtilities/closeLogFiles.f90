!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file closeLogFiles.f90
!> @brief Global subroutine (callable from C++) used by signal handlers, etc.

!!==============================================================================
!> @brief Global subroutine (callable from C++) used by signal handlers, etc.
!>
!> @details This subroutine explicitly closes the log-files before exiting,
!> such that non-flushed data are not lost in case of program exceptions.
!>
!> @author Knut Morten Okstad
!>
!> @date 13 Aug 2019

subroutine closeLogFiles ()

  use allocationModule , only : closeMemLogFile
  use reportErrorModule, only : getErrorFile, setErrorFile
  use progressModule   , only : lterm

  call closeMemLogFile()
  close(getErrorFile())
  call setErrorFile(6)
  if (lterm > 6) then
     close(lterm)
  else
     call flush(lterm)
  end if

end subroutine closeLogFiles
