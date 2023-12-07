!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file prescribedMotionModule.f90
!> @brief Prescribed motions from file.

!!==============================================================================
!> @brief Module with support for prescribed motions from file.
!> @details This module contains an array with prescribed displacement values
!> which are read from a binary file. Subroutines for opening and reading this
!> file is also provided.
!>
!> @author Knut Morten Okstad
!> @date 8 Apr 2020

module PrescribedMotionModule

  use kindModule, only : dp

  implicit none

  integer , save :: pDisFile = -11111 !< File handle for prescribed motion file
  real(dp), save :: nextTime = 0.0_dp !< Time for next update from file

  real(dp), allocatable, save :: pDispl(:) !< Prescribed displacements from file


contains

  !!============================================================================
  !> @brief Opens the binary file with prescribed nodal displacements.
  !>
  !> @param[in] fileName Name of binary file with prescribed displacements
  !> @param[out] pdNodes Node numbers for points with prescribed displacements
  !> @param[out] ncmp Number of displacement components per node on file
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph

  subroutine OpenMotionFile (fileName,pdNodes,ncmp,ierr)

    use BinaryDBInterface, only : read_p, openBinaryDB, closeBinaryDB
    use BinaryDBInterface, only : readTagDB, readIntDB, readDoubleDB
    use ProgressModule   , only : lterm
    use ReportErrorModule, only : warning_p, error_p, debugFileOnly_p
    use ReportErrorModule, only : reportError, AllocationError

    character(len=*), intent(in)  :: fileName
    integer, pointer, intent(out) :: pdNodes(:)
    integer         , intent(out) :: ncmp, ierr

    !! Local variables
    integer           :: istep, nnod
    character(len=32) :: ctag

    !! --- Logic section ---

    call openBinaryDB (fileName,read_p,pDisFile,ierr)
    if (ierr < 0) then
       call reportError (error_p,'Failed to open binary file '//fileName, &
            &            addString='OpenMotionFile')
       return
    end if

    call readTagDB (pDisFile,ctag,istep,ierr)
    if (ierr < 0) then
       call reportError (error_p,'Failed to read file tag from '//fileName, &
            &            addString='OpenMotionFile')
       return
    else if (ctag /= '#FEDEM nodal displacements') then
       call reportError (error_p,trim(fileName)// &
            ' is not a nodal displacement file','tag='//trim(ctag))
       ierr = -99
       return
    else
       nnod = ierr/10
       ncmp = mod(ierr,10)
       write(lterm,100) nnod,ncmp
100    format(5X,'Number of nodes/dofs in prescribed motion file:',I6,I3)
    end if

    allocate(pdNodes(nnod),pDispl(ncmp*nnod),STAT=ierr)
    if (ierr /= 0) then
       ierr = AllocationError('OpenMotionFile')
       return
    end if

    !! Read the node numbers
    call readIntDB (pDisFile,pdNodes(1),nnod,ierr)
    if (ierr > 0) then
       !! Check for time step data on the prescribed motion file
       call readIntDB (pDisFile,istep,-1,ierr)
       if (ierr > 0) then
          !! Yes, read time for first motion update
          call readDoubleDB (pDisFile,nextTime,1,ierr)
       else ! No valid data on this file
          call reportError (warning_p,'No data on file '//fileName)
          call closeBinaryDB (pDisFile,ierr)
          pDisFile = -12345
       end if
    end if

    if (ierr < 0) call reportError (debugFileOnly_p,'OpenMotionFile')
    if (ierr < 0 .or. pDisFile < 0) then
       deallocate(pdNodes,pDispl)
       nullify(pdNodes)
    else if (ierr > 0) then
       write(lterm,110) istep,nextTime
110    format(5X,'First time step:',I6,1PE12.5)
    end if

  end subroutine OpenMotionFile


  !!============================================================================
  !> @brief Reads prescribed displacements for next time step from binary file.
  !>
  !> @param[out] ierr Error flag
  !> @param[in] lpu File unit number for debug output of motions from file
  !>
  !> @callgraph @callergraph

  subroutine ReadMotionFile (ierr,lpu)

    use BinaryDBInterface, only : readIntDB, readDoubleDB, closeBinaryDB
    use ReportErrorModule, only : note_p, error_p, reportError

    integer,           intent(out) :: ierr
    integer, optional, intent(in)  :: lpu

    !! Local variables
    integer :: lerr, nextStep
    logical :: notice

    !! --- Logic section ---

    ierr = 0
    lerr = 0
    notice = .false.
    if (pDisFile >= 0 .and. allocated(pDispl)) then
       if (present(lpu)) notice = lpu > 0
       if (notice) then
          call reportError (note_p,'Update prescribed displacements from file')
       end if
       call readDoubleDB (pDisFile,pDispl(1),size(pDispl),ierr)
    end if
    if (ierr > 0) then

       !! Check if more data on prescribed motion file
       call readIntDB (pDisFile,nextStep,-1,lerr)
       if (lerr > 0) then
          !! Yes, read time for next motion update
          call readDoubleDB (pDisFile,nextTime,1,lerr)
       else
          !! No more prescribed motions on the file
          if (notice) then
             call reportError (note_p,'End of displacement file, closing it')
          end if
          call closeBinaryDB (pDisFile,lerr)
          pDisFile = -12345
       end if
       if (lerr < 0) ierr = lerr

    end if
    if (ierr < 0) then
       call reportError (error_p,'Failed to read prescribed motions file', &
            &            addString='ReadMotionFile')
    else if (lerr > 0 .and. notice) then
       write(lpu,110) nextStep,nextTime
110    format(10X,'Next time step on prescribed motion file:',I6,1PE12.5)
    end if

  end subroutine ReadMotionFile

end module PrescribedMotionModule
