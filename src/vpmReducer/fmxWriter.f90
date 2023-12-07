!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file fmxWriter.f90
!> @brief Shared library wrapper for reading and writing FEDEM fmx files.
!> @author Knut Morten Okstad, SAP SE
!> @date 22 Nov 2019

!> @brief Module encapsulation of the IO_FXM subroutine.

module FMXwriter

  implicit none

contains

  !> @brief Read/writes a double precision array from/to the specified fmx-file.
  !> @param[in] prefix File name prefix
  !> @param[in] itype Matrix type (1=stiffness, 2=mass, 3=gravity forces)
  !> @param data Array with matrix content
  !> @param[in] nval Size of the data array
  !> @param[in] chkSum Write to file if present, otherwise read from file
  !> @param[out] ierr Error flag (negative value indicates error, otherwise OK)

  subroutine io_FMX (prefix,itype,data,nval,chkSum,ierr)

    use kindModule       , only : dp, lfnam_p
    use binaryDBInterface, only : writeDoubleDB, readDoubleDB

    character(len=*) , intent(in)    :: prefix
    integer          , intent(in)    :: itype, nval
    real(dp)         , intent(inout) :: data(*)
    integer, optional, intent(in)    :: chkSum
    integer          , intent(out)   :: ierr

    character(len=lfnam_p) :: fileName
    character(len=32)      :: fileTag

    select case (itype)
    case (1)
       fileName = trim(prefix)//'_K.fmx'
       fileTag  = 'stiffness matrix'
    case (2)
       fileName = trim(prefix)//'_M.fmx'
       fileTag  = 'mass matrix'
    case (3)
       fileName = trim(prefix)//'_G.fmx'
       fileTag  = 'gravity force vectors'
    case default
       print *,'*** FMX: Unknown file type',itype
       ierr = -max(1,itype)
       return
    end select

    ierr = 0
    if (present(chkSum)) then
       print *,'  * Writing '//trim(fileTag)//' file: '//trim(fileName),nval
       call writeDoubleDB (fileName,fileTag,chkSum, &
            &              reshape(data(1:nval),(/nval,1/)),ierr)
    else
       print *,'  * Reading '//trim(fileTag)//' file: '//trim(fileName),nval
       call readDoubleDB (fileName,fileTag,nval,data(1),ierr)
    end if
    if (ierr < 0) print *,'*** Failure:',ierr

  end subroutine io_FMX

end module FMXwriter


!> @brief Writes a double precision array to the specified fmx-file.
!> @param[in] prefix File name prefix
!> @param[in] itype Matrix type (1=stiffness, 2=mass, 3=gravity forces)
!> @param[in] data Array with matrix content
!> @param[in] nval Size of the data array
!> @return Error status. Negative value indicates an error, otherwise OK

function writeFMX (prefix,itype,data,nval) result(ierr)
!DEC$ ATTRIBUTES DLLEXPORT :: WRITEFMX
  use FMXwriter, only : io_FMX
  character(len=*), intent(in)    :: prefix
  integer         , intent(in)    :: itype, nval
  double precision, intent(inout) :: data(*)
  integer         , parameter     :: chkSum = 0
  integer                         :: ierr
  call io_FMX (prefix,itype,data,nval,chkSum,ierr)
end function writeFMX


!> @brief Reads a double precision array from the specified fmx-file.
!> @param[in] prefix File name prefix
!> @param[in] itype Matrix type (1=stiffness, 2=mass, 3=gravity forces)
!> @param[out] data Array with matrix content
!> @param[in] nval Size of the data array
!> @return Error status. Negative value indicates an error, otherwise OK

function readFMX (prefix,itype,data,nval) result(ierr)
!DEC$ ATTRIBUTES DLLEXPORT :: READFMX
  use FMXwriter, only : io_FMX
  character(len=*), intent(in)    :: prefix
  integer         , intent(in)    :: itype, nval
  double precision, intent(inout) :: data(*)
  integer                         :: ierr
  call io_FMX (prefix,itype,data,nval,ierr=ierr)
end function readFMX
