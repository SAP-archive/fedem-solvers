!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!! Test program for writeBinaryDB.
!! Writes a file test.dat containing <nCount> doubles.
!!
!! $ writeBinaryTest <nCount> [<bufsize>]
!!
!! Written by: Knut M. Okstad 12.03.2002

program test

  use kindModule       , only : i4, i8, dp, nbd_p
  use binaryDBInterface, only : openBinaryDB, closeBinaryDB
  use binaryDBInterface, only : writeDoubleDB, readDoubleDB
  use binaryDBInterface, only : append_p, write_p, read_p

  implicit none

  real(dp), allocatable :: data(:)
  real(dp) :: dev, error
  integer(i8) :: i, j, nBytes, nCount, bufsize
  integer(i4) :: ierr, ifile
  character(len=16) :: cArg

  if (iargc() > 0) then
     call getarg(1,cArg)
     read(cArg,*) nCount
     print*,'Number of words to write',nCount
     if (nCount < 100) stop 'too few < 100'
  else
     call usage('<nCount> [<bufsize>]')
  end if

  if (iargc() > 1) then
     call getarg(2,cArg)
     read(cArg,*) bufSize
     if (bufSize < 1000_i8) then
        bufSize = 1000_i8
     else if (bufSize > nCount) then
        bufSize = nCount
     end if
  else
     bufSize = 1000_i8
  end if

  print*,'Using buffer size       ',bufSize
  allocate(data(bufSize),stat=ierr)
  if (ierr /= 0) stop 'allocation error'

  call openBinaryDB ('test.dat',append_p,ifile,ierr)
  if (ierr /= 0) call openBinaryDB ('test.dat',write_p,ifile,ierr)
  if (ierr /= 0) stop 'can not open test.dat'

  nBytes = 0_i8
  do i = 1_i8, bufSize
     data(i) = real(i,dp)
  end do
  do i = 1_i8, nCount, bufSize
     call printProgress (i)
     call writeDoubleDB (ifile,data(1),MIN(bufSize,nCount+1_i8-i),ierr)
     if (ierr < 0) then
        print*,'Failure after writing',nBytes,'bytes'
        print*,'I =',i,' IERR =',ierr
        stop
     else
        nBytes = nBytes + int(ierr,i8)
     end if
  end do

  print*,'Successfully wrote',nBytes,' bytes'
  call closeBinaryDB (ifile,ierr)
  if (ierr /= 0) stop 'can not close test.dat'

  print*,'Reading test.dat'
  call openBinaryDB ('test.dat',read_p,ifile,ierr)
  if (ierr /= 0) stop 'can not open test.dat'

  nBytes = 0_i8
  error = 0.0_dp
  do i = 1_i8, nCount, bufSize
     call printProgress (i)
     call readDoubleDB (ifile,data(1),MIN(bufSize,nCount+1_i8-i),ierr)
     if (ierr < 0) then
        print*,'Failure after reading',nBytes,'bytes'
        print*,'I =',i,' IERR =',ierr
        stop
     else
        nBytes = nBytes + int(ierr,i8)
     end if
     do j = 1, ierr/nbd_p
        dev = data(j) - real(j,dp)
        if (abs(dev) > 1.0e-16_dp) then
           error = error + abs(dev)
           print*,'Deviation at',j,dev
        end if
     end do
  end do

  print*,'Successfully read ',nBytes,' bytes'
  call closeBinaryDB (ifile,ierr)
  if (ierr /= 0) stop 'can not close test.dat'

  deallocate(data)

  if (error > 0.0_dp) then
     print*,'Accumulated/mean error',error,error/real(nCount,dp)
  end if

contains

  subroutine printProgress(count)
    integer(i8), intent(in) :: count
    integer(i8) :: nMod
    if (count <= 1_i8) then
       return
    else if (nCount > 100000_i8) then
       nMod = nCount/100_i8
    else
       nMod = nCount/10_i8
    end if
    if (MOD(count,nMod) < MOD(count-bufSize,nMod)) then
       print"(2PF5.1,'%')",REAL(count-1_i8,dp)/REAL(nCount,dp)
    end if
  end subroutine printProgress

end program test
