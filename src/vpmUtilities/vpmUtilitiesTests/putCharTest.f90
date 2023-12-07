!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!! Test program for putCharDB.
!! Inserts UsedTime = hh:mm:ss.dd into given RDB-file.
!!
!! $ putCharTest <filname>
!!
!! Written by: Knut M. Okstad 12.11.2002

program test

  use timerModule
  use binaryDBInterface

  implicit none

  integer           :: ierr, ifile, hours, minutes, seconds
  real              :: usedTime, hundreds
  character(len=12) :: chUsed
  real, external    :: cpusec
  character(len=64) :: fileName

  if (iargc() > 0) then
     call getarg(1,fileName)
  else
     call usage('<filename>')
  end if

  call openBinaryDB (fileName,reawri_p,ifile,ierr)
  if (ierr /= 0) stop 'can not open RBD file'

  usedTime = 1234567.89
  hours    = int(usedTime)/3600
  minutes  = int(usedTime)/60 - 60*hours
  seconds  = int(usedTime) - 60*minutes - 3600*hours
  hundreds = usedTime - int(usedTime)
  write(chUsed,"(i3,':',i2.2,':',i2.2,f3.2)") hours,minutes,seconds,hundreds
  call putCharDB (ifile,'UsedTime                =',chUsed,ierr)
  if (ierr < 0) then
     print*,'Could not insert used CPU time into RDB-file '//trim(fileName)
     stop 'putCharDB'
  end if
  call closeBinaryDB (ifile,ierr)
  if (ierr /= 0) stop 'can not close RDB-file'

end program test
