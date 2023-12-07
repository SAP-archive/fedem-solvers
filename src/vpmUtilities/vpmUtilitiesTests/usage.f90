!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!! Prints a usage message to console and exits.
!!
!! Written by: Knut M. Okstad 26.07.2016

subroutine usage (args)

  character(len=*), intent(in) :: args
  character(len=128) :: prgnam
  integer :: i

  call getarg(0,prgnam)
  i = index(prgnam,'/',.true.)
  if (i < 1) i = index(prgnam,'\',.true.)
  if (i > 0) prgnam = prgnam(i+1:)
  print*,'Usage: ',trim(prgnam),' ',trim(args)
  stop

end subroutine usage
