!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file pFUnitExtraArg.f90
!> @brief This file contains additional source code for the pFUnit driver.
!> @details This is an extension of the Fortran unit test system (pFUnit)
!> to enable passing additional command-line options to the unit test
!> executables when invoked via ctest.
!>
!> @author Knut Morten Okstad / SAP SE
!> @date 10 Apr 2020

!> @brief Module with additional source code for the pFUnit driver.
!> @details This module contains a function for extracting the directory
!> of the source code for the unit tests, which then can be used when
!> loading specified data files required by the tests. This way the unit
!> test executables can be invoked from an arbitrary working directory.
!>
!> The pf-files need to be compiled with -DPFUNIT_EXTRA_USAGE=pfUnitArgs
!> and -DPFUNIT_EXTRA_ARGS=get_srcdir to activate this code.

module pFUnitArgs

  !> Source directory of the unit tests
  character(len=:), allocatable, save :: srcdir

contains

  !> @brief Extracts the source directory from a command-line argument.
  function get_srcdir (argument)
    character(len=*), intent(in) :: argument
    logical :: get_srcdir
    integer :: len_arg
    len_arg = len_trim(argument)
    if (len_arg < 10) then
       get_srcdir = .false.
    else if (argument(1:9) == '--srcdir=') then
       get_srcdir = .true.
       allocate(character(len=len_arg-9) :: srcdir)
       srcdir = argument(10:)
    else
       get_srcdir = .false.
    end if
  end function get_srcdir

end module pFUnitArgs
