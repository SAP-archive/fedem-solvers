# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# Library setup

set ( LIB_ID vpmLinAlg_dummy )
set ( UNIT_ID ${DOMAIN_ID}_${PACKAGE_ID}_${LIB_ID} )

message ( STATUS "INFORMATION : Processing unit ${UNIT_ID}" )

# Avoid that the generated module files are moved to the common folder
unset ( CMAKE_Fortran_MODULE_DIRECTORY )

## Fortran90 source files
set ( F90_FILE_LIST dummyErrorFlag dummyFELinearSolver ../feDataModule )

foreach ( FILE ${F90_FILE_LIST} )
  list ( APPEND F90_SOURCE_FILES ${FILE}.f90 )
endforeach ( FILE ${F90_FILE_LIST} )

add_library ( ${LIB_ID} ${F90_SOURCE_FILES} )
target_link_libraries ( ${LIB_ID} vpmCommon_F90 )
