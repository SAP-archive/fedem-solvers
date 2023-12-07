# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# Library setup

set ( LIB_ID JDTire )
set ( UNIT_ID ${DOMAIN_ID}_${PACKAGE_ID}_${LIB_ID} )
 
message ( STATUS "INFORMATION : Processing unit ${UNIT_ID}" )

if ( FTENV_WARNINGS AND LINUX AND NOT USE_INTEL_FORTRAN )
  string ( APPEND CMAKE_Fortran_FLAGS " -Wno-compare-reals -Wno-conversion" )
  if ( CMAKE_C_COMPILER_VERSION VERSION_GREATER 6.0 )
    string ( APPEND CMAKE_Fortran_FLAGS " -Wno-unused-dummy-argument" )
  endif ( CMAKE_C_COMPILER_VERSION VERSION_GREATER 6.0 )
endif ( FTENV_WARNINGS AND LINUX AND NOT USE_INTEL_FORTRAN )

set ( CMAKE_Fortran_MODULE_DIRECTORY .. )
include_directories ( ${CMAKE_CURRENT_BINARY_DIR}/.. )

add_library ( ${LIB_ID}_F JD_dostire.f )
add_library ( ${LIB_ID}_F90 JD_dataModules.f90 JD_namelistModule.f90 )

target_link_libraries ( ${LIB_ID}_F ${LIB_ID}_F90 )
target_link_libraries ( ${LIB_ID}_F90 vpmUtilities_F90 )
