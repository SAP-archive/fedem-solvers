# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# Build setup

set ( LIB_ID_LIST vpmUtilities vpmCommon Femlib )
set ( PACKAGE_ID_LIST OP2converter SolutionMapper FMU
                      vpmReducer vpmSolver vpmStress vpmGraphExp
    )

set ( UNIT_ID ${DOMAIN_ID}_${PACKAGE_ID} )

message ( STATUS "INFORMATION : Processing unit ${UNIT_ID}" )

option ( BUILD_SOLVER_AS_DLL "Build solvers as shared libraries" false )
option ( USE_FFTPACK "Build dynamics solver with FFTPack usage" false )

if ( USE_FFTPACK )
  list ( APPEND LIB_ID_LIST FFTPack )
endif ( USE_FFTPACK )

get_filename_component ( PARENT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY )
include_directories ( "${PROJECT_BINARY_DIR}/f90_modules"
                      "${CMAKE_CURRENT_SOURCE_DIR}"
                      "${PARENT_SOURCE_DIR}/fedem-foundation/src" )

add_custom_target ( all_solvers )

if ( TARGET FFp )
  message ( STATUS "INFORMATION : Including FFp in all_solvers target" )
  add_dependencies ( all_solvers FFp )
endif ( TARGET FFp )

set ( CMAKE_Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/f90_modules" )
foreach ( FOLDER ${LIB_ID_LIST} )
  add_subdirectory ( ${FOLDER} )
endforeach ( FOLDER ${LIB_ID_LIST} )

unset ( CMAKE_Fortran_MODULE_DIRECTORY )
foreach ( FOLDER ${PACKAGE_ID_LIST} )
  add_subdirectory ( ${FOLDER} )
endforeach ( FOLDER ${PACKAGE_ID_LIST} )
