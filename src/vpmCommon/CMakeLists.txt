# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# Library setup

set ( LIB_ID vpmCommon )
set ( LIB_ID_LIST vpmLinAlg )
set ( UNIT_ID ${DOMAIN_ID}_${PACKAGE_ID}_${LIB_ID} )

message ( STATUS "INFORMATION : Processing unit ${UNIT_ID}" )

if ( WIN )
  set ( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Qdiag-disable:8290" )
endif ( WIN )

foreach ( FOLDER ${LIB_ID_LIST} )
  add_subdirectory ( ${FOLDER} )
endforeach ( FOLDER ${LIB_ID_LIST} )


## Files with header and source with same name
set ( COMPONENT_FILE_LIST SignalHandler )

## Pure header files, i.e., header files without a corresponding source file
set ( HEADER_FILE_LIST )

## Pure source files, i.e., without a corresponding header file
set ( SOURCE_FILE_LIST cmdLineArgInit cmdLineArgInitStd )

## Fortran90 source files
set ( F90_FILE_LIST environmentTypeModule
                    headingNamelistModule
                    idTypeModule
                    prescribedMotionModule
                    rdbModule
                    samModule
                    supElTypeModule
                    triadTypeModule
                    versionModule
    )


foreach ( FILE ${COMPONENT_FILE_LIST} )
  list ( APPEND CPP_SOURCE_FILES ${FILE}.C )
  list ( APPEND CPP_HEADER_FILES ${FILE}.H )
endforeach ( FILE ${COMPONENT_FILE_LIST} )

foreach ( FILE ${HEADER_FILE_LIST} )
  list ( APPEND CPP_HEADER_FILES ${FILE}.H )
endforeach ( FILE ${HEADER_FILE_LIST} )

foreach ( FILE ${SOURCE_FILE_LIST} )
  list ( APPEND CPP_SOURCE_FILES ${FILE}.C )
endforeach ( FILE ${SOURCE_FILE_LIST} )

foreach ( FILE ${F90_FILE_LIST} )
  list ( APPEND F90_SOURCE_FILES ${FILE}.f90 )
endforeach ( FILE ${F90_FILE_LIST} )


set ( DEPENDENCY_LIST_F vpmUtilities_F90 SAM )
set ( DEPENDENCY_LIST_C vpmUtilities_F90 FFaCmdLineArg Admin )

add_library ( ${LIB_ID}_F90 ${F90_SOURCE_FILES} )
add_library ( ${LIB_ID} ${CPP_SOURCE_FILES} ${CPP_HEADER_FILES} )
target_link_libraries ( ${LIB_ID} ${DEPENDENCY_LIST_C} )
target_link_libraries ( ${LIB_ID}_F90 ${DEPENDENCY_LIST_F} )
