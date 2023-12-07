# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# Build setup

set ( APPLICATION_SOURCE_FILE_LIST stress modes gage fpp )
set ( LIB_ID vpmStress )
set ( DLL_ID fedem_recovery_core )
set ( LIB_ID_LIST vpmStressTests )
set ( UNIT_ID ${DOMAIN_ID}_${PACKAGE_ID}_${LIB_ID} )

message ( STATUS "INFORMATION : Processing unit ${UNIT_ID}" )

foreach ( FOLDER ${LIB_ID_LIST} )
  add_subdirectory ( ${FOLDER} )
endforeach ( FOLDER ${LIB_ID_LIST} )


## Files with header and source with same name
set ( COMPONENT_FILE_LIST )
if ( BUILD_SOLVER_AS_DLL )
  set ( COMPONENT_FILE_LIST ../vpmCommon/FFCmdLineArgInterface )
endif ( BUILD_SOLVER_AS_DLL )

## Pure header files, i.e., header files without a corresponding source file
set ( HEADER_FILE_LIST )

## Pure source files, i.e., without a corresponding header file
set ( SOURCE_FILE_LIST stressInterface )

## Fortran 90 source files
set ( F90_FILE_LIST displacementModule
                    initiateTriadAndSupElTypeModule
                    samStressModule strainAndStressUtils
                    fpp fatigueModule strainCoatModule
                    gage strainRosetteModule strainGageModule
                    modes modesRoutines
                    stress stressRoutines elStressModule
                    saveNastranModule saveVTFModule
                    saveStrainCoatModule
                    saveStrainGageModule
                    saveStressModule
                    resStressModule
    )
if ( USE_FPPINTERFACE )
  string ( APPEND CMAKE_Fortran_FLAGS " -DFT_HAS_FPPINTERFACE" )
  list ( APPEND F90_FILE_LIST fppModule )
endif ( USE_FPPINTERFACE )

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


set ( FT_COMMON_LIBRARIES FemLib vpmUtilities_F90
                          FFlLib_F FFrLib_F FFaLib_F
    )
set ( FT_KERNEL_LIBRARIES vpmCommon_F90
                          FiDeviceFunctions_F90
                          FFpLib_F
    )

set ( DEPENDENCY_LIST ${FT_KERNEL_LIBRARIES} ${FT_COMMON_LIBRARIES} )


# Build and install

add_library ( ${LIB_ID}_F ${F90_SOURCE_FILES} )
target_link_libraries ( ${LIB_ID}_F ${DEPENDENCY_LIST} )
if ( BUILD_SOLVER_AS_DLL )
  string ( APPEND CMAKE_CXX_FLAGS " -DFT_USE_DLL" )
  add_library ( ${DLL_ID} SHARED ${CPP_SOURCE_FILES} ${CPP_HEADER_FILES} )
  target_link_libraries ( ${DLL_ID} LINK_PRIVATE ${LIB_ID}_F vpmCommon )
else ( BUILD_SOLVER_AS_DLL )
  add_library ( ${LIB_ID} STATIC ${CPP_SOURCE_FILES} ${CPP_HEADER_FILES} )
  target_link_libraries ( ${LIB_ID} LINK_PUBLIC ${LIB_ID}_F vpmCommon )
endif ( BUILD_SOLVER_AS_DLL )

set ( TARGET_LIST )
foreach ( APPL ${APPLICATION_SOURCE_FILE_LIST} )
  add_executable ( fedem_${APPL} ${APPL}main.C )
  if ( BUILD_SOLVER_AS_DLL )
    target_link_libraries ( fedem_${APPL} ${DLL_ID} )
  else ( BUILD_SOLVER_AS_DLL )
    target_link_libraries ( fedem_${APPL} ${LIB_ID} )
  endif ( BUILD_SOLVER_AS_DLL )
  list ( APPEND TARGET_LIST fedem_${APPL} )
  add_dependencies ( all_solvers fedem_${APPL} )
endforeach ( APPL ${APPLICATION_SOURCE_FILE_LIST} )
if ( BUILD_SOLVER_AS_DLL )
  list ( APPEND TARGET_LIST ${DLL_ID} )
endif ( BUILD_SOLVER_AS_DLL )
do_install ( "${TARGET_LIST}" ${CMAKE_INSTALL_PREFIX} )
