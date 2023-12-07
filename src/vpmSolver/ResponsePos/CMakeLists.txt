# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# Build setup

set ( LIB_ID ResponsePos )
set ( UNIT_ID ${DOMAIN_ID}_${PACKAGE_ID}_${LIB_ID} )

message ( STATUS "INFORMATION : Processing unit ${UNIT_ID}" )


## Fortran 90 source files
set ( F90_FILE_LIST  response_pos
                     ../bushingElementTypeModule
                     ../contactElementTypeModule
                     ../contactSurfaceModule
                     ../controlTypeModule
                     ../curvePointTypeModule
                     ../damperTypeModule
                     ../environmentTypeModule
                     ../extCtrlSysTypeModule
                     ../forceTypeModule
                     ../frictionTypeModule
                     ../functionTypeModule
                     ../massTypeModule
                     ../masterSlaveJointTypeModule
                     ../mechanismTypeModule
                     ../modesTypeModule
                     ../motionTypeModule
                     ../normTypeModule
                     ../profilerModule
                     ../roadTypeModule
                     ../saveModule
                     ../sensorTypeModule
                     ../springTypeModule
                     ../supElLoadTypeModule
                     ../systemTypeModule
                     ../tireTypeModule
                     ../userdefElTypeModule
                     ../windTurbineTypeModule
                     ../../vpmStress/saveStrainGageModule
                     ../../vpmStress/strainRosetteModule
                     ../../vpmStress/strainAndStressUtils
                     dummies
    )


foreach ( FILE ${F90_FILE_LIST} )
  list ( APPEND F90_SOURCE_FILES ${FILE}.f90 )
endforeach ( FILE ${F90_FILE_LIST} )


add_library ( ${LIB_ID} ${F90_SOURCE_FILES} )
target_link_libraries ( ${LIB_ID} vpmCommon_F90 FFaFunctionLib_F90 )


set ( APPL response_pos )
add_executable ( ${APPL} main.C )
target_link_libraries ( ${APPL} ${LIB_ID} )
add_dependencies ( all_solvers ${APPL} )
