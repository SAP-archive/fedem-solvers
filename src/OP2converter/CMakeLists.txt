# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# Build setup

set ( APPL fedem_op2fmx )
set ( LIB_ID OP2converter )
set ( UNIT_ID ${DOMAIN_ID}_${PACKAGE_ID}_${LIB_ID} )

message ( STATUS "INFORMATION : Processing unit ${UNIT_ID}" )

set ( DEPENDENCY_LIST vpmUtilities_F90 FFaCmdLineArg FFaOS Admin )


# Build and install

file ( GLOB F77_SOURCE_FILES *.f )
add_library ( ${LIB_ID} ${F77_SOURCE_FILES} )
target_link_libraries ( ${LIB_ID} ${DEPENDENCY_LIST} )

add_executable ( ${APPL} convertmain.C )
target_link_libraries ( ${APPL} ${LIB_ID} )
add_dependencies ( all_solvers ${APPL} )
install ( TARGETS ${APPL} RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" )
