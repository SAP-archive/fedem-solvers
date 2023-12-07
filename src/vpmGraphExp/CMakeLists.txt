# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# Build setup

set ( APPL fedem_graphexp )
set ( LIB_ID vpmGraphExpr )
set ( UNIT_ID ${DOMAIN_ID}_${PACKAGE_ID}_${LIB_ID} )

message ( STATUS "INFORMATION : Processing unit ${UNIT_ID}" )

set ( DEPENDENCY_LIST FFpExport FFaCmdLineArg FFaString Admin )


# Build and install

add_executable ( ${APPL} exportmain.C ../vpmCommon/cmdLineArgInitStd.C )
target_link_libraries ( ${APPL} ${DEPENDENCY_LIST} )
add_dependencies ( all_solvers ${APPL} )
do_install ( ${APPL} ${CMAKE_INSTALL_PREFIX} )
