# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# Build setup

set ( LIB_ID vpmUtilitiesTest )
set ( UNIT_ID ${DOMAIN_ID}_${PACKAGE_ID}_${LIB_ID} )

message ( STATUS "INFORMATION : Processing unit ${UNIT_ID}" )

# Old-style unit tests (executed manually from console)
set ( OLD_TESTS writeBinary putChar dmOpen )
if ( NOT pFUnit_FOUND )
  list ( APPEND OLD_TESTS rotation )
endif ( NOT pFUnit_FOUND )

foreach ( TEST ${OLD_TESTS} )
  add_executable ( test_${TEST} ${TEST}Test.f90 usage.f90 )
  target_link_libraries ( test_${TEST} vpmUtilities_F90 )
endforeach ( TEST ${OLD_TESTS} )

if ( pFUnit_FOUND )
# Add some unit tests using pFUnit (executed via ctest)
  enable_fortran_tests ()
  add_fortran_test ( test_rotation testRotation.pf vpmUtilities_F90 )
  add_fortran_test ( test_Bmatrix testBmatrix.pf vpmUtilities_F90 )
endif ( pFUnit_FOUND )
