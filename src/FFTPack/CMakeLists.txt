# SPDX-FileCopyrightText: 2024 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# Library setup

set ( LIB_ID fftpack51 )
set ( UNIT_ID ${DOMAIN_ID}_${PACKAGE_ID}_${LIB_ID} )
set ( SRC_DIR fftpack5.1/src )

message ( STATUS "INFORMATION : Processing unit ${UNIT_ID}" )

# Build the library with double precision reals
if ( WIN OR USE_INTEL_FORTRAN )
  string ( APPEND CMAKE_Fortran_FLAGS " -Qautodouble" )
else ( WIN OR USE_INTEL_FORTRAN )
  string ( APPEND CMAKE_Fortran_FLAGS " -fdefault-real-8" )
  if ( CMAKE_C_COMPILER_VERSION VERSION_GREATER 10.0 )
    string ( APPEND CMAKE_Fortran_FLAGS " -fallow-argument-mismatch" )
  endif ( CMAKE_C_COMPILER_VERSION VERSION_GREATER 10.0 )
endif ( WIN OR USE_INTEL_FORTRAN )

# Only including the subroutines that are in use in FEDEN
set ( SRC_FILES ${SRC_DIR}/cfft1i.f ${SRC_DIR}/cfft1f.f ${SRC_DIR}/cfft1b.f
                ${SRC_DIR}/mcfti1.f ${SRC_DIR}/c1fm1f.f ${SRC_DIR}/c1fm1b.f
                ${SRC_DIR}/c1f2kb.f ${SRC_DIR}/c1f2kf.f
                ${SRC_DIR}/c1f3kb.f ${SRC_DIR}/c1f3kf.f
                ${SRC_DIR}/c1f4kb.f ${SRC_DIR}/c1f4kf.f
                ${SRC_DIR}/c1f5kb.f ${SRC_DIR}/c1f5kf.f
                ${SRC_DIR}/c1fgkb.f ${SRC_DIR}/c1fgkf.f
                ${SRC_DIR}/factor.f ${SRC_DIR}/tables.f ${SRC_DIR}/xerfft.f )

add_library ( ${LIB_ID} ${SRC_FILES} )
