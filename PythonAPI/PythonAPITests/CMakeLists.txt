# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

set ( UNIT_ID "${DOMAIN_ID}_${PACKAGE_ID}_PythonAPITests" )

message ( STATUS "INFORMATION : Processing unit ${UNIT_ID}" )

################################################################################
# Set up regression tests for the Python wrapping of the FEDEM solver libraries.
# Each test uses a python driver (test_<name>.py) which is invoked by specifying
# the ${PYTHON_EXECUTABLE} as the command in the add_test() calls, and with a
# generated wrapper that invokes the actual test driver scripts as argument.

set ( TEST_DIR "${SRC_DIR}/vpmSolver/vpmSolverTests/solverTests" )
if ( PYTHON_DIR )
  set ( PYTHONPATH "${PYTHON_DIR}${CMAKE_PATH_SEP}${CMAKE_CURRENT_SOURCE_DIR}" )
else ( PYTHON_DIR )
  set ( PYTHONPATH "${CMAKE_CURRENT_SOURCE_DIR}" )
endif ( PYTHON_DIR )
if ( NOT CMAKE_BUILD_TYPE )
  # Build type not defined (Windows), assume Release build here
  message ( STATUS "NOTE : The Python API tests will always use Release build" )
  set ( CMAKE_BUILD_TYPE Release )
endif ( NOT CMAKE_BUILD_TYPE )

if ( TARGET fedem_solver_core )

  ##############################################################################
  # Generate a Python wrapper that loads the fedem solver shared object library
  # for current build configuration, while executing some solver/api tests.
  file ( GENERATE OUTPUT
         ${CMAKE_CURRENT_BINARY_DIR}/testPyAPI_$<CONFIGURATION>.py
         CONTENT "from os import environ
environ['FEDEM_SOLVER']='$<TARGET_FILE:fedem_solver_core>'
from sys import argv
do_restart = True if len(argv) > 1 else False
import test_api
import test_extfunc
if len(argv) > 1 and argv[1] == '--with-restart':
    import test_restart
import test_Microbatch
" )
  message ( STATUS "INFORMATION : Adding regression test PythonAPI" )
  if ( NOT TARGET AW-N80-80 AND FT_LARGE_MODELS )
    # Include the restart test only if a corresponding solver test is activated
    # (see .../solverTests/AW-N80-80/CMakeLists.txt), since this test will
    # re-use the reduced superelement matrices of that test
    set ( RESTART "--with-restart" )
  endif ( NOT TARGET AW-N80-80 AND FT_LARGE_MODELS )
  add_test ( PythonAPI ${PYTHON_EXECUTABLE} -B testPyAPI_${CMAKE_BUILD_TYPE}.py ${RESTART} )
  set_tests_properties ( PythonAPI PROPERTIES ENVIRONMENT
                         "PYTHONPATH=${PYTHONPATH};TEST_DIR=${TEST_DIR}" )

endif ( TARGET fedem_solver_core )
if ( TARGET fedem_reducer_core )

  ##############################################################################
  # Generate a Python wrapper that loads the fedem reducer shared object library
  # for current build configuration, while executing some reduction tests.
  file ( GENERATE OUTPUT
         ${CMAKE_CURRENT_BINARY_DIR}/testReducer_$<CONFIGURATION>.py
         CONTENT "from os import environ
environ['FEDEM_REDUCER']='$<TARGET_FILE:fedem_reducer_core>'
import test_reducer
" )
  message ( STATUS "INFORMATION : Adding regression test PythonReducer" )
  add_test ( PythonReducer ${PYTHON_EXECUTABLE} -B testReducer_${CMAKE_BUILD_TYPE}.py )
  set_tests_properties ( PythonReducer PROPERTIES ENVIRONMENT
                         "PYTHONPATH=${PYTHONPATH};TEST_DIR=${TEST_DIR}" )

endif ( TARGET fedem_reducer_core )
if ( TARGET fmxWriter )
  add_dependencies ( check fmxWriter )

  ##############################################################################
  # Generate a Python wrapper that loads the correct shared object library
  # for current build configuration, while executing some fmx-writing tests.
  file ( GENERATE OUTPUT
         ${CMAKE_CURRENT_BINARY_DIR}/testFMX_$<CONFIGURATION>.py
         CONTENT "from os import environ
environ['FMXWRITER']='$<TARGET_FILE:fmxWriter>'
import test_fmx
" )
  message ( STATUS "INFORMATION : Adding regression test FMXWrite" )
  add_test ( FMXwrite ${PYTHON_EXECUTABLE} -B testFMX_${CMAKE_BUILD_TYPE}.py )
  set_tests_properties ( FMXwrite PROPERTIES ENVIRONMENT "PYTHONPATH=${PYTHONPATH}" )

endif ( TARGET fmxWriter )
