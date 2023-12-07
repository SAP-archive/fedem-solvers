<!---
  SPDX-FileCopyrightText: 2023 SAP SE

  SPDX-License-Identifier: Apache-2.0

  This file is part of FEDEM - https://openfedem.org
--->

# FEDEM Solvers Build

Here is the details you need to perform a local build and test
of the FEDEM solvers

## Build instruction

The build system uses the cross-platform [CMake](https://cmake.org/) tool
to perform the build on Linux and Windows platforms.
The top level `CMakeLists.txt` file of this repository therefore contains
the necessary setup and configuration. It also controls the execution of
the unit- and regression test [see below](#regression-testing).

We use the packages [googletest](https://github.com/google/googletest) and
[pFUnit](https://github.com/Goddard-Fortran-Ecosystem/pFUnit) to implement
some unit tests for the C++ and Fortran code, respectively. Therefore,
you need to set the environment variables `GTEST_ROOT` and `PFUNIT` to point to
the root path of valid installations of the googletest and pFUnit frameworks,
respectively, before executing the cmake command shown below.

On Linux, the FEDEM solvers and tests can be built and executed by:

    mkdir Release; cd Release
    cmake .. -DCMAKE_BUILD_TYPE=Release
    make check
