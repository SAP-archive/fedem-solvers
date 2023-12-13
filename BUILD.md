<!---
  SPDX-FileCopyrightText: 2023 SAP SE

  SPDX-License-Identifier: Apache-2.0

  This file is part of FEDEM - https://openfedem.org
--->

# FEDEM Solvers Build

Currently, we support building on Windows
using Microsoft Visual Studio 2019 (or later) and Intel&reg; Fortran Compilers,
and on Linux using the GNU compilers (gcc and gfortran).
The build system is based on the cross-platform [CMake](https://cmake.org/) tool,
which is embedded in the Visual Studio 2019 installation.

## Build instruction

The top level `CMakeLists.txt` file of this repository contains
the necessary setup and configuration. It also controls the execution of
the unit- and regression test [see below](#regression-testing).

We use the packages [googletest](https://github.com/google/googletest) and
[pFUnit](https://github.com/Goddard-Fortran-Ecosystem/pFUnit) to implement
some unit tests for the C++ and Fortran code, respectively. Therefore,
you need to set the environment variables `GTEST_ROOT` and `PFUNIT` to point to
the root path of valid installations of the googletest and pFUnit frameworks,
respectively, before executing the `cmake` commands shown below.

- From a bash shell or command prompt, clone the sources of this repository:

      mkdir ~/Fedem-src
      cd ~/Fedem-src
      git clone --recurse-submodules git@github.com:SAP/fedem-solvers.git

- On Linux, the FEDEM solvers and tests can be built and executed by:

      cd ~/Fedem-src/fedem-solvers
      mkdir Release; cd Release
      cmake .. -DCMAKE_BUILD_TYPE=Release \
               -DBUILD_SOLVER_AS_DLL=ON -DBUILD_CONTROL_AS_DLL=ON -DBUILD_TEST_REPORTS=ON \
               -DUSE_CONCURRENT_RECOVERY=ON -DUSE_SP_RECOVERY=ON -DFT_LARGE_MODELS=ON -DFT_TOLERANCE=1.0e-12
      make check

- On Windows, use the following bat script to configure the build with
  Visual Studio 2019 and the Intel&reg; Fortran Compilers.

      @echo off
      title Open FEDEM solvers configuration for Windows
      call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 vs2019
      "%VS2019INSTALLDIR%\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe" ^
      -G "Visual Studio 16 2019" ^
      -DBUILD_SOLVER_AS_DLL=ON -DBUILD_CONTROL_AS_DLL=ON ^
      -DUSE_CONCURRENT_RECOVERY=ON -DUSE_SP_RECOVERY=ON ^
      -DFT_LARGE_MODELS=ON -DFT_TOLERANCE=1.0e-10 ^
      -S %USERPROFILE%\Fedem-src\fedem-solvers ^
      -B %USERPROFILE%\Fedem-build\solvers
      pause

- Then (still on Windows), open the generated solution file
  `%USERPROFILE%\Fedem-build\solvers\fedemKernel.sln`
  in Visual Studio and build the `all_solvers` target for `Release`
  configuration to compile all solver modules. Build the `check` target
  if you also want to execute the tests after building the solvers.
