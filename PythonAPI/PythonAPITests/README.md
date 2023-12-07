<!---
  SPDX-FileCopyrightText: 2023 SAP SE

  SPDX-License-Identifier: Apache-2.0

  This file is part of FEDEM - https://openfedem.org
--->

# Python API tests

This folder contains a set of regression tests for the Python API of Fedem.
The test execution is governed by `ctest` through the `CMakeLists.txt` file.

Some of the tests requires the dynamics solver binary only, and are therefore
invoked by the build process of the ftKernel repository.

The main bulk of the tests relies on accessing (or creating) Fedem model files
using the fedem_mdb module. They are therefore only invoked from the build
process of that module.
