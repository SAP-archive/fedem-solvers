<!---
  SPDX-FileCopyrightText: 2023 SAP SE

  SPDX-License-Identifier: Apache-2.0

  This file is part of FEDEM - https://openfedem.org
--->

# FEDEM solvers Changelog

## [fedem-8.0.1] (2024-01-25)

### :rocket: Added

- The frequency domain solution methods are reinstated,
  using the FFTPACK 5.1 library which now is embedded in this repository.

### :bug: Fixed

- (Cosmetic) A newline follows the object description in the frs-file header.
- Direct curve export from the dynamics solver does not work,
  due do the issue [SAP/fedem-gui#5](https://github.com/SAP/fedem-gui/issues/5).

## fedem-8.0.0 (2023-12-21)

### :rocket: Added

- First open source version of the FEDEM solvers, including a Python API

[fedem-8.0.1]: https://github.com/SAP/fedem-gui/compare/fedem-8.0.0...fedem-8.0.1
