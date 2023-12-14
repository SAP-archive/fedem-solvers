<!---
  SPDX-FileCopyrightText: 2023 SAP SE

  SPDX-License-Identifier: Apache-2.0

  This file is part of FEDEM - https://openfedem.org
--->

[![REUSE status](https://api.reuse.software/badge/github.com/SAP/fedem-solvers)](https://api.reuse.software/info/github.com/SAP/fedem-solvers)

# FEDEM Solvers

![Fedem Logo](https://github.com/SAP/fedem-foundation/blob/main/cfg/FedemLogo.png "Welcome to FEDEM")

**Welcome to FEDEM! - Finite Element Dynamics of Elastic Mechanisms.**

## About this project

This project contains the Fortran source code of the FEDEM solver modules.
In addition, the [PythonAPI](PythonAPI) folder contains the source code of
the Python package `fedempy`, which enables the use of python to control
execution of the FE Part Reducer and Dynamics Solver,
as well as creating new models and/or editing existing ones by scripting.

This project also uses elements from the
[fedem-foundation](https://github.com/SAP/fedem-foundation) project,
which is consumed as a submodule by this repository.
It also uses the [SAM library](https://github.com/SAP/sam-lib).

The FEDEM GUI application, which also can be used to control the solvers
and to post-process their results, is maintained in a parallel project
[fedem-gui](https://github.com/SAP/fedem-gui).

Refer to our web page [openfedem.org](https://openfedem.org/)
for overall information on the FEDEM project.

## Requirements and Setup

See [BUILD.md](BUILD.md) for detailed information on building and testing.
Notice the regression tests for the Dynamics Solver are maintained in a separate
repository [fedem-solver-tests](https://github.com/SAP/fedem-solver-tests),
which also is consumed as a submodule by this repository. Therefore,
make sure you get all submodules resolved when cloning this repository, e.g.:

    $ git clone --recurse-submodules git@github.com:SAP/fedem-solvers.git

## Contributing

This project is open to feature requests/suggestions, bug reports etc. via [GitHub issues](https://github.com/SAP/fedem-solvers/issues). Contribution and feedback are encouraged and always welcome. For more information about how to contribute, the project structure, as well as additional contribution information, see our [Contribution Guidelines](CONTRIBUTING.md).

## Security / Disclosure

If you find any bug that may be a security problem, please follow our instructions at [in our security policy](https://github.com/SAP/fedem-solvers/security/policy) on how to report it. Please do not create GitHub issues for security-related doubts or problems.

## Code of Conduct

We as members, contributors, and leaders pledge to make participation in our community a harassment-free experience for everyone. By participating in this project, you agree to abide by its [Code of Conduct](https://github.com/SAP/.github/blob/main/CODE_OF_CONDUCT.md) at all times.

## Licensing

Copyright 2023 SAP SE or an SAP affiliate company and fedem-solvers contributors. Please see our [LICENSE](LICENSE) for copyright and license information. Detailed information including third-party components and their licensing/copyright information is available [via the REUSE tool](https://api.reuse.software/info/github.com/SAP/fedem-solvers).
