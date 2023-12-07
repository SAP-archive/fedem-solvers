<!---
  SPDX-FileCopyrightText: 2023 SAP SE

  SPDX-License-Identifier: Apache-2.0

  This file is part of FEDEM - https://openfedem.org
--->

# The FEDEM python API

This folder contains the source code of the python wrapper of the FEDEM solver API.
The source code resides in the sub-folder [src/fedempy](src/fedempy) which is
compiled into the python package `fedempy` by the `setup.py` script.
In addition, some test drivers are placed in the sub-folder [PythonAPITests](PythonAPITests)
which are invoked as regression/integration tests by the cmake-based build system.

See [here](https://github.wdf.sap.corp/pages/FEDEM/digital-twin-fedempy/index.html),
for the extracted source code documentation for the python code, which is generated
using the [Sphinx](https://www.sphinx-doc.org) tool. That page also contains some
basic installation and end-user documentation.

## Local build of python package

Run the command

    python setup.py sdist

This will create the archive `fedempy-<VERSION>.tar.gz` where `<VERSION>` is the content of the
file `version.txt`. To install the built package, use the command:

    pip install dist/fedempy-<VERSION>.tar.gz

## Static code analysis

Static code analysis is done for each pull request and trending is done over time on the
master branch using SonarQube.

## Local checks on the code

Run `tox -e type` to type check the code with mypy.

Run `tox -e lint` to check coding standard and syntax.

Run `tox -e docs` to extract source code documentation in html-format.
The generated files are found in the sub-folder `.tox/docs/tmp/html`.
