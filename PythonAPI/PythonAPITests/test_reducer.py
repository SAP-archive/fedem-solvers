# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

"""
This Python test basically does the same as the reduction
of one FE part in the folder `../../../solverTests/CarSuspension`.
It therefore reuses the input-files from that folder.

The test uses two environment variables:
    FEDEM_REDUCER = Full path to the reducer shared object library
    TEST_DIR = Full path to the parent folder of the solver tests
"""

from os import environ
from fedempy.reducer import FedemReducer


# Start the fedem reducer in the $TEST_DIR/CarSuspension folder.
# From that folder, read solver options from the knuckle.fco file.
wrkdir = environ["TEST_DIR"] + "/CarSuspension"
slvopt = ["-cwd=" + wrkdir, "-fco=knuckle.fco", "-resfile=knuckle_py.res"]
solver = FedemReducer(environ["FEDEM_REDUCER"], slvopt)

print("\n#### Running FE part reducer in", wrkdir)
ierr = solver.run()
if ierr != 0:
    print(" *** Fedem reducer failed", ierr)
    print("     Check the output file knuckle_py.res")
    exit(ierr)

print("OK!")
