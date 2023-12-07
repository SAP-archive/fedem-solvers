# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

"""
This Python test basically does the same as the solver test
located in the folder `../../../solverTests/Cantilever-extFunc`.
It therefore reuses the input-files from that folder.

The test uses two environment variables:
    FEDEM_SOLVER = Full path to the solver shared object library
    TEST_DIR = Full path to the parent folder of the solver tests
"""

from os import environ
from math import sin, cos
from fedempy.solver import FedemSolver
from test_utils import compare_lists


def ext_func(func_id, x):
    """
    Hard-coded external function evaluations to test against.
    """
    if func_id == 1:
        return 1.0e7 * sin(x)
    elif func_id == 2:
        return 5.0e3 * (cos(2.5 * x) - 1.0)
    else:
        return 0.0


# List of function IDs to extract results for.
fId = [3, 4, 5]
n_out = len(fId)

# Start the fedem solver in the $TEST_DIR/Cantilever-extFunc folder.
# From that folder, read solver options from the Setup.fco file,
# and read model data from the Model.fsi file.
wrkdir = environ["TEST_DIR"] + "/Cantilever-extFunc"
slvopt = ["-cwd=" + wrkdir, "-fco=Setup.fco", "-fsifile=Model.fsi"]
solver = FedemSolver(environ["FEDEM_SOLVER"], slvopt, True)

print("\n#### Running dynamics solver in", wrkdir)
print("Comparing with response variables in TipPosition.asc")

# Read reference values to compare with from file
rfile = open(wrkdir + "/TipPosition.asc", "r")
for line in rfile:
    print(line.strip())
    if line[0:5] == "#DESC":
        next(rfile)  # skip the first data line (t=0.0)
        break

# Time step loop
ierr = 0
do_continue = True
references = []
while do_continue and solver.ierr.value == 0:
    # Evaluate the external functions
    t = solver.get_next_time()
    solver.set_ext_func(1, ext_func(1, t))
    solver.set_ext_func(2, ext_func(2, t))
    # Solve for next time step
    do_continue = solver.solve_next()
    # Extract the results
    outputs = [float(solver.get_function(func_id)) for func_id in fId]
    outputs.insert(0, float(t))
    # Read reference data from file
    reference = [float(x) for x in next(rfile).split()]
    # Compare response values
    ierr += compare_lists(t, outputs, reference, 1.0e-8)
    references.append(reference)
if ierr > 0:
    print(f" *** Detected {ierr} discrepancies, test not passed")
else:
    print(f"   * All response values match, checked {len(references)} steps")

# Simulation finished, terminate by closing down the result database, etc.
rfile.close()
ierr += abs(solver.solver_done())
if ierr == 0 and solver.ierr.value == 0:
    print("Time step loop OK, solver closed")
else:
    exit(ierr + abs(solver.ierr.value))

# Start over, read all input once again (from same working directory)
ierr = solver.solver_init(slvopt)
if ierr < 0:
    exit(ierr)

# Run 100 steps using the Window function
n_step = 100
inputs = [0.0] * (2 * n_step)
t = solver.get_current_time()
dt = solver.get_next_time() - t
for i in range(n_step):
    inputs[2 * i] = ext_func(1, t + dt * (i + 1))
    inputs[2 * i + 1] = ext_func(2, t + dt * (i + 1))
outputs, stat = solver.solve_window(n_step, inputs, fId)
if not stat:
    print("End of simulation was reached")

# Compare response values
ierr = 0
for i in range(n_step):
    t = references[i][0]
    output = [float(outputs[i * n_out + j]) for j in range(n_out)]
    output.insert(0, t)
    ierr += compare_lists(t, output, references[i], 1.0e-8)
if ierr > 0:
    print(f" *** Detected {ierr} discrepancies, test not passed")
else:
    print(f"   * All response values match, checked {n_step} steps")

# Simulation finished, terminate by closing down the result database, etc.
ierr += abs(solver.solver_done())
if ierr == 0 and solver.ierr.value == 0:
    print("Time window OK, solver closed")
else:
    exit(ierr + abs(solver.ierr.value))
