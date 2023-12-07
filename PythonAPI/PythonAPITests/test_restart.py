# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

"""
This Python test basically does the same as one of the solver tests
located in the folder `../../../solverTests/AW-N80-80-Windturbine`.
It therefore reuses the input-files from that folder.

The test uses two environment variables:
    FEDEM_SOLVER = Full path to the solver shared object library
    TEST_DIR = Full path to parent folder of the solver tests in the build tree
"""

from os import environ
from fedempy.solver import FedemSolver
from test_utils import compare_lists, read_reference_data

# List of function IDs to extract results for
fId = [
    155,
    156,
    153,
    154,
    135,
    136,
    139,
    140,
    143,
    144,
    147,
    148,
    149,
    150,
    145,
    146,
    141,
    142,
    137,
    138,
]
n_out = len(fId)

ierr = 0
n_step = 100  # Time window size

# Start the fedem solver in the $TEST_DIR/AW-N80-80-Windturbine folder.
wrkdir = environ["TEST_DIR"] + "/AW-N80-80-Windturbine"
# From that folder, read solver options from the Setup.fco file,
# and read model data from the Model.fsi and Gages.fsi files.
# Run with -recovery=2 to perform gage recovery during time integration.
slv_opts = [
    "-cwd=" + wrkdir,
    "-fco=Setup.fco",
    "-fsifile=Model.fsi",
    "-fsi2file=Gages.fsi",
    "-recovery=2",
]
solver = FedemSolver(environ["FEDEM_SOLVER"], slv_opts, True)

print("\n#### Running dynamics solver in", wrkdir)
print("Dimension of equation system:", solver.get_system_size())
print("State array size:", solver.state_size)
print("Gages array size:", solver.gauge_size)
# Read reference values to compare with from file
references = read_reference_data(wrkdir + "/Responses.asc")
print("Comparing with response variables in Responses.asc", len(references))

# Outer time window loop
step = 0
while ierr == 0:
    # Caution: We here assume the time step size in the model is constant
    t = solver.get_current_time()
    dt = solver.get_next_time() - t

    # Invoke the dynamics solver for this time window
    outputs, do_continue = solver.solve_window(n_step, None, fId)
    if not (do_continue and solver.ierr.value == 0):
        break

    # Save current state to core array
    if not solver.save_state():
        exit(-99)
    if not solver.save_gauges():
        exit(-99)

    # Compare the responses
    for i in range(n_step):
        output = [float(outputs[i * n_out + j]) for j in range(n_out)]
        output.insert(0, float(t + dt * (i + 1)))
        ierr += compare_lists(output[0], output, references[step], 1.0e-4)
        step += 1

    if ierr > 0:
        print(" *** Detected", ierr, "discrepancies, test not passed")
    else:
        print("All response values match, checked", n_step, "steps")

    # Time window finished, terminate by closing down the result database, etc.
    ierr += abs(solver.solver_done())
    if ierr == 0:
        print("Time window OK, solver closed")
    else:
        exit(ierr)

    # Restart the simulation, reading the state and gages from core arrays
    ierr = solver.solver_init(slv_opts, None, solver.state_data, solver.gauge_data)

# Simulation finished, terminate by closing down the result database, etc.
ierr += solver.solver_done()
if ierr == 0 and solver.ierr.value == 0:
    print("Simulation OK, solver closed")
else:
    exit(ierr + solver.ierr.value)

# Start the simulation again from the beginning
ierr = solver.solver_init(slv_opts)
if ierr < 0:
    exit(ierr)

print("\n#### Running dynamics solver in", wrkdir)
t = solver.get_current_time()
dt = solver.get_next_time() - t

# Invoke the dynamics solver for this time window
outputs, do_continue = solver.solve_window(n_step, None, fId)
if solver.ierr.value != 0:
    exit(solver.ierr.value)

# Save current state to core array
if not solver.save_state():
    exit(-99)

# Compare the responses
ierr = 0
for i in range(n_step):
    output = [float(outputs[i * n_out + j]) for j in range(n_out)]
    output.insert(0, float(t + dt * (i + 1)))
    ierr += compare_lists(output[0], output, references[i], 1.0e-4)

if ierr > 0:
    print(" *** Detected", ierr, "discrepancies, test not passed")
    exit(ierr)
else:
    print("All response values of first time window match")

# Restart the simulation, reading the state from from core array
# Notice, here we don't do a solverDone first, so the model reside in core.
# Only the solution state is reloaded.
print("\n#### Continuing with second time window")
ierr = solver.restart_from_state(solver.state_data)
if ierr < 0:
    exit(ierr)

# Invoke the dynamics solver for next time window
t = solver.get_current_time()
dt = solver.get_next_time() - t
outputs, do_continue = solver.solve_window(n_step, None, fId)
if solver.ierr.value != 0:
    exit(solver.ierr.value)

# Compare the responses
ierr = 0
for i in range(n_step):
    output = [float(outputs[i * n_out + j]) for j in range(n_out)]
    output.insert(0, float(t + dt * (i + 1)))
    ierr += compare_lists(output[0], output, references[n_step + i], 1.0e-4)

if ierr > 0:
    print(" *** Detected", ierr, "discrepancies, test not passed")
else:
    print("All response values of second time window match")

# Simulation finished, terminate by closing down the result database, etc.
ierr += abs(solver.solver_done())
if ierr == 0:
    print("Simulation OK, solver closed")
else:
    exit(ierr)
