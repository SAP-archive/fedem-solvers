# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

"""
This Python test is supposed to invoke all interfaces of the solver API, except
for those invoked through other tests, to ensure 100% coverage in a test run.
Consistency of the output data is checked when relevant.

The test uses input files from the solverTests repository,
and the two environment variables:
    FEDEM_SOLVER = Full path to the solver shared object library
    TEST_DIR = Full path to parent folder of the solver tests in the build tree

Here are all interfaces of the class FedemSolver,
see also the file ../src/fedempy/solver.py
(+: Tested in this driver, -: Tested by other drivers):

  + def solver_init(self, options, fsi=None, state_data=None, gauge_data=None):
  - def restart_from_state(self, state_data, write_to_rdb=2):
  - def solve_window(self, n_step, inputs=None, f_out=None, xtimes=None):
  - def get_state_size(self):
  - def get_gauge_size(self):
    def get_transformation_state_size(self):
    def get_part_deformation_state_size(self, base_id):
    def get_part_stress_state_size(self, base_id):
  - def save_state(self):
  - def save_gauges(self):
    def save_transformation_state(self, state_data):
    def save_part_state(self, base_id, def_state, str_state):
  + def solve_next(self, inp=None, inp_def=None, out_def=None, time_next=None):
  + def start_step(self, time_next=None):
  + def solve_iteration(self):
  + def finish_step(self):
  + def solve_modes(self, n_modes, dof_order=False, use_lapack=0):
  + def solve_inverse(self, x_val, x_def, g_def, out_def=None):
  + def solver_done(self, remove_singletons=None):
  - def solver_close(self):
    def run_all(self, options):
  - def set_ext_func(self, func_id, value=None):
  + def get_current_time(self):
  + def get_next_time(self):
    def check_times(self, xtimes, use_times=True):
  + def get_function(self, uid=0, tag=None, arg=None):
  - def get_functions(self, uids):
    def get_function_ids(self, tags):
  + def get_equations(self, bid):
  + def get_system_size(self):
  + def get_system_dofs(self):
  + def get_newton_matrix(self):
  + def get_stiffness_matrix(self):
  + def get_mass_matrix(self):
  + def get_damping_matrix(self):
  + def get_element_stiffness_matrix(self, bid):
  + def get_rhs_vector(self):
  + def get_external_force_vector(self):
  + def set_rhs_vector(self, r_vec):
  + def add_rhs_vector(self, r_vec):
    def compute_strains_from_displ(self, disp, gauge_ids):
    def get_current_strains(self, gauge_ids):
    def compute_rel_dist_from_displ(self, disp, ids):
    def compute_spring_var_from_displ(self, disp, ids):
    def compute_int_forces_from_displ(self, disp, ids):
    def get_joint_spring_stiffness(self, bid):
"""

from os import environ
from fedempy.solver import FedemSolver
from test_utils import compare_lists

################################################################################
# Start the fedem solver in the $TEST_DIR/Cantilever-inverse folder.
# From that folder, read solver options from the Setup.fco file,
# and read model data from the Model.fsi and Gavity.fsi files.
wrkdir = environ["TEST_DIR"] + "/Cantilever-inverse"
solver = FedemSolver(
    environ["FEDEM_SOLVER"],
    ["-cwd=" + wrkdir, "-fco=Setup.fco", "-fsi2file=Gravity.fsi"],  # NOSONAR
)

print("\n#### Running inverse solver in", wrkdir)
print("Comparing with response variables in Master-normalG.asc")

# Read reference values to compare with from file
rfile = open(wrkdir + "/Master-normalG.asc", "r")
for line in rfile:
    print(line.strip())
    if line[0:5] == "#DESC":
        next(rfile)  # skip the first data line (t=0.0)
        break

# List of functions for displacement input and force output
dId = [6, 5]
fId = [1, 2, 3, 4]

# List of equation numbers for the specified displacement DOFs,
# and for the unknown external force DOFs
meqn1 = solver.get_equations(113)
meqn2 = solver.get_equations(17)
meqn3 = solver.get_equations(114)
meqnd = [meqn1[1], meqn2[1]]  # d_y in Triads 113 and 17
meqnf = [meqn3[1], meqn3[5]]  # F_y and M_z in Triad 114
print("   * Specified displacements DOF", meqnd)
print("   * Unknown external force DOF", meqnf)

# Time step loop
IERR = 0
CONT = True
while CONT and solver.ierr.value == 0:
    # Evaluate the displacements (emulated sensor data)
    t = solver.get_next_time()
    dis = [solver.get_function(func_id, None, t) for func_id in dId]
    print("   * Solving inverse problem at t =", t)
    print("     dis =", dis)
    # Solve the inverse problem
    CONT = solver.solve_inverse(dis, meqnd, meqnf)
    # Extract and print the results
    outputs = [float(solver.get_function(func_id)) for func_id in fId]
    print("     out =", outputs)
    outputs.insert(0, float(t))
    # Read reference data from file
    reference = [float(x) for x in next(rfile).split()]
    # Compare response values
    IERR += compare_lists(t, outputs, reference, 1.0e-6)

if IERR > 0:
    print(" *** Detected", IERR, "discrepancies, test not passed")
else:
    print("All response values match")

EK, ok = solver.get_element_stiffness_matrix(19)  # base ID of first beam element
if ok:
    print("   * Element stiffness matrix:\n", EK)
else:
    print(" *** Failed to obtain element stiffness matrix")

# Simulation finished, terminate by closing down the result database, etc.
IERR += abs(solver.solver_done())
if IERR == 0 and solver.ierr.value == 0:
    print("Time step loop OK")
else:
    exit(IERR + abs(solver.ierr.value))

################################################################################
# Start the fedem solver again, now in the $TEST_DIR/DampedOscillatorer folder.
# From that folder, read solver options from the Setup.fco file,
# and read model data from the Model.fsi and linear-damper.fsi files.
# Note that we don't bother with the response values here, since the RHS vector
# is modified arbitrarily. The purpose is only the check that the methods work
# without exceptions.
wrkdir = environ["TEST_DIR"] + "/DampedOscillator"
IERR = solver.solver_init(
    ["-cwd=" + wrkdir, "-fco=Setup.fco", "-fsi2file=linear-damper.fsi"]
)
if IERR < 0:
    exit(IERR)

print("\n#### Running dynamics solver in", wrkdir)
print("Dimension of equation system:", solver.get_system_size())

# Time step loop, up to t=1.0
while solver.get_current_time() < 1.0 and solver.solve_next():
    if solver.ierr.value != 0:
        exit(solver.ierr.value)
print("Simulation paused at", solver.get_current_time())

# Start next time step
CONT = solver.start_step()
if solver.ierr.value != 0:
    exit(solver.ierr.value)

# Extract the system matrices
Kmat, ok = solver.get_stiffness_matrix()
if ok:
    print("Here is the stiffness matrix")
    print(Kmat)
else:
    exit(-99)

Cmat, ok = solver.get_damping_matrix()
if ok:
    print("Here is the damping matrix")
    print(Cmat)
else:
    exit(-99)

Mmat, ok = solver.get_mass_matrix()
if ok:
    print("Here is the mass matrix")
    print(Mmat)
else:
    exit(-99)

Nmat, ok = solver.get_newton_matrix()
if ok:
    print("Here is the Newton matrix")
    print(Nmat)
else:
    exit(-99)

Rvec, ok = solver.get_rhs_vector()
if ok:
    print("Here is the right-hand-side vector")
    print(Rvec)
else:
    exit(-99)

Rvec, ok = solver.get_external_force_vector()
if ok:
    print("Here is the external force vector")
    print(Rvec)
else:
    exit(-99)

# Modify the right-hand side vector
Radd = [0.0] * len(Rvec)
Radd[1] = 5.0
if not solver.add_rhs_vector(Radd):
    exit(-99)
Rvec, ok = solver.get_rhs_vector()
if ok:
    print("Here is the modified right-hand-side vector")
    print(Rvec)
else:
    exit(-99)

Radd[0] = -650
Radd[1] = -2
if not solver.set_rhs_vector(Radd):
    exit(-99)
Rvec, ok = solver.get_rhs_vector()
if ok:
    print("Here is the new right-hand-side vector")
    print(Rvec)
else:
    exit(-99)

# Iteration loop, finishing current time step
while CONT and solver.solve_iteration():
    pass  # NOSONAR

print("Simulation paused at", solver.get_current_time())

# Start next time step
if CONT:
    CONT = solver.start_step()

# Extract the system matrices
Nmat, ok = solver.get_newton_matrix()
if ok:
    print("Here is the Newton matrix")
    print(Nmat)
else:
    exit(-99)

# Finishing current time step (iterations)
if CONT:
    CONT = solver.finish_step()

print("Simulation continuing from", solver.get_current_time())

# Continue the time step loop until the end
while CONT and solver.ierr.value == 0:
    CONT = solver.solve_next()

# Simulation finished, terminate by closing down the result database, etc.
IERR = solver.solver_done()
if IERR == 0 and solver.ierr.value == 0:
    print("Time step loop OK")
else:
    exit(IERR + solver.ierr.value)

################################################################################
# Start the fedem solver again, now in the $TEST_DIR/BeamModes folder.
# From that folder, read solver options from the Setup.fco file,
# and read model data from the Model.fsi and Water.fsi files.
wrkdir = environ["TEST_DIR"] + "/BeamModes"
IERR = solver.solver_init(
    [
        "-cwd=" + wrkdir,
        "-fco=Setup.fco",
        "-eigInc=2.0",
        "-fsi2file=Water.fsi",
        "-timeInc=0.1",
        "-timeEnd=2.0",
    ]
)
if IERR < 0:
    exit(IERR)

print("\n#### Running dynamics solver in", wrkdir)
print("Dimension of equation system:", solver.get_system_size())
print("Number of degrees of freedom:", solver.get_system_dofs())


def get_indices(vec, target):
    """
    Get indices in vec, maching target value
    """
    return [i for i, x in enumerate(vec) if abs(x - target) < 1.0e-8 * abs(target)]


def print_modes(eig_vals, eig_vecs):
    """
    Utility function to print eigenmodes
    """
    if eig_vals is None:
        return
    print("Eigenvalues:", eig_vals)

    if eig_vecs is None:
        return
    print("Eigenvector dimension:", len(eig_vecs[0]))

    ivec = 0
    for evec in eig_vecs:
        ivec = ivec + 1
        vmax = max(evec)
        vmin = min(evec)
        print(
            "Eigenvector",
            ivec,
            " min =",
            vmin,
            get_indices(evec, vmin),
            " max =",
            vmax,
            get_indices(evec, vmax),
        )


# Do eigenvalue analysis at initial configuration
# Extract eigenvectors in DOF-order
e_val, e_vec, ok = solver.solve_modes(9, True)
print_modes(e_val, e_vec)

# Compare with reference values (from the file wet_modes.asc)
output = [float(e_val[i]) for i in (8, 6, 4, 2, 0)]
reference = [5.62979221, 4.47576094, 3.34034538, 2.2190156, 1.10713613]
IERR = compare_lists(0.0, output, reference, 0.001)
if IERR > 0:
    print(" *** Detected", IERR, "discrepancies, test not passed")
else:
    print("All eigenvalues match")

# Time step loop, up to t=1.0
while solver.get_current_time() < 1.0 and solver.solve_next():
    if solver.ierr.value != 0:
        exit(solver.ierr.value)
print("Simulation paused at", solver.get_current_time())

# Do eigenvalue analysis at t=1.0
# Extract eigenvectors in equation order
e_val, e_vec, ok = solver.solve_modes(12)
print_modes(e_val, e_vec)

# Simulation finished, terminate by closing down the result database, etc.
IERR += abs(solver.solver_done())
if IERR == 0 and solver.ierr.value == 0:
    print("Time step loop OK")
else:
    exit(IERR + abs(solver.ierr.value))
