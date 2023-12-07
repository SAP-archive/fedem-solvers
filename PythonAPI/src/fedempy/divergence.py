# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

"""
Convenience module for automatic simulation restart of diverging models.
"""
from ctypes import c_double
from json import dumps, load
from os import getcwd, mkdir, path


def ramp_dataframe(dfr, start_pos=0):
    """
    Creating a s-shape ramp, starting with 0 and ending with 1
    Second derivative at start/end will be zero
    Multiplying each column in dataframe with ramp-array
    """
    r_steps = len(dfr.index.tolist())  # ramp steps

    smooth = [1.0] * r_steps
    scale_step = 1.0 / (r_steps - start_pos - 1)
    for i in range(start_pos, r_steps):
        s_x = scale_step * (i - start_pos)
        smooth[i] = s_x * s_x * s_x * (s_x * (s_x * 6.0 - 15.0) + 10.0)

    return dfr.multiply(smooth, axis=0)


def jump_state(lib_dir):
    """
    Checks whether the current state should be jumped over or not.
    """
    jump_exist = False
    json_file = lib_dir + "/stream.json"
    if path.exists(json_file):
        with open(json_file, "r") as fd:
            frame_data = load(fd)
            jump_exist = frame_data.get("conv_type") == "JUMP_OVER"

    if jump_exist:
        # Reset conv_type to None
        with open(json_file, "w") as fd:
            fd.write(dumps({"conv_type": None}, indent=2))

    return jump_exist


def _opts_list(opts_dict):
    """
    Helper converting option dictionary into a list of solver options.
    """

    def __format_arg(key, value):
        if value is None:
            return None
        if isinstance(value, bool):
            return f"-{key}+" if value else f"-{key}-"
        return f"-{key}={value}"

    return [__format_arg(k, v) for k, v in opts_dict.items()]


def restart(solver, df, state, lib_dir, out_id, c_opts, x_times=None):  # NOSONAR
    """
    Restarts the simulation over a time window in case of divergence issues.
    When a divergent solution is detected, one (or more) of the following
    approaches may be attempted:

    * `TOL`:
      Increasing tolerance(s) for convergence checks
      (tolDispNorm, tolVecNorm, tolEnerSum, maxit).
    * `STATE_DYNAMIC`:
      Use dynamic equilibrium at restart.
      Restart can be performed at a specified number of steps before the divergent step,
      optionally with ramp-scaling activated (default activated).
    * `STATE_STATIC`:
      Use quasi-static equilibrium for a specified number of steps at restart.
      Restart can be performed at a specified number of steps before the divergent step,
      optionally with ramp-scaling activated (default activated) done n-steps backwards.
    * `NO_STATE_WITH_RAMP`:
      Restart of specified number of steps before the diverged step.
      The loads are ramped up from zero and the state is not used.
    * `JUMP_OVER`:
      Skip the rest if the current window and restart from the next,
      with the loads ramped up from zero.
    * `NO_CONV`:
      No convergence is obtained.
    * `FAILURE`:
      Other solver failure during re-initialization.

    Parameters
    ----------
    solver : FedemSolver
        Fedem dynamics solver instance
    df : DataFrame
        Input function values
    state : list of float
        State vector to restart simulation from
    lib_dir : str
        Path to the input/output files
    out_id : list of int
        List of user Ids identifying the output sensors in the model
    c_opts : dictionary
        Settings for the Fedem solver
    x_times : list of float, default=None
        Time list linked to the input

    Returns
    -------
    list of float
        Output sensor values for each time step
    str
        How the simulation was actually restarted (or not)
    """

    print("  ** The solver has diverged at t =", solver.get_current_time())

    # If file stream.json not existing, create it
    json_file = lib_dir + "/stream.json"
    if not path.isfile(json_file):
        with open(json_file, "w") as fd:
            fd.write(dumps({"conv_type": "TOL"}, indent=2))

    # Define state path for state file storage
    step_state_path = lib_dir + "/step_state"
    if not path.isdir(step_state_path):
        mkdir(step_state_path)

    # Input data
    input_data = df.values
    n_step = df.shape[0]

    # Default dictionary settings (solver)
    # Check that the option files exist before adding them
    opts = {"cwd": getcwd(), "terminal": -1}
    for ext in ("fco", "fop", "fao"):
        option_file = "fedem_solver." + ext
        if path.isfile(getcwd() + "/" + option_file):
            opts[ext] = option_file

    # Do not overwrite the original res-file
    ires = 0
    while ires <= 0:
        ires -= 1
        opts["resfile"] = "fedem_solver_r" + str(-ires) + ".res"
        if not path.isfile(getcwd() + "/" + opts["resfile"]):
            ires = -ires

    if c_opts is not None:
        opts.update(c_opts)

    # Convergence handling settings (solver)
    activate_ramp = opts.pop("use_ramping", True)
    state_pos = opts.pop("use_state_n_steps_behind", 4)
    nr_steps = opts.pop("nr_steps", 10)
    conv_types = opts.pop(
        "sequence",
        [
            "TOL",
            "STATE_DYNAMIC",
            "STATE_STATIC",
            "NO_STATE_WITH_RAMP",
            "JUMP_OVER",
            "NO_CONV",
        ],
    )
    # Ensure the trial loop exits with NO_CONV if nothing else succeeds
    if conv_types[-1] != "NO_CONV":
        conv_types.append("NO_CONV")

    use_ramp = False

    j = 0
    k = 0
    setoff = 0  # setoff for scaling function
    output = [0.0] * n_step * len(out_id)

    # Set convergence type to TOL
    conv_type = conv_types[0]

    # Initialise step state file size
    state_size = solver.get_state_size()
    step_state = state

    for k in range(len(conv_types)):
        # Shut down the solver
        solver.solver_done()

        # Reinitialize the solver
        ierr = solver.solver_init(_opts_list(opts), state_data=step_state)
        if ierr < 0:
            print(
                f" *** Solver initialization failure in divergence handling ({ierr})",
                flush=True,
            )
            return None, "FAILURE"

        print(f"     Trying restart {conv_type} at t = {solver.get_current_time()}")

        # Scaling down the input for the remaining period
        if use_ramp:
            input_data = ramp_dataframe(df, j - setoff).values
            use_ramp = False

        while j < n_step and solver.ierr.value == 0:
            # Solve for next time step
            t_next = None if x_times is None else x_times[j]
            if solver.solve_next(input_data[j], time_next=t_next):
                # Save state during stepping
                solver.save_state()
                if not solver.state_data:
                    print("  ** No state at step", j, flush=True)
                else:
                    step_state_file = step_state_path + "/step_" + str(j % state_pos)
                    with open(step_state_file, "wb") as step_file:
                        try:
                            step_file.write(solver.state_data)
                        except IOError:
                            print("  ** Failed to write state to", step_state_file)

                # Extract the results
                output[j * len(out_id) : j * len(out_id) + len(out_id)] = (
                    float(solver.get_function(out_id[i])) for i in range(len(out_id))
                )
                j += 1
            else:
                if solver.ierr.value != 0:
                    print("  ** The solver diverged at t =", solver.get_current_time())
                break

        # Leave loop in case of convergence/non-convergence and jump over
        conv_type = conv_types[k + 1]
        if conv_type == "JUMP_OVER":
            j = n_step
        if j == n_step:
            # Update data on json file
            with open(json_file, "w") as fd:
                fd.write(dumps({"conv_type": conv_type}, indent=2))
            break

        # Convergence issue still exists
        if conv_type == "NO_STATE_WITH_RAMP":
            step_state = None
            use_ramp = True
        elif conv_type == "NO_CONV":
            break  # No convergence reached, giving up
        elif conv_type in ("STATE_DYNAMIC", "STATE_STATIC"):
            if j < state_pos:
                step_state = state  # use state at the beginning
            else:
                step_state = (c_double * state_size)()
                step_state_file = step_state_path + "/step_" + str(j % state_pos)
                print(f"   * Read state from {step_state_file} [{state_size}]")
                with open(step_state_file, "rb") as step_file:
                    try:
                        step_file.readinto(step_state)
                    except FileNotFoundError:
                        print("  ** Failed to read state from", step_state_file)

            # Number of backwards steps
            j = 0 if j < state_pos else j - state_pos + 1

            # Activate ramp (external specification, user request, default activated)
            if activate_ramp:
                use_ramp = True
                setoff = 0  # means ramping starts with zero

            # Equilibrium fails at time, set solver parameter quasiStatic
            if conv_type == "STATE_STATIC":
                time = solver.get_next_time()
                dt = time - solver.get_current_time()
                opts["quasiStatic"] = str(time + nr_steps * dt)

        # Increment the restart res-file, to not overwrite previous ones
        ires += 1
        opts["resfile"] = "fedem_solver_r" + str(ires) + ".res"

    return output, conv_type
