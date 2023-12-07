# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

"""
Convenience module for running Fedem simulations as an operator.
Usage: Invoke `run(df)` with the input data in the DataFrame `df`.
"""
from os import environ, path
from sys import stdout

from numpy import zeros
from pandas import DataFrame

from fedempy.divergence import jump_state, ramp_dataframe, restart
from fedempy.fmm_solver import FmmSolver
from fedempy.solver import FedemException, FedemSolver


def _check_state(dts, use_state):
    """
    Checks whether solver state transfer should be used or not.
    """
    if not use_state or not bool(dts):
        return False, None

    if dts.window is not None and dts.window.overlap > 0:
        raise FedemException("State transfer with overlap not yet supported.")

    if dts.state is None:
        return True, None

    if int(dts.state.iat[0, 0]) < 1:
        # An invalid state array was provided, ignore it and do normal start instead.
        # The first three values are the step number, actual time, and increment size,
        # respectively. At least the first and third value should thus be positive.
        print(f"  ** Ignoring invalid state array {dts.state.head(3)}")
        return True, None

    return True, dts.state.to_numpy()


def start_solver(
    solver_options, cw_dir, use_state=False, old_state=None, ext_input=None, tstart=None
):
    """
    Starts the dynamics solver with the provided command-line options.

    Parameters
    ----------
    solver_options : list of str
        List of command-line arguments that are passed to the solver
    cw_dir : str
        Path to working directory of the solver process
    use_state : bool, default=False
        If True, internal state vectors will be allocated
    old_state : list of float, default=None
        State vector to restart simulation from
    ext_input : list of float, default=None
        Initial external function values, for initial equilibrium iterations
    tstart : float, default=None
        Optional start time of simulation, override setting in model file

    Returns
    -------
    FedemSolver
        The dynamics solver object
    int
        Zero on success, negative values indicate errors
    """

    def format_arg(key, value):
        if value is None:
            return None
        if isinstance(value, bool):
            return f"-{key}+" if value else f"-{key}-"
        return f"-{key}={value}"

    options = [format_arg(key, val) for key, val in solver_options.items()]
    if path.isdir(cw_dir):
        options.append("-cwd=" + cw_dir)
    if tstart is not None:
        options.append("-timeStart=" + str(tstart))

    if "FEDEM_SOLVER" not in environ:
        raise FedemException("Environment variable FEDEM_SOLVER not defined")

    solver = FedemSolver(environ["FEDEM_SOLVER"], None, use_state)
    status = solver.solver_init(options, None, old_state, extf_input=ext_input)
    if status >= 0 and solver.ierr.value != 0:
        status = -abs(solver.ierr.value)

    return solver, status


def start_fmm_solver(
    fmm_file,
    lib_dir,
    use_state=False,
    old_state=None,
    ext_input=None,
    t_start=None,
    keep_old_res=True,
):
    """
    Starts the dynamics solver on the provided Fedem model file.

    Parameters
    ----------
    fmm_file : str
        The Fedem model file to run simulation on
    lib_dir : str
        Path to where the model file is located,
        used only if `fmm_file` is a relative path
    use_state : bool, default=False
        If True, internal state vectors will be allocated
    old_state : list of float, default=None
        State vector to restart simulation from
    ext_input : list of float, default=None
        Initial external function values, for initial equilibrium iterations
    t_start : float, default=None
        Optional start time of simulation, override setting in model file
    keep_old_res : bool, default=True
        Option to not overwrite any existing res-file in the RDB directory

    Returns
    -------
    FmmSolver
        The dynamics solver object
    int
        Zero on success, negative values indicate errors
    """

    if fmm_file and not path.isabs(fmm_file) and path.isdir(lib_dir):
        fmm_file = lib_dir + "/" + fmm_file

    solver = FmmSolver(None, use_state)
    status = solver.start(
        fmm_file,
        keep_old_res,
        True,
        True,
        old_state,
        extf_input=ext_input,
        time_start=t_start,
    )

    if status >= 0 and solver.ierr.value != 0:
        status = -abs(solver.ierr.value)

    return solver, status


def _get_start_state(df, use_times):
    """
    Extracts start condition (time and input values) from given DataFrame.
    """
    if use_times:
        # Extract start time and associated input values from first row,
        # to be applied on the initial equilibirum iterations, if any
        return df.index[0], df.iloc[0].values

    if df.index[0] == 0:
        # First row has time stamp 0.0, assume the associated input values
        # then should be applied on the initial equilibrium iterations, if any
        return None, df.iloc[0].values

    # The first row will be applied as load on the first time increment
    return None, None


def _run_window(solver, n_step, df, out_id, xtimes, old_state, lib_dir, options):
    """
    Runs a new time window, with optional restart on divergence.
    """
    output = solver.solve_window(n_step, df.values.flatten(), out_id, xtimes)[0]
    if solver.ierr.value == 0:
        conv_type = "CONVERGED"
    else:
        output, conv_type = restart(
            solver,
            df,
            old_state,
            lib_dir,
            out_id,
            options,
            xtimes,
        )

    t_stop = solver.get_current_time()
    if conv_type == "FAILURE":
        ierr = -1  # Initialization failure, model is already deallocated
    else:
        ierr = solver.solver_done()
    if ierr == 0 and solver.ierr.value == 0:
        print(f"**** Time window finished at t={t_stop}", flush=True)
    elif conv_type in ("NO_CONV", "FAILURE"):
        print(" *** Dynamics solver failed, probably a divergent model.")
        print("     Check the fedem_solver.res file for details.", flush=True)
        raise FedemException(f"Failed to run solver ({solver.ierr.value}).")
    else:
        print(f"**** Time window stopped [{conv_type}]", flush=True)
        if ierr != 0 or solver.ierr.value != 0:
            print("     with error flags", ierr, solver.ierr.value)

    return output


def run(df, dts=None, **kwargs):
    """
    Run Fedem simulation over a time window with `df` as input.

    Parameters
    ----------
    df : DataFrame
        Input function values
    dts : DTSContext, default=None
        State data to be passed between each micro-batch
    kwargs : dict
        Dictionary containing output definitions and/or solver options

    Returns
    -------
    DataFrame
        Response values in output sensors
    """

    # Absolute path to app location
    lib_dir = kwargs.get("lib_dir", "/var/digitaltwin/app/lib")

    use_times = kwargs.get("use_times", False)

    if jump_state(lib_dir):
        df = ramp_dataframe(df)
        use_state, old_state = True, None
    else:
        use_state, old_state = _check_state(dts, kwargs.get("use_state", False))

    if old_state is None:
        print("\n**** Start solver on initial state", flush=True)
        t_start, ext_inp = _get_start_state(df, use_times)
    else:  # Restart from previous state
        istep = int(old_state[0][0])
        time0 = float(old_state[1][0])
        print(f"\n**** Start solver from step {istep} t={time0}", flush=True)
        t_start = None
        ext_inp = None

    if "solver_options" in kwargs:
        # Start the dynamics solver directly with the solver options provided
        solver, status = start_solver(
            kwargs["solver_options"], lib_dir, use_state, old_state, ext_inp, t_start
        )

    else:
        # Assume a Fedem model file is provided.
        # This will then create the solver input files and start the dynamics solver.
        # Any FE models will be reduced first, unless they already have been reduced.
        solver, status = start_fmm_solver(
            kwargs.get("fmm_file"),
            lib_dir,
            use_state,
            old_state,
            ext_inp,
            t_start,
            kwargs.get("keep_old_res", False),
        )

    if status < 0:
        stdout.flush()
        raise FedemException(f"Failed to start solver ({status}).")

    if status > 0:
        # The first row of the input dataframe was consumed by initial equilibrium iterations,
        # so drop it from the subsequent dynamics/quasi-static simulation
        df.drop(index=df.index[0], axis=0, inplace=True)
        t0 = "\n     The first input values are consumed by the initial equilibrium analysis."
    else:
        t0 = ""

    # Run the solver through the given time window
    n_step = df.shape[0]
    t_step = f"solving {n_step} steps [{df.index[0]},{df.index[-1]}]" + t0
    print("   * Solver successfully started,", t_step, flush=True)
    out_id = kwargs.get("output_ids", None)
    output = _run_window(
        solver,
        n_step,
        df,
        solver.get_function_ids(out_id),
        solver.check_times(df.index.values, use_times),
        old_state,
        lib_dir,
        kwargs.get("crash_options", None),
    )
    if out_id is None:
        return None

    # Create the output DataFrame
    n_out = len(output)
    n_col = len(out_id)
    n_row = int(n_out / n_col)
    print(f"   * Create output DataFrame [{n_col}x{n_row}].")
    np_outputs = zeros(n_out)
    np_outputs[:] = output[:]
    df_out = DataFrame(np_outputs.reshape((n_row, n_col), order="C"), index=df.index)
    df_out.columns = [str(oid) for oid in out_id]

    if use_state and solver.state_data:
        # Store the final state for next time window
        print(f"   * Store final state for next window ({solver.state_size.value})")
        new_state = DataFrame(index=range(solver.state_size.value))
        new_state.insert(0, "state", list(solver.state_data))
        dts.state = new_state

    return df_out
