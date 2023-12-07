# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

"""
Convenience module for running Fedem stress recovery with CUG export as an operator.
Usage: Invoke `run(df)` with the input data in the DataFrame `df`.
"""

from ctypes import c_double
from os import environ, makedirs, path
from sys import stdout

from fedempy.dts_operators.window import start_fmm_solver, start_solver
from fedempy.exporter import Exporter
from fedempy.solver import FedemException


def _get_out_dir():
    out_dir_env = environ.get("PEI_SERVICE_FILE_STORAGE")
    if out_dir_env is None:
        return "."
    out_dir = path.normcase(out_dir_env)
    if not path.isdir(out_dir):
        makedirs(out_dir)
    return out_dir


def run(df, **kwargs):
    """
    Run Fedem simulation with stress recovery and create CUG visualization.

    Parameters
    ----------
    df : Dataframe
        Input function values
    kwargs : dict
        Dictionary containing configuration parameters and/or solver options
    """

    # Absolute path to app location
    lib_dir = kwargs.get("lib_dir", "/var/digitaltwin/app/lib")
    # Input values for initial static equilibrium iterations
    ext_inp = df.iloc[0].values if df.index[0] == 0 else None

    if "solver_options" in kwargs:
        # Start the dynamics solver directly with the solver options provided
        solver, status = start_solver(
            kwargs["solver_options"], lib_dir, ext_input=ext_inp
        )

    else:
        # Assume a Fedem model file is provided.
        # This will then create the solver input files and start the dynamics solver.
        # Any FE models will be reduced first, unless they already have been reduced.
        solver, status = start_fmm_solver(
            kwargs.get("fmm_file", None), lib_dir, ext_input=ext_inp
        )

    if status < 0:
        stdout.flush()
        raise FedemException(f"Failed to start solver ({status}).")

    base_ids = []
    fe_parts = kwargs.get("fe_parts", {})
    vis_parts = kwargs.get("visualization_parts", {})

    for part in fe_parts:
        if part["recovery"]:
            base_ids.append(part["base_id"])

    if path.isdir(lib_dir):
        for part in [*fe_parts, *vis_parts]:
            part["path"] = lib_dir + "/" + part["path"]

    if "VIS_EXPORTER" not in environ:
        raise FedemException("Environment variable VIS_EXPORTER not defined")

    exporter = Exporter(
        fe_parts, vis_parts, environ["VIS_EXPORTER"], kwargs.get("vtfx_file", None)
    )

    if status > 0:
        # The first row of the input dataframe was consumed by initial equilibrium iterations,
        # so drop it from the subsequent dynamics/quasi-static simulation
        df.drop(index=df.index[0], axis=0, inplace=True)
        t0 = "\n     The first input values are consumed by the initial equilibrium analysis."
    else:
        t0 = ""

    t_step = f"solving {df.shape[0]} steps [{df.index[0]},{df.index[-1]}]" + t0
    print("   * Solver successfully started,", t_step, flush=True)
    xtimes = solver.check_times(df.index.values, kwargs.get("use_times", False))
    _solve(df, base_ids, exporter, solver, xtimes, lib_dir, "vtfx_file" in kwargs)


def _get_deform_state_array(bid, solver):
    """
    Allocates deformation state array for an FE part.
    """
    state_size = solver.get_part_deformation_state_size(bid)
    if state_size < 0:
        raise FedemException(f"No deformation state ({state_size}) for FE Part {bid}.")

    return (c_double * state_size)()


def _get_stress_state_array(bid, solver):
    """
    Allocates stress state array for an FE part.
    """
    state_size = solver.get_part_stress_state_size(bid)
    if state_size < 0:
        raise FedemException(f"No stress state ({state_size}) for FE Part {bid}.")

    return (c_double * state_size)()


def _solve(df, fe_parts, exporter, solver, xtimes, lib_dir, skip_cug):
    """
    Solve and recover deformations and von Mises stresses.
    """
    c_transformation = (c_double * solver.get_transformation_state_size())()
    c_deformation = {}
    c_stress = {}
    for bid in fe_parts:
        # Initialize the deformation and stress arrays
        c_deformation[bid] = _get_deform_state_array(bid, solver)
        c_stress[bid] = _get_stress_state_array(bid, solver)

    # Run the solver through the given time window
    n_steps = df.shape[0]
    for i in range(n_steps):
        next_time = xtimes[i] if xtimes is not None else None
        _continue = solver.solve_next(df.values[i], time_next=next_time)
        if solver.ierr.value != 0:
            raise FedemException(f"Failed to solve next step ({solver.ierr.value}).")

        # Get transformation and stress recovery results
        solver.save_transformation_state(c_transformation)
        for bid in fe_parts:
            solver.save_part_state(bid, c_deformation[bid], c_stress[bid])

        # Export the transformations and recovery results
        curr_time = solver.get_current_time()
        exporter.do_step(curr_time, c_transformation, c_deformation, c_stress)

        if not _continue:  # Reached simulation end time
            break

    endtime = solver.get_current_time()
    if solver.solver_done() == 0 and solver.ierr.value == 0:
        print(f"**** Time window finished at t={endtime}", flush=True)
    else:
        print(" *** Dynamics solver failed, probably a divergent model.")
        print("     Check the fedem_solver.res file for details.", flush=True)
        raise FedemException(f"Failed to run solver ({solver.ierr.value}).")

    if skip_cug:  # Don't write CUG database and retain the vtfx-file
        exporter.clean(endtime / n_steps, lib_dir)
    else:  # Write CUG database and delete the temporary vtfx-file
        exporter.clean(endtime / n_steps, lib_dir, _get_out_dir())


def stress_visualization_run(df, **kwargs):
    """
    For backwards compatibility.
    Consider remove and update existing apps using this entry accordingly.
    """
    print("  ** Warning: Calling depreciated method stress_visualization_run().")
    print("              Update your app replacing it by run() instead.")
    run(df, **kwargs)
