# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

"""
Convenience functions for running FEDEM sub-model simulations as an operator.
To use, call sub_model_run(df) with input data in df.
"""

import glob
import json
import logging
import os
import uuid
from subprocess import CalledProcessError, run

import numpy as np
import pandas as pd

from fedempy.solver import FedemSolver

log = logging.getLogger(__name__)


def _common_args(model_dir, out_dir, result_id):
    solver_args = [
        "-cwd=" + model_dir,
        "-allPrimaryVars-",
        "-allSecondaryVars-",
        "-terminal=-1",
    ]
    if os.path.isfile(model_dir + "/fedem_solver.fsi"):
        solver_args.append("-fsifile=" + model_dir + "/fedem_solver.fsi")
    if os.path.isfile(model_dir + "/fedem_solver.fco"):
        solver_args.append("-fco=" + model_dir + "/fedem_solver.fco")
    if os.path.isfile(model_dir + "/fedem_solver.fop"):
        solver_args.append("-fop=" + model_dir + "/fedem_solver.fop")
    if os.path.isfile(model_dir + "/fedem_solver.fao"):
        solver_args.append("-fao=" + model_dir + "/fedem_solver.fao")
    solver_args.append("-resfile=" + out_dir + "/" + result_id + ".res")

    return solver_args


def _recovery_args(out_dir, frs_file):
    return [
        "-frs3file=" + out_dir + "/" + frs_file,
        "-recovery=1",
        "-partDeformation=2",
        "-partVMStress=0",
    ]


def _solve_global(model_dir, out_dir, input_data, n_steps, n_funcs):
    """
    Solve top-level global FEDEM model.
    input_data must contain n_funcs columns,
    mapped by index to external functions in the fedem model.
    """
    input_def = range(1, n_funcs + 1)
    result_id = str(uuid.uuid4())
    solver_args = _common_args(model_dir, out_dir, result_id)
    solver_args.extend(_recovery_args(out_dir, result_id + "_Recovery_glob.frs"))

    solver = FedemSolver(os.environ["FEDEM_SOLVER"], solver_args)
    for i in range(n_steps):
        solver.solve_next(input_data[i], input_def)
    solver.solver_done()
    del solver

    return result_id


def _solve_sub_level(model_dir, out_dir, parent_frsfile, n_steps, conf, output):
    """
    Solve sub-level model. Recover strain-gauges if bottom-level
    """
    result_id = str(uuid.uuid4())

    # Run solmap
    map_file = glob.glob(model_dir + "/*.map")[0]
    fnd_file = out_dir + "/" + result_id + "_interp.fnd"
    f_solver = os.path.normpath(os.environ["FEDEM_SOLVER"])
    f_solmap = os.path.dirname(f_solver) + "/fedem_solmap"

    try:
        process = run(
            [
                f_solmap,
                "-frsFile=" + parent_frsfile,
                "-outFile=" + fnd_file,
                "-mapFile=" + map_file,
            ],
            capture_output=True,
            check=True,
        )
    except CalledProcessError as err:
        log.exception("fedem_solmap exited with return code: %d", err.returncode)
        log.debug("fedem_solmap stdout: %r", err.stdout)
        log.debug("fedem_solmap stderr: %r", err.stderr)
        raise
    else:
        log.debug("fedem_solmap stdout: %r", process.stdout)
        log.debug("fedem_solmap stderr: %r", process.stderr)

    subfolders = [f.path for f in os.scandir(model_dir) if f.is_dir()]
    # If subfolder(s), recover stress and run sub-levels, if not recover gauges.
    if len(subfolders) > 0:
        solver_args = _common_args(model_dir, out_dir, result_id)
        solver_args.append("-displacementfile=" + fnd_file)
        solver_args.extend(_recovery_args(out_dir, result_id + "_Recovery_sub.frs"))

        solver = FedemSolver(os.environ["FEDEM_SOLVER"], solver_args)
        for i in range(n_steps):
            solver.solve_next()
        solver.solver_done()
        del solver

        recovery_file_name = glob.glob(out_dir + "/" + result_id + "*Recovery*")[0]

        for sub_folder in subfolders:
            _solve_sub_level(
                sub_folder, out_dir, recovery_file_name, n_steps, conf, output
            )
    else:
        # If a bottom level of submodels is reached, solve strain gauge recovery
        folder_name = os.path.basename(os.path.normpath(model_dir))
        solver_args = _common_args(model_dir, out_dir, result_id)
        solver_args.extend(
            [
                "-displacementfile=" + fnd_file,
                "-recovery=2",
                "-frs3file=" + out_dir + "/" + result_id + "_Recovery_sub_gage.frs",
            ]
        )

        solver = FedemSolver(os.environ["FEDEM_SOLVER"], solver_args)
        for i in range(n_steps):
            solver.solve_next()
            for vals in conf[folder_name]:
                output[i, vals[0]] = solver.get_function(vals[1])
        solver.solver_done()
        del solver


def sub_model_run(df):
    """
    Run sub-model simulation.
    First solves the global model, then recursively solves the submodels.
    Reads external submodel_config.json file located in lib-folder.

    Parameters
    ----------
    df : Pandas dataframe
        Dataframe containing input data
    """
    lib_dir = "/var/digitaltwin/app/lib"
    if not os.path.isdir(lib_dir):
        lib_dir = os.getcwd()

    # setup input
    n_steps, n_inputs = df.shape
    global_input = df.values

    # setup output
    conf = {}
    with open(lib_dir + "/submodel_config.json") as conf_file:
        conf = json.load(conf_file)

    # Loop over list of output-functions for the submodels.
    # Create output_defs as dict with output index as key, and function name as value
    output_defs = {}
    for output_functions in conf.values():
        for func in output_functions:
            output_defs[func[0]] = func[2]

    output = np.zeros((n_steps, len(output_defs)))

    # Set-up folder for temporary frs-files
    fedem_out_dir = "/var/digitaltwin/app/lib/tmp"
    if not os.path.exists(fedem_out_dir):
        os.mkdir(fedem_out_dir)

    # Clean
    files = glob.glob(fedem_out_dir + "/*")
    for temp_file in files:
        os.remove(temp_file)

    # Solve global model
    result_id_global = _solve_global(
        lib_dir + "/model", fedem_out_dir, global_input, n_steps, n_inputs
    )
    recovery_file_name = glob.glob(
        fedem_out_dir + "/" + result_id_global + "*Recovery*"
    )[0]

    # Solve sublevels recursively
    subfolders = [f.path for f in os.scandir(lib_dir + "/model") if f.is_dir()]
    for sub_folder in subfolders:
        _solve_sub_level(
            sub_folder, fedem_out_dir, recovery_file_name, n_steps, conf, output
        )

    # Create output dataframe
    column_names = []
    for key in sorted(output_defs):
        column_names.append(output_defs[key])

    df_out = pd.DataFrame(output, index=df.index, columns=column_names)

    return df_out
