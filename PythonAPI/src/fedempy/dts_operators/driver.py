# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

"""
Convenience drivers for running FEDEM simulations as an operator.

These drivers can also take the model input as low-code (python) or
no-code (yaml) files, and the corresponding fmm-file of the model
will then be generated the first time the operator is invoked.
"""

import fedempy.dts_operators.fmm_creator as lownocode
import fedempy.dts_operators.stress_visualization as fpy_stress
import fedempy.dts_operators.window as fpy_window


def run_window(df, dts=None, **kwargs):
    """
    Run fedem simulation with `df` as input.

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
    Dataframe
        Response values in output sensors
    """

    # Check if solver_options is specified, if yes run fedem_solver directly
    if "solver_options" in kwargs:
        return fpy_window.run(df, dts, **kwargs)

    # Find or create the fmm-file
    fmm_path = lownocode.create_fmm(
        kwargs.get("lib_dir", "/var/digitaltwin/app/lib"),
        kwargs.get("fmm_file", None),
        kwargs.get("model_file", None),
    )

    # Run the dynamics solver
    return fpy_window.run(df, dts, fmm_file=fmm_path, **kwargs)


def run_stress_visualization(df, **kwargs):
    """
    Run batch fedem simulation with `df` as input,
    to create CUG stress visualization of the FE-parts.

    Parameters
    ----------
    df : DataFrame
        Input function values
    kwargs : dict
        Dictionary containing output definitions and/or solver options
    """

    # Check if solver_options is specified, if yes run fedem_solver directly
    if "solver_options" in kwargs:
        return fpy_stress.run(df, **kwargs)

    # Find or create the fmm-file
    fmm_path = lownocode.create_fmm(
        kwargs.get("lib_dir", "/var/digitaltwin/app/lib"),
        kwargs.get("fmm_file", None),
        kwargs.get("model_file", None),
    )

    # Run the dynamics solver using the given/generated fmm-file
    return fpy_stress.run(df, fmm_file=fmm_path, **kwargs)
