# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

"""
Utility module, for generating fmm-files from low/no-code model files.
"""

from glob import glob
from importlib import import_module
from inspect import isfunction
from os import getcwd, path
from pathlib import Path

from fedempy.yaml_parser import ModelYAML


def _abs_path(file_name, dir_name):
    """Returns the absolute path of a file."""
    if file_name:
        if not path.isabs(file_name):
            file_name = dir_name + "/" + file_name
        if not path.isfile(file_name):
            file_name = None
    return file_name


def _call_func(full_module_name, func_name, *argv):
    """Generic function invoking a low-code modeling script."""
    module = import_module(full_module_name)
    for attribute_name in dir(module):
        attribute = getattr(module, attribute_name)
        if isfunction(attribute) and attribute_name == func_name:
            attribute(*argv)


def create_fmm(lib_dir, fmm_file=None, mod_file=None):
    """
    Creates a Fedem model file (fmm) based on a low-code (python)
    or no-code (yaml) model file. If the `fmm_file` is specified
    directly, that file will be used instead. If the `lib_dir` folder
    already contains fmm-files, the first one found will be used.
    This will be the case when starting a new micro-batch in a streaming app.

    Parameters
    ----------
    lib_dir : str
        Absolute path to directory where the fmm-file should be stored
    fmm_file : str, default=None
        Name of fmm-file, if provided explicitly
    mod_file : str, default=None
        Name of low/no-code model file to generate fmm-file from

    Returns
    -------
    str
        Absolute path to (generated) fmm-file
    """

    if not path.isdir(lib_dir):
        lib_dir = getcwd()  # use current working directory if no lib_dir

    # If fmm_file is specified explicitly, use it if the file exists
    fmm_file = _abs_path(fmm_file, lib_dir)
    if fmm_file:
        return fmm_file

    # Check if lib_dir already contains fmm-files
    fmm_files = glob(lib_dir + "/*.fmm")
    if not fmm_files and mod_file:
        # No fmm-files found, but a low/no-code model file was specified
        mod_file = _abs_path(mod_file, lib_dir)
        if not mod_file:
            return None

        # Extract the file extension
        ext = mod_file[mod_file.rfind(".") :]

        if ext == ".py":
            # Create fmm-file from low-code (py-file)
            _call_func(Path(mod_file).stem, "model_builder", lib_dir)
        elif ext in {".yaml", ".yml"}:
            # Create fmm-file from no-code (yaml-file)
            model = ModelYAML(mod_file)
            model.build()
            model.model.close(True)

        fmm_files = glob(lib_dir + "/*.fmm")

    # Return the (first) fmm-file found (should only be one)
    return path.abspath(fmm_files[0]) if fmm_files else None
