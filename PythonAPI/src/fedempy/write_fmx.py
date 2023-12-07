# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

"""
Python wrapper for native FMX-writer library.
"""

import os
from ctypes import byref, c_char_p, c_double, c_int, cdll

if "FMXWRITER" in os.environ:
    _lib = cdll.LoadLibrary(os.environ["FMXWRITER"])
    if os.name == "nt":  # Windows
        writeFMX = _lib.WRITEFMX
        readFMX = _lib.READFMX
    else:  # Linux
        writeFMX = _lib.writefmx_
        readFMX = _lib.readfmx_


def write(fnam, ityp, data):
    """
    This function writes a rectangular matrix as a binary FMX-file for FEDEM.

    Parameters
    ----------
    fnam : str
        Absolute path to the fmx-file to be written
    ityp : int
        Type of matrix to write,
        1=stiffness matrix, 2=mass matrix, 3=gravity force vectors
    data : list of float
        Matrix content, column-wise storage

    Returns
    -------
    int
        Zero on success, otherwise negative
    """

    cfil = c_char_p(fnam)
    ctyp = c_int(ityp)
    cdat = (c_double * len(data))()
    cdat[:] = data
    clen = c_int(len(data))
    return writeFMX(cfil, byref(ctyp), cdat, byref(clen), len(fnam))


def read(fnam, ityp, data):
    """
    This function reads a rectangular matrix from a binary FMX-file for FEDEM.

    Parameters
    ----------
    fnam : str
        Absolute path to the fmx-file to be read
    ityp : int
        Type of matrix to read,
        1=stiffness matrix, 2=mass matrix, 3=gravity force vectors
    data : list of float
        Matrix content, column-wise storage

    Returns
    -------
    int
        Zero on success, otherwise negative
    """

    cfil = c_char_p(fnam)
    ctyp = c_int(ityp)
    cdat = (c_double * len(data))()
    clen = c_int(len(data))
    ierr = readFMX(cfil, byref(ctyp), cdat, byref(clen), len(fnam))
    data[:] = cdat
    return ierr
