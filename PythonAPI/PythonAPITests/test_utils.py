# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

"""
Utility functions to facilitate regression testing of the Python wrapper.
"""

import argparse as ap
import numpy as np


def diff_values(v1, v2, eps):
    """
    Compares two float values with some tolerance.
    """
    diff = abs(v1 - v2)
    refv = abs(v1 + v2) * eps * 0.5
    if diff > refv and diff > eps:
        return diff
    return 0.0


# pylint: disable=invalid-name, consider-using-enumerate
def compare_lists(t, l1, l2, eps=1.0e-12):
    """
    Compares two lists of real values with some tolerance.
    """
    ndiff = 0
    for i in range(len(l1)):
        diff = diff_values(l1[i], l2[i], eps)
        if diff > 0.0:
            print(
                "t =",
                t,
                ": Value",
                i + 1,
                "does not match",
                l1[i],
                l2[i],
                "diff =",
                diff,
            )
            ndiff += 1

    return ndiff


def read_reference_data(file_name):
    """
    Reads time history reference data from a file.
    """
    refsd = []
    rfile = open(file_name, "r")
    for line in rfile:
        if line[0:5] == "#DESC":
            count = 0
            for head in line.split():
                count += 1
                if count == 1:
                    print("#Reference quantities:")
                else:
                    print(count, head)
        elif line[0] == "#":
            print(line.strip())  # comment line
        else:
            refsd.append([float(x) for x in line.split()])

    rfile.close()
    return refsd
