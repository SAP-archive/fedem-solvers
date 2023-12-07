# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

"""
Tests for fmx-writing
"""

from fedempy.write_fmx import write, read

fnam = b"testfil"
data = [1.0, 2.0, 3.0, 4.0]

# Write stiffness matrix
status = write(fnam, 1, data)
# Write masse matrix
status += write(fnam, 2, data)
# Write gravitation forces
status += write(fnam, 3, data)
if status < 0:
    exit(status)

# Read in the stiffness matrix
new_data = [0.0, 0.0, 0.0, 0.0]
status += read(fnam, 1, new_data)
for i in range(0, 4):
    if data[i] != new_data[i]:
        status -= 1
if status < 0:
    exit(status)

# Read in the mass matrix
new_data = [0.0, 0.0, 0.0, 0.0]
status += read(fnam, 2, new_data)
for i in range(0, 4):
    if data[i] != new_data[i]:
        status -= 1
if status < 0:
    exit(status)

# Read in the gravitation forces
new_data = [0.0, 0.0, 0.0, 0.0]
status += read(fnam, 3, new_data)
for i in range(0, 4):
    if data[i] != new_data[i]:
        status -= 1

exit(status)
