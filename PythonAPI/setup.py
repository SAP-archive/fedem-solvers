# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

"""
Setuptools package definition
"""
from setuptools import find_packages, setup

# build pipeline updates version.txt with build timestamp
VERSION = open("version.txt").read().rstrip()
# "A Better Pip Workflow": https://www.kennethreitz.org/essays/a-better-pip-workflow
REQUIRES = list(open("requirements-to-freeze.txt")) + list(
    open("requirements-internal.txt")
)
DEV_REQUIRES = list(open("requirements-dev.txt"))


setup(
    name="fedempy",
    version=VERSION,
    description="Python wrapper for the FEDEM solvers",
    author="Knut Morten Okstad, SAP SE",
    author_email="Knut.Morten.Okstad@sap.com",
    url="www.sap.com",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    include_package_data=True,
    license="SAP License",
    install_requires=REQUIRES,
    extras_require={"dev": DEV_REQUIRES},
    classifiers=["Intended Audience :: Internal", "Programming Language :: Python"],
    test_suite="test",
)
