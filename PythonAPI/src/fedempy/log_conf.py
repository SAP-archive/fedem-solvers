# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

"""
Logger
"""
# -*-coding:utf-8-*-

import logging

use_logger = False  # Set to True to activate this logger


class DummyLogger:
    """
    Dummy logger class doing nothing.
    """

    def info(self, message):
        """
        info
        """


def get_logger(module_name="log.txt"):
    """
    Get file logger
    """

    if not use_logger:
        return DummyLogger()

    # logger
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)

    # create console handler and set level to debug
    f_handler = logging.FileHandler(module_name)
    f_handler.setLevel(logging.DEBUG)

    # create formatter
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    # add formatter to f_handler
    f_handler.setFormatter(formatter)

    # add ch to logger
    logger.addHandler(f_handler)

    return logger


# end of file
