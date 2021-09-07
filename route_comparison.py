#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Comparison of green and short routes using openrouteservice"""

__author__ = "Christina Ludwig, GIScience Research Group, Heidelberg University"
__email__ = "christina.ludwig@uni-heidelberg.de"


import argparse
import sys
from modules.mapmatching import mapmatching
from modules.utils import load_config, init_logger
from modules.routing import generate_routes


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Comparison of short and green routes using openrouteservice"
    )
    parser.add_argument(
        "--config",
        "-c",
        required=True,
        dest="config_file",
        type=str,
        help="Configuration file",
    )
    args = parser.parse_args()
    logger = init_logger("comparison")
    config = load_config(args.config_file)

    try:
        generate_routes(config)
    except Exception:
        logger.exception("Error during simulation of routes:")
        sys.exit(1)

    try:
        mapmatching(config)
    except Exception:
        logger.exception("Error during map matching:")
        sys.exit(1)
