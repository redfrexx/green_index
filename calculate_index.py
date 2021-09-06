#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Calculates index of highways for usage openrouteservice as extended storage"""

__author__ = "Christina Ludwig, GIScience Research Group, Heidelberg University"
__email__ = "christina.ludwig@uni-heidelberg.de"

import os
import sys
import argparse
from modules.download import download_features
from modules.green_index import calc_green_index

from modules.utils import load_config, init_logger


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Calculates the index of each OSM highway based on provided raster or "
        "vector file. The highways are downloaded using the ohsome API."
    )
    parser.add_argument(
        "--bbox",
        "-b",
        required=True,
        dest="bbox",
        type=str,
        help="Bounding box in geographic coordinates as string without whitespace e.g. 'minx,miny,maxx,maxy')",
    )
    parser.add_argument(
        "--timestamp",
        "-t",
        required=False,
        dest="timestamp",
        type=str,
        default=None,
        help="ISO formatted timestamp for download of highways from OSM.",
    )
    parser.add_argument(
        "--width",
        "-w",
        required=True,
        dest="width",
        type=float,
        help="Width of the buffer around highway segments in meters",
    )
    parser.add_argument(
        "--vector",
        "-v",
        required=False,
        dest="vector_file",
        type=str,
        help="Path to vector file containing features counted nearby highways",
    )
    parser.add_argument(
        "--raster",
        "-r",
        required=False,
        dest="raster_file",
        type=str,
        help="Path to raster file used to calculate mean value within area nearby highway",
    )
    parser.add_argument(
        "--outputdirectory",
        "-o",
        required=True,
        dest="output_dir",
        type=str,
        help="Path to existing output directory.",
    )
    args = parser.parse_args()

    logger = init_logger("calculate index")

    try:
        assert (
            args.vector_file is not None or args.raster_file is not None
        ), "Either raster or vector file must be given."
        assert not (
            (args.vector_file is not None) and (args.raster_file is not None)
        ), "Either raster or vector file must be given."
    except AssertionError as e:
        logger.critical(e)
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    if not os.path.exists(os.path.join(args.output_dir, "highways.geojson")):
        logger.info("Downloading highway features...")
        config_file_tags = "./config/config_tags.json"
        config_tags = load_config(config_file_tags)
        try:
            download_features(
                bbox=args.bbox,
                timestamp=args.timestamp,
                layers=config_tags["highways"],
                outdir=args.output_dir,
            )
        except Exception:
            logger.exception("Error during download of highways:")
            sys.exit(1)

    try:
        logger.info("Calculating green index...")
        calc_green_index(args.width, args.output_dir, raster_file=args.raster_file)
    except Exception:
        logger.exception("Error during green index calculation:")
        sys.exit(1)
