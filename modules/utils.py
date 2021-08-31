#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""General utility functions"""

__author__ = "Christina Ludwig, GIScience Research Group, Heidelberg University"
__email__ = "christina.ludwig@uni-heidelberg.de"


import logging
import os
import json
import pygeos
import numpy as np
import math


def create_subfolder(out_dir: str, name: str):
    """
    Creates a subfolder in out_dir with given name
    :param out_dir: Output directory
    :param name: Name of new subfolder
    :return:
    """
    new_subfolder = os.path.join(out_dir, name)
    os.makedirs(new_subfolder, exist_ok=True)
    return new_subfolder


def check_config(config, aoi_name):
    """
    Check for missing parameters in config file
    :return:
    """
    # Check general parameters
    general_parameters = ["output_dir", "init_timestamp"]
    for par in general_parameters:
        assert par in config.keys(), f"Parameter '{par}' missing in config file."
    assert aoi_name in config["aois"], f"AOI '{aoi_name}' not found in config file"

    # Check AOI parameters
    aoi_parameters = ["target_geoms", "timestamp"]
    for par in aoi_parameters:
        assert (
            par in config["aois"][aoi_name].keys()
        ), f"Parameter 'aois/{aoi_name}/{par}' missing in config file."


def init_logger(name, log_file_name=None):
    """
    Set up a logger instance with stream and file logger
    :return:
    """
    assert os.path.exists(
        os.path.dirname(log_file_name)
    ), "Error during logger setup: Directory of log file does not exist."
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S",
    )
    # Add stream handler
    streamhandler = logging.StreamHandler()
    streamhandler.setLevel(logging.INFO)
    streamhandler.setFormatter(formatter)
    logger.addHandler(streamhandler)
    # Log file handler
    if log_file_name:
        filehandler = logging.FileHandler(filename=log_file_name)
        filehandler.setLevel(logging.INFO)
        filehandler.setFormatter(formatter)
        logger.addHandler(filehandler)
    return logger


def load_config(config_file):
    """
    Load config file
    :return:
    """
    assert os.path.exists(
        config_file
    ), f"Configuration file does not exist: {os.path.abspath(config_file)}"
    with open(config_file, "r") as src:
        config = json.load(src)
    return config


def calc_compactness(features):
    """
    Selects sliver polygons
    :return:
    """
    pygeos_geoms = features.apply(
        lambda x: pygeos.from_shapely(x["geometry"]), axis=1
    ).values
    area = pygeos.area(pygeos_geoms)
    perimeter = pygeos.length(pygeos_geoms)
    return (area * 4 * np.pi) / (perimeter * perimeter)


def build_ohsome_filters(filter_items: list):
    """
    Builds the filter string based a dictionary of tags
    :return:
    """
    filters = {}
    for item in filter_items:
        geom_filter = " or ".join(["geometry:{0}".format(i) for i in item["geoms"]])
        for key, values in item["tags"].items():
            if isinstance(values, list):
                tag_filter = "{0} in ({1})".format(key, ", ".join(values))
            elif isinstance(values, str):
                tag_filter = f"{key}={values}"
            else:
                raise TypeError(f"{type(values)} not supported for tag values")
            filters[item["name"]] = f"({tag_filter}) and ({geom_filter})"
    return filters


def split_bbox(bbox, tile_size):
    """
    Split bounding box in tiles

    :param bbox: bounding box in format [xmin, ymin, xmax, ymax]
    :param tile_size: the size of the tiles in x and y direction in crs coordinates in format (x_tilesize, y_tilesize)
    :return:
    """
    x_min, y_min, x_max, y_max = bbox
    dx, dy = tile_size

    # Number of full tiles in x and y direction
    x_tiles = math.floor(round(x_max - x_min, 6) / dx)
    y_tiles = math.floor(round(y_max - y_min, 6) / dy)

    # Remainder of bbox in x and y direction
    x_rest = round(x_max - x_min, 6) % dx
    y_rest = round(y_max - y_min, 6) % dy

    for y in range(0, y_tiles):
        for x in range(0, x_tiles):
            yield tuple(
                np.array(
                    [
                        x_min + dx * x,
                        y_min + dy * y,
                        x_min + dx * (x + 1),
                        y_min + dy * (y + 1),
                    ]
                ).round(6)
            )
        if x_rest != 0:
            yield tuple(
                np.array(
                    [x_min + dx * (x + 1), y_min + dy * y, x_max, y_min + dy * (y + 1)]
                ).round(6)
            )

    # Last row
    if y_rest != 0:
        for x in range(0, x_tiles):
            yield tuple(
                np.array(
                    [x_min + dx * x, y_min + dy * y_tiles, x_min + dx * (x + 1), y_max]
                ).round(6)
            )
        yield tuple(
            np.array([x_min + dx * x_tiles, y_min + dy * y_tiles, x_max, y_max]).round(
                6
            )
        )
