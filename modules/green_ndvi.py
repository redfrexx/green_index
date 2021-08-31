#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Preprocessing of imagery for greenness model"""

__author__ = "Christina Ludwig, GIScience Research Group, Heidelberg University"
__email__ = "christina.ludwig@uni-heidelberg.de"


import pandas as pd
import geopandas as gpd
import logging
import os
import numpy as np
from shapely.geometry import box
import json
from rasterio import features
import glob
import rasterio as rio

import modules.utils as utils
from modules.green_osm import calculate_probability_mass


def read_ndvi(in_dir=None, aoi_name=None, ndvi_file=None):
    """
    Searches and reads in ndvi file
    :return:
    """
    if ndvi_file is None:
        # Search NDVI file
        ndvi_files = glob.glob(os.path.join(in_dir, "ndvi", aoi_name + "_ndvi_s2.tif"))
        if len(ndvi_files) == 0:
            return 1, "No NDVI file found."
        elif len(ndvi_files) > 1:
            return 1, "Too many matching NDVI files found."
        ndvi_file = ndvi_files[0]

    # Open NDVI file
    with rio.open(ndvi_file) as src:
        ndvi = src.read(1)
        affine = src.transform
        crs = src.crs
        nodata = src.nodatavals[0]

    return 0, [ndvi, affine, crs, nodata]


def calc_ndvi_stats(targets_df, ndvi, affine):
    """
    Get NDVI mean value for each target part
    :param targets_df:
    :return:
    """

    # targets_df_reproj = targets_df.copy().to_crs(crs)

    # targets_df_reproj.geometry = targets_df_reproj.geometry.buffer(buffer)

    # EXTRACT NDVI VALUES FOR POLYGONS
    shapes = (
        (geom, id)
        for geom, id in zip(targets_df.geometry.values, targets_df.index.to_numpy())
    )
    image = features.rasterize(shapes, out_shape=ndvi.shape, transform=affine, fill=-1)

    ids = targets_df.index.to_numpy()
    ndvi_values = [ndvi[image == id] for id in ids]
    # ndvi_mean, ndvi_median = zip(*[(bn.nanmean(vals), bn.nanmedian(vals)) for vals in ndvi_values])
    # calculate mean and std of NDVI
    # ndvi_df = pd.DataFrame({"ndvi_mean": ndvi_mean, "ndvi_median": ndvi_median}, index=target_parts_ndvi_df.index, dtype="float16")
    # ndvi_df = {id:vals for id, vals in zip(ids, ndvi_values)}

    return ids, ndvi_values


def green_from_ndvi(aoi_name, config, sensor="s2"):
    """
    Fuses evidence from multiple sources (OSM, NDVI, OSM context)
    :param config:
    :return:
    """

    logger = logging.getLogger("root." + __name__)

    # Input parameters from config file ------------------
    polygon_file = os.path.abspath(config["aois"][aoi_name]["target_geoms"])
    out_dir_green = os.path.join(config["output_dir"], aoi_name, "green")

    green_c = config["aois"][aoi_name]["fuzzy_centers"][sensor]["green"]
    mixed_c = config["aois"][aoi_name]["fuzzy_centers"][sensor]["mixed"]
    grey_c = config["aois"][aoi_name]["fuzzy_centers"][sensor]["grey"]
    d = config["aois"][aoi_name]["fuzzy_centers"][sensor]["d"]

    green_old_min = mixed_c - d
    green_old_max = green_c
    grey_old_min = grey_c
    grey_old_max = mixed_c + d

    # Read NDVI file
    if sensor == "s2":
        success, result = read_ndvi(in_dir=out_dir_green, aoi_name=aoi_name)
        buffer = -5
    elif sensor == "vhr":
        success, result = read_ndvi(ndvi_file=config["aois"][aoi_name]["vhr"])
        buffer = -0.5
    else:
        print("unknown sensor")

    if success == 1:
        logger.critical(result)
        exit(1)
    ndvi, affine, crs, nodata = result

    # Bounding box for reading polygons
    bbox_geom = box(*config["aois"][aoi_name]["bbox"])

    # Read Targets
    logger.info("Reading data...")
    targets_df = gpd.read_file(polygon_file, bbox=bbox_geom)
    targets_df = targets_df.loc[:, ("TARGET_ID", "geometry")]

    # Clip to bounding box
    bbox_df = gpd.GeoDataFrame({"geometry": [bbox_geom]}, crs={"init": "epsg:4326"})
    targets_df = gpd.overlay(bbox_df, targets_df, how="intersection")

    targets_df_small = targets_df.copy().to_crs({"init": str(crs)})
    targets_df_small["geometry"] = targets_df_small.buffer(buffer)
    targets_df_small = targets_df_small[
        targets_df_small.geometry.is_valid & ~targets_df_small.geometry.is_empty
    ]

    # Calculate NDVI for each target part
    logger.info("Extracting NDVI values...")
    ids, ndvi_values = calc_ndvi_stats(targets_df_small, ndvi, affine)

    logger.info("Calculating probability masses...")
    prob_masses = {
        id: calculate_probability_mass(
            vals,
            green_old_min=green_old_min,
            green_old_max=green_old_max,
            grey_old_min=grey_old_min,
            grey_old_max=grey_old_max,
        )
        for id, vals in zip(ids, ndvi_values)
    }
    targets_ndvi = pd.DataFrame(prob_masses).T
    targets_ndvi.columns = ["green", "grey", "green_grey"]

    targets_out = targets_df.loc[:, ["TARGET_ID"]].join(targets_ndvi, how="left")
    targets_out = pd.DataFrame(targets_out).loc[~targets_out.index.isna()]
    targets_out["green"].replace(np.nan, 0, inplace=True)
    targets_out["grey"].replace(np.nan, 0, inplace=True)
    targets_out["green_grey"].replace(np.nan, 1, inplace=True)

    targets_out.set_index("TARGET_ID", inplace=True)
    targets_out = targets_out.round(2)

    # Export to file
    logger.info("Writing to file...")
    outfile = os.path.join(
        out_dir_green,
        aoi_name + "_green_{0}.csv".format(os.path.basename(sensor).split(".")[0]),
    )
    targets_out.to_csv(outfile)


if __name__ == "__main__":
    aoi_name = "heidelberg"
    config_file = "../config/config_meingruen.json"

    # Read configuration parameters
    with open(config_file, "r") as src:
        config = json.load(src)

    bbox = config[aoi_name]["bbox"]

    # SET UP LOGGING
    # ---------------
    logger = utils.setup_logger(config["log_dir"], aoi_name)

    green_from_ndvi(aoi_name, config)
