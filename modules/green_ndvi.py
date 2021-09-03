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
from modules.green_osm import calculate_probability_mass, read_ndvi


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


def green_from_ndvi(config):
    """
    Fuses evidence from multiple sources (OSM, NDVI, OSM context)
    :param config:
    :return:
    """

    aoi_name = config["name"]
    lu_polygons_file = os.path.join(
        config["output_dir"], aoi_name, f"{aoi_name}_lu_polygons_green.shp"
    )
    ndvi_dir = os.path.join(config["output_dir"], aoi_name, "ndvi")

    logger = logging.getLogger("root." + __name__)

    # Input parameters from config file ------------------
    green_c = config["fuzzy_centers"]["green"]
    mixed_c = config["fuzzy_centers"]["mixed"]
    grey_c = config["fuzzy_centers"]["grey"]
    d = config["fuzzy_centers"]["d"]

    green_old_min = mixed_c - d
    green_old_max = green_c
    grey_old_min = grey_c
    grey_old_max = mixed_c + d

    # Read NDVI file
    ndvi, affine, crs, nodata = read_ndvi(ndvi_dir=ndvi_dir, aoi_name=aoi_name)
    buffer = -5

    # Read Targets
    logger.info("Reading data...")
    lu_polygons = gpd.read_file(lu_polygons_file)

    # Delete columns for beliefs if they already exist
    if "g_ndvi" in lu_polygons.columns:
        lu_polygons.drop(["g_ndvi", "n_ndvi", "gn_ndvi"], axis=1, inplace=True)

    targets_df_small = lu_polygons.copy().to_crs({"init": str(crs)})
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
    green_ndvi = pd.DataFrame(prob_masses).T
    green_ndvi.columns = ["g_ndvi", "n_ndvi", "gn_ndvi"]

    lu_polygons_ndvi = lu_polygons.join(green_ndvi, how="left")
    lu_polygons_ndvi = pd.DataFrame(lu_polygons_ndvi).loc[
        ~lu_polygons_ndvi.index.isna()
    ]

    # Set belief of polygons which do not contain any reliable information about greenness from the sentinel-2 image
    lu_polygons_ndvi["g_ndvi"].replace(np.nan, 0, inplace=True)
    lu_polygons_ndvi["n_ndvi"].replace(np.nan, 0, inplace=True)
    lu_polygons_ndvi["gn_ndvi"].replace(np.nan, 1, inplace=True)
    lu_polygons_ndvi = lu_polygons_ndvi.round(2)

    gpd.GeoDataFrame(lu_polygons_ndvi).to_file(lu_polygons_file)


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
