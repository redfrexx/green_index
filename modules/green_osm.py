#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Preprocessing of OSM data for greenness model"""

__author__ = "Christina Ludwig, GIScience Research Group, Heidelberg University"
__email__ = "christina.ludwig@uni-heidelberg.de"


import pandas as pd
import numpy as np
import geopandas as gpd
import glob
import os
import rasterio as rio
from rasterio import features
import logging


def get_ndvi_values(features_df, ndvi, nodata, affine):
    """
    Extract NDVI values within polygons from raster
    :param features_df:
    :param ndvi:
    :param nodata:
    :param affine:
    :return:
    """

    # features_df = featuresdf[~features_df.geometry.isnull()]
    features_df = features_df[
        features_df.geometry.is_valid & ~features_df.geometry.is_empty
    ]
    if len(features_df) == 0:
        return []

    # Rasterize features
    shapes = ((geom, 2) for geom in features_df.geometry)
    image = rio.features.rasterize(
        shapes, out_shape=ndvi.shape, transform=affine, fill=-1
    )
    return ndvi[(image == 2) & (ndvi != nodata)]


def calculate_probability_mass(
    ndvi_values,
    green_old_min=0.3,
    green_old_max=0.9,
    grey_old_min=0.1,
    grey_old_max=0.7,
):
    """
    Calculate probabilities
    :param tags_df:
    :param ndvi_values:
    :return:
    """
    if len(ndvi_values) == 0:
        return 0, 0, 1
    green = round(
        np.vectorize(
            lambda x: prob_mass_green_from_ndvi(x, green_old_min, green_old_max),
            otypes=["float32"],
        )(ndvi_values).mean(),
        2,
    )
    not_green = round(
        np.vectorize(
            lambda x: prob_mass_grey_from_ndvi(x, grey_old_min, grey_old_max),
            otypes=["float32"],
        )(ndvi_values).mean(),
        2,
    )
    either = round((1 - green - not_green), 2)
    return green, not_green, either


def prob_mass_green_from_ndvi(ndvi, old_min=0.3, old_max=0.9):
    """
    Calculate probability mass for greenness from NDVI values
    :param ndvi:
    :param old_min:
    :param old_max:
    :return:
    """
    if ndvi < old_min:
        return 0
    elif ndvi >= old_max:
        return 1
    else:
        new_max = 1
        new_min = 0
        old_range = old_max - old_min
        new_range = new_max - new_min
        return (((ndvi - old_min) * new_range) / old_range) + new_min


def prob_mass_grey_from_ndvi(ndvi, old_min=0.1, old_max=0.7):
    """
    Calculates probability masses for grey from NDVI values
    :param ndvi:
    :param old_min:
    :param old_max:
    :return:
    """

    # Not Green belief
    if ndvi > old_max:
        return 0
    elif ndvi < old_min:
        return 1
    else:
        new_max = 1
        new_min = 0
        old_range = old_max - old_min
        new_range = new_max - new_min
        return 1 - (((ndvi - old_min) * new_range) / old_range) + new_min


def read_ndvi(ndvi_dir, aoi_name):
    """
    Searches and reads in ndvi file
    :return:
    """
    # Search NDVI file
    ndvi_files = glob.glob(os.path.join(ndvi_dir, aoi_name + "*ndvi.tif"))
    if len(ndvi_files) == 0:
        raise FileNotFoundError("No NDVI file found.")
    elif len(ndvi_files) > 1:
        raise FileNotFoundError(
            "Too many NDVI files found: {0}".format(",".join(ndvi_files))
        )

    # Open NDVI file
    with rio.open(ndvi_files[0]) as src:
        ndvi = src.read(1)
        affine = src.transform
        crs = src.crs
        nodata = src.nodatavals[0]

    return [ndvi, affine, crs, nodata]


def greenness_of_osm_tags(config):
    """
    Compute the greenness of OSM tags

    :param aoi_name:
    :param config:
    :return:
    """
    buffer = -5

    logger = logging.getLogger("root." + __name__)
    logger.info("Deriving belief about greenness of OSM tags ...")

    # Read config parameters
    green_c = config["fuzzy_centers"]["green"]
    mixed_c = config["fuzzy_centers"]["mixed"]
    grey_c = config["fuzzy_centers"]["grey"]
    d = config["fuzzy_centers"]["d"]
    aoi_name = config["name"]
    ndvi_dir = os.path.join(config["output_dir"], config["name"], "ndvi")
    lu_polygons_file = os.path.join(
        config["output_dir"], config["name"], f"{aoi_name}_lu_polygons.shp"
    )
    out_file = os.path.join(
        config["output_dir"], config["name"], aoi_name + "_green_tags.csv"
    )

    # Read NDVI file
    ndvi, affine, crs, nodata = read_ndvi(ndvi_dir, aoi_name)

    # Read landuse polygons
    targets_df = gpd.read_file(lu_polygons_file).to_crs({"init": str(crs)})

    # Unique land use tags
    targets_df["tags"].replace(np.nan, "None", inplace=True)
    osm_tags = targets_df["tags"].unique()

    # Buffer targets inwards
    targets_df.geometry = targets_df.geometry.buffer(buffer)

    # Create empty dataframe which holds greenness of each tag
    tag_greenness = {}

    # Calculate parameters for membership functions
    green_old_min = mixed_c - d
    green_old_max = green_c
    grey_old_min = grey_c
    grey_old_max = mixed_c + d

    # iterate over tags
    for tag in osm_tags:

        print(tag)
        # Filter features by tag and remove empty geometries
        targets_with_tag = targets_df.loc[targets_df["tags"] == tag]

        # Extract NDVI values for tag
        ndvi_values = get_ndvi_values(targets_with_tag, ndvi, nodata, affine)
        if len(ndvi_values) == 0:
            tag_greenness[tag] = [None, None, None, 0]
            logger.info(f"{tag}: No NDVI samples")
            print(f"{tag}: No NDVI samples")
            continue

        # Calculate green/non-green probabilities
        green, not_green, green_or_not_green = calculate_probability_mass(
            ndvi_values,
            green_old_min=green_old_min,
            green_old_max=green_old_max,
            grey_old_min=grey_old_min,
            grey_old_max=grey_old_max,
        )

        tag_greenness[tag] = [green, not_green, green_or_not_green, len(ndvi_values)]

    tag_greenness_df = pd.DataFrame(tag_greenness).T.round(2)
    tag_greenness_df.columns = ["green", "grey", "green_grey", "npixels"]
    tag_greenness_df.index.name = "tag"

    tag_greenness_df.to_csv(out_file)

    return tag_greenness_df


def green_from_osm(config):
    """
    Assigns evidence for green, non_green, either to each polygon
    :param config:
    :return:
    """

    lu_polygons_file = os.path.join(
        config["output_dir"], config["name"], f"{config['name']}_lu_polygons.shp"
    )
    green_tags_file = os.path.join(
        config["output_dir"], config["name"], f"{config['name']}_green_tags.csv"
    )
    out_file = os.path.join(
        config["output_dir"], config["name"], f"{config['name']}_lu_polygons_green.shp"
    )
    # Load lookup table for osm tags and greenness
    green_tags = pd.read_csv(green_tags_file).set_index("tag")

    # Read target parts and clip to bbox
    lu_polygons = gpd.read_file(lu_polygons_file)

    # Dataframe to store beliefs
    lu_polygons["g_osm"] = 0
    lu_polygons["n_osm"] = 0
    lu_polygons["gn_osm"] = 1
    lu_polygons["nsamples"] = 0

    for tag in green_tags.index:
        print(tag)
        green, grey, green_grey, nsamples = green_tags.loc[tag]
        lu_polygons.loc[
            (lu_polygons["tags"] == tag),
            ("g_osm", "n_osm", "gn_osm", "nsamples"),
        ] = (green, grey, green_grey, nsamples)

    lu_polygons.to_file(out_file)


if __name__ == "__main__":

    pass
