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
import json
from shapely.geometry import box


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


def read_ndvi(out_dir, aoi_name):
    """
    Searches and reads in ndvi file
    :return:
    """
    # Search NDVI file
    ndvi_files = glob.glob(os.path.join(out_dir, "ndvi", aoi_name + "*ndvi_s2.tif"))
    if len(ndvi_files) == 0:
        return 1, "No NDVI file found."
    elif len(ndvi_files) > 1:
        return 1, "Too many matching NDVI files found."

    # Open NDVI file
    with rio.open(ndvi_files[0]) as src:
        ndvi = src.read(1)
        affine = src.transform
        crs = src.crs
        nodata = src.nodatavals[0]

    return 0, [ndvi, affine, crs, nodata]


def greenness_of_osm_tags(aoi_name, config):
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
    out_dir_green = os.path.join(config["output_dir"], aoi_name, "green")
    green_c = config["aois"][aoi_name]["fuzzy_centers"]["s2"]["green"]
    mixed_c = config["aois"][aoi_name]["fuzzy_centers"]["s2"]["mixed"]
    grey_c = config["aois"][aoi_name]["fuzzy_centers"]["s2"]["grey"]
    d = config["aois"][aoi_name]["fuzzy_centers"]["s2"]["d"]
    targets_file = config["aois"][aoi_name]["target_geoms"]

    # Read NDVI file
    success, result = read_ndvi(out_dir_green, aoi_name)
    if success == 1:
        logger.critical(result)
        return 1
    ndvi, affine, crs, nodata = result

    # Read landuse polygons
    assert os.path.exists(targets_file), "{0} not found.".format(targets_file)
    targets_df = gpd.read_file(targets_file).to_crs({"init": str(crs)})

    # Unique land use tags
    targets_df["TARGET_type"].replace(np.nan, "None", inplace=True)
    osm_tags = targets_df.TARGET_type.unique()

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
        targets_with_tag = targets_df.loc[targets_df["TARGET_type"] == tag]

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

    tag_greenness_df.to_csv(
        os.path.join(out_dir_green, aoi_name + "_greenness_osm_tags.csv")
    )


def green_from_osm(aoi_name, config):
    """
    Assigns evidence for green, non_green, either to each polygon
    :param config:
    :return:
    """
    out_dir_green = os.path.join(config["output_dir"], aoi_name, "green")

    # Load lookup table for osm tags and greenness
    tag_doc_file = os.path.join(out_dir_green, aoi_name + "_greenness_osm_tags.csv")
    assert os.path.exists(
        tag_doc_file
    ), "OSM - Greenness Lookup file does not exist: {0}".format(
        os.path.abspath(tag_doc_file)
    )
    tag_evidence = pd.read_csv(tag_doc_file).set_index("tag")

    # Read target parts and clip to bbox
    polygon_file = config["aois"][aoi_name]["target_geoms"]
    bbox_geom = box(*config["aois"][aoi_name]["bbox"])
    targets_df = gpd.read_file(polygon_file, bbox=bbox_geom)

    # Clip to bounding box
    bbox_df = gpd.GeoDataFrame({"geometry": [bbox_geom]}, crs={"init": "epsg:4326"})
    targets_df = gpd.overlay(bbox_df, targets_df, how="intersection")

    targets_df = targets_df.loc[:, ("TARGET_ID", "TARGET_type", "geometry")]

    # Add empty columns
    empty_green_df = pd.DataFrame(
        {"green": 0, "grey": 0, "green_grey": 1.0, "nsamples": 0},
        dtype="float16",
        index=targets_df.index,
    )
    targets_df = pd.concat([targets_df, empty_green_df], axis=1, join="inner")
    targets_df["nsamples"] = targets_df["nsamples"].astype("int32")

    for tag in tag_evidence.index:
        print(tag)
        green, grey, green_grey, nsamples = tag_evidence.loc[tag]
        targets_df.loc[
            (targets_df["TARGET_type"] == tag),
            ("green", "grey", "green_grey", "nsamples"),
        ] = (green, grey, green_grey, int(nsamples))

    # Write to geojson
    # features_json = targets_df.to_json(na="drop")
    # outfile = os.path.join(out_dir_green, aoi_name + "_greenness_osm.geojson")
    # with open(outfile, "w") as src:
    #    src.write(features_json)

    # Write to csv
    outfile = os.path.join(out_dir_green, aoi_name + "_green_osm.csv")
    pd.DataFrame(
        targets_df.loc[:, ["TARGET_ID", "green", "grey", "green_grey"]]
    ).to_csv(outfile, index=False)


if __name__ == "__main__":

    pass
