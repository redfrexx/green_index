#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Integrate trees in greenness map"""

__author__ = "Christina Ludwig, GIScience Research Group, Heidelberg University"
__email__ = "christina.ludwig@uni-heidelberg.de"


import geopandas as gpd
import pandas as pd
import os


def integrate_trees(config):
    """
    Integrates trees from OSM into the greenness map
    :param aoi_name:
    :param config:
    :return:
    """

    aoi_name = config["name"]

    lu_polygons_file = os.path.join(
        config["output_dir"], aoi_name, f"{aoi_name}_greenness.shp"
    )
    trees_file = os.path.join(config["output_dir"], aoi_name, "osm/trees/trees.geojson")

    # Import data
    lu_polygons = gpd.read_file(lu_polygons_file)
    trees = gpd.read_file(trees_file)
    trees = trees[["geometry"]].to_crs(lu_polygons.crs)

    # Buffer trees by default size of 3 meters
    trees["geometry"] = trees.buffer(3)

    # Clip trees out of lu polygons
    lu_polygons_no_trees = gpd.overlay(lu_polygons, trees, how="symmetric_difference")

    trees["green"] = 0.90
    trees["grey"] = 0.05
    trees["green_grey"] = 0.05

    lu_polygons_trees = pd.concat([lu_polygons_no_trees, trees], axis=0)

    output_file = os.path.join(
        config["output_dir"], aoi_name, f"{aoi_name}_greenness.shp"
    )
    lu_polygons_trees.to_file(output_file)

    return 0
