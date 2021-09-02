#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Calculates green index"""

__author__ = "Christina Ludwig, GIScience Research Group, Heidelberg University"
__email__ = "christina.ludwig@uni-heidelberg.de"

import fiona
import geopandas as gpd
from rasterstats import zonal_stats
import rasterio as rio
import pandas as pd
import os


def calc_green_index(width, output_dir, vector_file=None, raster_file=None):
    """
    Calculates the green index for each highway feature
    :return:
    """
    green_index_file = os.path.join(output_dir, "green_index.geojson")
    green_index_csv_file = os.path.join(output_dir, "green_index.csv")
    highway_file = os.path.join(output_dir, "highways.geojson")

    # Check if coordinate system is metric
    if raster_file:
        with rio.open(raster_file) as src:
            crs = src.crs
        assert (
            crs.is_projected
        ), f"Coordinate references system of raster file is not projected (epsg:{crs.epsg})."
    if vector_file:
        with fiona.open(vector_file) as src:
            crs = src.crs
        assert (
            crs.is_projected
        ), f"Coordinate references system of vector file is not projected (epsg:{crs.epsg})."

    # Buffer highways
    highways = gpd.read_file(highway_file).loc[:, ["@osmId", "geometry"]]
    highways_buffered = highways.copy().to_crs(crs)
    highways_buffered["geometry"] = highways_buffered.buffer(width)

    # Extract green index for each highway
    if raster_file:
        res = zonal_stats(
            highways_buffered["geometry"].values, raster_file, stats="mean"
        )
        highways["index"] = pd.DataFrame(res)["mean"]
        del highways_buffered
        highways["index"] = highways["index"] * 100
        highways["index"] = highways["index"].round(0).astype(int)

    if vector_file:
        raise NotImplementedError("Processing of vector file not implemented yet.")

    # Export to file
    highways["id"] = highways["@osmId"].map(lambda x: x.split("/")[1])
    highways.to_file(green_index_file, driver="GeoJSON")

    # Create csv for openrouteservice
    # highways_csv = pd.DataFrame(highways.loc[:, ["id", "index"]])
    highways.loc[:, ["id", "index"]].to_csv(
        green_index_csv_file,
        sep=",",
        decimal=".",
        header=False,
        index=False,
    )
