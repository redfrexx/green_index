#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Calculates green index"""

__author__ = "Christina Ludwig, GIScience Research Group, Heidelberg University"
__email__ = "christina.ludwig@uni-heidelberg.de"

import geopandas as gpd
from rasterstats import zonal_stats
from shapely.geometry import box
import rasterio as rio
import pandas as pd
import os
import subprocess
import numpy as np
from modules.utils import load_config


def calc_green_index(aoi_name, config):
    """
    Calculates the green index for each highway feature
    :return:
    """

    highway_file = os.path.join(
        config["output_dir"], aoi_name, f"{aoi_name}_highways.geojson"
    )
    green_spaces_file = os.path.join(
        config["output_dir"], aoi_name, "green", "greenness.shp"
    )

    # Convert green spaces from geojson to tif
    green_spaces_file_tif_tmp = os.path.splitext(green_spaces_file)[0] + "_tmp.tif"
    cmd = [
        "gdal_rasterize",
        "-a",
        "green",
        "-of",
        "GTiff",
        "-tr",
        "5",
        "5",
        "-ot",
        "Float32",
        green_spaces_file,
        green_spaces_file_tif_tmp,
    ]
    subprocess.call(cmd)
    green_spaces_file_tif = os.path.splitext(green_spaces_file)[0] + ".tif"
    cmd = [
        "gdalwarp",
        "-t_srs",
        "EPSG:" + config["aois"][aoi_name]["epsg"],
        green_spaces_file_tif_tmp,
        green_spaces_file_tif,
    ]
    subprocess.call(cmd)
    os.unlink(green_spaces_file_tif_tmp)

    # Load and buffer highways
    bbox = box(*config["aois"][aoi_name]["bbox"])
    highways = gpd.read_file(highway_file, bbox=bbox).to_crs(
        {"init": "epsg:{0}".format(config["aois"][aoi_name]["epsg"])}
    )
    highways_buffered = highways.copy()
    highways_buffered["geometry"] = highways_buffered.buffer(30)

    # Extract green index for each highway
    res = zonal_stats(
        highways_buffered["geometry"].values, green_spaces_file_tif, stats="mean"
    )
    highways["green"] = pd.DataFrame(res)["mean"]
    del highways_buffered

    # Export to file
    highways["green"] = highways["green"] * 100
    highways["green"] = highways["green"].round(0).astype(int)
    highways.to_file(
        os.path.join(config["output_dir"], aoi_name, "highways_green_index.shp")
    )
    highways_csv = pd.DataFrame(highways.loc[:, ["@osmId", "green"]])
    highways_csv["id"] = highways_csv["@osmId"].map(lambda x: x.split("/")[1])
    highways_csv.loc[:, ["id", "green"]].to_csv(
        os.path.join(config["output_dir"], aoi_name, "highways_green_index.csv"),
        sep=",",
        decimal=".",
        header=False,
        index=False,
    )

    return 0


if __name__ == "__main__":

    aoi_name = "dresden"
    config_file = "/Users/chludwig/Development/meinGruen/code/ors_green_index/config/config_meingruen.json"
    green_spaces_file = (
        "/Users/chludwig/Data/meingruen/stadt/DD/GV2017/VegClasses_neueMethodik.tif"
    )

    config = load_config(config_file)

    with rio.open(green_spaces_file) as src:
        gv_classes = src.read(1)

        gv_binary = np.where((gv_classes > 0) & (gv_classes < 10), 1, 0)
        gv_binary = np.where(np.isnan(gv_binary), -1, gv_binary).astype("int8")

        green_spaces_file_tif = os.path.join(
            config["out_dir"], aoi_name, "gv_binary.tif"
        )
        with rio.open(
            green_spaces_file_tif,
            "w",
            driver="GTiff",
            height=gv_binary.shape[0],
            width=gv_binary.shape[1],
            count=1,
            dtype="int8",
            nodata=-1,
            crs=src.crs,
            transform=src.transform,
        ) as dst:
            dst.write(gv_binary, 1)

    del gv_classes, gv_binary

    calc_green_index(aoi_name, config, green_spaces_file_tif)

    os.unlink(green_spaces_file_tif)
