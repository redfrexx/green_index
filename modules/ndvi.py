#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Calculation of NDVI using Earth Engine API"""

__author__ = "Christina Ludwig, GIScience Research Group, Heidelberg University"
__email__ = "christina.ludwig@uni-heidelberg.de"


import ee
import requests
import os
import zipfile
import re
import logging
import json
import subprocess

from modules.utils import split_bbox, load_config, create_subfolder


def get_filename_from_cd(cd):
    """
    Get filename from content-disposition
    :returns
    """
    if not cd:
        return None
    fname = re.findall("filename=(.+)", cd)
    if len(fname) == 0:
        return None
    return fname[0]


def create_ee_polygon(bbox):
    """
    Create a polygon geometry from corner coordinates
    :param bbox: list containing xmin, ymin, xmax, ymax
    :return:
    """
    xmin, ymin, xmax, ymax = bbox
    coords = [[xmin, ymax], [xmin, ymin], [xmax, ymin], [xmax, ymax], [xmin, ymax]]
    return coords, ee.Algorithms.GeometryConstructors.Polygon(
        [[xmin, ymax], [xmin, ymin], [xmax, ymin], [xmax, ymax], [xmin, ymax]],
        "EPSG:4326",
    )


def ndvi(config, google_credentials):
    """
    Calculates a NDVI max composite using Google Earth Engine

    :param aoi_name:
    :param config:
    :return:
    """

    logger = logging.getLogger("root." + __name__)
    logger.info("Calculating NDVI ...")

    # Read config file
    aoi_name = config["name"]
    bbox = config["bbox"]
    year = config["ndvi_year"]
    cloudcov = config["cloud_coverage"]
    start_date = f"{year}-01-01"
    end_date = f"{year}-12-31"
    target_crs = "epsg:{0}".format(config["epsg"])
    output_dir = create_subfolder(config["output_dir"], f"{aoi_name}/ndvi")

    # bbox in tiles
    bbox_tiles = split_bbox(bbox, (0.25, 0.25))
    ndvi_file_names = []
    for i, bbox in enumerate(bbox_tiles):
        print(i)

        # Create polygon for AOI
        credentials = ee.ServiceAccountCredentials(
            google_credentials["service_account"],
            google_credentials["service_account_json"],
        )
        ee.Initialize(credentials=credentials)
        coords, region = create_ee_polygon(bbox)

        # Find suitable Sentinel-2 scenes
        collection = (
            ee.ImageCollection("COPERNICUS/S2")
            .filterDate(start_date, end_date)
            .filterMetadata("CLOUDY_PIXEL_PERCENTAGE", "less_than", cloudcov)
            .filterBounds(region)
        )
        images = collection.toList(500)
        logger.info("Number of Sentinel-2 scenes: {0}".format(len(images.getInfo())))

        # Write metadata of scenes to file
        scene_info_file = os.path.join(
            output_dir, "{0}_{1}_scenes.json".format(aoi_name, i)
        )
        with open(scene_info_file, "w") as dst:
            info = json.dumps(images.getInfo(), indent=4)
            dst.write(info)

        # Calculate NDVI max composite --------------------------------------

        def addNDVI(image):
            ndvi = image.normalizedDifference(["B8", "B4"]).rename("NDVI")
            return image.addBands(ndvi)

        ndvi4scenes = collection.map(addNDVI)
        max_ndvi = ndvi4scenes.qualityMosaic("NDVI").select("NDVI")

        # Download NDVI ----------------------------------------------------
        out = max_ndvi.getDownloadUrl(
            {
                "name": "_".join([aoi_name, start_date, end_date]),
                "scale": 10,
                "crs": target_crs,
                "region": coords,
                "fileFormat": "GeoTIFF",
                "maxPixels": 1e12,
            }
        )
        r = requests.get(out, allow_redirects=True)
        zip_file_name = get_filename_from_cd(r.headers.get("content-disposition"))
        zip_file_path = os.path.join(output_dir, zip_file_name)
        open(zip_file_path, "wb").write(r.content)

        # Unzip
        with zipfile.ZipFile(zip_file_path, "r") as zip_ref:
            zip_ref.extractall(output_dir)

        # Delete zip file
        if os.path.exists(zip_file_path):
            os.unlink(zip_file_path)
        tfw_file = zip_file_path[:-4] + ".NDVI.tfw"
        if os.path.exists(tfw_file):
            os.unlink(tfw_file)

        # Rename ndvi file
        ndvi_file = os.path.splitext(zip_file_path)[0] + ".NDVI.tif"
        ndvi_file_new = os.path.join(output_dir, f"{aoi_name}_ndvi_{i}.tif")
        os.rename(ndvi_file, ndvi_file_new)
        ndvi_file_names.append(ndvi_file_new)

    # merge ndvi files
    ndvi_file = os.path.join(output_dir, f"{aoi_name}_ndvi.tif")
    cmd = ["gdal_merge.py", "-o", ndvi_file] + ndvi_file_names
    subprocess.call(cmd)
    [os.unlink(f) for f in ndvi_file_names]
    logger.info("Calculating NDVI - done.")


if __name__ == "__main__":
    aoi = "saopaulo"
    config_path = "../config/config_saopaulo.json"
    config = load_config(config_path)
    ndvi(aoi, config)
