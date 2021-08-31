#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Calculates green index from OSM and Sentinal-2 data to generate green routes with openrouteservice"""

__author__ = "Christina Ludwig, GIScience Research Group, Heidelberg University"
__email__ = "christina.ludwig@uni-heidelberg.de"

import os
import datetime as dt
import sys
import traceback

from modules.download import download_features
from modules.green_osm import greenness_of_osm_tags, green_from_osm
from modules.green_ndvi import green_from_ndvi
from modules.lu_polygons import generate_landuse_polygons
from modules.ndvi import ndvi
from modules.fusion import fuse
from modules.green_index import calc_green_index

from modules.utils import (
    create_subfolder,
    check_config,
    load_config,
    init_logger,
)


if __name__ == "__main__":

    # Input parameters
    aoi_name = "test2"
    config_file = "./config/config_aois.json"
    config_file_tags = "./config/config_tags.json"
    credentials_file = "./config/google_credentials.json"

    # MODULES TO BE EXECUTED
    run_download_traffic = False
    run_download_landuse = False
    run_download_buildings = False
    run_landuse_polygons = True
    run_ndvi = False
    run_green_from_osm = False
    run_green_from_s2 = False
    run_fusion = False
    run_green_index = False

    # Preparation
    # ----------------------------

    # Load configuration parameters
    config = load_config(config_file)
    check_config(config, aoi_name)
    config_tags = load_config(config_file_tags)

    # Create output and logs folders
    out_dir_aoi = os.path.join(config["output_dir"], aoi_name)
    os.makedirs(out_dir_aoi, exist_ok=True)
    log_dir = os.path.join(out_dir_aoi, "logs")
    os.makedirs(log_dir, exist_ok=True)

    # Set up Logger
    log_file_name = os.path.join(
        log_dir, "{0}_{1}.log".format(dt.datetime.now().strftime("%Y%m%d"), aoi_name)
    )
    if os.path.exists(log_file_name):
        os.unlink(log_file_name)
    logger = init_logger("root", log_file_name)
    logger.info("Start processing")

    # Preprocessing for Greenness Model
    # -----------------------------------

    if run_download_landuse:
        try:
            landuse_feature_dir = create_subfolder(out_dir_aoi, "ohsome/landuse")
            download_features(
                bbox=config["aois"][aoi_name]["bbox"],
                timestamp=config["aois"][aoi_name]["timestamp"],
                tags=config_tags["landuse_features"],
                outdir=landuse_feature_dir,
            )
        except Exception as e:
            logger.critical(e.with_traceback())
            sys.exit(1)

    if run_download_traffic:
        try:
            traffic_feature_dir = create_subfolder(out_dir_aoi, "ohsome/traffic")
            download_features(
                bbox=config["aois"][aoi_name]["bbox"],
                timestamp=config["aois"][aoi_name]["timestamp"],
                tags=config_tags["traffic_features"],
                outdir=traffic_feature_dir,
            )
        except Exception as e:
            logger.critical(e.with_traceback())
            sys.exit(1)

    if run_download_buildings:
        try:
            building_feature_dir = create_subfolder(out_dir_aoi, "ohsome/buildings")
            download_features(
                bbox=config["aois"][aoi_name]["bbox"],
                timestamp=config["aois"][aoi_name]["timestamp"],
                tags=config_tags["building_features"],
                outdir=building_feature_dir,
            )
        except Exception as e:
            logger.critical(e.with_traceback())
            sys.exit(1)

    if run_landuse_polygons:
        try:
            feature_dir = os.path.join(out_dir_aoi, "ohsome")
            lu_polygons_dir = create_subfolder(out_dir_aoi, "polygons")
            lu_polygons = generate_landuse_polygons(
                feature_dir, epsg=config["aois"][aoi_name]["epsg"]
            )
            lu_polygons.to_file(
                os.path.join(lu_polygons_dir, "lu_polygons.geojson"), driver="GeoJSON"
            )
        except Exception as e:
            logger.critical(e)
            logger.critical(traceback.print_exc())
            sys.exit(1)

    if run_ndvi:
        try:
            credentials = load_config(credentials_file)
            ndvi(aoi_name, config, credentials)
        except Exception as err:
            logger.critical(err)
            sys.exit(3)

    if run_green_from_osm:
        try:
            greenness_of_osm_tags(aoi_name, config)
            green_from_osm(aoi_name, config)
        except Exception as err:
            logger.critical(err)
            sys.exit(3)

    if run_green_from_s2:
        try:
            green_from_ndvi(aoi_name, config, "s2")
        except Exception as err:
            logger.critical(err)
            sys.exit(3)

    if run_fusion:
        try:
            fuse(aoi_name, config)
        except Exception as err:
            logger.critical(err)
            sys.exit(3)

    if run_green_index:
        try:
            calc_green_index(aoi_name, config)
        except Exception as err:
            logger.critical(err)
            sys.exit(3)
