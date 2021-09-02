#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Calculates green index from OSM and Sentinal-2 data to generate green routes with openrouteservice"""

__author__ = "Christina Ludwig, GIScience Research Group, Heidelberg University"
__email__ = "christina.ludwig@uni-heidelberg.de"

import os
import sys
import argparse
from modules.download import download_features
from modules.green_osm import greenness_of_osm_tags, green_from_osm
from modules.green_ndvi import green_from_ndvi
from modules.lu_polygons import generate_landuse_polygons
from modules.ndvi import ndvi
from modules.fusion import fuse

from modules.utils import (
    create_subfolder,
    check_config,
    load_config,
    init_logger,
    rasterize,
)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Calculates the greenness based on OSM and Sentinel-2 data."
    )
    parser.add_argument(
        "--aoi",
        "-a",
        required=True,
        dest="config_file",
        type=str,
        help="Path to configuration file (.json) describing the area of interest (aoi)",
    )
    parser.add_argument(
        "--google_cred",
        "-g",
        required=True,
        dest="google_cred_file",
        type=str,
        help="Path to json file containing credentials for using google earth engine",
    )
    args = parser.parse_args()

    # Set modules to be executed to True
    run_download_traffic = False
    run_download_landuse = False
    run_download_buildings = False
    run_download_trees = False
    run_landuse_polygons = False
    run_ndvi = False
    run_green_osm_tags = False
    run_green_osm = False
    run_green_ndvi = False
    run_fusion = True
    run_rasterization = True

    # Load configuration parameters
    config = load_config(args.config_file)
    check_config(config)

    # Create output folders
    aoi_name = config["name"]
    out_dir_aoi = os.path.join(config["output_dir"], aoi_name)
    os.makedirs(out_dir_aoi, exist_ok=True)

    log_file = os.path.join(out_dir_aoi, "log.log")
    logger = init_logger("calculate green", log_file)

    config_tags = load_config("./config/config_tags.json")

    if run_download_landuse:
        logger.info("Download land use features...")
        try:
            landuse_feature_dir = create_subfolder(out_dir_aoi, "osm/landuse")
            download_features(
                bbox=config["bbox"],
                timestamp=config["timestamp"],
                layers=config_tags["landuse"],
                outdir=landuse_feature_dir,
            )
        except Exception:
            logger.exception("Error during download of landuse features:")
            sys.exit(1)

    if run_download_traffic:
        logger.info("Download traffic features...")
        try:
            traffic_feature_dir = create_subfolder(out_dir_aoi, "osm/traffic")
            download_features(
                bbox=config["bbox"],
                timestamp=config["timestamp"],
                layers=config_tags["traffic"],
                outdir=traffic_feature_dir,
            )
        except Exception:
            logger.exception("Error during download of traffic features:")
            sys.exit(1)

    if run_download_buildings:
        logger.info("Download building features...")
        try:
            building_feature_dir = create_subfolder(out_dir_aoi, "osm/buildings")
            download_features(
                bbox=config["bbox"],
                timestamp=config["timestamp"],
                layers=config_tags["buildings"],
                outdir=building_feature_dir,
            )
        except Exception:
            logger.exception("Error during download of buildings:")
            sys.exit(1)

    if run_download_trees:
        logger.info("Download tree features...")
        try:
            trees_dir = create_subfolder(out_dir_aoi, "osm/trees")
            download_features(
                bbox=config["bbox"],
                timestamp=config["timestamp"],
                layers=config_tags["trees"],
                outdir=trees_dir,
            )
        except Exception:
            logger.exception("Error during download of trees:")
            sys.exit(1)

    if run_landuse_polygons:
        logger.info("Generate land use polygons ...")
        try:
            generate_landuse_polygons(config)
        except Exception:
            logger.exception("Error during generation of landuse polygons:")
            sys.exit(1)

    if run_ndvi:
        logger.info("Calculate NDVI ...")
        try:
            credentials = load_config(args.google_cred_file)
            ndvi(config, credentials)
        except Exception:
            logger.exception("Error during generation of NDVI calculation:")
            sys.exit(1)

    if run_green_osm_tags:
        logger.info("Calculate greenness of OSM tags...")
        try:
            greenness_of_osm_tags(config)
        except Exception:
            logger.exception("Error during calculating greenness of OSM tags:")
            sys.exit(1)

    if run_green_osm:
        logger.info("Calculate greenness of polygons from OSM data...")
        try:
            green_from_osm(config)
        except Exception:
            logger.exception("Error during greenness calculation from OSM data:")
            sys.exit(1)

    if run_green_ndvi:
        logger.info("Calculate greenness of polygons from NDVI...")
        try:
            green_from_ndvi(config)
        except Exception:
            logger.exception("Error during greenness calculation from Sentinel-2 data:")
            sys.exit(3)

    if run_fusion:
        logger.info("Fuse greenness from OSM and NDVI...")
        try:
            class_ndvi_osm_geo = fuse(config)
        except Exception:
            logger.exception("Error during greenness fusion:")
            sys.exit(3)

    if run_rasterization:
        logger.info("Rasterize greenness of polygons...")
        try:
            lu_polygons_file = os.path.join(
                out_dir_aoi, f"{aoi_name}_lu_polygons_final.shp"
            )
            rasterize(lu_polygons_file, column="green")
        except Exception:
            logger.exception("Error during green index calculation:")
            sys.exit(3)
