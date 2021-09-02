#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Calculates green index from OSM and Sentinal-2 data to generate green routes with openrouteservice"""

__author__ = "Christina Ludwig, GIScience Research Group, Heidelberg University"
__email__ = "christina.ludwig@uni-heidelberg.de"

import os
import datetime as dt
import sys
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
    rasterize,
)


if __name__ == "__main__":

    # Input parameters
    aoi_name = "test2"
    config_file = "./config/config_aois.json"
    config_file_tags = "./config/config_tags.json"
    credentials_file = "./config/google_credentials.json"

    # Set modules to be executed to True
    run_download_traffic = False
    run_download_landuse = False
    run_download_buildings = False
    run_download_ors_highways = False
    run_download_trees = False
    run_landuse_polygons = False
    run_ndvi = False
    run_green_osm_tags = False
    run_green_osm = False
    run_green_s2 = False
    run_fusion = False
    run_rasterization = True

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
                layers=config_tags["landuse"],
                outdir=landuse_feature_dir,
            )
        except Exception:
            logger.exception("Error during download of landuse features:")
            sys.exit(1)

    if run_download_traffic:
        try:
            traffic_feature_dir = create_subfolder(out_dir_aoi, "ohsome/traffic")
            download_features(
                bbox=config["aois"][aoi_name]["bbox"],
                timestamp=config["aois"][aoi_name]["timestamp"],
                layers=config_tags["traffic"],
                outdir=traffic_feature_dir,
            )
        except Exception:
            logger.exception("Error during download of traffic features:")
            sys.exit(1)

    if run_download_buildings:
        try:
            building_feature_dir = create_subfolder(out_dir_aoi, "ohsome/buildings")
            download_features(
                bbox=config["aois"][aoi_name]["bbox"],
                timestamp=config["aois"][aoi_name]["timestamp"],
                layers=config_tags["buildings"],
                outdir=building_feature_dir,
            )
        except Exception:
            logger.exception("Error during download of buildings:")
            sys.exit(1)

    if run_download_ors_highways:
        try:
            building_feature_dir = create_subfolder(out_dir_aoi, "ohsome/highways")
            download_features(
                bbox=config["aois"][aoi_name]["bbox"],
                timestamp=config["aois"][aoi_name]["timestamp"],
                layers=config_tags["highways"],
                outdir=building_feature_dir,
            )
        except Exception:
            logger.exception("Error during download of highways:")
            sys.exit(1)

    if run_download_trees:
        try:
            building_feature_dir = create_subfolder(out_dir_aoi, "ohsome/trees")
            download_features(
                bbox=config["aois"][aoi_name]["bbox"],
                timestamp=config["aois"][aoi_name]["timestamp"],
                layers=config_tags["trees"],
                outdir=building_feature_dir,
            )
        except Exception:
            logger.exception("Error during download of buildings:")
            sys.exit(1)

    if run_landuse_polygons:
        try:
            feature_dir = os.path.join(out_dir_aoi, "ohsome")
            lu_polygons_dir = create_subfolder(out_dir_aoi, "polygons")
            lu_polygons = generate_landuse_polygons(
                in_dir=feature_dir, epsg=config["aois"][aoi_name]["epsg"]
            )
            lu_polygons.to_file(
                os.path.join(lu_polygons_dir, f"{aoi_name}_lu_polygons.shp")
            )
        except Exception:
            logger.exception("Error during generation of landuse polygons:")
            sys.exit(1)

    if run_ndvi:
        ndvi_dir = create_subfolder(out_dir_aoi, "ndvi")
        try:
            credentials = load_config(credentials_file)
            ndvi(aoi_name, config, credentials, ndvi_dir)
        except Exception:
            logger.exception("Error during generation of NDVI calculation:")
            sys.exit(1)

    if run_green_osm_tags:
        ndvi_dir = create_subfolder(out_dir_aoi, "ndvi")
        lu_polygons_file = os.path.join(out_dir_aoi, f"{aoi_name}_lu_polygons.shp")
        try:
            green_tags = greenness_of_osm_tags(
                aoi_name, config, ndvi_dir, lu_polygons_file
            )
            green_tags.to_csv(os.path.join(out_dir_aoi, aoi_name + "_green_tags.csv"))
        except Exception:
            logger.exception("Error during calculating greenness of OSM tags:")
            sys.exit(1)

    if run_green_osm:
        ndvi_dir = create_subfolder(out_dir_aoi, "ndvi")
        lu_polygons_file = os.path.join(out_dir_aoi, f"{aoi_name}_lu_polygons.shp")
        green_tags_file = os.path.join(out_dir_aoi, aoi_name + "_green_tags.csv")
        try:
            lu_polygons_osm = green_from_osm(lu_polygons_file, green_tags_file)
            lu_polygons_osm.to_file(lu_polygons_file)
        except Exception:
            logger.exception("Error during greenness calculation from OSM data:")
            sys.exit(1)

    if run_green_s2:
        lu_polygons_file = os.path.join(out_dir_aoi, f"{aoi_name}_lu_polygons.shp")
        ndvi_dir = create_subfolder(out_dir_aoi, "ndvi")
        try:
            lu_polygons_ndvi = green_from_ndvi(
                aoi_name, config, ndvi_dir, lu_polygons_file
            )
            lu_polygons_ndvi.to_file(lu_polygons_file)
        except Exception:
            logger.exception("Error during greenness calculation from Sentinel-2 data:")
            sys.exit(3)

    if run_fusion:
        beliefs_dir = create_subfolder(out_dir_aoi, "beliefs")
        lu_polygons_file = os.path.join(out_dir_aoi, f"{aoi_name}_lu_polygons.shp")
        lu_polygons_file2 = os.path.join(
            out_dir_aoi, f"{aoi_name}_lu_polygons_final.shp"
        )
        try:
            class_ndvi_osm_geo = fuse(aoi_name, config, lu_polygons_file)
            class_ndvi_osm_geo.to_file(lu_polygons_file2)
        except Exception:
            logger.exception("Error during greenness fusion:")
            sys.exit(3)

    if run_rasterization:
        try:
            lu_polygons_file = os.path.join(
                out_dir_aoi, f"{aoi_name}_lu_polygons_final.shp"
            )
            rasterize(lu_polygons_file, column="green")
        except Exception:
            logger.exception("Error during green index calculation:")
            sys.exit(3)
