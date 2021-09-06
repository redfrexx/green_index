#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""__description__"""

__author__ = "Christina Ludwig, GIScience Research Group, Heidelberg University"
__email__ = "christina.ludwig@uni-heidelberg.de"

import logging

import pygeos
import geopandas as gpd
from modules.utils import calc_compactness
from shapely.geometry import box
import glob
import os
import pandas as pd
import numpy as np

gpd.options.use_pygeos = True


def generate_street_blocks_gpd(in_dir):
    """
    Generate street blcoks using geopandas (this is 8 times slower than generate_street_blocks())
    :param in_dir:
    :param out_dir:
    :return:
    """
    # Load traffic features
    files = glob.glob(os.path.join(in_dir, "*.geojson"))
    all_features = []
    for f in files:
        features = gpd.read_file(f)
        features.drop("@snapshotTimestamp", axis=1, inplace=True)
        features["type"] = os.path.basename(f).split(".")[0]
        all_features.append(features)
    features_df = pd.concat(all_features)

    # Buffer line features
    line_features = features_df.loc[
        features_df.geometry.map(lambda x: x.geom_type != "Polygon")
    ]
    line_features.loc[:, "geometry"] = line_features.buffer(0.00001)

    poly_features = features_df.loc[
        features_df.geometry.map(lambda x: x.geom_type == "Polygon")
    ]
    all_features = pd.concat([line_features, poly_features])

    # Bounding box as polygon
    bbox_geom = box(*features_df.total_bounds)
    bbox_df = gpd.GeoDataFrame({"geometry": [bbox_geom]}, crs="epsg:4326")

    # Calculate symmetric difference
    geoms_diff = gpd.overlay(all_features, bbox_df, how="symmetric_difference")

    return geoms_diff


def polygons_from_traffic(in_dir):
    """
    Generates street blocks by intersecting traffic features with each other. Sliver polygons with width < 10m
    are merged with neighbouring polygons
    :return:
    """
    # Load traffic features
    files = glob.glob(os.path.join(in_dir, "traffic", "*.geojson"))
    assert len(files) > 0, f"No OSM features not found in {in_dir}"
    all_features = []
    for f in files:
        features = gpd.read_file(f, names=["geometry"])
        features = features.loc[:, ["geometry"]]
        features["type"] = os.path.basename(f).split(".")[0]
        all_features.append(features)
    features_df = pd.concat(all_features)

    # Bounding box as polygon
    bbox_geom = pygeos.from_shapely(box(*features_df.total_bounds))

    # Buffer line Features
    # todo: adjust buffers based on traffic feature type
    line_features = features_df.loc[
        features_df.geometry.map(
            lambda x: x.geom_type not in ("Polygon", "MultiPolygon")
        )
    ]
    line_geoms = line_features.apply(
        lambda x: pygeos.from_shapely(x["geometry"]), axis=1
    )
    line_geoms_buf = pygeos.buffer(line_geoms, 0.00005)

    # Merge buffered line features with polygon features
    poly_features = features_df.loc[
        features_df.geometry.map(lambda x: x.geom_type in ("Polygon", "MultiPolygon"))
    ]
    poly_geoms = poly_features.apply(
        lambda x: pygeos.from_shapely(x["geometry"]), axis=1
    )
    all_geoms = np.concatenate(
        (np.array(poly_geoms).ravel(), np.array(line_geoms_buf).ravel())
    )
    all_geoms_union = pygeos.union_all(all_geoms)

    # Calculate symmetric difference
    geoms_diff = pygeos.symmetric_difference(bbox_geom, all_geoms_union)
    geom_diff_df = gpd.GeoDataFrame(
        {"geometry": pygeos.get_parts(geoms_diff)}, crs="epsg:4326"
    )
    geom_diff_df = geom_diff_df.loc[
        geom_diff_df.geometry.map(lambda x: x.geom_type in ("Polygon", "Multipolygon"))
    ]
    geom_diff_df = geom_diff_df.explode()

    return geom_diff_df


def polygons_from_landuse(in_dir, street_blocks, epsg):
    """
    Generates land use polygons
    Features contained in bigger feature: Clip small feature from big feature
    Features overlapping two features: Split features in three parts
    :return
    """
    street_blocks["tags"] = np.empty((len(street_blocks), 0)).tolist()
    street_blocks = street_blocks.to_crs(epsg=epsg)
    lu_polygons = street_blocks.copy()

    # Load traffic features
    files = glob.glob(os.path.join(in_dir, "landuse", "*.geojson"))
    for f in files:
        key = os.path.basename(f).split(".")[0]
        features = gpd.read_file(f)
        if len(features) == 0:
            continue
        features = features.loc[:, ["geometry", key]]
        features = features.to_crs(epsg=epsg)

        values = features[key].unique()
        for val in values:
            selected_features = features.loc[features[key] == val]
            # Land use blocks which intersect selected land use features
            intersected = gpd.overlay(
                lu_polygons,
                selected_features,
                how="intersection",
                keep_geom_type=True,
            )
            intersected["tags"] = intersected["tags"].map(lambda x: x + [(key, val)])

            # Land use blocks which don't intersect the selected features
            difference = gpd.overlay(lu_polygons, selected_features, how="difference")

            lu_polygons = pd.concat(
                [intersected[["geometry", "tags"]], difference[["geometry", "tags"]]],
                axis=0,
                ignore_index=True,
            )
    lu_polygons["tags"] = lu_polygons["tags"].map(
        lambda tags: ";".join([f"{t[0]}={t[1]}" for t in tags])
    )
    lu_polygons = lu_polygons.loc[
        lu_polygons.geometry.map(lambda x: x.geom_type in ("Polygon", "MultiPolygon"))
    ]
    lu_polygons = lu_polygons.explode()

    return lu_polygons


def clean_polygons(lu_polygons):
    """
    Cleans land use polygons by removing small sliver polygons
    :return:
    """
    lu_polygons["geometry"] = lu_polygons.buffer(-3, join_style=2, resolution=2).buffer(
        3, join_style=2, resolution=2
    )
    lu_polygons.reset_index(inplace=True, drop=True)
    lu_polygons = lu_polygons.explode()
    lu_polygons["compactness"] = calc_compactness(lu_polygons)
    lu_polygons["area"] = lu_polygons.area
    lu_polygons = lu_polygons.loc[
        ~((lu_polygons.area < 2500) & (lu_polygons.compactness < 0.05))
    ]
    return lu_polygons.drop("compactness", axis=1)


def clip_buildings(in_dir, lu_polygons):
    """
    Clip buildings out of land use polygons. convert all geometries to polygons. Point and line geometries are dropped.
    :param aoi_name:
    :param config:
    :return:
    """
    buildings_file = os.path.join(in_dir, "buildings", "building.geojson")
    buildings = gpd.read_file(buildings_file).to_crs(lu_polygons.crs)

    # Clip out buildings from land use polygons
    lu_polygons_clip = gpd.overlay(lu_polygons, buildings, how="difference")
    lu_polygons_clip.reset_index(drop=True, inplace=True)

    # Convert geometry collections and multipolygons to polygons
    geom_collections = lu_polygons_clip.loc[
        lu_polygons_clip.geometry.map(
            lambda x: x.geom_type in ("GeometryCollection", "MultiPolygon")
        )
    ]
    geom_collections = geom_collections.explode()
    geom_collections.reset_index(drop=True, inplace=True)
    geom_collections = geom_collections.loc[
        geom_collections.geometry.map(lambda x: x.geom_type == "Polygon")
    ]
    # Select all polygon features
    polygons = lu_polygons_clip.loc[
        lu_polygons_clip.geometry.map(lambda x: x.geom_type == "Polygon")
    ]
    # Merge all features
    lu_polygons_cleaned = pd.concat([polygons, geom_collections], axis=0)
    return lu_polygons_cleaned.loc[~lu_polygons_cleaned.is_empty]


def generate_landuse_polygons(config):
    """
    Generates landuse polygons based on traffic network and landuse features from OSM
    :return:
    """

    osm_dir = os.path.join(config["output_dir"], config["name"], "osm")
    out_file = os.path.join(
        config["output_dir"], config["name"], f"{config['name']}_lu_polygons.shp"
    )

    street_blocks = polygons_from_traffic(osm_dir)
    # street_blocks.to_file(
    #    os.path.join(
    #        config["output_dir"], config["name"], f"{config['name']}_street_blocks.shp"
    #    )
    # )

    lu_polygons = polygons_from_landuse(osm_dir, street_blocks, config["epsg"])
    # lu_polygons.to_file(
    #    os.path.join(
    #        config["output_dir"],
    #        config["name"],
    #        f"{config['name']}_lu_polygons_raw.shp",
    #    )
    # )

    lu_polygons_clean = clean_polygons(lu_polygons)
    # lu_polygons_clean.to_file(
    #    os.path.join(
    #        config["output_dir"],
    #        config["name"],
    #        f"{config['name']}_lu_polygons_clean.shp",
    #    )
    # )

    lu_polygons_no_building = clip_buildings(osm_dir, lu_polygons_clean)
    lu_polygons_no_building = lu_polygons_no_building.loc[
        lu_polygons_no_building.geometry.map(
            lambda x: x.geom_type in ("Polygon", "MultiPolygon")
        )
    ]
    lu_polygons_no_building.reset_index(drop=True)
    lu_polygons_no_building = lu_polygons_no_building.explode()
    lu_polygons_no_building.reset_index(drop=True, inplace=True)
    lu_polygons_no_building.to_file(out_file)
