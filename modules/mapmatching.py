#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Match routes to highway network"""

__author__ = "Christina Ludwig, GIScience Research Group, Heidelberg University"
__email__ = "christina.ludwig@uni-heidelberg.de"

import logging
import os

import fiona
from shapely.geometry import LineString, MultiLineString, Point, box
import geopandas as gpd
import pygeos
import numpy as np


def polygon2line(geom):
    """
    Converts Polygon and Multipolygon objects to LineStrings. For Multipolyongs only the exterior ring is used.
    :param geom: shapely.geometry object
    :return:
    """
    if geom.geom_type in ("Polygon", "MultiPolygon"):
        return LineString(geom.exterior)
    else:
        return geom


def mapmatching(config):
    """
    Match routes to highway network
    :param config:
    :return:
    """

    logger = logging.getLogger("simulate routes")

    out_dir = config["out_dir"]
    steep_level = config["steep_level"]
    profile = config["profile"]
    n_routes = config["n_routes"]
    epsg = config["epsg"]
    aoi_name = config["name"]

    # Create output directories -------------------------
    out_dir_job = os.path.join(
        out_dir, aoi_name, f"{aoi_name}_{profile}_{steep_level}_{n_routes}"
    )
    gpkg_file = os.path.join(out_dir_job, "routes.gpkg")
    outfile = os.path.join(out_dir_job, "short_green_comparison.shp")

    # Highways file
    logger.info("Preparing highway data...")
    splitted_highways_file = os.path.join(out_dir_job, "splitted_highways.shp")
    if os.path.exists(splitted_highways_file):
        highways_splitted = gpd.read_file(splitted_highways_file)
    else:
        highways_file = os.path.join(
            config["out_dir"], aoi_name, "green_index", "highways.geojson"
        )
        highways = gpd.read_file(highways_file)
        highways_splitted = split_highways(highways)
        highways_splitted = highways_splitted.to_crs(epsg=epsg)
        highways_splitted = highways_splitted.explode()
        highways_splitted["geometry"] = highways_splitted.geometry.map(
            lambda x: polygon2line(x)
        )
        highways_splitted.to_file(splitted_highways_file)

    highways_buffered = highways_splitted.copy()
    highways_buffered.geometry = highways_buffered.geometry.map(lambda x: x.buffer(0.1))
    highways_buffered.reset_index(inplace=True, drop=True)

    for layer in fiona.listlayers(gpkg_file):
        logger.info(f"Matching {layer} routes...")
        routes = gpd.read_file(gpkg_file, layer=layer)
        routes = routes.to_crs(epsg=epsg)
        highways_splitted[layer] = 0
        for r in routes.iterrows():
            matched_highway_ids = map_matching(
                r[1].geometry, highways_buffered, crs=routes.crs
            )
            highways_splitted.loc[matched_highway_ids.index, layer] += 1

    highways_splitted.to_file(outfile)


def split_highways(highways):
    """
    Split highway features at intersections
    :param highways:
    :return:
    """
    geoms = highways.apply(lambda x: pygeos.from_shapely(x["geometry"]), axis=1)
    splitted_geoms = pygeos.union_all(np.array(geoms))
    return gpd.GeoDataFrame(
        {"geometry": pygeos.get_parts(splitted_geoms)}, crs="epsg:4326"
    )


def map_matching(route_geometry, highways_buffered, crs):
    """
    Matches route to highways
    :param route:
    :param highways_buffered:
    :return:
    """
    coords = convert_line_to_points(route_geometry, 5)  # extract_coords(route_geometry)
    coords = filter_duplicate_points(coords)
    point_df = gpd.GeoDataFrame({"geometry": [Point(c) for c in coords]}, crs=crs)
    matching_highways = gpd.tools.sjoin(
        point_df, highways_buffered, how="left", op="within"
    )

    highway_counts = matching_highways.groupby("index_right").count()
    if len(highway_counts) == 0:
        return highway_counts
    else:
        highway_counts = highway_counts.loc[highway_counts["geometry"] >= 2]
    return highway_counts


def filter_duplicate_points(coords):
    """
    Filter out duplicate points
    :param coords:
    :return:
    """
    coords_t = []
    [coords_t.extend(list(c.coords)) for c in coords]
    coords_t = set(coords_t)
    return [Point(c) for c in coords_t]


def convert_line_to_points(geom, dist):
    """
    Converts a line string into a list of points with distance
    :param dist:
    :return:
    """
    geom = remove_third_dimension(geom)
    if isinstance(geom, LineString):
        geom = [geom]

    points = []
    for ls in geom:
        dist_all = dist
        if len(points) == 0 or (not points[-1].equals(Point(ls.coords[0]))):
            points.append(Point(ls.coords[0]))
        points.append(ls.interpolate(dist_all))
        while (ls.distance(points[-1]) < 1e-8) and (
            not points[-1].equals(Point(ls.coords[-1]))
        ):
            dist_all += dist
            points.append(ls.interpolate(dist_all))
    return points


def remove_third_dimension(geom):
    """
    Remove z coordinate of geometries
    :param geom:
    :return:
    """
    if geom.is_empty:
        return geom

    elif isinstance(geom, LineString):
        return LineString([xy[0:2] for xy in list(geom.coords)])

    elif isinstance(geom, MultiLineString):
        lines = list(geom.geoms)
        new_lines = []
        for line in lines:
            new_lines.append(remove_third_dimension(line))

        return MultiLineString(new_lines)

    else:
        print("not supported")


def extract_coords(geom):
    """
    Returns the coordinates of a line string or MultiLineString
    :return:
    """
    if isinstance(geom, LineString):
        return list(set(geom.coords))
    elif isinstance(geom, MultiLineString):
        result = []
        [result.extend(list(g.coords)) for g in geom.geoms]
        return list(set(result))
    else:
        raise ValueError(f"{type(geom)} not supported.")
