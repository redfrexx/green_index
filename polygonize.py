#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""__description__
"""

__author__ = "Christina Ludwig, GIScience Research Group, Heidelberg University"
__email__ = "christina.ludwig@uni-heidelberg.de"


import pygeos
import geopandas as gpd
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from utils import split_highways
gpd.options.use_pygeos = True
from shapely.geometry import box, LineString
import glob
import os
import pandas as pd

in_dir = "./data/dresden/raw"
out_dir = "./data/dresden/polygons"


def generate_street_blcoks():
    """
    Generates streets blocks
    :return:
    """
    files = glob.glob(os.path.join(in_dir, "*.geojson"))

    all_features = []
    for f in files:
        features = gpd.read_file(f)
        features.drop('@snapshotTimestamp', axis=1, inplace=True)
        features["type"] = os.path.basename(f).split(".")[0]
        all_features.append(features)

    features_df = pd.concat(all_features)
    features_df.plot()
    plt.show()



# Add bounding box as linestring
#features_df = features_df.append({"@osmId": "border",
#                 "geometry": LineString(box(*features_df.total_bounds).exterior)},
#                ignore_index=True)

bbox_geom = pygeos.from_shapely(box(*features_df.total_bounds))
#bbox_df = gpd.GeoDataFrame({"geometry": [bbox_geom]})
#bbox_df.to_file(os.path.join(out_dir, "bbox.geojson"), driver="GeoJSON")

features_df.geometry.map(lambda x: x.geom_type).unique()
#features_df["geometry"] = features_df.geometry.map(lambda x: x if x.geom_type != "Polygon" else LineString(x.exterior))

poly_features = features_df.loc[features_df.geometry.map(lambda x: x.geom_type == "Polygon")]
poly_geoms = poly_features.apply(lambda x: pygeos.from_shapely(x["geometry"]), axis=1)

# Line Features
line_features = features_df.loc[features_df.geometry.map(lambda x: x.geom_type != "Polygon")]
line_geoms = line_features.apply(lambda x: pygeos.from_shapely(x["geometry"]), axis=1)
line_geoms_buf = pygeos.buffer(line_geoms, 0.00001)
#line_geoms_buf_df = gpd.GeoDataFrame({"geometry": pygeos.get_parts(line_geoms_buf)})
#line_geoms_buf_df.plot(color="blue")
#plt.show()
#line_geoms_buf_df.to_file(os.path.join(out_dir, "test.geojson"), driver="GeoJSON")

# Difference
all_geoms = pd.concat([line_geoms_buf, poly_geoms], axis=0)
all_geoms_union = pygeos.union_all(all_geoms)
#all_geoms_union_df = gpd.GeoDataFrame({"geometry": [all_geoms_union]})
#all_geoms_union_df.to_file(os.path.join(out_dir, "all_poly_union.geojson"), driver="GeoJSON")


geoms_diff = pygeos.symmetric_difference(bbox_geom, all_geoms_union)
geom_diff_df = gpd.GeoDataFrame({"geometry": pygeos.get_parts(geoms_diff)})
#geom_diff_df.plot()
#plt.show()
geom_diff_df.to_file(os.path.join(out_dir, "dresden_polygons.geojson"), driver="GeoJSON")

# Split lines at intersections
#geoms = features_df.apply(lambda x: pygeos.from_shapely(x["geometry"]), axis=1)
#geoms_split = pygeos.union_all(np.array(geoms))

# Polygonize
#geoms = features.apply(lambda x: pygeos.from_shapely(x["geometry"]), axis=1)
#polygons = pygeos.polygonize(np.array(geoms))
#polygons = pygeos.polygonize(pygeos.get_parts(geoms_split))

# Convert to geopandas dataframe
#polygons_shapely = pygeos.to_shapely(polygons)
#polygons_df = gpd.GeoDataFrame({"geometry": pygeos.get_parts(geoms_split)})

#polygons_df.plot()
#plt.show()

#polygons_df.to_file(os.path.join(out_dir, "polygons.geojson"), driver="GeoJSON")

