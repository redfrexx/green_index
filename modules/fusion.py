#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Fuse green information from OSM and Sentinel-2"""

__author__ = "Christina Ludwig, GIScience Research Group, Heidelberg University"
__email__ = "christina.ludwig@uni-heidelberg.de"

import os
import geopandas as gpd
import pandas as pd
from pyds import MassFunction


def convert_to_mass(data: pd.DataFrame, adapt_masses: bool = False):
    """
    Converts data frame to series of probability masses
    :param data: Dataframe containing mass values for g_xxx (green), n_xxx (not green), gn_xxx (green / not green)
    :param adapt_masses: If true, adapt masses so the maximum mass value is 0.95
    :return: pd.DataFrame containing MassFunction objects
    """
    label_dict = {x: x.split("_")[0] for x in data.columns}
    return data.rename(columns=label_dict).apply(
        lambda row: as_MassFunction(row, adapt_masses=adapt_masses), axis=1
    )


def as_MassFunction(row, adapt_masses=True):
    """
    Adaptes the masses to a maximum value of 0.95. Data manipulation is performed inplace.
    :param row: pd.Series with fields g (green), n (not green) and gn (green / not green)
    :param adapt_masses: If true, adapt masses so the maximum mass value is 0.95
    :return:
    """
    row_cp = row.copy()
    if adapt_masses:
        if row_cp["g"] > 0.95:
            row_cp["gn"] = row_cp["gn"] + row_cp["g"] - 0.95
            row_cp["g"] = 0.95
        elif row_cp["n"] > 0.95:
            row_cp["gn"] = row_cp["gn"] + row_cp["n"] - 0.95
            row_cp["n"] = 0.95
    return MassFunction(dict(row_cp))


def classify_row(row: MassFunction):
    """
    Classifies beliefs into green, grey or uncertain based on pignistic probability
    :param row: pd.Series with fields g (green), n (not green) and gn (green / not green)
    :returns
    """
    prob = row.pignistic()
    if prob["g"] > prob["n"]:
        label = "green"
    elif prob["g"] < prob["n"]:
        label = "grey"
    else:
        label = "uncertain"
    return prob["g"], prob["n"], label


def classify_mass(mass):
    """
    Derives a classification based on probability masses
    :returns Dataframe with columns green, grey, green_grey,
    """
    mass_df = pd.DataFrame({"mass": mass})
    classified = mass_df.apply(
        lambda x: classify_row(x["mass"]), axis=1, result_type="expand"
    )
    classified.columns = [
        "g_prob",
        "n_prob",
        "class",
    ]

    mass_values = mass_df.apply(lambda r: dict(r["mass"]), axis=1, result_type="expand")
    mass_values.columns = ["green", "grey", "green_grey"]

    return mass_values.join(classified)


def fuse_masses(mass1, mass2):
    """
    Fuses two probability masses
    :returns
    """
    # Convert to masses
    both_masses = pd.DataFrame({"m1": mass1, "m2": mass2})

    # Fuse masses
    both_masses["fused"] = both_masses.apply(
        lambda row: row["m1"].combine_conjunctive(row["m2"], importance_sampling=False),
        axis=1,
    )

    return both_masses["fused"]


def fuse(config):
    """
    Fuses information from OSM and Sentinel-2 on greenness
    :param aoi_name:
    :param config:
    :return:
    """

    aoi_name = config["name"]
    lu_polygons_file = os.path.join(
        config["output_dir"], aoi_name, f"{aoi_name}_lu_polygons_green.shp"
    )

    lu_polygons = gpd.read_file(lu_polygons_file)
    lu_polygons = lu_polygons
    lu_polygons["area"] = lu_polygons.area

    ndvi_beliefs = lu_polygons[["g_ndvi", "n_ndvi", "gn_ndvi"]]
    osm_beliefs = lu_polygons[["g_osm", "n_osm", "gn_osm"]]

    # Convert data to mass functions
    ndvi_mass = convert_to_mass(ndvi_beliefs)
    osm_mass = convert_to_mass(osm_beliefs)

    ndvi_osm_mass = fuse_masses(ndvi_mass, osm_mass)

    # Classify based on beliefs
    class_ndvi_osm = classify_mass(ndvi_osm_mass)

    class_ndvi_osm_geo = gpd.GeoDataFrame(
        class_ndvi_osm.join(lu_polygons.loc[:, "geometry"])
    )
    class_ndvi_osm_geo.crs = lu_polygons.crs
    class_ndvi_osm_geo = class_ndvi_osm_geo[["green", "grey", "green_grey", "geometry"]]
    class_ndvi_osm_geo.to_file(lu_polygons_file)

    return class_ndvi_osm_geo
