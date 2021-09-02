#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Fuse green information from OSM and Sentinel-2"""

__author__ = "Christina Ludwig, GIScience Research Group, Heidelberg University"
__email__ = "christina.ludwig@uni-heidelberg.de"

import os
import geopandas as gpd
import numpy as np

from pyds import MassFunction
import pandas as pd


def convert_to_mass(data, adapt_masses=False):
    """
    Converts data frame to series of probability masses
    :return
    """
    label_dict = {x: x.split("_")[0] for x in data.columns}
    return data.rename(columns=label_dict).apply(
        lambda row: as_MassFunction(row, adapt_masses=adapt_masses), axis=1
    )


def as_MassFunction(row, adapt_masses=True):
    """
    Adaptes the masses to a maximum value of 0.95. Data manipulation is performed inplace.
    :param row:
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


def classify_green(row):
    """
    Classifies beliefs into green, grey or uncertain based on majority vote
    :returns
    """
    uc_green = row.pl("g") - row.bel("g")
    uc_grey = row.pl("n") - row.bel("n")
    prob = row.pignistic()
    if prob["g"] > prob["n"]:
        label = "green"
    elif prob["g"] < prob["n"]:
        label = "grey"
    else:
        label = "uncertain"
    return uc_green, uc_grey, prob["g"], prob["n"], label


def mass_to_class(mass):
    """
    Dervies a classification based on probability masses
    :returns
    """
    mass_df = pd.DataFrame({"mass": mass})
    classified = mass_df.apply(
        lambda x: classify_green(x["mass"]), axis=1, result_type="expand"
    )
    classified.columns = [
        "g_unc",
        "n_unc",
        "g_prob",
        "n_prob",
        "class",
    ]

    fused_mass = mass_df.apply(lambda r: dict(r["mass"]), axis=1, result_type="expand")
    fused_mass.columns = ["green", "grey", "green_grey"]

    return fused_mass.join(classified)


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
        config["output_dir"], aoi_name, f"{aoi_name}_lu_polygons.shp"
    )
    lu_polygons_file_fused = os.path.join(
        config["output_dir"], aoi_name, f"{aoi_name}_lu_polygons_final.shp"
    )

    lu_polygons = gpd.read_file(lu_polygons_file)
    lu_polygons = lu_polygons
    lu_polygons["area"] = lu_polygons.area
    # lu_polygons["tags"].replace(np.nan, "None", inplace=True)

    ndvi_beliefs = lu_polygons[["g_ndvi", "n_ndvi", "gn_ndvi"]]
    osm_beliefs = lu_polygons[["g_osm", "n_osm", "gn_osm"]]

    # Convert data to mass functions

    ndvi_mass = convert_to_mass(ndvi_beliefs)
    osm_mass = convert_to_mass(osm_beliefs)

    ndvi_osm_mass = fuse_masses(ndvi_mass, osm_mass)

    # Classify based on beliefs
    # class_ndvi = mass_to_class(ndvi_mass)
    # class_osm = mass_to_class(osm_mass)
    class_ndvi_osm = mass_to_class(ndvi_osm_mass)

    class_ndvi_osm_geo = gpd.GeoDataFrame(
        class_ndvi_osm.join(lu_polygons.loc[:, "geometry"])
    )
    # class_ndvi_osm_geo = class_ndvi_osm_geo.join(class_osm, rsuffix="_osm")
    # class_ndvi_osm_geo = class_ndvi_osm_geo.join(class_ndvi, rsuffix="_ndvi")

    class_ndvi_osm_geo.crs = lu_polygons.crs

    class_ndvi_osm_geo.to_file(lu_polygons_file_fused)

    return class_ndvi_osm_geo
