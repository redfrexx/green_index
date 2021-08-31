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

GREEN_LABEL_DICT = {"green": "a", "grey": "b", "green_grey": "ab"}
PUBLIC_LABEL_DICT = {"public": "a", "private": "b", "public_private": "ab"}


def convert_to_mass(data, label_dict, adapt_masses=False):
    """
    Converts data frame to series of probability masses
    :return
    """
    # Convert to masses
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
        if row_cp["a"] > 0.95:
            row_cp["ab"] = row_cp["ab"] + row_cp["a"] - 0.95
            row_cp["a"] = 0.95
        elif row_cp["b"] > 0.95:
            row_cp["ab"] = row_cp["ab"] + row_cp["b"] - 0.95
            row_cp["b"] = 0.95

    return MassFunction(dict(row_cp))


def classify_green(row):
    """
    Classifies beliefs into green or grey
    :returns
    """
    uc_green = row.pl("a") - row.bel("a")
    uc_grey = row.pl("b") - row.bel("b")
    prob = row.pignistic()
    if prob["a"] > prob["b"]:
        label = "green"
    elif prob["a"] < prob["b"]:
        label = "grey"
    else:
        label = "uncertain"
    return uc_green, uc_grey, prob["a"], prob["b"], label


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
        "uc_green",
        "uc_grey",
        "prob_green",
        "prob_grey",
        "pred_class",
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


def fuse(aoi_name, config):
    """
    Fuses information from OSM and Sentinel-2 on greenness
    :param aoi_name:
    :param config:
    :return:
    """

    green_results_dir = os.path.join(config["output_dir"], aoi_name, "green")
    ndvi_s2_file = os.path.join(green_results_dir, f"{aoi_name}_green_s2.csv")
    ndvi_osm_file = os.path.join(green_results_dir, f"{aoi_name}_green_osm.csv")
    target_file = os.path.join(config["aois"][aoi_name]["target_geoms"])

    s2_data = pd.read_csv(ndvi_s2_file).set_index("TARGET_ID")
    osm_data = pd.read_csv(ndvi_osm_file).set_index("TARGET_ID")

    targets = gpd.read_file(target_file).set_index("TARGET_ID")
    targets = targets.to_crs(epsg=config["aois"][aoi_name]["epsg"])
    targets["area"] = targets.area
    targets["TARGET_type"].replace(np.nan, "None", inplace=True)

    green_label_dict = {"green": "a", "grey": "b", "green_grey": "ab"}

    s2_mass = convert_to_mass(s2_data, label_dict=green_label_dict)
    osm_mass = convert_to_mass(osm_data, label_dict=green_label_dict)

    s2_osm_mass = fuse_masses(s2_mass, osm_mass)

    class_s2 = mass_to_class(s2_mass)
    class_osm = mass_to_class(osm_mass)
    class_s2_osm = mass_to_class(s2_osm_mass)

    class_s2_osm_geo = gpd.GeoDataFrame(class_s2_osm.join(targets.loc[:, "geometry"]))
    class_s2_osm_geo = class_s2_osm_geo.join(class_osm, rsuffix="_osm")
    class_s2_osm_geo = class_s2_osm_geo.join(class_s2, rsuffix="_s2")

    class_s2_osm_geo.crs = targets.crs
    class_s2_osm_geo.to_file(
        os.path.join(config["output_dir"], aoi_name, "green", "greenness.shp"),
        driver="ESRI Shapefile",
    )
