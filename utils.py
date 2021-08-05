#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""__description__
"""

__author__ = "Christina Ludwig, GIScience Research Group, Heidelberg University"
__email__ = "christina.ludwig@uni-heidelberg.de"

import os
import json
import geopandas as gpd


def load_config(config_file):
    """
    Load config file
    :return:
    """
    assert os.path.exists(config_file), "Configuration file does not exist: {}".format(os.path.abspath(config_file))
    with open(config_file, 'r') as src:
        config = json.load(src)
    return config


def build_ohsome_filters(filter_items: list):
    """
    Builds the filter string based a dictionary of tags
    :return:
    """
    filters = {}
    for item in filter_items:
        geom_filter = " or ".join(["geometry:{}".format(i) for i in item["geoms"]])
        for key, values in item["tags"].items():
            if isinstance(values, list):
                tag_filter = "{} in ({})".format(key, ", ".join(values))
            elif isinstance(values, str):
                tag_filter = "{}={}".format(key, values)
            else:
                raise TypeError("{} not supported for tag values".format(type(values)))
            filters[item["name"]] = ("({}) and ({})".format(tag_filter, geom_filter))
    return filters



def split_highways(highways):
    """
    Split highway features at intersections
    :param highways:
    :return:
    """
    new_highways = gpd.GeoDataFrame()
    for splitted_geom in highways.geometry.unary_union:
        part = gpd.GeoDataFrame([[splitted_geom]], columns=['geometry'])
        new_highways = new_highways.append(part)
    return new_highways