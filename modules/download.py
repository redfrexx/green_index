#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Download OSM data using ohsome API"""

__author__ = "Christina Ludwig, GIScience Research Group, Heidelberg University"
__email__ = "christina.ludwig@uni-heidelberg.de"

import os
from ohsome import OhsomeClient, OhsomeException
from modules.utils import load_config


def build_ohsome_filters(layer: dict):
    """
    Builds the filter string based a dictionary of tags
    :return:
    """
    tag_filters = []
    for key, values in layer["tags"].items():
        if isinstance(values, list):
            tag_filter = "{0} in ({1})".format(key, ", ".join(values))
        elif isinstance(values, str):
            tag_filter = f"{key}={values}"
        else:
            raise TypeError(f"{type(values)} not supported for tag values")
        tag_filters.append(tag_filter)
    ohsome_filter = "({0})".format(" or ".join(tag_filters))
    if "geoms" in layer.keys():
        geom_filter = " or ".join(["geometry:{0}".format(i) for i in layer["geoms"]])
        ohsome_filter += f" and ({geom_filter})"
    if "types" in layer.keys():
        type_filter = " or ".join(["type:{0}".format(i) for i in layer["types"]])
        ohsome_filter += f" and ({type_filter})"

    return ohsome_filter


def download_features(bbox, layers, outdir, timestamp=None):
    """
    Downloads all OSM highways for the specified timestamp and bbox

    :param bbox:
    :param aoi_name:
    :param config:
    :param db:
    :return:
    """
    client = OhsomeClient()

    for layer in layers:
        ohsome_filter_str = build_ohsome_filters(layer)

        response = client.elements.geometry.post(
            bboxes=bbox, time=timestamp, filter=ohsome_filter_str, properties="tags"
        )
        response.to_json(os.path.join(outdir, f"{layer['name']}.geojson"))


if __name__ == "__main__":

    config_file = "./config/config.json"
    config = load_config(config_file)

    landuse_feature_dir = os.path.join(config["outdir"], "ohsome", "landuse")
    if not os.path.exists(landuse_feature_dir):
        os.mkdir(landuse_feature_dir)
    download_features(
        bbox=config["bbox"],
        timestamp=config["timestamp"],
        layers=config["landuse_features"],
        outdir=landuse_feature_dir,
    )

    traffic_feature_dir = os.path.join(config["outdir"], "ohsome", "traffic")
    if not os.path.exists(traffic_feature_dir):
        os.mkdir(traffic_feature_dir)
    download_features(
        bbox=config["bbox"],
        timestamp=config["timestamp"],
        ohsome_filter=config["traffic_features"],
        outdir=traffic_feature_dir,
    )
