#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Download OSM data using ohsome API"""

__author__ = "Christina Ludwig, GIScience Research Group, Heidelberg University"
__email__ = "christina.ludwig@uni-heidelberg.de"

import os
from ohsome import OhsomeClient, OhsomeException
from modules.utils import load_config


def build_ohsome_filters(filter_items: list):
    """
    Builds the filter string based a dictionary of tags
    :return:
    """
    filters = {}
    for item in filter_items:
        if "geoms" in item.keys():
            geom_filter = " or ".join(["geometry:{0}".format(i) for i in item["geoms"]])
        if "types" in item.keys():
            type_filter = " or ".join(["type:{0}".format(i) for i in item["types"]])
        for key, values in item["tags"].items():
            if isinstance(values, list):
                tag_filter = "{0} in ({1})".format(key, ", ".join(values))
            elif isinstance(values, str):
                tag_filter = f"{key}={values}"
            else:
                raise TypeError(f"{type(values)} not supported for tag values")
            filter = tag_filter
            if "geoms" in item.keys():
                filter += f" and ({geom_filter})"
            if "types" in item.keys():
                filter += f" and ({type_filter})"

            filters[item["name"]] = filter
    return filters


def download_features(bbox, timestamp, tags, outdir):
    """
    Downloads all OSM highways for the specified timestamp and bbox

    :param bbox:
    :param aoi_name:
    :param config:
    :param db:
    :return:
    """
    client = OhsomeClient()
    ohsome_filters = build_ohsome_filters(tags)

    for name, fltr in ohsome_filters.items():
        print(name)
        try:
            response = client.elements.geometry.post(
                bboxes=bbox, time=timestamp, filter=fltr, properties="tags"
            )
        except OhsomeException as e:
            print(e)
            continue
        response.to_json(os.path.join(outdir, f"{name}.geojson"))


if __name__ == "__main__":

    config_file = "./config/config.json"
    config = load_config(config_file)

    landuse_feature_dir = os.path.join(config["outdir"], "ohsome", "landuse")
    if not os.path.exists(landuse_feature_dir):
        os.mkdir(landuse_feature_dir)
    download_features(
        bbox=config["bbox"],
        timestamp=config["timestamp"],
        tags=config["landuse_features"],
        outdir=landuse_feature_dir,
    )

    traffic_feature_dir = os.path.join(config["outdir"], "ohsome", "traffic")
    if not os.path.exists(traffic_feature_dir):
        os.mkdir(traffic_feature_dir)
    download_features(
        bbox=config["bbox"],
        timestamp=config["timestamp"],
        tags=config["traffic_features"],
        outdir=traffic_feature_dir,
    )
