#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""__description__
"""

__author__ = "Christina Ludwig, GIScience Research Group, Heidelberg University"
__email__ = "christina.ludwig@uni-heidelberg.de"

import os
from ohsome import OhsomeClient, OhsomeException
from importlib import reload
from utils import build_ohsome_filters, load_config
import logging



def download_traffic_features(bbox, timestamp, traffic_tags, outdir):
    """
    Downloads all OSM highways for the specified timestamp and bbox

    :param bbox:
    :param aoi_name:
    :param config:
    :param db:
    :return:
    """

    logger = logging.getLogger("root." + __name__)

    # Get highway features form ohsome api -----------------------------------------
    client = OhsomeClient()
    ohsome_filters = build_ohsome_filters(traffic_tags)

    for name, fltr in ohsome_filters.items():
        print(name)
        try:
            response = client.elements.geometry.post(bboxes=bbox, time=timestamp, filter=fltr)
        except OhsomeException as e:
            print(e)
            continue
        response.to_json(os.path.join(outdir, "{}.geojson".format(name)))


if __name__ == "__main__":

    config_file = "./config/config_dresden.json"
    config = load_config(config_file)

    download_traffic_features(bbox=config["bbox"],
                              timestamp=config["timestamp"],
                              traffic_tags=config["traffic_features"],
                              outdir=os.path.join(config["outdir"], "ohsome"))