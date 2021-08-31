#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for utils"""

__author__ = "Christina Ludwig, GIScience Research Group, Heidelberg University"
__email__ = "christina.ludwig@uni-heidelberg.de"

from modules.utils import build_ohsome_filters, calc_compactness
from shapely.geometry import box, Polygon
import geopandas as gpd
import numpy as np


def test_calc_compactness():
    """
    Calculate the compactness of a rectangle
    :return:
    """
    minx = 0
    maxx = 2
    miny = 0
    maxy = 2
    expected_compactness = 0.78
    polygon_df = gpd.GeoDataFrame({"geometry": [Polygon(box(minx, miny, maxx, maxy))]})

    actual_compactness = calc_compactness(polygon_df)

    np.testing.assert_almost_equal(actual_compactness, expected_compactness, decimal=2)


def test_build_filters():
    """
    Test if filter string is created correctly
    :return:
    """
    tags = [
        {"tags": {"highway": ["primary", "secondary"]}, "geoms": ["polygon", "line"]},
        {"tags": {"railway": "*"}, "geoms": ["line"]},
    ]

    expected_filters = [
        "(highway in (primary, secondary)) and (geometry:polygon or geometry:line)",
        "(railway=*) and (geometry:line)",
    ]

    actual_filters = build_ohsome_filters(tags)

    assert expected_filters == actual_filters
