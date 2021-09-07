#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Simulates routes using openrouteservice"""

__author__ = "Christina Ludwig, GIScience Research Group, Heidelberg University"
__email__ = "christina.ludwig@uni-heidelberg.de"


import logging
import os
import random
import itertools
import subprocess
from shapely.geometry import Point
import requests
import copy
import fiona
import shutil
import geopandas as gpd
import openrouteservice as ors
from shapely.geometry import LineString, MultiLineString
import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt

LOW_GREEN_THRESHOLD = 10
HIGH_GREEN_THRESHOLD = 60
GREEN_WEIGHT = 1.0


class Route(object):
    """Route calculated using openrouteservice"""

    def __init__(
        self,
        params=None,
        base_url=None,
        key=None,
        profile="foot-walking",
        fmt="geojson",
        file=None,
    ):
        """
        Initializes parameters and sends request to ORS server

        :param params: dict
        :param base_url: string
        """

        self.__dataframe = None
        self.__headers = {
            "headers": {
                "Accept": "application/json, application/geo+json, application/gpx+xml, img/png; charset=utf-8",
                "Content-Type": "application/json; charset=utf-8",
            }
        }

        self.params = params
        if base_url:
            self.base_url = base_url
        else:
            self.base_url = "https://api.openrouteservice.org/"
        self.key = key
        self.profile = profile
        self.fmt = fmt

        if file is None:
            self.__json_response = self.request_route()
        else:
            self.load(file)

    def load(self, file):
        """
        Loads route from file
        :param file:
        :return:
        """
        with open(file) as src:
            self.__json_response = json.load(src)

    def request_route(self):
        """
        Send route request to ORS server

        :param params: dict containing request parameters
        :param base_url: string of ORS base url
        :return: dict of ORS response
        """
        if self.base_url is not None:
            client = ors.Client(base_url=self.base_url)
        else:
            client = ors.Client(base_url=self.base_url, key=self.key)

        # Send request
        try:
            return client.request(
                url="v2/directions/{0}/{1}".format(self.profile, self.fmt),
                post_json=self.params,
                requests_kwargs=self.__headers,
                get_params=[],
            )
        except ors.exceptions.ApiError as e:
            raise ValueError(e)
        finally:
            del client

    @property
    def json_response(self):
        """
        Returns the ORS response as a dictionary
        :return: dict
        """
        return self.__json_response

    @property
    def coordinates(self):
        """
        Returns the coordinates of the route from the ORS response
        :return: list of coordinates
        """
        return self.json_response["features"][0]["geometry"]["coordinates"]

    @property
    def extras(self):
        """
        Returns the extra information from the ORS response
        :return:
        """
        try:
            return self.json_response["features"][0]["properties"]["extras"]
        except Exception:
            return None

    def values(self, criteria):
        """
        Returns the values for a certain criterion
        :param criterion: 'green', 'noise' or 'steepness'
        :return: values of criterion along route
        """
        return np.concatenate(
            [np.repeat(v[2], v[1] - v[0]) for v in self.extras[criteria]["values"]]
        )

    @property
    def green_exposure(self):
        """
        Returns the overall exposure to greenness of the route
        :return: green exposure
        """
        summary = self.summary("green")
        return sum(summary["value"] * summary["distance"]) / summary["distance"].sum()

    @property
    def solar_exposure(self):
        """
        Returns the overall exposure to greenness of the route
        :return: green exposure
        """
        summary = self.summary("shadow")
        return sum(summary["value"] * summary["distance"]) / summary["distance"].sum()

    @property
    def steepness_exposure(self):
        """
        Returns the overall exposure to positive and negative steepness of the route
        :return: steepness exposure for negative and positive values
        """
        summary = self.summary("steepness")
        pos = []
        neg = []
        dist_Neg = []
        dist_Pos = []
        for o in range(len(summary["value"])):
            if summary["value"][o] > 0:
                pos.append(summary["value"][o])
                dist_Pos.append(summary["distance"][o])
            else:
                neg.append(summary["value"][o])
                dist_Neg.append(summary["distance"][o])

        if sum(dist_Neg) != 0:
            res2 = sum(np.array(neg) * np.array(dist_Neg)) / sum(dist_Neg)
        else:
            res2 = np.nan
        if sum(dist_Pos) != 0:
            res1 = sum(np.array(pos) * np.array(dist_Pos)) / sum(dist_Pos)
        else:
            res1 = np.nan
        return [res1, res2]

    @property
    def noise_exposure(self):
        """
        Returns the overall exposure to noise of the route
        :return: noise exposure
        """
        summary = self.summary("noise")
        return sum(summary["value"] * summary["distance"]) / summary["distance"].sum()

    @property
    def duration(self):
        """
        Returns the overall duration of the route
        :return: Duration
        """
        return self.json_response["features"][0]["properties"]["summary"]["duration"]

    @property
    def distance(self):
        """
        Returns the overall distance of the route
        :return: Distance
        """
        return self.json_response["features"][0]["properties"]["summary"]["distance"]

    @property
    def descent(self):
        """
        Returns the overall distance of the route
        :return: Distance
        """
        return self.json_response["features"][0]["properties"]["descent"]

    @property
    def ascent(self):
        """
        Returns the overall distance of the route
        :return: Distance
        """
        return self.json_response["features"][0]["properties"]["ascent"]

    def summary(self, criterion):
        """
        Returns the summary for a certain criterion of the ORS response as a pandas dataframe
        :param criterion: 'green', 'noise' or 'steepness'
        :return: Dataframe with summary
        """
        if criterion in self.extras.keys():
            return pd.DataFrame(self.extras[criterion]["summary"])
        else:
            raise ValueError("criterion '%s' does not exist.")

    def plot_summary(self, criterion):
        """
        Returns a bar plot of the summary for a certain criterion
        :param criterion: 'green', 'noise' or 'steepness'
        :return: Bar plot showing summary
        """
        summary = self.summary(criterion)
        return plt.bar(x=summary["value"], height=summary["amount"], color="green")

    @property
    def route_segments(self):
        """
        Returns segments of the route
        :return: list of LineStrings
        """
        n_segments = len(self.coordinates) - 1
        segments = []
        for i in range(0, n_segments):
            segments.append(LineString(self.coordinates[i : i + 2]))
        return segments

    # todo write test for this function
    def as_dataframe(self):
        """
        Converts the route and its extra information into a geopandas dataframe
        :return: GeoDataFrame with route information
        """
        if self.__dataframe is not None:
            return self.__dataframe
        else:
            df = gpd.GeoDataFrame({"geometry": self.route_segments}, crs="epsg:4326")
            if self.extras:
                for k in self.extras.keys():
                    df[k] = self.values(k)
                # Dissolve line strings
                columns = list(df.columns.drop("geometry"))
                df = df.dissolve(by=columns, as_index=True).reset_index()
                df = df.loc[~df.is_empty]
                df.geometry = df.geometry.apply(
                    lambda x: MultiLineString([x]) if isinstance(x, LineString) else x
                )
        self.__dataframe = df
        return self.__dataframe

    def to_geojson(self, outfile, driver):
        """
        Writes the route to a geojson file
        :param outfile: Path to output file as string
        :return:
        """
        self.as_dataframe().to_file(outfile, driver=driver)

    def plot(self, *args, **kwargs):
        """
        Plots the route on a map
        :param args:
        :param kwargs:
        :return: plotted route
        """
        return self.as_dataframe().plot(*args, **kwargs)

    def to_file(self, outfile):
        """
        Writes the whole response to file.
        :param outfile:
        :return:
        """
        with open(outfile, "w") as dst:
            json.dump(self.json_response, dst, indent=4)


def generate_random(polygon):
    """
    Generates a random points within a polygon geometry
    :param polygon:
    :return:
    """
    minx, miny, maxx, maxy = polygon.bounds
    while True:
        pnt = Point(random.uniform(minx, maxx), random.uniform(miny, maxy))
        if polygon.contains(pnt):
            yield pnt


def merge_tmp_files(out_gpkg, tmp_gpkg, crs="epsg:4326"):
    """
    Merge all temporary data chunks
    :param outfile:
    :param name:
    :return:
    """
    name = os.path.basename(tmp_gpkg).split(".")[0]
    cmd = [
        "ogrmerge.py",
        "-single",
        "-t_srs",
        f"epsg={crs}",
        "-nln",
        name,
        "-overwrite_layer",
        "-o",
        out_gpkg,
        tmp_gpkg,
    ]
    try:
        subprocess.check_call(cmd)
    except Exception:
        cmd = [
            "ogrmerge.py",
            "-single",
            "-t_srs",
            f"epsg={crs}",
            "-nln",
            name,
            "-f",
            "GPKG",
            "-o",
            out_gpkg,
            tmp_gpkg,
        ]
        subprocess.check_call(cmd)
    os.unlink(tmp_gpkg)


def generate_routes(config, debug=False):
    """
    Simulate pedestiran and bike routes using openrouteservice instance
    :param config:
    :param debug:
    :return:
    """

    # Config parameters
    random.seed(10)

    # Config file params
    ors_url = config["ors_url"]
    steep_level = config["steep_level"]
    profile = config["profile"]
    n_routes = config["n_routes"]
    districts_file = config["districts_file"]
    col_name = config["col_name"]
    out_dir = config["out_dir"]
    aoi_name = config["name"]
    min_length = config["min_length"]
    max_length = config["max_length"]
    preference = "fastest"

    logger = logging.getLogger("comparison")

    # Create output directories -------------------------
    out_dir_job = os.path.join(
        out_dir, aoi_name, f"{aoi_name}_{profile}_{steep_level}_{n_routes}"
    )
    out_gpkg = os.path.join(out_dir_job, "routes.gpkg")

    if os.path.exists(out_dir_job):
        shutil.rmtree(out_dir_job)
    os.makedirs(out_dir_job, exist_ok=True)

    try:
        requests.get(ors_url + "health")
    except requests.exceptions.ConnectionError:
        logger.critical(f"Connection error: {ors_url} is not available.")

    # Load districts
    districts = gpd.read_file(districts_file)
    districts[col_name] = districts[col_name].map(lambda x: x.replace(" ", "_"))
    districts.set_index(col_name, inplace=True)

    pairs = list(itertools.combinations_with_replacement(list(districts.index), 2))

    # Request parameters
    body_normal = {
        "instructions": "false",
        "preference": preference,
        "extra_info": ["green"],
        "elevation": False,
        "continue_straight": True,
        "options": {
            "avoid_features": ["ferries"],
            "profile_params": {
                "weightings": {"steepness_difficulty": steep_level, "green": 0.0}
            },
        },
    }
    body_green = copy.deepcopy(body_normal)
    body_green["options"]["profile_params"]["weightings"]["green"] = GREEN_WEIGHT

    route_number = 1
    for start_idx, end_idx in pairs:

        start_district = districts.loc[start_idx]
        end_district = districts.loc[end_idx]

        start_generator = generate_random(start_district.geometry)
        end_generator = generate_random(end_district.geometry)

        normal_route_geoms = []
        green_route_geoms = []
        normal_routes_stats = []
        green_routes_stats = []
        same_low_green = []
        same_high_green = []
        calculated_routes_for_pair = 0
        tries = 0
        while calculated_routes_for_pair < n_routes and tries < n_routes * 3:

            tries += 1
            start = next(start_generator)
            end = next(end_generator)

            try:
                body_normal["coordinates"] = list(start.coords) + list(end.coords)
                normal_route = Route(
                    params=body_normal, base_url=ors_url, profile=profile, fmt="geojson"
                )

                body_green["coordinates"] = list(start.coords) + list(end.coords)
                green_route = Route(
                    params=body_green, base_url=ors_url, profile=profile, fmt="geojson"
                )

                invalid = normal_route.duration > green_route.duration * 1.02

                if normal_route.distance < min_length:
                    logger.info(
                        f"Route is too short: {normal_route.distance} m. Generate new route."
                    )
                    invalid = True

                if normal_route.distance > max_length:
                    logger.info(
                        f"Route is too long: {normal_route.distance} m. Generate new route."
                    )
                    invalid = True

            except (TimeoutError, ValueError) as e:
                logger.exception(e)
                continue
            except requests.exceptions.ConnectionError:
                logger.exception(f"Connection error: {ors_url} is not available.")
                continue
            except Exception:
                logger.exception("Error during ORS request: ")
                continue

            if invalid:
                continue
            else:
                normal_route_geoms.append(normal_route.as_dataframe().cascaded_union)
                normal_routes_stats.append(
                    [
                        route_number,
                        normal_route.green_exposure,
                        normal_route.distance,
                        normal_route.duration,
                    ]
                )
                green_route_geoms.append(green_route.as_dataframe().cascaded_union)
                green_routes_stats.append(
                    [
                        route_number,
                        green_route.green_exposure,
                        green_route.distance,
                        green_route.duration,
                    ]
                )

                normal_route_df = normal_route.as_dataframe()
                green_route_df = green_route.as_dataframe()

                # Get low greenness part of routes
                normal_low_green = normal_route_df.loc[
                    normal_route_df.green <= LOW_GREEN_THRESHOLD
                ].cascaded_union
                green_low_green = green_route_df.loc[
                    green_route_df.green <= LOW_GREEN_THRESHOLD
                ].cascaded_union
                same_geom = green_low_green.intersection(normal_low_green)
                if (
                    isinstance(same_geom, MultiLineString)
                    or isinstance(same_geom, LineString)
                ) and not same_geom.is_empty:
                    same_low_green.append(same_geom)
                else:
                    same_low_green.append(None)

                # Get high greenness part of routes
                normal_high_green = normal_route_df.loc[
                    normal_route_df.green >= HIGH_GREEN_THRESHOLD
                ].cascaded_union
                green_high_green = green_route_df.loc[
                    green_route_df.green >= HIGH_GREEN_THRESHOLD
                ].cascaded_union
                same_geom = green_high_green.intersection(normal_high_green)
                if (
                    isinstance(same_geom, MultiLineString)
                    or isinstance(same_geom, LineString)
                ) and not same_geom.is_empty:
                    same_high_green.append(same_geom)
                else:
                    same_high_green.append(None)

                route_number += 1
                calculated_routes_for_pair += 1

        if len(normal_routes_stats) == 0:
            logging.critical("No routes generated.")
            continue

        # Writes chunks of data to file
        route_ids, green_factor, distance, duration = zip(*normal_routes_stats)
        tmp_normal_routes_df = gpd.GeoDataFrame(
            {
                "geometry": normal_route_geoms,
                "start": np.repeat(start_idx, n_routes),
                "end": np.repeat(end_idx, n_routes),
                "route_id": route_ids,
                "green": green_factor,
                "distance": distance,
                "duration": duration,
            },
            crs="epsg:4326",
        )
        route_ids, green_factor, distance, duration = zip(*green_routes_stats)
        tmp_green_routes_df = gpd.GeoDataFrame(
            {
                "geometry": green_route_geoms,
                "start": np.repeat(start_idx, n_routes),
                "end": np.repeat(end_idx, n_routes),
                "route_id": route_ids,
                "green": green_factor,
                "distance": distance,
                "duration": duration,
            },
            crs="epsg:4326",
        )

        tmp_low_green_df = gpd.GeoDataFrame(
            {"geometry": same_low_green, "route_id": route_ids}, crs="epsg:4326"
        ).dropna()
        tmp_high_green_df = gpd.GeoDataFrame(
            {"geometry": same_high_green, "route_id": route_ids}, crs="epsg:4326"
        ).dropna()

        # Write routes to temporary files
        try:
            tmp_normal_routes_df.to_file(
                out_gpkg, layer="normal", mode="a", driver="GPKG"
            )
        except OSError:
            tmp_normal_routes_df.to_file(out_gpkg, layer="normal", driver="GPKG")
        except ValueError as e:
            logging.info(e)
            continue

        try:
            tmp_green_routes_df.to_file(
                out_gpkg, layer="green", mode="a", driver="GPKG"
            )
        except fiona.errors.DriverError:
            tmp_green_routes_df.to_file(out_gpkg, layer="green", driver="GPKG")
        except ValueError as e:
            logging.info(e)
            continue

        try:
            tmp_low_green_df.to_file(
                out_gpkg, layer="low_green", mode="a", driver="GPKG"
            )
        except fiona.errors.DriverError:
            tmp_low_green_df.to_file(out_gpkg, layer="low_green", driver="GPKG")
        except ValueError as e:
            logging.info(e)
            continue

        try:
            tmp_high_green_df.to_file(
                out_gpkg, layer="high_green", mode="a", driver="GPKG"
            )
        except fiona.errors.DriverError:
            tmp_high_green_df.to_file(out_gpkg, layer="high_green", driver="GPKG")
        except ValueError as e:
            logging.info(e)
            continue

        del (
            green_factor,
            distance,
            duration,
            green_routes_stats,
            normal_routes_stats,
            normal_route_geoms,
            green_route_geoms,
            tmp_low_green_df,
        )  # , tmp_high_green_df,

        avoided_segments = []
        prefered_segments = []
        shared = []
        for r_id in range(len(tmp_normal_routes_df)):
            green_geom = tmp_green_routes_df.loc[r_id, "geometry"]
            normal_geom = tmp_normal_routes_df.loc[r_id, "geometry"]
            avoided = normal_geom.difference(green_geom)
            prefered = green_geom.difference(normal_geom)
            shared.append((normal_geom.length - avoided.length) / normal_geom.length)
            avoided_segments.append(avoided)
            prefered_segments.append(prefered)

        tmp_green_routes_df["shared"] = shared
        avoided_df = gpd.GeoDataFrame({"geometry": avoided_segments}, crs="epsg:4326")
        prefered_df = gpd.GeoDataFrame({"geometry": prefered_segments}, crs="epsg:4326")
        avoided_df = avoided_df.loc[~avoided_df.geometry.map(lambda x: x.is_empty)]
        prefered_df = prefered_df.loc[~prefered_df.geometry.map(lambda x: x.is_empty)]

        try:
            avoided_df.to_file(out_gpkg, layer="avoided", mode="a", driver="GPKG")
        except fiona.errors.DriverError:
            avoided_df.to_file(out_gpkg, layer="avoided", driver="GPKG")
        except ValueError as e:
            logging.info(e)
            continue

        try:
            prefered_df.to_file(out_gpkg, layer="preferred", mode="a", driver="GPKG")
        except fiona.errors.DriverError:
            prefered_df.to_file(out_gpkg, layer="preferred", driver="GPKG")
        except ValueError as e:
            logging.info(e)
            continue

        del tmp_normal_routes_df, tmp_green_routes_df
