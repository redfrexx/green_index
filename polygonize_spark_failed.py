#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""__description__
"""

__author__ = "Christina Ludwig, GIScience Research Group, Heidelberg University"
__email__ = "christina.ludwig@uni-heidelberg.de"


import geopandas as gpd
from pyspark.sql import SparkSession
from sedona.register import SedonaRegistrator
from sedona.utils import KryoSerializer, SedonaKryoRegistrator

spark = SparkSession.builder. \
    master("local[*]"). \
    appName("Sedona App"). \
    config("spark.serializer", KryoSerializer.getName). \
    config("spark.kryo.registrator", SedonaKryoRegistrator.getName). \
    config("spark.jars.packages",
           "org.apache.sedona:sedona-python-adapter-3.0_2.12:1.0.0-incubating,org.datasyslab:geotools-wrapper:geotools-24.1"). \
    getOrCreate()

SedonaRegistrator.registerAll(spark)

highways_file = "/Users/chludwig/Development/meinGruen/code/green_index/data/(highway_in_(motorway,_trunk,_primary,_secondary,_tertiary,_residential,_unclassified,_motorway_link,_trunk_link,_primary_link,_secondary_link,_tertiary_link,_living_street))_and_(geometry:polygon_or_geometry:line).geojson"
highways = gpd.read_file(highways_file)

highways_spark = spark.createDataFrame(highways)
type(highways_spark)
highways_spark.printSchema()
highways_spark.show()
highways_spark.createOrReplaceTempView("highways")

# Node line network
spark.sql("SELECT ST_Node(geometry) FROM highways;").show()


res = spark.sql("SELECT st_polygonize(geometry) As geomtextrep FROM (SELECT geometry FROM highways_spark) As foo;")

point_csv_df = spark.read.format("csv").\
    option("delimiter", ",").\
    option("header", "false").\
    load("./data/testpoint.csv")

point_csv_df.createOrReplaceTempView("pointtable")
type(point_csv_df)
point_df = spark.sql("select ST_Point(cast(pointtable._c0 as Decimal(24,20)), cast(pointtable._c1 as Decimal(24,20))) as arealandmark from pointtable")
point_df.show(5)