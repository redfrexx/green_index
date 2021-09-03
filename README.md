# (Green) Index Calculation for OpenRouteService

This repository contains source code to

1. Calculate the greenness (i.e. presence of vegetation) of street blocks based on OpenStreetMap and Sentinel-2 data
2. Calculate the green index for each OSM highway feature to be used in an openrouteservice instance to calculate green routes.

## Installation

Python 3.x is required and the packages listed in `requirements.txt`. You can set up a new python environment with all dependencies using `pip`:

```
$ python3 -m venv env
$ source env/bin/activate
$ python3 -m pip install -r requirements.txt
```

## Usage

### 1. Greenness calculation

The first step is to calculate the greenness of individual street blocks within the area of interest using the `calculate_greenness.py` script which takes two required parameters:

```
$ calculate_greenness.py -c config_sample.json -g google_credentials_sample.json
```

#### Configuration file: -c / --config

All parameters need to specify the area of interest must be given in a configuration file.

``` json
{
  "name": "MA_LU_debug",
  "bbox": [8.46874,49.4971,8.49213,49.50995],
  "epsg": "32632",
  "timestamp": "2021-08-15",
  "cloud_coverage": 5,
  "ndvi_year": 2020,
  "output_dir": "./data"
  "fuzzy_centers": {
    "green": 0.71,
    "mixed": 0.43,
    "grey": 0.15,
    "d": 0.094
  },
}
```

| Parameter | Explanation                                           |
|-----------|-------------------------------------------------------|
| name | The name of the area of interest. Used in output file names. |
| bbox | Bounding box of area of interest in geographic coordinates, format: (minx, miny, max, maxy)|
| epsg | The epsg of a projected coordinate reference system suitable for the are of interest.|
| timestamp | Timestamp of the OSM data used for processing. |
| ndvi_year| Year used to calculate the annual maximum NDVI. |
| output_dir | Path to a existing directory in which output data will be stored. |
| fuzzy centers | NDVI values used as centers for the classes |


#### Google Credentials: -g / --google_cred

A json file containing your personal credentials need to use Google Earth Engine. You need to [create a service account with google](https://developers.google.com/earth-engine/guides/service_account) and generate a key.


### 2. Green Index Calculation

After the greenness is calculated, the green index of each OSM highway feature can be calculated.


```
$ python calculate_index.py -h

Usage: calculate_index.py [-h] --bbox BBOX [--timestamp TIMESTAMP] --width WIDTH [--vector VECTOR_FILE] [--raster RASTER_FILE] --outputdirectory OUTPUT_DIR

Calculates the index of each OSM highway based on provided raster or vector file. The highways are downloaded using the ohsome API.

optional arguments:
  -h, --help            show this help message and exit
  --bbox BBOX, -b BBOX  Bounding box in geographic coordinates as string without whitespace e.g. 'minx,miny,maxx,maxy')
  --timestamp TIMESTAMP, -t TIMESTAMP
                        ISO formatted timestamp for download of highways from OSM. By default the latest timestamp available in the ohsome API will be used.
  --width WIDTH, -w WIDTH
                        Width of the buffer around highway segments in meters
  --vector VECTOR_FILE, -v VECTOR_FILE
                        Path to vector file containing features counted nearby highways
  --raster RASTER_FILE, -r RASTER_FILE
                        Path to raster file used to calculate mean value within area nearby highway
  --outputdirectory OUTPUT_DIR, -o OUTPUT_DIR
                        Path to existing output directory.

```

#### Examples:

```
$ calculate_index.py -b 13.7176,51.0298,13.7957,51.0731 -r green.tif -w 20 -o index
```

```
$ calculate_index.py -b 13.7176,51.0298,13.7957,51.0731 -v benches.geojson -w 20
```

### 3. Set up OpenRouteService instance

Last step is to set up an instance of this [openrouteservice repository](https://github.com/redfrexx/openrouteservice/tree/shadow-trees).
Put the green index csv file and an OSM file covering your area of interest in the folder `openrouteservice/docker/data`. Replace the files in the  openrouteservice/docker/docker-compose.yml accordingly.
