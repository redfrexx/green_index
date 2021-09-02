# (Green) Index Calculation for OpenRouteService


## Usage

### 1. Greenness calculation

```
$ calculate_greenness.py -a config_sample.json -g google_credentials_sample.json
```


### 2. Index calculation

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

```
$ calculate_index.py -b 13.7176,51.0298,13.7957,51.0731 -r green.tif -w 20 -o index
```

```
$ calculate_index.py -b 13.7176,51.0298,13.7957,51.0731 -v benches.geojson -w 20
```
