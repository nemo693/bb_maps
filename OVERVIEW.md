# Forest Damage Detection

Service to the forest service of South Tyrol

## 1. Overall Application Overview

This application is a multi-stage pipeline designed for detecting forest damage/mortality using satellite imagery (Sentinel-2 BOA data, preprocessed in FORCE). The core of the detection process is the *fordead* python package. Data preprocessing and postprocessing are handled by R scripts instead.

In the preprocessing part and the *fordead* part, the pipeline processes data separately for each grid cell (FORCE tiles). In the postprocessing sections, it mosaics the data and works on the whole area. During the postprocessing, "historical" data is integrated, and final yearly and monthly damage products in both raster and vector formats are generated.
