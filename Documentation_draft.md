# Forest Damage Detection

Service to the forest service of South Tyrol

## 1. Overall Application Overview

This application is a multi-stage pipeline designed for detecting forest damage/mortality using satellite imagery (Sentinel-2 BOA data, preprocessed in FORCE). The core of the detection process is the *fordead* python package. Data preprocessing and postprocessing are handled by R scripts instead.

In the preprocessing part and the *fordead* part, the pipeline processes data separately for each grid cell (FORCE tiles). In the postprocessing sections, it mosaics the data and works on the whole area. During the postprocessing, "historical" data is integrated, and final yearly and monthly damage products in both raster and vector formats are generated.

## 2. Scripts composing the pipeline

The application consists of several R scripts and a Python package (`fordead_plain`) with specific modules. The python package script is saved in a separate folder, since it was modified and it requires custom installation (see installation notes).

### R Scripts (Orchestration and Post-processing)

-   **`update_province.R`**
    -   **Purpose:** The top-level orchestrator of the entire pipeline. It sets the working directory, defines output paths, manages backups of previous runs, and sequentially calls all other R scripts. It also includes commands for generating raster overviews on the final layers and uploading them to the maps portal.
    -   **Inputs:**
        -   Previous run files (from `/mnt/CEPH_PROJECTS/WALDSCHAEDEN/working_folder/outputs/fordead_15`). Used for backing these up.
        -   Hardcoded paths to previous yearly/monthly shapefiles for area comparison.
    -   **Outputs:**
        -   Backed-up previous run files.
        -   Final yearly and monthly damage GeoTIFFs and Shapefiles (via `gdaladdo` and `scp` commands).
    -   **Key Actions:**
        -   Sets `setwd()`.
        -   Defines `update_name` and `outfold`.
        -   Moves files from `fordead_15` to a backup directory.
        -   Sources `B_RUN_simple.R`, `2B_merge_tiled_outputs.R`, `3B_index_merging.R`, `4-5B_revised.R`, `6-7B_revised.R`.
        -   Executes `gdaladdo` commands.
        -   Provides `scp` commands for file transfer.
        -   Performs area comparisons using `terra::vect` and `expanse`.
-   **`B_RUN_simple.R`**
    -   **Purpose:** Manages the acquisition and organization of Sentinel-2 BOA imagery for processing. It defines spatial grids and temporal parameters, and iterates through them, calling `1B_fordead_loop_simple.R` for each grid and detection period.
    -   **Inputs:**
        -   Configuration parameters (`get_new_images`, `sentinel_only`).
        -   Spatial grid definitions (`grids`).
        -   Temporal parameters (`train_period_max`, `detection_start_dates`, `detection_end_dates`, `ignored_doy`).
    -   **Outputs:**
        -   Organized BOA imagery in `_boa` and `_update` directories for each grid.
    -   **Key Actions:**
        -   Defines global parameters for the FORDEAD process.
        -   Loops through defined `detection_start_dates` and `detection_end_dates`.
        -   For each grid, manages BOA image files, moving "future" images to an `_update` directory and restoring "past" images to the main `_boa` directory.
        -   Sources `1B_fordead_loop_simple.R` within its main loop.
-   **`1B_fordead_loop_simple.R`**
    -   **Purpose:** This script is the bridge between the R orchestration and the Python FORDEAD core. It handles initial data preprocessing, QAI mask generation (in R), and prepares and executes the Python FORDEAD scripts for each grid cell and vegetation index.
    -   **Inputs:**
        -   Parameters from `B_RUN_simple.R` (e.g., `grids`, `train_period_min`, `train_period_max`, `ignored_doy`).
        -   Raw Sentinel-2 BOA and QAI TIFFs from `in_path`.
        -   Forest mask shapefile (`forest_mask`).
    -   **Outputs:**
        -   Organized BOA data and generated QAI masks in `_boa` directories.
        -   Parameter files (`param.py`) for Python scripts.
        -   Intermediate mask files (`mask_nodata.tif`, `valid_obs.tif`).
        -   Outputs from Python FORDEAD scripts (e.g., `coeff_model.tif`, `sufficient_coverage_mask.tif`).
    -   **Key Actions:**
        -   Defines numerous settings for the FORDEAD process (paths, thresholds, indices).
        -   Manages the `qai_version` and generates/activates QAI masks based on Sentinel-2 QAI data.
        -   Prepares a `param.py` file with various parameters for the Python scripts.
        -   Calls external Python scripts (`fordead_1.py`, `fordead_2.py`, `fordead_3.py` -- likely corresponding to `step1_compute_masked_vegetationindex.py`, `step2_train_model.py`, `step3_dieback_detection.py` in the `fordead_plain/steps` directory) using `system()` calls.
        -   Performs a "dummy mask" computation to ensure algorithm stability.
-   **`2B_merge_tiled_outputs.R`**
    -   **Purpose:** Mosaics the tiled outputs (both raster and vector) generated by the FORDEAD process for each vegetation index into full-area maps.
    -   **Inputs:**
        -   Tiled FORDEAD outputs (e.g., `confidence_index.tif`, `dieback.shp`, `stress_periods.shp`, `mask_nodata.tif`, `valid_obs.tif`, `first_date_dieback.tif`, etc.) from `for_out_path`.
        -   `index.list` (implicitly from `1B_fordead_loop_simple.R`).
    -   **Outputs:**
        -   Merged GeoTIFFs and Shapefiles for various outputs (e.g., `output_NDWI_merged.shp`, `output_confidence_NDVI_merged.tif`).
    -   **Key Actions:**
        -   Loops through each vegetation index in `index.list`.
        -   Identifies and merges corresponding tiled vector files using `terra::vect` and `rbind`.
        -   Identifies and merges corresponding tiled raster files using `terra::vrt` and `writeRaster`.
-   **`3B_index_merging.R`**
    -   **Purpose:** Refines the merged outputs by applying additional masking, adjusting temporal periods, and deriving annual products. It also performs area calculations and comparisons.
    -   **Inputs:**
        -   Merged NDVI and NDWI shapefiles (e.g., `output_NDWI_merged.shp`, `output_NDVI_merged.shp`) from `fold` (which is `for_out_path`).
        -   Province boundaries shapefile (`province.shp`).
        -   LAFIS (agricultural measures) shapefile (`lafis.shp`).
        -   A base Sentinel-2 mosaic GeoTIFF (`base`) for rasterization.
        -   A reference bark beetle damage shapefile (`bark_beetle_damages_province.shp`) for comparison.
    -   **Outputs:**
        -   Masked and refined NDVI shapefiles and GeoTIFFs (e.g., `NDVI_merged_masked_with_NDWI.shp`, `NDVI_merged_masked_with_NDWI_province_lafis.shp`, `NDVI_merged_masked_with_NDWI_province_lafis.tif`).
        -   Annual damage GeoTIFF (`NDVI_merged_masked_with_NDWI_province_only_annual.tif`) and shapefile (`bark_beetle_damages_province_annual.shp`).
    -   **Key Actions:**
        -   Loads NDWI and NDVI vector data and rasterizes them.
        -   Masks the NDVI raster with the NDWI raster.
        -   Applies masking using province and LAFIS data.
        -   Adjusts date labels for November and winter periods.
        -   Derives annual damage products.
        -   Calculates and compares damage areas with a reference dataset.
-   **`4-5B_revised.R`**
    -   **Purpose:** Prepares the yearly and monthly damage products for final output and visualization by applying color palettes, projecting data, and converting formats.
    -   **Inputs:**
        -   Yearly and monthly damage rasters (e.g., `NDVI_merged_masked_with_NDWI_province_only_annual.tif`, `NDVI_merged_masked_with_NDWI_province_lafis.tif`) from `prod_fol` (which is `for_out_path`).
        -   `update_name` (from `update_province.R`).
    -   **Outputs:**
        -   Final yearly and monthly damage GeoTIFFs and Shapefiles (e.g., `changes_yearly_damages_july_2025_25832.tif`, `changes_monthly_damages_july_2025_25832.shp`).
    -   **Key Actions:**
        -   Loads yearly and monthly damage rasters.
        -   Defines and applies color palettes and labels.
        -   Projects rasters to EPSG:25832.
        -   Converts rasters to vector format.
        -   Writes final GeoTIFF and Shapefile outputs.
        -   Includes post-processing adjustments for specific date labels (November, winter periods) and recalculates `last_detectable_date`.
-   **`6-7B_revised.R`**
    -   **Purpose:** Performs the final post-processing steps, integrating previous results, applying additional filters (stress periods, shadows, single pixels), and generating the definitive output products.
    -   **Inputs:**
        -   Current yearly and monthly damage rasters (from `4-5B_revised.R`).
        -   Previous monthly and yearly damage rasters (hardcoded paths).
        -   NDVI and NDWI stress period rasters (e.g., `output_nb_periods_stress_NDVI_merged.tif`).
        -   QAI data for shadow masking (hardcoded path).
        -   Forest inspectorates shapefile (`ForestInspectorates_polygon.shp`).
        -   A base raster for projection (`base`).
    -   **Outputs:**
        -   Final yearly and monthly damage GeoTIFFs and Shapefiles (e.g., `changes_monthly_damages_july_2025_25832_final.tif`).
        -   Intermediate post-processing files moved to `post_proc_intermediates` subdirectory.
    -   **Key Actions:**
        -   Loads current and previous damage rasters.
        -   Integrates old detections into current results.
        -   Applies stress period filtering based on combined NDVI/NDWI stress.
        -   Applies shadow masking using QAI data and a specific geographic area (Vipiteno).
        -   Removes single-pixel detections based on area.
        -   Applies all final masks and filters.
        -   Sets color tables and levels for final outputs.
        -   Writes final GeoTIFF and Shapefile outputs with user confirmation.
        -   Moves intermediate files to a dedicated subdirectory.

### Python Package (`fordead_plain`)

The `fordead_plain` package contains the core logic for the damage detection. The R scripts call specific modules within this package.

-   **`fordead_plain/steps/step1_compute_masked_vegetationindex.py`**
    -   **Purpose:** Computes masks and masked vegetation indices for each Sentinel date, filtering by cloudiness and applying source/soil/user-defined masks.
    -   **Inputs:**
        -   Sentinel-2 data directory.
        -   Cloudiness threshold, interpolation order, Sentinel source, mask application flags, soil detection flag, formula for user mask, vegetation index type, compression flag, ignored periods, extent shapefile, path to VI dictionary.
    -   **Outputs:**
        -   Vegetation index data (NetCDF).
        -   Mask data (GeoTIFF).
        -   Soil data (GeoTIFFs for state, first date, count).
        -   Updates `TileInfo` object.
    -   **Key Actions:**
        -   Initializes/imports `TileInfo` object.
        -   Adds parameters to `TileInfo`.
        -   Detects and adds Sentinel data paths to `TileInfo`.
        -   Computes cloudiness percentage.
        -   Imports and resamples Sentinel data.
        -   Computes vegetation index.
        -   Computes and applies various masks (source, soil, user-defined).
        -   Writes vegetation index and mask outputs.
        -   Saves `TileInfo` object.
-   **`fordead_plain/steps/step2_train_model.py`**
    -   **Purpose:** Uses Sentinel dates to train a periodic vegetation index model capable of predicting the vegetation index at any date.
    -   **Inputs:**
        -   Data directory (containing vegetation indices and masks from Step 1).
        -   Minimum number of valid dates for model computation.
        -   Minimum and maximum last dates for training.
        -   Flag for correcting VI using large-scale median VI.
    -   **Outputs:**
        -   Model coefficients (`coeff_model.tif`).
        -   First detection date index (`first_detection_date_index.tif`).
        -   Sufficient coverage mask (`sufficient_coverage_mask.tif`).
        -   Updates `TileInfo` object.
    -   **Key Actions:**
        -   Imports `TileInfo` object.
        -   Adds parameters to `TileInfo`.
        -   Imports stacked masked vegetation indices and masks.
        -   Determines detection dates.
        -   Computes `sufficient_coverage_mask`.
        -   Optionally corrects VI using large-scale median VI.
        -   Models the vegetation index.
        -   Writes model outputs.
        -   Saves `TileInfo` object.
-   **`fordead_plain/steps/step3_dieback_detection.py`**
    -   **Purpose:** Detects anomalies by comparing the vegetation index with its model prediction. Identifies pixels suffering from dieback based on successive anomalies and saves information on stress periods.
    -   **Inputs:**
        -   Data directory (containing model outputs from Step 2).
        -   Anomaly threshold.
        -   Maximum number of stress periods.
        -   Stress index mode.
        -   Vegetation index type and path to VI dictionary (if step1 was skipped).
    -   **Outputs:**
        -   Anomaly data (GeoTIFFs for each date).
        -   Dieback data (GeoTIFFs for state, first date, unconfirmed first date, count).
        -   Stress data (GeoTIFFs for dates, number of periods, cumulative difference, number of dates, stress index).
        -   `too_many_stress_periods_mask.tif`.
        -   Updates `TileInfo` object.
    -   **Key Actions:**
        -   Imports `TileInfo` object.
        -   Adds parameters to `TileInfo`.
        -   Imports necessary data (first detection date index, model coefficients, dieback/stress data).
        -   Imports masked vegetation indices for new dates.
        -   Predicts vegetation index using the model.
        -   Detects anomalies.
        -   Detects dieback based on anomalies.
        -   Saves stress period information.
        -   Writes anomaly, dieback, and stress data outputs.
        -   Saves `TileInfo` object.

## 3. Data Flow

The data flows through the application in a sequential manner:

1.  **Raw Satellite Data (Sentinel-2 BOA & QAI):** Stored in `in_path` (e.g., `/mnt/CEPH_PROJECTS/sao/SENTINEL-2/SentinelVegetationProducts/FORCE/level2/`).
2.  **`B_RUN_simple.R`:** Organizes raw data into grid-specific `_boa` and `_update` directories.
3.  **`1B_fordead_loop_simple.R`:**
    -   Reads organized BOA and QAI data.
    -   Generates R-based QAI masks.
    -   Passes parameters to Python scripts via `param.py`.
    -   Python scripts (Step 1, 2, 3) process data, generating:
        -   Vegetation indices and masks (`fordead_output_*/VegetationIndex`, `fordead_output_*/Mask`).
        -   Model coefficients (`fordead_output_*/DataModel/coeff_model.tif`).
        -   Anomaly, dieback, and stress data (`fordead_output_*/DataAnomalies`, `fordead_output_*/DataDieback`, `fordead_output_*/DataStress`).
4.  **`2B_merge_tiled_outputs.R`:** Reads tiled outputs from the Python steps and mosaics them into full-area raster and vector files (e.g., `output_NDWI_merged.shp`, `output_confidence_NDVI_merged.tif`).
5.  **`3B_index_merging.R`:** Reads merged outputs, applies further masking (NDWI, province, LAFIS), and generates refined monthly and annual products (e.g., `NDVI_merged_masked_with_NDWI_province_lafis.tif`).
6.  **`4-5B_revised.R`:** Reads the refined monthly and annual products, applies styling, and generates initial final GeoTIFFs and Shapefiles (e.g., `changes_yearly_damages_july_2025_25832.tif`).
7.  **`6-7B_revised.R`:** Reads the initial final products, integrates historical data, applies additional filters (stress, shadow, single pixel), and produces the definitive final GeoTIFFs and Shapefiles (e.g., `changes_monthly_damages_july_2025_25832_final.tif`).
8.  **`update_province.R`:** Uses the final GeoTIFFs for `gdaladdo` and `scp` commands to generate overviews and upload data.

## 4. Identified Inconsistencies/Areas for Improvement

During the analysis, several inconsistencies and areas for improvement were noted:

-   **Python Script Locations:** The R scripts (`1B_fordead_loop_simple.R`) refer to Python scripts (`fordead_1.py`, `fordead_2.py`, `fordead_3.py`) located in `/mnt/CEPH_PROJECTS/WALDSCHAEDEN/working_folder/fordead_git/forest_outputs`. However, the core logic appears to reside in `fordead_plain/steps/stepX_*.py` within the provided directory structure. It needs to be clarified which set of Python scripts is the authoritative one for documentation and execution. If the `fordead_1.py` etc. are wrappers, their content should be documented.

-   **Hardcoded Paths:** This is a pervasive issue across multiple R scripts (`3B_index_merging.R`, `4-5B_revised.R`, `6-7B_revised.R`). Many input and output file paths are hardcoded, making the pipeline inflexible and difficult to manage if the directory structure or data sources change. Examples include:

    -   `province.shp`, `lafis.shp`, `base` raster in `3B_index_merging.R`.
    -   `old_monthly` and `old_yearly` paths in `6-7B_revised.R` pointing to a specific date.
    -   QAI input file in `6-7B_revised.R` being hardcoded to a single tile and date.

-   **Redundant Logic:**

    -   The calculation of `last_detectable_date` is duplicated in `3B_index_merging.R` and `4-5B_revised.R`. This logic should be centralized.
    -   QAI processing logic (converting bit values, filtering) is present in both `1B_fordead_loop_simple.R` (R) and `6-7B_revised.R` (R). This redundancy could lead to inconsistencies if not carefully managed.

-   **Parameter Management:** The comment in `6-7B_revised.R` (`# todo: write parameters separately as a file (e.g. outfold keeps reverting to an older setting)`) indicates a known issue with parameter persistence. While Python scripts use `param.py`, a more robust and centralized system for managing parameters across all R scripts would improve maintainability and prevent unexpected behavior.

-   **Fragile Date/Naming Logic:**

    -   The "fix november" and "fix april/winter" adjustments in `3B_index_merging.R` and `4-5B_revised.R` suggest specific workarounds that might indicate underlying issues with date handling or data consistency.
    -   Hardcoded years in output filenames (e.g., `_2025_25832`) in `4-5B_revised.R` and `6-7B_revised.R` require manual updates for each new year.

-   **Implicit Dependencies:** Several scripts rely on variables defined in previously sourced scripts (e.g., `index.list` in `2B_merge_tiled_outputs.R` from `1B_fordead_loop_simple.R`), which can make it harder to understand individual script functionality in isolation.

-   **Hardcoded Visuals:** Color palettes and labels are directly embedded in `4-5B_revised.R`, making them less dynamic and harder to update without modifying the script itself.
