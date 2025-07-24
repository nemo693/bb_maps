# ------------------------------------------------------------------------------
# Script:       01_Import_S2_data.R
#
# Level of interaction: None.
#
# Purpose:      Sets key parameters for the following steps.
#               Defines the tiles to be analyzed, as well as temporal parameters.
#               Iterates through the tiles, calling the core FORDEAD
#               processing script for each grid and detection period.
#
# Date:         26_05_2025
#
# Inputs:
#   - Raw Sentinel-2 BOA data (from `in_path`).
#
# Parameters (fixed):
#   - Configuration parameters (`get_new_images`, `sentinel_only`).
#   - Selection of the tiles (grids) to analyze (`grids`).
#   - Temporal parameters (`train_period_max`, `detection_start_dates`, `detection_end_dates`, `ignored_doy`).
#   - End of the detection period (`detection_end_dates`). Can be left open (now set to end of 2030).
#
# Outputs:
#   - Organized BOA imagery in `_boa` and `_update` directories for each grid.
#   All available observations are put in the _boa folders when this script is
#   run.
#
# Operations:
#   1. Defines global parameters for the FORDEAD process.
#   2. Loops through defined detection periods.
#   3. For each grid, manages BOA image files, moving "future" images to an `_update` directory
#      and restoring "past" images to the main `_boa` directory.
#   4. Sources `02_Execute_Core_FORDEAD_Processing.R` within its main loop to execute the core detection algorithm.
# ------------------------------------------------------------------------------

cat ("Setting paths/n")

# Path to the level2 FORCE folder
in_path = "/mnt/CEPH_PROJECTS/sao/SENTINEL-2/SentinelVegetationProducts/FORCE/level2/" 
# Path to the forest mask (shp)
forest_mask = "/mnt/CEPH_PROJECTS/WALDSCHAEDEN/GIS/land_cover/forest_map/forest_mask/forest_mask_sarah_2019_mmu_1000m.shp"
# Where the BOA exploded data will be stored
out_path = "/mnt/CEPH_PROJECTS/WALDSCHAEDEN/working_folder/outputs/fordead/" # where the boa data will be stored
# Where the analyses will be carried out. Needs to be the same folder where the
# previous updates were elaborated.
for_out_path = "/mnt/CEPH_PROJECTS/WALDSCHAEDEN/working_folder/outputs/fordead_15/"
# Path of the fordead_git folder (where the fordead python scripts are kept)
script_path = "/mnt/CEPH_PROJECTS/WALDSCHAEDEN/working_folder/fordead_git/forest_outputs" # path where the scripts are located


cat("Setting parameters/n")

# --- Configuration Parameters ---

# Flag to determine whether to download new satellite images. Set to TRUE to enable new image acquisition.
get_new_images <- TRUE
# Restrict analysis to Sentinel satellite data only. Currently, the script only supports Sentinel data.
sentinel_only <- TRUE

# --- Analysis Grid Definition ---
# Define spatial grid cells for processing. Each element represents a unique grid ID (X,Y coordinates).

# These grid IDs are used to organize and process data for specific geographic areas.
grids <- c(
  "X0000_Y0003", "X0000_Y0004", "X0000_Y0005",
  "X0001_Y0003", "X0001_Y0004", "X0001_Y0005", 
  "X0002_Y0003", "X0002_Y0004", "X0002_Y0005",
  "X0003_Y0002", "X0003_Y0003", "X0003_Y0004", "X0003_Y0005",
  "X0004_Y0002", "X0004_Y0003", "X0004_Y0004", "X0004_Y0005",
  "X0005_Y0003", "X0005_Y0004"
)

# --- Temporal Parameters ---
# Defines the end date of the training period, used to establish baseline forest conditions.
# This corresponds to the start date of detection period. Do not edit. This parameter must
# not be changed, or fordead will run from scratch.
train_period_min <- c("2019-12-31")
# For areas where the training data is insufficient within the training period,
# train_period_max allows the model to exploit data until this date to train (limited number
# of areas affected).
train_period_max = "2021-11-30"

# Defines the end date for the damage detection period. This also influences
# which images are imported.
# Images after this date are still used for reverting detections, if available.
# This is not relevant if updating the maps. The parameter needs to be set
# anyways.
end_detect <- c("2030-12-31")

# Exclude specific periods when NDVI is naturally low. This prevents false
# detections/reversions during seasonal low-vegetation periods.
# This parameter must not be changed, or fordead will run from scratch.
ignored_doy <- '["11-15","04-15"]'

# Minimum number of training observations required to train the model
min_train_data = 30

# indices to be used
index.list = c("NDVI", "NDWI")



# --- Execute Detection Algorithm ---
# This section sources the core FORDEAD processing script
# (`02_Execute_Core_FORDEAD_Processing.R`), which performs the main damage
# detection logic.

cat("Executing detection algorithm\n")
source("02_Execute_Core_FORDEAD_Processing.R")

# Print a final success message after all detection periods have been processed.
cat("Forest damage detection (fordead part) update completed successfully!\n")

