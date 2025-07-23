# Set the working directory for the script. This ensures all relative paths are resolved correctly.
setwd("/mnt/CEPH_PROJECTS/WALDSCHAEDEN/Scripts_fordead_forest_service/")

# Define the name for the current update, typically representing the month.
update_name = "july" 
# Define the base output folder for the final products of this run.
outfold <- "/mnt/CEPH_PROJECTS/WALDSCHAEDEN/Products/FORDEAD_01_07_2025/" 

# Before starting, manually make a copy of the folder "fordead_15", just in case something went wrong with the update (~55 mins)
# since the updates always build on the latest one, always keep at least one copy
# the copy is in fordead_15_updates/fordead_15

# --- Backup Previous Run Outputs ---
# This section handles backing up the results from the previous pipeline execution.
# It's crucial for maintaining a history of outputs and for recovery in case of errors.

# List all files (excluding directories) from the previous run's output folder.
# FLAG: The path for previous run files is hardcoded. Consider making this configurable.
prev_run_files = list.files("/mnt/CEPH_PROJECTS/WALDSCHAEDEN/working_folder/outputs/fordead_15", full.names = T)
prev_run_files = prev_run_files[!file.info(prev_run_files)$isdir]

# If there are files from the previous run, move them to a designated backup directory.
if(length(prev_run_files) > 0) {
  # Define the target backup directory. This should ideally be dynamic (e.g., include a timestamp).
  # FLAG: The backup directory path is hardcoded. Consider making this configurable and dynamic (e.g., using `update_name` or current date).
  new_dir <- "/mnt/CEPH_PROJECTS/WALDSCHAEDEN/working_folder/outputs/fordead_15_forest_service_updates/jun_2025/" 
  # Create the new backup directory if it doesn't already exist.
  dir.create(new_dir)
  # Construct the new paths for the files in the backup directory.
  new_paths <- file.path(new_dir, basename(prev_run_files))
  # Move the files from the previous output location to the backup directory.
  file.rename(from = prev_run_files, to = new_paths)
}

# --- Pipeline Execution: Script Calls ---
# This section sequentially sources the R scripts that constitute the main processing pipeline.
# Each script performs a specific stage of the forest damage detection and mapping.
source("01_Import_S2_data.R") # Manages satellite imagery acquisition and organization.
source("03_Mosaic_FORDEAD_Outputs.R") # Merges tiled FORDEAD outputs.
source("04_Refine_And_Mask_Damage_Products.R") # Refines and masks damage products.
source("05_Style_And_Project_Damage_Maps.R") # Prepares and styles damage maps for output.
source("06_Integrate_And_Refine_Damage_Products.R") # Integrates historical data and applies final filters.

# --- Finalization Steps ---

# Generate Overviews (Pyramids) for the final raster products.
# This improves rendering performance in GIS applications.
# `raster_filename_month` and `raster_filename_year` are expected to be defined
# in the last sourced script (06_Integrate_And_Refine_Damage_Products.R).
system(paste0('gdaladdo -r mode ', raster_filename_month, ' 2 4 8 16 32'))
system(paste0('gdaladdo -r mode ', raster_filename_year, ' 2 4 8 16 32'))

# --- Data Upload Commands (Manual Execution) ---
# These commands are provided for manual execution to upload the final products
# to a remote server using SCP (Secure Copy Protocol).
# The password for 'eo_forest' user is required.
paste0("scp ", raster_filename_year," eo_forest@193.106.181.22:changes_yearly_damages_latest.tif")
#password: AqyowLyoNkqRGMTmJK8Z
paste0("scp ", raster_filename_month," eo_forest@193.106.181.22:changes_monthly_damages_latest.tif")








# --- Area Comparison with Historical Data ---
# This section loads historical damage shapefiles and calculates their areas
# for comparison with the current run's results.
library(terra) # Load the 'terra' package for spatial data handling.

# Load vector data for previous runs. Paths are hardcoded.
# FLAG: Hardcoded paths for historical data. Consider making these configurable or deriving them dynamically.
# previous_run <- vect("/mnt/CEPH_PROJECTS/WALDSCHAEDEN/Products/FORDEAD_25_08_2024/changes_yearly_damages_aug_2024_25832.shp")
previous_run <- vect("/mnt/CEPH_PROJECTS/WALDSCHAEDEN/Products/FORDEAD_03_11_2024/changes_yearly_damages_oct_2024_25832.shp")
sept_run     <- vect("/mnt/CEPH_PROJECTS/WALDSCHAEDEN/Products/FORDEAD_30_09_2024/changes_yearly_damages_sept_2024_25832.shp")
current_run  <- vect("/mnt/CEPH_PROJECTS/WALDSCHAEDEN/Products/FORDEAD_09_06_2025/changes_yearly_damages_june_2025_25832_final.shp")

# Calculate and print the area (in hectares) for each loaded shapefile.
expanse(current_run, unit = "ha")
expanse(sept_run, unit = "ha")
expanse(previous_run, unit = "ha")

