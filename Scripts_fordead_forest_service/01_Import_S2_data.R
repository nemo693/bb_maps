# ------------------------------------------------------------------------------
# Script:       01_Import_S2_data.R
# Purpose:      Manages the acquisition and organization of Sentinel-2 BOA imagery for processing.
#               It defines spatial grids and temporal parameters, and iterates through them,
#               calling the core FORDEAD processing script for each grid and detection period.
# Date:         26_05_2025
#
# Inputs:
#   - Configuration parameters (`get_new_images`, `sentinel_only`).
#   - Spatial grid definitions (`grids`).
#   - Temporal parameters (`train_period_max`, `detection_start_dates`, `detection_end_dates`, `ignored_doy`).
#   - Raw Sentinel-2 BOA data (from `in_path`).
#
# Outputs:
#   - Organized BOA imagery in `_boa` and `_update` directories for each grid.
#
# Operations:
#   1. Defines global parameters for the FORDEAD process.
#   2. Loops through defined detection periods.
#   3. For each grid, manages BOA image files, moving "future" images to an `_update` directory
#      and restoring "past" images to the main `_boa` directory.
#   4. Sources `02_Execute_Core_FORDEAD_Processing.R` within its main loop to execute the core detection algorithm.
# ------------------------------------------------------------------------------

# --- Configuration Parameters ---

# Flag to determine whether to download new satellite images. Set to TRUE to enable new image acquisition.
get_new_images <- TRUE    
# Restrict analysis to Sentinel satellite data only. Currently, the script only supports Sentinel data.
sentinel_only <- TRUE     

# Output directory for analysis results (example commented out).
# for_out_path_base <- "/mnt/CEPH_PROJECTS/WALDSCHAEDEN/working_folder/outputs/fordead_30_b/"

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
train_period_max <- "2022-11-30"  

# Defines the start dates for each damage detection period.
detection_start_dates <- c("2019-12-31")  
# Defines the end dates for each damage detection period. This also influences which images are imported.
# Images after this date are still used for reverting detections, if available.
detection_end_dates <- c("2025-12-31")    

# Optional: Exclude specific periods when NDVI is naturally low.
# This prevents false detections/reversions during seasonal low-vegetation periods.
ignored_doy <- '["11-15","04-15"]'  


# ==============================================================================
# MAIN PROCESSING LOOP
# ==============================================================================

# Iterate through each defined detection period. The loop runs once for each element in `detection_start_dates`.
for (period_idx in seq_along(detection_start_dates)) {
  
  # Set the current training period minimum and detection end date based on the current iteration.
  train_period_min <- detection_start_dates[period_idx]
  end_detect <- detection_end_dates[period_idx]
  
  # Print a message indicating the current detection period being processed.
  cat("Processing detection period:", train_period_min, "to", end_detect, "\n")

  
  # --- File Management for Each Grid Cell ---
  # This section organizes BOA (Bottom-Of-Atmosphere) imagery files for each grid cell.
  # It moves images outside the current detection period to a temporary `_update` directory
  # and restores images from `_update` that fall within the current period.
  
  for (grid_id in grids) {
    
    # Print a message indicating the current grid being processed.
    cat("Processing grid:", grid_id, "\n")
    
    # Define directory paths for the current grid's BOA data.
    # `update_dir` is used to temporarily store images that are outside the current detection period.
    # `main_boa_dir` is the primary directory for BOA data relevant to the current detection period.
    update_dir <- paste0("/mnt/CEPH_PROJECTS/WALDSCHAEDEN/working_folder/outputs/fordead/", 
                         grid_id, "_boa_update/")
    main_boa_dir <- paste0("/mnt/CEPH_PROJECTS/WALDSCHAEDEN/working_folder/outputs/fordead/", 
                           grid_id, "_boa/")
    
    # --- Handle Future Images (beyond detection period) ---
    # This block identifies and moves BOA images that are dated after the current
    # detection period's end date to a temporary update directory (`_boa_update`).
    # This ensures that only relevant images for the current analysis period are processed.
    # The condition `as.Date(end_detect) < as.Date("2025-12-31")` suggests this block is primarily
    # active for periods before the final end date, allowing for incremental updates.
    if (as.Date(end_detect) < as.Date("2025-12-31")) {
      
      # Get all files currently present in the main BOA directory for the grid.
      all_files <- list.files(main_boa_dir, full.names = TRUE)
      
      # Create the update directory if it does not already exist.
      if (!dir.exists(update_dir)) {
        dir.create(update_dir, recursive = TRUE)
      }
      
      # Identify files dated after detection end period
      # The `length(all_files) > 2` condition is a heuristic to ensure there are enough files
      # to safely extract dates without errors (assuming a specific naming convention).
      # FLAG: The `[1:(length(all_files) - 2)]` indexing for `file_dates` might be fragile
      # if the naming convention or number of non-image files changes.
      if (length(all_files) > 2) {  
        file_dates <- as.Date(basename(all_files)[1:(length(all_files) - 2)], 
                              format = "%Y%m%d")
        future_files <- all_files[file_dates > as.Date(end_detect)]
        
        # Move identified future files to the update directory.
        if (length(future_files) > 0) {
          success <- file.rename(from = future_files, 
                                 to = paste0(update_dir, basename(future_files)))
          cat("Moved", sum(success), "future files to update directory\n")
        } else {
          # This branch is expected when no images dated after the detection end are found,
          # which is typically the case for incremental updates where all relevant images
          # for the current period are already in the main BOA directory.
          cat("No future files to move for grid", grid_id, "\n") 
        }
      }
    }
    
    # --- Handle Past Images (before detection period) ---
    # This block identifies and moves images from the `_boa_update` directory back to the
    # main BOA directory (`main_boa_dir`) if their dates fall within the current detection period.
    # This is crucial for restoring images that were previously moved out during a different (later) detection period.
    if (dir.exists(update_dir)) {
      # List all files currently present in the update directory.
      update_files <- list.files(update_dir, full.names = TRUE)
      
      if (length(update_files) > 0) {
        # Extract dates from filenames and identify files whose dates are before the `end_detect` date.
        # These are files that were previously moved out but are now relevant for the current period.
        update_file_dates <- as.Date(basename(update_files), format = "%Y%m%d")
        files_to_restore <- update_files[update_file_dates < as.Date(end_detect)]
        
        # If files need to be restored, move them back to the main BOA directory.
        if (length(files_to_restore) > 0) {
          success <- file.rename(from = files_to_restore,
                                 to = paste0(main_boa_dir, basename(files_to_restore)))
          cat("Restored", sum(success), "files to main BOA directory\n")
        }
      }
    }
  }
  
  # --- Execute Detection Algorithm ---
  # This section sources the core FORDEAD processing script (`02_Execute_Core_FORDEAD_Processing.R`)
  # for the current detection period and grid. This script performs the main damage detection logic.
  cat("Executing detection algorithm for period", period_idx, "\n")
  source("02_Execute_Core_FORDEAD_Processing.R")
  
  # Print a message indicating the completion of processing for the current detection period.
  cat("Completed processing for detection period", period_idx, "\n\n")
}

# Print a final success message after all detection periods have been processed.
cat("Forest damage detection update completed successfully!\n")

