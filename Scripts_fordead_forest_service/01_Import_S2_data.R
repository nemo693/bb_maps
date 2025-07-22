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

# Data acquisition settings
get_new_images <- TRUE    # Flag to download new satellite images
sentinel_only <- TRUE     # Restrict analysis to Sentinel satellite data only/ CURRENTLY ONLY WORKS IF TRUE

# Output directory for analysis results
# for_out_path_base <- "/mnt/CEPH_PROJECTS/WALDSCHAEDEN/working_folder/outputs/fordead_30_b/"

# --- Analysis Grid Definition ---
# Define spatial grid cells for processing (organized by X,Y coordinates)
grids <- c(
  "X0000_Y0003", "X0000_Y0004", "X0000_Y0005",
  "X0001_Y0003", "X0001_Y0004", "X0001_Y0005", 
  "X0002_Y0003", "X0002_Y0004", "X0002_Y0005",
  "X0003_Y0002", "X0003_Y0003", "X0003_Y0004", "X0003_Y0005",
  "X0004_Y0002", "X0004_Y0003", "X0004_Y0004", "X0004_Y0005",
  "X0005_Y0003", "X0005_Y0004"
)

# --- Temporal Parameters ---
# Training period: used to establish baseline forest conditions
train_period_max <- "2022-11-30"  # End of training period

# Detection periods: define start and end dates for damage detection
detection_start_dates <- c("2019-12-31")  # Start of detection period(s)
detection_end_dates <- c("2025-12-31")    # End of detection period(s) - on its own, this only determines the end of the exported period, if dates after this are available they are still used for reverting detections. In the context of the script it is also used to select the images to import

# Optional: Exclude specific periods when NDVI is naturally low
# This prevents false detections/reversions during seasonal low-vegetation periods
ignored_doy <- '["11-15","04-15"]'  # Day-of-year exclusions (MM-DD format)


# ==============================================================================
# MAIN PROCESSING LOOP
# ==============================================================================

# Process each detection period
for (period_idx in seq_along(detection_start_dates)) {
  
  # Set current period parameters
  train_period_min <- detection_start_dates[period_idx]
  end_detect <- detection_end_dates[period_idx]
  
  cat("Processing detection period:", train_period_min, "to", end_detect, "\n")
  
  # --- File Management for Each Grid Cell ---
  # Organize BOA imagery files based on detection period
  
  for (grid_id in grids) {
    
    cat("Processing grid:", grid_id, "\n")
    
    # Define directory paths for current grid
    update_dir <- paste0("/mnt/CEPH_PROJECTS/WALDSCHAEDEN/working_folder/outputs/fordead/", 
                         grid_id, "_boa_update/")
    main_boa_dir <- paste0("/mnt/CEPH_PROJECTS/WALDSCHAEDEN/working_folder/outputs/fordead/", 
                           grid_id, "_boa/")
    
    # --- Handle Future Images (beyond detection period) ---
    # Move images dated after detection end to update directory
    if (as.Date(end_detect) < as.Date("2025-12-31")) {
      
      # Get all files in main BOA directory
      all_files <- list.files(main_boa_dir, full.names = TRUE)
      
      # Create update directory if it doesn't exist
      if (!dir.exists(update_dir)) {
        dir.create(update_dir, recursive = TRUE)
      }
      
      # Identify files dated after detection end period
      if (length(all_files) > 2) {  # Ensure we have files to process
        file_dates <- as.Date(basename(all_files)[1:(length(all_files) - 2)], 
                              format = "%Y%m%d")
        future_files <- all_files[file_dates > as.Date(end_detect)]
        
        # Move future files to update directory
        if (length(future_files) > 0) {
          success <- file.rename(from = future_files, 
                                 to = paste0(update_dir, basename(future_files)))
          cat("Moved", sum(success), "future files to update directory\n")
        } else {
          cat("No future files to move for grid", grid_id, "\n") # should always go to this case for updating the maps, used to run retrospective analyses
        }
      }
    }
    
    # --- Handle Past Images (before detection period) ---
    # Move images from update directory back to main directory if they belong to detection period
    if (dir.exists(update_dir)) {
      update_files <- list.files(update_dir, full.names = TRUE)
      
      if (length(update_files) > 0) {
        # Check which files have dates within the detection period
        update_file_dates <- as.Date(basename(update_files), format = "%Y%m%d")
        files_to_restore <- update_files[update_file_dates < as.Date(end_detect)]
        
        if (length(files_to_restore) > 0) {
          success <- file.rename(from = files_to_restore,
                                 to = paste0(main_boa_dir, basename(files_to_restore)))
          cat("Restored", sum(success), "files to main BOA directory\n")
        }
      }
    }
  }
  
  # --- Execute Detection Algorithm ---
  # Run the main forest damage detection script with current parameters
  cat("Executing detection algorithm for period", period_idx, "\n")
  source("02_Execute_Core_FORDEAD_Processing.R")
  
  cat("Completed processing for detection period", period_idx, "\n\n")
}

cat("Forest damage detection update completed successfully!\n")

# ==============================================================================
# NOTES:
# - Ensure all required dependencies are installed and accessible
# - Check that file paths exist and are accessible before running
# - Monitor disk space as satellite imagery processing can be storage-intensive
# - Consider running in parallel for large grid sets if computational resources allow
# ==============================================================================