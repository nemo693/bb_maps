for (grid_id in grids) {
  
  # Print a message indicating the current grid being processed.
  cat("Processing grid:", grid_id, "\n")
  
  # Define directory paths for the current grid's BOA data.
  # `update_dir` is used to temporarily store images that are outside the
  # current detection period in the paper script. Here, it is just used to bring
  # back images that might have been moved by the paper's script. Kept for
  # safety.
  # `main_boa_dir` is the primary directory for BOA data relevant to the current detection period.
  update_dir <- paste0("/mnt/CEPH_PROJECTS/WALDSCHAEDEN/working_folder/outputs/fordead/", 
                       grid_id, "_boa_update/")
  
  main_boa_dir <- paste0("/mnt/CEPH_PROJECTS/WALDSCHAEDEN/working_folder/outputs/fordead/", 
                         grid_id, "_boa/")
  
  # --- Handle Past Images in boa_update folder (before detection period) ---
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