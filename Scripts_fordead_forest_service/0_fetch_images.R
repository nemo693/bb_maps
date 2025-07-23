# ==============================================================================
# SENTINEL-2 DATA PROCESSING - Unzip new acquisitions and generate FORCE queue
# Author: Emilio Dorigatti | Created: 26_06_2025
# ==============================================================================

# This script automates the process of unzipping new Sentinel-2 acquisitions
# and generating a queue file for the FORCE (Forest Observation and Resource
# Classification Engine) processing system. It identifies newly acquired Sentinel-2
# data within a specified date range, unzips them to a target directory, and then
# creates a queue file that FORCE can use to process these new acquisitions.

# PARAMETERS ===================================================================
# Define the base directory where raw Sentinel-2 L1C data is stored.
basedata_path <- "/mnt/CEPH_BASEDATA/SATELLITE/SENTINEL/SENTINEL-2/L1C/"
# Define the target directory where unzipped Sentinel-2 data will be copied for FORCE processing.
copy_to <- "/mnt/CEPH_PROJECTS/sao/SENTINEL-2/SentinelVegetationProducts/FORCE/level1_test/"

# Define the year for which data is being processed.
year <- "2025"

# DATE RANGE -------------------------------------------------------------------
# Get the start date for data fetching from previously generated FORCE queue files.
# It looks for the latest date in existing queue files and starts from the day after.
prev_queue <- substr(
  list.files(
    "/mnt/CEPH_PROJECTS/sao/SENTINEL-2/SentinelVegetationProducts/FORCE/param/level1/sentinel/",
    pattern = "S2_"), 
  start = 15, stop = 24
)
# Sort dates in descending order and take the first one, then add one day to get the new start date.
from_date <- sort(as.Date(prev_queue, format = "%d_%m_%Y") +1, decreasing = TRUE)[1]

# Define the end date for data fetching as yesterday. This is primarily used for naming the output queue file.
until_date <- as.Date(Sys.time()) - 1

# Construct the full path for the output FORCE queue file.
queue_file_out <- paste0(
  "/mnt/CEPH_PROJECTS/sao/SENTINEL-2/SentinelVegetationProducts/FORCE/param/level1/sentinel/S2_",
  format(from_date, "%d_%m_%Y"), "-", format(until_date, "%d_%m_%Y"), ".txt"
)

# DATA PROCESSING ==============================================================

# Get lists of source Sentinel-2 files, categorized by tile identifier (TPS, TQS, TPT, TQT).
# FLAG: This approach assumes fixed tile identifiers. If new tiles are introduced, this section would need updating.
basedata <- list.files(paste0(basedata_path, year, "/"))
basedata_TPS <- basedata[grep("TPS", basedata)]
basedata_TQS <- basedata[grep("TQS", basedata)]
basedata_TPT <- basedata[grep("TPT", basedata)]
basedata_TQT <- basedata[grep("TQT", basedata)]
# Combine all tile-specific lists into a single list for iteration.
basedata <- list(basedata_TPS, basedata_TPT, basedata_TQS, basedata_TQT)

# Get lists of already processed files in the target FORCE directories, categorized by tile.
# This is used to avoid re-processing already unzipped files.
force_TPS <- list.files(paste0(copy_to, "T32TPS/"), recursive = FALSE, pattern = year)
force_TQS <- list.files(paste0(copy_to, "T32TQS/"), recursive = FALSE, pattern = year)
force_TPT <- list.files(paste0(copy_to, "T32TPT/"), recursive = FALSE, pattern = year)
force_TQT <- list.files(paste0(copy_to, "T32TQT/"), recursive = FALSE, pattern = year)
# Combine all tile-specific lists into a single list for iteration.
force <- list(force_TPS, force_TPT, force_TQS, force_TQT)

# Initialize an empty list to store paths of successfully processed (unzipped) files.
processed_files <- list()

# Loop through each tile's data to identify and unzip new acquisitions.
for(i in 1:4) {
  # Get the base data files and already processed files for the current tile.
  b <- basedata[[i]]
  f <- force[[i]]
  # Extract the tile identifier from the filename of an already processed file.
  # FLAG: Assumes a fixed filename structure for extracting tile ID (characters 39 to 44).
  tile <- substr(f[1], 39, 44) 
  
  # Compare filenames without extensions to identify files that need to be copied.
  # FLAG: Assumes fixed length for filename and extension for both base and force files.
  b_no_ex <- substr(b, 1, (nchar(b[1]) - 4)) # Remove .zip extension
  f_no_ex <- substr(f, 1, (nchar(f[1]) - 5)) # Remove .SAFE extension
  
  # Identify files in the base data that are not yet in the processed (force) directory.
  to_copy_log <- !(b_no_ex %in% f_no_ex)
  b_to_copy <- b[to_copy_log]
  
  # Define the target directory for the current tile and create it if it doesn't exist.
  target <- paste0(copy_to, tile)
  if(!dir.exists(target)) dir.create(target)
  
  # Loop through the files identified for copying and unzip them if they fall within the date range.
  for(n in 1:length(b_to_copy)) {
    curr <- b_to_copy[n]
    # Extract the date from the current filename.
    # FLAG: Assumes a fixed date string position (characters 12 to 19) and format in the filename.
    date <- as.Date(substr(curr, 12, 19), "%Y%m%d")
    # If the file's date is on or after the 'from_date', unzip it.
    if(date >= from_date) {
      # Execute the unzip command. This will extract the contents of the .zip file to the target directory.
      system(paste0("unzip ", basedata_path, year, "/", curr," -d ", target))
      # Add the name of the successfully processed file to the list.
      processed_files[(length(processed_files)+1)] <- curr
    }
  }
}

# GENERATE QUEUE FILE ==========================================================
# This section prepares a queue file for the FORCE processing system.

# List all directories within the 'copy_to' path that end with ".SAFE" (unzipped Sentinel-2 product directories).
dirs <- list.dirs(copy_to, recursive = TRUE)
dirs_safe <- dirs[grep(".SAFE$", dirs)]
# Extract just the base names of these .SAFE directories.
dirs_safe_basename <- basename(dirs_safe)

# Filter the list of .SAFE directories to include only those within the specified date range.
# FLAG: Assumes a fixed date string position (characters 46 to 53) and format within the .SAFE directory name.
dirs_safe_basename <- dirs_safe_basename[as.Date(substr(dirs_safe_basename, 46, 53), "%Y%m%d") >= from_date]

# Create queue entries in the format required by FORCE: "/data/level1/TILE/FILENAME.SAFE QUEUED".
# FLAG: Assumes a fixed tile ID position (characters 39 to 44) in the .SAFE directory name.
queue <- paste0("/data/level1_test/", substr(dirs_safe_basename, 39, 44), "/", dirs_safe_basename, " QUEUED")

# Write the generated queue entries to the output queue file.
write.table(queue, queue_file_out, quote = FALSE, col.names = FALSE, row.names = FALSE)

# Print a summary of the processing.
cat("Processed", length(processed_files), "files\n")
cat("Queue file:", queue_file_out, "with", length(queue), "entries\n")

