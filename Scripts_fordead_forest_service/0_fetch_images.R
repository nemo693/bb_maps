# ==============================================================================
# SENTINEL-2 DATA PROCESSING - Unzip new acquisitions and generate FORCE queue
# Author: Emilio Dorigatti | Created: 26_06_2025
# ==============================================================================

# PARAMETERS ===================================================================
basedata_path <- "/mnt/CEPH_BASEDATA/SATELLITE/SENTINEL/SENTINEL-2/L1C/"
copy_to <- "/mnt/CEPH_PROJECTS/sao/SENTINEL-2/SentinelVegetationProducts/FORCE/level1_test/"

year <- "2025"

# DATE RANGE -------------------------------------------------------------------
# Get start date from previous queue files
prev_queue <- substr(
  list.files(
    "/mnt/CEPH_PROJECTS/sao/SENTINEL-2/SentinelVegetationProducts/FORCE/param/level1/sentinel/",
    pattern = "S2_"), 
  start = 15, stop = 24
)
from_date <- sort(as.Date(prev_queue, format = "%d_%m_%Y") +1, decreasing = TRUE)[1]

# End date: yesterday. Used only to name the queue file.
until_date <- as.Date(Sys.time()) - 1

# Queue file path
queue_file_out <- paste0(
  "/mnt/CEPH_PROJECTS/sao/SENTINEL-2/SentinelVegetationProducts/FORCE/param/level1/sentinel/S2_",
  format(from_date, "%d_%m_%Y"), "-", format(until_date, "%d_%m_%Y"), ".txt"
)

# DATA PROCESSING ==============================================================

# Get source files by tile
basedata <- list.files(paste0(basedata_path, year, "/"))
basedata_TPS <- basedata[grep("TPS", basedata)]
basedata_TQS <- basedata[grep("TQS", basedata)]
basedata_TPT <- basedata[grep("TPT", basedata)]
basedata_TQT <- basedata[grep("TQT", basedata)]
basedata <- list(basedata_TPS, basedata_TPT, basedata_TQS, basedata_TQT)

# Get existing processed files by tile
force_TPS <- list.files(paste0(copy_to, "T32TPS/"), recursive = FALSE, pattern = year)
force_TQS <- list.files(paste0(copy_to, "T32TQS/"), recursive = FALSE, pattern = year)
force_TPT <- list.files(paste0(copy_to, "T32TPT/"), recursive = FALSE, pattern = year)
force_TQT <- list.files(paste0(copy_to, "T32TQT/"), recursive = FALSE, pattern = year)
force <- list(force_TPS, force_TPT, force_TQS, force_TQT)

# Process each tile - unzip new files
processed_files <- list()

for(i in 1:4) {
  # get files in the 
  b <- basedata[[i]]
  f <- force[[i]]
  tile <- substr(f[1], 39, 44) 
  
  # Compare filenames without extensions
  b_no_ex <- substr(b, 1, (nchar(b[1]) - 4))
  f_no_ex <- substr(f, 1, (nchar(f[1]) - 5))
  
  # Find files to copy
  to_copy_log <- !(b_no_ex %in% f_no_ex)
  b_to_copy <- b[to_copy_log]
  
  target <- paste0(copy_to, tile)
  if(!dir.exists(target)) dir.create(target)
  
  # Unzip files within date range
  for(n in 1:length(b_to_copy)) {
    curr <- b_to_copy[n]
    date <- as.Date(substr(curr, 12, 19), "%Y%m%d")
    if(date >= from_date) {
      system(paste0("unzip ", basedata_path, year, "/", curr," -d ", target))
      processed_files[(length(processed_files)+1)] <- curr
    }
  }
}

# GENERATE QUEUE FILE ==========================================================
# Format: "/data/level1/TILE/FILENAME.SAFE DONE"
dirs <- list.dirs(copy_to, recursive = TRUE)
dirs_safe <- dirs[grep(".SAFE$", dirs)]
dirs_safe_basename <- basename(dirs_safe)

# Filter by date range
dirs_safe_basename <- dirs_safe_basename[as.Date(substr(dirs_safe_basename, 46, 53), "%Y%m%d") >= from_date]

# Create queue entries
queue <- paste0("/data/level1_test/", substr(dirs_safe_basename, 39, 44), "/", dirs_safe_basename, " QUEUED")

# Write queue file
write.table(queue, queue_file_out, quote = FALSE, col.names = FALSE, row.names = FALSE)

cat("Processed", length(processed_files), "files\n")
cat("Queue file:", queue_file_out, "with", length(queue), "entries\n")
