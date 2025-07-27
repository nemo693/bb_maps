# ------------------------------------------------------------------------------
# Script:       02_Execute_Core_FORDEAD_Processing.R
#
# Level of interaction: None.
#
# Purpose:      This script acts as the bridge between the R orchestration and the Python FORDEAD core.
#               It handles initial data preprocessing (S2 BOA data ingestion), R-based QAI mask generation, and prepares
#               and executes the Python FORDEAD scripts for each grid cell and vegetation index.
#
# Inputs:
#   - Parameters from 01_Import_S2_data.R.
#   - Raw Sentinel-2 BOA and QAI TIFFs from `in_path`.
#   - Forest mask shapefile.
#
# Outputs:
#   - Organized BOA data and generated QAI masks in `_boa` directories.
#   - Parameter files (`param.py`) for Python scripts.
#   - Outputs from Python FORDEAD scripts.
#
# Operations:
#   1. Defines several settings for the FORDEAD process (thresholds, python paths).
#   2. Imports BOA files and explodes them (each band in a separate tif file), and imports and processes QAI files into masks (masking cloud cover and snow).
#   3. Prepares a `param.py` file with various parameters for the Python scripts.
#   4. Calls external Python scripts (`fordead_1.py`, `fordead_2.py`, `fordead_3.py`) using `system()` calls.

library(terra)
library(reticulate)

####### SETTINGS ########

# Control flags for script execution and data handling.
run.fordead = T       # Flag to enable/disable the FORDEAD algorithm execution.
get_new_images = T    # Flag to enable/disable fetching new satellite images (needed if forcing mask computation).
overwrite = T         # Flag to enable/disable overwriting existing output files.

# Conda environment and Python executable paths for `reticulate`.
conda = "/opt/conda/bin/conda"
myenv_python = "/opt/conda/envs/fordead_plain/bin/python "
myenv_bin = "/opt/conda/envs/fordead_plain/bin"

# paths of the fordead scripts - derive from previous
fordead_path_1 = paste0(script_path, "/fordead_1.py")# 
fordead_path_2 = paste0(script_path, "/fordead_2.py")# 
fordead_path_3 = paste0(script_path, "/fordead_3.py")# 

# threshold of anomaly, 0.16 default. Don't change this or fordead will run from scratch. 
# deviation of observed index value from modeled index value which determines an anomaly
thr_anomaly = 0.16

# Custom formula for additional local masking. Improves the handling of missed
# clouds and deep shadows in dark forest areas. One of the possible points of improvement.
formula = "(B2 > 600) | (B3 <= 0) | (B4 <=0) | (B4 > 1000)| (((B8-B4)/(B8+B4)) > 1) | (((B8-B4)/(B8+B4)) < -1)"

# AOI (does not seem to work)
area = "/mnt/CEPH_PROJECTS/WALDSCHAEDEN/GIS/boundaries/province/province_from_forest_inspectorates_b1000_32632.shp"

# Soil detection (not working well with our data - can be tuned to extract soil,
# but bb detection at the same time is not working fine) The weakest point of
# the analysis. Some conflicts occur between soil detection and bark beetle
# mapping. Areas flagged as soil are excluded from the analysis. This means no
# detections or reversions are possible in these areas after that date. This was
# working fine but in the latest updates is shows more and more.
soil = "True"

# Province mask
prov_mask = "True"

# Detrend (fordead parameters, needs to be set)
general_trend_det = "False" # leave as is

# Months to be considered (for fetching new data and copying it to the boa folders)
# could be reduced, to reduce the duplication of S2 data. To test. might interfere with data ingestion below.
months = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12") 

# How to label the changes (monthly, biweekly, every acquisition)
out.freq = "M" # e.g. "sentinel" for all acquisition dates. Other options are e.g.: 'M' (every month), '3M' (every three months), '15D' (every 15 days)"
# prints a future warning, should be investigated.

# Output a single file for every period? - not working last time I checked, leave set to False
mult.files = "False"

# Get start time of analysis
mask_time = Sys.time()

# feature to kill
qai_version = "scx" # snow cloud shadow (filters for snow and cloud, not cloud shadow)
# this is a parameter which I kept fixed for the maps. The script would allow to
# generate other versions corresponding to different quality bits combinations
# or to include external masks too (e.g. other fmask implementations as tested
# in the past) Also an interaction with the paper's script. Can be moved on a
# script on its own.

# feature to kill
# Force the generation of the masks. Leave as is.
force_masks = F

# feature to kill
sentinel_only = T

####### RUN ########

# Loop through each grid (tile) for processing.
for(g in grids) {
  
  gc() # Perform garbage collection to free up memory.
  
  updating = F # Initialize a flag to track if the current process is an update.
  
  print(paste0("processing grid ", g)) # Display the current grid being processed for user information.
  
  # Function to rename bands of a raster according to spectral bands and dates.
  # Used throughout the script
  # This is crucial for consistent naming conventions across different BOA images.
  # Reads a file or a stack, and names the bands based on the dates in the paths.
    myf = function(boa) {
      # Load the BOA TIFF file as a raster image.
      img = rast(boa)
      
      # Extract the date from the filename (assuming it's part of the filename).
      dd = strsplit(boa, "/")
      dd = as.data.frame(dd)
      dd = dd[nrow(dd),] # Get the last element (filename).
      dd = substr(dd, 1, nchar(dd)-4) # Remove file extension.
      
      # Rename the raster bands to include the date and band name (e.g., "YYYYMMDD_B02").
      names(img) <- paste0(rep(dd, each = 10), c("_B02", "_B03", "_B04", "_B05", "_B06", "_B07", "_B08", "_B8A", "_B11", "_B12"))
      return(img)
    }
  
  #### PREPROCESSING ####
  
  ### FIRST TIME ####
  
  # Check if the BOA output directory for the current grid exists.
  # If it does not exist, it indicates an initial run for this grid, and the script will stop.
  if(!dir.exists(paste0(out_path, g, "_boa"))) {
    print(paste0("Directory ", out_path, g, "_boa ", "does not exist.")) # Inform the user that the directory is missing.
    stop() # Halt execution as the required input directory is not found.
    
  }    else { 
    
    # If the BOA output directory already exists, this indicates an update scenario.
    updating = T # Set the 'updating' flag to TRUE, indicating an update operation.
    
    # List all original BOA TIFF files from the Sentinel-2 input path for the current grid.
    or_f = list.files(paste0(in_path, g), full.names = T)
    or_f = or_f[grep("BOA.tif$", or_f)]
    # If `sentinel_only` is TRUE, filter the list to include only Sentinel-2 files.
    if(sentinel_only) {
      or_f = or_f[grep("SEN2", or_f)]  
    } 
    
    # Load images for the current tile and rename the bands according to bands and dates
    img = myf(boa = or_f)
    
    # Get the list of dates from the names of the bands
    dates = unique(substr(names(img), 1, 8))
    dates = dates[substr(dates, 5,6) %in% months]
    
    # Get the subfolders already present in the boa output folder - pattern to avoid counting non-BOA folder(s)
    p_f = list.files(paste0(out_path, g, "_boa"), full.names = T, pattern = "20")
    
    
    ##### Get new images from the FORCE folder #####
    
    # Check if there are more original images than already processed ones and if `get_new_images` flag is true.
    if((length(p_f)-1) < length(or_f) & get_new_images) { 
      
      # If there more original images than processed one, then we want to process the ones that are missing
      print(paste0("Adding missing elements. Tile: ", g)) # Inform the user that new elements are being added.
      
      # Update the QAI version file. This might be needed to ensure consistency in mask generation.
      # don't worry about this one
      saveRDS(qai_version, paste0(script_path, "/qai_version"))
      
      # Create a temporary directory for storing exploded (unpacked) images.
      out_img = paste0(out_path, "tmp/")
      # Empty temporary directory from possible leftovers.
      file.remove(list.files(out_img, full.names = T, recursive = T), recursive = T)
      file.remove(list.dirs(out_img, full.names = T, recursive = T), recursive = T)
      file.remove(list.dirs(out_img, full.names = T, recursive = T), recursive = T)
      dir.create(out_img)
      
      # Get the paths to the original Sentinel-2 BOA files for the current grid.
      f = list.files(paste0(in_path, g), full.names = T)
      f = f[grep("BOA.tif$", f)]
      f = f[grep("SEN2", f)]
      
      # Rename bands of the newly fetched images for consistency.
      img = myf(boa = f)
      
      # Extract unique dates from the band names and filter by specified months.
      dates = unique(substr(names(img), 1, 8))
      dates = dates[substr(dates, 5,6) %in% months]
      
      ##### Process BOA Data and QAI Masks #####
      # This section processes Sentinel-2 BOA files and their corresponding QAI masks,
      # organizing them into a structured temporary folder before copying to the main BOA directory.
      
      # Get the paths to the QAI files for the current grid.
      fm = list.files(paste0(in_path, g), full.names = T)
      fm = fm[grep("QAI.tif$", fm)]
      fm = fm[grep("SEN2", fm)]
      
      # Check for a mismatch between the number of QAI files and BOA files.
      # If there's a mismatch, it indicates a data integrity issue, and the script will stop.
      if(length(fm) != length(f)) {
        print("Number of masks and images not matching. Stopping.")
        stop()
      }
      
      # Process each image by date. Add images that are missing from the current list. 
      for(d in dates) {
        # If a folder for the given date (assumed to contain processed BOA and
        # QAI data, meaning the image was already ingested) does not exist,
        # create it and copy the relevant data. Otherwise do nothing and process
        # the next one.
        if(! d %in% basename(p_f)) {
          file.remove(paste0(out_img, d)) # Ensure the temporary directory for this date is clean.
          dir.create(paste0(out_img, d)) # Create the temporary directory for the current date.
          
          # Write each band of the image to a separate TIFF file within the temporary directory.
          for(i in 1:10) {
            r = img[[grep(d, names(img))[i]]]
            writeRaster(r, paste0(out_img, d, "/", names(img)[grep(d, names(img))[i]], ".tif"), overwrite = F)
          }
          print(paste0("Image ", d)) # Indicate which image is being processed.
        }
      }     
      
      print("Processing masks") # Placeholder comment, masks are handled separately.
      
      ##### Masks folder preparation and move exploded BOA to final folders #####
      
      # Re-get the paths to the QAI files for the current grid.
      f = list.files(paste0(in_path, g), full.names = T)
      f = f[grep("QAI.tif$", f)]
      f = f[grep("SEN2", f)]
      
      # Create empty mask folders for dates that haven't been processed yet.
      for(d in dates[! dates %in% basename(p_f)]) {
        mdir = paste0(out_img, d, "/masks")
        file.remove(mdir) # Ensure the mask directory is clean.
        dir.create(mdir) # Create the mask directory.
      }
      
      # Copy the processed temporary BOA and mask folders to the main BOA output folder.
      file.copy(from = list.dirs(out_img, full.names = T)[-1], 
                to = paste0(out_path, g, "_boa/"), 
                recursive = T)
      # End of the update process for new images.
    }
    
    #### QAI ####
    # This section manages the QAI (Quality Assessment Index) masks, ensuring the correct version is active.
    
    # If the BOA output folder exists, this section handles renaming masks according to the version ID.
    # This is part of managing different QAI mask versions and ensuring consistency.
  
    # List original BOA files for the current grid.
    f = list.files(paste0(in_path, g), full.names = T)
    f = f[grep("BOA.tif$", f)]
    f = f[grep("SEN2", f)]
    
    # Load images and rename bands for consistency.
    img = myf(boa = f)
    
    # Extract unique dates from the band names and filter by specified months.
    dates = unique(substr(names(img), 1, 8))
    dates = dates[substr(dates, 5,6) %in% months]
    
    # Read the last QAI mask version used to create masks from a saved RData file.
    curr_qai_version = readRDS(paste0(script_path, "/qai_version")) 

    # If a new QAI version is detected, or if an update/force mask generation is
    # triggered, this block generates or activates the corresponding masks and
    # deactivates previous ones.
    if(curr_qai_version != qai_version | updating == T | force_masks == T) {
      # select the qai files in the source folder - ALSO HERE, THIS IS A DUPLICATE
      f = list.files(paste0(in_path, g), full.names = T)
      f = f[grep("QAI.tif$", f)]
      f = f[grep("SEN2", f)]
      
      # destination folder for current tile
      boa_d = paste0(out_path, g, "_boa")
      
      # get/update the masks
      if(get_new_images == T) {
        # If `force_masks` is FALSE, the script checks for existing masks and activates them if they match the target QAI version.
        if(force_masks == F) {
          avail_masks = list.files(boa_d, recursive = T, full.names = T, pattern = qai_version)
          active_masks = list.files(boa_d, recursive = T, full.names = T, pattern = "CLM")
          avail_masks_dates = substr(basename(avail_masks), 5, 12)
          missing_masks = dates[!(dates %in% avail_masks_dates)]
          
          # If the currently active masks do not match the target QAI version, proceed to update them.
          if(!(length(grep(qai_version, active_masks)) == length(dates))) { # if the active masks match the ones to be activated, do nothing
            
            # inactivate active masks not matching the target ones if any
            # If there are any active masks that do not match the target QAI version, inactivate them.
            if(!identical(active_masks, character(0))) {
              file.rename(from = active_masks[!grepl(qai_version, active_masks)], 
                          to = gsub(pattern = "CLM", "INA", x = active_masks[!grepl(qai_version, active_masks)]))
            }
            
            print("bulk activation")
            # Activate the target masks by renaming them from 'INA' (inactive) to 'CLM' (active).
            new_names = gsub(pattern = "INA", replacement = "CLM", x = avail_masks)
            file.rename(from = avail_masks, 
                        to = new_names)
          } 
          
        } else dates_sub = dates
                
        print("Missing masks:")
        print(missing_masks)
        
        # Iterate through each missing mask date to either deactivate old masks or create new ones.
        for(d in missing_masks) {
          
          # just keep the creation option here
          # get the masks for this date (if any)
          bbbb = paste0(boa_d, "/", d, "/masks/")
          # If a mask already exists for this date, deactivate it by renaming it from 'CLM' (active) to 'INA' (inactive).
          if (length(list.files(bbbb, pattern = "_CLM_")) != 0) {
            file.rename(from = paste0(bbbb, list.files(bbbb, pattern = "_CLM_")), 
                        to = paste0(bbbb, substr(list.files(bbbb, pattern = "_CLM_"), 1, 13),
                                    "INA", 
                                    substr(list.files(bbbb, pattern = "_CLM_"), 17, 25)))
          }
          
          # If a mask matching the desired QAI version exists and `force_masks` is FALSE, activate it.
          if (length(list.files(bbbb, pattern = qai_version)) != 0 & force_masks == F) {
            file.rename(from = paste0(bbbb, 
                                      list.files(bbbb, pattern = qai_version)), 
                        to = paste0(bbbb, substr(list.files(bbbb, pattern = qai_version), 1, 13),
                                    "CLM", 
                                    substr(list.files(bbbb, pattern = qai_version), 17, 25)))
            
          } else { # If the desired mask does not exist, create it.
            file.remove(list.files(bbbb, pattern = qai_version, full.names = T), recursive = T)
            if(!fmask_py_only) {
              # Create the mask from the QAI file. In case of multiple observations on the same day, the first one is selected.
              mask = rast(f[grep(d, f)][1]) # in case of 2 observations on the same day, take the first one (the same is done above, when writing only the first 10 bands)
              mask_out = mask
              
              # Extract unique integer values from the QAI mask.
              qai_values = freq(mask)[,2]
              
              # Generate a dataframe with the unique QAI values observed in the bitmask
              qai_bit_values = as.data.frame(qai_values)
              
              # Convert QAI integer values to 16-bit binary strings for interpretation based on FORCE documentation.
              for(qai in 1:nrow(qai_bit_values)) qai_bit_values$bit[qai] = paste(sapply(strsplit(paste(rev(intToBits(qai_bit_values$qai_values[qai]))),""),`[[`,2),collapse="")
              # select all where last 4 positions are 0 (https://force-eo.readthedocs.io/en/latest/howto/qai.html?highlight=qai)
              # valid data = 0 / cloud state = 00 / cloud shadow flag = 0 or 1 / snow flag = 0
              
              # Select cloud-free pixel values based on specific bit patterns (01000$ or 00000$) as per FORCE QAI documentation.
              qai_bit_values_cloudfree = qai_bit_values[grepl("01000$|00000$", qai_bit_values$bit),] # select here the valid pixel values
              qai_bit_values_cloudfree$bit = substr(qai_bit_values_cloudfree$bit, nchar(qai_bit_values_cloudfree$bit)-15, nchar(qai_bit_values_cloudfree$bit))
              
              # write the numbers referring to the bit combination we want into a character 
              no_cloud_bit = as.vector(qai_bit_values_cloudfree$qai_values)
              
              # Reclassify the QAI layer to mask out clouds and cloud shadows.
              if(length(no_cloud_bit) > 0) mask_out[mask %in% no_cloud_bit] <- 0
            }
            
            names(mask) <- "Mask"
            # Assign a name to the mask layer and write the processed mask to a TIFF file.
            writeRaster(mask_out, paste0(bbbb, qai_version, "_", d, "_CLM_mask.tif"), datatype = "INT2S", overwrite = T)  
            gc()
          }
          print(paste0("mask ", d))
          
        }
      }
      
      saveRDS(qai_version, paste0(script_path, "/qai_version"))
    }
    
  }
  
  # added
  gc()
  
  #### FORDEAD RUN ####
  
  # Initiate FORDEAD processing if the BOA output folder exists and the `run.fordead` flag is TRUE.
  if(length(list.files(paste0(out_path, g, "_boa"))) != 0  & run.fordead == T) {
    
    # Iterate through each vegetation index to generate parameter files and execute the FORDEAD processing.
    for (ii in index.list) {
      
      # Initialize a data frame to store parameters for the FORDEAD Python scripts.
      param = data.frame()
      
      # Populate the parameter data frame with configuration settings for the FORDEAD Python scripts.
      # grid
      param[1,1] = paste0("g =", "'", g, "'")
      
      # index to be used
      param[2,1] = paste0("index = ", "'", ii, "'")
      
      # input  folder
      in_dir = paste0(out_path, g, "_boa")
      param[3,1] = paste0("in_dir = ", "'", in_dir, "'")
      
      # data folder
      data_fol = paste0(for_out_path, "fordead_output_", ii, "_", g)
      param[4,1] = paste0("data_dir = ", "'", data_fol, "'")
      
      # min_train_period
      param[5,1] = paste0("min_last_date_training = ", "'", train_period_min, "'")
      
      # max
      param[6,1] = paste0("max_last_date_training = ", "'", train_period_max, "'")
      
      # formula
      param[7,1] = paste0("formula = ", "'", formula, "'")
      
      # min_train_data
      param[8,1] = paste0("min_dates= ", min_train_data)
      
      # thr anomaly
      param[9,1] = paste0("thr_anomaly= ", thr_anomaly)
      
      # doy_range
      param[10,1] = paste0("ignored_doy= ", ignored_doy)
      
      # area
      param[11,1] = paste0("area= ", "'", area, "'")
      
      # soi detection
      param[12,1] = paste0("soil= ", soil)
      
      # general_trend_detection
      param[13,1] = paste0("trend= ", general_trend_det)
      
      # output frequency
      param[14,1] = paste0("out_freq= ", "'", out.freq, "'")
      
      # multiple files
      param[15,1] = paste0("mult_files= ", mult.files)
      
      # multiple files
      param[16,1] = paste0("p_mask= ", prov_mask)
      
      # multiple files
      param[17,1] = paste0("mask_version= ", "'", qai_version, "'")
      
      # multiple files
      param[18,1] = paste0("end_detect= ", "'", end_detect, "'")
      
      # forest mask file
      param[19,1] = paste0("forest_path= ", "'", forest_mask, "'")
      
      # Write the generated parameters to a `param.py` file, which will be used by the Python FORDEAD scripts.
      write.table(param, paste0(dirname(fordead_path_1), "/param.py"), quote = F, row.names = F, col.names = F)
      
      # Configure the Python environment for `reticulate`.
      Sys.getenv()
      Sys.setenv(RETICULATE_PYTHON = myenv_bin)
      #use_condaenv(condaenv = "fordead_plain", conda = conda)
      
      print("running p1")
      # Execute the first FORDEAD Python script for preprocessing.
      system(paste0(myenv_python, " ",fordead_path_1))
      
      # Save a copy of the parameter file to the output directory for the current FORDEAD run.
      write.table(param, paste0(for_out_path, "fordead_output_", ii, "_", g, "/param.py"), quote = F, row.names = F, col.names = F)
      
      fordead_time_1 = Sys.time()
      
      print("Time elapsed for 1st fordead script")
      print(mask_time - fordead_time_1)
      
      ##### Compute Dummy Mask for Algorithm Stability #####
      # This section computes a "dummy mask" which is necessary for the algorithm to run stably.
      # It sets to valid the mask for a single date for pixels where all observations are otherwise masked.
      # This works in combination with a maximum training date set to the end of the acquisition period only.
      
      # Construct the path to the mask files, retrieve them, and then extract their base names and dates.
      # Only masks within the defined training period are selected for further processing.
      mask_path = paste0(data_fol, "/Mask")
      mask_files = list.files(mask_path, full.names =  T)
      bn = basename(mask_files)
      bns = substr(bn, 6, 15)
      #train_dates = bns[as.Date(bns) < train_period_min]
      
      # Import the masks, filtering them to include only those within the training period.
      mask = rast(mask_files[as.Date(bns) <= train_period_max])
      
      # Import masks specifically for the minimum training period.
      mask_min = rast(mask_files[as.Date(bns) <= train_period_min])
      
      # Calculate the sum of valid (non-NA) observations for each pixel.
      smask = sum(mask, na.rm = T)
      train_dates_number = abs((nlyr(mask) - smask))
      
      # Identify locations where the number of valid observations is below the minimum training data threshold.
      mask_nodata = smask >= (nlyr(mask) - min_train_data) # 5 is the number of min obs
      
      # Facilitate reruns: if `mask_nodata.tif` exists, read it; otherwise, create it along with `valid_obs.tif`.
      if(file.exists(paste0(data_fol, "/mask_nodata.tif"))) {
        print("reading existing mask_nodata.tif file")
        mask_nodata = rast(paste0(data_fol, "/mask_nodata.tif"))
      } else { # If `mask_nodata.tif` does not exist, create and write it along with `valid_obs.tif`.
        print("writing mask_nodata.tif file")
        writeRaster(train_dates_number , paste0(data_fol, "/valid_obs.tif"), overwrite = overwrite)
        writeRaster(mask_nodata, paste0(data_fol, "/mask_nodata.tif"), overwrite = overwrite)
        # Edit the mask for the first date and overwrite it. This helps resolve singular matrix errors.
        for(i in mask_files[1:min_train_data]) {
          mask_1 = i
          edited_mask = rast(mask_1)
          edited_mask[mask_nodata == 1] <- 0
          edited_mask[is.na(edited_mask)] <- 0
          writeRaster(edited_mask, mask_1, overwrite = T)
        }
      }
      
      # Run the second part of fordead

      # Check if the coefficient model file exists to determine whether to fill holes in the data.
      fill = !file.exists(paste0(data_fol, "/DataModel/coeff_model.tif"))
      
      print("running p2")
      # Execute the second FORDEAD Python script.
      system(paste0(myenv_python, " ", fordead_path_2))
      
      # Mask out from the model coefficient raster the areas that do not match
      # the min number of observations during the max training period (with the
      # file previously generated)
      
      coeff = rast(paste0(data_fol, "/DataModel/coeff_model.tif"))
      
      # If `fill` is TRUE, fill holes in the `coeff` raster using a median filter.
      if(fill) {
        w = 3
        
        while(w < 7) {
          coeff <- terra::focal(coeff, w = w, fun = "median", na.policy = "only", na.rm = T)
          w <- w + 2
        }
      }
      
      # Mask out areas from the coefficient model where there are insufficient valid observations.
      coeff[mask_nodata] <- NA
      writeRaster(coeff, paste0(data_fol, "/DataModel/coeff_model.tif"), overwrite)
      # Write the processed coefficient model to a TIFF file.
      # Load the sufficient coverage mask and then mask out areas with insufficient valid observations.
      scm = rast(paste0(data_fol, "/TimelessMasks/sufficient_coverage_mask.tif"))
      scm[mask_nodata] <- NA
      writeRaster(scm, paste0(data_fol, "/TimelessMasks/sufficient_coverage_mask.tif"), overwrite)
      
      print("running p3")
      # Execute the third FORDEAD Python script.
      system(paste0(myenv_python, " ", fordead_path_3))
      
      gc() # Perform garbage collection to free up memory after processing each vegetation index.
    }
  }
  gc() # Perform garbage collection to free up memory after processing each grid.
}