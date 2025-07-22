# ------------------------------------------------------------------------------
# Script:       02_Execute_Core_FORDEAD_Processing.R
# Purpose:      This script acts as the bridge between the R orchestration and the Python FORDEAD core.
#               It handles initial data preprocessing, R-based QAI mask generation, and prepares
#               and executes the Python FORDEAD scripts for each grid cell and vegetation index.
# Date:         (Original Date from script, e.g., 26_05_2025)
#
# Inputs:
#   - Parameters from 01_Import_S2_data.R (e.g., `grids`, `train_period_min`, `train_period_max`, `ignored_doy`).
#   - Raw Sentinel-2 BOA and QAI TIFFs from `in_path`.
#   - Forest mask shapefile (`forest_mask`).
#
# Outputs:
#   - Organized BOA data and generated QAI masks in `_boa` directories.
#   - Parameter files (`param.py`) for Python scripts.
#   - Intermediate mask files (`mask_nodata.tif`, `valid_obs.tif`).
#   - Outputs from Python FORDEAD scripts (e.g., `coeff_model.tif`, `sufficient_coverage_mask.tif`).
#
# Operations:
#   1. Defines numerous settings for the FORDEAD process (paths, thresholds, indices).
#   2. Manages the `qai_version` and generates/activates QAI masks based on Sentinel-2 QAI data.
#   3. Prepares a `param.py` file with various parameters for the Python scripts.
#   4. Calls external Python scripts (`fordead_1.py`, `fordead_2.py`, `fordead_3.py`) using `system()` calls.
#   5. Performs a "dummy mask" computation to ensure algorithm stability.
# ------------------------------------------------------------------------------
# to add: check if the tile is already processed with the latest data (skip)

##### DATA PREPROCESSING NEEDS TO BE SEPARATED AS IT IS NOW IN COMMON WITH THE OTHER SCRIPT#####
# for duplicates, the first image is taken (e.g. among sensors ABC)

# install.packages("remotes")
# remotes::install_github("nathan-russell/hashmap")
# remotes::install_github("jkruppa/dataTools")

library(terra)
library(reticulate)
#library(dataTools)

####### SETTINGS ########

run.fordead = T
get_new_images = T # needed on if forcing mask computation
overwrite = T
combine_fmask = F
fmask_py_only = F

conda = "/opt/conda/bin/conda"
myenv_python = "/opt/conda/envs/fordead_plain/bin/python "
myenv_bin = "/opt/conda/envs/fordead_plain/bin"

# path to the level2 coreg tif folder which contains the boa and qai files
in_path = "/mnt/CEPH_PROJECTS/sao/SENTINEL-2/SentinelVegetationProducts/FORCE/level2/" 
forest_mask = "/mnt/CEPH_PROJECTS/WALDSCHAEDEN/GIS/land_cover/forest_map/forest_mask/forest_mask_sarah_2019_mmu_1000m.shp"
# Forest mask options (choose appropriate mask file)
# Option 1: Sarah's 2019 forest mask with 1000 sqm mmu
# forest_mask <- "/mnt/CEPH_PROJECTS/WALDSCHAEDEN/GIS/land_cover/forest_map/forest_mask/forest_mask_sarah_2019_mmu_1000m.shp"
# Option 2: Sarah's 2018 forest mask with 1000 sqm mmu  
# forest_mask <- "/mnt/CEPH_PROJECTS/WALDSCHAEDEN/working_folder/outputs/fordead_17/formask_2018_1000sqm_mmu.shp"


# output path
out_path = "/mnt/CEPH_PROJECTS/WALDSCHAEDEN/working_folder/outputs/fordead/" # where the boa data will be stored
for_out_path = "/mnt/CEPH_PROJECTS/WALDSCHAEDEN/working_folder/outputs/fordead_15/" # where the analyses output will be stored / path of the analysis to update

# path of the fordead_git folder
script_path = "/mnt/CEPH_PROJECTS/WALDSCHAEDEN/working_folder/fordead_git/forest_outputs" # path where the scripts are located

# paths of the fordead scripts - derive from previous
fordead_path_1 = paste0(script_path, "/fordead_1.py")# 
fordead_path_2 = paste0(script_path, "/fordead_2.py")# 
fordead_path_3 = paste0(script_path, "/fordead_3.py")# 

# ignored time frame
#ignored_doy = '["11-15","04-15"]' # defined in calling script


# threshold of anomaly, 0.16 default

thr_anomaly = 0.16 # 

# custom formula for additional local masking

formula = "(B2 > 600) | (B3 <= 0) | (B4 <=0) | (B4 > 1000)| (((B8-B4)/(B8+B4)) > 1) | (((B8-B4)/(B8+B4)) < -1)"
#formula = "(((B8-B4)/(B8+B4)) > 1) | (((B8-B4)/(B8+B4)) < -1) | (((B8A-B11)/(B8A+B11)) > 1) | (((B8A-B11)/(B8A+B11)) < -1)"

# minimum number of training observations (overrides the parameter of fordead?)

min_train_data = 30


# AOI (does not seem to work)

area = "/mnt/CEPH_PROJECTS/WALDSCHAEDEN/GIS/boundaries/province/province_from_forest_inspectorates_b1000_32632.shp"


# Training period

# train_period_min = "2018-12-31"
train_period_min = "2019-12-31"
train_period_max = "2021-11-30"

# indices to be used

index.list = c("NDVI", "NDWI")
#index.list = c("NDVI")
#index.list = c("NDWI")
#index.list = c("brightness")


# soil detection (not working well with our data - can be tuned to extract soil, but bb detection at the same time is not working fine)

soil = "True"


# prov mask - UNUSED

prov_mask = "True"


# detrend - last time I checked it was not working

general_trend_det = "False"


# grids to be analysed

# grids <- c(
#   "X0000_Y0003",
#   "X0000_Y0004",
#   "X0000_Y0005",
#   "X0001_Y0003",
#   "X0001_Y0004",
#   "X0001_Y0005",
#   "X0002_Y0003",
#   "X0002_Y0004",
#   "X0002_Y0005",
#   "X0003_Y0002",
#   "X0003_Y0003",
#   "X0003_Y0004",
#   "X0003_Y0005",
#   "X0004_Y0002",
#   "X0004_Y0003",
#   "X0004_Y0004",
#   "X0004_Y0005",
#   "X0005_Y0003",
#   "X0005_Y0004"
# )


#grids = c("X0001_Y0004", "X0000_Y0004", "X0004_Y0004", "X0004_Y0003")
#grids = c("X0004_Y0005")

# months to be considered (for fetching new data and copying it to the boa folders - only once)

months = c("01","02","03","04", "05", "06", "07", "08", "09", "10", "11","12")


# how to label the changes (monthly, biweekly, every acquisition)

out.freq = "M" # e.g. "sentinel" for all acquisition dates. Other options are e.g.: 'M' (every month), '3M' (every three months), '15D' (every 15 days)"


# output a single file for every period? - not working last time I checked

mult.files = "False"


# get start time of analysis

mask_time = Sys.time()


# qai settings

# QAI VERSION SETTINGS STILL NEED TO BE CHANGED MANUALLY WHEN A NEW VERSION IS ADDED, IN THE MIDDLE OF THE SCRIPT
qai_version = "scx" # snow cloud shadow (filters for snow and cloud, not cloud cover)
# THIS NEEDS STORING OF THE INFO AND A DICTIONARY

# force the generation of the masks

force_masks = F

# illumination state bits for mask computation
# 00 - good (incidence angle < 55°, best quality for top. correction)
# 01 - medium (incidence angle 55°–80°, good quality for top. correction)
# 10 - poor (incidence angle > 80°, low quality for top. correction)
# 11 - shadow (incidence angle > 90°, no top. correction applied)

i_state_bits = NA

sentinel_only = T


####### RUN ########

# for every grid
for(g in grids) { #what is this?? shouldn't it just be tiles?
  
  gc()
  
  updating = F
  
  print(paste0("processing grid ", g))
  
  #### if boa folder !exists ####
  
  # function for renaming bands of a raster according to spectral bands and dates
  myf = function(boa) {
    #print(boa)
    
    img = rast(boa)
    
    dd = strsplit(boa, "/")
    dd = as.data.frame(dd)
    dd = dd[nrow(dd),] # get the last element
    dd = substr(dd, 1, nchar(dd)-4)
    names(img) <- paste0(rep(dd, each = 10), c("_B02", "_B03", "_B04", "_B05", "_B06", "_B07", "_B08", "_B8A", "_B11", "_B12")
    )
    return(img)
  }
  
  #### PREPROCESSING ####
  
  ### FIRST TIME ####
  
  if(!dir.exists(paste0(out_path, g, "_boa"))) {
    
    print("directory does not exist") # if the directory does not exist, return an error
    stop()
    
  }    else { # UPDATE: if the directory exists, check that all the images have been copied there, add the missing ones
    
    updating = T # set flag
    
    # list the original BOA files from SEN2
    or_f = list.files(paste0(in_path, g), full.names = T)
    or_f = or_f[grep("BOA.tif$", or_f)]
    if(sentinel_only) {
      or_f = or_f[grep("SEN2", or_f)]  
    } 
    
    # load images for the current tile and rename the bands according to spectral bands and dates
    img = myf(boa = or_f)
    
    # get the list of dates from the names of the bands
    dates = unique(substr(names(img), 1, 8))
    dates = dates[substr(dates, 5,6) %in% months]
    
    # get the subfolders of the data already present in the boa output folder - pattern to avoid counting non-BOA folder(s)
    p_f = list.files(paste0(out_path, g, "_boa"), full.names = T, pattern = "20")
    
    ##### get new images #####
    
    if((length(p_f)-1) < length(or_f) & get_new_images) { # if the images in the source folder are more than those that were processed, get the new ones
      
      print(paste0("Adding missing elements. Tile: ", g))
      
      # update the mask qai version file - needed here?
      saveRDS(qai_version, paste0(script_path, "/qai_version"))
      
      # create temporary directory for storing the exploded images
      out_img = paste0(out_path, "tmp/")
      file.remove(list.files(out_img, full.names = T, recursive = T), recursive = T)
      file.remove(list.dirs(out_img, full.names = T, recursive = T), recursive = T)
      file.remove(list.dirs(out_img, full.names = T, recursive = T), recursive = T)
      dir.create(out_img)
      
      # get the paths to the original boa files
      f = list.files(paste0(in_path, g), full.names = T)
      f = f[grep("BOA.tif$", f)]
      f = f[grep("SEN2", f)]
      
      # rename bands according to spectral bands and dates
      img = myf(boa = f)
      
      # get the list of dates from the paths
      # THE DATES COULD PROBABLY BE EXTRACTED AT THE VERY BEGINNING, ONCE FOR ALL (in source and in target)
      dates = unique(substr(names(img), 1, 8))
      dates = dates[substr(dates, 5,6) %in% months]
      
      ##### boa ####
      # get only sentinel 2 boa files and transform them according to the needed folder structure in the /tmp folder #
      
      # Run through the dates and the bands. Create a folder for each date and 
      # write each band separately
      
      # get the qai paths
      fm = list.files(paste0(in_path, g), full.names = T)
      fm = fm[grep("QAI.tif$", fm)]
      fm = fm[grep("SEN2", fm)]
      
      # if there is a mismatch, quit. (errors were present at some point, this
      # raises the issue)
      if(length(fm) != length(f)) {
        print("Number of masks and images not matching. Stopping.")
        stop()
      }
      
      # process images
      for(d in dates) {
        # if a folder for the given date (assumed to contain the processed boa and
        # qai) does not exist yet, create it and copy there the data
        if(! d %in% basename(p_f)) {
          
          file.remove(paste0(out_img, d))
          dir.create(paste0(out_img, d))
          
          for(i in 1:10) {
            r = img[[grep(d, names(img))[i]]]
            writeRaster(r, paste0(out_img, d, "/", names(img)[grep(d, names(img))[i]], ".tif"), overwrite = F)
          }
          print(paste0("image ", d))
        }
      }     
      
      print("skipping masks")
      
      ##### !!! masks to be defined uniquely in a function/two #####
      
      # get the paths
      f = list.files(paste0(in_path, g), full.names = T)
      f = f[grep("QAI.tif$", f)]
      f = f[grep("SEN2", f)]
      
      # create the EMPTY (for now) masks folder
      for(d in dates[! dates %in% basename(p_f)]) {
        mdir = paste0(out_img, d, "/masks")
        file.remove(mdir)
        dir.create(mdir)
      }
      
      # copy to the boa folder
      file.copy(from = list.dirs(out_img, full.names = T)[-1], 
                to = paste0(out_path, g, "_boa/"), 
                recursive = T)
      # end of update
    }
    
    #### QAI VERSION CHECK ####
    
    # if the folder exists, rename masks according to the version id (save a flag)
    
    # get the dates - to revise and make faster, for now just copied from above
    print(g)
    
    # get the paths to the boa files - ALSO HERE, TO MAKE ONCE AT THE BEGINNING
    f = list.files(paste0(in_path, g), full.names = T)
    f = f[grep("BOA.tif$", f)]
    f = f[grep("SEN2", f)]
    
    # rename according to bands and dates
    img = myf(boa = f)
    
    # get the list of dates from the paths
    dates = unique(substr(names(img), 1, 8))
    dates = dates[substr(dates, 5,6) %in% months]
    
    # read the last mask version used to create masks
    curr_qai_version = readRDS(paste0(script_path, "/qai_version")) 
    
    print(paste0("current qai version is: ", curr_qai_version))
    
    # if a new qai version is detected, the corresponding masks are generated or
    # activated and the previous ones are renamed and inactivated
    if(curr_qai_version != qai_version | updating == T | force_masks == T) {
      
      print("checkpoint 1")
      
      # select the qai files in the source folder - ALSO HERE, THIS IS A DUPLICATE
      f = list.files(paste0(in_path, g), full.names = T)
      f = f[grep("QAI.tif$", f)]
      f = f[grep("SEN2", f)]
      
      # destination folder for current tile
      boa_d = paste0(out_path, g, "_boa")
      
      # get/update the masks
      if(get_new_images == T) {
        
        if(force_masks == F) {
          
          avail_masks = list.files(boa_d, recursive = T, full.names = T, pattern = qai_version)
          active_masks = list.files(boa_d, recursive = T, full.names = T, pattern = "CLM")
          
          avail_masks_dates = substr(basename(avail_masks), 5, 12)
          missing_masks = dates[!(dates %in% avail_masks_dates)]
          
          if(!(length(grep(qai_version, active_masks)) == length(dates))) { # if the active masks match the ones to be activated, do nothing
            
            # inactivate active masks not matching the target ones if any
            if(!identical(active_masks, character(0))) {
              file.rename(from = active_masks[!grepl(qai_version, active_masks)], 
                          to = gsub(pattern = "CLM", "INA", x = active_masks[!grepl(qai_version, active_masks)]))
            }
            
            
            print("bulk activation")
            # activate target masks
            
            new_names = gsub(pattern = "INA", replacement = "CLM", x = avail_masks)
            
            file.rename(from = avail_masks, 
                        to = new_names)
          } 
          
        } else dates_sub = dates
                
        print("missing masks:")
        print(missing_masks)
        
        for(d in missing_masks) {
          
          # inactivate previous mask(s) (CLM to INA)
          bbbb = paste0(boa_d, "/", d, "/masks/")
          
          # if a mask already exist for this date, deactivate it
          if (length(list.files(bbbb, pattern = "_CLM_")) != 0) {
            file.rename(from = paste0(bbbb, list.files(bbbb, pattern = "_CLM_")), 
                        to = paste0(bbbb, substr(list.files(bbbb, pattern = "_CLM_"), 1, 13),
                                    "INA", 
                                    substr(list.files(bbbb, pattern = "_CLM_"), 17, 25)))
          }
          
          # if a mask matching the desired qai version exists, activate it
          if (length(list.files(bbbb, pattern = qai_version)) != 0 & force_masks == F) {
            file.rename(from = paste0(bbbb, 
                                      list.files(bbbb, pattern = qai_version)), 
                        to = paste0(bbbb, substr(list.files(bbbb, pattern = qai_version), 1, 13),
                                    "CLM", 
                                    substr(list.files(bbbb, pattern = qai_version), 17, 25)))
            
          } else {
            
            file.remove(list.files(bbbb, pattern = qai_version, full.names = T), recursive = T)
            
            if(!fmask_py_only) {
              
              # if the desired mask does not exist, then create it
              mask = rast(f[grep(d, f)][1]) # in case of 2 observations on the same day, take the first one (the same is done above, when writing only the first 10 bands)
              mask_out = mask
              
              # get the integer values that occur in your image (this may not be complete
              # though and you may have to find a way to include all missing combinations if you want to automatize this)
              qai_values = freq(mask)[,2]
              
              qai_bit_values = as.data.frame(qai_values)
              
              # this is an ugly line of code. However, it should convert the integers in bit "strings"
              # so that you can interprete them according to the table in the force website
              for(qai in 1:nrow(qai_bit_values)) qai_bit_values$bit[qai] = paste(sapply(strsplit(paste(rev(intToBits(qai_bit_values$qai_values[qai]))),""),`[[`,2),collapse="")
              #for(i in 1:nrow(qai_bit_values)) qai_bit_values$bit[i] = substr(paste(sapply(strsplit(paste(rev(intToBits(qai_bit_values$qai_values[i]))),""),`[[`,2),collapse=""), 17, 32)
              # select all where last 4 positions are 0 (https://force-eo.readthedocs.io/en/latest/howto/qai.html?highlight=qai)
              # valid data = 0 / cloud state = 00 / cloud shadow flag = 0 or 1 / snow flag = 0
              
              # select the numbers with 0000 bit value at the beginning
              qai_bit_values_cloudfree = qai_bit_values[grepl("01000$|00000$", qai_bit_values$bit),] # select here the valid pixel values
              qai_bit_values_cloudfree$bit = substr(qai_bit_values_cloudfree$bit, nchar(qai_bit_values_cloudfree$bit)-15, nchar(qai_bit_values_cloudfree$bit))
              
              ####### illumination state bits ######
              if(!is.na(i_state_bits)) {
                qai_bit_values_cloudfree = qai_bit_values_cloudfree[grepl(i_state_bits, substr(qai_bit_values_cloudfree$bit, 4, 5)), ]
              }
              
              #qai_bit_values_cloudfree_snow = qai_bit_values[grepl("10000$", qai_bit_values$bit),] # select here the valid pixel values
              
              # write the numbers referring to the bit combination we want into a character 
              no_cloud_bit = as.vector(qai_bit_values_cloudfree$qai_values)
              #no_cloud_bit_snow = as.vector(qai_bit_values_cloudfree_snow$qai_values)
              
              # reclassify QAI layer for clouds and clouds + shadow
              if(length(no_cloud_bit) > 0) mask_out[mask %in% no_cloud_bit] <- 0
              #if(length(no_cloud_bit_snow) > 0) mask_out[mask %in% no_cloud_bit_snow] <- 1
              #mask_out[mask_out > 1] <- 1 # for theia, all pixels with values != from 0 are considered as cloudy
              
              
            }
            
            if(combine_fmask == T) {
              fma = list.files(paste0("/mnt/CEPH_PROJECTS/sao/SENTINEL-2/SentinelVegetationProducts/FORCE/masks/fmask/", g), pattern = d, full.names = T)
              if(!identical(fma, character(0))) {
                fma = rast(fma)
                mask_out[(mask_out == 0) & fma == 2] <- 1 
              }
            }
            
            if(fmask_py_only == T) {
              fma = list.files(paste0("/mnt/CEPH_PROJECTS/sao/SENTINEL-2/SentinelVegetationProducts/FORCE/masks/fmask/", g), pattern = d, full.names = T)
              if(!identical(fma, character(0))) {
                fma = rast(fma)
                mask_out = fma
                mask_out[fma != 2] <- 0 
              } 
            }
            
            names(mask) <- "Mask"
            
            writeRaster(mask_out, paste0(bbbb, qai_version, "_", d, "_CLM_mask.tif"), datatype = "INT2S", overwrite = T)  
            gc()
          }
          print(paste0("mask ", d))
          
        }
      }
      
      saveRDS(qai_version, paste0(script_path, "/qai_version"))
    }
    
  }
  
  #### FORDEAD RUN ####
  
  #break
  
  if(length(list.files(paste0(out_path, g, "_boa"))) != 0  & run.fordead == T) {
    
    for (ii in index.list) {
      
      # generate parameter file
      param = data.frame()
      
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
      
      write.table(param, paste0(dirname(fordead_path_1), "/param.py"), quote = F, row.names = F, col.names = F)
      
      #call fordead
      Sys.getenv()
      Sys.setenv(RETICULATE_PYTHON = myenv_bin)
      #use_condaenv(condaenv = "fordead_plain", conda = conda)
      
      print("running p1")
      # launch the fordead preprocessing
      system(paste0(myenv_python, " ",fordead_path_1))
      
      # save the parameter file
      write.table(param, paste0(for_out_path, "fordead_output_", ii, "_", g, "/param.py"), quote = F, row.names = F, col.names = F)
      
      fordead_time_1 = Sys.time()
      
      print("Time elapsed for 1st fordead script")
      print(mask_time - fordead_time_1)
      
      ##### compute the dummy mask #####
      #(necessary for the algorithm to run - sets to valid the mask for 1 date only for those pixels where all observations are masked)
      # works in combination with a max train date set to the end of the acquisition period only.
      
      # get the masks paths and select the first one
      mask_path = paste0(data_fol, "/Mask")
      mask_files = list.files(mask_path, full.names =  T)
      
      # work only with the masks during the training period
      bn = basename(mask_files)
      bns = substr(bn, 6, 15)
      #train_dates = bns[as.Date(bns) < train_period_min]
      
      # import
      mask = rast(mask_files[as.Date(bns) <= train_period_max])
      
      mask_min = rast(mask_files[as.Date(bns) <= train_period_min])
      
      # fill NAs
      # mask_filled = mask
      # ff = function(i, ...) {i[is.na(i)] <- 0}
      # mask_filled = ff(mask_filled)
      # mask_filled[is.na(mask)] <- 0
      
      smask = sum(mask, na.rm = T)
      train_dates_number = abs((nlyr(mask) - smask))
      
      # get the locations where all observations are masked (==1)
      mask_nodata = smask >= (nlyr(mask) - min_train_data) # 5 is the number of min obs
      
      # facilitate reruns (masks are edited at the first one, needing a rerun)
      if(file.exists(paste0(data_fol, "/mask_nodata.tif"))) {
        print("reading existing mask_nodata.tif file")
        mask_nodata = rast(paste0(data_fol, "/mask_nodata.tif"))
      } else {
        print("writing mask_nodata.tif file")
        writeRaster(train_dates_number , paste0(data_fol, "/valid_obs.tif"), overwrite = overwrite)
        writeRaster(mask_nodata, paste0(data_fol, "/mask_nodata.tif"), overwrite = overwrite)
        # edit the mask for the first date and overwrite it
        # in combination with lowered min obs seems to solve the singular matrix error.
        
        for(i in mask_files[1:min_train_data]) {
          mask_1 = i
          edited_mask = rast(mask_1)
          edited_mask[mask_nodata == 1] <- 0
          edited_mask[is.na(edited_mask)] <- 0
          writeRaster(edited_mask, mask_1, overwrite = T)
        }
      }
      
      # run the second part of fordead
      #
      
      fill = !file.exists(paste0(data_fol, "/DataModel/coeff_model.tif"))
      
      print("running p2")
      system(paste0(myenv_python, " ", fordead_path_2))
      
      # mask out fromt the model the areas that do not match the min number of observations during the max training period (with the file previously generated)
      
      coeff = rast(paste0(data_fol, "/DataModel/coeff_model.tif"))
      
      # fill holes
      if(fill) {
        w = 3
        
        while(w < 7) {
          coeff <- terra::focal(coeff, w = w, fun = "median", na.policy = "only", na.rm = T)
          w <- w + 2
        }
      }
      
      
      coeff[mask_nodata] <- NA
      writeRaster(coeff, paste0(data_fol, "/DataModel/coeff_model.tif"), overwrite)
      scm = rast(paste0(data_fol, "/TimelessMasks/sufficient_coverage_mask.tif"))
      scm[mask_nodata] <- NA
      writeRaster(scm, paste0(data_fol, "/TimelessMasks/sufficient_coverage_mask.tif"), overwrite)
      
      print("running p3")
      system(paste0(myenv_python, " ", fordead_path_3))
      
      gc()
    }
  }
  ########
  gc()
}
# 
# # plots
# fordead_path_4 = "/mnt/CEPH_PROJECTS/WALDSCHAEDEN/working_folder/fordead_git/forest_outputs/fordead_4_fix.py"
# # 
# system(paste0(myenv_python, " ", fordead_path_4))
# system(paste0(myenv_python, " ", "/mnt/CEPH_PROJECTS/WALDSCHAEDEN/working_folder/fordead_git/forest_outputs/fordead_4_fix.py"))
# 
# should probably loop through indices



