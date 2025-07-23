# ------------------------------------------------------------------------------
# Script:       04_Refine_And_Mask_Damage_Products.R
# Purpose:      Refines the merged outputs from the FORDEAD process by applying additional
#               masking steps, adjusting temporal periods, and deriving annual products.
#               It also includes functionality for area calculations and comparisons.
# Date:         (Original Date from script, if available)
#
# Inputs:
#   - Merged NDVI and NDWI shapefiles (`output_NDWI_merged.shp`, `output_NDVI_merged.shp`) from `for_out_path`.
#   - Province boundaries shapefile (`/mnt/CEPH_PROJECTS/WALDSCHAEDEN/GIS/forest_administration/forest_inspectorates/ForestInspectorates_polygon.shp`).
#   - LAFIS (agricultural measures) shapefile (`/mnt/CEPH_PROJECTS/WALDSCHAEDEN/GIS/land_cover/LAFIS_20210413/AGRI_MEASURES_20210413154430_polygon.shp`).
#   - Base Sentinel-2 mosaic GeoTIFF (`/mnt/CEPH_PROJECTS/WALDSCHAEDEN/Products/MOSAICS/SENTINEL2_MOSAIC_4BANDS_LAST/final/SENTINEL2_MOSAIC_20220601_20220930_7L_COUNT.tif`) for rasterization.
#   - Reference bark beetle damage shapefile (`/mnt/CEPH_PROJECTS/WALDSCHAEDEN/GIS/forest_impacts/bark_beetle_damages/damages_full_province/bark_beetle_damages_province.shp`) for comparison.
#
# Outputs:
#   - Masked and refined NDVI shapefiles and GeoTIFFs (e.g., `NDVI_merged_masked_with_NDWI.shp`,
#     `NDVI_merged_masked_with_NDWI_province_lafis.shp`, `NDVI_merged_masked_with_NDWI_province_lafis.tif`) in `for_out_path`.
#   - Annual damage GeoTIFF (`NDVI_merged_masked_with_NDWI_province_only_annual.tif`) and shapefile
#     (`bark_beetle_damages_province_annual.shp`) in `for_out_path`.
#
# Operations:
#   1. Loads NDWI and NDVI vector data and rasterizes them.
#   2. Masks the NDVI raster with the NDWI raster to reduce commission errors.
#   3. Applies further spatial masking using province boundaries and LAFIS data.
#   4. Adjusts date labels for specific periods (November, April/winter) and determines the `last_detectable_date`.
#   5. Derives annual damage products from the monthly data.
#   6. Calculates and compares damage areas with a reference bark beetle damage dataset.
#   7. Writes intermediate and final masked/refined products to disk.
# ------------------------------------------------------------------------------
# This script combines the information from the runs with the two indices 
# The output deriving from the NDVI run is masked based on the output from the
# NDWI run to minimize commission errors.
# 
# LAFIS is also used to perform some extra masking steps.
#
# The periods are adjusted to reflect the winter gap and to consider the last
# possible date of detection.
#
# The yearly product is also derived from the monthly one here.

library(terra)

#fold = "/mnt/CEPH_PROJECTS/WALDSCHAEDEN/working_folder/outputs/fordead_15"
fold = for_out_path # Set the working folder to the output path from previous steps.

province = vect("/mnt/CEPH_PROJECTS/WALDSCHAEDEN/GIS/forest_administration/forest_inspectorates/ForestInspectorates_polygon.shp") # Load the province boundaries shapefile.
lafis = vect("/mnt/CEPH_PROJECTS/WALDSCHAEDEN/GIS/land_cover/LAFIS_20210413/AGRI_MEASURES_20210413154430_polygon.shp") # Load the LAFIS (agricultural measures) shapefile.
lafis_nf = lafis[!grepl("FO|AL1|AL2|AL3|AL4|AL5|AL6|AL7|AL8|AL9|ANA|PA2|PA3", lafis$CODE),] # Filter out specific agricultural measures from the LAFIS data.

overwrite = T

######## get last detection date

images = list.files(paste0(fold, "/fordead_output_NDWI_X0001_Y0004/VegetationIndex")) # List all image files from a specific NDWI output directory to extract dates.

dates = substr(images, 17, 26) # Extract date strings from the image filenames.

last_detectable_date = dates[(length(dates)-2)] # Determine the last detectable date from the extracted dates.

######## load ndwi & rasterize

ndwi = vect(paste0(fold, "output_NDWI_merged.shp")) # Load the merged NDWI shapefile.

# wrong conversion, just need this as a mask

p = unique(ndwi$period) # Extract unique period values from the NDWI data.

#indx = which(p == ndvi$period)

fun = function(x) {which(p == x)}

#sapply(as.data.frame(ndvi$period), FUN = fun)
indx = sapply(ndwi$period, FUN = fun) # Apply the function to find the index for each period in the NDWI data.
#l = lapply(ndvi$period, FUN = fun)

#i = rep(p, length(ndvi$period))

ndwi$p = indx # Assign the period indices to the NDWI data.

base = rast("/mnt/CEPH_PROJECTS/WALDSCHAEDEN/Products/MOSAICS/SENTINEL2_MOSAIC_4BANDS_LAST/final/SENTINEL2_MOSAIC_20220601_20220930_7L_COUNT.tif")

ndwi_rast= rasterize(ndwi, base, field = "p") # Rasterize the NDWI vector data using the base raster and period indices.


gc()
######### load ndvi and rasterize

ndvi = vect(paste0(fold, "output_NDVI_merged.shp")) # Load the merged NDVI shapefile. # Load the merged NDVI shapefile.

p = unique(ndvi$period) # Extract unique period values from the NDVI data. # Extract unique period values from the NDVI data.

#indx = which(p == ndvi$period)

fun = function(x) {which(p == x)}

#sapply(as.data.frame(ndvi$period), FUN = fun)
indx = sapply(ndvi$period, FUN = fun) # Apply the function to find the index for each period in the NDVI data.
#l = lapply(ndvi$period, FUN = fun)

#i = rep(p, length(ndvi$period))

ndvi$p = indx # Assign the period indices to the NDVI data.

base = rast("/mnt/CEPH_PROJECTS/WALDSCHAEDEN/Products/MOSAICS/SENTINEL2_MOSAIC_4BANDS_LAST/final/SENTINEL2_MOSAIC_20220601_20220930_7L_COUNT.tif")

ndvi_rast= rasterize(ndvi, base, field = "p") # Rasterize the NDVI vector data using the base raster and period indices. # Rasterize the NDVI vector data using the base raster and period indices. # Rasterize the NDVI vector data using the base raster and period indices. # Rasterize the NDVI vector data using the base raster and period indices.

gc()

# mask ndvi with ndwi

ndvi_rast_masked = mask(ndvi_rast, ndwi_rast) # Mask the NDVI raster with the NDWI raster to reduce commission errors.

ndvi_masked = as.polygons(ndvi_rast_masked, values = T)

ndvi_masked$first_anomaly = p[ndvi_masked$p]

# fix november
november = ndvi_masked[grep("11-", ndvi_masked$first_anomaly),]$first_anomaly # Extract November periods from the masked NDVI data.
november_updated = paste0(substr(november, 1, 21), "15") # Update the November period to a specific date (15th of the month).
ndvi_masked[grep("11-", ndvi_masked$first_anomaly),]$first_anomaly <- november_updated

# fix the last period with the last detectable date for west of province
last_period = ndvi_masked[grep(substr(last_detectable_date, 1, 7), ndvi_masked$first_anomaly),]$first_anomaly # Extract the last detectable period from the masked NDVI data.
last_period_updated = paste0(substr(last_period, 1, 13), last_detectable_date) # Update the last detectable period with the actual last detectable date.
ndvi_masked[grep(substr(last_detectable_date, 1, 7), ndvi_masked$first_anomaly),]$first_anomaly <- last_period_updated # Update the `first_anomaly` attribute for the last detectable period.

# fix april/winter
winter = ndvi_masked[grep("04-", ndvi_masked$first_anomaly),]$first_anomaly # Extract April/winter periods from the masked NDVI data.
winter_updated = paste0((as.numeric(substr(winter, 1, 4)) -1), "-11-16 - ", substr(winter, 1, 4), "-04-30") # Update the April/winter period to span from November 16th of the previous year to April 30th of the current year.
ndvi_masked[grep("11-", ndvi_masked$first_anomaly),]$first_anomaly <- november_updated

writeVector(ndvi_masked, paste0(fold, "NDVI_merged_masked_with_NDWI.shp"), overwrite = overwrite) # Write the masked NDVI data to a shapefile.

# add masking for the province and lafis

province_rast = rasterize(province, base) # Rasterize the province boundaries.

lafis_rast = rasterize(lafis_nf, base) # Rasterize the filtered LAFIS data.

ndvi_rast_masked = mask(ndvi_rast_masked, province_rast) # Mask the NDVI raster with the province boundaries.

ndvi_rast_masked = mask(ndvi_rast_masked, lafis_rast, inverse = T) # Mask the NDVI raster with the inverse of the LAFIS raster to exclude agricultural areas.

ndvi_masked = as.polygons(ndvi_rast_masked, values = T) # Convert the masked NDVI raster (after province and LAFIS masking) to polygons, retaining the raster values as attributes.

ndvi_masked$first_anomaly = p[ndvi_masked$p] # Assign the corresponding period labels to the `first_anomaly` attribute of the newly masked NDVI data.

writeVector(ndvi_masked, paste0(fold, "NDVI_merged_masked_with_NDWI_province_lafis.shp"), overwrite = overwrite) # Write the NDVI data, masked by province and LAFIS, to a shapefile.

as.integer(expanse(ndvi_masked, unit = "ha")) # Calculate the area of the masked NDVI data in hectares.

ndvi_masked_rast_out = rasterize(ndvi_masked, base, field = "first_anomaly") # Rasterize the masked NDVI data using the base raster and 'first_anomaly' field.

writeRaster(
  ndvi_masked_rast_out,
  paste0(fold, "NDVI_merged_masked_with_NDWI_province_lafis.tif"),
  overwrite = overwrite
) # Write the rasterized, masked NDVI data to a TIFF file.

#annual

ndvi_masked_annual = ndvi_masked
# Convert monthly/period-based first_anomaly dates to annual labels for yearly product generation.
ndvi_masked_annual$first_anomaly[grep("2019", ndvi_masked_annual$first_anomaly)] <- "2019"
ndvi_masked_annual$first_anomaly[grep("2020", ndvi_masked_annual$first_anomaly)] <- "2020"
ndvi_masked_annual$first_anomaly[grep("2021", ndvi_masked_annual$first_anomaly)] <- "2021"
ndvi_masked_annual$first_anomaly[grep("2022", ndvi_masked_annual$first_anomaly)] <- "2022"
ndvi_masked_annual$first_anomaly[grep("2023", ndvi_masked_annual$first_anomaly)] <- "2023"
ndvi_masked_annual$first_anomaly[grep("2024", ndvi_masked_annual$first_anomaly)] <- "2024"
ndvi_masked_annual$first_anomaly[grep("2025", ndvi_masked_annual$first_anomaly)] <- "2025"

#aggregate(ndvi_masked_annual, field = "first_anomaly")

ndvi_masked_rast_out = rasterize(ndvi_masked_annual, base, field = "first_anomaly") # Rasterize the annual masked NDVI data using the base raster and 'first_anomaly' field. 

writeRaster(
  ndvi_masked_rast_out,
  paste0(fold, "NDVI_merged_masked_with_NDWI_province_only_annual.tif"),
  overwrite = overwrite
) # Write the rasterized annual masked NDVI data to a TIFF file.

###

ch_p = vect("/mnt/CEPH_PROJECTS/WALDSCHAEDEN/GIS/forest_impacts/bark_beetle_damages/damages_full_province/bark_beetle_damages_province.shp") # Load the reference bark beetle damages shapefile.
#ch_p = vect("/mnt/CEPH_PROJECTS/WALDSCHAEDEN/GIS/forest_impacts/bark_beetle_damages/damages_full_province/bark_beetle_damages_province_annual.shp")
ch = ndvi_masked # Assign the masked NDVI data to a new variable for area calculations.
  
# area 2020
sum(expanse(ch[grep("2020", ch$first_anomaly),]))/10000 # Calculate the area of detected damages for 2020 in hectares.
sum(expanse(ch_p[grep("2020", ch_p$YEAR_EVENT)]))/10000 # Calculate the area of reference damages for 2020 in hectares.

# area 2021
sum(expanse(ch[grep("2021", ch$first_anomaly),]))/10000 # Calculate the area of detected damages for 2021 in hectares.
sum(expanse(ch_p[grep("2021", ch_p$YEAR_EVENT)]))/10000 # Calculate the area of reference damages for 2021 in hectares.

# area 2022
sum(expanse(ch[grep("2022", ch$first_anomaly),]))/10000 # Calculate the area of detected damages for 2022 in hectares.
sum(expanse(ch_p[grep("2022", ch_p$YEAR_EVENT)]))/10000 # Calculate the area of reference damages for 2022 in hectares.

# area 2023
sum(expanse(ch[grep("2023", ch$first_anomaly),]))/10000 # Calculate the area of detected damages for 2023 in hectares.
sum(expanse(ch_p[grep("2023", ch_p$YEAR_EVENT)]))/10000 # Calculate the area of reference damages for 2023 in hectares.

# area 2024
sum(expanse(ch[grep("2024", ch$first_anomaly),]))/10000 # Calculate the area of detected damages for 2024 in hectares.
sum(expanse(ch_p[grep("2024", ch_p$YEAR_EVENT)]))/10000 # Calculate the area of reference damages for 2024 in hectares.

# area 2025
sum(expanse(ch[grep("2025", ch$first_anomaly),]))/10000 # Calculate the area of detected damages for 2025 in hectares.
sum(expanse(ch_p[grep("2025", ch_p$YEAR_EVENT)]))/10000 # Calculate the area of reference damages for 2025 in hectares.

### yearly shp

r = rast(paste0(fold, "NDVI_merged_masked_with_NDWI_province_only_annual.tif")) # Load the annual masked NDVI raster.

v = as.polygons(r) # Convert the annual raster to polygons.

writeVector(v, paste0(fold, "bark_beetle_damages_province_annual.shp"), overwrite = overwrite) # Write the annual bark beetle damages to a shapefile.

gc()