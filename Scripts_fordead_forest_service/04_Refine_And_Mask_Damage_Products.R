# ------------------------------------------------------------------------------
# Script:       04_Refine_And_Mask_Damage_Products.R
# Purpose:      Refines the merged outputs from the FORDEAD process by applying additional
#               masking steps (only keep agreement areas between NDVI and NDWI,
#               mask out agricultural areas, exclude areas outside of the
#               province), adjusting temporal periods (correct labeling for
#               winter), and deriving annual products. It also includes
#               functionality for area calculations and comparisons.
#
# Inputs:
#   - Merged NDVI and NDWI shapefiles (`output_NDWI_merged.shp`, `output_NDVI_merged.shp`) from `for_out_path`.
#   - Province boundaries shapefile.
#   - LAFIS.
#   - Base Sentinel-2 mosaic GeoTIFF for rasterization.
#
# Outputs:
#   - Masked and refined NDVI shapefiles and GeoTIFFs in `for_out_path`.
#   - Annual and monthly damage GeoTIFF in `for_out_path`.
#
# Operations:
#   1. Loads NDWI and NDVI vector data and rasterizes them.
#   2. Masks the NDVI raster with the NDWI raster to reduce commission errors.
#   3. Applies further spatial masking using province boundaries and LAFIS data.
#   4. Adjusts date labels for specific periods (November, April/winter) and determines the `last_detectable_date` (third to last S2 acquisition).
#   5. Derives annual damage products from the monthly data.
#   6. Writes intermediate and final masked/refined products to disk.
# ------------------------------------------------------------------------------

library(terra)

fold = for_out_path # Set the working folder to the output path from previous steps.

province = vect("/mnt/CEPH_PROJECTS/WALDSCHAEDEN/GIS/forest_administration/forest_inspectorates/ForestInspectorates_polygon.shp") # Load the province boundaries shapefile.
lafis = vect("/mnt/CEPH_PROJECTS/WALDSCHAEDEN/GIS/land_cover/LAFIS_20210413/AGRI_MEASURES_20210413154430_polygon.shp") # Load the LAFIS (agricultural measures) shapefile.
lafis_nf = lafis[!grepl("FO|AL1|AL2|AL3|AL4|AL5|AL6|AL7|AL8|AL9|ANA|PA2|PA3", lafis$CODE),] # Filter out specific agricultural measures from the LAFIS data.
base = rast("/mnt/CEPH_PROJECTS/WALDSCHAEDEN/Products/MOSAICS/SENTINEL2_MOSAIC_4BANDS_LAST/final/SENTINEL2_MOSAIC_20220601_20220930_7L_COUNT.tif")

overwrite = T

######## Get last detectable date

images = list.files(paste0(fold, "/fordead_output_NDWI_X0001_Y0004/VegetationIndex")) # List all image files from a specific NDWI output directory to extract dates.
dates = substr(images, 17, 26) # Extract date strings from the image filenames.
last_detectable_date = dates[(length(dates)-2)] # Determine the last detectable date from the extracted dates.

# Load ndwi & rasterize (used as a mask later)

ndwi = vect(paste0(fold, "output_NDWI_merged.shp")) # Load the merged NDWI shapefile.
# wrong conversion, just need this as a mask
p = unique(ndwi$period) # Extract unique period values from the NDWI data.

fun = function(x) {which(p == x)}

indx = sapply(ndwi$period, FUN = fun) # Apply the function to find the index for each period in the NDWI data.

ndwi$p = indx # Assign the period indices to the NDWI data.

ndwi_rast= rasterize(ndwi, base, field = "p") # Rasterize the NDWI vector data using the base raster and period indices.

gc()

# Load ndvi and rasterize

ndvi = vect(paste0(fold, "output_NDVI_merged.shp")) # Load the merged NDVI shapefile.

p = unique(ndvi$period) # Extract unique period values from the NDVI data.

fun = function(x) {which(p == x)}

indx = sapply(ndvi$period, FUN = fun) # Apply the function to find the index for each period in the NDVI data.

ndvi$p = indx # Assign the period indices to the NDVI data.

base = rast("/mnt/CEPH_PROJECTS/WALDSCHAEDEN/Products/MOSAICS/SENTINEL2_MOSAIC_4BANDS_LAST/final/SENTINEL2_MOSAIC_20220601_20220930_7L_COUNT.tif")

ndvi_rast= rasterize(ndvi, base, field = "p") # Rasterize the NDVI vector data using the base raster and period indices. # Rasterize the NDVI vector data using the base raster and period indices. # Rasterize the NDVI vector data using the base raster and period indices. # Rasterize the NDVI vector data using the base raster and period indices.

gc()

# mask ndvi with ndwi

ndvi_rast_masked = mask(ndvi_rast, ndwi_rast) # Mask the NDVI raster with the NDWI raster to reduce commission errors.

ndvi_masked = as.polygons(ndvi_rast_masked, values = T)

ndvi_masked$first_anomaly = p[ndvi_masked$p]

# Fix the label for November (1-15 Nov, rather than 1-30 Nov, as the doy window ends on Nov 15th)
november = ndvi_masked[grep("11-", ndvi_masked$first_anomaly),]$first_anomaly # Extract November periods from the masked NDVI data.
november_updated = paste0(substr(november, 1, 21), "15") # Update the November period to correct end (15th of the month).
ndvi_masked[grep("11-", ndvi_masked$first_anomaly),]$first_anomaly <- november_updated

# Fix the label for the last period with the overall last detectable date (based on west of province)
last_period = ndvi_masked[grep(substr(last_detectable_date, 1, 7), ndvi_masked$first_anomaly),]$first_anomaly # Extract the last detectable period from the masked NDVI data.
last_period_updated = paste0(substr(last_period, 1, 13), last_detectable_date) # Update the last detectable period with the actual last detectable date.
ndvi_masked[grep(substr(last_detectable_date, 1, 7), ndvi_masked$first_anomaly),]$first_anomaly <- last_period_updated # Update the `first_anomaly` attribute for the last detectable period.

# Fix the label for winter/april
# All that happens between the 15th of November and the 30th of April inevitably falls under the same label.
winter = ndvi_masked[grep("04-", ndvi_masked$first_anomaly),]$first_anomaly # Extract April/winter periods from the masked NDVI data.
winter_updated = paste0((as.numeric(substr(winter, 1, 4)) -1), "-11-16 - ", substr(winter, 1, 4), "-04-30") # Update the April/winter period to span from November 16th of the previous year to April 30th of the current year.
ndvi_masked[grep("11-", ndvi_masked$first_anomaly),]$first_anomaly <- november_updated

writeVector(ndvi_masked, paste0(fold, "NDVI_merged_masked_with_NDWI.shp"), overwrite = overwrite) # Write the masked NDVI data to a shapefile.

# Mask areas outside of the province, or covered by agri areas.

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

# Generate the annual versions
ndvi_masked_annual = ndvi_masked
# Convert monthly/period-based first_anomaly dates to annual labels for yearly product generation.
# FLAG: Fragile logic. Hardcoded year values will require manual updates. Delete after the new logic has been tested.
# ndvi_masked_annual$first_anomaly[grep("2019", ndvi_masked_annual$first_anomaly)] <- "2019"
# ndvi_masked_annual$first_anomaly[grep("2020", ndvi_masked_annual$first_anomaly)] <- "2020"
# ndvi_masked_annual$first_anomaly[grep("2021", ndvi_masked_annual$first_anomaly)] <- "2021"
# ndvi_masked_annual$first_anomaly[grep("2022", ndvi_masked_annual$first_anomaly)] <- "2022"
# ndvi_masked_annual$first_anomaly[grep("2023", ndvi_masked_annual$first_anomaly)] <- "2023"
# ndvi_masked_annual$first_anomaly[grep("2024", ndvi_masked_annual$first_anomaly)] <- "2024"
# ndvi_masked_annual$first_anomaly[grep("2025", ndvi_masked_annual$first_anomaly)] <- "2025"
# Improved logic: Extract 4-digit year from first_anomaly and assign it. 
ndvi_masked_annual$first_anomaly <- stringr::str_extract(ndvi_masked_annual$first_anomaly, "\\b\\d{4}\\b")

# Rasterize the annual masked NDVI data using the base raster and 'first_anomaly' field. 
ndvi_masked_rast_out = rasterize(ndvi_masked_annual, base, field = "first_anomaly") 

writeRaster(
  ndvi_masked_rast_out,
  paste0(fold, "NDVI_merged_masked_with_NDWI_province_only_annual.tif"),
  overwrite = overwrite
) # Write the rasterized annual masked NDVI data to a TIFF file.

# Killed, to delete after verifying that it's not needed in the next run
# Generate yearly .shp
# r = rast(paste0(fold, "NDVI_merged_masked_with_NDWI_province_only_annual.tif")) # Load the annual masked NDVI raster.
# v = as.polygons(r) # Convert the annual raster to polygons.
# writeVector(v, paste0(fold, "bark_beetle_damages_province_annual.shp"), overwrite = overwrite) # Write the annual bark beetle damages to a shapefile.
# 
gc()