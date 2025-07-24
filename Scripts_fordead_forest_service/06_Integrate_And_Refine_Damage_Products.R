
# ------------------------------------------------------------------------------
# Script:       06_Integrate_And_Refine_Damage_Products.R
# Purpose:      Performs the final post-processing steps, integrating previous results,
#               applying additional filters (stress periods, shadows, single pixels),
#               and generating the definitive output products.
# Date:         (Original Date from script, if available)
#
# Inputs:
#   - Current yearly and monthly damage rasters (from `05_Style_And_Project_Damage_Maps.R`).
#   - Previous monthly and yearly damage rasters (hardcoded paths in the script).
#   - NDVI and NDWI stress period rasters (e.g., `/mnt/CEPH_PROJECTS/WALDSCHAEDEN/working_folder/outputs/fordead_15/output_nb_periods_stress_NDVI_merged.tif`).
#   - QAI data for shadow masking (hardcoded path: `/mnt/CEPH_PROJECTS/sao/SENTINEL-2/SentinelVegetationProducts/FORCE/level2_misc/level2_tiles/T32TPT/20241103_LEVEL2_SEN2A_QAI.tif`).
#   - Forest inspectorates shapefile (`/mnt/CEPH_PROJECTS/WALDSCHAEDEN/GIS/forest_administration/forest_inspectorates/ForestInspectorates_polygon.shp`).
#   - A base raster for projection (hardcoded path: `/mnt/CEPH_PROJECTS/WALDSCHAEDEN/Products/MOSAICS/LANDSAT_MOSAIC_4BANDS_CENTER/final/SENTINEL2_MOSAIC_19840601_19900930_19870715_7C_CIR.tif`).
#
# Outputs:
#   - Final yearly and monthly damage GeoTIFFs and Shapefiles (e.g., `changes_monthly_damages_july_2025_25832_final.tif`).
#   - Intermediate post-processing files moved to `post_proc_intermediates` subdirectory.
#
# Operations:
#   1. Loads current and previous damage rasters.
#   2. Integrates old detections into current results.
#   3. Applies stress period filtering based on combined NDVI/NDWI stress.
#   4. Applies shadow masking using QAI data and a specific geographic area (Vipiteno).
#   5. Removes single-pixel detections based on area.
#   6. Applies all final masks and filters.
#   7. Sets color tables and levels for final outputs.
#   8. Writes final GeoTIFF and Shapefile outputs with user confirmation.
#   9. Moves intermediate files to a dedicated subdirectory.
# ------------------------------------------------------------------------------

# todo: set permanent changes to period with at least 2 months of check
# todo: revise stress period filtering
# todo: write parameters separately as a file (e.g. outfold keeps reverting to an older setting)

# =============================================================================
# CONFIRMATION FUNCTION
# =============================================================================
confirm_write <- function(filename) {
  response <- readline(prompt = paste("Do you want to write or overwrite", filename, "? (yes/no): "))
  tolower(response) == "yes"
}


library(terra)

# Clear temporary files (optional)
# terra::tmpFiles(T, T, T, T)

# ==============================================================================
# CONFIGURATION AND INPUT LOADING
# ==============================================================================

# Load base raster

### CHECK IN THE OTHER SCRITPS AS THIS LINE IS WRONG AND IT WAS COPIED FROM SOMEWHERE ELSE
base = rast("/mnt/CEPH_PROJECTS/WALDSCHAEDEN/Products/MOSAICS/LANDSAT_MOSAIC_4BANDS_CENTER/final/SENTINEL2_MOSAIC_19840601_19900930_19870715_7C_CIR.tif")[[1]]

# Load current detection results
yearly_original = rast(paste0(outfold, prefix, "_", "yearly_damages_", update_name, "_2025_25832.tif"))
yearly = yearly_original
monthly = rast(paste0(outfold, prefix, "_", "monthly_damages_", update_name, "_2025_25832.tif"))

# Load previous detection results for integration
old_monthly = rast("/mnt/CEPH_PROJECTS/WALDSCHAEDEN/Products/FORDEAD_09_06_2025/changes_monthly_damages_june_2025_25832_final.tif") # The latest update that was delivered (using this pipeline). should run at regular intervals and come up with a rule here
old_yearly = rast("/mnt/CEPH_PROJECTS/WALDSCHAEDEN/Products/FORDEAD_09_06_2025/changes_yearly_damages_june_2025_25832_final.tif")

# Initial frequency check
cat("Initial frequencies:\n")
print("Old monthly:")
print(freq(old_monthly))
print("Current monthly:")
print(freq(monthly))
print("Current yearly:")
print(freq(yearly))


# glpat-pvyAi7X6uQqG5yf2_85b
# =============================================================================
# INTEGRATE OLD DETECTIONS
# =============================================================================

cat("\n=== Integrating old detections ===\n")

# Filter old monthly detections (adjust threshold as needed)
# Remove values > 36 (was 28 in original script - adjust based on your needs)
old_monthly_filtered = old_monthly
old_monthly_filtered[old_monthly_filtered > 39] <- NA

# Create condition mask for where old detections should be preserved
condition1 = !is.na(old_monthly_filtered)

# Integrate old detections into current results
monthly[condition1] <- old_monthly_filtered[condition1]
yearly[condition1] <- old_yearly[condition1]

# Calculate difference for monitoring
difference = yearly - yearly_original

cat("Frequencies after integration:\n")
print("Monthly:")
print(freq(monthly))
print("Yearly:")
print(freq(yearly))

# Save intermediate outputs
if (confirm_write(paste0(outfold, prefix, "_monthly_integration.tif"))) {
  writeRaster(monthly, paste0(outfold, prefix, "_monthly_integration.tif"), overwrite=TRUE, datatype = "INT1U")
}
if (confirm_write(paste0(outfold, prefix, "_yearly_integration.tif"))) {
  writeRaster(yearly, paste0(outfold, prefix, "_yearly_integration.tif"), overwrite=TRUE, datatype = "INT1U")
}

freq(yearly)

# =============================================================================
# STRESS PERIOD FILTERING (ADVANCED FILTERING)
# =============================================================================

cat("\n=== Applying stress period filtering ===\n")

# Load stress period masks
ndvi_stress = rast("/mnt/CEPH_PROJECTS/WALDSCHAEDEN/working_folder/outputs/fordead_15/output_nb_periods_stress_NDVI_merged.tif")
ndwi_stress = rast("/mnt/CEPH_PROJECTS/WALDSCHAEDEN/working_folder/outputs/fordead_15/output_nb_periods_stress_NDWI_merged.tif")

# Create combined stress filter (2 or more stress periods)
filter_2or2 = (ndvi_stress > 1) | (ndwi_stress > 1)
filter_2or2 = project(filter_2or2, monthly)
filter_2or2[filter_2or2 == 0] <- NA

# Create mask for specific conditions (monthly == 38 and stress periods)
mask = (filter_2or2 == 1 & (monthly == 38))
mask[!mask] <- NA

# Apply stress period filtering
monthly_filtered = mask(monthly, mask, inverse = T)
yearly_filtered = mask(yearly, monthly_filtered)

cat("Stress filtering results:\n")
print("Before filtering (rows 30-39):")
print(freq(monthly)[30:39,])
print("After filtering (rows 30-39):")
print(freq(monthly_filtered)[30:39,])

# Save intermediate outputs
if (confirm_write(paste0(outfold, prefix, "_monthly_stress.tif"))) {
  writeRaster(monthly_filtered, paste0(outfold, prefix, "_monthly_stress.tif"), overwrite=TRUE)
}
if (confirm_write(paste0(outfold, prefix, "_yearly_stress.tif"))) {
  writeRaster(yearly_filtered, paste0(outfold, prefix, "_yearly_stress.tif"), overwrite=TRUE)
}

# =============================================================================
# SHADOW MASKING (VIPITENO/STERZING AREA)
# =============================================================================

cat("\n=== Applying shadow masking ===\n")

# Load QAI (Quality Assessment Index) data
qai = rast("/mnt/CEPH_PROJECTS/sao/SENTINEL-2/SentinelVegetationProducts/FORCE/level2_misc/level2_tiles/T32TPT/20241103_LEVEL2_SEN2A_QAI.tif")
mask_out = qai

# Process QAI bit values
qai_values = freq(mask_out)[,2]
qai_bit_values = as.data.frame(qai_values)

# Convert to binary representation
for(qai in 1:nrow(qai_bit_values)) {
  qai_bit_values$bit[qai] = paste(sapply(strsplit(paste(rev(intToBits(qai_bit_values$qai_values[qai]))),""),`[[`,2),collapse="")
}

# Extract relevant bit positions
qai_bit_values$bit = substr(qai_bit_values$bit, nchar(qai_bit_values$bit)-15, nchar(qai_bit_values$bit))

# Select cloud-free values (adjust bit pattern as needed)
qai_sel = grep("11", substr(qai_bit_values$bit, 4,5))
qai_bit_values_cloudfree = qai_bit_values[qai_sel,]

# Apply illumination state filtering if specified
if(exists("i_state_bits") && !is.na(i_state_bits)) {
  qai_bit_values_cloudfree = qai_bit_values_cloudfree[grepl(i_state_bits, substr(qai_bit_values_cloudfree$bit, 4, 5)), ]
}

# Create cloud-free mask
no_cloud_bit = as.vector(qai_bit_values_cloudfree$qai_values)
if(length(no_cloud_bit) > 0) mask_out[mask_out %in% no_cloud_bit] <- 0

# Load forest inspectorate boundaries and filter for Vipiteno
i = vect("/mnt/CEPH_PROJECTS/WALDSCHAEDEN/GIS/forest_administration/forest_inspectorates/ForestInspectorates_polygon.shp")
i = i[i$NAME_IT == "Vipiteno",]

# Apply spatial masking
mask_out[mask_out != 0] <- NA
mask_out = project(mask_out, yearly_filtered)
mask_out = mask(mask_out, i)
mask_out = crop(mask_out, yearly_filtered)

# Additional filtering for specific monthly values
mask_out[monthly_filtered != 38] <- NA

# Apply shadow mask
yearly_filtered = mask(yearly_filtered, mask_out, inverse = T)
monthly_filtered = mask(monthly_filtered, mask_out, inverse = T)

# Save intermediate outputs
if (confirm_write(paste0(outfold, prefix, "_monthly_shadow.tif"))) {
  writeRaster(monthly_filtered, paste0(outfold, prefix, "_monthly_shadow.tif"), overwrite=TRUE)
}
if (confirm_write(paste0(outfold, prefix, "_yearly_shadow.tif"))) {
  writeRaster(yearly_filtered, paste0(outfold, prefix, "_yearly_shadow.tif"), overwrite=TRUE)
}

# =============================================================================
# SINGLE PIXEL REMOVAL
# =============================================================================

cat("\n=== Removing single pixel detections ===\n")

# Create mask for connected components
onepx_mask = yearly_filtered
onepx_mask[!is.na(onepx_mask)] <- 0

# Convert to polygons and calculate areas
onepx_mask_shp = as.polygons(onepx_mask, dissolve = T)
onepx_mask_shp = disagg(onepx_mask_shp)
onepx_mask_shp$area = expanse(onepx_mask_shp, unit = "ha")

cat("Total area before filtering:", sum(onepx_mask_shp$area), "ha\n")

# Filter out small patches (>0.011 ha threshold)
onepx_mask_shp_filtered = onepx_mask_shp[onepx_mask_shp$area > 0.011]

cat("Total area after filtering:", sum(onepx_mask_shp_filtered$area), "ha\n")

# Convert back to raster
onepx_mask = rasterize(onepx_mask_shp_filtered, monthly)

# =============================================================================
# FINAL PROCESSING AND OUTPUT
# =============================================================================

cat("\n=== Finalizing outputs ===\n")

# Apply final masks
monthly_final = mask(monthly_filtered, onepx_mask)
yearly_final = mask(yearly_filtered, onepx_mask)

# Create vector versions
monthly_vect_25832 = as.polygons(monthly_final)
yearly_vect_25832 = as.polygons(yearly_final)

# Set up color tables and levels (ensure these variables are defined)
if(exists("periods") && exists("colors")) {
  levels(monthly_final) = periods_df
  coltab(monthly_final) = colors_df
}

if(exists("lev_year_df") && exists("col_year_df")) {
  levels(yearly_final) <- lev_year_df
  coltab(yearly_final) <- col_year_df
}

# Final frequency check
cat("Final frequencies:\n")
print("Monthly (rows 30-39):")
print(freq(monthly_final)[30:39,])
print("Yearly:")
print(freq(yearly_final))

# =============================================================================
# SAVE OUTPUTS (with confirmation)
# =============================================================================

cat("\n=== Saving final outputs ===\n")

# Define output filenames
raster_filename_month = paste0(outfold, prefix, "_", "monthly_damages_", update_name, "_2025_25832_final.tif")
raster_filename_year = paste0(outfold, prefix, "_", "yearly_damages_", update_name, "_2025_25832_final.tif")
vector_filename_month = paste0(outfold, prefix, "_", "monthly_damages_", update_name, "_2025_25832_final.shp")
vector_filename_year = paste0(outfold, prefix, "_", "yearly_damages_", update_name, "_2025_25832_final.shp")

# Write raster outputs with confirmation
if (confirm_write(raster_filename_month)) {
  writeRaster(monthly_final, raster_filename_month, datatype = "INT1U", overwrite = TRUE)
}
if (confirm_write(raster_filename_year)) {
  writeRaster(yearly_final, raster_filename_year, datatype = "INT1U", overwrite = TRUE)
}

# Write vector outputs with confirmation
if (confirm_write(vector_filename_month)) {
  writeVector(monthly_vect_25832, vector_filename_month, overwrite = TRUE)
}
if (confirm_write(vector_filename_year)) {
  writeVector(yearly_vect_25832, vector_filename_year, overwrite = TRUE)
}

print("Step 4 was successful") # Indicate successful completion of the processing step.
print(paste("Outputs written to:", outfold)) # Print the path where the outputs are saved.

dir.create(outfold, showWarnings = FALSE) # Create the output folder if it doesn't already exist, suppressing warnings.
files <- list.files(outfold, full.names = TRUE)[!grepl("final", list.files(outfold))] # List all files in the output folder that do not contain "final" in their name. These are considered intermediate files to be moved.
file.rename(files, file.path(outfold, "post_proc_intermediates", basename(files))) # Move the intermediate files to the 'post_proc_intermediates' subdirectory.
