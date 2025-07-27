# Level of interaction: every run.
#               - update path to previous outputs (used to fix past detections).
#               - run this once, then evaluate if optional masking is needed
#
# Script:       06_Integrate_And_Refine_Damage_Products.R
#
# Purpose:      Performs the final post-processing steps, applying additional 
#               filters (stress periods, single pixels), integrating previous 
#               results and generating the final output products.
#
# Inputs:
#   - Yearly and monthly damage rasters (from `05_Style_And_Project_Damage_Maps.R`).
#   - Monthly and yearly damage rasters from previous updates.
#   - NDVI and NDWI stress period rasters (for optional masking).
#
# Outputs:
#   - Final yearly and monthly damage GeoTIFFs and Shapefiles.
#
# Operations:
#   1. Loads current and previous damage rasters.
#   3. Applies stress period filtering based on combined NDVI/NDWI stress.
#   5. Removes single-pixel detections based on area.
#   2. Integrates old detections into current results.
#   7. Sets color tables and levels for final outputs.
#   8. Writes final GeoTIFF and Shapefile outputs with user confirmation.
#   9. Moves intermediate files to a dedicated subdirectory.

 
stress_filter = F

# Check here the index of periods to choose which ones to pick for the following steps
print(periods_df) # refer to the value


####
compare_outputs <- function(layer, old_layer, view = F) {
  f1 = freq(layer)
  f2 = freq(old_layer)
  t = f1$value
  count1 = f1$count
  count2 = f2$count
  length_diff = length(count1) - length(count2)
  if(length_diff > 0) count2_padded = c(count2, rep(NA, length_diff)) else count2_padded = count2
  px_diff = count1 - count2_padded
  percentage_diff = (px_diff) / count2_padded * 100
  df <- data.frame(
    t = t,
    count_old = count2_padded,
    count_new = count1,
    perc_diff = percentage_diff,
    px_diff = px_diff
  )
  print(df)
  if(view) View(df)
}

library(terra)

# Clear temporary files (optional)
# terra::tmpFiles(T, T, T, T)

# ==============================================================================
# CONFIGURATION AND INPUT LOADING
# ==============================================================================

# Load current detection results
yearly_original = rast(paste0(outfold, prefix, "_", "yearly_damages_", update_name, "_2025_25832.tif"))
yearly = yearly_original
monthly = rast(paste0(outfold, prefix, "_", "monthly_damages_", update_name, "_2025_25832.tif"))

# Load previous detection results for integration
# 1.INTERACTION ----
# Update the paths to the update from which the previous detections to be fixed
# are taken (usually the last one before the one currently being generated).
old_monthly = rast("/mnt/CEPH_PROJECTS/WALDSCHAEDEN/Products/FORDEAD_09_06_2025/changes_monthly_damages_june_2025_25832_final.tif") # The latest update that was delivered (using this pipeline). should run at regular intervals and come up with a rule here
old_yearly = rast("/mnt/CEPH_PROJECTS/WALDSCHAEDEN/Products/FORDEAD_09_06_2025/changes_yearly_damages_june_2025_25832_final.tif")

# Optional, can be useful for debugging
# Compare the final output from the previous run to the current output (pre post-processing)
# Can be deployed at any point of the script adjusting inputs to compare two different versions.
#compare_outputs(monthly, old_monthly)

# =============================================================================
# STRESS PERIOD FILTERING (OPTIONAL ADVANCED FILTERING)
# =============================================================================

# skip on first run
# eventually apply after output evaluation

if(stress_filter) {
  cat("\n=== Applying stress period filtering ===\n")
  
  # Load number of stress periods file
  # These files contain the number of periods when a pixel was detected as
  # disturbed, but then reverted to healty (stress periods as defined by fordead).
  # These areas, switching back and forth, can be considered unreliable detections 
  # especially during periods neighboring winter, where data quality is lower.
  # These areas can be masked out for specific periods, a decision to take 
  # following visual inspection of the final product.
  # FLAG: Hardcoded paths. These should be made dynamic to avoid issues when running the pipeline in different environments.
  ndvi_stress = rast("/mnt/CEPH_PROJECTS/WALDSCHAEDEN/working_folder/outputs/fordead_15/output_nb_periods_stress_NDVI_merged.tif")
  ndwi_stress = rast("/mnt/CEPH_PROJECTS/WALDSCHAEDEN/working_folder/outputs/fordead_15/output_nb_periods_stress_NDWI_merged.tif")
  
  # Create combined stress filter: if 2 or more stress periods for at least one of
  # the two indices occur, exclude the unreliably detected pixel
  filter_2or2 = (ndvi_stress > 1) | (ndwi_stress > 1)
  filter_2or2 = project(filter_2or2, monthly)
  filter_2or2[filter_2or2 == 0] <- NA
  
  # Create the mask limited to the specific period which needs to be targeted (monthly == 38 and stress periods flag)
  # these are typically one or more periods between October and April
  
  # 2(optional).INTERACTION ----
  mask = (filter_2or2 == 1 & (monthly == 38)) # set here the period(s) for which excessive noise is observed
  mask[!mask] <- NA
  
  # Apply mask
  monthly_filtered = mask(monthly, mask, inverse = T)
  yearly_filtered = mask(yearly, monthly_filtered)
  
  # Compare the final output from the previous run to the current output (pre post-processing)
  compare_outputs(monthly_filtered, old_monthly)
  
  # Save intermediate outputs
  writeRaster(monthly_filtered, paste0(outfold, prefix, "_monthly_stress.tif"), overwrite=TRUE, datatype = "INT1U")
  writeRaster(yearly_filtered, paste0(outfold, prefix, "_yearly_stress.tif"), overwrite=TRUE, datatype = "INT1U")
} else {
  monthly_filtered = monthly
  yearly_filtered = yearly
}

##

monthly_final = monthly_filtered
yearly_final = yearly_filtered

# Set up color tables and levels (ensure these variables are defined)
if(exists("periods") && exists("colors")) {
  levels(monthly_final) = periods_df
  coltab(monthly_final) = colors_df
}

if(exists("lev_year_df") && exists("col_year_df")) {
  levels(yearly_final) <- lev_year_df
  coltab(yearly_final) <- col_year_df
}

# =============================================================================
# INTEGRATE OLD DETECTIONS
# =============================================================================

cat("\n=== Integrating old detections ===\n")

compare_outputs(monthly_final, old_monthly)

# 3.INTERACTION ----
# define here the period before which detections of the previous run are to be
# considered permanent the comparison above should helo, together with a visual
# inspection of the results of the current and previous runs
#### Hard Fix Detections Before a Given Date in the Previous Update ####
fix_until = 37
old_monthly_filtered = old_monthly
old_monthly_filtered[old_monthly > fix_until] <- NA

# Create condition mask for where old detections should be preserved
condition1 = (monthly_final <= fix_until) # where the detections pre-date the fix date - delete
condition2 = !is.na(old_monthly_filtered) # where the old layer had some detections - write
condition12 = condition1 | condition2

# Integrate old detections into current results
monthly_final[condition12] <- old_monthly_filtered[condition12]
yearly_final[condition12] <- old_yearly[condition12]
yearly_final = mask(yearly_final, monthly_final)

# Save intermediate outputs
writeRaster(monthly_final, paste0(outfold, prefix, "_monthly_integration.tif"), overwrite=TRUE, datatype = "INT1U")
writeRaster(yearly_final, paste0(outfold, prefix, "_yearly_integration.tif"), overwrite=TRUE, datatype = "INT1U")

# =============================================================================
# SINGLE PIXEL REMOVAL
# =============================================================================

cat("\n=== Removing single pixel detections ===\n")

# Initiate mask
onepx_mask = yearly_final
onepx_mask[!is.na(onepx_mask)] <- 0

# Convert to polygons and calculate areas
onepx_mask_shp = as.polygons(onepx_mask, dissolve = T)
onepx_mask_shp = disagg(onepx_mask_shp)
onepx_mask_shp$area = expanse(onepx_mask_shp, unit = "ha")

cat("Total area before single-pixels filtering:", sum(onepx_mask_shp$area), "ha\n")

# Filter out small patches (>0.011 ha threshold)
onepx_mask_shp_filtered = onepx_mask_shp[onepx_mask_shp$area > 0.011]

cat("Total area after single-pixels filtering:", sum(onepx_mask_shp_filtered$area), "ha\n")

# Convert back to raster
onepx_mask = rasterize(onepx_mask_shp_filtered, monthly)
condition_1px = !condition12 & is.na(onepx_mask)

# Apply onepx mask only to non-fixed dates
monthly_final[condition_1px] = NA
yearly_final = mask(yearly_final, monthly_final)

compare_outputs(layer = monthly_final, old_layer = old_monthly)
compare_outputs(layer = yearly_final, old_layer = old_yearly)

# =============================================================================
# FINAL PROCESSING AND OUTPUT
# =============================================================================

cat("\n=== Finalizing outputs ===\n")

# Create vector versions
monthly_vect_25832 = as.polygons(monthly_final)
yearly_vect_25832 = as.polygons(yearly_final)

# Plot Optional ----
# Generate plot comparing evolution throughout years
# Final frequency check
f = freq(monthly_final)
f[(nrow(f)+1):(nrow(f)+10), ] = 0
# Convert the date column to Date type
f$Date <- as.Date(substr(f[,2], 1, 10))
# Extract year and day of year
f$Year <- format(f$Date, "%Y")
f$DayOfYear <- as.numeric(format(f$Date, "%j"))
plot(NULL, xlim = c(80, 320), ylim = range(f[,3], na.rm = TRUE),
     xlab = "Day of Year", ylab = "Value", main = "Yearly Series Overlaid")
# Assign colors
years <- unique(na.omit(f$Year))
colors <- rainbow(length(years))
# Plot each year
for (i in seq_along(years)) {
  subset_data <- f[f$Year == years[i], ]
  lines(subset_data$DayOfYear, subset_data[,3], col = colors[i], lwd = 2)
}
legend("topright", legend = years, col = colors, lty = 1, lwd = 2)


# =============================================================================
# SAVE OUTPUTS
# =============================================================================

cat("\n=== Saving final outputs ===\n")

# Define output filenames
raster_filename_month = paste0(outfold, prefix, "_", "monthly_damages_", update_name, "_2025_25832_final.tif")
raster_filename_year = paste0(outfold, prefix, "_", "yearly_damages_", update_name, "_2025_25832_final.tif")
vector_filename_month = paste0(outfold, prefix, "_", "monthly_damages_", update_name, "_2025_25832_final.shp")
vector_filename_year = paste0(outfold, prefix, "_", "yearly_damages_", update_name, "_2025_25832_final.shp")

# Write raster outputs with confirmation
writeRaster(monthly_final, raster_filename_month, datatype = "INT1U", overwrite = TRUE)
writeRaster(yearly_final, raster_filename_year, datatype = "INT1U", overwrite = TRUE)
writeVector(monthly_vect_25832, vector_filename_month, overwrite = TRUE)
writeVector(yearly_vect_25832, vector_filename_year, overwrite = TRUE)

cat("Processing completed. Inspect the outputs. \n")