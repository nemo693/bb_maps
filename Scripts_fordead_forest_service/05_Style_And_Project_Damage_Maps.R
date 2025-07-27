# Script:       05_Style_And_Project_Damage_Maps.R
#
# Level of interaction: update (once a year)
#               - update colors and labels with the first update of every year.
#               - update year in naming (at the first update of each year).
#
# Purpose:      Prepares the yearly and monthly damage products for final output
#               and visualization by applying color palettes, projecting data,
#               and converting formats.
#
# Inputs:
#   - Yearly and monthly damage rasters from `prod_fol`.
#   - `update_name` (from `00_RUN_Province_Damage_Update.R`).
#
# Outputs:
#   - Yearly and monthly damage GeoTIFFs and Shapefiles - semi-final versions.
#
# Operations:
#   1. Loads the yearly and monthly damage rasters.
#   2. Defines and applies predefined color palettes and labels to these rasters.
#   3. Projects the data to EPSG:25832.
#   4. Converts the raster data to vector (polygon) format.
#   5. Includes post-processing adjustments for specific date labels (November, winter periods) and recalculates the `last_detectable_date`.
#   6. Writes the final GeoTIFF and Shapefile outputs for both yearly and monthly products.


library(terra)

# .0 Interaction: enter current year ----
year_filename = 2025

# Configuration
overwrite <- TRUE  # Set this as needed
prefix <- "changes"
prod_fol <- for_out_path
dir.create(outfold, showWarnings = FALSE) #  defined in 00_...

# ============================================================================
# PART 1: PROCESS YEARLY DAMAGES
# ============================================================================

print("Processing yearly damages...")

# Load the yearly damage raster, which is the output from the previous step.
r_yearly <- rast(paste0(prod_fol, "NDVI_merged_masked_with_NDWI_province_only_annual.tif"))

#.1 Interaction: once a year pick the color for the new year ----
# Define the color palette for yearly damage visualization.
col_year <- c(
  "#6f8592", #2020
  "#a6cee3", 
  "#0000FF", 
  "#FF9900",
  "#d80077",
  "#500074" # 2025
)
# Create a data frame for the color table, mapping values to colors.
col_year_df <- data.frame(
  value = 0:(length(col_year)-1), 
  color = col_year
)

# .2 Interaction: once a year, define the new label ----
# Define labels for each year of change.
lev_year <- c(
  'Change 2020',
  'Change 2021',
  'Change 2022',
  'Change 2023',
  'Change 2024',
  'Change 2025'
)
# Create a data frame for the levels, mapping values to labels.
lev_year_df <- data.frame(
  value = 0:(length(lev_year)-1), 
  label = lev_year
)

# Apply the defined color table and levels to the yearly raster.
coltab(r_yearly) <- col_year_df
levels(r_yearly) <- lev_year_df

# Project the yearly data to EPSG:25832 for standardized spatial reference.
print("Projecting yearly data to EPSG:25832...")
yearly_25832 <- project(r_yearly, "EPSG:25832", threads = TRUE)

# Convert the projected yearly raster data to a vector (polygon) format.
yearly_vect_25832 <- as.polygons(yearly_25832)

# Write the processed yearly raster output to a GeoTIFF file.
print("Writing yearly outputs...")
writeRaster(yearly_25832, 
            paste0(outfold, prefix, "_yearly_damages_", update_name, "_",  year_filename, "_25832.tif"), 
            datatype = "INT1U",
            overwrite = overwrite)

# Write the processed yearly vector output to a Shapefile.
writeVector(yearly_vect_25832, 
            paste0(outfold, prefix, "_yearly_damages_", update_name, "_", year_filename, "_25832.shp"),
            overwrite = overwrite)

# Clean up unused objects from memory to free up resources.
rm(r_yearly, yearly_25832, yearly_vect_25832)
gc()

# ============================================================================
# PART 2: PROCESS MONTHLY DAMAGES
# ============================================================================

print("Processing monthly damages...")

# Load the monthly damage raster, which is the output from the previous step.
r_monthly <- rast(paste0(prod_fol, "NDVI_merged_masked_with_NDWI_province_lafis.tif"))

# Extract unique period values from the monthly raster to define monthly periods.
periods <- unique(r_monthly)[,1]

# Create a data frame for the monthly periods, mapping values to labels.
periods_df <- data.frame(
  value = 0:(length(periods)-1), 
  label = periods
)

# .3 Interaction: once a year, define the new colors for the season ----
# Define a comprehensive color palette for monthly damage visualization.
colors <- c(
  "#f0f0f0", "#d8d8d8", "#c1c1c1", "#a9a9a9", "#919191", "#797979", "#626262", "#4a4a4a",
  "#a6cee3", "#8fb1e4", "#7993e4", "#6276e5", "#4c58e5", "#353be6", "#1f1de6", "#0800e7",
  "#99ff00", "#85ed04", "#71dc08", "#5dca0c", "#4ab811", "#36a615", "#229519", "#0e831d",
  "#ff9900", "#ed8b06", "#dc7e0b", "#ca7011", "#b86216", "#a6541c", "#954721", "#833927",
  "#E6E6FA", "#D8BFD8", "#E0B0FF", "#DA70D6", "#9370DB", "#9966CC", "#B19CD9", "#A69CD9",
  "#CB0000", "#D01F1F", "#D53E3E", "#DA5D5D", "#DF7B7B", "#E49A9A", "#E59F9F", "#E9B9B9"
)

# Trim the color palette to match the number of unique periods.
colors = colors[1:length(periods)]

# Create a data frame for the color table, mapping values to colors for monthly data.
colors_df <- data.frame(
  value = 0:(length(colors)-1), 
  color = colors
)

# Apply the defined color table and levels to the monthly raster.
coltab(r_monthly) <- colors_df
levels(r_monthly) <- periods_df

# Project the monthly data to EPSG:25832 for standardized spatial reference.
print("Projecting monthly data to EPSG:25832...")
monthly_25832 <- project(r_monthly, "EPSG:25832", threads = TRUE)

# ============================================================================
# PART 2b: POST-PROCESSING ADJUSTMENTS
# ============================================================================

print("Applying post-processing adjustments...")

# Get current labels for modification from the projected monthly data.
current_levels <- levels(monthly_25832)[[1]]

# Fix November labels: change the day to the 15th of the month for consistency.
november_mask <- grepl("-11-", current_levels$label)
if(any(november_mask)) {
  current_levels$label[november_mask] <- paste0(substr(current_levels$label[november_mask], 1, 21), "15")
}

# Determine the last detectable date by checking image paths from two orbits.
images_paths <- c(
  paste0(prod_fol, "/fordead_output_NDVI_X0001_Y0004/VegetationIndex"),
  paste0(prod_fol, "/fordead_output_NDVI_X0003_Y0004/VegetationIndex")
)

last_detectable_date <- ""
for(path in images_paths) {
  if(dir.exists(path)) {
    images <- list.files(path)
    if(length(images) > 2) {
      # Extract date strings from image filenames (assuming format "YYYYMMDD_HHMMSS_DATE").
      dates <- substr(images, 17, 26)
      # Get the second to last date, which is considered the last detectable date.
      current_last <- dates[length(dates) - 2]
      # Update last_detectable_date if the current one is more recent.
      if(current_last > last_detectable_date) {
        last_detectable_date <- current_last
      }
    }
  }
}

# Update the last label in the monthly periods with the actual last detectable date.
if(nchar(last_detectable_date) > 0) {
  last_idx <- nrow(current_levels)
  current_levels$label[last_idx] <- paste0(substr(current_levels$label[last_idx], 1, 13), last_detectable_date)
}

# Fix winter periods: April entries are adjusted to represent the winter period spanning from November 16th of the previous year to April 30th of the current year.
winter_mask <- grepl("04-01", current_levels$label)
if(any(winter_mask)) {
  winter_labels <- current_levels$label[winter_mask]
  year <- as.numeric(substr(winter_labels, 1, 4))
  current_levels$label[winter_mask] <- paste0((year - 1), "-11-16 - ", year, "-04-30")
}

# Apply the updated labels to the monthly raster.
levels(monthly_25832) <- current_levels

# Convert the projected monthly raster data to a vector (polygon) format.
monthly_vect_25832 <- as.polygons(monthly_25832)

# ============================================================================
# PART 2c: WRITE FINAL OUTPUTS
# ============================================================================

print("Writing monthly outputs...")
# Write the processed monthly raster output to a GeoTIFF file.
writeRaster(monthly_25832, 
            paste0(outfold, prefix, "_monthly_damages_", update_name, "_", year_filename, "_25832.tif"), 
            datatype = "INT1U", 
            overwrite = overwrite)
# Write the processed monthly vector output to a Shapefile.
writeVector(monthly_vect_25832, 
            paste0(outfold, prefix, "_monthly_damages_", update_name, "_", year_filename, "_25832.shp"),
            overwrite = overwrite)

# Clean up unused objects from memory to free up resources.
rm(r_monthly, monthly_25832, monthly_vect_25832)

print("Step 4 was successful")
print(paste("Outputs written to:", outfold))

# Clean up unused objects from memory to free up resources.
gc()

