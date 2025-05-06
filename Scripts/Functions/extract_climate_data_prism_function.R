## ---------------------------
##
## Script name: Function to parse and extract multiple PRISM data variables
##
## Based on original by: Dr. Joan Dudney
## Modified by: [Your Name]
##
## Date Created: 2024-10-22
## Date Modified: 2025-04-22
##
## Copyright (c) Joan Dudney, 2024
## Email: dudney@ucsb.edu
##
## ---------------------------
##
## Notes: Works with multiple PRISM variables (tmax, ppt, vpdmax, etc.)
##   
##
## ---------------------------


## packages
librarian::shelf(terra, sf, tidyverse, lubridate)


# Function to parse dates and variable types from PRISM folder names
parse_prism_info <- function(filepath) {
  # Get the folder name
  folder_name <- basename(dirname(filepath))
  
  # Extract variable type (tmax, ppt, vpdmax, etc.)
  var_type <- str_extract(folder_name, "(?<=PRISM_)[a-z]+(?=_)")
  
  # Extract date from folder name (matches patterns like PRISM_tmax_stable_4kmM3_198209_bil
  # or PRISM_ppt_provisional_4kmM3_202005_bil)
  date_str <- str_extract(folder_name, "\\d{6}(?=_bil)")
  
  if(is.na(date_str) || is.na(var_type)) {
    warning("Could not extract date or variable type from folder: ", folder_name)
    return(NULL)
  }
  
  # Parse the YYYYMM format
  year <- substr(date_str, 1, 4)
  month <- substr(date_str, 5, 6)
  
  # Since we're dealing with monthly data, set day to 1
  day <- "01"
  
  # Validate components
  year_num <- as.numeric(year)
  month_num <- as.numeric(month)
  
  if(is.na(year_num) || is.na(month_num)) {
    warning("Invalid date components in: ", date_str)
    return(NULL)
  }
  
  if(year_num < 1800 || year_num > 2100 || month_num < 1 || month_num > 12) {
    warning("Date out of valid range: ", date_str)
    return(NULL)
  }
  
  tryCatch({
    date <- ymd(paste(year, month, day, sep="-"))
    return(list(
      variable = var_type,
      year = year_num,
      month = month_num,
      day = 1,  # Always 1 for monthly data
      date = date,
      is_provisional = grepl("provisional", folder_name),
      is_stable = grepl("stable", folder_name),
      resolution = str_extract(folder_name, "4km\\w+")
    ))
  }, error = function(e) {
    warning("Error creating date from: ", date_str)
    return(NULL)
  })
}

# Function to find and process PRISM files
find_prism_files <- function(base_dir, variables = NULL) {
  # List all subdirectories that match PRISM pattern
  all_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
  
  # Create pattern to match any PRISM variable or specific variables if provided
  if(is.null(variables)) {
    var_pattern <- "[a-z]+"
  } else {
    var_pattern <- paste(variables, collapse = "|")
  }
  
  # Pattern to match PRISM directories
  dir_pattern <- paste0("PRISM_", var_pattern, "_(?:stable|provisional)_4km\\w+_\\d{6}_bil")
  prism_dirs <- all_dirs[grepl(dir_pattern, basename(all_dirs))]
  
  if(length(prism_dirs) == 0) {
    stop("No PRISM directories found matching expected pattern")
  }
  
  # For each directory, find the .bil file
  all_files <- vector()
  for(dir in prism_dirs) {
    bil_files <- list.files(
      path = dir,
      pattern = "\\.bil$",
      full.names = TRUE
    )
    if(length(bil_files) > 0) {
      # Take the first .bil file if multiple exist
      all_files <- c(all_files, bil_files[1])
    }
  }
  
  return(all_files)
}

# Main function to extract PRISM data
extract_prism_data <- function(base_dir, waypoints_sf, variables = NULL) {
  cat("Starting PRISM data extraction...\n")
  
  # Find all PRISM files, optionally filtering by variable type
  cat("Searching for PRISM files...\n")
  prism_files <- find_prism_files(base_dir, variables)
  
  if(length(prism_files) == 0) {
    stop("No PRISM files found")
  }
  
  cat(sprintf("Found %d PRISM files\n", length(prism_files)))
  
  # Get CRS from first PRISM file
  first_raster <- try(rast(prism_files[1]))
  if(inherits(first_raster, "try-error")) {
    stop("Cannot read first PRISM file: ", prism_files[1])
  }
  prism_crs <- crs(first_raster)
  
  # Transform waypoints to match PRISM CRS
  cat("Transforming waypoints to match PRISM CRS...\n")
  waypoints_transformed <- st_transform(waypoints_sf, prism_crs)
  waypoints_vect <- vect(waypoints_transformed)
  
  # Initialize results dataframe
  results <- data.frame()
  
  # Process files with a progress bar
  cat("Processing PRISM files...\n")
  pb <- txtProgressBar(min = 0, max = length(prism_files), style = 3)
  
  for(i in seq_along(prism_files)) {
    file <- prism_files[i]
    tryCatch({
      # Read the PRISM raster
      prism_raster <- rast(file)
      
      # Extract values for waypoints
      extracted_values <- terra::extract(prism_raster, waypoints_vect)
      
      # Get variable type, date and stability information
      info <- parse_prism_info(file)
      
      if(!is.null(info)) {
        # Create temporary dataframe for this file
        temp_df <- data.frame(
          waypoint_id = waypoints_sf$id,
          longitude = st_coordinates(waypoints_transformed)[,1],
          latitude = st_coordinates(waypoints_transformed)[,2],
          value = extracted_values[,2],
          variable = info$variable,
          year = info$year,
          month = info$month,
          date = info$date,
          is_provisional = info$is_provisional,
          is_stable = info$is_stable,
          resolution = info$resolution,
          filename = basename(file),
          filepath = file
        )
        
        # Append to results
        results <- bind_rows(results, temp_df)
      }
    }, error = function(e) {
      warning("Error processing file: ", file, "\n", e$message)
    })
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Convert to wider format if desired
  # results_wide <- results %>%
  #   pivot_wider(
  #     id_cols = c(waypoint_id, longitude, latitude, year, month, date),
  #     names_from = variable,
  #     values_from = value
  #   )
  
  return(results)
}

# Example usage:
# 1. Extract all available PRISM variables
# all_data <- extract_prism_data("path/to/prism/data", waypoints)
#
# 2. Extract specific variables
# climate_data <- extract_prism_data("path/to/prism/data", waypoints, 
#                                   variables = c("tmax", "ppt", "vpdmax"))
