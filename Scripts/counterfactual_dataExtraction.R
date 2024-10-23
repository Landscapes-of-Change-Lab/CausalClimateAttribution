## ---------------------------
##
## Script name: Extacting counterfactual temp data
##
## Author: Dr. Joan Dudney
##
## Date Created: 2024-10-21
##
## Copyright (c) Joan Dudney, 2024
## Email: dudney@ucsb.edu
##
## ---------------------------
##
## Notes: this code extracts estimated counterfactual tmax nc files
##   to random waypoints across CA
## source data: https://zenodo.org/records/5036364
##
## ---------------------------

## Packages

librarian::shelf(terra, here, tidyverse, sf, tmap, lubridate)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1. Reading in files
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## sourcing needed function
source(here("Scripts", "Functions", "extract_counterfactual_function.R"))

## CA boundary
ca <- vect(here("Data","SpatialFiles", "CA_polygon_nad83zone11.shp"))
plot(ca)

## Random points
random_points <- read_csv(here("Data","SpatialFiles", "ca_random_points.csv"))

# Convert points to sf object
points_sf <- st_as_sf(random_points, 
                      coords = c("lon", "lat"), 
                      crs = 4326)

# Visualize the data
tmap_mode("view")

tm_shape(points_sf) +
  tm_dots(col = "blue", size = 0.5, alpha = 0.6) +
  tm_layout(title = "Random Sampling Points in California",
            title.position = c("center", "top"),
            legend.position = c("right", "bottom")) +
  tm_basemap("OpenStreetMap")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2. Extracting climate data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## counterfactual data
nc_files <- list.files(here("Data",  "Counterfactuals"), 
                       pattern = "\\.nc$", 
                       full.names = TRUE)
raster_list <- rast(nc_files)

extracted_counter <- terra::extract(raster_list, points_sf)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## 3. Process all the files
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Organize the data
all_data <- lapply(nc_files, function(file) {
  process_nc_file(file, points_sf)
}) %>%
  bind_rows()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## 4. Clean and organize the final dataset
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## easy extraction using terra
terra_vect <- vect(points_sf)

final_data <- all_data %>%
  arrange(point_id, date) %>%
  # Convert temperature to Celsius 
  mutate(temperature = temperature - 273.15) %>%
  # Add spatial coordinates from waypoints
  left_join(
    data.frame(
      point_id = 1:nrow(terra_vect),
      longitude = crds(terra_vect)[,1],
      latitude = crds(terra_vect)[,2]
    ),
    by = "point_id"
  )

# Visualize the results
summarydat <- final_data %>% 
  group_by(year, point_id, longitude, latitude) %>% 
  summarize(tmax=mean(temperature))

summarydat %>% 
  ggplot(aes(x=year, y=tmax))+
  geom_point()

# Save the results 
#write.csv(summarydat, "Data/extracted_counterfactual_data.csv", row.names = FALSE)
