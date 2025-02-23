
## ---------------------------
##
## Script name: Extracts PRISM data to waypoints
##
## Author: Dr. Joan Dudney
##
## Date Created: 2024-10-22
##
## Copyright (c) Joan Dudney, 2024
## Email: dudney@ucsb.edu
##
## ---------------------------
##
## Notes: Sources the "extract_prism_function.R" script
##   to extract prism data to waypoints of choice
##
## ---------------------------

## pacakges
librarian::shelf(terra,sf, tidyverse, lubridate, here)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1. Reading and cleaning data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source(here("Scripts", "Functions", "extract_prism_function.R"))

## writing the directory
base_dir <- here("Data", "PRISM_tmax")

## reading in waypoints
waypoints <- read_csv(here("Data", "SpatialFiles","ca_random_points.csv"))

# Convert points to sf object
waypoints_sf <- st_as_sf(waypoints,
                         coords = c("lon", "lat"),
                         crs = 4326)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2. Extracting prism data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## extract prism data
prism_data <- extract_prism_data(base_dir, waypoints_sf)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3. Summarizing and exporting data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## summary dataset

summary_dat <- prism_data %>% 
  group_by(waypoint_id, latitude, longitude, year) %>% 
  summarize(m.tmax = mean(tmax, na.omit=T)) %>% 
  rename(tmax = m.tmax)

summary_dat %>% 
  ggplot(aes(x=year, y=tmax))+
  geom_point()

#write_csv(summary_dat, here("Data", "extracted_prism_data.csv"))
