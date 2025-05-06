
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
librarian::shelf(terra, sf, tidyverse, lubridate, here)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1. Read in data, project, and visualize
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## data to run all analyses
#itrdbdat <- read_csv(here("Data", "Primary_data", "ITRDB_PIAL_latlon.csv"))

itrdbdat <- read_csv(here("Data", "ITRB_species_latlon.csv"))  
  
domsp <- itrdbdat %>% 
  group_by(species_id) %>% 
  summarise(unique_plots = n_distinct(collection_id))

itrdbdat_sp <- itrdbdat %>% 
  filter(species_id %in% c("psme", "pico", "pied", "tsme", "pcgl","pcen", "pifl"))

## filter itrdb and change column names to work with function
itrdb_simp <- itrdbdat_sp %>% 
  dplyr::select(collection_id, latitude, longitude) %>% 
  as.data.frame() %>%
  unique() %>% 
  rename(lon = longitude, lat = latitude) %>% 
  rename(id = collection_id) %>% 
  filter(!grepl("CAN", id)) %>% 
  filter(!grepl("MEX", id))

# Convert points to sf object
waypoints_itrdb <- st_as_sf(itrdb_simp, 
                         coords = c("lon", "lat"), 
                         crs = 4326) # WGS84 coordinate system

# interactive viewing mode
tmap_mode("view")

# Create map
tm_map <- tm_shape(waypoints_itrdb) +
  tm_dots(col = "id", 
          palette = "viridis",
          title = "Collection ID",
          size = 1,
          popup.vars = c("id")) +
  tm_basemap(server = c("OpenStreetMap", "Esri.WorldImagery"))

# Display the map
tm_map

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2. Extract waypoints
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

source(here("Scripts", "Functions", "extract_climate_data_prism_function.R"))

## writing the directory
base_dir <- here("Data", "PRISM_tmax_ppt")

climate_data <- extract_prism_data(base_dir, waypoints_itrdb, 
                                   variables = c("tmax", "ppt")) %>% 
  filter(!is.na(value))


#write_csv(climate_data, here("Data", "raw_prism_itrdb_sevenDomSp.csv"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3. Summarizing and exporting data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## summary dataset
summary_dat <- climate_data %>%
  dplyr::select(waypoint_id, variable,value, month, year, longitude, latitude) %>% 
  rename(collection_id = waypoint_id) %>% 
  mutate(growing=ifelse(month%in%c(10:12), year+1, year)) %>% 
  filter(!month%in%c(7:9)) %>% ## don't include July, August, September
  filter(growing!=1900) %>% 
  pivot_wider(names_from=variable, values_from=value) %>% 
  group_by(collection_id, growing) %>% 
  summarize(tmax = mean(tmax, na.rm=T), ppt = sum(ppt, na.rm=T))%>%
  rename(year=growing) %>% 
  ungroup() 

#write_csv(summary_dat, here("Data", "extracted_prism_data_itrdb.csv"))




# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# #
# #
# #                            Extracting prism data for SN
# #
# #
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
# ## data to run all analyses
# rwldat <- read_csv(here("Data", "cleaned_RWL_climdat.csv")) %>% 
#   dplyr::select(plot_id_needle, lat, long) %>% 
#   rename(lon = long, id = plot_id_needle) %>% 
#   distinct()
# 
# ## combining datasets
# comb_latlong <- itrdbdat %>%
#   full_join(rwldat)
# 
# # Convert points to sf object
# waypoints_comb <- st_as_sf(comb_latlong, 
#                             coords = c("lon", "lat"), 
#                             crs = 4326) # WGS84 coordinate system
# 
# # interactive viewing mode
# tmap_mode("view")
# 
# # Create map
# tm_map <- tm_shape(waypoints_comb) +
#   tm_dots(col = "id", 
#           palette = "viridis",
#           title = "Plot ID",
#           size = .5,
#           popup.vars = c("id")) +
#   tm_basemap(server = "CartoDB.Positron")
# 
# tm_map
# 
# ## extracting prism data
# # Convert points to sf object
# waypoints_rwl <- st_as_sf(rwldat, 
#                            coords = c("lon", "lat"), 
#                            crs = 4326) # WGS84 coordinate system
# 
# source(here("Scripts", "Functions", "extract_climate_data_prism_function.R"))
# 
# ## writing the directory
# base_dir <- here("Data", "PRISM_tmax_ppt")
# 
# climate_data <- extract_prism_data(base_dir, waypoints_rwl, 
#                                    variables = c("tmax", "ppt"))
# 
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # 3. Summarizing and exporting data
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# ## summary dataset
# summary_dat <- climate_data %>%
#   dplyr::select(waypoint_id, variable,value, month, year, longitude, latitude) %>% 
#   rename(plot_id_needle = waypoint_id) %>% 
#   mutate(growing=ifelse(month%in%c(10:12), year+1, year)) %>% 
#   filter(!month%in%c(7:9)) %>% ## don't include July, August, September
#   filter(growing!=1900) %>% 
#   pivot_wider(names_from=variable, values_from=value) %>% 
#   group_by(plot_id_needle, growing) %>% 
#   summarize(tmax = mean(tmax, na.rm=T), ppt = sum(ppt, na.rm=T))%>%
#   rename(year=growing) %>% 
#   ungroup() 
# 
# #write_csv(summary_dat, here("Data", "prism_data_SN_new.csv"))
# rwldat <- read_csv(here("Data", "cleaned_RWL_climdat.csv"))
# 

