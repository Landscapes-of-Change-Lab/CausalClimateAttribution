
## ---------------------------
##
## Script name: Random CA waypoints
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
## Notes: Samples points within CA
##   for counterfactual figure
##
## ---------------------------



# Load required libraries

librarian::shelf(sf, dplyr, tigris)

# Function to generate random points in California
generate_random_ca_points <- function(n_points = 100) {
  # Get California state boundary
  ca <- states() %>%
    filter(STUSPS == "CA") %>%
    st_transform(4326)  # Transform to WGS84 (lat/long)
  
  # Get the bounding box
  bbox <- st_bbox(ca)
  
  # Initialize empty list for points
  points_list <- list()
  point_count <- 0
  
  # Generate points until we have enough
  while(point_count < n_points) {
    # Generate a random point within the bounding box
    random_point <- st_point(c(
      runif(1, bbox["xmin"], bbox["xmax"]),
      runif(1, bbox["ymin"], bbox["ymax"])
    ))
    
    # Convert to sf object
    point_sf <- st_sfc(random_point, crs = 4326)
    
    # Check if point is within California
    if(st_intersects(point_sf, ca, sparse = FALSE)[1]) {
      point_count <- point_count + 1
      points_list[[point_count]] <- random_point
    }
  }
  
  # Convert points list to sf object
  points_sf <- st_sfc(points_list, crs = 4326)
  
  # Create sf data frame
  points_df <- st_sf(
    geometry = points_sf,
    id = 1:n_points,
    lat = st_coordinates(points_sf)[,2],
    lon = st_coordinates(points_sf)[,1]
  )
  
  return(points_df)
}

# Generate 100 random points
set.seed(123)  # for reproducibility
random_points <- generate_random_ca_points(100)

# View the first few points
head(random_points)

# Export to CSV (optional)
write.csv(st_drop_geometry(random_points), "Shapefiles/ca_random_points.csv", row.names = FALSE)

# Plot to verify (optional)
plot(st_geometry(states()[states()$STUSPS == "CA",]))
plot(st_geometry(random_points), add = TRUE, col = "red", pch = 20)
