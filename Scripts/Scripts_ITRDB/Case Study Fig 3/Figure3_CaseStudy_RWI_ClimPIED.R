## ---------------------------
##
## Script name: 
##
## Author: Dr. Joan Dudney
##
## Date Created: 2024-09-07
##
## Copyright (c) Joan Dudney, 2024
## Email: dudney@ucsb.edu
##
## ---------------------------
##
## Notes: Creates RWI and climate variable
##   Figures in Figure 3
##
## ---------------------------

## packages
librarian::shelf(sjPlot, ggeffects, patchwork, tidyverse, here,terra, sf,tmap,
                 lme4, plotrix, ggpubr, mgcv, nlme, fixest, plotrix, egg, ggpmisc,
                 mvtnorm, clubSandwich, rasterVis, broom.mixed, scales,RColorBrewer)

theme_set(
  theme_bw(base_size = 15)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)

select <- dplyr::select


## itrdb data
itrdbdat <- read_csv(here("Data", "ITRB_species_latlon.csv")) %>%
  filter(species_id == "pied") 

## isolating sites for PIED
pied <- itrdbdat %>% 
  select(species_id, collection_id) %>% 
  distinct()
  

## prism data
#prism_itrdb <- read_csv(here("Data", "raw_prism_data_itrdb.csv"))
prism_itrdb <- read_csv(here("Data", "raw_prism_itrdb_sevenDomSp.csv"))

pied_clim <- prism_itrdb %>%
  select(waypoint_id, variable,value, month, year, longitude, latitude) %>% 
  rename(collection_id = waypoint_id) %>%
  left_join(pied) %>% 
  filter(!is.na(species_id))

ggplot(pied_clim, aes(x=longitude, y=latitude, color = collection_id))+
  geom_point()

climdat <- pied_clim %>% 
  mutate(growing=ifelse(month%in%c(10:12), year+1, year)) %>% 
  filter(!month%in%c(7:9)) %>% ## don't include July, August, September
  filter(growing!=1900) %>% 
  pivot_wider(names_from=variable, values_from=value) %>% 
  group_by(collection_id, growing) %>% 
  summarize(tmax = mean(tmax, na.rm=T), ppt = sum(ppt, na.rm=T))%>%
  rename(year=growing) %>% 
  ungroup() %>% 
  rename(plot = collection_id)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#                           FIGURES of RWI and Climate for Figure 3
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Temperature and precip trends
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
meantemp <- climdat %>% 
  group_by(year) %>% 
  summarize(m.temp = mean(tmax, na.omit=T))


tempfig <- meantemp %>% 
  filter(year<2019) %>% 
  ggplot(aes(x = year, y=m.temp))+
  geom_point(size=.2)+
  geom_line(linewidth=.1)+
  geom_smooth(size=.7, linetype=2, color="#91376F", fill = "#91376F", alpha=.05)+
  ylab("Temperature (°C)")+
  xlab("Year")+
  scale_x_continuous(breaks = seq(1900, 2020, by = 20))


meanprecip <- climdat %>% 
  group_by(year) %>% 
  summarize(m.ppt = mean(ppt, na.omit=T))

summary(lm(m.ppt~year, data=meanprecip))

pptfig <- meanprecip %>% 
  filter(year<2019) %>% 
  ggplot(aes(x = year, y=m.ppt))+
  geom_point(size=.2)+
  geom_line(linewidth=.1)+
  geom_smooth(size=.7, linetype=2, color="#303077", fill = "#303077", alpha=.05)+
  ylab("Precipitation (mm)")+
  xlab("Year")+
  scale_x_continuous(breaks = seq(1900, 2020, by = 20))


pptfig / tempfig


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# RWI FIGURE
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# RWI FIGURE
mean_rwidat <- itrdbdat %>% 
  filter(year>1899) %>% 
  group_by(collection_id, year) %>% 
  summarize(rwi = mean(rwi, na.omit=T))

meanrwi <- itrdbdat %>% 
  group_by(year) %>% 
  summarize(meanrwi = mean(rwi, na.omit=T))

## combine datasets
combrwi <- mean_rwidat %>% 
  left_join(meanrwi)


# Function to darken colors
darken_color <- function(color, factor = 0.3) {
  rgb_val <- col2rgb(color)
  darkened_rgb <- rgb_val * (1 - factor)
  darkened_rgb <- pmax(pmin(darkened_rgb, 255), 0) # Ensure values are within 0-255
  rgb(darkened_rgb[1], darkened_rgb[2], darkened_rgb[3], maxColorValue = 255)
}

# Generate a blue-green palette with 27 colors from RColorBrewer
base_palette <- brewer.pal(n = 9, name = "Blues")

# Extend the palette to 27 colors by interpolating between the base colors
extended_palette <- colorRampPalette(base_palette)(108)

# Darken the colors using the custom darken_color function
darkened_palette <- sapply(extended_palette, darken_color, factor = 0.1)

theme_set(
  theme_bw(base_size = 22)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)

rwifig <- combrwi %>% 
  ggplot(aes(x = year, y=rwi, color=collection_id))+
  geom_point(size=.2)+
  geom_line(linewidth=.1)+
  scale_color_manual(values = extended_palette) +
  geom_line(aes(x=year, y=meanrwi), color = "black", size = 1)+
  #scale_color_manual(values="darkgrey")+
  ylim(0,3)+
  theme(legend.position = "none")+
  ylab("RWI")+
  xlab("Year")+
  scale_x_continuous(breaks = seq(1900, 2020, by = 20))

rwifig

combrwi %>% 
  filter(year>1980) %>% 
  ggplot(aes(x = year, y=rwi, color=collection_id))+
  geom_point(size=.2)+
  geom_line(linewidth=.1)+
  guides(color=F)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# shapefile of pied plots
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pied_itrdb <- itrdbdat %>% 
  filter(species_id == "pied") %>% 
  rename(lon = longitude, lat = latitude)

waypoints_pied <- st_as_sf(pied_itrdb, 
                            coords = c("lon", "lat"), 
                            crs = 4326) # WGS84 coordinate system
#write_sf(waypoints_pied, here("Data", "pied_itrdb_plots.shp"))
# interactive viewing mode
tmap_mode("view")

# Create map
tm_map <- tm_shape(waypoints_pied) +
  tm_dots(col = "collection_id", 
          palette = "viridis",
          title = "Collection ID",
          size = 1) +
          #popup.vars = c("id")) +
  tm_basemap(server = c("OpenStreetMap", "Esri.WorldImagery"))

# Display the map
tm_map



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

meantemp <- paneldat %>% 
  group_by(year) %>% 
  summarize(m.temp = mean(tmax, na.omit=T))


tempfig <- meantemp %>% 
  filter(year<2019) %>% 
  ggplot(aes(x = year, y=m.temp))+
  geom_point(size=.2)+
  geom_line(linewidth=.1)+
  geom_smooth(size=.7, linetype=2, color="#91376F", fill = "#91376F", alpha=.05)+
  ylab("Temperature (°C)")+
  xlab("Year")+
  scale_x_continuous(breaks = seq(1900, 2020, by = 20))


meanprecip <- climdat %>% 
  group_by(year) %>% 
  summarize(m.ppt = mean(ppt, na.omit=T))

summary(lm(m.ppt~year, data=meanprecip))

pptfig <- meanprecip %>% 
  filter(year<2019) %>% 
  ggplot(aes(x = year, y=m.ppt))+
  geom_point(size=.2)+
  geom_line(linewidth=.1)+
  geom_smooth(size=.7, linetype=2, color="#303077", fill = "#303077", alpha=.05)+
  ylab("Precipitation (mm)")+
  xlab("Year")+
  scale_x_continuous(breaks = seq(1900, 2020, by = 20))


pptfig / tempfig
