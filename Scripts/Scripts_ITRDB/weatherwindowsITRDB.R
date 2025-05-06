## ---------------------------
##
## Script name: ITRDB climate windows
##
## Author: Dr. Joan Dudney
##
## Date Created: 2025-04-23
##
## Copyright (c) Joan Dudney, 2025
## Email: dudney@ucsb.edu
##
## ---------------------------
##
## Notes: creates different climate windows 
##   for ITRDB sites with PIAL
##
## ---------------------------

select <- dplyr::select

## PIAL data
## itrdb data
# itrdbdat <- read_csv(here("Data", "Primary_data", "ITRDB_PIAL_latlon.csv")) %>% 
#   select(-c(species_id, family, genus, gymno_angio))

## PIPO data
itrdbdat <- read_csv(here("Data", "ITRB_species_latlon.csv")) %>%
  filter(species_id %in% c("psme", "pico", "pied", "tsme", "pcgl","pcen", "pifl"))
  

## prism data
#prism_itrdb <- read_csv(here("Data", "raw_prism_data_itrdb.csv"))
prism_itrdb <- read_csv(here("Data", "raw_prism_itrdb_sevenDomSp.csv"))

growing_dat <- prism_itrdb %>%
  select(waypoint_id, variable,value, month, year, longitude, latitude) %>% 
  rename(collection_id = waypoint_id) %>% 
  mutate(growing=ifelse(month%in%c(10:12), year+1, year)) %>% 
  filter(!month%in%c(7:9)) %>% ## don't include July, August, September
  filter(growing!=1900) %>% 
  pivot_wider(names_from=variable, values_from=value) %>% 
  group_by(collection_id, growing) %>% 
  summarize(tmax = mean(tmax, na.rm=T), ppt = sum(ppt, na.rm=T))%>%
  rename(year=growing) %>% 
  ungroup() 


prism_seasonSummerJune_S <- prism_itrdb %>%
  select(waypoint_id, variable,value, month, year, longitude, latitude) %>% 
  rename(collection_id = waypoint_id) %>% 
  mutate(month=as.numeric(month),
         growing=ifelse(month%in%c(10:12), year+1, year)) %>%
  filter(month%in%c(6:9)) %>%
  filter(growing!=1900) %>%
  pivot_wider(names_from=variable, values_from=value) %>% 
  group_by(collection_id, growing) %>% 
  summarize(tmaxSummer = mean(tmax, na.rm=T),
            pptSummer = sum(ppt, na.rm=T))%>%
  rename(year=growing) %>%
  ungroup()

prism_season_ann <- prism_itrdb %>%
  select(waypoint_id, variable,value, month, year, longitude, latitude) %>% 
  rename(collection_id = waypoint_id) %>% 
  filter(year!=1900) %>%
  pivot_wider(names_from=variable, values_from=value) %>% 
  group_by(collection_id, year) %>% 
  summarize(tmax_an = mean(tmax, na.rm=T),
            ppt_an = sum(ppt, na.rm=T))%>%
  ungroup()


prism_lagged <- prism_itrdb %>%
  select(waypoint_id, variable,value, month, year, longitude, latitude) %>% 
  rename(collection_id = waypoint_id) %>% 
  mutate(month=as.numeric(month),
         growing=ifelse(month%in%c(10:12), year+1, year)) %>% 
  filter(!month%in%c(7:9)) %>% ## don't include July, August, September
  filter(growing!=1900) %>% 
  pivot_wider(names_from=variable, values_from=value) %>% 
  group_by(collection_id, growing) %>% 
  summarize(tmax = mean(tmax, na.rm=T),
            ppt = sum(ppt, na.rm=T))%>%
  rename(year=growing) %>% 
  ungroup() %>% 
  # Add lag variables for tmax and ppt by plot
  group_by(collection_id) %>%
  arrange(collection_id, year) %>%
  mutate(laggedtmax = lag(tmax, 1),
         laggedprecip = lag(ppt, 1)) %>%
  ungroup()


## tree age
pipo_age <- itrdbdat %>% 
  group_by(tree, collection_id) %>% 
  mutate(year = as.numeric(year)) %>% 
  mutate(maxyear= max(year), minyear=min(year)) %>% ## calculate age
  mutate(age = year - minyear + 1) %>% 
  filter(age<1000) %>% 
  filter(year > 1900) %>% 
  filter(!grepl("CAN", collection_id)) %>%
  filter(!grepl("MEX", collection_id)) %>% 
  left_join(growing_dat) %>% 
  left_join(prism_seasonSummerJune_S) %>% 
  left_join(prism_season_ann) %>% 
  left_join(prism_lagged)

hist(pipo_age$age, col = "blue")

#x = filter(pipo_age, collection_id == "NM559")

#write_csv(pipo_age, here("Data", "SevenSpecies_ITRDB_climatewindows.csv"))
