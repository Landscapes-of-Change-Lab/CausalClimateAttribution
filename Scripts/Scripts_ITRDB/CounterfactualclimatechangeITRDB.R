## ---------------------------
##
## Script name: Counterfactual estimate
##
## Author: Dr. Joan Dudney
##
## Date Created: 2024-08-02
##
## Copyright (c) Joan Dudney, 2024
## Email: dudney@ucsb.edu
##
## ---------------------------
##
## Notes:
##   This script estimates the counterfactual climate and growth scenarios
##    using a MC simulation
## ---------------------------


## packages
librarian::shelf(sjPlot, ggeffects, patchwork, tidyverse, broom,progress, here,
                 lme4, plotrix, ggpubr, mgcv, nlme, fixest, plotrix, egg, ggpmisc,
                 mvtnorm, clubSandwich, rasterVis, broom.mixed, scales,RColorBrewer, splines, zoo)


theme_set(
  theme_bw(base_size = 20)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)

select <- dplyr::select


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Importing and cleaning data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## importing the panel ITRDB data
paneldat_itrdb <- read_csv(here("Data", "SevenSpecies_ITRDB_climatewindows.csv"))

paneldat <- paneldat_itrdb %>%
  mutate(tree_id = paste0(collection_id, "_", tree)) %>% 
  filter(species_id == "pied") %>%
  select(-tree) %>% ## remove duplicated tree ids 
  rename(tree = tree_id, plot = collection_id)

## creating a dataset with unique tree and plot
treedat <- paneldat %>% 
  select(tree, plot) %>% 
  distinct()

## importing climate data for all itrdb data
## panel data above does not have all climate data for predictions
itrdbdat <- read_csv(here("Data", "ITRDB_species_latlon.csv")) %>%
  filter(species_id == "pied") 

## isolating sites for PIED
pied <- itrdbdat %>% 
  select(species_id, collection_id) %>% 
  distinct()

## prism data
prism_itrdb <- read_csv(here("Data", "raw_prism_itrdb_sevenDomSp.csv"))

pied_clim <- prism_itrdb %>%
  select(waypoint_id, variable,value, month, year, longitude, latitude) %>% 
  rename(collection_id = waypoint_id) %>%
  left_join(pied) %>% 
  filter(!is.na(species_id))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#                    COUNTERFACTUAL SCENARIOS
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## all climate data
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

## temperature model
temp_mod <- lm(tmax ~ ns(year, df = 4) + factor(plot), data=climdat)
summary(temp_mod)

## precipitation model
precip_mod <- lm(ppt ~ ns(year, df = 1) + factor(plot), data=climdat)
summary(precip_mod)

## create a new dataset and add predicted temperature
preddat <- climdat
preddat$predtmax <- predict(temp_mod)
preddat$predppt <- predict(precip_mod)

# ## visualize the predicted data compared to observed data
# 
# ## temp
# preddat %>% 
#   ggplot(aes(x=year, y=tmax, color="tmax"))+
#   geom_line()+
#   geom_line(aes(y=predtmax, color="predtmax"))+
#   geom_line()
# 
# ## precip
# preddat %>% 
#   ggplot(aes(x=year, y=ppt, color="precip"))+
#   geom_line()+
#   geom_line(aes(y=predppt, color="predprecip"))+
#   geom_line()

## calculate plot-level mean tmax
meanpre80 <- preddat %>% 
  filter(year<1980) %>% 
  group_by(plot) %>% 
  summarize(meantmax = mean(tmax), meanppt = mean(ppt)) 

## join plot means with main dataset
mean_counterdat <- preddat %>% 
  left_join(meanpre80)

## add year-level heterogeneity to the counterfactual scenario
## formula for each predicted tmax: (tmax - predtmax) + meantmax
## time period 1980-2018
counterdat80 <- mean_counterdat %>% 
  mutate(tmax_counter = (tmax - predtmax) + meantmax) %>% ## adding heterogeneity
  mutate(tmax_counter_historic = ifelse(year>1979, tmax_counter, tmax)) %>%  ## include counterfactual values after 1979 only
  mutate(ppt_counter = (ppt - predppt) + meanppt) %>% 
  mutate(ppt_counter_historic = ifelse(year>1979, ppt_counter, ppt))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#     MC SIMULATION TO ESTIMATE THE EFFECT OF CHANGES IN TEMPERATURE ON GROWTH
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## estimating the model and extracting coefficient matrix
fe_mod <-  feols(rwi ~ tmax * ppt | tree + year,
                 data= paneldat, cluster = ~ plot)

summary(fe_mod)

V_CR1 = vcov(fe_mod)
V_CR1=as.matrix(V_CR1)
coef_vector = fe_mod$coefficients


## randomly drawing coefficients from distribution
draw = rmvnorm(n = 1000, mean = coef_vector, sigma = V_CR1)


## creating counterfactual and actual dataframes

## include precip counter
counterdat <- counterdat80 %>%
  select(plot, tmax_counter_historic, ppt, year) %>%
  rename(tmax = tmax_counter_historic) %>%
  left_join(treedat)%>% 
  filter(year<2014)

actualdat <- counterdat80 %>%
  left_join(treedat) %>% 
  select(plot, tmax, ppt, year, tree) %>% 
  filter(year<2014)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# MONTE CARLO SIMULATION
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## create an empty dataframe for results
results=data.frame()

pb <- progress_bar$new(
  format = "  Processing [:bar] :percent eta: :eta",
  total = 1000,
  clear = FALSE,
  width = 60
)


for (i in 1:10){
  
  pb$tick()
  
  d <- draw[i,]
  modified_fe_mod <-  fe_mod
  modified_fe_mod$coefficients <- d
  actualdat$vals_actual <-  predict(modified_fe_mod, newdata = actualdat)
  actualdat$vals_counter <-  predict(modified_fe_mod, newdata = counterdat)
  
  results_dat_1980 <- actualdat %>% 
    filter(year >= 1981 & year <= 1991) %>%  
    mutate(tree_diff = vals_actual - vals_counter) %>% 
    mutate(iteration = i, period = "1981-1991")
  
  results_dat_2000 <- actualdat %>% 
    filter(year >= 1992 & year <= 2002) %>% 
    mutate(tree_diff = vals_actual - vals_counter) %>% 
    mutate(iteration = i, period = "1992-2002")
  
  results_dat_2010 <- actualdat %>% 
    filter(year >= 2003 & year <= 2013) %>% 
    mutate(tree_diff = vals_actual - vals_counter) %>% 
    mutate(iteration = i, period = "2003-2013")
  
  results <- rbind(results, results_dat_1980, results_dat_2000, results_dat_2010)
}

head(results)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Main results
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mainresults <- results

meanchange <- mainresults %>% 
  group_by(period, iteration) %>% 
  summarise(mean_diff = mean(tree_diff)) %>%
  group_by(period) %>% 
  summarize(lower_ci = round(quantile(mean_diff, 0.025), digits = 3),
            rwi_mean = round(mean(mean_diff), digits=3),
            upper_ci = round(quantile(mean_diff, 0.975), digits = 3))


perchange <- mainresults %>% 
  group_by(period, iteration) %>% 
  summarise(actual = mean(vals_actual), counter = mean(vals_counter), 
            perchange = (actual - counter)/counter*100) %>% 
  group_by(period) %>% 
  summarize(lower_ci = round(quantile(perchange, 0.025), digits = 3),
            rwi_mean = round(mean(perchange), digits=3),
            upper_ci = round(quantile(perchange, 0.975), digits = 3))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# THE FIGURES
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# FIGURE 6B-C
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## new dataframe to create figures from
pred_results <- results 

## calculate historical mean precip for each plot
tercdat <- climdat %>% 
  filter(year < 1979) %>% 
  group_by(plot) %>% 
  summarize(meanppt = mean(ppt))

tercile_boundaries <- quantile(tercdat$meanppt, probs = c(1/3, 2/3))

tercdat <- tercdat %>%
  mutate(precip_terciles = case_when(
    meanppt < tercile_boundaries[1] ~ "0-33.3%",
    meanppt >= tercile_boundaries[1] & meanppt < tercile_boundaries[2] ~ "33.4-66.7%",
    TRUE ~ "66.8-100%"
  ))

# Sample the dataset to be able to plot it
calc_terc <- pred_results %>%
  left_join(tercdat)


# Calculate the median and 95% CI for plotting 
summarized_data <- calc_terc%>%
  group_by(year, precip_terciles, iteration) %>%
  summarise(mean_diff = mean(tree_diff)) %>%
  group_by(year, precip_terciles) %>% 
  summarize(lower_ci = quantile(mean_diff, 0.025),
            rwi_mean = mean(mean_diff),
            upper_ci = quantile(mean_diff, 0.975))

# Figure
ribbonplot = summarized_data %>% 
  filter(precip_terciles!="33.4-66.7%") %>% 
  ggplot( aes(x = year, y = rwi_mean, color = precip_terciles, fill = precip_terciles)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #geom_point(alpha = 0.01) +
  geom_ribbon(
              aes(y = rwi_mean, ymin = lower_ci, ymax = upper_ci),
              alpha = 0.2, color = NA) +
  geom_line( aes(y = rwi_mean), size = 1) +
  scale_color_manual(values = c("0-33.3%" = "#303077", "33.4-66.7%" = "#A4BDD7", "66.8-100%" = "#8c6e87")) +
  scale_fill_manual(values = c("0-33.3%" = "#303077", "33.4-66.7%" = "#A4BDD7", "66.8-100%" = "#8c6e87")) +
  labs(y = "Δ RWI", x = "Year", color = "Precip. terciles", fill = "Precip. terciles") +
  theme(legend.position = c(-.00000001,.15))+
  theme(
    legend.position = c(0.22, 0.38),  
    legend.background = element_rect(fill = "transparent"), 
    legend.key = element_rect(fill = "transparent", color = NA)
  )

ribbonplot


precipterc <- pred_results %>% 
  left_join(tercdat) %>% 
  group_by(period, precip_terciles, iteration) %>% 
  summarise(mean_diff = mean(tree_diff)) %>%
  summarize(lower_ci = quantile(mean_diff, 0.025),
            rwi_mean = mean(mean_diff),
            upper_ci = quantile(mean_diff, 0.975)) %>% 
  ungroup()

# Define the color palette
color_palette <- c("0-33.3%" = "#303077", "33.4-66.7%" = "#A4BDD7", "66.8-100%" = "#8c6e87")

diff_fig <- precipterc %>% 
  ggplot(aes(x = factor(period), y = rwi_mean, color = precip_terciles, fill = precip_terciles)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci, width = 0.1, color = "black"))+
  geom_point(size = 1) +
  facet_grid(~precip_terciles) +
  labs(title = "Estimated climate change impacts",
       x = "Time period",
       y = "Δ RWI") +
  #ylim(-.07,.03)+
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values = color_palette) +
  guides(fill=F, color=F)+
  #theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
        #axis.text.x = element_text(angle = 45, hjust = .4, vjust = 0.5))

diff_fig


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# FIGURE 6A
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


## visualize difference between actual and counterfactual scenarios
counter_clim_fig <- counterdat80 %>% 
  filter(year<2015) %>% 
  group_by(year) %>% 
  summarize(meantmax = mean(tmax), meancounter = mean(tmax_counter_historic)) %>% 
  select(year, meantmax, meancounter) %>% 
  pivot_longer(-year) %>% 
  ggplot(aes(x=year, y=value, color=name))+
  geom_line(linewidth = .5, alpha=.9)+
  geom_point(size = .9, alpha=.9)+
  geom_smooth(se=F, span=1)+
  labs(color = "",
       x = "Year",
       y = "Temperature (°C)")+
  scale_color_manual(
    values = c( "meantmax" = "#C75B77" , "meancounter" = "#303077"),
    labels = c("Counterfactual", "Actual"))+
  theme(legend.position=c(.80,.98),
        legend.background = element_blank())

counter_clim_fig


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Combining panels
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theme_set(
  theme_bw(base_size = 20)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)


counter_clim_fig + (ribbonplot/diff_fig) & plot_annotation(tag_levels = "A")








