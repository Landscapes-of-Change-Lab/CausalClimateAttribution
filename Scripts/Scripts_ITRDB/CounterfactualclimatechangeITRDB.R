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
                 lme4, plotrix, ggpubr, nlme, fixest, plotrix, egg, ggpmisc,
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

## isolating sites of PIED
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
# Counterfactual climate scenarios
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#========== Subsetting and further cleaning data ==========

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

climdat_short <- climdat %>% 
  filter(year>1979)

clim_means <- climdat %>%
  filter(year < 1979) %>%
  group_by(plot) %>%
  summarize(hist_mtmax = mean(tmax),
            hist_mppt = mean(ppt))

climdat_short <- climdat_short %>%
  left_join(clim_means, by = "plot") %>%
  mutate(tmax_anom = tmax - hist_mtmax,
         ppt_anom = ppt - hist_mppt,
         plus_year = year - 1980)

#========== Climate models ==========

## temperature model
temp_mod <- lm(tmax_anom ~ plus_year + 0, data=climdat_short)
summary(temp_mod)

## precipitation model
precip_mod <- lm(ppt_anom ~ plus_year + 0, data=climdat_short)
summary(precip_mod)

## counterfactual climate data
climdat_short <- climdat_short %>%
  mutate(tmax_pred_anom = predict(temp_mod),
         ppt_pred_anom = predict(precip_mod),
         tmax_cf = tmax - tmax_pred_anom,
         ppt_cf = ppt - ppt_pred_anom)


#========== Counterfactual datasets ==========

## creating counterfactual and actual dataframes
counterdat <- climdat_short %>%
  select(plot, tmax_cf, ppt, year) %>%
  rename(tmax = tmax_cf) %>%
  left_join(treedat)%>% 
  filter(year<2014)

actualdat <- climdat_short%>%
  left_join(treedat) %>% 
  select(plot, tmax, ppt, year, tree) %>% 
  filter(year<2014)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
# MC simulation to estimate the effect of climate change on tree growth
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



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# MC simulation
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## create an empty dataframe for results
results=data.frame()

pb <- progress_bar$new(
  format = "  Processing [:bar] :percent eta: :eta",
  total = 1000,
  clear = FALSE,
  width = 60
)


for (i in 1:1000){
  
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
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
# Main results
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
#
#
# The figures
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#========== Ribbon plot figure  ==========

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

summarized_data <-  pred_results %>% 
  #filter(year>1990) %>% 
  left_join(tercdat) %>% 
  group_by(year, precip_terciles, iteration) %>%
  summarize(mean_diff = mean(tree_diff)) %>%
  group_by(year, precip_terciles) %>% 
  summarize(lower_ci = quantile(mean_diff, 0.025),
            rwi_mean = mean(mean_diff),
            upper_ci = quantile(mean_diff, 0.975))

# Figure
ribbonplot = ggplot(summarized_data, aes(x = year, y = tree_diff, color = precip_terciles, fill = precip_terciles)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #geom_point(alpha = 0.01) +
  geom_ribbon(data = summarized_data, 
              aes(y = rwi_mean, ymin = lower_ci, ymax = upper_ci),
              alpha = 0.2, color = NA) +
  geom_line(data = summarized_data, aes(y = rwi_mean), size = 1) +
  scale_color_manual(values = c("0-33.3%" = "#303077", "33.4-66.7%" = "#A4BDD7", "66.8-100%" = "#8c6e87")) +
  scale_fill_manual(values = c("0-33.3%" = "#303077", "33.4-66.7%" = "#A4BDD7", "66.8-100%" = "#8c6e87")) +
  labs(y = "Δ RWI", x = "Year", color = "Precip. terciles", fill = "Precip. terciles") +
  #theme(legend.position = c(.2,.65))+
  theme(legend.position=c(.2,.35),
    legend.background = element_rect(fill = "transparent"), 
    legend.key = element_rect(fill = "transparent", color = NA) 
  )

ribbonplot


#========== Mean change figure  ==========

## data
precipterc <- pred_results %>% 
  left_join(tercdat) %>% 
  group_by(period, precip_terciles, iteration) %>% 
  summarise(mean_diff = mean(tree_diff)) %>%
  summarize(lower_ci = quantile(mean_diff, 0.025),
            rwi_mean = mean(mean_diff),
            upper_ci = quantile(mean_diff, 0.975)) %>% 
  ungroup()

## Define the color palette
color_palette <- c("0-33.3%" = "#303077", "33.4-66.7%" = "#A4BDD7", "66.8-100%" = "#8c6e87")


## figure
diff_fig <- precipterc %>% 
  ggplot(aes(x = factor(period), y = rwi_mean, color = precip_terciles, fill = precip_terciles)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci, width = 0.07), color = "black")+
  geom_point(size = 1.5) +
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


#========== Temperature counterfactual figure ==========

## fig temp trends
## visualize difference between actual and counterfactual scenarios
allclimdat <- climdat_short %>% 
  select(plot, year, tmax_cf) %>% 
  right_join(climdat) %>% 
  mutate(tmax_cf_all = ifelse(year < 1980, tmax, tmax_cf)) %>% 
  filter(year<2017)

counter_clim_fig <- allclimdat %>% 
  group_by(year) %>% 
  summarize(meantmax = mean(tmax), meancounter = mean(tmax_cf_all)) %>% 
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
  theme(legend.position=c(.74,.95),
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



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#  Robustness check for precip counterfactual
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


## creating counterfactual and actual dataframes
counterdat <- climdat_short %>%
  select(plot, tmax_cf, ppt_cf, year) %>%
  rename(tmax = tmax_cf,
         ppt = ppt_cf) %>%
  left_join(treedat)%>% 
  filter(year<2014)

actualdat <- climdat_short%>%
  left_join(treedat) %>% 
  select(plot, tmax, ppt, year, tree) %>% 
  filter(year<2014)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# MC Simulation with precipitation
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## create an empty dataframe for results
results_precip=data.frame()

pb <- progress_bar$new(
  format = "  Processing [:bar] :percent eta: :eta",
  total = 1000,
  clear = FALSE,
  width = 60
)


for (i in 1:100){
  
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
  
  results_precip <- rbind(results_precip, results_dat_1980, results_dat_2000, results_dat_2010)
}

head(results_precip)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Calculations
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


meanchange_noprecip <- meanchange %>% 
  mutate(Scenario = "Without precip counterfactual")
  
  
precip_main <-  results_precip %>% 
  group_by(period, iteration) %>% 
  summarise(mean_diff = mean(tree_diff)) %>%
  group_by(period) %>% 
  summarize(lower_ci = round(quantile(mean_diff, 0.025), digits = 3),
            rwi_mean = round(mean(mean_diff), digits=3),
            upper_ci = round(quantile(mean_diff, 0.975), digits = 3)) %>% 
  mutate(Scenario = "With precip counterfactual")


comb_main_res <- precip_main %>% 
  full_join(meanchange_noprecip)


## figure
comb_main_res %>% 
  ggplot(aes(x = factor(period), y = rwi_mean, color = Scenario, fill = Scenario, shape = Scenario)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci, width = 0.07), 
                position = position_dodge(width = 0.5), color = "black") +
  geom_point(size = 3, position = position_dodge(width = 0.5))+ 
  labs(title = "Estimated climate change impacts",
       x = "Time period",
       y = "Δ RWI") +
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values = color_palette) 



color_palette <- c("Without precip counterfactual" = "#303077", "With precip counterfactual" = "#8c6e87")
