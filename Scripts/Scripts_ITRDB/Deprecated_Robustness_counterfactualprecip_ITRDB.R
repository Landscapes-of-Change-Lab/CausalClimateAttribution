## ---------------------------
##
## Script name: Robustness check: counterfactual precip
##
## Author: Dr. Joan Dudney
##
## Date Created: 2025-05-06
##
## Copyright (c) Joan Dudney, 2025
## Email: dudney@ucsb.edu
##
## ---------------------------
##
## Notes: This code runs a MC simulation to compare
##   counterfactual scenarios with and without precipitation
##
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
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#                    COUNTERFACTUAL SCENARIOS
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(tidyverse)
library(mgcv)  # For GAM models
library(splines)  # For spline models

# Prepare data similar to your example
climdat_short <- climdat %>% 
  filter(year > 1979)

clim_means <- climdat %>%
  filter(year < 1979) %>%
  group_by(plot) %>%
  summarize(hist_mtmax = mean(tmax),
            hist_mppt = mean(ppt))

climdat_short <- climdat_short %>%
  left_join(clim_means, by = "plot") %>%
  mutate(tmax_anom = tmax - hist_mtmax,
         plus_year = year - 1980)

# Create empty columns to store predictions
climdat_short <- climdat_short %>%
  mutate(
    tmax_pred_anom_linear = NA_real_,
    tmax_pred_anom_spline = NA_real_,
    tmax_pred_anom_gam = NA_real_,
    tmax_pred_anom_quad = NA_real_
  )

#========== TEMPERATURE MODELS ==========

# 1. Linear model (as in your example)
temp_mod_linear <- lm(tmax_anom ~ plus_year + 0, data = climdat_short)
summary(temp_mod_linear)
climdat_short$tmax_pred_anom_linear <- predict(temp_mod_linear)

# 2. Spline model
temp_mod_spline <- lm(tmax_anom ~ bs(plus_year, knots = c(10, 20, 30)) + 0, data = climdat_short)
summary(temp_mod_spline)
climdat_short$tmax_pred_anom_spline <- predict(temp_mod_spline)

# 3. GAM model
temp_mod_gam <- gam(tmax_anom ~ s(plus_year) + 0, data = climdat_short)
summary(temp_mod_gam)
climdat_short$tmax_pred_anom_gam <- predict(temp_mod_gam)

# 4. Quadratic model
temp_mod_quad <- lm(tmax_anom ~ plus_year + I(plus_year^2) + 0, data = climdat_short)
summary(temp_mod_quad)
climdat_short$tmax_pred_anom_quad <- predict(temp_mod_quad)


#========== CREATE COUNTERFACTUALS ==========

# Generate counterfactual climate variables for each model type
climdat_short <- climdat_short %>%
  mutate(
    tmax_cf_linear = tmax - tmax_pred_anom_linear,
    tmax_cf_spline = tmax - tmax_pred_anom_spline,
    tmax_cf_gam = tmax - tmax_pred_anom_gam,
    tmax_cf_quad = tmax - tmax_pred_anom_quad,
  )

#========== VISUALIZE RESULTS ==========

# Plot temperature anomaly models
ggplot(climdat_short, aes(x = year)) +
  geom_point(aes(y = tmax_anom), alpha = 0.3) +
  geom_line(aes(y = tmax_pred_anom_linear, color = "Linear"), size = 1) +
  geom_line(aes(y = tmax_pred_anom_spline, color = "Spline"), size = 1) + 
  geom_line(aes(y = tmax_pred_anom_gam, color = "GAM"), size = 1) +
  geom_line(aes(y = tmax_pred_anom_quad, color = "Quadratic"), size = 1) +
  labs(title = "Temperature Anomaly Models",
       y = "Temperature Anomaly (°C)",
       color = "Model Type") +
  theme_minimal()

# Compare counterfactual vs actual
ggplot(climdat_short, aes(x = year)) +
  geom_line(aes(y = tmax, color = "Observed"), size = 1) +
  geom_line(aes(y = tmax_cf_linear, color = "Linear CF"), size = 1, linetype = "dashed") +
  geom_line(aes(y = tmax_cf_spline, color = "Spline CF"), size = 1, linetype = "dashed") +
  geom_line(aes(y = tmax_cf_gam, color = "GAM CF"), size = 1, linetype = "dashed") +
  geom_line(aes(y = tmax_cf_quad, color = "Quad CF"), size = 1, linetype = "dashed") +
  labs(title = "Observed vs Counterfactual Temperature",
       y = "Temperature (°C)",
       color = "Series") +
  theme_minimal()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#                   ESTIMATING THE DIFFERENT COUNTERFACTUAL GROWTH SCENARIOS
#
#
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
actualdat <- climdat_short %>% 
  filter(year>1979 & year < 2014) %>% 
  select(plot, tmax, ppt, year) %>% 
  left_join(treedat)

counter_gam <- climdat_short %>% 
  filter(year>1979 & year < 2014) %>% 
  select(plot, tmax_cf_gam, ppt, year) %>% 
  rename(tmax = tmax_cf_gam) %>% 
  left_join(treedat)

counter_linear <- climdat_short %>% 
  filter(year>1979 & year < 2014) %>% 
  select(plot, tmax_cf_linear, ppt, year) %>% 
  rename(tmax = tmax_cf_linear) %>% 
  left_join(treedat)

counter_spline <- climdat_short %>% 
  filter(year>1979 & year < 2014) %>% 
  select(plot, tmax_cf_spline, ppt, year) %>% 
  rename(tmax = tmax_cf_spline) %>% 
  left_join(treedat)

counter_qu <- climdat_short %>% 
  filter(year>1979 & year < 2014) %>% 
  select(plot, tmax_cf_quad, ppt, year) %>% 
  rename(tmax = tmax_cf_quad) %>% 
  left_join(treedat)

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

for (i in 1:5){
  
  pb$tick()
  
  d <- draw[i,]
  modified_fe_mod <-  fe_mod
  modified_fe_mod$coefficients <- d
  actualdat$vals_actual <-  predict(modified_fe_mod, newdata = actualdat)
  actualdat$vals_gam <-  predict(modified_fe_mod, newdata = counter_gam)
  actualdat$vals_spline <-  predict(modified_fe_mod, newdata = counter_spline)
  actualdat$vals_qu <-  predict(modified_fe_mod, newdata = counter_qu)
  actualdat$vals_linear <-  predict(modified_fe_mod, newdata = counter_linear)
  
  results_dat_1980 <- actualdat %>% 
    filter(year >= 1981 & year <= 1991) %>%  
    mutate(tree_diff_gam = vals_actual - vals_gam) %>% 
    mutate(tree_diff_spline = vals_actual - vals_spline) %>% 
    mutate(tree_diff_qu = vals_actual - vals_qu) %>% 
    mutate(tree_diff_linear = vals_actual - vals_linear) %>% 
    mutate(iteration = i, period = "1981-1991")
  
  results_dat_2000 <- actualdat %>% 
    filter(year >= 1992 & year <= 2002) %>% 
    mutate(tree_diff_gam = vals_actual - vals_gam) %>% 
    mutate(tree_diff_spline = vals_actual - vals_spline) %>% 
    mutate(tree_diff_qu = vals_actual - vals_qu) %>% 
    mutate(tree_diff_linear = vals_actual - vals_linear) %>% 
    mutate(iteration = i, period = "1992-2002")
  
  results_dat_2010 <- actualdat %>% 
    filter(year >= 2003 & year <= 2013) %>% 
    mutate(tree_diff_gam = vals_actual - vals_gam) %>% 
    mutate(tree_diff_spline = vals_actual - vals_spline) %>% 
    mutate(tree_diff_qu = vals_actual - vals_qu) %>% 
    mutate(tree_diff_linear = vals_actual - vals_linear) %>% 
    mutate(iteration = i, period = "2003-2013")
  
  results <- rbind(results, results_dat_1980, results_dat_2000, results_dat_2010)
}

head(results)







#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# THE FIGURES
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


# Calculate the mean and 95% CI for each year and tercile
average_data <- pred_results %>% 
  left_join(tercdat) %>% 
  group_by(period, iteration, precip_terciles) %>%
  summarise(
    linear_diff = mean(tree_diff_linear),
    qu_diff = mean(tree_diff_qu),
    spline_diff = mean(tree_diff_spline),
    gam_diff = mean(tree_diff_gam)) %>% 
  group_by(period, precip_terciles) %>%
  summarise(
    lower_ci_linear = quantile(linear_diff, probs=0.025),
    upper_ci_linear = quantile(linear_diff, probs=0.975),
    mean_linear = mean(linear_diff),
    lower_ci_gam = quantile(gam_diff, probs=0.025),
    upper_ci_gam = quantile(gam_diff, probs=0.975),
    mean_gam = mean(gam_diff),
    lower_ci_spline = quantile(spline_diff, probs=0.025),
    upper_ci_spline = quantile(spline_diff, probs=0.975),
    mean_spline = mean(spline_diff),
    lower_ci_qu = quantile(qu_diff, probs=0.025),
    upper_ci_qu = quantile(qu_diff, probs=0.975),
    mean_qu = mean(qu_diff)
  ) %>%
  pivot_longer(
    cols = starts_with("lower_") | starts_with("upper_") | starts_with("mean_"),  
    names_to = c(".value", "model"),
    names_pattern = "(lower_ci|upper_ci|mean)_(.*)"
  )



average_data %>% 
  filter(period != "1981-1991") %>% 
  mutate(model = factor(model, levels = c("linear", "qu", "gam", "spline"))) %>% 
  ggplot(aes(x = factor(period), y = mean, color = model, fill = model)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), 
                width = 0.1, 
                position = position_dodge(width = 0.5)) +  # Match dodge width with points
  facet_grid(~precip_terciles)+
  labs(title = "Estimated climate change impacts",
       x = "Time period",
       y = "Δ log(RWL)", color="Method", fill="Method") +
  scale_color_manual(values = c("#91376F","#303077","#44A894","#EAB94B", "#5e6988"),
                     labels = c("Linear","Quadratic", "GAM", 
                                "Spline"))+
  scale_fill_manual(values = c("#91376F","#303077","#44A894","#EAB94B", "#5e6988"),
                    labels = c("Linear","Quadratic", "GAM", 
                               "Spline"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 45, hjust = .4, vjust = 0.5),
        plot.title = element_text(hjust = 0.5))



























