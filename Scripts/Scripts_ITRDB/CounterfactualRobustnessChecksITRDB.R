## ---------------------------
##
## Script name: Robustness checks of counterfactual scenarios
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
##   This script tests whether the main results are robust to different 
##    estimates of the counterfactual temperature scenarios
## ---------------------------


## packages
librarian::shelf(sjPlot, ggeffects, patchwork, tidyverse, broom, progress,
                 lme4, plotrix, ggpubr, mgcv, nlme, fixest, plotrix, egg, ggpmisc,
                 mvtnorm, clubSandwich, rasterVis, broom.mixed, scales, RColorBrewer, splines, zoo)

theme_set(
  theme_bw(base_size = 20)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1. Reading in data and cleaning data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

select <- dplyr::select

## importing the panel ITRDB data
paneldat_itrdb <- read_csv(here("Data", "SevenSpecies_ITRDB_climatewindows.csv"))

paneldat <- paneldat_itrdb %>%
  mutate(tree_id = paste0(collection_id, "_", tree)) %>% 
  filter(species_id == "pied") %>%
  select(-tree) %>% ## remove duplicated tree ids 
  rename(tree = tree_id, plot = collection_id)

## creating a dataset with tree and plot
treedat <- paneldat %>% 
  select(tree, plot) %>% 
  distinct()

## importing climate data for ITRDB sites
itrdbdat <- read_csv(here("Data", "ITRDB_species_latlon.csv")) %>%
  filter(species_id == "pied") 

## isolating sites for PIED
pied <- itrdbdat %>% 
  select(species_id, collection_id) %>% 
  distinct() %>% 
  rename(plot = collection_id)

## prism data
prism_itrdb <- read_csv(here("Data", "raw_prism_itrdb_sevenDomSp.csv"))

pied_clim <- prism_itrdb %>%
  select(waypoint_id, variable,value, month, year, longitude, latitude) %>% 
  rename(plot = waypoint_id) %>%
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

## estimate changes in maximum temperature as a function of year for each plot
## use a flexible spline to improve predictive accuracy

climdat <- pied_clim %>% 
  mutate(growing=ifelse(month%in%c(10:12), year+1, year)) %>% 
  filter(!month%in%c(7:9)) %>% ## don't include July, August, September
  filter(growing!=1900) %>% 
  pivot_wider(names_from=variable, values_from=value) %>% 
  group_by(plot, growing) %>% 
  summarize(tmax = mean(tmax, na.rm=T), ppt = sum(ppt, na.rm=T))%>%
  rename(year=growing) %>% 
  ungroup() 

## temperature 
temp_mod <- lm(tmax ~ ns(year, df = 4) + factor(plot), data=climdat)
summary(temp_mod)

## create a new dataset and add predicted temperature
preddat <- climdat
preddat$predtmax <- predict(temp_mod)

## calculate plot-level mean tmax
meanpre80 <- preddat %>% 
  filter(year<1980) %>% 
  group_by(plot) %>% 
  summarize(meantmax = mean(tmax)) 

meanpre80 <- preddat %>% 
  filter(year<1980) %>% 
  group_by(plot) %>% 
  summarize(meantmax = mean(tmax)) 

## join plot means with main dataset
mean_counterdat <- preddat %>% 
  left_join(meanpre80)

## add year-level heterogeneity to the counterfactual scenario
## formula for each predicted tmax: (tmax - predtmax) + meantmax
## time period 1980-2018
counterdat80 <- mean_counterdat %>% 
  mutate(tmax_counter = (tmax - predtmax) + meantmax) %>% ## adding heterogeneity
  mutate(tmax_counter_historic = ifelse(year>1979, tmax_counter, tmax))

check <- counterdat80 %>% 
  select(plot, year, tmax, tmax_counter_historic) 

check1 <- counterdat80 %>% 
  select(plot, year, tmax, tmax_counter_historic) %>% 
  rename(tmax_meanall = tmax_counter_historic)

check_full <- check %>% 
  left_join(check1) %>% 
  filter(year>1980)

check_full %>% 
  filter(year>1980) %>% 
  ggplot(aes(x = year, y = tmax_counter_historic))+
  #geom_point()+
  #geom_point(aes(x = year, y = tmax_meanall), color="blue")
  geom_smooth(aes(x = year, y = tmax_meanall), color = "blue")+
  geom_smooth(aes(x = year, y = tmax_counter_historic), color = "green")+
  geom_smooth(aes(x = year, y = tmax), color = "red")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## estimating different splines
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rwi = climdat

# GAM detrended
gam_model <- gam(tmax ~ s(year) + plot, data = rwi)
rwi$gam_trend <- predict(gam_model)
rwi$gam_detrended <- rwi$tmax - rwi$gam_trend + mean(rwi$tmax)

# 2. Linear Detrending method
linear_model <- lm(tmax ~ year + factor(plot), data = rwi)
rwi$linear_trend <- predict(linear_model)
rwi$linear_detrended <- rwi$tmax - rwi$linear_trend + mean(rwi$tmax)

# 3. Cubic Spline method
spline_fit <- smooth.spline(rwi$year, rwi$tmax, df = 20)
rwi$spline_trend <- predict(spline_fit, rwi$year)$y
rwi$spline_detrended <- rwi$tmax - rwi$spline_trend + mean(rwi$tmax)

# 4. Quadratic detrended
qu_model <- lm(tmax ~ year + I(year^2) + factor(plot), data = rwi)
rwi$qu_trend <- predict(qu_model)
rwi$qu_detrended <- rwi$tmax - rwi$qu_trend + mean(rwi$tmax)

# 5. Target method
temp_mod <- lm(tmax ~ ns(year, df = 4) + factor(plot), data=rwi)
rwi$target_trend <- predict(temp_mod)
rwi$target_detrend <- rwi$tmax - rwi$target_trend + mean(rwi$tmax)
rwi$target_detrended <- ifelse(rwi$year <= 1979, 
                            rwi$tmax, 
                            rwi$target_detrend)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#                                    Text
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Base counterfactual scenario function
create_counterfactual <- function(model, data_df, model_name) {
  
  # Create prediction dataset
  preddat <- data_df
  preddat$predtmax <- predict(model, newdata = preddat)
  
  # Calculate plot-level mean tmax
  meanpre80 <- preddat %>%
    filter(year<1980) %>% 
    group_by(plot) %>% 
    summarize(meantmax = mean(tmax))
  
  # Join plot means with main dataset
  mean_counterdat <- preddat %>% 
    left_join(meanpre80)
  
  # Add year-level heterogeneity to the counterfactual scenario
  # Formula for each predicted tmax: (tmax - predtmax) + meantmax
  counterdat <- mean_counterdat %>% 
    mutate(tmax_counter = (tmax - predtmax) + meantmax) %>%
    mutate(tmax_counter_historic = ifelse(year > 1979, tmax_counter, tmax),
           model_type = model_name)
  
  return(counterdat)
}

# 1. Natural Spline Model (as in original code)
ns_model <- lm(tmax ~ ns(year, df = 4) + factor(plot), data = rwi)
ns_counter <- create_counterfactual(ns_model, rwi, "Natural Spline")

# 2. GAM Model
gam_model <- gam(tmax ~ s(year) + factor(plot), data = rwi)
gam_counter <- create_counterfactual(gam_model, rwi, "GAM")

# 3. Linear Model
linear_model <- lm(tmax ~ year + factor(plot), data = rwi)
linear_counter <- create_counterfactual(linear_model, rwi, "Linear")

# 4. Quadratic Model
qu_model <- lm(tmax ~ year + I(year^2) + factor(plot), data = rwi)
qu_counter <- create_counterfactual(qu_model, rwi, "Quadratic")

# 5. Smooth Spline (requires different handling due to different predict interface)
spline_fit <- smooth.spline(rwi$year, rwi$tmax, df = 20)

# For spline_fit we need a different approach
preddat_spline <- rwi
preddat_spline$predtmax <- predict(spline_fit, x = rwi$year)$y

# Calculate plot-level mean tmax
meanpre80_spline <- preddat_spline %>%
  filter(year<1980) %>% 
  group_by(plot) %>% 
  summarize(meantmax = mean(tmax))

# Join plot means with main dataset
mean_counterdat_spline <- preddat_spline %>% 
  left_join(meanpre80_spline)

# Add year-level heterogeneity to the counterfactual scenario
spline_counter <- mean_counterdat_spline %>% 
  mutate(tmax_counter = (tmax - predtmax) + meantmax) %>%
  mutate(tmax_counter_historic = ifelse(year > 1979, tmax_counter, tmax),
         model_type = "Smooth Spline")

# Combine all counterfactual datasets
all_counters <- bind_rows(
  ns_counter,
  gam_counter,
  linear_counter,
  qu_counter,
  spline_counter
)

# Visualize the different counterfactual scenarios
# Calculate means by year and model type for clearer visualization
yearly_means <- all_counters %>%
  group_by(year, model_type) %>%
  summarize(
    mean_tmax = mean(tmax),
    mean_counter = mean(tmax_counter_historic)
  )

# Plot the results
ggplot(yearly_means, aes(x = year)) +
  geom_line(aes(y = mean_tmax), color = "black", linetype = "solid", size = 1) +
  geom_line(aes(y = mean_counter, color = model_type), linetype = "dashed", size = 1) +
  geom_vline(xintercept = 1979.5, linetype = "dotted", color = "gray50") +
  geom_smooth(aes(y = mean_counter, color = model_type), se=F)+
  geom_smooth(aes(y = mean_tmax), color = "black", linetype = "solid", size = 1) +
  labs(
    title = "Comparison of Temperature Counterfactuals by Model Type",
    subtitle = "Solid black line = actual temperatures; Dashed lines = counterfactual scenarios",
    x = "Year",
    y = "Temperature (tmax)",
    color = "Model Type"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Summary statistics of counterfactual differences (post-1979)
counter_summary <- all_counters %>%
  filter(year > 1979) %>%
  group_by(model_type) %>%
  summarize(
    mean_diff = mean(tmax - tmax_counter_historic),
    max_diff = max(tmax - tmax_counter_historic),
    min_diff = min(tmax - tmax_counter_historic),
    sd_diff = sd(tmax - tmax_counter_historic)
  )

print(counter_summary)

# Return the counterfactual datasets
return(all_counters)




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#                   ESTIMATING THE DIFFERENT COUNTERFACTUAL GROWTH SCENARIOS
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fe_mod <-  feols(rwi ~ tmax * ppt | tree + year,
                 data= paneldat, cluster = ~ plot)

summary(fe_mod)

V_CR1 = vcov(fe_mod)
V_CR1=as.matrix(V_CR1)
coef_vector = fe_mod$coefficients


## randomly drawing coefficients from distribution
draw = rmvnorm(n = 1000, mean = coef_vector, sigma = V_CR1)


## creating counterfactual and actual dataframes
actualdat <- rwi %>% 
  left_join(treedat)%>% 
  filter(year<2014)

counter_gam <- rwi %>% 
  select(-tmax) %>% 
  rename(tmax = gam_detrended) %>% 
  select(plot, year, ppt, tmax) %>% 
  left_join(treedat)%>% 
  filter(year<2014)

counter_linear <- rwi %>% 
  select(-tmax) %>% 
  rename(tmax = linear_detrended) %>% 
  select(plot, year, tmax, ppt)%>% 
  left_join(treedat)%>% 
  filter(year<2014)

# counter_spline <- counterdat80 %>%
#   select(-tmax) %>% 
#   rename(tmax =  tmax_counter_historic) %>%
#   select(plot, year, ppt, tmax)%>%
#   left_join(treedat)%>%
#   filter(year<2014)

# counter_spline <- rwi %>%
#   select(-tmax) %>%
#   rename(tmax =  target_detrended) %>%
#   select(plot, year, ppt, tmax)%>%
#   left_join(treedat)%>%
#   filter(year<2014)

counter_spline <- rwi %>%
  select(-tmax) %>%
  rename(tmax =  target_detrended) %>%
  select(plot, year, ppt, tmax)%>%
  left_join(treedat)%>%
  filter(year<2014)


counter_qu <- rwi %>% 
  select(-tmax) %>% 
  rename(tmax = qu_detrended) %>% 
  select(plot, year, ppt, tmax)%>% 
  left_join(treedat)%>% 
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


# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # MONTE CARLO SIMULATION
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# ## create an empty dataframe for results
# results=data.frame()
# 
# pb <- progress_bar$new(
#   format = "  Processing [:bar] :percent eta: :eta",
#   total = 1000,
#   clear = FALSE,
#   width = 60
# )
# 
# 
# for (i in 1:1000){
#   
#   pb$tick()
#   
#   d <- draw[i,]
#   modified_fe_mod <-  fe_mod
#   modified_fe_mod$coefficients <- d
#   actualdat$vals_actual <-  predict(modified_fe_mod, newdata = actualdat)
#   actualdat$vals_counter <-  predict(modified_fe_mod, newdata = counterdat)
#   
#   results_dat_1980 <- actualdat %>% 
#     filter(year >= 1982 & year <= 1992) %>%  
#     mutate(tree_diff = vals_actual - vals_counter) %>% 
#     mutate(iteration = i, period = "1982-1992")
#   
#   results_dat_2000 <- actualdat %>% 
#     filter(year >= 1993 & year <= 2003) %>% 
#     mutate(tree_diff = vals_actual - vals_counter) %>% 
#     mutate(iteration = i, period = "1993-2003")
#   
#   results_dat_2010 <- actualdat %>% 
#     filter(year >= 2004 & year <= 2014) %>% 
#     mutate(tree_diff = vals_actual - vals_counter) %>% 
#     mutate(iteration = i, period = "2004-2014")
#   
#   results <- rbind(results, results_dat_1980, results_dat_2000, results_dat_2010)
# }
# 
# head(results)







#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# THE FIGURES
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## new dataframe to create figures from
pred_results <- results
#write_csv(pred_results, here("Data", "Primary_data", "counterfactual_robustness.csv"))

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


# Calculate the mean and 95% CI for each year and tercile
average_data <- calc_terc %>%
  group_by(period, precip_terciles) %>%
  summarise(
    linear_diff = mean(tree_diff_linear),
    qu_diff = mean(tree_diff_qu),
    spline_diff = mean(tree_diff_spline),
    gam_diff = mean(tree_diff_gam)) %>% 
  pivot_longer(
    cols = -c(precip_terciles, period),   
    names_to = "model",                   
    values_to = "mean"                   
  ) %>% 
  mutate(model = str_remove(model, "_diff"))

confdata <- calc_terc%>%
  group_by(period, precip_terciles, iteration) %>%
  summarise(
    mean_diff_linear = mean(tree_diff_linear),
    mean_diff_gam = mean(tree_diff_gam), 
    mean_diff_spline = mean(tree_diff_spline),
    mean_diff_qu = mean(tree_diff_qu)
    ) %>% 
  group_by(period, precip_terciles) %>% 
  summarise(
    lower_ci_linear = quantile(mean_diff_linear, probs=0.025),
    upper_ci_linear = quantile(mean_diff_linear, probs=0.975),
    lower_ci_gam = quantile(mean_diff_gam, probs=0.025),
    upper_ci_gam = quantile(mean_diff_gam, probs=0.975),
    lower_ci_spline = quantile(mean_diff_spline, probs=0.025),
    upper_ci_spline = quantile(mean_diff_spline, probs=0.975),
    lower_ci_qu = quantile(mean_diff_qu, probs=0.025),
    upper_ci_qu = quantile(mean_diff_qu, probs=0.975)
  ) %>%
  pivot_longer(
    cols = starts_with("lower_ci_") | starts_with("upper_ci_"),  
    names_to = c(".value", "model"),
    names_pattern = "(lower|upper)_ci_(.*)")


results_means_ci <- average_data %>% 
  left_join(confdata) 

results_means_ci %>% 
  mutate(model = factor(model, levels = c("linear", "qu", "gam", "spline"))) %>% 
  ggplot(aes(x = factor(period), y = mean, color = model, fill = model)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = 0.1, 
                position = position_dodge(width = 0.5)) +  # Match dodge width with points
  facet_grid(~precip_terciles)+
  labs(title = "Estimated climate change impacts",
       x = "Time period",
       y = "Δ RWI", color="Method", fill="Method") +
  #scale_color_manual(values = color_palette) +
  #scale_fill_manual(values = color_palette) +
  #guides(fill=F, color=F)+
  scale_color_manual(values = c("#91376F","#303077","#44A894","#EAB94B", "#5e6988"),
                     labels = c("Linear","Quadratic", "GAM", 
                                "Spline"))+
  scale_fill_manual(values = c("#91376F","#303077","#44A894","#EAB94B", "#5e6988"),
                    labels = c("Linear","Quadratic", "GAM", 
                               "Spline"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5))



