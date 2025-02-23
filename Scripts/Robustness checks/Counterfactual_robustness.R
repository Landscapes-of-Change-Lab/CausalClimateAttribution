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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1. Reading in data and extracting counterfactual data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## RWL data
rwldat <- read_csv(here("Data", "Primary_data" ,"paneldat_RWL.csv"))

## sourcing function to process nc files
source(here("Scripts", "Functions", "extract_counterfactual_function.R"))

## whitebark pine plots
wb_plot <- read_csv(here("Data", "Primary_data", "plotlatlong.csv"))

## projecting plots
wb_sf <- st_as_sf(wb_plot , 
                      coords = c("longitude", "latitude"), 
                      crs = 4326)
wb_sf_vect <- vect(wb_sf)

## reading in counterfactual data
nc_files <- list.files(here("Data",  "Counterfactual"), 
                       pattern = "\\.nc$", 
                       full.names = TRUE)
raster_list <- rast(nc_files)

## extracting to wb plots
extracted_counter <- terra::extract(raster_list, wb_sf_vect)


# Organizing the data
all_data <- lapply(nc_files, function(file) {
  process_nc_file(file, wb_sf_vect)
}) %>%
  bind_rows()


## easy extraction using terra

final_data <- all_data %>%
  arrange(point_id, date) %>%
  # Convert temperature to Celsius 
  mutate(temperature = temperature - 273.15) %>%
  # Add spatial coordinates from waypoints
  left_join(
    data.frame(
      point_id = 1:nrow(wb_sf_vect),
      longitude = crds(wb_sf_vect)[,1],
      latitude = crds(wb_sf_vect)[,2]
    ),
    by = "point_id"
  )


## add in plot names
plotdat <- wb_plot %>% 
  mutate(point_id = 1:nrow(wb_sf_vect))

# summarize data
counterdat <- final_data %>% 
  group_by(year, point_id, longitude, latitude) %>% 
  summarize(counter_tmax=mean(temperature)) %>% 
  left_join(plotdat) %>% 
  select(-point_id) %>% 
  rename(plot_id_needle = plot)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## estimating different splines
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rwl = rwldat %>% 
  left_join(counterdat)


# GAM detrended
gam_model <- gam(tmax ~ s(year), data = rwl)
rwl$gam_trend <- predict(gam_model)
rwl$gam_detrended <- rwl$tmax - rwl$gam_trend + mean(rwl$tmax)

# 2. Linear Detrending method
linear_model <- lm(tmax ~ year, data = rwl)
rwl$linear_trend <- predict(linear_model)
rwl$linear_detrended <- rwl$tmax - rwl$linear_trend + mean(rwl$tmax)

# 3. Cubic Spline method
spline_fit <- smooth.spline(rwl$year, rwl$tmax, df = 20)
rwl$spline_trend <- predict(spline_fit, rwl$year)$y
rwl$spline_detrended <- rwl$tmax - rwl$spline_trend + mean(rwl$tmax)

# 4. Quadratic detrended
qu_model <- lm(tmax ~ year + I(year^2), data = rwl)
rwl$qu_trend <- predict(qu_model)
rwl$qu_detrended <- rwl$tmax - rwl$qu_trend + mean(rwl$tmax)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#                   ESTIMATING THE DIFFERENT COUNTERFACTUAL GROWTH SCENARIOS
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paneldat <- rwldat 


## extracting coefficient matrix
fe_mod <-  feols(log(value+0.01) ~ tmax * ppt + year | tree_id,
                 data= paneldat, cluster = ~ plot_id_needle)

summary(fe_mod)

V_CR1 = vcov(fe_mod)
V_CR1=as.matrix(V_CR1)
coef_vector = fe_mod$coefficients


## randomly drawing coefficients from distribution
draw = rmvnorm(n = 1000, mean = coef_vector, sigma = V_CR1)


## creating counterfactual and actual dataframes
actualdat <- paneldat

counter_gam <- rwl %>% 
  select(-tmax) %>% 
  rename(tmax = gam_detrended)

counter_linear <- rwl %>% 
  select(-tmax) %>% 
  rename(tmax = linear_detrended)

counter_spline <- rwl %>% 
  select(-tmax) %>% 
  rename(tmax =  spline_detrended)

counter_qu <- rwl %>% 
  select(-tmax) %>% 
  rename(tmax = qu_detrended)

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

for (i in 1:1000){
  
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
    filter(year>1979) %>% 
    mutate(tree_diff_gam = vals_actual - vals_gam) %>% 
    mutate(tree_diff_spline = vals_actual - vals_spline) %>% 
    mutate(tree_diff_qu = vals_actual - vals_qu) %>% 
    mutate(tree_diff_linear = vals_actual - vals_linear) %>% 
    mutate(iteration = i, period = "1980-2018")
  
  results_dat_2000 <- actualdat %>% 
    filter(year>1999) %>% 
    mutate(tree_diff_gam = vals_actual - vals_gam) %>% 
    mutate(tree_diff_spline = vals_actual - vals_spline) %>% 
    mutate(tree_diff_qu = vals_actual - vals_qu) %>% 
    mutate(tree_diff_linear = vals_actual - vals_linear) %>% 
    mutate(iteration = i, period = "2000-2018")
  
  results_dat_2010 <- actualdat %>% 
    filter(year>2009) %>% 
    mutate(tree_diff_gam = vals_actual - vals_gam) %>% 
    mutate(tree_diff_spline = vals_actual - vals_spline) %>% 
    mutate(tree_diff_qu = vals_actual - vals_qu) %>% 
    mutate(tree_diff_linear = vals_actual - vals_linear) %>% 
    mutate(iteration = i, period = "2010-2018")
  
  results <- rbind(results, results_dat_1980, results_dat_2000, results_dat_2010)
}

head(results)







#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# THE FIGURES
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## new dataframe to create figures from
pred_results <- results
#write_csv(pred_results, here("Data", "Primary_data", "counterfactual_robustness.csv"))


set.seed(123)  # for reproducibility
tercile_boundaries <- quantile(pred_results$ppt, probs = c(1/3, 2/3))

# First, let's prepare the data
calc_terc <- pred_results %>%
  #sample_n(10000) %>% 
  mutate(precip_terciles = case_when(
    ppt < 0.667222 ~ "0-33.3%",
    ppt >= 0.667222 & ppt < 1.045141 ~ "33.4-66.7%",
    TRUE ~ "67.1-100%"
  ))


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
  group_by(period, precip_terciles) %>%
  summarise(
    lower_ci_linear = quantile(tree_diff_linear, probs=0.025),
    upper_ci_linear = quantile(tree_diff_linear, probs=0.975),
    lower_ci_gam = quantile(tree_diff_gam, probs=0.025),
    upper_ci_gam = quantile(tree_diff_gam, probs=0.975),
    lower_ci_spline = quantile(tree_diff_spline, probs=0.025),
    upper_ci_spline = quantile(tree_diff_spline, probs=0.975),
    lower_ci_qu = quantile(tree_diff_qu, probs=0.025),
    upper_ci_qu = quantile(tree_diff_qu, probs=0.975)
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
       y = "Δ log(RWL)", color="Method", fill="Method") +
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
  theme(axis.text.x = element_text(angle = 45, hjust = .4, vjust = 0.5),
        plot.title = element_text(hjust = 0.5))


# results_means_ci %>% 
#   mutate(model = factor(model, levels = c("linear", "gam", "spline"))) %>% 
#   ggplot(aes(x = factor(period), y = mean, color = model, fill = model)) +
#   geom_point(size = 3,position = position_dodge(width = 0.5)) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, position = "dodge") +
#   facet_grid(~precip_terciles) +
#   labs(title = "Estimated climate change impacts",
#        x = "Time period",
#        y = "Δ log(RWL)") +
#   scale_color_manual(values = color_palette) +
#   scale_fill_manual(values = color_palette) +
#   guides(fill=F, color=F)+
#   theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
#   theme(axis.text.x = element_text(angle = 45, hjust = .4, vjust = 0.5),
#         plot.title = element_text(hjust = 0.5))


