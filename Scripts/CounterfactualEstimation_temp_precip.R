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

## data
# rwldat <- read_csv(here("Data", "Primary_data", "paneldat_RWL.csv"))[c(1:7)]
# climdat <- read_csv(here("Data", "Primary_data","climatewindows.csv"))
# water_year <- read_csv(here("Data", "Primary_data", "PANEL_data_rwi.csv"))[-c(3,9,10)]
# 
# paneldat <- rwldat %>% 
#   left_join(water_year) %>% 
#   na.omit()

paneldat <- read_csv(here("Data", "cleaned_RWL_climdat.csv"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#                    COUNTERFACTUAL TEMPERATURE SCENARIO
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## estimate changes in maximum temperature as a function of year for each plot
## use a flexible spline to improve predictive accuracy

## temperature (hereafter tmax)
temp_mod <- lm(tmax ~ ns(year, df = 4) + factor(plot_id_needle), data=paneldat)
summary(temp_mod)

## precipitation (hereafter precip)
precip_mod <- lm(ppt ~ ns(year, df = 4) + factor(plot_id_needle), data=paneldat)
summary(precip_mod)


## create a new dataset and add predicted temperature
preddat <- paneldat
preddat$predtmax <- predict(temp_mod)
preddat$predppt <- predict(precip_mod)

## visualize the predicted data compared to observed data

## temp
preddat %>% 
  ggplot(aes(x=year, y=tmax, color="tmax"))+
  geom_line()+
  geom_line(aes(y=predtmax, color="predtmax"))+
  geom_line()

## precip
preddat %>% 
  ggplot(aes(x=year, y=ppt, color="precip"))+
  geom_line()+
  geom_line(aes(y=predppt, color="predprecip"))+
  geom_line()

## calculate plot-level mean tmax
meanpre80 <- preddat %>% 
  filter(year<1980) %>% 
  group_by(plot_id_needle) %>% 
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
  mutate(ppt_counter_historic = ifelse(year>1999, ppt_counter, ppt))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# FIGURE 6A
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## visualize difference between actual and counterfactual scenarios
counter_clim_fig <- counterdat80 %>% 
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
  theme(legend.position=c(.70,.90),
        legend.background = element_blank())

counter_clim_fig



## precipitation figure
counterdat80 %>% 
  group_by(year) %>% 
  summarize(meanppt = mean(ppt), meancounter = mean(ppt_counter_historic)) %>% 
  select(year, meanppt, meancounter) %>% 
  pivot_longer(-year) %>% 
  ggplot(aes(x=year, y=value, color=name))+
  geom_line(linewidth = .5, alpha=.9)+
  geom_point(size = .9, alpha=.9)+
  geom_smooth(se=F, span=1)+
  labs(color = "",
       x = "Year",
       y = "Precipitation (mm)")+
  scale_color_manual(
    values = c( "meanppt" = "#C75B77" , "meancounter" = "#303077"),
    labels = c("Counterfactual", "Actual"))+
  theme(legend.position=c(.70,.90),
        legend.background = element_blank())



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#                   ESTIMATING THE DIFFERENT COUNTERFACTUAL GROWTH SCENARIOS
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## rename data
paneldat <- counterdat80

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
counterdat <- counterdat80 %>% 
  select(plot_id_needle, tmax_counter_historic, ppt_counter_historic, tree_id, year) %>% 
  rename(tmax = tmax_counter_historic,
         ppt = ppt_counter_historic)

actualdat <- counterdat80 %>% 
  select(plot_id_needle, tmax, ppt, tree_id, year)



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
  actualdat$vals_counter <-  predict(modified_fe_mod, newdata = counterdat)
  
  results_dat_1980 <- actualdat %>% 
    filter(year>1979) %>% 
    mutate(tree_diff = vals_actual - vals_counter) %>% 
    mutate(iteration = i, period = "1980-2018")
  
  results_dat_2000 <- actualdat %>% 
    filter(year>1999) %>% 
    mutate(tree_diff = vals_actual - vals_counter) %>% 
    mutate(iteration = i, period = "2000-2018")
  
  results_dat_2010 <- actualdat %>% 
    filter(year>2009) %>% 
    mutate(tree_diff = vals_actual - vals_counter) %>% 
    mutate(iteration = i, period = "2010-2018")
  
  results <- rbind(results, results_dat_1980, results_dat_2000, results_dat_2010)
}

head(results)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# THE FIGURES
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## new dataframe to create figures from
pred_results <- results
#write_csv(pred_results, here("Data", "Primary_data", "MC_simulation_data.csv"))

set.seed(123)  # for reproducibility
tercile_boundaries <- quantile(pred_results$ppt, probs = c(1/3, 2/3))

# First, let's prepare the data
# calc_terc <- pred_results %>%
#   sample_n(10000) %>% 
#   mutate(precip_terciles = case_when(
#     ppt < 0.667222 ~ "0-33.3%",
#     ppt >= 0.667222 & ppt < 1.045141 ~ "33.4-66.7%",
#     TRUE ~ "67.1-100%"
#   ))

##   667.222  1045.141 
calc_terc <- pred_results %>%
  sample_n(10000) %>% 
  mutate(precip_terciles = case_when(
    ppt < 667.222 ~ "0-33.3%",
    ppt >= 667.222 & ppt < 1045.141 ~ "33.4-66.7%",
    TRUE ~ "67.1-100%"
  ))

# Calculate the mean and 95% CI for each year and tercile
summarized_data <- calc_terc%>%
  group_by(year, precip_terciles) %>%
  summarise(
    mean_diff = mean(tree_diff),
    se = sd(tree_diff) / sqrt(n()),
    lower_ci = mean_diff - 1.96 * se,
    upper_ci = mean_diff + 1.96 * se,
    .groups = "drop"
  )

# Figure
ribbonplot = ggplot(summarized_data, aes(x = year, y = tree_diff, color = precip_terciles, fill = precip_terciles)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #geom_point(alpha = 0.01) +
  geom_ribbon(data = summarized_data, 
              aes(y = mean_diff, ymin = lower_ci, ymax = upper_ci),
              alpha = 0.2, color = NA) +
  geom_line(data = summarized_data, aes(y = mean_diff), size = 1) +
  scale_color_manual(values = c("0-33.3%" = "#303077", "33.4-66.7%" = "#A4BDD7", "67.1-100%" = "#8c6e87")) +
  scale_fill_manual(values = c("0-33.3%" = "#303077", "33.4-66.7%" = "#A4BDD7", "67.1-100%" = "#8c6e87")) +
  labs(y = "Δ log(RWL)", x = "Year", color = "Precip. terciles", fill = "Precip. terciles") +
  #theme(legend.position = c(.2,.65))+
  theme(
    legend.background = element_rect(fill = "transparent"), 
    legend.key = element_rect(fill = "transparent", color = NA) 
  )

ribbonplot



# precipterc <- pred_results %>% 
#   mutate(precip_terciles = case_when(
#     ppt < 0.667222 ~ "0-33.3%",
#     ppt >= 0.667222 & ppt < 1.045141 ~ "33.4-66.7%",
#     TRUE ~ "67.1-100%")) %>% 
#   group_by(period, precip_terciles) %>% 
#   summarize(meandiff = mean(tree_diff),lowCI = quantile(tree_diff, probs = .025), highCI = quantile(tree_diff, probs = 0.95)) %>% 
#   ungroup()

precipterc <- pred_results %>% 
  mutate(precip_terciles = case_when(
    ppt < 667.222 ~ "0-33.3%",
    ppt >= 667.222 & ppt < 1045.141 ~ "33.4-66.7%",
    TRUE ~ "67.1-100%")) %>% 
  group_by(period, precip_terciles) %>% 
  summarize(meandiff = mean(tree_diff),lowCI = quantile(tree_diff, probs = .025), highCI = quantile(tree_diff, probs = 0.95)) %>% 
  ungroup()



# Define the color palette
color_palette <- c("0-33.3%" = "#303077", "33.4-66.7%" = "#A4BDD7", "67.1-100%" = "#8c6e87")

diff_fig <- precipterc %>% 
  ggplot(aes(x = factor(period), y = meandiff, color = precip_terciles, fill = precip_terciles)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = lowCI, ymax = highCI), width = 0.1) +
  facet_grid(~precip_terciles) +
  labs(title = "Estimated climate change impacts",
       x = "Time period",
       y = "Δ log(RWL)") +
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values = color_palette) +
  guides(fill=F, color=F)+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 45, hjust = .4, vjust = 0.5),
        plot.title = element_text(hjust = 0.5))

diff_fig



theme_set(
  theme_bw(base_size = 20)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)


counter_clim_fig + (ribbonplot/diff_fig) & plot_annotation(tag_levels = "A")




# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# #
# #
# #                                   Without precip counterfactual
# #
# #
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
# 
# ## rename data
# paneldat <- counterdat80
# 
# ## extracting coefficient matrix
# fe_mod <-  feols(log(value+0.01) ~ tmax * ppt + year | tree_id,
#                  data= paneldat, cluster = ~ plot_id_needle)
# 
# summary(fe_mod)
# 
# V_CR1 = vcov(fe_mod)
# V_CR1=as.matrix(V_CR1)
# coef_vector = fe_mod$coefficients
# 
# 
# ## randomly drawing coefficients from distribution
# draw = rmvnorm(n = 1000, mean = coef_vector, sigma = V_CR1)
# 
# 
# ## creating counterfactual and actual dataframes
# counterdat_new <- counterdat80 %>% 
#   select(plot_id_needle, tmax_counter_historic, ppt, tree_id, year) %>% 
#   rename(tmax = tmax_counter_historic)
# 
# actualdat_new <- counterdat80 %>% 
#   select(plot_id_needle, tmax, ppt, tree_id, year)
# 
# 
# 
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # MONTE CARLO SIMULATION
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# ## create an empty dataframe for results
# results_new=data.frame()
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
#   actualdat_new$vals_actual <-  predict(modified_fe_mod, newdata = actualdat_new)
#   actualdat_new$vals_counter <-  predict(modified_fe_mod, newdata = counterdat_new)
#   
#   results_dat_1980 <- actualdat_new %>% 
#     filter(year>1979) %>% 
#     mutate(tree_diff = vals_actual - vals_counter) %>% 
#     mutate(iteration = i, period = "1980-2018")
#   
#   results_dat_2000 <- actualdat_new %>% 
#     filter(year>1999) %>% 
#     mutate(tree_diff = vals_actual - vals_counter) %>% 
#     mutate(iteration = i, period = "2000-2018")
#   
#   results_dat_2010 <- actualdat_new %>% 
#     filter(year>2009) %>% 
#     mutate(tree_diff = vals_actual - vals_counter) %>% 
#     mutate(iteration = i, period = "2010-2018")
#   
#   results_new <- rbind(results_new, results_dat_1980, results_dat_2000, results_dat_2010)
# }
# 
# head(results_new)
# 
# 
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # THE FIGURES
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# ## new dataframe to create figures from
# pred_results_new <- results_new
# 
# 
# set.seed(123)  # for reproducibility
# tercile_boundaries <- quantile(pred_results_new$ppt, probs = c(1/3, 2/3))
# ##   667.222  1045.141 
# 
# 
# calc_terc_new <- pred_results_new %>%
#   sample_n(10000) %>% 
#   mutate(precip_terciles = case_when(
#     ppt < 667.222 ~ "0-33.3%",
#     ppt >= 667.222 & ppt < 1045.141 ~ "33.4-66.7%",
#     TRUE ~ "67.1-100%"
#   ))
# 
# # Calculate the mean and 95% CI for each year and tercile
# summarized_data_new <- calc_terc_new%>%
#   group_by(year, precip_terciles) %>%
#   summarise(
#     mean_diff = mean(tree_diff),
#     se = sd(tree_diff) / sqrt(n()),
#     lower_ci = mean_diff - 1.96 * se,
#     upper_ci = mean_diff + 1.96 * se,
#     .groups = "drop"
#   )
# 
# # Figure
# ribbonplot_new = ggplot(summarized_data_new, aes(x = year, y = tree_diff, color = precip_terciles, fill = precip_terciles)) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   #geom_point(alpha = 0.01) +
#   geom_ribbon(data = summarized_data, 
#               aes(y = mean_diff, ymin = lower_ci, ymax = upper_ci),
#               alpha = 0.2, color = NA) +
#   geom_line(data = summarized_data, aes(y = mean_diff), size = 1) +
#   scale_color_manual(values = c("0-33.3%" = "#303077", "33.4-66.7%" = "#A4BDD7", "67.1-100%" = "#8c6e87")) +
#   scale_fill_manual(values = c("0-33.3%" = "#303077", "33.4-66.7%" = "#A4BDD7", "67.1-100%" = "#8c6e87")) +
#   labs(y = "Δ log(RWL)", x = "Year", color = "Precip. terciles", fill = "Precip. terciles") +
#   #theme(legend.position = c(.2,.65))+
#   theme(
#     legend.background = element_rect(fill = "transparent"), 
#     legend.key = element_rect(fill = "transparent", color = NA) 
#   )
# 
# ribbonplot_new
# 
# 
# precipterc_new <- pred_results_new %>% 
#   mutate(precip_terciles = case_when(
#     ppt < 667.222 ~ "0-33.3%",
#     ppt >= 667.222 & ppt < 1045.141 ~ "33.4-66.7%",
#     TRUE ~ "67.1-100%")) %>% 
#   group_by(period, precip_terciles) %>% 
#   summarize(meandiff = mean(tree_diff),lowCI = quantile(tree_diff, probs = .025), highCI = quantile(tree_diff, probs = 0.95)) %>%
#   mutate(Data = "Observed precip") %>% 
#   ungroup()
# 
# 
# ## combine both counterfactual estimates
# precipterc_all <- precipterc %>% 
#   mutate(Data = "Detrended precip") %>% 
#   rbind(precipterc_new)
# 
# precipterc_all %>% 
#   ggplot(aes(x = factor(period), y = meandiff, color = Data, fill = Data)) +
#   geom_point(position = position_dodge(width = 0.5), size = 3) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   geom_errorbar(aes(ymin = lowCI, ymax = highCI), width = 0.1, 
#                 position = position_dodge(width = 0.5)) +
#   facet_grid(~precip_terciles) +
#   labs(
#     x = "Time period",
#     y = "Δ log(RWL)") +
#   scale_color_manual(values = c("#303077", "#8c6e87")) +
#   scale_fill_manual(values = c("#303077", "#8c6e87")) +
#   #guides(fill=F, color=F)+
#   theme(axis.text.x = element_text(angle = 45, hjust = 0.4, vjust = 0.5),
#         plot.title = element_text(hjust = 0.5))
# 
# 
# 






