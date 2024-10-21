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
librarian::shelf(sjPlot, ggeffects, patchwork, tidyverse, broom,
                 lme4, plotrix, ggpubr, mgcv, nlme, fixest, plotrix, egg, ggpmisc,
                 mvtnorm, clubSandwich, rasterVis, broom.mixed, scales,RColorBrewer, splines, zoo)


## data
rwldat <- read_csv(here("Data", "paneldat_RWL.csv"))



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#                    COUNTERFACTUAL TEMPERATURE SCENARIO
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


temp_mod <- lm(tmax ~ ns(year, df = 4) + factor(plot_id_needle), data=rwldat)
summary(temp_mod)


## Visualizing the predicted data
preddat <- rwldat
preddat$predtmax <- predict(temp_mod)

preddat %>% 
  ggplot(aes(x=year, y=tmax, color="tmax"))+
  geom_line()+
  geom_line(aes(y=predtmax, color="predtmax"))+
  geom_line()

## Adding year-level heterogeneity back into the counterfactual values
meanpre80 <- preddat %>% 
  filter(year<1980) %>% 
  group_by(plot_id_needle) %>% 
  summarize(meantmax = mean(tmax)) 

mean_counterdat <- preddat %>% 
  left_join(meanpre80)

## time period 1980-2018
counterdat80 <- mean_counterdat %>% 
  mutate(tmax_counter = (tmax - predtmax) + meantmax) %>% 
  mutate(tmax_counter_historic = ifelse(year>1979, tmax_counter, tmax))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# FIGURE 6A
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#                   ESTIMATING THE DIFFERENT COUNTERFACTUAL GROWTH SCENARIOS
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
  select(plot_id_needle, tmax_counter_historic, ppt, tree_id, year) %>% 
  rename(tmax = tmax_counter_historic)

actualdat <- counterdat80 %>% 
  select(plot_id_needle, tmax, ppt, tree_id, year)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# MONTE CARLO SIMULATION
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## create an empty dataframe for results
results=data.frame()

for (i in 1:1000){
  
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


set.seed(123)  # for reproducibility
tercile_boundaries <- quantile(pred_results$ppt, probs = c(1/3, 2/3))

# First, let's prepare the data
calc_terc <- pred_results %>%
  sample_n(10000) %>% 
  mutate(precip_terciles = case_when(
    ppt < 0.667222 ~ "0-33.3%",
    ppt >= 0.667222 & ppt < 1.045141 ~ "33.4-66.7%",
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
ribbonplot = ggplot(prepared_data, aes(x = year, y = tree_diff, color = precip_terciles, fill = precip_terciles)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #geom_point(alpha = 0.01) +
  geom_ribbon(data = summarized_data, 
              aes(y = mean_diff, ymin = lower_ci, ymax = upper_ci),
              alpha = 0.2, color = NA) +
  geom_line(data = summarized_data, aes(y = mean_diff), size = 1) +
  scale_color_manual(values = c("0-33.3%" = "#303077", "33.4-66.7%" = "#A4BDD7", "67.1-100%" = "#8c6e87")) +
  scale_fill_manual(values = c("0-33.3%" = "#303077", "33.4-66.7%" = "#A4BDD7", "67.1-100%" = "#8c6e87")) +
  labs(y = "Δ log(RWL)", x = "Year", color = "Precip. terciles", fill = "Precip. terciles") +
  theme(legend.position = c(.2,.65))+
  theme(
    legend.background = element_rect(fill = "transparent"), 
    legend.key = element_rect(fill = "transparent", color = NA) 
  )

ribbonplot



precipterc <- pred_results %>% 
  mutate(precip_terciles = case_when(
    ppt < 0.667222 ~ "0-33.3%",
    ppt >= 0.667222 & ppt < 1.045141 ~ "33.4-66.7%",
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

## density figs
# dens_fig <- pred_results %>%
#   sample_n(1000) %>% 
#   mutate(precip_terciles = ifelse(ppt<0.662672, "0-33",
#                                   ifelse(ppt>=0.662672&ppt<1.034786 , "34-66", "67-100"))) %>%  
#   select(vals_actual, vals_counter, tree_id, iteration, plot_id_needle, period, precip_terciles) %>%
#   pivot_longer(-c(tree_id, iteration, plot_id_needle, period, precip_terciles)) %>%
#   filter(period == "2010-2018") %>% 
#   ggplot(aes(x=value, color=name, fill = name))+
#   geom_density(alpha=.3)+
#   labs(color = "", fill = "",
#        title = "Predicted RWL across precip. terciles",
#        x = "Predicted log(RWL)",
#        y = "Density")+
#   guides(fill = F, color = F)+
#   scale_fill_manual(
#     values = c("vals_counter" = "#303077", "vals_actual" = "#C75B77"),
#     labels = c("Actual", "Counterfactual"))+
#   scale_color_manual(
#     values = c("vals_counter" = "#303077", "vals_actual" = "#C75B77"),
#     labels = c("Actual", "Counterfactual"))+
#   theme(legend.position=c(.25,.90),
#         legend.background = element_blank())+
#   facet_grid(~precip_terciles)
# 
# 
# dens_fig



library(ggplot2)
library(grid)

# Create a simple ggplot
p <- ggplot() + geom_point() + ggtitle("Test Plot")

# Build the plot
built_plot <- ggplot_build(p)

# Extract the font family
default_font <- built_plot$plot$theme$text$family

# Print the result
print(paste("The default font family is:", default_font))


