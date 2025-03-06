## ---------------------------
##
## Script name: Counterfactual figure
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
## Notes: Reads in data to create the counterfactual figure
##   
##
## ---------------------------

librarian::shelf(here, tidyverse, mgcv, zoo)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1. Reading in and cleaning data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## counterfactual data
counter <- read_csv(here("Data","Primary_data", "extracted_counterfactual_data.csv")) %>% 
  rename(counter_tmax = tmax)

## observed prism data
prismdat <- read_csv(here("Data","Primary_data", "extracted_prism_data.csv")) %>% 
  rename(prism_tmax = tmax) %>% 
  rename(point_id = waypoint_id)

## combining datasets

allclimdat <- counter %>% 
  left_join(prismdat)
  

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2. Creating the counterfactual figure
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theme_set(
  theme_bw(base_size = 20)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)

counterfig <- allclimdat %>% 
  na.omit() %>% 
  group_by(year) %>% 
  summarize(meantmax = mean(prism_tmax), meancounter = mean(counter_tmax)) %>% 
  pivot_longer(-(year)) %>% 
  ggplot(aes(x=year, y=value, color=name))+
  geom_line(linewidth = .3, alpha=.6)+
  geom_point(size = .5, alpha=.6)+
  geom_smooth(se=F, span=1)+
  labs(title = "Counterfactual data", color = "",
       x = "Year",
       y = "Temperature (°C)")+
  scale_color_manual(
    values = c( "meantmax" = "#C75B77" , "meancounter" = "#303077"),
    labels = c("Counterfactual", "Observed"))+
  theme(legend.position=c(.3,.15),
        legend.background = element_blank())+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks = seq(1900, 2020, by = 20))

counterfig

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3. Using different splines to create the comparison figure
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

splinedat <- prismdat %>% 
  rename(tmax = prism_tmax) %>% 
  na.omit()

# GAM detrended
gam_model <- gam(tmax ~ s(year), data = splinedat)
splinedat$gam_trend <- predict(gam_model)
splinedat$gam_detrended <- splinedat$tmax - splinedat$gam_trend + mean(splinedat$tmax)

# 2. Linear Detrending method
linear_model <- lm(tmax ~ year, data = splinedat)
splinedat$linear_trend <- predict(linear_model)
splinedat$linear_detrended <- splinedat$tmax - splinedat$linear_trend + mean(splinedat$tmax)

# 3. Cubic Spline method
spline_fit <- smooth.spline(splinedat$year, splinedat$tmax, df = 20)
splinedat$spline_trend <- predict(spline_fit, splinedat$year)$y
splinedat$spline_detrended <- splinedat$tmax - splinedat$spline_trend + mean(splinedat$tmax)

# 4. Quadratic detrended
qu_model <- lm(tmax ~ year + I(year^2), data = splinedat)
splinedat$qu_trend <- predict(qu_model)
splinedat$qu_detrended <- splinedat$tmax - splinedat$qu_trend + mean(splinedat$tmax)

# Cleaning data
data_long <- splinedat %>%
  group_by(year) %>% 
  summarize(meant=mean(tmax), qu_trend = mean(qu_detrended), gamdtrend = mean(gam_detrended),
            meanspline = mean(spline_detrended), meanLtrend = mean(linear_detrended)) %>% 
  pivot_longer(-year)

# Create plot
splinefig <- data_long %>% 
  filter(year>1960&year<2018) %>% 
  ggplot(aes(x = year, y = value, color = name)) +
  geom_line(linewidth = .3, alpha=.6)+
  geom_point(size = .5, alpha=.9)+
  geom_smooth(se=F, span=1)+
  labs(title = "Detrending",
       x = "Year",
       y = "Temperature (°C)",
       color = "Method") +
  scale_color_manual(values = c("#91376F","#303077","#44A894","#EAB94B", "#5e6988"),
                     labels = c("GAM", "Linear", 
                                "Quadratic", 
                                "Observed","Spline"))+
  #guides(fill=F, color=F)+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks = seq(1960, 2020, by = 10))

splinefig


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4. Figure showing time periods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

temp_data <- prismdat %>% 
  na.omit() %>% 
  group_by(year) %>% 
  summarize(tmax = mean(prism_tmax)) %>% 
  filter(year<2019)

# Calculate means 
mean_1900_1980 <- temp_data %>%
  filter(year >= 1900, year<2019) %>%
  summarise(mean_temp = mean(tmax)) %>%
  pull(mean_temp)

mean_1980_2018 <- temp_data %>%
  filter(year >= 1980, year<2019) %>%
  summarise(mean_temp = mean(tmax)) %>%
  pull(mean_temp)

mean_1990_2018 <- temp_data %>%
  filter(year >= 1990, year<2019) %>%
  summarise(mean_temp = mean(tmax)) %>%
  pull(mean_temp)

mean_2000_2018 <- temp_data %>%
  filter(year >= 2000, year<2019) %>%
  summarise(mean_temp = mean(tmax)) %>%
  pull(mean_temp)


# Create the plot
windows <- temp_data %>%
  filter(year<2019) %>% 
  ggplot(aes(x = year, y = tmax)) +
    geom_line(color =  "#5e6988") +
    geom_point(size = .5, alpha=.9, color="#5e6988")+
    annotate("segment", x = 1900, xend = 1980, 
           y = mean_1900_1980, yend = mean_1900_1980,
           color = "#303077", linewidth = 1) +
    annotate("segment", x = 1980, xend = 2020, 
           y = mean_1980_2018, yend = mean_1980_2018,
           color = "#44A894", linewidth = 1) +
    annotate("segment", x = 1990, xend = 2020, 
           y = mean_1990_2018, yend = mean_1990_2018,
           color = "#91376F", linewidth = 1) +
    annotate("segment", x = 2000, xend = 2020, 
           y = mean_2000_2018, yend = mean_2000_2018,
           color = "#EAB94B", linewidth = 1) +
    labs(title = "Time periods",
      x = "Year",
      y = "Observed temperature (°C)")+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_x_continuous(breaks = seq(1900, 2020, by = 20))

windows
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5. Combining figures
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

counterfig + windows + splinefig  & plot_annotation(tag_levels = "A") 

