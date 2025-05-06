## ---------------------------
##
## Script name: Effects of temp on RWI
##
## Author: Dr. Joan Dudney
##
## Date Created: 2024-10-23
##
## Copyright (c) Joan Dudney, 2024
## Email: dudney@ucsb.edu
##
## ---------------------------
##
## Notes: Estimates the effects of temperature
##   on RWI as a robustness check
##
## ---------------------------


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1. Packages and data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


## packages
librarian::shelf(sjPlot, ggeffects, patchwork, tidyverse, here,
                 lme4, plotrix, ggpubr, mgcv, nlme, fixest, plotrix, egg, ggpmisc,
                 mvtnorm, clubSandwich, rasterVis, broom.mixed, scales,RColorBrewer)

theme_set(
  theme_bw(base_size = 15)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)


rwidat <- read_csv(here("Data", "Primary_data", "PANEL_data_rwi.csv")) %>% 
  mutate(ppt = ppt/1000)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2. The model
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fe_mod <-  feols(value ~ tmax + ppt | tree_id + year,
                 data= rwidat, cluster = ~ plot_id_needle)

summary(fe_mod)
tab_model(fe_mod, digits = 4, show.ci = F, show.se = T)

