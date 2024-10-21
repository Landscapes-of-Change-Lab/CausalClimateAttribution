## ---------------------------
##
## Script name: Marginal Effects
##
## Author: Dr. Joan Dudney
##
## Date Created: 2024-10-20
##
## Copyright (c) Joan Dudney, 2024
## Email: dudney@ucsb.edu
##
## ---------------------------
##
## Notes: Estimates the fixed effects panel model of growth and
##   plots the marginal effects of climate variables
##
## ---------------------------


## packages
librarian::shelf(sjPlot, lme4, patchwork, tidyverse, ggeffects,
                 fixest, clubSandwich, here, lmtest, sandwich, marginaleffects)

theme_set(
  theme_bw(base_size = 15)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)


## data
rwldat <- read_csv(here("Data", "paneldat_RWL.csv"))
climdat <- read_csv(here("Data", "climatewindows.csv"))

paneldat <- rwldat %>% 
  left_join(climdat)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# The panel model 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fe_mod <-  feols(log(value+0.01) ~ tmax * ppt + year | tree_id,
                 data= paneldat, cluster = ~ plot_id_needle)

summary(fe_mod)
modelsummary(fe_mod)

tab_model(fe_mod, digits=4, show.ci = F, show.se = T)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Margins plot
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Compute temperature marginal effects
mtemp <- slopes(fe_mod, 
              variables = "tmax",
              by = "ppt",
              newdata = datagrid(ppt = seq(min(paneldat$ppt), max(paneldat$ppt), length.out = 100)))

# the plot
marg_temp = ggplot(mtemp, aes(x = ppt*1000, y = estimate, color)) +
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_line(color="#91376F") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color="#91376F", fill="#91376F") +
  labs(x = "Precipitation (mm)", 
       y = "Marginal effect of temperature")+
  scale_x_continuous(breaks = seq(0, 3800, by = 500))


# Compute marginal effects of precip
mprecip <- slopes(fe_mod, 
                variables = "ppt",
                by = "tmax",
                newdata = datagrid(tmax = seq(min(paneldat$tmax), max(paneldat$tmax), length.out = 100)))

# The plot
marg_precip = ggplot(mprecip, aes(x = tmax, y = estimate)) +
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_line(color= "#303077") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color="#303077", fill = "#303077") +
  labs(x = "Temperature (°C)", 
       y = "Marginal effect of precipitation")



marg_temp + marg_precip +   plot_annotation(tag_levels = "A")
 

