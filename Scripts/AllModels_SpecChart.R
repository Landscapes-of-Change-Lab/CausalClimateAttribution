## ---------------------------
##
## Script name: CC Panel Model
##
## Author: Dr. Joan Dudney
##
## Date Created: 2024-09-15
##
## Copyright (c) Joan Dudney, 2024
## Email: dudney@ucsb.edu
##
## ---------------------------
##
## Notes: This code estimates the panel model
##    which is used to estimate the causal effect
##    of temperature on growth
## ---------------------------



## packages
librarian::shelf(sjPlot, lme4, patchwork, tidyverse, 
                 fixest, clubSandwich, here, lmtest, sandwich, mgcv,
                 marginaleffects)

## note: need to fix lagged weather data

## data
rwldat <- read_csv(here("Data", "paneldat_RWL.csv"))
climdat <- read_csv(here("Data", "climatewindows.csv"))

paneldat <- rwldat %>% 
  left_join(climdat)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1. The panel model 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fe_mod <-  feols(log(value + 0.01) ~ tmax * ppt + year | tree_id,
                 data= paneldat, cluster = ~ plot_id_needle)

fe_mod <-  feols(log(value+1) ~ tmax * ppt + year| tree_id,
                      data= paneldat, cluster = ~ plot_id_needle)
summary(fe_mod)
# avg_slopes(fe_mod)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2. Without clustered standard errors 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fe_mod_noCE <-  feols(log(value+1) ~ tmax * ppt + year| tree_id,
                      data= paneldat, vcov = "iid")

summary(fe_mod_noCE)
tab_model(fe_mod, digits = 4, show.se = T, show.ci = F)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3. heteroskedasticity-robust standard errors
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fe_mod_hetero <-  feols(log(value+0.01) ~ tmax * ppt + year | tree_id,
                        data= paneldat, vcov = "hetero")

summary(fe_mod_hetero)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4. Panel model without year effects
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fe_mod_noyr <-  feols(log(value+0.01) ~ tmax * ppt | tree_id,
                      data= paneldat)

summary(fe_mod_noyr)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5. Panel model with quadratic terms
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fe_mod_quad <-  feols(log(value+0.01) ~ tmax + ppt + I(tmax^2) + I(ppt^2) + year | tree_id,
                      data= paneldat, cluster = ~ plot_id_needle)
summary(fe_mod_quad)
tab_model(fe_mod_quad, digits = 4, show.se = T, show.ci = F)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 6. Panel model without interaction effects 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fe_mod_no_int<-  feols(log(value+0.01) ~ tmax + ppt + year | tree_id,
                     data= paneldat, cluster = ~ plot_id_needle)

summary(fe_mod_no_int)
tab_model(fe_mod_no_int)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 7. Panel model with lagged effects
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fe_mod_lag <-  feols(log(value+0.01) ~ tmax*ppt + laggedtmax + laggedprecip + year | tree_id,
                     data= paneldat, cluster = ~ plot_id_needle)

summary(fe_mod_lag)
tab_model(fe_mod_lag, digits = 4, show.se = T, show.ci = F)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 8. Panel model with annual weather
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fe_mod_an <-  feols(log(value+0.01) ~ tmax_an * ppt_an + year | tree_id,
                    data= paneldat, cluster = ~ plot_id_needle)

summary(fe_mod_an)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 9. Panel model with summer weather
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fe_mod_summer <-  feols(log(value+0.01) ~ tmaxSummer * pptSummer + year | tree_id,
                        data= paneldat, cluster = ~ plot_id_needle)

summary(fe_mod_summer)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 10. LM model
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lm_mod <- lm(log(value+0.01) ~ tmax * ppt + year + factor(plot_id_needle),
             data = paneldat)

summary(lm_mod)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 11. LMM with random intercepts
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lmm_mod_int = lmer(log(value+0.01) ~ tmax *  ppt + year +  (1|tree_id/plot_id_needle), 
                   control = lmerControl(optimizer = "bobyqa",
                                         optCtrl = list(maxfun = 100000)), 
                   data = paneldat)

summary(lmm_mod_int)
tab_model(lmm_mod_int, show.se = T, show.ci = F, digits = 4)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 12. LMM with random slopes and intercepts
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lmm_mod_rslopes = lmer(log(value+0.01) ~ tmax +   ppt  + year +
                         (1 + tmax | plot_id_needle)+
                         (1 + ppt | plot_id_needle)+
                         (1 + year | plot_id_needle), 
                       control = lmerControl(optimizer = "bobyqa",
                                             optCtrl = list(maxfun = 100000)), 
                       data = paneldat)


summary(lmm_mod_rslopes)
tab_model(lmm_mod_rslopes, show.se = T, show.ci = F, digits = 4)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 13. Just random slopes
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## just random slopes
# lmm_mod_rslopes_only = lmer(log(value+0.01) ~ tmax +   ppt  + year +
#                               (0 + tmax | tree_id)+
#                               (0 + ppt | tree_id)+
#                               (0 + year | tree_id), 
#                             control = lmerControl(optimizer = "bobyqa"), 
#                             data = paneldat)


# Updated model with scaled variables and increased iterations
lmm_mod_rslopes_only <- lmer(log(value+0.01) ~ tmax + ppt + year +
                               (0 + tmax | tree_id) +
                               (0 + ppt | tree_id) +
                               (0 + year | tree_id), 
                             control = lmerControl(optimizer = "bobyqa",
                                                   optCtrl = list(maxfun = 100000)), 
                             data = paneldat)


tab_model(lmm_mod_rslopes_only, show.se = T, show.ci = F, digits = 4)




rm(rwldat, climdat)


length(ls(envir = .GlobalEnv))




