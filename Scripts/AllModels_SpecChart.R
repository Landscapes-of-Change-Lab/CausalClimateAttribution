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

source(here("Scripts", "Functions", "spec_chart_function.R"))


## note: need to fix lagged weather data


## data
# rwldat <- read_csv(here("Data", "Primary_data", "paneldat_RWL.csv"))
# climdat <- read_csv(here("Data", "Primary_data","climatewindows.csv"))
wdir <- 'G:/My Drive/collaborations/ecology/CausalInferenceClimateAttribution/'
rwldat <- read_csv(paste0(wdir, "paneldat_RWL.csv"))
climdat <- read_csv(paste0(wdir, "climatewindows.csv"))

paneldat <- rwldat %>% 
  left_join(climdat)


winsorize_value <- paneldat %>% filter(value>0) %>% pull(value) %>% min()
paneldat <- paneldat %>% 
  mutate(value_w = ifelse(value==0, winsorize_value, value))


# Function to extract coefficients and SE
extract_coef_se <- function(model) {
  slopes_df <- avg_slopes(model) %>% 
    as_tibble()
  slopes_df %>% 
    filter(term %in% c("tmax", "tmaxSummer", "tmax_an")) %>% 
    select(term, estimate, std.error)
}


specs <- data.frame(coef=NaN, 
                    se=NaN,
                    stder_clust = NaN,
                    stder_het = NaN,
                    clim_lin = NaN,
                    clim_quad = NaN,
                    clim_int = NaN,
                    ww_wy = NaN,
                    ww_an = NaN,
                    ww_sum = NaN,
                    ww_lag = NaN,
                    struc_plotfe = NaN,
                    struc_ri = NaN) 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1. The panel model 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fe_mod <-  feols(log(value_w) ~ tmax * ppt + year | tree_id,
                 data= paneldat, cluster = ~ plot_id_needle)
summary(fe_mod)
avg_slopes(fe_mod)

fe_coef <- extract_coef_se(fe_mod)
new_row <-  data_frame(coef=fe_coef$estimate, 
                       se=fe_coef$std.error,
                       stder_clust = TRUE,
                       stder_het = FALSE,
                       clim_lin = TRUE, 
                       clim_quad = FALSE,
                       clim_int = FALSE,
                       ww_wy = TRUE,
                       ww_an = FALSE,
                       ww_sum = FALSE,
                       ww_lag = FALSE,
                       struc_plotfe = TRUE,
                       struc_ri = FALSE)
specs <- rbind(specs, new_row)
# avg_slopes(fe_mod)

# 
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # 2. Without clustered standard errors 
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# fe_mod_noCE <-  feols(log(value_w) ~ tmax * ppt + year| tree_id,
#                       data= paneldat, vcov = "iid")
# 
# summary(fe_mod_noCE)
# tab_model(fe_mod, digits = 4, show.se = T, show.ci = F)
# 
# fe_coef <- extract_coef_se(fe_mod)
# new_row <-  data_frame(coef=fe_coef$estimate, 
#                        se=fe_coef$std.error,
#                        stder_clust = TRUE,
#                        stder_het = FALSE,
#                        clim_lin = TRUE, 
#                        clim_quad = FALSE,
#                        clim_int = FALSE,
#                        ww_wy = TRUE,
#                        ww_an = FALSE,
#                        ww_sum = FALSE,
#                        ww_lag = FALSE,
#                        struc_plotfe = TRUE,
#                        struc_ri = FALSE)
# specs <- rbind(specs, new_row)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3. heteroskedasticity-robust standard errors
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fe_mod_hetero <-  feols(log(value_w) ~ tmax * ppt + year | tree_id,
                        data= paneldat, vcov = "hetero")

summary(fe_mod_hetero)

fe_hetero_coef <- extract_coef_se(fe_mod_hetero)
new_row <-  data_frame(coef=fe_hetero_coef$estimate, 
                       se=fe_hetero_coef$std.error,
                       stder_clust = FALSE,
                       stder_het = TRUE,
                       clim_lin = TRUE, 
                       clim_quad = FALSE,
                       clim_int = FALSE,
                       ww_wy = TRUE,
                       ww_an = FALSE,
                       ww_sum = FALSE,
                       ww_lag = FALSE,
                       struc_plotfe = TRUE,
                       struc_ri = FALSE)
specs <- rbind(specs, new_row)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4. Panel model without year effects
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fe_mod_noyr <-  feols(log(value + 0.01) ~ tmax * ppt | tree_id,
                      data= paneldat)

summary(fe_mod_noyr)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5. Panel model with quadratic terms
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fe_mod_quad <-  feols(log(value + 0.01) ~ tmax + ppt + I(tmax^2) + I(ppt^2) + year | tree_id,
                      data= paneldat, cluster = ~ plot_id_needle)
summary(fe_mod_quad)
tab_model(fe_mod_quad, digits = 4, show.se = T, show.ci = F)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 6. Panel model without interaction effects 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fe_mod_no_int<-  feols(log(value + 0.01) ~ tmax + ppt + year | tree_id,
                     data= paneldat, cluster = ~ plot_id_needle)

summary(fe_mod_no_int)
tab_model(fe_mod_no_int)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 7. Panel model with lagged effects
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fe_mod_lag <-  feols(log(value + 0.01) ~ tmax*ppt + laggedtmax + laggedprecip + year | tree_id,
                     data= paneldat, cluster = ~ plot_id_needle)

summary(fe_mod_lag)
tab_model(fe_mod_lag, digits = 4, show.se = T, show.ci = F)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 8. Panel model with annual weather
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fe_mod_an <-  feols(log(value + 0.01) ~ tmax_an * ppt_an + year | tree_id,
                    data= paneldat, cluster = ~ plot_id_needle)

summary(fe_mod_an)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 9. Panel model with summer weather
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fe_mod_summer <-  feols(log(value + 0.01) ~ tmaxSummer * pptSummer + year | tree_id,
                        data= paneldat, cluster = ~ plot_id_needle)

summary(fe_mod_summer)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 10. LM model
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lm_mod <- lm(log(value + 0.01) ~ tmax * ppt + year + factor(plot_id_needle),
             data = paneldat)

summary(lm_mod)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 11. LMM with random intercepts
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lmm_mod_int = lmer(log(value + 0.01) ~ tmax *  ppt + year +  (1|tree_id/plot_id_needle), 
                   control = lmerControl(optimizer = "bobyqa",
                                         optCtrl = list(maxfun = 100000)), 
                   data = paneldat)

summary(lmm_mod_int)
tab_model(lmm_mod_int, show.se = T, show.ci = F, digits = 4)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 12. LMM with random slopes and intercepts
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lmm_mod_rslopes = lmer(log(value + 0.01) ~ tmax +   ppt  + year +
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
lmm_mod_rslopes_only <- lmer(log(value + 0.01) ~ tmax + ppt + year +
                               (0 + tmax | tree_id) +
                               (0 + ppt | tree_id) +
                               (0 + year | tree_id), 
                             control = lmerControl(optimizer = "bobyqa",
                                                   optCtrl = list(maxfun = 100000)), 
                             data = paneldat)


tab_model(lmm_mod_rslopes_only, show.se = T, show.ci = F, digits = 4)




rm(rwldat, climdat)


length(ls(envir = .GlobalEnv))





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create spec chart --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Create figure
highlight_n <- 1

## Define label structure
labels <- list("Standard errors" = c("Clustered", "Heteroskedasticity\nrobust"),
               "Climate relationships" = c("Linear", "Quadratic", "Interaction"),
               "Weather window" = c("Water year", "Annual", "Summer", "Lagged"),
               "Model structure" = c("Fixed effects", "Random intercept"))

specs <- specs %>% drop_na()
svg(paste0(wdir, 'robustness.svg'), width = 9, height = 14)
par(oma=c(1,0,1,1))

robustness_fig <- schart(specs, labels, highlight=highlight_n, order = "asis", 
                         # cex=(1.2),fonts=c(2,3),
                         # heights = c(.4,.6),
                         n=c(10), ci = c(.95),
                         # n=c(1, 2, 1, 1, 4, 3), ci = c(.95),
                         ylab = "Average marginal effect of\ntemperature on growth",
                         col.est=c("grey80", "dodgerblue4"),
                         col.dot=c("grey60","grey95","grey95","dodgerblue4"),
                         bg.dot=c("grey60","grey95","grey95","dodgerblue4"),
                         lwd.symbol=1)

# text(x = 1, y = -.008, label = "1",
#      col = "black",   # Color of the text
#      font = 2,      # Bold face
#      cex = 1)     # Size
# 
# text(x = 3, y = -.008, label = "2",
#      col = "black",   # Color of the text
#      font = 2,      # Bold face
#      cex = 1)     # Size
# 
# text(x = 4, y = -.008, label = "3",
#      col = "black",   # Color of the text
#      font = 2,      # Bold face
#      cex = 1)     # Size
# 
# text(x = 6, y = -.008, label = "4",
#      col = "black",   # Color of the text
#      font = 2,      # Bold face
#      cex = 1)     # Size
# 
# text(x = 8, y = -.008, label = "5",
#      col = "black",   # Color of the text
#      font = 2,      # Bold face
#      cex = 1)     # Size
# 
# text(x = 10, y = -.008, label = "6",
#      col = "black",   # Color of the text
#      font = 2,      # Bold face
#      cex = 1)     # Size
# 
# text(x = 11, y = -.008, label = "7",
#      col = "black",   # Color of the text
#      font = 2,      # Bold face
#      cex = 1)     # Size
# 
# text(x = 12, y = -.008, label = "8",
#      col = "black",   # Color of the text
#      font = 2,      # Bold face
#      cex = 1)     # Size
# 
# text(x = 13, y = -.008, label = "9",
#      col = "black",   # Color of the text
#      font = 2,      # Bold face
#      cex = 1)     # Size
# 
# text(x = 15, y = -.008, label = "10",
#      col = "black",   # Color of the text
#      font = 2,      # Bold face
#      cex = 1)     # Size
# 
# text(x = 16, y = -.008, label = "11",
#      col = "black",   # Color of the text
#      font = 2,      # Bold face
#      cex = 1)     # Size
# 
# text(x = 17, y = -.008, label = "12",
#      col = "black",   # Color of the text
#      font = 2,      # Bold face
#      cex = 1)     # Size

dev.off()

