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
# rwldat <- read_csv(paste0(wdir, "paneldat_RWL.csv"))
# climdat <- read_csv(paste0(wdir, "climatewindows.csv"))
# latlong <- read_csv(paste0(wdir, "plotlatlong.csv")) %>% rename(plot_id_needle = plot)
# paneldat <- rwldat %>% 
#   left_join(climdat) %>% 
#   left_join(latlong, by = "plot_id_needle")

wdir <- 'G:/My Drive/collaborations/ecology/CausalInferenceClimateAttribution/'
paneldat <- read_csv(paste0(wdir, "cleaned_RWL_climdat.csv"))





winsorize_value <- paneldat %>% filter(value>0) %>% pull(value) %>% min()
paneldat <- paneldat %>% 
  mutate(value_w = ifelse(value==0, winsorize_value, value))


# Function to extract coefficients and SE
extract_coef_se <- function(model) {
  slopes_df <- avg_slopes(model) %>% 
    as_tibble() %>% 
    print()
  slopes_df %>% 
    filter(term %in% c("tmax", "tmaxSummer", "tmax_an")) %>% 
    select(term, estimate, std.error)
}


specs <- data.frame(coef=NaN, 
                    se=NaN,
                    clim_lin = NaN,
                    clim_quad = NaN,
                    clim_int = NaN,
                    clim_lag = FALSE,
                    ww_wy = NaN,
                    ww_an = NaN,
                    ww_sum = NaN,
                    struc_yrtrend = NaN,
                    struc_plotfe = NaN,
                    struc_treefe = NaN,
                    struc_ri = NaN,
                    struc_rs = NaN,
                    drop_max = NaN,
                    drop_min = NaN,
                    stder_clust = NaN,
                    stder_het = NaN) 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1. The panel model 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fe_mod <-  feols(log(value + 0.01) ~ tmax * ppt + year | tree_id,
                 data= paneldat, cluster = ~ plot_id_needle)

# fe_mod <-  feols(log(value_w) ~ tmax * ppt + year | tree_id,
#                  data= paneldat, vcov = "conley")

# summary(fe_mod)
# avg_slopes(fe_mod)

fe_coef <- extract_coef_se(fe_mod)
new_row <-  data_frame(coef=fe_coef$estimate, 
                       se=fe_coef$std.error,
                       clim_lin = TRUE,
                       clim_quad = FALSE,
                       clim_int = TRUE,
                       clim_lag = FALSE,
                       ww_wy = TRUE,
                       ww_an = FALSE,
                       ww_sum = FALSE,
                       struc_yrtrend = TRUE,
                       struc_treefe = TRUE,
                       struc_plotfe = FALSE,
                       struc_ri = FALSE,
                       struc_rs = FALSE,
                       drop_max = FALSE,
                       drop_min = FALSE,
                       stder_clust = TRUE,
                       stder_het = FALSE)
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

fe_mod_hetero <-  feols(log(value + 0.01) ~ tmax * ppt + year | tree_id,
                        data= paneldat, vcov = "hetero")

# summary(fe_mod_hetero)
# avg_slopes(fe_mod_hetero)

fe_hetero_coef <- extract_coef_se(fe_mod_hetero)
new_row <-  data_frame(coef=fe_hetero_coef$estimate, 
                       se=fe_hetero_coef$std.error,
                       clim_lin = TRUE,
                       clim_quad = FALSE,
                       clim_int = TRUE,
                       clim_lag = FALSE,
                       ww_wy = TRUE,
                       ww_an = FALSE,
                       ww_sum = FALSE,
                       struc_yrtrend = TRUE,
                       struc_treefe = TRUE,
                       struc_plotfe = FALSE,
                       struc_ri = FALSE,
                       struc_rs = FALSE,
                       drop_max = FALSE,
                       drop_min = FALSE,
                       stder_clust = FALSE,
                       stder_het = TRUE)
specs <- rbind(specs, new_row)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5. Panel model with quadratic terms
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fe_mod_quad <-  feols(log(value + 0.01) ~ tmax + ppt + I(tmax^2) + I(ppt^2) + year | tree_id,
                      data= paneldat, cluster = ~ plot_id_needle)
summary(fe_mod_quad)
fe_quad_coef <- extract_coef_se(fe_mod_quad)
new_row <-  data_frame(coef=fe_quad_coef$estimate, 
                       se=fe_quad_coef$std.error,
                       clim_lin = FALSE,
                       clim_quad = TRUE,
                       clim_int = FALSE,
                       clim_lag = FALSE,
                       ww_wy = TRUE,
                       ww_an = FALSE,
                       ww_sum = FALSE,
                       struc_yrtrend = TRUE,
                       struc_treefe = TRUE,
                       struc_plotfe = FALSE,
                       struc_ri = FALSE,
                       struc_rs = FALSE,
                       drop_max = FALSE,
                       drop_min = FALSE,
                       stder_clust = TRUE,
                       stder_het = FALSE)
specs <- rbind(specs, new_row)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 6. Panel model without interaction effects 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fe_mod_no_int<-  feols(log(value + 0.01) ~ tmax + ppt + year | tree_id,
                     data= paneldat, cluster = ~ plot_id_needle)

summary(fe_mod_no_int)
fe_no_int_coef <- extract_coef_se(fe_mod_no_int)
new_row <-  data_frame(coef=fe_no_int_coef$estimate, 
                       se=fe_no_int_coef$std.error,
                       clim_lin = TRUE,
                       clim_quad = FALSE,
                       clim_int = FALSE,
                       clim_lag = FALSE,
                       ww_wy = TRUE,
                       ww_an = FALSE,
                       ww_sum = FALSE,
                       struc_yrtrend = TRUE,
                       struc_treefe = TRUE,
                       struc_plotfe = FALSE,
                       struc_ri = FALSE,
                       struc_rs = FALSE,
                       drop_max = FALSE,
                       drop_min = FALSE,
                       stder_clust = TRUE,
                       stder_het = FALSE)
specs <- rbind(specs, new_row)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 7. Panel model with lagged effects
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fe_mod_lag <-  feols(log(value + 0.01) ~ tmax*ppt + laggedtmax + laggedprecip + year | tree_id,
                     data= paneldat, cluster = ~ plot_id_needle)

summary(fe_mod_lag)
fe_lag_coef <- extract_coef_se(fe_mod_lag)
new_row <-  data_frame(coef=fe_lag_coef$estimate, 
                       se=fe_lag_coef$std.error,
                       clim_lin = TRUE,
                       clim_quad = FALSE,
                       clim_int = TRUE,
                       clim_lag = TRUE,
                       ww_wy = TRUE,
                       ww_an = FALSE,
                       ww_sum = FALSE,
                       struc_yrtrend = TRUE,
                       struc_treefe = TRUE,
                       struc_plotfe = FALSE,
                       struc_ri = FALSE,
                       struc_rs = FALSE,
                       drop_max = FALSE,
                       drop_min = FALSE,
                       stder_clust = TRUE,
                       stder_het = FALSE)
specs <- rbind(specs, new_row)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 8. Panel model with annual weather
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fe_mod_an <-  feols(log(value + 0.01) ~ tmax_an * ppt_an + year | tree_id,
                    data= paneldat, cluster = ~ plot_id_needle)
summary(fe_mod_an)
fe_an_coef <- extract_coef_se(fe_mod_an)
new_row <-  data_frame(coef=fe_an_coef$estimate, 
                       se=fe_an_coef$std.error,
                       clim_lin = TRUE,
                       clim_quad = FALSE,
                       clim_int = TRUE,
                       clim_lag = FALSE,
                       ww_wy = FALSE,
                       ww_an = TRUE,
                       ww_sum = FALSE,
                       struc_yrtrend = TRUE,
                       struc_treefe = TRUE,
                       struc_plotfe = FALSE,
                       struc_ri = FALSE,
                       struc_rs = FALSE,
                       drop_max = FALSE,
                       drop_min = FALSE,
                       stder_clust = TRUE,
                       stder_het = FALSE)
specs <- rbind(specs, new_row)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 9. Panel model with summer weather
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fe_mod_summer <-  feols(log(value + 0.01) ~ tmaxSummer * pptSummer + year | tree_id,
                        data= paneldat, cluster = ~ plot_id_needle)

summary(fe_mod_summer)
fe_summer_coef <- extract_coef_se(fe_mod_summer)
new_row <-  data_frame(coef=fe_summer_coef$estimate, 
                       se=fe_summer_coef$std.error,
                       clim_lin = TRUE,
                       clim_quad = FALSE,
                       clim_int = TRUE,
                       clim_lag = FALSE,
                       ww_wy = FALSE,
                       ww_an = FALSE,
                       ww_sum = TRUE,
                       struc_yrtrend = TRUE,
                       struc_treefe = TRUE,
                       struc_plotfe = FALSE,
                       struc_ri = FALSE,
                       struc_rs = FALSE,
                       drop_max = FALSE,
                       drop_min = FALSE,
                       stder_clust = TRUE,
                       stder_het = FALSE)
specs <- rbind(specs, new_row)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 10. LM model
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# lm_mod <- lm(log(value + 0.01) ~ tmax * ppt + year + factor(plot_id_needle),
#              data = paneldat)
# 
# summary(lm_mod)

plotfe_mod <- feols(log(value + 0.01) ~ tmax * ppt + year | plot_id_needle,
                     data = paneldat, cluster = ~ plot_id_needle)
plotfe_coef <- extract_coef_se(plotfe_mod)
new_row <-  data_frame(coef=plotfe_coef$estimate, 
                       se=plotfe_coef$std.error,
                       clim_lin = TRUE,
                       clim_quad = FALSE,
                       clim_int = TRUE,
                       clim_lag = FALSE,
                       ww_wy = TRUE,
                       ww_an = FALSE,
                       ww_sum = FALSE,
                       struc_yrtrend = TRUE,
                       struc_treefe = FALSE,
                       struc_plotfe = TRUE,
                       struc_ri = FALSE,
                       struc_rs = FALSE,
                       drop_max = FALSE,
                       drop_min = FALSE,
                       stder_clust = TRUE,
                       stder_het = FALSE)
specs <- rbind(specs, new_row)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 11. LMM with random intercepts
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rei_mod = lmer(log(value + 0.01) ~ tmax *  ppt + year +  (1|tree_id/plot_id_needle), 
                   control = lmerControl(optimizer = "bobyqa",
                                         optCtrl = list(maxfun = 100000)), 
                   data = paneldat)

# summary(lmm_mod_int)
# tab_model(lmm_mod_int, show.se = T, show.ci = F, digits = 4)

rei_coef <- extract_coef_se(rei_mod)
new_row <-  data_frame(coef=rei_coef$estimate, 
                       se=rei_coef$std.error,
                       clim_lin = TRUE,
                       clim_quad = FALSE,
                       clim_int = TRUE,
                       clim_lag = FALSE,
                       ww_wy = TRUE,
                       ww_an = FALSE,
                       ww_sum = FALSE,
                       struc_yrtrend = TRUE,
                       struc_treefe = FALSE,
                       struc_plotfe = FALSE,
                       struc_ri = TRUE, # NOTE: Not sure if i fully understand this lmer model - may need to improve labels here
                       struc_rs = FALSE,
                       drop_max = FALSE,
                       drop_min = FALSE,
                       stder_clust = FALSE, # NOTE: Need to do extra work to get these lmer models with clustered standard errors
                       stder_het = FALSE)
specs <- rbind(specs, new_row)

# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # 12. LMM with random slopes and intercepts
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# reis_mod = lmer(log(value + 0.01) ~ tmax +   ppt  + year +
#                          (1 + tmax | plot_id_needle)+
#                          (1 + ppt | plot_id_needle)+
#                          (1 + year | plot_id_needle), 
#                        control = lmerControl(optimizer = "bobyqa",
#                                              optCtrl = list(maxfun = 100000)), 
#                        data = paneldat)
# 
# 
# # summary(lmm_mod_rslopes)
# # tab_model(lmm_mod_rslopes, show.se = T, show.ci = F, digits = 4)
# 
# 
# reis_coef <- extract_coef_se(reis_mod)
# new_row <-  data_frame(coef=reis_coef$estimate, 
#                        se=reis_coef$std.error,
#                        clim_lin = TRUE,
#                        clim_quad = FALSE,
#                        clim_int = FALSE, # Note: Why did you turn off the interaction in this model?
#                        ww_wy = TRUE,
#                        ww_an = FALSE,
#                        ww_sum = FALSE,
#                        ww_lag = FALSE,
#                        struc_yrtrend = TRUE,
#                        struc_treefe = FALSE,
#                        struc_plotfe = FALSE,
#                        struc_ri = TRUE, # NOTE: Not sure if i fully understand this lmer model - may need to improve labels here
#                        struc_rs = TRUE,
#                        stder_clust = FALSE, # NOTE: Need to do extra work to get these lmer models with clustered standard errors
#                        stder_het = FALSE)
# specs <- rbind(specs, new_row)

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
res_mod <- lmer(log(value + 0.01) ~ tmax * ppt + year +
                               (0 + tmax | tree_id) +
                               (0 + ppt | tree_id) +
                               (0 + year | tree_id), 
                             control = lmerControl(optimizer = "bobyqa",
                                                   optCtrl = list(maxfun = 100000)), 
                             data = paneldat)

res_coef <- extract_coef_se(res_mod)
new_row <-  data_frame(coef=res_coef$estimate, 
                       se=res_coef$std.error,
                       clim_lin = TRUE,
                       clim_quad = FALSE,
                       clim_int = TRUE,
                       clim_lag = FALSE,
                       ww_wy = TRUE,
                       ww_an = FALSE,
                       ww_sum = FALSE,
                       struc_yrtrend = TRUE,
                       struc_treefe = FALSE,
                       struc_plotfe = FALSE,
                       struc_ri = FALSE, # NOTE: Not sure if i fully understand this lmer model - may need to improve labels here
                       struc_rs = TRUE,
                       drop_max = FALSE,
                       drop_min = FALSE,
                       stder_clust = FALSE, # NOTE: Need to do extra work to get these lmer models with clustered standard errors
                       stder_het = FALSE)
specs <- rbind(specs, new_row)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4. Panel model without year effects
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fe_mod_noyr <-  feols(log(value + 0.01) ~ tmax * ppt | tree_id,
                      data= paneldat, cluster = ~ plot_id_needle)

summary(fe_mod_noyr)
fe_noyr_coef <- extract_coef_se(fe_mod_noyr)
new_row <-  data_frame(coef=fe_noyr_coef$estimate, 
                       se=fe_noyr_coef$std.error,
                       clim_lin = TRUE,
                       clim_quad = FALSE,
                       clim_int = TRUE,
                       clim_lag = FALSE,
                       ww_wy = TRUE,
                       ww_an = FALSE,
                       ww_sum = FALSE,
                       struc_yrtrend = FALSE,
                       struc_treefe = TRUE,
                       struc_plotfe = FALSE,
                       struc_ri = FALSE,
                       struc_rs = FALSE,
                       drop_max = FALSE,
                       drop_min = FALSE,
                       stder_clust = TRUE,
                       stder_het = FALSE)
specs <- rbind(specs, new_row)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Jackknife --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plots <- paneldat %>% pull(plot_id_needle) %>% unique()
drop_coefs <- tibble()
for (p in plots) {
  drop_dat <- paneldat %>% filter(plot_id_needle != p)
  
  drop_mod <-  feols(log(value + 0.01) ~ tmax * ppt + year | tree_id,
                     data= drop_dat, cluster = ~ plot_id_needle)
  drop_coef <- extract_coef_se(drop_mod)
  drop_coefs <- rbind(drop_coefs, drop_coef)
}

min_mod <- drop_coefs %>%  slice(which.min(estimate))
max_mod <- drop_coefs %>%  slice(which.max(estimate))

new_row <-  data_frame(coef=max_mod$estimate, 
                       se=max_mod$std.error,
                       clim_lin = TRUE,
                       clim_quad = FALSE,
                       clim_int = TRUE,
                       clim_lag = FALSE,
                       ww_wy = TRUE,
                       ww_an = FALSE,
                       ww_sum = FALSE,
                       struc_yrtrend = TRUE,
                       struc_treefe = TRUE,
                       struc_plotfe = FALSE,
                       struc_ri = FALSE,
                       struc_rs = FALSE,
                       drop_max = TRUE,
                       drop_min = FALSE,
                       stder_clust = TRUE,
                       stder_het = FALSE)
specs <- rbind(specs, new_row)


new_row <-  data_frame(coef=min_mod$estimate, 
                       se=min_mod$std.error,
                       clim_lin = TRUE,
                       clim_quad = FALSE,
                       clim_int = TRUE,
                       clim_lag = FALSE,
                       ww_wy = TRUE,
                       ww_an = FALSE,
                       ww_sum = FALSE,
                       struc_yrtrend = TRUE,
                       struc_treefe = TRUE,
                       struc_plotfe = FALSE,
                       struc_ri = FALSE,
                       struc_rs = FALSE,
                       drop_max = FALSE,
                       drop_min = TRUE,
                       stder_clust = TRUE,
                       stder_het = FALSE)
specs <- rbind(specs, new_row)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create spec chart --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Create figure
highlight_n <- 1

## Define label structure
labels <- list("Climate relationships" = c("Linear", "Quadratic", "Interaction", "Lagged"),
               "Weather window" = c("Water year", "Annual", "Summer"),
               "Model structure" = c("Year trend", "Tree FEs", "Plot FEs", "Nested random intercepts", "Nested random slopes"),
               "Jackknife by plot" = c("Minimum effect", "Maximum effect"),
               "Standard errors" = c("Clustered at plot", "Heteroskedasticity robust"))

specs <- specs %>% drop_na()
png(paste0(wdir, 'robustness.png'), width = 12, height = 10, units = "in", res = 350)
par(oma=c(1,0,1,1))

robustness_fig <- schart(specs, labels, highlight=highlight_n, order = "asis", 
                         # cex=(1.2),fonts=c(2,3),
                         # heights = c(.4,.6),
                         # n=c(13), ci = c(.95),
                         n=c(2, 3,  2, 4, 2), ci = c(.95),
                         ylab = "Average marginal effect of\ntemperature on growth",
                         col.est=c("grey80", "dodgerblue4"),
                         col.dot=c("grey60","grey95","grey95","dodgerblue4"),
                         bg.dot=c("grey60","grey95","grey95","dodgerblue4"),
                         lwd.symbol=1)

text(x = 1, y = .025, label = "A",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 2, y = .025, label = "B",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 4, y = .025, label = "C",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 5, y = .025, label = "D",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 6, y = .025, label = "E",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 8, y = .025, label = "F",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 9, y = .025, label = "G",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 11, y = .025, label = "H",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 12, y = .025, label = "I",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 13, y = .025, label = "J",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 14, y = .025, label = "K",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 16, y = .025, label = "L",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

text(x = 17, y = .025, label = "M",
     col = "black",   # Color of the text
     font = 2,      # Bold face
     cex = 1)     # Size

dev.off()




