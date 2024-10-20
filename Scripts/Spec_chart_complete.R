## ---------------------------
##
## Script name:
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
## Notes: This code creates the specification chart in Figure X
##  that uses code from: https://github.com/ArielOrtizBobea/spec_chart
##
## ---------------------------

librarian::shelf(here, broom, broom.mixed, here, tidyverse)

# Load specification chart function and model results

source(here("Scripts", "spec_chart_function.R"))
source(here("Scripts", "Models_forSpecChart.R"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Function to pull coefficients
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


extract_and_combine_coefs <- function(models, model_names) {
  # Check if models and model_names have the same length
  if (length(models) != length(model_names)) {
    stop("The number of models and model names must be the same.")
  }

  # Function to extract coefficients and SE
  extract_coef_se <- function(model, model_name) {
    if (inherits(model, "lmerMod")) {
      tidy(model) %>%
        filter(effect == "fixed", term %in% c("tmax", "ppt", "year")) %>%
        select(term, estimate, std.error) %>%
        mutate(model = model_name)
    } else {
      tidy(model) %>%
        filter(term %in% c("tmax", "ppt", "year", "tmaxSummer", "tmax_an")) %>%
        select(term, estimate, std.error) %>%
        mutate(model = model_name)
    }
  }

  # Apply extract_coef_se to each model and combine results
  all_coefs <- map2_dfr(models, model_names, extract_coef_se)

  return(all_coefs)
}

#ls(envir = .GlobalEnv)

models_list <- list(fe_mod_noCE, fe_mod_hetero, fe_mod_noyr, fe_mod, fe_mod_quad, fe_mod_no_int, fe_mod_an, 
                    fe_mod_summer, fe_mod_lag, lmm_mod_int, lmm_mod_rslopes, lmm_mod_rslopes_only)


# Corresponding names for each model
model_names <- c("No clustering", "Heteroskedasticity", "W/out year", "Target model",
                 "Quadratic", "No interaction", "Annual weather", "Summer weather", "Lagged weather",
                 "Intercepts", "Intercepts & slopes", "Slopes only")

# Extract and combine coefficients
all_coefs <- extract_and_combine_coefs(models_list, model_names)


spec_data <- all_coefs %>%
  mutate(term = case_when(
    term == "tmaxSummer" ~ "tmax",
    term == "tmax_an" ~ "tmax",
    TRUE ~ term
  )) %>% 
  pivot_wider(names_from = term,
              values_from = c(estimate, std.error),
              names_glue = "{term}_{.value}") %>%
  mutate(across(where(is.numeric), ~round(., 4))) %>%
  mutate(
    fe_mod_noCE = model == "No clustering",
    fe_mod_hetero = model == "Heteroskedasticity",
    fe_mod_noyr = model == "W/out year",
    fe_mod = model == "Target model",
    fe_mod_quad = model == "Quadratic",
    fe_mod_cube = model == "No interaction",
    fe_mod_an = model == "Annual weather",
    fe_mod_summer = model == "Summer weather",
    fe_mod_lag = model == "Lagged weather",
    lmm_mod_int = model == "Intercepts",
    lmm_mod_rslopes = model == "Intercepts & slopes",
    lmm_mod_rslopes_only = model == "Slopes only"
  ) %>%
  #select(-model)
  select(-c(model, ppt_estimate, ppt_std.error, year_estimate, year_std.error))


data <- as.data.frame(spec_data)

labels <- list("FE models:" = c("No clustering", "Heteroskedasticity", "W/out year", "Target model", "No interaction", "Cubic"),
                 "Weather windows:" = c("Annual weather", "Summer weather", "Lagged weather"),
                 "Mixed models:" = c("Intercepts", "Intercepts & slopes", "Slopes only"))

#===============================================================================
# THE SPEC CHART
#===============================================================================
# Looks better when there is an outer margins
par(oma=c(1,0,1,1))

## use schart function to create spec chart
schart(data, labels, cex=(1.2),fonts=c(2,3), n=c(6,3,3), highlight=4, ci=c(.95,.99),
       col.est=c("grey60","#91376F"), col.dot=c("grey60","grey95","grey95","#91376F"),
       col.est2=c("grey","grey"))



