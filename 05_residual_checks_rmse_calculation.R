# Clear environment and perform garbage collection
rm(list = ls())
gc()

# Load required libraries
library(tidyverse)
library(sf)
library(sdmTMB)
library(Metrics)

# Read in model objects
eeev_mod <- read_rds("output/eeev_model_primary_ver01.rds")
eeev_mod_alt <- read_rds("output/eeev_model_alt_ver01.rds")

# Sanity checks
sanity(eeev_mod)#all good
sanity(eeev_mod_alt)#all good

# Residual checks
sim_eeev_mod <- simulate(eeev_mod, nsim = 500, type = "mle-mvn")
sim_eeev_mod_alt <- simulate(eeev_mod_alt, nsim = 500, type = "mle-mvn")

resid_eeev_mod <- dharma_residuals(sim_eeev_mod, eeev_mod, return_DHARMa = TRUE)
resid_eeev_mod_alt <- dharma_residuals(sim_eeev_mod_alt, eeev_mod_alt, return_DHARMa = TRUE)

plot(resid_eeev_mod)
plot(resid_eeev_mod_alt)

## Calculating RMSE
# Read in testing data
dat_test <- read_rds("data/eeev_env_covs_test_scaled.rds") %>% 
  st_drop_geometry()

# Predict on testing data
eeev_mod_testing_preds <- eeev_mod %>% 
  predict(newdata = dat_test, 
          se_fit = F, # avoid calculating standard errors for computational effeciency 
          # random effect structures:
          re_form = NA, # mean population level estimates for site/county intercepts
          re_form_iid = NA) %>% # include spatial random effects (SPDE random fields)
  mutate(date = make_date(month = as.character(month), year = as.character(year)))

eeev_mod_alt_testing_preds <- eeev_mod_alt %>% 
  predict(newdata = dat_test, 
          se_fit = F, # avoid calculating standard errors for computational effeciency 
          # random effect structures:
          re_form = NA, # mean population level estimates for site/county intercepts
          re_form_iid = NA) %>% # include spatial random effects (SPDE random fields)
  mutate(date = make_date(month = as.character(month), year = as.character(year)))

# Convert prediction estimates to response scale
eeev_mod_testing_preds <- eeev_mod_testing_preds %>% 
  mutate(predcited_prop = 1 - exp(-exp(est)))

eeev_mod_alt_testing_preds <- eeev_mod_alt_testing_preds %>% 
  mutate(predcited_prop = 1 - exp(-exp(est)))

# Calculate rmse
# Total
rmse_eeev_mod_total <- rmse(eeev_mod_testing_preds$eeev_prop, eeev_mod_testing_preds$predcited_prop)
rmse_eeev_mod_alt_total <- rmse(eeev_mod_alt_testing_preds$eeev_prop, eeev_mod_alt_testing_preds$predcited_prop)

# Yearly
rmse_eeev_mod_yearly <- eeev_mod_testing_preds %>% 
  group_by(year) %>% 
  summarise(rmse = rmse(eeev_prop, predcited_prop))

rmse_eeev_mod_alt_yearly <- eeev_mod_alt_testing_preds %>% 
  group_by(year) %>% 
  summarise(rmse = rmse(eeev_prop, predcited_prop))

# Monthly
rmse_eeev_mod_monthly <- eeev_mod_testing_preds %>% 
  group_by(year, month) %>% 
  summarise(rmse = rmse(eeev_prop, predcited_prop))

rmse_eeev_mod_alt_monthly <- eeev_mod_alt_testing_preds %>% 
  group_by(year, month) %>% 
  summarise(rmse = rmse(eeev_prop, predcited_prop))

# Add "model" columns to identify each set
rmse_yearly <- bind_rows(
  rmse_eeev_mod_yearly  %>% mutate(model = 'Primary'),
  rmse_eeev_mod_alt_yearly %>% mutate(model = 'Alternative')
) %>% select(model, year, rmse)

rmse_monthly <- bind_rows(
  rmse_eeev_mod_monthly %>% mutate(model = 'Primary'),
  rmse_eeev_mod_alt_monthly %>% mutate(model = 'Alternative')
) %>% select(model, year, month, rmse)

# Add total summary rows
rmse_total <- tibble(
  model = c('Primary', 'Alternative'),
  year = NA,
  month = NA,
  rmse = c(rmse_eeev_mod_total, rmse_eeev_mod_alt_total)
)

# Combine all into one tidy table
rmse_table <- bind_rows(
  rmse_total,
  rmse_yearly %>% mutate(month = NA),
  rmse_monthly
) %>%
  arrange(model, year, month)

print(rmse_table, n = Inf)

# Extract fixed effects parameter estimates
eeev_mod_fixed_effects <- tidy(eeev_mod, effects = "fixed", conf.int = TRUE)
eeev_mod_alt_fixed_effects <- tidy(eeev_mod_alt, effects = "fixed", conf.int = TRUE)

# Creating a tidy table of fixed effects
var_lookup <- c(
  "prcp_lag1" = "Precipitation Lag 1 Month",
  "prcp_lag5" = "Precipitation Lag 5 Months",
  "prcp_lag12" = "Precipitation Lag 12 Months",
  "tmax_lag6" = "Tmax Lag 6 Months",
  "tmax_lag12" = "Tmax Lag 12 Months",
  "tmin_lag1" = "Tmin Lag 1 Month",
  "forest" = "Forest Cover",
  "wetlands" = "Wetland Cover"
)

# Clean and transform
tidy_table_eeev_mod <- eeev_mod_fixed_effects %>%
  mutate(
    variable = case_when(
      term == "(Intercept)" ~ "Intercept",
      TRUE ~ str_extract(term, "(?<=poly\\().+?(?=,)")  
    ),
    degree = case_when(
      term == "(Intercept)" ~ NA_character_,
      str_detect(term, "1$") ~ "Linear",
      str_detect(term, "2$") ~ "Quadratic"
    ),
    variable_clean = if_else(
      variable == "Intercept", "Intercept",
      var_lookup[variable]
    )
  ) %>%
  select(
    Variable = variable_clean,
    `Linear / Quad` = degree,
    Estimate = estimate,
    `Std. Error` = std.error,
    `95% CI (Lower)` = conf.low,
    `95% CI (Upper)` = conf.high
  )

# Alternative model
var_lookup_alt <- c(
  "tmin" = "Tmin (Current)",
  "tmin_lag6" = "Tmin Lag 6 Months",
  "tmin_lag12" = "Tmin Lag 12 Months",
  "tmax_lag1" = "Tmax Lag 1 Month",
  "forest" = "Forest Cover",
  "wetlands" = "Wetland Cover"
)

# Clean and tidy table
tidy_table_eeev_mod_alt <- eeev_mod_alt_fixed_effects %>%
  mutate(
    variable = case_when(
      term == "(Intercept)" ~ "Intercept",
      TRUE ~ str_extract(term, "(?<=poly\\().+?(?=,)")  
    ),
    degree = case_when(
      term == "(Intercept)" ~ NA_character_,
      str_detect(term, "1$") ~ "Linear",
      str_detect(term, "2$") ~ "Quadratic"
    ),
    variable_clean = if_else(
      variable == "Intercept", "Intercept",
      var_lookup_alt[variable]
    )
  ) %>%
  select(
    Variable = variable_clean,
    `Linear / Quad` = degree,
    Estimate = estimate,
    `Std. Error` = std.error,
    `95% CI (Lower)` = conf.low,
    `95% CI (Upper)` = conf.high
  )


# Save data
write_csv(rmse_table, "output/rmse_table_ver01.csv")
write_csv(tidy_table_eeev_mod, "output/fixed_effects_primary_model.csv")
write_csv(tidy_table_eeev_mod_alt, "output/fixed_effects_alt_model.csv")
write_rds(resid_eeev_mod, "output/eeev_model_primary_residuals_ver01.rds")
write_rds(resid_eeev_mod_alt, "output/eeeev_model_alt_residuals_ver01.rds")
write_rds(eeev_mod_testing_preds, "output/eeev_model_primary_testing_preds_ver01.rds")
write_rds(eeev_mod_alt_testing_preds, "output/eeev_model_alt_testing_preds_ver01.rds")
save.image("06_residual_checks_rmse_calculation.RData")