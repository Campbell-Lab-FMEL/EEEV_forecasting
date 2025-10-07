# Load required libraries
library(sdmTMB)
library(tidyverse)
library(sf)
library(qpcR)
library(fmesher)

# Read in dat
dat_train <- read_rds("data/eeev_env_covs_train_scaled.rds")
dat_test <- read_rds("data/eeev_env_covs_test_scaled.rds")

# Read in model objects and mesh
eeev_mod <- read_rds("output/eeev_model_primary_ver01.rds")
eeev_mod_alt <- read_rds("output/eeev_model_alt_ver01.rds")
mesh_sdmtmb <- read_rds("output/mesh_sdmTMB_750.rds")

# Combine data for predictions
dat <- bind_rows(dat_train, dat_test) %>%
  arrange(year, month, site_id) %>%
  mutate(fold = ifelse(year %in% 2018:2019, 1, 2))

# Dropping unused factor levels and geometry from training and testing data
# Training
dat_train <- dat_train %>% 
  st_drop_geometry() %>% 
  droplevels()

dat_train <- as.data.frame(dat_train)

# Testing
dat_test <- dat_test %>%
  st_drop_geometry() %>%
  droplevels()

dat_test <- as.data.frame(dat_test)

## Fitting model variants for model selection
## EEEV Model variants (based on eeev_mod)
# eeev_mod is already fitted - this is the full spatiotemporal model with fixed + random effects

# Random effects + spatiotemporal
eeev_mod_ranef <- sdmTMB(eeev_prop ~ (1 | county) + (1 | site_id),
                         data = dat_train,
                         weights = dat_train$weights,
                         extra_time = 157:180,
                         mesh = mesh_sdmtmb,
                         spatial = "off",
                         spatiotemporal = "ar1",
                         family = binomial(link = "cloglog"),
                         time = "time_step")

# Fixed effects + spatiotemporal
eeev_mod_fixed <- sdmTMB(eeev_prop ~ poly(prcp_lag1, 2) + poly(prcp_lag5, 2) +
                           poly(prcp_lag12, 2) + poly(tmax_lag6, 2) +
                           poly(tmax_lag12, 2) + poly(tmin_lag1, 2) +
                           poly(forest, 2) + poly(wetlands, 2),
                         data = dat_train,
                         weights = dat_train$weights,
                         extra_time = 157:180,
                         mesh = mesh_sdmtmb,
                         spatial = "off",
                         spatiotemporal = "ar1",
                         family = binomial(link = "cloglog"),
                         time = "time_step")

# Intercept + spatiotemporal
eeev_mod_intercept <- sdmTMB(eeev_prop ~ 1,
                             data = dat_train,
                             weights = dat_train$weights,
                             extra_time = 157:180,
                             mesh = mesh_sdmtmb,
                             spatial = "off",
                             spatiotemporal = "ar1",
                             family = binomial(link = "cloglog"),
                             time = "time_step")

# Non-spatial models (spatiotemporal = "off")
# Fixed + random effects, non-spatial
glmmtmb_eeev <- sdmTMB(eeev_prop ~ poly(prcp_lag1, 2) + poly(prcp_lag5, 2) +
                         poly(prcp_lag12, 2) + poly(tmax_lag6, 2) +
                         poly(tmax_lag12, 2) + poly(tmin_lag1, 2) +
                         poly(forest, 2) + poly(wetlands, 2) +
                         (1 | county) + (1 | site_id),
                       data = dat_train,
                       weights = dat_train$weights,
                       extra_time = 157:180,
                       mesh = mesh_sdmtmb,
                       spatial = "off",
                       spatiotemporal = "off",
                       family = binomial(link = "cloglog"),
                       time = "time_step")

# Random effects non-spatial
glmmtmb_eeev_ranef <- sdmTMB(eeev_prop ~ (1 | county) + (1 | site_id),
                             data = dat_train,
                             weights = dat_train$weights,
                             extra_time = 157:180,
                             mesh = mesh_sdmtmb,
                             spatial = "off",
                             spatiotemporal = "off",
                             family = binomial(link = "cloglog"),
                             time = "time_step")

# Intercept non-spatial
glmmtmb_eeev_intercept <- sdmTMB(eeev_prop ~ 1,
                                 data = dat_train,
                                 weights = dat_train$weights,
                                 extra_time = 157:180,
                                 mesh = mesh_sdmtmb,
                                 spatial = "off",
                                 spatiotemporal = "off",
                                 family = binomial(link = "cloglog"),
                                 time = "time_step")

## EEEV Model Alt variants (based on eeev_mod_alt)
# eeev_mod_alt is already fitted - this is the full spatiotemporal model with fixed + random effects

# Random effects + spatiotemporal
eeev_mod_alt_ranef <- sdmTMB(eeev_prop ~ (1 | county) + (1 | site_id),
                             data = dat_train,
                             weights = dat_train$weights,
                             extra_time = 157:180,
                             mesh = mesh_sdmtmb,
                             spatial = "off",
                             spatiotemporal = "ar1",
                             family = binomial(link = "cloglog"),
                             time = "time_step")

# Fixed effects + spatiotemporal
eeev_mod_alt_fixed <- sdmTMB(eeev_prop ~ poly(tmin, 2) + poly(tmin_lag6, 2) +
                               poly(tmin_lag12, 2) + poly(tmax_lag1, 2) +
                               poly(forest, 2) + poly(wetlands, 2),
                             data = dat_train,
                             weights = dat_train$weights,
                             extra_time = 157:180,
                             mesh = mesh_sdmtmb,
                             spatial = "off",
                             spatiotemporal = "ar1",
                             family = binomial(link = "cloglog"),
                             time = "time_step")

# Intercept + spatiotemporal
eeev_mod_alt_intercept <- sdmTMB(eeev_prop ~ 1,
                                 data = dat_train,
                                 weights = dat_train$weights,
                                 extra_time = 157:180,
                                 mesh = mesh_sdmtmb,
                                 spatial = "off",
                                 spatiotemporal = "ar1",
                                 family = binomial(link = "cloglog"),
                                 time = "time_step")

# Non-spatial models (spatiotemporal = "off")
# Fixed + random effects, non-spatial
glmmtmb_eeev_alt <- sdmTMB(eeev_prop ~ poly(tmin, 2) + poly(tmin_lag6, 2) +
                             poly(tmin_lag12, 2) + poly(tmax_lag1, 2) +
                             poly(forest, 2) + poly(wetlands, 2) +
                             (1 | county) + (1 | site_id),
                           data = dat_train,
                           weights = dat_train$weights,
                           extra_time = 157:180,
                           mesh = mesh_sdmtmb,
                           spatial = "off",
                           spatiotemporal = "off",
                           family = binomial(link = "cloglog"),
                           time = "time_step")

# Random effects non-spatial
glmmtmb_eeev_alt_ranef <- sdmTMB(eeev_prop ~ (1 | county) + (1 | site_id),
                                 data = dat_train,
                                 weights = dat_train$weights,
                                 extra_time = 157:180,
                                 mesh = mesh_sdmtmb,
                                 spatial = "off",
                                 spatiotemporal = "off",
                                 family = binomial(link = "cloglog"),
                                 time = "time_step")

# Intercept non-spatial
glmmtmb_eeev_alt_intercept <- sdmTMB(eeev_prop ~ 1,
                                     data = dat_train,
                                     weights = dat_train$weights,
                                     extra_time = 157:180,
                                     mesh = mesh_sdmtmb,
                                     spatial = "off",
                                     spatiotemporal = "off",
                                     family = binomial(link = "cloglog"),
                                     time = "time_step")

## Create model lists
eeev_mods <- list(
  eeev_mod,  
  eeev_mod_ranef,
  eeev_mod_fixed,
  eeev_mod_intercept,
  glmmtmb_eeev,
  glmmtmb_eeev_ranef,
  glmmtmb_eeev_intercept)

eeev_mods_alt <- list(
  eeev_mod_alt,   
  eeev_mod_alt_ranef,
  eeev_mod_alt_fixed,
  eeev_mod_alt_intercept,
  glmmtmb_eeev_alt,
  glmmtmb_eeev_alt_ranef,
  glmmtmb_eeev_alt_intercept)

# Save models
write_rds(eeev_mods, "output/eeev_mods_selection.rds")
write_rds(eeev_mods_alt, "output/eeev_mods_alt_selection.rds")

# Read in model lists
eeev_mods <- read_rds("output/eeev_mods_selection.rds")
eeev_mods_alt <- read_rds("output/eeev_mods_alt_selection.rds")

## Model comparison metrics
# Conditional percent deviance explained using log-likelihood
eeev_mods_devexp <- data.frame(
  model = c(
    "GMRF, fixed effects, random effects",
    "GMRF, random effects", 
    "GMRF, fixed effects",
    "GMRF, intercept",
    "Non-spatial, fixed effects, random effects",
    "Non-spatial, random effects",
    "Non-spatial, intercept"),
  devexp = lapply(eeev_mods, FUN = function(x){
    # Calculate deviance explained using log-likelihood
    # Deviance = -2 * log-likelihood
    null_deviance <- -2 * logLik(eeev_mods[[7]])  # intercept-only model
    model_deviance <- -2 * logLik(x)
    1 - (model_deviance / null_deviance)
  }) %>%
    unlist())

eeev_mods_alt_devexp <- data.frame(
  model = c(
    "GMRF, fixed effects, random effects",
    "GMRF, random effects",
    "GMRF, fixed effects", 
    "GMRF, intercept",
    "Non-spatial, fixed effects, random effects",
    "Non-spatial, random effects",
    "Non-spatial, intercept"),
  devexp = lapply(eeev_mods_alt, FUN = function(x){
    # Calculate deviance explained using log-likelihood
    null_deviance <- -2 * logLik(eeev_mods_alt[[7]])  # intercept-only model
    model_deviance <- -2 * logLik(x)
    1 - (model_deviance / null_deviance)
  }) %>%
    unlist())


# Conditional AIC
aic_eeev <- data.frame(
  model = c(
    "GMRF, fixed effects, random effects",
    "GMRF, random effects",
    "GMRF, fixed effects",
    "GMRF, intercept",
    "Non-spatial, fixed effects, random effects",
    "Non-spatial, random effects",
    "Non-spatial, intercept"),
  aic = lapply(eeev_mods, FUN = function(x){
    AIC(x)}) %>%
    unlist()) %>%
  mutate(delta_aic = aic - min(aic),
         aic_weights = round(qpcR::akaike.weights(aic)$weights, 3))

aic_eeev_alt <- data.frame(
  model = c(
    "GMRF, fixed effects, random effects",
    "GMRF, random effects",
    "GMRF, fixed effects",
    "GMRF, intercept",
    "Non-spatial, fixed effects, random effects",
    "Non-spatial, random effects", 
    "Non-spatial, intercept"),
  aic = lapply(eeev_mods_alt, FUN = function(x){
    AIC(x)}) %>%
    unlist()) %>%
  mutate(delta_aic = aic - min(aic),
         aic_weights = round(qpcR::akaike.weights(aic)$weights, 3))

print(aic_eeev)
print(aic_eeev_alt)

# Primary model comparison table
primary_model_comparison <- data.frame(
  Model = c(
    "GMRF, fixed effects, random effects",
    "GMRF, random effects",
    "GMRF, fixed effects",
    "GMRF, intercept",
    "Non-spatial, fixed effects, random effects",
    "Non-spatial, random effects",
    "Non-spatial, intercept"
  ),
  `Deviance Explained` = round(eeev_mods_devexp$devexp, 3),
  AIC = round(aic_eeev$aic, 1),
  `ΔAIC` = round(aic_eeev$delta_aic, 1),
  `AIC Weight` = aic_eeev$aic_weights,
  check.names = FALSE  
)

# Alternative model comparison table
alt_model_comparison <- data.frame(
  Model = c(
    "GMRF, fixed effects, random effects",
    "GMRF, random effects", 
    "GMRF, fixed effects",
    "GMRF, intercept",
    "Non-spatial, fixed effects, random effects",
    "Non-spatial, random effects",
    "Non-spatial, intercept"
  ),
  `Deviance Explained` = round(eeev_mods_alt_devexp$devexp, 3),
  AIC = round(aic_eeev_alt$aic, 1),
  `ΔAIC` = round(aic_eeev_alt$delta_aic, 1),
  `AIC Weight` = aic_eeev_alt$aic_weights,
  check.names = FALSE
)


# Save tables
write_csv(primary_model_comparison, "output/primary_model_comparison.csv")
write_csv(alt_model_comparison, "output/alt_model_comparison.csv")