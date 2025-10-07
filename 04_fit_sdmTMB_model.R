# Load required libraries
library(sdmTMB)
library(sf)
library(tidyverse)
library(DHARMa)

# Read in data
dat_train <- read_rds("data/eeev_env_covs_train_scaled.rds")
mesh_sdmtmb <- read_rds("output/mesh_sdmTMB_750.rds")

# Changing to data.frame and dropping unused factor levels and geometry
dat_train <- dat_train %>% 
  st_drop_geometry() %>% 
  droplevels()

dat_train <- as.data.frame(dat_train)


# Fit primary model with selected predictors
eeev_mod <- sdmTMB(eeev_prop ~ poly(prcp_lag1, 2) + poly(prcp_lag5, 2) + 
                     poly(prcp_lag12, 2) + poly(tmax_lag6, 2) + 
                     poly(tmax_lag12, 2) + poly(tmin_lag1, 2) + 
                     poly(forest, 2) + poly(wetlands, 2) + 
                     (1 | county) + (1 | site_id), 
                   data = dat_train,
                   weights = dat_train$weights,
                   extra_time = 157:180,
                   mesh = mesh_sdmtmb,
                   spatial = "off",
                   spatiotemporal = "ar1",
                   family = binomial(link = "cloglog"),
                   time = "time_step")

# Save model output
write_rds(eeev_mod, "output/eeev_model_primary_ver00.rds")

# Fit alternate model with secondary candidate set
eeev_mod_alt <- sdmTMB(eeev_prop ~ poly(tmin, 2) + poly(tmin_lag6, 2) + 
                         poly(tmin_lag12, 2) + poly(tmax_lag1, 2) + 
                         poly(forest, 2) + poly(wetlands, 2) + 
                         (1 | county) + (1 | site_id), 
                       data = dat_train,
                       weights = dat_train$weights,
                       extra_time = 157:180,
                       mesh = mesh_sdmtmb,
                       spatial = "off",
                       spatiotemporal = "ar1",
                       family = binomial(link = "cloglog"),
                       time = "time_step")

# Save model output
write_rds(eeev_mod_alt, "output/eeev_model_alt_ver00.rds")