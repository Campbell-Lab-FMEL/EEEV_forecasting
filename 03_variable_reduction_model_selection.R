# Clear environment and perform garbage collection
rm(list = ls())
gc()

# Load required libraries
library(tidyverse)
library(terra)
library(sf)
library(fmesher)
library(sdmTMB)
library(blockCV)
library(mgcv)
library(olsrr)
library(car)
library(broom)
library(performance)
library(glmmTMB)
library(Metrics)
library(DHARMa)


# Read in data
dat <- read_rds("data/eeev_env_covs.rds")

# Read and process spatial domain
domain <- read_rds("data/fl_poly_trim.rds")
domain <- st_make_valid(domain)
plot(domain)
#domain <- st_geometry(domain)

# Extracting sites that fall within domain, the domain is a non-convex hull around sites that have at least
# 1 positive case during the sampling period. We do this to reduce the noise in the data.
# Create sites object
sites_sf <- dat %>%
  group_by(site_id) %>%
  summarise(geometry = st_geometry(first(geometry)))

# Transform crs
sites_sf <- st_transform(sites_sf, st_crs(domain))

plot(sites_sf)
plot(domain)

# Identify which sites fall within the domain
sites_in_domain <- sites_sf[st_within(sites_sf, st_geometry(domain), sparse = FALSE), ]

# Extract site_ids that fall within the domain
site_ids <- unique(sites_in_domain$site_id)

# Filter data
dat <- dat %>% 
  filter(site_id %in% site_ids)

# Filter to remove years earlier than 2005, before which EEEV couldn't be differentiated
# from Highlands J virus
dat <- dat %>% 
  filter(year >= 2005)

# Check
range(dat$year)#2005 - 2019

# Preparing data for model
min_year <- min(dat$year)

dat_monthly <- dat %>%
  st_transform(32617) %>%
  mutate(
    date = make_date(month = as.character(month), year = as.character(year)),
    disease = "eeev",
    time_step = ((dat$year - min_year) * 12 + dat$month), 
    weights = testing * 6
  ) %>%
  mutate(geometry = geometry / 1000) %>%
  bind_cols(as.data.frame(st_coordinates(.))) %>%
  arrange(year, month)

# Subset into training and testing data
# Training 
dat_train <- dat_monthly %>% 
  filter(year %in% 2005:2017)

range(dat_train$year)#2005 - 2017
range(dat_train$time_step)#1 - 156

# Testing
dat_test <- dat_monthly %>% 
  filter(!year %in% 2005:2017)

range(dat_test$year)#2018 - 2019
range(dat_test$time_step)#157 - 180

# Scale predictor variables
# Identify vars to scale
vars_to_scale <- names(dat_train)[
  grepl("^(prcp|tmin|tmax)", names(dat_train)) |
    names(dat_train) %in% c("developed", "cropland", "natural", "forest", "wetlands")
]

# Compute scaling parameters from training data only
scaling_params <- dat_train %>%
  st_drop_geometry() %>%
  summarise(across(all_of(vars_to_scale),
                   list(mean = ~mean(.x, na.rm = TRUE),
                        sd = ~sd(.x, na.rm = TRUE)))) %>%
  pivot_longer(cols = everything(),
               names_to = c("variable", "stat"),
               names_pattern = "^(.*)_(mean|sd)$",
               values_to = "value") %>%
  pivot_wider(names_from = stat, values_from = value)

# Scale training data using its own stats
dat_train <- dat_train %>%
  mutate(across(all_of(vars_to_scale),
                ~ (. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE)))

# Scale test data using training means and sds
dat_test <- dat_test %>%
  mutate(across(all_of(vars_to_scale), 
                ~ {
                  varname <- cur_column()
                  mu <- scaling_params$mean[scaling_params$variable == varname]
                  sigma <- scaling_params$sd[scaling_params$variable == varname]
                  (. - mu) / sigma
                }))

# Convert factor columns
dat_train <- dat_train %>% mutate(across(c("year", "site_id", "county"), as.factor))
dat_test  <- dat_test  %>% mutate(across(c("year", "site_id", "county"), as.factor))

# Save data
write_csv(dat_train, "data/eeev_env_covs_train_scaled.csv")
write_rds(dat_train, "data/eeev_env_covs_train_scaled.rds")
write_csv(dat_test, "data/eeev_env_covs_test_scaled.csv")
write_rds(dat_test, "data/eeev_env_covs_test_scaled.rds")
write_csv(scaling_params, "data/var_scaling_params.csv")
write_csv(dat, "data/eeev_env_covs_filter.csv")
write_rds(dat, "data/eeev_env_covs_filter.rds")

# Read in data
# dat_train <- readRDS("data/eeev_env_covs_train_scaled.rds")
# dat_test <- readRDS("data/eeev_env_covs_test_scaled.rds")

# Checking collinearity
full_model_glmmtmb <- glmmTMB(eeev_prop ~ poly(prcp, 2) + poly(prcp_lag1, 2) + 
                                poly(prcp_lag5, 2) + poly(prcp_lag12, 2) + 
                                poly(tmax, 2) + poly(tmax_lag1, 2) + 
                                poly(tmax_lag6, 2) + poly(tmax_lag12, 2) + 
                                poly(tmin, 2) + poly(tmin_lag1, 2) + poly(tmin_lag6, 2) +
                                poly(tmin_lag12, 2) + poly(developed, 2) + poly(forest,2) + 
                                poly(wetlands, 2) + (1 | county/site_id), 
                              family = binomial(link = "cloglog"),
                              data    = dat_train,
                              weights = weights)

check_collinearity(full_model_glmmtmb)

# Removing tmin, highest VIF
glmmtmb_no_tmin <- glmmTMB(eeev_prop ~ poly(prcp, 2) + poly(prcp_lag1, 2) + 
                                   poly(prcp_lag5, 2) + poly(prcp_lag12, 2) + 
                                   poly(tmax, 2) + poly(tmax_lag1, 2) + 
                                   poly(tmax_lag6, 2) + poly(tmax_lag12, 2) + 
                                   poly(tmin_lag1, 2) + poly(tmin_lag6, 2) + poly(tmin_lag12, 2) +
                                   poly(developed, 2) + poly(forest,2) + 
                                   poly(wetlands, 2) + (1 | county/site_id), 
                                 family = binomial(link = "cloglog"),
                                 data    = dat_train,
                                 weights = weights)

check_collinearity(glmmtmb_no_tmin)

# Removing tmin_lag12, highest VIF
glmmtmb_no_tmin_lag12 <- glmmTMB(eeev_prop ~ poly(prcp, 2) + poly(prcp_lag1, 2) + 
                             poly(prcp_lag5, 2) + poly(prcp_lag12, 2) + 
                             poly(tmax, 2) + poly(tmax_lag1, 2) + 
                             poly(tmax_lag6, 2) + poly(tmax_lag12, 2) + 
                             poly(tmin_lag1, 2) + poly(tmin_lag6, 2) + 
                             poly(developed, 2) + poly(forest,2) + 
                             poly(wetlands, 2) + (1 | county/site_id), 
                           family = binomial(link = "cloglog"),
                           data    = dat_train,
                           weights = weights)

check_collinearity(glmmtmb_no_tmin_lag12)

# Removing tmin_lag6, highest VIF
glmmtmb_no_tmin_lag6 <- glmmTMB(eeev_prop ~ poly(prcp, 2) + poly(prcp_lag1, 2) + 
                                   poly(prcp_lag5, 2) + poly(prcp_lag12, 2) + 
                                   poly(tmax, 2) + poly(tmax_lag1, 2) + 
                                   poly(tmax_lag6, 2) + poly(tmax_lag12, 2) + 
                                   poly(tmin_lag1, 2) + 
                                   poly(developed, 2) + poly(forest,2) + 
                                   poly(wetlands, 2) + (1 | county/site_id), 
                                 family = binomial(link = "cloglog"),
                                 data    = dat_train,
                                 weights = weights)

check_collinearity(glmmtmb_no_tmin_lag6)

# Removing tmax_lag1, highest VIF
glmmtmb_no_tmax_lag1 <- glmmTMB(eeev_prop ~ poly(prcp, 2) + poly(prcp_lag1, 2) + 
                                  poly(prcp_lag5, 2) + poly(prcp_lag12, 2) + 
                                  poly(tmax, 2) + 
                                  poly(tmax_lag6, 2) + poly(tmax_lag12, 2) + 
                                  poly(tmin_lag1, 2) + 
                                  poly(developed, 2) + poly(forest,2) + 
                                  poly(wetlands, 2) + (1 | county/site_id), 
                                family = binomial(link = "cloglog"),
                                data    = dat_train,
                                weights = weights)

check_collinearity(glmmtmb_no_tmax_lag1)


# Fitting a model with dropped variables + landcover
glmmtmb_removed <- glmmTMB(eeev_prop ~ poly(tmin, 2) + poly(tmin_lag6, 2) + 
                             poly(tmin_lag12, 2) + poly(tmax_lag1) +
                             poly(developed, 2) + poly(forest,2) + 
                             poly(wetlands, 2) + (1 | county/site_id), 
                           family = binomial(link = "cloglog"),
                           data    = dat_train,
                           weights = weights)

check_collinearity(glmmtmb_removed)# no multicollinearity

# Save models
glmmtmb_objects <- ls(pattern = "glmmtmb")

save(list = glmmtmb_objects, file = "04_VIF_checks.RData")

# Create a spatial mesh
# Extract locations from training data
locs <- dat_train %>% 
  st_coordinates() %>% 
  as.data.frame()

# Create domain
# Transform crs
domain <- domain %>% 
  st_transform(32617)

st_geometry(domain) <- st_geometry(domain)/1000

# fl_coords <- domain %>%
#   st_buffer(2) %>%
#   st_sample(8000) %>%
#   st_coordinates() %>%
#   as.data.frame() %>% 
#   bind_rows(locs)
# 
# fl_coords_sf <- st_as_sf(fl_coords, coords = c("X", "Y"), crs = 32617, remove = FALSE)
# 
# plot(st_geometry(domain))
# plot(st_geometry(fl_coords_sf))
# 
# domain <- fm_nonconvex_hull(
#   fl_coords_sf,                
#   convex = -0.02,           
#   resolution = c(56, 51))
# 
# #plot to verify
# plot(st_geometry(domain), border = "red", lwd = 2)
# plot(st_geometry(fl_coords_sf), add = TRUE, pch = 16, col = "blue")

# Make spatial mesh
# Calculate max edge parameter
max_edge <- max(
  diff(range(locs$X)),
  diff(range(locs$Y))
)

mesh <- fm_mesh_2d(loc = locs, boundary = domain, max.edge = c(5, 10) * (max_edge / 15),
                   offset = c(max_edge / 1, max_edge / 20),
                   cutoff = 22)

# Check
mesh$n#253 vertices

plot(mesh)
plot(st_geometry(dat_train), add = TRUE)

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

# Make mesh with sdmTMB
mesh_sdmtmb <- make_mesh(data = dat_train, c("X", "Y"), n_knots = 100, mesh = mesh)

# Check
mesh_sdmtmb$mesh$n#253 vertices


# Clearing environment of large objects
rm(dat, dat_monthly, scaling_params, vars_to_scale, min_year, sites_in_domain, 
   sites_sf, site_ids)

# Checking for multicollinearity using conditional index
# Select predictors
predictors <- dat_train %>%
  select(prcp, prcp_lag1, prcp_lag5, prcp_lag12,
         tmax, tmax_lag1, tmax_lag6, tmax_lag12,
         tmin, tmin_lag1, tmin_lag6, tmin_lag12,
         developed, forest, wetlands)


# Create fake outcome and dummy model
lm_fake <- lm(runif(nrow(predictors)) ~ ., data = predictors)

# Run condition index diagnostics
ci_diag <- ols_coll_diag(lm_fake)

# Print CI table
ci_diag$eig_cindex

# Print VIF table
ci_diag$vif_t

# Flagging problematic vars
ci_df <- ci_diag$eig_cindex
high_ci_threshold <- 30
high_vdp_threshold <- 0.5

# Find rows with high condition index
high_ci_rows <- which(ci_df$`Condition Index` > high_ci_threshold)

# For each such row, find variables with VDP >= threshold
problem_vars <- unique(unlist(
  apply(ci_df[high_ci_rows, -(1:2)], 1, function(vdp_row) {
    names(vdp_row)[which(vdp_row >= high_vdp_threshold)]
  })
))

problem_vars#null

## Stepwise backward selection for variable selection
# Define predictor variables these were selected based on VIF checks above, tmin, tmin_lag6, 
# tmin_lag12, and tmax_lag1 were dropped sequentially until all remaining covs had adj.VIF < 5
current_predictors <- c(
  "prcp", "prcp_lag1", "prcp_lag5", "prcp_lag12",
  "tmax", "tmax_lag6", "tmax_lag12",
  "tmin_lag1", "developed", "forest", "wetlands"
)

# Build formula with 2nd order polynomial terms and random intercepts
poly_terms <- paste0("poly(", current_predictors, ", 2)")
full_formula <- as.formula(
  paste("eeev_prop ~", paste(poly_terms, collapse = " + "),
        "+ (1 | county) + (1 | site_id)")
)

# Fit full model
full_model <- sdmTMB(
  formula = full_formula,
  data = dat_train,
  weights = dat_train$weights,
  mesh = mesh_sdmtmb,
  spatial = "on",
  spatiotemporal = "off",
  family = binomial(link = "cloglog"),
  reml = FALSE
)


# Set current model to start backward selection
current_model <- full_model
best_AIC <- AIC(current_model)

# Initialize storage
step_models <- list(full_model = full_model)
step_AICs <- c(full_model = best_AIC)
removal_log <- data.frame(
  step = 0,
  removed = NA,
  AIC = best_AIC,
  delta_AIC = 0,
  predictors = paste(current_predictors, collapse = ","))

# Backward selection loop
improvement <- TRUE
step <- 1

system.time({
  while (improvement && length(current_predictors) > 1) {
    message("Step ", step, ": testing removal of ", length(current_predictors), " predictors")
    improvement <- FALSE
    best_local_AIC <- best_AIC
    best_local_model <- NULL
    best_var_to_remove <- NULL
    
    for (var in current_predictors) {
      message("  Testing removal of: ", var)
      temp_vars <- setdiff(current_predictors, var)
      temp_poly_terms <- paste0("poly(", temp_vars, ", 2)")
      temp_formula <- as.formula(
        paste("eeev_prop ~", paste(temp_poly_terms, collapse = " + "),
              "+ (1 | county) + (1 | site_id)")
      )
      
      temp_model <- tryCatch(
        sdmTMB(
          formula = temp_formula,
          data = dat_train,
          weights = dat_train$weights,
          mesh = mesh_sdmtmb,
          spatial = "on",
          spatiotemporal = "off",
          family = binomial(link = "cloglog"),
          reml = FALSE
        ),
        error = function(e) NULL
      )
      
      if (!is.null(temp_model)) {
        temp_AIC <- AIC(temp_model)
        message("    AIC: ", round(temp_AIC, 2))
        if (temp_AIC < best_local_AIC) {
          best_local_AIC <- temp_AIC
          best_local_model <- temp_model
          best_var_to_remove <- var
        }
      } else {
        message("    model failed")
      }
    }
    
    if (!is.null(best_local_model) && best_local_AIC < best_AIC) {
      message("  Removing ", best_var_to_remove, " improved AIC to ", round(best_local_AIC, 2))
      current_model <- best_local_model
      delta_AIC <- best_AIC - best_local_AIC
      best_AIC <- best_local_AIC
      current_predictors <- setdiff(current_predictors, best_var_to_remove)
      step_name <- paste0("minus_", best_var_to_remove)
      step_models[[step_name]] <- current_model
      step_AICs[step_name] <- best_AIC
      removal_log <- rbind(removal_log, data.frame(
        step = step,
        removed = best_var_to_remove,
        AIC = best_AIC,
        delta_AIC = round(delta_AIC, 2),
        predictors = paste(current_predictors, collapse = ",")
      ))
      improvement <- TRUE
      step <- step + 1
    } else {
      message("  No improvement found at step ", step)
    }
  }
})

# Running the backward selection on the alternate set of predictors. These were the vars that were dropped
# based on VIF. We now fit those variables as an alternate model with landcover vars
current_predictors_alt <- c(
  "tmin", "tmin_lag6", "tmin_lag12", "tmax_lag1", 
  "developed", "forest", "wetlands"
)

# Build formula with 2nd order polynomial terms and random intercepts
poly_terms_alt <- paste0("poly(", current_predictors_alt, ", 2)")
full_formula_alt <- as.formula(
  paste("eeev_prop ~", paste(poly_terms_alt, collapse = " + "),
        "+ (1 | county) + (1 | site_id)")
)

# Fit full model
full_model_alt <- sdmTMB(
  formula = full_formula_alt,
  data = dat_train,
  weights = dat_train$weights,
  mesh = mesh_sdmtmb,
  spatial = "on",
  spatiotemporal = "off",
  family = binomial(link = "cloglog"),
  reml = FALSE
)

# Set current model to start backward selection
current_model_alt <- full_model_alt
best_AIC_alt <- AIC(current_model_alt)

# Initialize storage
step_models_alt <- list(full_model = full_model_alt)
step_AICs_alt <- c(full_model = best_AIC_alt)
removal_log_alt <- data.frame(
  step = 0,
  removed = NA,
  AIC = best_AIC_alt,
  delta_AIC = 0,
  predictors = paste(current_predictors_alt, collapse = ",")
)

# Backward selection loop
improvement_alt <- TRUE
step_alt <- 1

system.time({
  while (improvement_alt && length(current_predictors_alt) > 1) {
    message("Step ", step_alt, ": testing removal of ", length(current_predictors_alt), " predictors")
    improvement_alt <- FALSE
    best_local_AIC_alt <- best_AIC_alt
    best_local_model_alt <- NULL
    best_var_to_remove_alt <- NULL
    
    for (var_alt in current_predictors_alt) {
      message("  Testing removal of: ", var_alt)
      temp_vars_alt <- setdiff(current_predictors_alt, var_alt)
      temp_poly_terms_alt <- paste0("poly(", temp_vars_alt, ", 2)")
      temp_formula_alt <- as.formula(
        paste("eeev_prop ~", paste(temp_poly_terms_alt, collapse = " + "),
              "+ (1 | county) + (1 | site_id)")
      )
      
      temp_model_alt <- tryCatch(
        sdmTMB(
          formula = temp_formula_alt,
          data = dat_train,
          weights = dat_train$weights,
          mesh = mesh_sdmtmb,
          spatial = "on",
          spatiotemporal = "off",
          family = binomial(link = "cloglog"),
          reml = FALSE
        ),
        error = function(e) NULL
      )
      
      if (!is.null(temp_model_alt)) {
        temp_AIC_alt <- AIC(temp_model_alt)
        message("    AIC: ", round(temp_AIC_alt, 2))
        if (temp_AIC_alt < best_local_AIC_alt) {
          best_local_AIC_alt <- temp_AIC_alt
          best_local_model_alt <- temp_model_alt
          best_var_to_remove_alt <- var_alt
        }
      } else {
        message("    model failed")
      }
    }
    
    if (!is.null(best_local_model_alt) && best_local_AIC_alt < best_AIC_alt) {
      message("  Removing ", best_var_to_remove_alt, " improved AIC to ", round(best_local_AIC_alt, 2))
      current_model_alt <- best_local_model_alt
      delta_AIC_alt <- best_AIC_alt - best_local_AIC_alt
      best_AIC_alt <- best_local_AIC_alt
      current_predictors_alt <- setdiff(current_predictors_alt, best_var_to_remove_alt)
      step_name_alt <- paste0("minus_", best_var_to_remove_alt)
      step_models_alt[[step_name_alt]] <- current_model_alt
      step_AICs_alt[step_name_alt] <- best_AIC_alt
      removal_log_alt <- rbind(removal_log_alt, data.frame(
        step = step_alt,
        removed = best_var_to_remove_alt,
        AIC = best_AIC_alt,
        delta_AIC = round(delta_AIC_alt, 2),
        predictors = paste(current_predictors_alt, collapse = ",")
      ))
      improvement_alt <- TRUE
      step_alt <- step_alt + 1
    } else {
      message("  No improvement found at step ", step_alt)
    }
  }
})


# Print variable removal log
removal_log

# Extract top models by AIC
top_model_name <- names(sort(step_AICs))[1]
top_model <- step_models[[top_model_name]]
top_model_name_alt <- names(sort(step_AICs_alt))[1]
top_model_alt <- step_models_alt[[top_model_name_alt]]

summary(top_model)
sanity(top_model)
top_model$model$convergence


summary(top_model_alt)
sanity(top_model_alt)
top_model_alt$model$convergence


# Fitting intercept only models
sdmtmb_intercept <- sdmTMB(
  formula = eeev_prop ~ 1,
  data = dat_train,
  weights = dat_train$weights,
  mesh = mesh_sdmtmb,
  spatial = "on",
  spatiotemporal = "off",
  family = binomial(link = "cloglog"),
  reml = FALSE
)

sdmtmb_ranef <- sdmTMB(
  formula = eeev_prop ~ (1 | county) + (1 | site_id),
  data = dat_train,
  weights = dat_train$weights,
  mesh = mesh_sdmtmb,
  spatial = "on",
  spatiotemporal = "off",
  family = binomial(link = "cloglog"),
  reml = FALSE
)

glmmtmb_intercept <- sdmTMB(
  formula = eeev_prop ~ 1,
  data = dat_train,
  weights = dat_train$weights,
  mesh = mesh_sdmtmb,
  spatial = "off",
  spatiotemporal = "off",
  family = binomial(link = "cloglog"),
  reml = FALSE
)

glmmtmb_ranef <- sdmTMB(
  formula = eeev_prop ~ (1 | county) + (1 | site_id),
  data = dat_train,
  weights = dat_train$weights,
  mesh = mesh_sdmtmb,
  spatial = "off",
  spatiotemporal = "off",
  family = binomial(link = "cloglog"),
  reml = FALSE
)

AIC(full_model, top_model_objects[[1]], top_model_objects[[2]], sdmtmb_intercept, 
    sdmtmb_ranef, glmmtmb_intercept, glmmtmb_ranef)

# Checking multicollinearity
top_model_glmm <- glmmTMB(eeev_prop ~ poly(prcp_lag1, 2) + 
                            poly(prcp_lag5, 2) + poly(prcp_lag12, 2) + 
                            poly(tmax_lag6, 2) + poly(tmax_lag12, 2) + 
                            poly(tmin_lag1, 2) + poly(forest, 2) + 
                            poly(wetlands, 2) + (1 | county/site_id),
                            family = binomial(link = "cloglog"),
                            data    = dat_train,
                            weights = weights)

check_collinearity(top_model_glmm)#all good

# Simulating residuals and checking for temporal and spatial autocorrelation
resid_top_model_glmm <- simulateResiduals(top_model_glmm)

# Spatial autocorrelation
testSpatialAutocorrelation(resid_top_model_glmm, x = dat_train$X, y = dat_train$Y)

# Temporal autocorrelation
check_autocorrelation(top_model_glmm)

# Dropping collinear variables
top_model_noncol_glmm <- glmmTMB(eeev_prop ~ poly(prcp_lag12, 2) + poly(tmax_lag1, 2) +
                                    poly(tmin_lag6, 2) + poly(forest, 2) +
                                    poly(wetlands, 2) + (1 | county/site_id),
                                  family = binomial(link = "cloglog"),
                                  data    = dat_train,
                                  weights = testing)

# Check again
check_collinearity(top_model_noncol_glmm)#all good

# Fitting this model in sdmTMB
 system.time(
 eeev_spatial_mod <- sdmTMB(
   formula = eeev_prop ~ poly(prcp_lag12, 2) + poly(tmax_lag1, 2) +
     poly(tmin_lag6, 2) + poly(forest, 2) +
     poly(wetlands, 2) + (1 | county) + (1 | site_id),
   data = dat_train,
   weights = dat_train$testing,
   mesh = mesh_sdmtmb,
   spatial = "on",
   spatiotemporal = "off",
   family = binomial(link = "cloglog"),
   time = "time_step",
   reml = FALSE
 )
 )

# Extracting top model from backward selection that includes precipitation
formula_mod_back <- formula(top_model)

# There were convergence issues with the top model, trying without the spatial term, essentially a glmmtmb
top_model_non_spatial <- sdmTMB(
  formula = formula_mod_back,
  data = dat_train,
  weights = dat_train$weights,
  mesh = mesh_sdmtmb,
  spatial = "off",
  spatiotemporal = "off",
  family = binomial(link = "cloglog"),
  reml = FALSE
)

sanity(top_model_non_spatial)#better, no warnings
top_model_non_spatial$model$convergence
summary(top_model_non_spatial)

AIC(top_model, top_model_non_spatial)#non spatial is slightly better

# Trying to adjust the mesh sequentially from 250 - 1000 to assess convergence
# 500 vert mesh
mesh_500 <- fm_mesh_2d(loc = locs, boundary = domain, max.edge = c(5, 10) * (max_edge / 15),
                   offset = c(max_edge / 1, max_edge / 20),
                   cutoff = 12)

mesh_500$n

# Make mesh with sdmTMB
mesh_sdmtmb_500 <- make_mesh(data = dat_train, c("X", "Y"), n_knots = 100, mesh = mesh_500)

# Check
mesh_sdmtmb_500$mesh$n#543 vertices

plot(mesh_sdmtmb_500)

# Refit model with new mesh
system.time(
top_model_500 <- sdmTMB(
  formula = formula_mod_back,
  data = dat_train,
  weights = dat_train$weights,
  mesh = mesh_sdmtmb_500,
  spatial = "on",
  spatiotemporal = "off",
  family = binomial(link = "cloglog"),
  reml = FALSE
)
)#5.2 min

sanity(top_model_500)#seems better
top_model_500$model$convergence
summary(top_model_500)

# 750 vert mesh
mesh_750 <- fm_mesh_2d(loc = locs, boundary = domain, max.edge = c(5, 10) * (max_edge / 15),
                       offset = c(max_edge / 1, max_edge / 20),
                       cutoff = 8.5)

mesh_750$n

# Make mesh with sdmTMB
mesh_sdmtmb_750 <- make_mesh(data = dat_train, c("X", "Y"), n_knots = 100, mesh = mesh_750)

# Check
mesh_sdmtmb_750$mesh$n#761

plot(mesh_sdmtmb_750)

# Refit model with new mesh
system.time(
  top_model_750 <- sdmTMB(
    formula = formula_mod_back,
    data = dat_train,
    weights = dat_train$weights,
    mesh = mesh_sdmtmb_750,
    spatial = "on",
    spatiotemporal = "off",
    family = binomial(link = "cloglog"),
    reml = FALSE
  )
)#5 min

sanity(top_model_750)
top_model_750$model$convergence
summary(top_model_750)

# 1000 vert mesh
mesh_1000 <- fm_mesh_2d(loc = locs, boundary = domain, max.edge = c(5, 10) * (max_edge / 15),
                       offset = c(max_edge / 1, max_edge / 20),
                       cutoff = 6)

mesh_1000$n

# Make mesh with sdmTMB
mesh_sdmtmb_1000 <- make_mesh(data = dat_train, c("X", "Y"), n_knots = 100, mesh = mesh_1000)

# Check
mesh_sdmtmb_1000$mesh$n#1083

plot(mesh_sdmtmb_1000)

# Refit model with new mesh
system.time(
  top_model_1000 <- sdmTMB(
    formula = formula_mod_back,
    data = dat_train,
    weights = dat_train$weights,
    mesh = mesh_sdmtmb_1000,
    spatial = "on",
    spatiotemporal = "off",
    family = binomial(link = "cloglog"),
    reml = FALSE
  )
)#5 min

sanity(top_model_1000)
top_model_1000$model$convergence
summary(top_model_1000)

# From the above steps, the best mesh size to use is 750. We will rerun a backward selection on the 
# top model from the previous coarse mesh backward selection to further reduce predictors
current_predictors_750 <- c(
  "prcp_lag1", "prcp_lag5", "prcp_lag12",
  "tmax_lag6", "tmax_lag12",
  "tmin_lag1", "forest", "wetlands"
)

# Build formula with 2nd order polynomial terms and random intercepts
poly_terms_750 <- paste0("poly(", current_predictors_750, ", 2)")
full_formula_750 <- as.formula(
  paste("eeev_prop ~", paste(poly_terms_750, collapse = " + "),
        "+ (1 | county) + (1 | site_id)")
)

# Fit full model
system.time(
full_model_750 <- sdmTMB(
  formula = full_formula_750,
  data = dat_train,
  weights = dat_train$weights,
  mesh = mesh_sdmtmb_750,
  spatial = "on",
  spatiotemporal = "off",
  family = binomial(link = "cloglog"),
  reml = FALSE
)
)#5 min


# Set current model to start backward selection
current_model_750 <- full_model_750
best_AIC_750 <- AIC(current_model_750)

# Initialize storage
step_models_750 <- list(full_model_750 = full_model_750)
step_AICs_750 <- c(full_model_750 = best_AIC_750)
removal_log_750 <- data.frame(
  step = 0,
  removed = NA,
  AIC = best_AIC_750,
  delta_AIC = 0,
  predictors = paste(current_predictors_750, collapse = ",")
)

# Backward selection loop
improvement_750 <- TRUE
step_750 <- 1

system.time({
  while (improvement_750 && length(current_predictors_750) > 1) {
    message("Step ", step_750, ": testing removal of ", length(current_predictors_750), " predictors")
    improvement_750 <- FALSE
    best_local_AIC_750 <- best_AIC_750
    best_local_model_750 <- NULL
    best_var_to_remove_750 <- NULL
    
    for (var in current_predictors_750) {
      message("  Testing removal of: ", var)
      temp_vars_750 <- setdiff(current_predictors_750, var)
      temp_poly_terms_750 <- paste0("poly(", temp_vars_750, ", 2)")
      temp_formula_750 <- as.formula(
        paste("eeev_prop ~", paste(temp_poly_terms_750, collapse = " + "),
              "+ (1 | county) + (1 | site_id)")
      )
      
      temp_model_750 <- tryCatch(
        sdmTMB(
          formula = temp_formula_750,
          data = dat_train,
          weights = dat_train$weights,
          mesh = mesh_sdmtmb_750,
          spatial = "on",
          spatiotemporal = "off",
          family = binomial(link = "cloglog"),
          reml = FALSE
        ),
        error = function(e) NULL
      )
      
      if (!is.null(temp_model_750)) {
        temp_AIC_750 <- AIC(temp_model_750)
        message("    AIC: ", round(temp_AIC_750, 2))
        if (temp_AIC_750 < best_local_AIC_750) {
          best_local_AIC_750 <- temp_AIC_750
          best_local_model_750 <- temp_model_750
          best_var_to_remove_750 <- var
        }
      } else {
        message("    model failed")
      }
    }
    
    if (!is.null(best_local_model_750) && best_local_AIC_750 < best_AIC_750) {
      message("  Removing ", best_var_to_remove_750, " improved AIC to ", round(best_local_AIC_750, 2))
      current_model_750 <- best_local_model_750
      delta_AIC_750 <- best_AIC_750 - best_local_AIC_750
      best_AIC_750 <- best_local_AIC_750
      current_predictors_750 <- setdiff(current_predictors_750, best_var_to_remove_750)
      step_name_750 <- paste0("minus_", best_var_to_remove_750)
      step_models_750[[step_name_750]] <- current_model_750
      step_AICs_750[step_name_750] <- best_AIC_750
      removal_log_750 <- rbind(removal_log_750, data.frame(
        step = step_750,
        removed = best_var_to_remove_750,
        AIC = best_AIC_750,
        delta_AIC = round(delta_AIC_750, 2),
        predictors = paste(current_predictors_750, collapse = ",")
      ))
      improvement_750 <- TRUE
      step_750 <- step_750 + 1
    } else {
      message("  No improvement found at step ", step_750)
    }
  }
})

# Rerunning the backward selection on the alternate set of predictors (ones dropped based on VIF),
# the first round dropped developed with a mesh of 250 vertices, rerunning the backward selection on
# these predictors with a mesh of 750 vertices
current_predictors_alt_750 <- c(
  "tmin", "tmin_lag6", "tmin_lag12", "tmax_lag1", 
  "forest", "wetlands"
)

# Build formula with 2nd order polynomial terms and random intercepts
poly_terms_alt_750 <- paste0("poly(", current_predictors_alt_750, ", 2)")
full_formula_alt_750 <- as.formula(
  paste("eeev_prop ~", paste(poly_terms_alt_750, collapse = " + "),
        "+ (1 | county) + (1 | site_id)")
)

# Fit full model
full_model_alt_750 <- sdmTMB(
  formula = full_formula_alt_750,
  data = dat_train,
  weights = dat_train$weights,
  mesh = mesh_sdmtmb_750,
  spatial = "on",
  spatiotemporal = "off",
  family = binomial(link = "cloglog"),
  reml = FALSE
)

# Set current model to start backward selection
current_model_alt_750 <- full_model_alt_750
best_AIC_alt_750 <- AIC(current_model_alt_750)

# Initialize storage
step_models_alt_750 <- list(full_model = full_model_alt_750)
step_AICs_alt_750 <- c(full_model = best_AIC_alt_750)
removal_log_alt_750 <- data.frame(
  step = 0,
  removed = NA,
  AIC = best_AIC_alt_750,
  delta_AIC = 0,
  predictors = paste(current_predictors_alt_750, collapse = ",")
)

# Backward selection loop
improvement_alt_750 <- TRUE
step_alt_750 <- 1

system.time({
  while (improvement_alt_750 && length(current_predictors_alt_750) > 1) {
    message("Step ", step_alt_750, ": testing removal of ", length(current_predictors_alt_750), " predictors")
    improvement_alt_750 <- FALSE
    best_local_AIC_alt_750 <- best_AIC_alt_750
    best_local_model_alt_750 <- NULL
    best_var_to_remove_alt_750 <- NULL
    
    for (var_alt_750 in current_predictors_alt_750) {
      message("  Testing removal of: ", var_alt_750)
      temp_vars_alt_750 <- setdiff(current_predictors_alt_750, var_alt_750)
      temp_poly_terms_alt_750 <- paste0("poly(", temp_vars_alt_750, ", 2)")
      temp_formula_alt_750 <- as.formula(
        paste("eeev_prop ~", paste(temp_poly_terms_alt_750, collapse = " + "),
              "+ (1 | county) + (1 | site_id)")
      )
      
      temp_model_alt_750 <- tryCatch(
        sdmTMB(
          formula = temp_formula_alt_750,
          data = dat_train,
          weights = dat_train$weights,
          mesh = mesh_sdmtmb_750,
          spatial = "on",
          spatiotemporal = "off",
          family = binomial(link = "cloglog"),
          reml = FALSE
        ),
        error = function(e) NULL
      )
      
      if (!is.null(temp_model_alt_750)) {
        temp_AIC_alt_750 <- AIC(temp_model_alt_750)
        message("    AIC: ", round(temp_AIC_alt_750, 2))
        if (temp_AIC_alt_750 < best_local_AIC_alt_750) {
          best_local_AIC_alt_750 <- temp_AIC_alt_750
          best_local_model_alt_750 <- temp_model_alt_750
          best_var_to_remove_alt_750 <- var_alt_750
        }
      } else {
        message("    model failed")
      }
    }
    
    if (!is.null(best_local_model_alt_750) && best_local_AIC_alt_750 < best_AIC_alt_750) {
      message("  Removing ", best_var_to_remove_alt_750, " improved AIC to ", round(best_local_AIC_alt_750, 2))
      current_model_alt_750 <- best_local_model_alt_750
      delta_AIC_alt_750 <- best_AIC_alt_750 - best_local_AIC_alt_750
      best_AIC_alt_750 <- best_local_AIC_alt_750
      current_predictors_alt_750 <- setdiff(current_predictors_alt_750, best_var_to_remove_alt_750)
      step_name_alt_750 <- paste0("minus_", best_var_to_remove_alt_750)
      step_models_alt_750[[step_name_alt_750]] <- current_model_alt_750
      step_AICs_alt_750[step_name_alt_750] <- best_AIC_alt_750
      removal_log_alt_750 <- rbind(removal_log_alt_750, data.frame(
        step = step_alt_750,
        removed = best_var_to_remove_alt_750,
        AIC = best_AIC_alt_750,
        delta_AIC = round(delta_AIC_alt_750, 2),
        predictors = paste(current_predictors_alt_750, collapse = ",")
      ))
      improvement_alt_750 <- TRUE
      step_alt_750 <- step_alt_750 + 1
    } else {
      message("  No improvement found at step ", step_alt_750)
    }
  }
})

# No variables dropped for both the primary and alternate models during the second round of backward
# selection with a 750 vert mesh. Moving forward with these candidate sets for final model fitting.

# Save important objects
save.image("04_variable_reduction_model_selection.RData")
write_rds(domain, "output/mesh_domain.rds")
write_rds(mesh_sdmtmb_750, "output/mesh_sdmTMB_750.rds")