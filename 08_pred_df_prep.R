# Clear environment and perform garbage collection
rm(list = ls())
gc()

# Load required libraries
library(tidyverse)
library(sf)
library(terra)
library(tigris)
library(rmapshaper)
library(rnaturalearth)
library(sdmTMB)

# Read in scaling parameters
scaling_params <- read_csv("data/var_scaling_params.csv")

# Read and process spatial domain
domain <- read_rds("data/fl_poly_trim.rds")
domain <- st_make_valid(domain)
domain <- domain %>% 
  st_transform(32617)

#st_geometry(domain) <- st_geometry(domain)/1000

blank <- rast("data/tmin_stack_monthly.tif")$`1-2016` > Inf

blank <- blank %>%
  crop(st_geometry(domain), mask = T)

# Read in monthly climate data and mask to florida boundary
tmin_monthly_raw <- rast("data/tmin_stack_monthly.tif") %>%
  crop(blank, mask = T)

prcp_monthly_raw <- rast("data/prcp_stack_monthly.tif") %>%
  crop(blank, mask = T)

tmax_monthly_raw <- rast("data/tmax_stack_monthly.tif") %>%
  crop(blank, mask = T)

# Create climate dataframe from rasters
climate_monthly_df <- prcp_monthly_raw %>%
  as.data.frame(xy = T) %>%
  pivot_longer(-c(x, y), names_to = "layer") %>%
  rename("prcp" = "value",
         "date" = "layer") %>%
  bind_cols(
    tmax_monthly_raw %>%
      as.data.frame(xy = T) %>%
      pivot_longer(-c(x, y), names_to = "layer") %>%
      select(value) %>%
      rename("tmax" = "value"),
    tmin_monthly_raw %>%
      as.data.frame(xy = T) %>%
      pivot_longer(-c(x, y), names_to = "layer") %>%
      select(value) %>%
      rename("tmin" = "value")) %>%
  mutate(prcp_lag1 = lag(prcp, 1), 
         prcp_lag5 = lag(prcp, 5),
         prcp_lag6 = lag(prcp, 6),
         prcp_lag12 = lag(prcp, 12),
         tmax_lag1 = lag(tmax, 1), 
         tmax_lag5 = lag(tmax, 5),
         tmax_lag6 = lag(tmax, 6),
         tmax_lag12 = lag(tmax, 12),
         tmin_lag1 = lag(tmin, 1), 
         tmin_lag5 = lag(tmin, 5),
         tmin_lag6 = lag(tmin, 6),
         tmin_lag12 = lag(tmin, 12),
         month = str_split_fixed(date, "-", 2)[,1],
         year = str_split_fixed(date, "-", 2)[,2],
         month = as.character(month),
         year = as.character(year)) %>%
  filter(year %in% c(2001:2021)) %>%
  rename("lon" = "x", "lat" = "y") %>%
  mutate(lat = lat / 1000,
         lon = lon / 1000) %>%
  select(lon, lat, year, month, 
         prcp, prcp_lag1, prcp_lag5, prcp_lag6, prcp_lag12, 
         tmax, tmax_lag1, tmax_lag5, tmax_lag6, tmax_lag12, 
         tmin, tmin_lag1, tmin_lag5, tmin_lag6, tmin_lag12)

# Read in landcover data
forest <- rast("data/forest_annual.tif") %>% 
  crop(blank)

developed <- rast("data/developed_annual.tif") %>% 
  crop(blank)

wetlands <- rast("data/wetlands_annual.tif") %>% 
  crop(blank)

natural <- rast("data/natural_annual.tif") %>% 
  crop(blank)

# Create data frame
landcover_df <- as.data.frame(developed, xy = T) %>%
  pivot_longer(-c(x, y), names_to = "layer") %>%
  rename("developed" = "value") %>%
  bind_cols(
    wetlands %>%
      as.data.frame(xy = T) %>%
      pivot_longer(-c(x, y), names_to = "layer") %>%
      select(value) %>%
      rename("wetlands" = "value"),
    natural %>%
      as.data.frame(xy = T) %>%
      pivot_longer(-c(x, y), names_to = "layer") %>%
      select(value) %>%
      rename("natural" = "value"),
    forest %>%
      as.data.frame(xy = T) %>%
      pivot_longer(-c(x, y), names_to = "layer") %>%
      select(value) %>%
      rename("forest" = "value")) %>%
  rename("year" = "layer", "lon" = "x", "lat" = "y") %>%
  mutate(lat = lat / 1000,
         lon = lon / 1000) %>%
  filter(year %in% c(2001:2021)) %>%
  select(lon, lat, year, 
         forest, developed, wetlands, natural) 

landcover_monthly_df <- landcover_df %>%
  sdmTMB::replicate_df(time_name = "month", time_values = 1:12) %>%
  mutate(month = as.character(month),
         year = as.character(year)) %>%
  select(lon, lat, year, month, forest, developed, wetlands, natural) 

# Join with climate
env_monthly_df <- climate_monthly_df %>% 
  left_join(landcover_monthly_df) %>% 
  na.omit() %>% 
  select(lon, lat, year, month, prcp, prcp_lag1, prcp_lag5, prcp_lag6, prcp_lag12, 
  tmax, tmax_lag1, tmax_lag5, tmax_lag6, tmax_lag12, tmin, tmin_lag1, tmin_lag5, 
  tmin_lag6, tmin_lag12, forest, developed, wetlands, natural)

# Scale based on training data mean and sd
env_monthly_scaled <- env_monthly_df %>%
  mutate(
    across(all_of(scaling_params$variable[scaling_params$variable %in% names(.)]), 
           ~ (. - scaling_params$mean[scaling_params$variable == cur_column()]) / 
             scaling_params$sd[scaling_params$variable == cur_column()])
  )

# Checks
scaled_vars <- scaling_params$variable[scaling_params$variable %in% names(env_monthly_scaled)]

env_monthly_scaled %>%
  summarise(across(all_of(scaled_vars), list(mean = mean, sd = sd), .names = "{.col}_{.fn}"))

original_prcp <- env_monthly_df$prcp[1:3]
scaled_prcp <- env_monthly_scaled$prcp[1:3]
manual_scaled <- (original_prcp - scaling_params$mean[scaling_params$variable == "prcp"]) / 
  scaling_params$sd[scaling_params$variable == "prcp"]

data.frame(
  original = original_prcp,
  scaled = scaled_prcp,
  manual = manual_scaled,
  match = abs(scaled_prcp - manual_scaled) < 1e-10
)#all good

# Save data
write_rds(env_monthly_scaled, "data/pred_df.rds")
write_rds(climate_monthly_df, "data/climate_pred_df.rds")
write_rds(landcover_monthly_df, "data/landcover_pred_df.rds")
save.image("09_pred_df_prep.RData")
