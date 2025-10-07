#c Clear environment and perform garbage collection 
rm(list = ls())
gc()

# Load required libraries
library(tidyverse)
library(data.table)
library(terra)
library(sf)


# Read in data
dat <- read_csv("data/compiled_eeev_data.csv")

# Convert to data.table
dat <- as.data.table(dat)

# Cleaning weekly data
dat_weekly <- dat[, .(
  testing = max(testing, na.rm = TRUE),
  n_positive_eeev = sum(n_positive_eeev, na.rm = TRUE),
  lon = unique(lon),
  lat = unique(lat)
), by = .(county, site_id, year, month, week)]


# Aggregate data by month
dat_eeev_monthly <- dat_weekly[, .(
  testing = sum(testing, na.rm = TRUE), 
  eeev = sum(n_positive_eeev, na.rm = TRUE)
), by = .(year, month, county, site_id, lon, lat)]

# Read in stacked monthly rasters
precip <- rast("environ/prcp_stack_monthly.tif")
tmax <- rast("environ/tmax_stack_monthly.tif")
tmin <- rast("environ/tmin_stack_monthly.tif")

# Clean EEEV data
dat_eeev_clean <- dat_eeev_monthly[!is.na(lon) & !is.na(lat)]
sites <- unique(dat_eeev_clean[, .(site_id, lon, lat)])

# Create spatial vector
coords <- vect(sites[, .(lon, lat)], crs = "EPSG:4326") %>% 
  project("EPSG:32617")

# Extract climate data
prcp_extract <- extract(precip, coords, method = "bilinear", xy = TRUE, ID = FALSE) %>%
  bind_cols(sites) %>%
  pivot_longer(
    cols = matches("\\d+-\\d{4}"),
    names_to = c("month", "year"),
    names_sep = "-",
    names_transform = list(month = as.integer, year = as.integer),
    values_to = "prcp")

tmax_extract <- extract(tmax, coords, method = "bilinear", xy = TRUE, ID = FALSE) %>%
  bind_cols(sites) %>%
  pivot_longer(
    cols = matches("\\d+-\\d{4}"),
    names_to = c("month", "year"),
    names_sep = "-",
    names_transform = list(month = as.integer, year = as.integer),
    values_to = "tmax")

tmin_extract <- extract(tmin, coords, method = "bilinear", xy = TRUE, ID = FALSE) %>%
  bind_cols(sites) %>%
  pivot_longer(
    cols = matches("\\d+-\\d{4}"),
    names_to = c("month", "year"),
    names_sep = "-",
    names_transform = list(month = as.integer, year = as.integer),
    values_to = "tmin")

# Combine climate data
clim_extract <- prcp_extract %>%
  left_join(tmax_extract[, c("site_id", "month", "year", "tmax")], by = c("site_id", "month", "year")) %>%
  left_join(tmin_extract[, c("site_id", "month", "year", "tmin")], by = c("site_id", "month", "year"))

# Add lagged dates
data_monthly <- dat_eeev_clean

data_monthly <- data_monthly %>%
  mutate(across(month, as.integer)) %>%
  mutate(date = make_date(year, month)) %>%
  mutate(across(date, list(
    lag1 = ~ .x %m-% months(1),
    lag5 = ~ .x %m-% months(5),
    lag6 = ~ .x %m-% months(6),
    lag12 = ~ .x %m-% months(12)
  ))) %>%
  mutate(across(starts_with("date_lag"), list(
    month = ~ month(.x),
    year = ~ year(.x)
  ), .names = "{.col}_{.fn}"))

# Join current and lagged climate data
data_monthly_clim <- data_monthly %>%
  left_join(clim_extract, by = c("site_id", "month", "year")) %>%
  rename(prcp = prcp, tmax = tmax, tmin = tmin)

data_monthly_clim <- data_monthly_clim %>%
  left_join(clim_extract, by = c(
    "site_id" = "site_id",
    "date_lag1_month" = "month",
    "date_lag1_year" = "year"
  )) %>% 
  rename(
    prcp_lag1 = prcp.y, 
    tmax_lag1 = tmax.y,
    tmin_lag1 = tmin.y
  )

data_monthly_clim <- data_monthly_clim %>%
  left_join(clim_extract, by = c(
    "site_id" = "site_id",
    "date_lag5_month" = "month",
    "date_lag5_year" = "year"
  )) %>%
  rename(
    prcp_lag5 = prcp,
    tmax_lag5 = tmax,
    tmin_lag5 = tmin
  )

data_monthly_clim <- data_monthly_clim %>%
  left_join(clim_extract, by = c(
    "site_id" = "site_id",
    "date_lag6_month" = "month",
    "date_lag6_year" = "year"
  )) %>%
  rename(
    prcp_lag6 = prcp,
    tmax_lag6 = tmax,
    tmin_lag6 = tmin
  )

data_monthly_clim <- data_monthly_clim %>%
  left_join(clim_extract, by = c(
    "site_id" = "site_id",
    "date_lag12_month" = "month",
    "date_lag12_year" = "year"
  )) %>%
  rename(
    prcp_lag12 = prcp,
    tmax_lag12 = tmax,
    tmin_lag12 = tmin
  )

# Cleaning
data_monthly_clim_clean <- data_monthly_clim %>%
  select(year, month, county, site_id, lon = lon.x, lat = lat.x,
         testing, eeev, 
         prcp = prcp.x, tmax = tmax.x, tmin = tmin.x,
         prcp_lag1, tmax_lag1, tmin_lag1,
         prcp_lag5, tmax_lag5, tmin_lag5,
         prcp_lag6, tmax_lag6, tmin_lag6,
         prcp_lag12, tmax_lag12, tmin_lag12)

data_monthly_clim_clean <- data_monthly_clim_clean %>%
  arrange(site_id, year, month)

# Sanity check
alachua_1_sanity_check <- data_monthly_clim_clean %>%
  filter(site_id == "Alachua_1") %>%
  arrange(year, month) %>%
  select(year, month, prcp, prcp_lag1, prcp_lag5, tmax, tmax_lag5) %>%
  mutate(
    prcp_lag1_expected = lag(prcp, 1),
    prcp_lag5_expected = lag(prcp, 5),
    tmax_lag5_expected = lag(tmax, 5)
  ) %>%
  select(year, month, prcp_lag1, prcp_lag1_expected,
         prcp_lag5, prcp_lag5_expected,
         tmax_lag5, tmax_lag5_expected)#seems good

# Read in nlcd data
nlcd_prop <- read_csv("environ/site_landcover_prop_2.csv")

# Join with climate data
data_monthly_nlcd <- data_monthly_clim_clean %>%
  left_join(nlcd_prop, by = c("site_id" = "id", "year" = "year"))

# Create a column that calculates the proportion of positive eeev cases, asuming 6 birds were tested per testing event
data_monthly_nlcd <- data_monthly_nlcd[, eeev_prop := fifelse(testing > 0, 
                                                              eeev / (6 * testing), 0)]
# Convert to sf
data_monthly_nlcd_sf <- st_as_sf(
  data_monthly_nlcd,
  coords = c("lon", "lat"),
  crs = 4326)

# Write data
write_csv(data_monthly_nlcd_sf, "data/eeev_env_covs.csv")
write_rds(data_monthly_nlcd_sf, "data/eeev_env_covs.rds")