# Load required packages
library(sf)
library(terra)
library(spdep)
library(tidyverse)
library(ebirdst)
library(patchwork)

#read in statewide predictions
eeev_mod_statewide_preds <- read_rds("output/eeev_model_primary_state_preds.rds")

# Read in prediction data frame
pred_df <- read_rds("data/pred_df.rds")

# Assemble pred_df 
pred_df_prepared <- pred_df %>%
  filter(year %in% as.character(2005:2019)) %>%
  mutate(
    county = NA,
    site_id = NA) %>%
  left_join(
    expand_grid(
      year = as.character(2005:2019),
      month = as.character(1:12)
    ) %>%
      arrange(year, month) %>%
      mutate(time_step = row_number()),  
    by = c("year", "month")
  ) %>%
  rename(X = lon, Y = lat) %>%
  mutate(year = as.factor(year)); rm(pred_df)

# Process statewide predictions
eeev_mod_monthly_statewide <- pred_df_prepared %>%
  select(X, Y, year, month) %>%
  bind_cols(
    data.frame(
      mean = rowMeans(eeev_mod_statewide_preds),
      se = apply(eeev_mod_statewide_preds, 1, FUN = function(x){sd(x)/sqrt(100)}))) %>%
  mutate(
    # Convert from cloglog link space to response space
    mean = 1 - exp(-exp(mean))); rm(eeev_mod_statewide_preds)

# Create a mean monthly raster to compare with bird data
# Create mean monthly data (averaged across 2005-2019)
eeev_mean_monthly <- eeev_mod_monthly_statewide %>%
  group_by(X, Y, month) %>%
  summarise(
    mean_eeev = mean(mean, na.rm = TRUE),
    se_eeev = mean(se, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(month = as.numeric(month), 
         X_utm = X * 1000,
         Y_utm = Y * 1000)

# Create raster for each month
create_monthly_raster <- function(month_num) {
  month_data <- eeev_mean_monthly %>%
    filter(month == month_num) %>%
    select(X_utm, Y_utm, mean_eeev) %>%
    rename(X = X_utm, Y = Y_utm)
  
  # Convert to terra raster
  r <- rast(month_data, type = "xyz", crs = "EPSG:32617")  
  names(r) <- paste0("EEEV_month_", sprintf("%02d", month_num))
  return(r)
}

# Create all 12 monthly rasters
eeev_monthly_rasters <- map(1:12, create_monthly_raster)
names(eeev_monthly_rasters) <- paste0("month_", sprintf("%02d", 1:12))

# Convert to geographic coordinates (lat/lon)
eeev_monthly_rasters_latlon <- map(eeev_monthly_rasters, ~{
  project(.x, "EPSG:4326")  
})

# Check the results
eeev_monthly_rasters[[1]]  
plot(eeev_monthly_rasters[[1]], main = "EEEV Month 1 (UTM)")

eeev_monthly_rasters_latlon[[1]]    
plot(eeev_monthly_rasters_latlon[[1]], main = "EEEV Month 1 (Lat/Lon)")

# Save both versions
# UTM version
walk2(eeev_monthly_rasters, 1:12, ~{
  writeRaster(.x, 
              filename = paste0("output/rasters/eeev_mean_month_", sprintf("%02d", .y), "_utm.tif"),
              overwrite = TRUE)
})

# Lat/lon version 
walk2(eeev_monthly_rasters_latlon, 1:12, ~{
  writeRaster(.x, 
              filename = paste0("output/rasters/eeev_mean_month_", sprintf("%02d", .y), "_latlon.tif"),
              overwrite = TRUE)
})


#download ebird status data
ebirdst_runs %>% 
       filter(str_detect(common_name, 
                         "Yellow-rumped Warbler|Eastern Phoebe|Hermit Thrush|Pine Warbler|Green Heron|Wood Thrush|Common Yellowthroat|White-eyed Vireo|Red-eyed Vireo|Black-crowned Night Heron|American Robin|Northern Cardinal")) %>%
       select(species_code, common_name, scientific_name)

species_info <- tibble(
  common_name = c("Yellow-rumped Warbler", "Eastern Phoebe", "Hermit Thrush", 
                  "Black-crowned Night Heron", "Pine Warbler", "Green Heron",
                  "Wood Thrush", "Common Yellowthroat", "White-eyed Vireo", "Red-eyed Vireo", "American Robin", "Northern Cardinal"),
  species_code = c("yerwar", "easpho", "herthr", "bcnher", "pinwar", 
                   "grnher", "woothr", "comyel", "whevir", "reevir1", "amerob", "norcar")  
)

# Set up data directory
if (!dir.exists("data/ebirdst")) dir.create("data/ebirdst", recursive = TRUE)

# Function to download bird data
download_species_data <- function(species_code, species_name) {
  cat("Downloading annual abundance for", species_name, "(", species_code, ")...\n")
  try({
    download_path <- ebirdst_download_status(
      species = species_code,
      path = "data/ebirdst/",
      pattern = "abundance_(median|upper|lower)_3km",
      download_abundance = TRUE,
      show_progress = TRUE
    )
    cat("Successfully downloaded", species_name, "\n\n")
    return(download_path)
  }, silent = FALSE)
}

download_results <- map2(species_info$species_code, species_info$common_name,
                         ~download_species_data(.x, .y))

# Load data
yerwar_data <- load_raster(product = "abundance", species = "yerwar", path = "data/ebirdst/", metric = "median")

easpho_data <- load_raster(product = "abundance", species = "easpho", path = "data/ebirdst/", metric = "median")

herthr_data <- load_raster(product = "abundance", species = "herthr", path = "data/ebirdst/", metric = "median")

bcnher_data <- load_raster(product = "abundance", species = "bcnher", path = "data/ebirdst/", metric = "median")

pinwar_data <- load_raster(product = "abundance", species = "pinwar", path = "data/ebirdst/", metric = "median")

grnher_data <- load_raster(product = "abundance", species = "grnher", path = "data/ebirdst/", metric = "median")

woothr_data <- load_raster(product = "abundance", species = "woothr", path = "data/ebirdst/", metric = "median")

comyel_data <- load_raster(product = "abundance", species = "comyel", path = "data/ebirdst/", metric = "median")

whevir_data <- load_raster(product = "abundance", species = "whevir", path = "data/ebirdst/", metric = "median")

reevir1_data <- load_raster(product = "abundance", species = "reevir1", path = "data/ebirdst/", metric = "median")

amerob_data <- load_raster(product = "abundance", species = "amerob", path = "data/ebirdst/", metric = "median")

norcar_data <- load_raster(product = "abundance", species = "norcar", path = "data/ebirdst/", metric = "median")

bird_annual_abundance <- list(
  yerwar = yerwar_data,
  easpho = easpho_data,
  herthr = herthr_data,
  bcnher = bcnher_data,
  pinwar = pinwar_data,
  grnher = grnher_data,
  woothr = woothr_data,
  comyel = comyel_data,
  whevir = whevir_data,
  reevir1 = reevir1_data,
  amerob = amerob_data, 
  norcar = norcar_data
)

# Define week-to-month mapping (52 weeks -> 12 months)
week_to_month <- c(
  rep(1, 4),   # Jan: weeks 1-4
  rep(2, 4),   # Feb: weeks 5-8
  rep(3, 5),   # Mar: weeks 9-13
  rep(4, 4),   # Apr: weeks 14-17
  rep(5, 4),   # May: weeks 18-21
  rep(6, 5),   # Jun: weeks 22-26
  rep(7, 4),   # Jul: weeks 27-30
  rep(8, 4),   # Aug: weeks 31-34
  rep(9, 5),   # Sep: weeks 35-39
  rep(10, 4),  # Oct: weeks 40-43
  rep(11, 4),  # Nov: weeks 44-47
  rep(12, 5)   # Dec: weeks 48-52
)

# Get Florida extent for cropping
# Read EEEV monthly rasters back from saved files
eeev_monthly_rasters <- list()

for(i in 1:12) {
  filename <- paste0("output/rasters/eeev_mean_month_", sprintf("%02d", i), "_utm.tif")
  eeev_monthly_rasters[[i]] <- rast(filename)
}

names(eeev_monthly_rasters) <- paste0("month_", sprintf("%02d", 1:12))


eeev_temp <- project(eeev_monthly_rasters[[1]], "EPSG:8857")
fl_extent <- ext(eeev_temp)

# Initialize storage
bird_monthly_data <- list()

# Process yerwar (Yellow-rumped Warbler)
yerwar_cropped <- crop(bird_annual_abundance[["yerwar"]], fl_extent)
bird_monthly_data[["yerwar"]] <- list()
for(month in 1:12) {
  week_indices <- which(week_to_month == month)
  bird_monthly_data[["yerwar"]][[month]] <- mean(yerwar_cropped[[week_indices]], na.rm = TRUE)
}

# Process easpho (Eastern Phoebe)
easpho_cropped <- crop(bird_annual_abundance[["easpho"]], fl_extent)
bird_monthly_data[["easpho"]] <- list()
for(month in 1:12) {
  week_indices <- which(week_to_month == month)
  bird_monthly_data[["easpho"]][[month]] <- mean(easpho_cropped[[week_indices]], na.rm = TRUE)
}

# Process herthr (Hermit Thrush)
herthr_cropped <- crop(bird_annual_abundance[["herthr"]], fl_extent)
bird_monthly_data[["herthr"]] <- list()
for(month in 1:12) {
  week_indices <- which(week_to_month == month)
  bird_monthly_data[["herthr"]][[month]] <- mean(herthr_cropped[[week_indices]], na.rm = TRUE)
}

# Process bcnher (Black-crowned Night Heron)
bcnher_cropped <- crop(bird_annual_abundance[["bcnher"]], fl_extent)
bird_monthly_data[["bcnher"]] <- list()
for(month in 1:12) {
  week_indices <- which(week_to_month == month)
  bird_monthly_data[["bcnher"]][[month]] <- mean(bcnher_cropped[[week_indices]], na.rm = TRUE)
}

# Process pinwar (Pine Warbler)
pinwar_cropped <- crop(bird_annual_abundance[["pinwar"]], fl_extent)
bird_monthly_data[["pinwar"]] <- list()
for(month in 1:12) {
  week_indices <- which(week_to_month == month)
  bird_monthly_data[["pinwar"]][[month]] <- mean(pinwar_cropped[[week_indices]], na.rm = TRUE)
}

# Process grnher (Green Heron)
grnher_cropped <- crop(bird_annual_abundance[["grnher"]], fl_extent)
bird_monthly_data[["grnher"]] <- list()
for(month in 1:12) {
  week_indices <- which(week_to_month == month)
  bird_monthly_data[["grnher"]][[month]] <- mean(grnher_cropped[[week_indices]], na.rm = TRUE)
}

# Process woothr (Wood Thrush)
woothr_cropped <- crop(bird_annual_abundance[["woothr"]], fl_extent)
bird_monthly_data[["woothr"]] <- list()
for(month in 1:12) {
  week_indices <- which(week_to_month == month)
  bird_monthly_data[["woothr"]][[month]] <- mean(woothr_cropped[[week_indices]], na.rm = TRUE)
}

# Process comyel (Common Yellowthroat)
comyel_cropped <- crop(bird_annual_abundance[["comyel"]], fl_extent)
bird_monthly_data[["comyel"]] <- list()
for(month in 1:12) {
  week_indices <- which(week_to_month == month)
  bird_monthly_data[["comyel"]][[month]] <- mean(comyel_cropped[[week_indices]], na.rm = TRUE)
}

# Process whevir (White-eyed Vireo)
whevir_cropped <- crop(bird_annual_abundance[["whevir"]], fl_extent)
bird_monthly_data[["whevir"]] <- list()
for(month in 1:12) {
  week_indices <- which(week_to_month == month)
  bird_monthly_data[["whevir"]][[month]] <- mean(whevir_cropped[[week_indices]], na.rm = TRUE)
}

# Process reevir1 (Red-eyed Vireo)
reevir1_cropped <- crop(bird_annual_abundance[["reevir1"]], fl_extent)
bird_monthly_data[["reevir1"]] <- list()
for(month in 1:12) {
  week_indices <- which(week_to_month == month)
  bird_monthly_data[["reevir1"]][[month]] <- mean(reevir1_cropped[[week_indices]], na.rm = TRUE)
}

# Process amerob (American Robin)
amerob_cropped <- crop(bird_annual_abundance[["amerob"]], fl_extent)
bird_monthly_data[["amerob"]] <- list()
for(month in 1:12) {
  week_indices <- which(week_to_month == month)
  bird_monthly_data[["amerob"]][[month]] <- mean(amerob_cropped[[week_indices]], na.rm = TRUE)
}

# Process norcar (Northern Cardinal)
norcar_cropped <- crop(bird_annual_abundance[["norcar"]], fl_extent)
bird_monthly_data[["norcar"]] <- list()
for(month in 1:12) {
  week_indices <- which(week_to_month == month)
  bird_monthly_data[["norcar"]][[month]] <- mean(norcar_cropped[[week_indices]], na.rm = TRUE)
}

# Save processed monthly bird data
saveRDS(bird_monthly_data, "data/ebirdst/bird_monthly_abundance.rds")

# Read in bird monthly data
#bird_monthly_data <- read_rds("data/ebirdst/bird_monthly_abundance.rds")


# Reproject all bird data to UTM 32617
for(species in names(bird_monthly_data)) {
  for(month in 1:12) {
    bird_monthly_data[[species]][[month]] <- project(bird_monthly_data[[species]][[month]], "EPSG:32617")
  }
}

# Resampling EEEV rasters to match resolution of bird data
# Use the first bird raster as template for resampling
template_raster <- bird_monthly_data[["yerwar"]][[1]]

# Resample all EEEV monthly rasters to match bird resolution
eeev_monthly_resampled <- list()

for(i in 1:12) {
  eeev_monthly_resampled[[i]] <- resample(eeev_monthly_rasters[[i]], 
                                          template_raster, 
                                          method = "bilinear")
  cat("Resampled EEEV month", i, "\n")
}

names(eeev_monthly_resampled) <- paste0("month_", sprintf("%02d", 1:12))

# Check the result
eeev_monthly_resampled[[1]]
bird_monthly_data[["yerwar"]][[1]]

# Verify everything matches
identical(ext(eeev_monthly_resampled[[1]]), ext(bird_monthly_data[["yerwar"]][[1]]))
identical(res(eeev_monthly_resampled[[1]]), res(bird_monthly_data[["yerwar"]][[1]]))

# Pixel-wise correlations

calculate_pixel_correlations <- function(eeev_rasters, bird_monthly_data, species_code) {
  
  results <- expand_grid(
    species = species_code,
    eeev_month = 1:12,
    bird_lag = 0:3  # 0=contemporary, 1-3=lags
  ) %>%
    mutate(
      bird_month = eeev_month - bird_lag,
      correlation = NA_real_,
      p_value = NA_real_,
      n_pixels = NA_integer_
    ) %>%
    filter(bird_month >= 1 & bird_month <= 12)  
  
  for(i in 1:nrow(results)) {
    eeev_month <- results$eeev_month[i]
    bird_month <- results$bird_month[i]
    
    # Extract pixel values
    eeev_vals <- values(eeev_rasters[[eeev_month]], na.rm = FALSE)
    bird_vals <- values(bird_monthly_data[[species_code]][[bird_month]], na.rm = FALSE)
    
    # Remove NA pairs
    complete_pairs <- complete.cases(eeev_vals, bird_vals)
    
    if(sum(complete_pairs) > 20) {  
      cor_test <- cor.test(eeev_vals[complete_pairs], bird_vals[complete_pairs])
      results$correlation[i] <- cor_test$estimate
      results$p_value[i] <- cor_test$p.value
      results$n_pixels[i] <- sum(complete_pairs)
    }
  }
  
  return(results)
}

# Species to analyze
species_codes <- names(bird_monthly_data)

# Calculate pixel-wise correlations
pixel_correlations <- map_dfr(species_codes, ~{
  calculate_pixel_correlations(eeev_monthly_resampled, bird_monthly_data, .x)
})

# Plot pixel-wise correlations
pixel_plot <- pixel_correlations %>%
  ggplot(aes(x = eeev_month, y = correlation, color = factor(bird_lag))) +
  geom_line() + geom_point() +
  facet_wrap(~species, scales = "free_y") +
  scale_color_viridis_d(name = "Bird lag\n(months)") +
  labs(x = "EEEV Month", y = "Pixel-wise Correlation", 
       title = "EEEV vs Bird Abundance: Pixel-wise Correlations") +
  theme_minimal()

# Prepare data with migrant/resident classification
plot_data <- pixel_correlations %>%
  # Add species type classification with seasonal timing
  mutate(
    species_type = case_when(
      # Migrants: amerob, easpho, herthr, reevir1, woothr, yerwar
      species %in% c("amerob", "easpho", "herthr", "reevir1", "woothr", "yerwar") ~ "Migrants",
      # Residents: norcar, comyel, whevir, grnher, pinwar, bcnher
      species %in% c("norcar", "comyel", "whevir", "grnher", "pinwar", "bcnher") ~ "Residents",
      TRUE ~ "Other"
    ),
    # Seasonal migration timing
    migration_timing = case_when(
      species == "reevir1" ~ "Spring Migrant",
      species %in% c("amerob", "easpho", "woothr", "yerwar", "herthr") ~ "Fall Migrant",
      species %in% c("norcar", "comyel", "whevir", "grnher", "pinwar", "bcnher") ~ "Year-round Resident",
      TRUE ~ "Other"
    ),
    # Create cleaner species names
    species_clean = case_when(
      species == "reevir1" ~ "Red-eyed Vireo",
      species == "norcar" ~ "Northern Cardinal",
      species == "amerob" ~ "American Robin",
      species == "comyel" ~ "Common Yellowthroat",
      species == "whevir" ~ "White-eyed Vireo",
      species == "yerwar" ~ "Yellow-rumped Warbler",
      species == "grnher" ~ "Green Heron",
      species == "woothr" ~ "Wood Thrush",
      species == "easpho" ~ "Eastern Phoebe",
      species == "herthr" ~ "Hermit Thrush",
      species == "pinwar" ~ "Pine Warbler",
      species == "bcnher" ~ "Black-crowned Night Heron",
      TRUE ~ species
    ),
    # Filter lag data
    include_species = case_when(
      # Show all lags for key migrants
      species %in% c("reevir1", "amerob", "woothr", "yerwar", "easpho", "herthr") ~ TRUE,
      # Show only contemporary for residents
      species %in% c("norcar", "comyel", "whevir", "grnher", "pinwar", "bcnher") & bird_lag == 0 ~ TRUE,
      TRUE ~ FALSE
    ),
    # Create ordering factor for seasonal arrangement
    species_order = case_when(
      # Spring migrant
      species == "reevir1" ~ 1,
      # Residents
      species == "pinwar" ~ 2,
      species == "comyel" ~ 3,
      species == "whevir" ~ 4,
      species == "norcar" ~ 5,
      species == "bcnher" ~ 6,
      species == "grnher" ~ 7,
      # Fall migrants
      species == "amerob" ~ 8,
      species == "yerwar"~ 9,
      species == "herthr" ~ 10,
      species == "easpho" ~ 11,
      species == "woothr"~ 12,
      TRUE ~ 99
    )
  ) %>%
  filter(include_species, !is.na(correlation)) %>%
  # Create lag labels
  mutate(
    lag_label = case_when(
      bird_lag == 0 ~ "Contemporary",
      bird_lag == 1 ~ "1-month lag",
      bird_lag == 2 ~ "2-month lag",
      bird_lag == 3 ~ "3-month lag"
    ),
    # Add migration timing to species name
    species_display = case_when(
      species == "reevir1" ~ paste0(species_clean, " (SM)"),
      species %in% c("norcar", "comyel", "whevir", "grnher", "pinwar", "bcnher") ~ paste0(species_clean, " (R)"),
      species %in% c("amerob", "easpho", "woothr", "yerwar", "herthr") ~ paste0(species_clean, " (FM)"),
      TRUE ~ species_clean
    )
  )

# Create levels
species_levels <- plot_data %>%
  select(species_display, species_order) %>%
  distinct() %>%
  arrange(species_order) %>%
  pull(species_display)

# Add the ordered factor
plot_data <- plot_data %>%
  mutate(species_display_ordered = factor(species_display, levels = species_levels))

# Create plot with ordered facets
pixel_plot_clean <- plot_data %>%
  ggplot(aes(x = eeev_month, y = correlation)) +
  # Add horizontal reference lines
  geom_hline(yintercept = 0, color = "gray50", linetype = "dashed", alpha = 0.7) +
  geom_hline(yintercept = c(-0.3, 0.3), color = "gray70", linetype = "dotted", alpha = 0.5) +
  # Add correlation lines and points
  geom_line(aes(color = lag_label, linetype = species_type),
            linewidth = 0.8, alpha = 0.8) +
  geom_point(aes(color = lag_label, shape = species_type),
             size = 2.5, alpha = 0.9) +
  # Facet by ordered species
  facet_wrap(~species_display_ordered, nrow = 3, ncol = 4, strip.position = "top") +
  # Color scheme for lag effects
  scale_color_manual(
    name = "Temporal Pattern",
    values = c("Contemporary" = "#2E86AB",      
               "1-month lag" = "#A23B72",       
               "2-month lag" = "#F18F01",         
               "3-month lag" = "#C73E1D"),      
    labels = c("Contemporary" = "Contemporary (same month)",
               "1-month lag" = "1-month lag",
               "2-month lag" = "2-month lag",
               "3-month lag" = "3-month lag")
  ) +
  # Shape by migration status
  scale_shape_manual(
    name = "Migration Status",
    values = c("Migrants" = 16,     
               "Residents" = 17)
  ) +
  # Line type by migration status  
  scale_linetype_manual(
    name = "Migration Status",
    values = c("Migrants" = "solid",
               "Residents" = "dashed")
  ) +
  # Clear axis labels
  scale_x_continuous(
    name = "Month",
    breaks = c(1, 3, 5, 7, 9, 11),
    labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov")
  ) +
  scale_y_continuous(
    name = "Correlation with Predicted EEEV Seroprevalence",
    breaks = seq(-0.4, 0.6, 0.2),
    limits = c(-0.4, 0.6)
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "gray95", color = "white"),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(face = "bold", size = 11),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30"),
    plot.caption = element_text(size = 8, color = "gray50")
  ) +
  labs(caption = "SM = Spring Migrant, R = Resident, FM = Fall Migrant")

pixel_plot_clean

# Save the plot
ggsave("figs/eeev/eeev_mod_bird_correlations_lag_focus.jpeg", pixel_plot_clean, 
       width = 12, height = 12, dpi = 600)

# Create filtered correlation plot
plot_data_filt <- pixel_correlations %>%
  filter(species %in% c("amerob", "reevir1", "pinwar")) %>%  
  mutate(
    species_type = case_when(
      species %in% c("amerob", "easpho", "herthr", "reevir1", "woothr", "yerwar") ~ "Migrants",
      species %in% c("norcar", "comyel", "whevir", "grnher", "pinwar", "bcnher") ~ "Residents",
      TRUE ~ "Other"
    ),
    # Seasonal migration timing
    migration_timing = case_when(
      species == "reevir1" ~ "Spring Migrant",
      species %in% c("amerob", "easpho", "woothr", "yerwar", "herthr") ~ "Fall Migrant",
      species %in% c("norcar", "comyel", "whevir", "grnher", "pinwar", "bcnher") ~ "Year-round Resident",
      TRUE ~ "Other"
    ),
    species_clean = case_when(
      species == "reevir1" ~ "Red-eyed Vireo",
      species == "norcar" ~ "Northern Cardinal",
      species == "amerob" ~ "American Robin",
      species == "comyel" ~ "Common Yellowthroat",
      species == "whevir" ~ "White-eyed Vireo",
      species == "yerwar" ~ "Yellow-rumped Warbler",
      species == "grnher" ~ "Green Heron",
      species == "woothr" ~ "Wood Thrush",
      species == "easpho" ~ "Eastern Phoebe",
      species == "herthr" ~ "Hermit Thrush",
      species == "pinwar" ~ "Pine Warbler",
      species == "bcnher" ~ "Black-crowned Night Heron",
      TRUE ~ species
    ),
    include_species = case_when(
      species %in% c("reevir1", "amerob", "woothr", "yerwar", "easpho", "herthr") ~ TRUE,
      species %in% c("norcar", "comyel", "whevir", "grnher", "pinwar", "bcnher") & bird_lag == 0 ~ TRUE,
      TRUE ~ FALSE
    ),
    # Create ordering factor
    species_order = case_when(
      species == "reevir1" ~ 1,  
      species == "pinwar" ~ 2,     
      species == "amerob" ~ 3,   
      TRUE ~ 99
    )
  ) %>%
  filter(include_species) %>%
  mutate(
    lag_label = case_when(
      bird_lag == 0 ~ "Contemporary",
      bird_lag == 1 ~ "1-month lag",
      bird_lag == 2 ~ "2-month lag",
      bird_lag == 3 ~ "3-month lag"
    ),
    # Add migration timing to species name
    species_display = case_when(
      species == "reevir1" ~ paste0(species_clean, " (SM)"),
      species == "pinwar" ~ paste0(species_clean, " (R)"),
      species == "amerob" ~ paste0(species_clean, " (FM)"),
      TRUE ~ species_clean
    )
  )

# Create levels
species_levels_filt <- plot_data_filt %>%
  select(species_display, species_order) %>%
  distinct() %>%
  arrange(species_order) %>%
  pull(species_display)

# Add the ordered factor
plot_data_filt <- plot_data_filt %>%
  mutate(species_display_ordered = factor(species_display, levels = species_levels_filt))

# Create filtered pixel correlation plot with consistent axes
pixel_plot_filt <- plot_data_filt %>%
  ggplot(aes(x = eeev_month, y = correlation)) +
  geom_hline(yintercept = 0, color = "gray50", linetype = "dashed", alpha = 0.7) +
  geom_hline(yintercept = c(-0.3, 0.3), color = "gray70", linetype = "dotted", alpha = 0.5) +
  geom_line(aes(color = lag_label, linetype = species_type), linewidth = 0.8, alpha = 0.8) +
  geom_point(aes(color = lag_label, shape = species_type), size = 2.5, alpha = 0.9) +
  facet_wrap(~species_display_ordered, nrow = 1, ncol = 3) +  # Same as main plot structure
  scale_color_manual(
    name = "Temporal Pattern",
    values = c("Contemporary" = "#2E86AB", "1-month lag" = "#A23B72",
               "2-month lag" = "#F18F01", "3-month lag" = "#C73E1D"),
    labels = c("Contemporary" = "Contemporary (same month)",
               "1-month lag" = "1-month lag", 
               "2-month lag" = "2-month lag",
               "3-month lag" = "3-month lag")
  ) +
  scale_shape_manual(name = "Migration Status", values = c("Migrants" = 16, "Residents" = 17)) +
  scale_linetype_manual(name = "Migration Status", values = c("Migrants" = "solid", "Residents" = "dashed")) +
  scale_x_continuous(
    name = "Month", 
    breaks = c(1, 3, 5, 7, 9, 11),
    labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov")
  ) +
  scale_y_continuous(
    name = "Correlation with Predicted EEEV Seroprevalence",
    breaks = seq(-0.4, 0.6, 0.2),
    limits = c(-0.4, 0.6)  # Same limits as main plot
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "gray95", color = "white"),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(face = "bold", size = 11),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3)
  ) +
  labs(caption = "SM = Spring Migrant, R = Resident, FM = Fall Migrant")

pixel_plot_filt

ggsave("figs/eeev/eeev_bird_correlations_filtered.jpeg", pixel_plot_filt,
       width = 10, height = 6, dpi = 600)

write_rds(pixel_correlations, "output/bird_pixel_correlations.rds")

