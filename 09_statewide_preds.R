# Load required libraries
library(sdmTMB)
library(tidyverse)
library(sf)
library(terra)
library(evaluate)
library(wesanderson)
library(ggpubr)
library(cowplot)

# Read in prediction data frame
pred_df <- read_rds("data/pred_df.rds")

# Read in models
eeev_mod <- read_rds("output/eeev_model_primary_ver01.rds")
eeev_mod_alt <- read_rds("output/eeev_model_alt_ver01.rds")

# Assemble pred_df for predictions
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
  mutate(year = as.factor(year))

rm(pred_df)

# Get statewide predictions
# Primary model
eeev_mod_state_preds <- eeev_mod %>% 
  predict(newdata = pred_df_prepared, 
           se_fit = FALSE, 
           re_form_iid = NA, 
           nsim = 100)

# Save data
write_rds(eeev_mod_state_preds, "output/eeev_model_primary_state_preds.rds")

# Alternative model
eeev_mod_alt_state_preds <- eeev_mod_alt %>% 
  simulate(newdata = pred_df_prepared, 
           se_fit = FALSE, 
           re_form_iid = NA, 
           nsim = 100)

# Save data
write_rds(eeev_mod_alt_state_preds, "output/eeev_model_alt_state_preds.rds")

# Read in prediction outputs
eeev_mod_statewide_raw <- read_rds("output/eeev_model_primary_state_preds.rds")

# Define custom truncating function for plotting
truncate <- function(obj, lower, upper){
  if(missing(upper)){
    obj[obj > stats::quantile(obj, lower, na.rm = T)] <- stats::quantile(obj, lower, na.rm = T)
  }
  if(missing(lower)){
    obj[obj > stats::quantile(obj, upper, na.rm = T)] <- stats::quantile(obj, upper, na.rm = T)
  }
  else(
    obj[obj > stats::quantile(obj, c(lower, upper, na.rm = T))] <- stats::quantile(obj, c(lower, upper, na.rm = T))
  )
  return(obj)
}

# Process monthly statewide predictions
monthly_statewide <- pred_df_prepared %>%
  select(X, Y, year, month) %>%
  bind_cols(
    data.frame(
      mean = rowMeans(eeev_mod_statewide_raw),
      se = apply(eeev_mod_statewide_raw, 1, FUN = function(x){sd(x)/sqrt(100)}))) %>%
  mutate(
    # Convert from cloglog link space to response space
    mean = 1 - exp(-exp(mean)),
    mean_0999 = truncate(mean, upper = 0.999),
    mean_099 = truncate(mean, upper = 0.99),
    mean_098 = truncate(mean, upper = 0.98),
    mean_095 = truncate(mean, upper = 0.95),
    mean_090 = truncate(mean, upper = 0.90),
    # SE transformation for cloglog
    se = 1 - exp(-exp(se)),
    se_0999 = truncate(se, upper = 0.999),
    se_099 = truncate(se, upper = 0.99),
    se_098 = truncate(se, upper = 0.98),
    se_095 = truncate(se, upper = 0.95)
  ); rm(eeev_mod_statewide_raw)

# Add month names
monthly_statewide <- monthly_statewide %>%
  left_join(
    data.frame(
      month = as.character(1:12),
      name = c("January", "February", "March", "April", "May", "June", 
               "July", "August", "September", "October", "November", "December")
    )
  )

monthly_statewide$names <- factor(monthly_statewide$name, 
                                  levels = c("January", "February", "March", "April", "May", "June",
                                             "July", "August", "September", "October", "November", "December"))

# Create annual aggregated data (collapse months within years)
annual_statewide <- monthly_statewide %>%
  group_by(X, Y, year) %>%
  summarize(
    mean = mean(mean),
    mean_0999 = mean(mean_0999),
    mean_099 = mean(mean_099),
    mean_098 = mean(mean_098),
    mean_095 = mean(mean_095),
    mean_090 = mean(mean_090),
    se = mean(se),
    se_0999 = mean(se_0999),
    se_099 = mean(se_099),
    se_098 = mean(se_098),
    se_095 = mean(se_095),
    .groups = 'drop'
  )

# Create summary by location (across all years and months)
monthly_statewide_summ <- monthly_statewide %>%
  group_by(X, Y) %>%
  summarize(
    mean = mean(mean),
    sum = sum(mean),
    se = mean(se),
    mean_0999 = mean(mean_0999),
    mean_099 = mean(mean_099),
    mean_098 = mean(mean_098),
    mean_095 = mean(mean_095),
    mean_090 = mean(mean_090),
    se_0999 = mean(se_0999),
    se_099 = mean(se_099),
    se_098 = mean(se_098),
    se_095 = mean(se_095), 
    .groups = 'drop'
  ) %>%
  mutate(
    sum_0999 = truncate(sum, upper = 0.999),
    sum_099 = truncate(sum, upper = 0.99),
    sum_098 = truncate(sum, upper = 0.98),
    sum_095 = truncate(sum, upper = 0.95),
    sum_090 = truncate(sum, upper = 0.90)
  )

# Define color palette and read Florida shapefile
pal <- wes_palette("Zissou1", 256, type = "continuous")
fl <- read_rds("data/fl_polygon_crop.rds")

# Annual aggregated predictions (2005-2019)
annual_panel <- ggplot() +
  geom_raster(data = annual_statewide, aes(x = X, y = Y, fill = mean_098)) +
  geom_sf(data = fl, fill = NA, col = "black") +
  scale_fill_gradientn("EEEV\nSeroprevalence",
                       colors = pal,
                       limits = c(0, max(annual_statewide$mean_098)),
                       breaks = seq(0, max(annual_statewide$mean_098), length.out = 3),
                       labels = c("Low", "Mid", "High"),
                       na.value = NA) +
  facet_wrap(~year, nrow = 5, ncol = 3) +
  theme_void() +
  theme(
    legend.position = "right",
    #legend.position.inside = c(0.95, 0.5),
    #legend.title.position = "left",
    legend.title = element_text(size = 10, face = "bold", hjust = 0.5, vjust = 0.5),
    legend.key.height = unit(0.45, "cm"),
    legend.key.width = unit(0.4, "cm"),
    legend.ticks = element_blank(),
    legend.frame = element_rect(colour = "black", linewidth = 0.3),
    legend.text = element_text(size = 7, colour = "black", face = "bold"),
    strip.text = element_text(face = "bold", size = 10), 
    plot.margin = margin(0, -1, 0, -1, unit = "cm")
  )

annual_panel

# Monthly predictions for 2005 and 2010
monthly_panel_2005_2010 <- ggplot() +
  geom_raster(data = monthly_statewide %>%
                filter(year %in% as.character(c(2005, 2010))), 
              aes(x = X, y = Y, fill = mean_098)) +
  geom_sf(data = fl, fill = NA, col = "black") +
  scale_fill_gradientn("Monthly\nEEEV\nSeroprevalence",
                       colors = pal, na.value = NA) +
  facet_grid(rows = vars(year), cols = vars(names), switch = "both") +
  theme_void() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 8),
    strip.text.y.left = element_text(angle = 90), 
    plot.margin = margin(-1, 1, -1, 1, unit = "cm")
  )

monthly_panel_2005_2010

year_month_panel <- plot_grid(
  annual_panel, 
  monthly_panel_2005_2010, 
  nrow = 2, 
  rel_heights = c(1, 0.25)
)

year_month_panel

ggsave("figs/eeev/year_month_preds_panel.jpeg", year_month_panel, 
       height = 10, width = 8, units = "in", dpi = 600)

# All monthly predictions 2005-2019
monthly_panel_all <- ggplot() +
  geom_raster(data = monthly_statewide, aes(x = X, y = Y, fill = mean_098)) +
  geom_sf(data = fl, fill = NA, col = "black") +
  scale_fill_gradientn("Monthly\nEEEV\nSeroprevalence",
                       colors = pal,
                       limits = c(0, max(monthly_statewide$mean_098)),
                       breaks = seq(0, max(monthly_statewide$mean_098), length.out = 3),
                       labels = c("Low", "Mid", "High"),
                       na.value = NA) +
  facet_grid(rows = vars(year), cols = vars(names), switch = "both") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.key.height = unit(0.45, "cm"),
    legend.key.width = unit(0.4, "cm"),
    legend.text = element_text(size = 7, colour = "black", face = "bold"),
    legend.frame = element_rect(colour = "black", linewidth = 0.3),
    strip.text = element_text(face = "bold", size = 6),
    strip.text.y.left = element_text(angle = 90), 
    plot.margin = margin(-1, 0, -1, 0, unit = "cm")
  )

monthly_panel_all

# Save
ggsave("figs/eeev/monthly_all_years.png", monthly_panel_all, 
       width = 8, height = 11, units = "in", dpi = 600)

# Mean and Uncertainty plots
base <- ggplot() +
  theme_void() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.3, 0.45),
    legend.key.height = unit(0.55, "cm"),
    legend.key.width = unit(0.4, "cm"),
    legend.ticks = element_blank(),
    legend.frame = element_rect(colour = "black", linewidth = 0.3),
    legend.text = element_text(size = 7, colour = "black", face = "bold")
  )

monthly_uncertainty <- base +
  geom_raster(data = monthly_statewide_summ, aes(x = X, y = Y, fill = se)) +
  geom_sf(data = fl, fill = NA, col = "black") +
  scale_fill_gradientn("Uncertainty",
                       colors = pal,
                       limits = c(min(monthly_statewide_summ$se), max(monthly_statewide_summ$se)),
                       breaks = seq(min(monthly_statewide_summ$se), max(monthly_statewide_summ$se), length.out = 5),
                       labels = seq(min(monthly_statewide_summ$se), max(monthly_statewide_summ$se), length.out = 5) %>% round(3),
                       na.value = NA)

monthly_mean <- base +
  geom_raster(data = monthly_statewide_summ, aes(x = X, y = Y, fill = mean_090)) +
  geom_sf(data = fl, fill = NA, col = "black") +
  scale_fill_gradientn("Mean EEEV Seroprevalence",
                       colors = pal,
                       limits = c(0, max(monthly_statewide_summ$mean_090)),
                       breaks = seq(0, max(monthly_statewide_summ$mean_090), length.out = 5),
                       labels = seq(0, max(monthly_statewide_summ$mean_090), length.out = 5) %>% round(5),
                       na.value = NA)

monthly_mean

monthly_uncertainty

# Combine
combined_plot <- ggarrange(
  monthly_mean,
  monthly_uncertainty,
  nrow = 1, ncol = 2, align = "hv"
)

mean_uncertainty <- plot_grid(
  monthly_mean + 
    ggtitle("Prediction Mean") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12)),
  monthly_uncertainty + 
    ggtitle("Uncertainty") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12)),
  nrow = 1, ncol = 2, align = "hv"
)

mean_uncertainty

ggsave("figs/eeev/mean_uncertainty.png", mean_uncertainty, width = 8, height = 4, dpi = 600)

## Repeat for alternate model:

# Read in predictions
eeev_mod_alt_statewide_raw <- read_rds("output/eeev_model_alt_state_preds.rds")

# Process monthly statewide predictions
monthly_statewide_alt <- pred_df_prepared %>%
  select(X, Y, year, month) %>%
  bind_cols(
    data.frame(
      mean = rowMeans(eeev_mod_alt_statewide_raw),
      se = apply(eeev_mod_alt_statewide_raw, 1, FUN = function(x){sd(x)/sqrt(100)}))) %>%
  mutate(
    # Convert from cloglog link space to response space
    mean = 1 - exp(-exp(mean)),
    mean_0999 = truncate(mean, upper = 0.999),
    mean_099 = truncate(mean, upper = 0.99),
    mean_098 = truncate(mean, upper = 0.98),
    mean_095 = truncate(mean, upper = 0.95),
    mean_090 = truncate(mean, upper = 0.90),
    # SE transformation for cloglog
    se = 1 - exp(-exp(se)),
    se_0999 = truncate(se, upper = 0.999),
    se_099 = truncate(se, upper = 0.99),
    se_098 = truncate(se, upper = 0.98),
    se_095 = truncate(se, upper = 0.95)
  ); rm(eeev_mod_alt_statewide_raw)

# Add month names
monthly_statewide_alt <- monthly_statewide_alt %>%
  left_join(
    data.frame(
      month = as.character(1:12),
      name = c("January", "February", "March", "April", "May", "June",
               "July", "August", "September", "October", "November", "December")
    )
  )

monthly_statewide_alt$names <- factor(monthly_statewide_alt$name,
                                      levels = c("January", "February", "March", "April", "May", "June",
                                                 "July", "August", "September", "October", "November", "December"))

# Create annual aggregated data (collapse months within years)
annual_statewide_alt <- monthly_statewide_alt %>%
  group_by(X, Y, year) %>%
  summarize(
    mean = mean(mean),
    mean_0999 = mean(mean_0999),
    mean_099 = mean(mean_099),
    mean_098 = mean(mean_098),
    mean_095 = mean(mean_095),
    mean_090 = mean(mean_090),
    se = mean(se),
    se_0999 = mean(se_0999),
    se_099 = mean(se_099),
    se_098 = mean(se_098),
    se_095 = mean(se_095),
    .groups = 'drop'
  )

# Create summary by location (across all years and months)
monthly_statewide_summ_alt <- monthly_statewide_alt %>%
  group_by(X, Y) %>%
  summarize(
    mean = mean(mean),
    sum = sum(mean),
    se = mean(se),
    mean_0999 = mean(mean_0999),
    mean_099 = mean(mean_099),
    mean_098 = mean(mean_098),
    mean_095 = mean(mean_095),
    mean_090 = mean(mean_090),
    se_0999 = mean(se_0999),
    se_099 = mean(se_099),
    se_098 = mean(se_098),
    se_095 = mean(se_095),
    .groups = 'drop'
  ) %>%
  mutate(
    sum_0999 = truncate(sum, upper = 0.999),
    sum_099 = truncate(sum, upper = 0.99),
    sum_098 = truncate(sum, upper = 0.98),
    sum_095 = truncate(sum, upper = 0.95),
    sum_090 = truncate(sum, upper = 0.90)
  )

# Annual aggregated predictions (2005-2019)
annual_panel_alt <- ggplot() +
  geom_raster(data = annual_statewide_alt, aes(x = X, y = Y, fill = mean_098)) +
  geom_sf(data = fl, fill = NA, col = "black") +
  scale_fill_gradientn("EEEV\nSeroprevalence",
                       colors = pal,
                       limits = c(0, max(annual_statewide_alt$mean_098)),
                       breaks = seq(0, max(annual_statewide_alt$mean_098), length.out = 3),
                       labels = c("Low", "Mid", "High"),
                       na.value = NA) +
  facet_wrap(~year, nrow = 5, ncol = 3) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold", hjust = 0.5, vjust = 0.5),
    legend.key.height = unit(0.45, "cm"),
    legend.key.width = unit(0.4, "cm"),
    legend.ticks = element_blank(),
    legend.frame = element_rect(colour = "black", linewidth = 0.3),
    legend.text = element_text(size = 7, colour = "black", face = "bold"),
    strip.text = element_text(face = "bold", size = 10),
    plot.margin = margin(0, -1, 0, -1, unit = "cm")
  )

annual_panel_alt

# Monthly predictions for 2005 and 2010
monthly_panel_2005_2010_alt <- ggplot() +
  geom_raster(data = monthly_statewide_alt %>%
                filter(year %in% as.character(c(2005, 2010))),
              aes(x = X, y = Y, fill = mean_098)) +
  geom_sf(data = fl, fill = NA, col = "black") +
  scale_fill_gradientn("Monthly\nEEEV\nSeroprevalence",
                       colors = pal, na.value = NA) +
  facet_grid(rows = vars(year), cols = vars(names), switch = "both") +
  theme_void() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 8),
    strip.text.y.left = element_text(angle = 90),
    plot.margin = margin(-1, 1, -1, 1, unit = "cm")
  )

year_month_panel_alt <- plot_grid(
  annual_panel_alt,
  monthly_panel_2005_2010_alt,
  nrow = 2,
  rel_heights = c(1, 0.25)
)

year_month_panel_alt

ggsave("figs/eeev/year_month_preds_panel_alt.jpeg", year_month_panel_alt,
       height = 10, width = 8, units = "in", dpi = 600)

# All monthly predictions 2005-2019 
monthly_panel_all_alt <- ggplot() +
  geom_raster(data = monthly_statewide_alt, aes(x = X, y = Y, fill = mean_098)) +
  geom_sf(data = fl, fill = NA, col = "black") +
  scale_fill_gradientn("Monthly\nEEEV\nSeroprevalence",
                       colors = pal,
                       limits = c(0, max(monthly_statewide_alt$mean_098)),
                       breaks = seq(0, max(monthly_statewide_alt$mean_098), length.out = 3),
                       labels = c("Low", "Mid", "High"),
                       na.value = NA) +
  facet_grid(rows = vars(year), cols = vars(names), switch = "both") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.key.height = unit(0.45, "cm"),
    legend.key.width = unit(0.4, "cm"),
    legend.text = element_text(size = 7, colour = "black", face = "bold"),
    legend.frame = element_rect(colour = "black", linewidth = 0.3),
    strip.text = element_text(face = "bold", size = 6),
    strip.text.y.left = element_text(angle = 90),
    plot.margin = margin(-1, 0, -1, 0, unit = "cm")
  )

ggsave("figs/eeev/monthly_all_years_alt.png", monthly_panel_all_alt,
       width = 8, height = 11, units = "in", dpi = 600)

# Mean and Uncertainty plots
monthly_uncertainty_alt <- base +
  geom_raster(data = monthly_statewide_summ_alt, aes(x = X, y = Y, fill = se)) +
  geom_sf(data = fl, fill = NA, col = "black") +
  scale_fill_gradientn("Uncertainty",
                       colors = pal,
                       limits = c(min(monthly_statewide_summ_alt$se), max(monthly_statewide_summ_alt$se)),
                       breaks = seq(min(monthly_statewide_summ_alt$se), max(monthly_statewide_summ_alt$se), length.out = 5),
                       labels = seq(min(monthly_statewide_summ_alt$se), max(monthly_statewide_summ_alt$se), length.out = 5) %>% round(3),
                       na.value = NA)

monthly_mean_alt <- base +
  geom_raster(data = monthly_statewide_summ_alt, aes(x = X, y = Y, fill = mean_090)) +
  geom_sf(data = fl, fill = NA, col = "black") +
  scale_fill_gradientn("Mean EEEV Seroprevalence",
                       colors = pal,
                       limits = c(0, max(monthly_statewide_summ_alt$mean_090)),
                       breaks = seq(0, max(monthly_statewide_summ_alt$mean_090), length.out = 5),
                       labels = seq(0, max(monthly_statewide_summ_alt$mean_090), length.out = 5) %>% round(5),
                       na.value = NA)

# Combine
mean_uncertainty_alt <- plot_grid(
  monthly_mean_alt +
    ggtitle("Prediction Mean") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12)),
  monthly_uncertainty_alt +
    ggtitle("Uncertainty") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12)),
  nrow = 1, ncol = 2, align = "hv"
)

ggsave("figs/eeev/mean_uncertainty_alt.png", mean_uncertainty_alt, width = 8, height = 4, dpi = 600)