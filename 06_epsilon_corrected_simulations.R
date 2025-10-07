# Load required libraries
library(tidyverse)
library(sdmTMB)
library(sf)
library(ggpubr)

# Read in data
dat_train <- read_rds("data/eeev_env_covs_train_scaled.rds")
dat_test <- read_rds("data/eeev_env_covs_test_scaled.rds")

# Read in model objects
eeev_mod <- read_rds("output/eeev_model_primary_ver01.rds")
eeev_mod_alt <- read_rds("output/eeev_model_alt_ver01.rds")

# Combine data for simulation
dat <- bind_rows(dat_train, dat_test) %>% 
  arrange(year, month, site_id) %>% 
  st_drop_geometry()

# Epsilon corrected simulations
# Primary model
eeev_mod_sim <- simulate(eeev_mod, 
                         newdata = dat, 
                         size = dat$weights,
                         type = "mle-mvn", 
                         mle_mvn_samples = "multiple", 
                         nsim = 500, 
                         observation_error = FALSE)

# Alternate model
eeev_mod_alt_sim <- simulate(eeev_mod_alt, 
                             newdata = dat, 
                             size = dat$weights,
                             type = "mle-mvn", 
                             mle_mvn_samples = "multiple", 
                             nsim = 500, 
                             observation_error = FALSE)

# Save simulations
write_rds(eeev_mod_sim, "output/eeev_model_primary_eps_correct_sim_ver01.rds")
write_rds(eeev_mod_alt_sim, "output/eeev_model_alt_eps_correct_sim_ver01.rds")

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

# Change names
primary_sims <- eeev_mod_sim
alt_sims <- eeev_mod_alt_sim

# Set attributes for get_index_sims()
attr(primary_sims, "time") <- "time_step"
attr(primary_sims, "link") <- "response"
row.names(primary_sims) <- 1:nrow(dat)

attr(alt_sims, "time") <- "time_step"
attr(alt_sims, "link") <- "response"
row.names(alt_sims) <- 1:nrow(dat)

# Test with one confidence level to see dimensions
test_result <- get_index_sims(primary_sims, 
                              level = 0.95,
                              agg_function = function(x) sum(x),
                              area_function = function(x, area) x * area)

#check the dimensions
nrow(dat)
dim(primary_sims) 
dim(test_result)
colnames(test_result)
head(test_result)
length(unique(dat$time_step))
range(dat$time_step)
head(table(dat$time_step))

# Trying with aggregation
# Convert simulations to data frame with time_step
sim_df <- as.data.frame(primary_sims) %>%
  mutate(time_step = dat$time_step) %>%
  # Aggregate across sites within each time_step (sum the proportions)
  group_by(time_step) %>%
  summarise(across(starts_with("V"), sum, na.rm = TRUE), .groups = "drop")

# Convert back to matrix format
aggregated_sims <- as.matrix(sim_df[, -1])  
row.names(aggregated_sims) <- sim_df$time_step

# Set attributes
attr(aggregated_sims, "time") <- "time_step"
attr(aggregated_sims, "link") <- "response"

# Now test get_index_sims on the aggregated data
test_aggregated <- get_index_sims(aggregated_sims, 
                                  level = 0.95,
                                  agg_function = function(x) sum(x),
                                  area_function = function(x, area) x * area)

#check results
dim(aggregated_sims)
dim(test_aggregated)
head(test_aggregated)

# Convert simulations from counts back to proportions
eeev_mod_sim_props <- eeev_mod_sim / dat$weights
eeev_mod_alt_sim_props <- eeev_mod_alt_sim / dat$weights

# Handle division by zero (where weights = 0)
eeev_mod_sim_props[dat$weights == 0] <- 0
eeev_mod_alt_sim_props[dat$weights == 0] <- 0

# Create time lookup 
time_lookup <- dat %>%
  select(time_step, date) %>%
  distinct() %>%
  mutate(
    year = year(date),    
    month = month(date)   
  ) %>%
  arrange(time_step)

# Check the time lookup
head(time_lookup, 10)
range(time_lookup$date)

# Primary model - weighted average
primary_sim_df <- as.data.frame(eeev_mod_sim_props) %>%
  mutate(time_step = dat$time_step,
         weights = dat$weights) %>%
  group_by(time_step) %>%
  summarise(
    across(starts_with("V"), \(x) weighted.mean(x, w = weights, na.rm = TRUE)),
    .groups = "drop"
  )

# Aggregate
primary_sims_agg <- as.matrix(primary_sim_df[, -1])
row.names(primary_sims_agg) <- primary_sim_df$time_step
attr(primary_sims_agg, "time") <- "time_step"
attr(primary_sims_agg, "link") <- "response"

# Alternative model - weighted average 
alt_sim_df <- as.data.frame(eeev_mod_alt_sim_props) %>%
  mutate(time_step = dat$time_step,
         weights = dat$weights) %>%
  group_by(time_step) %>%
  summarise(
    across(starts_with("V"), \(x) weighted.mean(x, w = weights, na.rm = TRUE)),
    .groups = "drop"
  )

# Aggregate
alt_sims_agg <- as.matrix(alt_sim_df[, -1])
row.names(alt_sims_agg) <- alt_sim_df$time_step
attr(alt_sims_agg, "time") <- "time_step"
attr(alt_sims_agg, "link") <- "response"

# Check aggregation worked
dim(primary_sims_agg)
dim(alt_sims_agg)

# Test one result to check scale

test_result <- suppressWarnings(get_index_sims(primary_sims_agg, 
                                               level = 0.95,
                                               agg_function = function(x) mean(x),
                                               area_function = function(x, area) x * area))

head(test_result, 10)

round(min(test_result$est), 6)
round(max(test_result$est), 6)


# Get index simulations with confidence intervals

# Primary model
eeev_primary_index <- bind_rows(
  suppressWarnings(get_index_sims(primary_sims_agg, 
                                  level = 0.95,
                                  agg_function = function(x) mean(x),
                                  area_function = function(x, area) x * area)) %>%
    mutate(ci = "95", upr = truncate(upr, upper = 0.95)),
  suppressWarnings(get_index_sims(primary_sims_agg, 
                                  level = 0.80,
                                  agg_function = function(x) mean(x),
                                  area_function = function(x, area) x * area)) %>%
    mutate(ci = "80", upr = truncate(upr, upper = 0.95)),
  suppressWarnings(get_index_sims(primary_sims_agg, 
                                  level = 0.50,
                                  agg_function = function(x) mean(x),
                                  area_function = function(x, area) x * area)) %>%
    mutate(ci = "50", upr = truncate(upr, upper = 0.95))
) %>%
  mutate(model = "Primary")

# Alternative model  
eeev_alt_index <- bind_rows(
  suppressWarnings(get_index_sims(alt_sims_agg, 
                                  level = 0.95,
                                  agg_function = function(x) mean(x),
                                  area_function = function(x, area) x * area)) %>%
    mutate(ci = "95", upr = truncate(upr, upper = 0.95)),
  suppressWarnings(get_index_sims(alt_sims_agg, 
                                  level = 0.80,
                                  agg_function = function(x) mean(x),
                                  area_function = function(x, area) x * area)) %>%
    mutate(ci = "80", upr = truncate(upr, upper = 0.95)),
  suppressWarnings(get_index_sims(alt_sims_agg, 
                                  level = 0.50,
                                  agg_function = function(x) mean(x),
                                  area_function = function(x, area) x * area)) %>%
    mutate(ci = "50", upr = truncate(upr, upper = 0.95))
) %>%
  mutate(model = "Alternative")


# Join with index results using the fixed time_lookup
all_sims <- bind_rows(eeev_primary_index, eeev_alt_index) %>%
  left_join(time_lookup, by = "time_step") %>%
  mutate(data_type = case_when(
    year %in% 2005:2017 ~ "Training",
    year %in% 2018:2019 ~ "Testing",
    TRUE ~ "Other"
  ))

# Check the joined results
head(all_sims %>% select(time_step, est, date, year, month, model, data_type), 10)


# Prepare observed data
# Training data - weighted average
eeev_observed_train <- dat_train %>%
  st_drop_geometry() %>%
  group_by(date) %>%
  summarise(eeev_prop = weighted.mean(eeev_prop, w = weights, na.rm = TRUE), .groups = "drop") %>%
  mutate(data_type = "Training")

# Testing data - weighted average
eeev_observed_test <- dat_test %>%
  st_drop_geometry() %>%
  group_by(date) %>%
  summarise(eeev_prop = weighted.mean(eeev_prop, w = weights, na.rm = TRUE), .groups = "drop") %>%
  mutate(data_type = "Testing")

# Combine observed data
eeev_observed <- bind_rows(eeev_observed_train, eeev_observed_test)

# Check observed data
head(eeev_observed, 10)

round(min(eeev_observed$eeev_prop), 6)
round(max(eeev_observed$eeev_prop), 6)

# Compare scales
round(min(all_sims$est), 6)
round(max(all_sims$est), 6)
round(min(eeev_observed$eeev_prop), 6)
round(max(eeev_observed$eeev_prop), 6)

# Create plot
temporal_preds_plot <- ggplot() +
  geom_area(data = eeev_observed,
            aes(x = date, y = eeev_prop, fill = paste(data_type, "Observed")),
            col = "grey20", alpha = 0.7) +
  ggdist::geom_lineribbon(data = all_sims,
                          aes(x = date, y = est, ymin = lwr, ymax = upr,
                              fill = "Model Prediction", alpha = ci),
                          col = "#2d5016", linewidth = 1) +
  geom_vline(xintercept = as.Date("2018-01-01"),
             linetype = "dashed", color = "grey40", linewidth = 1) +
  annotate("text", x = as.Date("2017-11-01"), y = max(eeev_observed$eeev_prop, na.rm = TRUE),
           label = "Prediction horizon", color = "grey40", fontface = "bold", angle = 90) +
  facet_wrap(~factor(model, levels = c("Primary", "Alternative")), ncol = 1) +
  scale_fill_manual("Estimate",
                    values = c("Training Observed" = "#F69C73FF",
                               "Testing Observed" = "#E31A1C",
                               "Model Prediction" = "#2d5016"),
                    breaks = c("Training Observed", "Testing Observed", "Model Prediction")) +  
  scale_alpha_manual("Prediction CI (%)",
                     values = c("95" = 0.15, "80" = 0.25, "50" = 0.30)) +
  guides(fill = guide_legend(order = 1),        
         alpha = guide_legend(order = 2)) +
  scale_x_date(expand = expansion(mult = c(0.02, 0.02))) +  
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.05))) +  
  labs(y = "EEEV Seropositive Proportions", x = "Date") +
  theme_pubclean() +
  theme(legend.position = "bottom",
        legend.background = element_blank(),
        legend.title = element_text(size = 10, face = "bold", color = "grey10"),
        legend.text = element_text(size = 10, color = "grey10"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.text = element_text(color = "black", face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_blank())  

temporal_preds_plot

#save plot
ggsave("figs/eeev/temporal_predictions_eps_corrected_ver01.png", temporal_preds_plot,
       height = 8, width = 9, units = "in", dpi = 600)

# Create individual plots

model_blue <- "#1f77b4"  

# Create Primary Model Plot
primary_plot <- ggplot() +
  geom_area(data = eeev_observed,
            aes(x = date, y = eeev_prop, fill = paste(data_type, "Observed")),
            col = "grey20", alpha = 0.7) +
  ggdist::geom_lineribbon(data = filter(all_sims, model == "Primary"),
                          aes(x = date, y = est, ymin = lwr, ymax = upr,
                              fill = "Model Prediction", alpha = ci),
                          col = model_blue, linewidth = 1) +
  geom_vline(xintercept = as.Date("2018-01-01"),
             linetype = "dashed", color = "grey40", linewidth = 1) +
  annotate("text", x = as.Date("2017-11-01"), y = max(eeev_observed$eeev_prop, na.rm = TRUE),
           label = "Prediction horizon", color = "grey40", fontface = "bold", angle = 90) +
  scale_fill_manual("Estimate",
                    values = c("Training Observed" = "#F69C73FF",
                               "Testing Observed" = "#E31A1C",
                               "Model Prediction" = model_blue),
                    breaks = c("Training Observed", "Testing Observed", "Model Prediction")) +
  scale_alpha_manual("Prediction CI (%)",
                     values = c("95" = 0.15, "80" = 0.25, "50" = 0.30)) +
  guides(fill = guide_legend(order = 1),        
         alpha = guide_legend(order = 2)) +
  scale_x_date(expand = expansion(mult = c(0.02, 0.02))) +  
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.05))) +  
  labs(y = "EEEV Seroprevalence", x = "Date") +
  theme_pubclean() +
  theme(legend.position = "bottom",
        legend.background = element_blank(),
        legend.title = element_text(size = 10, face = "bold", color = "grey10"),
        legend.text = element_text(size = 10, color = "grey10"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.text = element_text(color = "black", face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 12, face = "bold", hjust = 0.5))

# Create Alternative Model Plot
alternative_plot <- ggplot() +
  geom_area(data = eeev_observed,
            aes(x = date, y = eeev_prop, fill = paste(data_type, "Observed")),
            col = "grey20", alpha = 0.7) +
  ggdist::geom_lineribbon(data = filter(all_sims, model == "Alternative"),
                          aes(x = date, y = est, ymin = lwr, ymax = upr,
                              fill = "Model Prediction", alpha = ci),
                          col = model_blue, linewidth = 1) +
  geom_vline(xintercept = as.Date("2018-01-01"),
             linetype = "dashed", color = "grey40", linewidth = 1) +
  annotate("text", x = as.Date("2017-11-01"), y = max(eeev_observed$eeev_prop, na.rm = TRUE),
           label = "Prediction horizon", color = "grey40", fontface = "bold", angle = 90) +
  scale_fill_manual("Estimate",
                    values = c("Training Observed" = "#F69C73FF",
                               "Testing Observed" = "#E31A1C",
                               "Model Prediction" = model_blue),
                    breaks = c("Training Observed", "Testing Observed", "Model Prediction")) +
  scale_alpha_manual("Prediction CI (%)",
                     values = c("95" = 0.15, "80" = 0.25, "50" = 0.30)) +
  guides(fill = guide_legend(order = 1),        
         alpha = guide_legend(order = 2)) +
  scale_x_date(expand = expansion(mult = c(0.02, 0.02))) +  
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.05))) +  
  labs(y = "EEEV Seroprevalence", x = "Date") +
  theme_pubclean() +
  theme(legend.position = "bottom",
        legend.background = element_blank(),
        legend.title = element_text(size = 10, face = "bold", color = "grey10"),
        legend.text = element_text(size = 10, color = "grey10"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.text = element_text(color = "black", face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 12, face = "bold", hjust = 0.5))

# Display the plots
primary_plot
alternative_plot

# Save the plots
ggsave("figs/eeev/primary_model_temporal_predictions.jpeg", primary_plot,
       width = 10, height = 6, dpi = 600)

ggsave("figs/eeev/alternative_model_temporal_predictions.jpeg", alternative_plot,
       width = 10, height = 6, dpi = 600)

#save workspace image
save.image("07_epsilon_corrected_simulations.RData")
