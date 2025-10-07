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

# Combine data for predictions
dat <- bind_rows(dat_train, dat_test) %>% 
  arrange(year, month, site_id) %>% 
  st_drop_geometry()


# Create new data for response curves
n <- 10 #number of points along response curve

# Primary model - significant variables:
# prcp_lag12 (significant effect)
eeev_mod_newdata_prcp_lag12 <- eeev_mod %>%
  predict(se_fit = TRUE, re_form = NA, re_form_iid = NA, type = "response",
          newdata = data.frame(
            prcp_lag1 = rep(mean(dat$prcp_lag1, na.rm = TRUE), n),
            prcp_lag5 = rep(mean(dat$prcp_lag5, na.rm = TRUE), n),
            prcp_lag12 = seq(min(dat$prcp_lag12, na.rm = TRUE),
                             max(dat$prcp_lag12, na.rm = TRUE), length.out = n),
            tmax_lag6 = rep(mean(dat$tmax_lag6, na.rm = TRUE), n),
            tmax_lag12 = rep(mean(dat$tmax_lag12, na.rm = TRUE), n),
            tmin_lag1 = rep(mean(dat$tmin_lag1, na.rm = TRUE), n),
            forest = rep(mean(dat$forest, na.rm = TRUE), n),
            wetlands = rep(mean(dat$wetlands, na.rm = TRUE), n),
            county = NA, site_id = NA, time_step = 180)
  )

# tmax_lag6 (marginally significant)
eeev_mod_newdata_tmax_lag6 <- eeev_mod %>%
  predict(se_fit = TRUE, re_form = NA, re_form_iid = NA, type = "response",
          newdata = data.frame(
            prcp_lag1 = rep(mean(dat$prcp_lag1, na.rm = TRUE), n),
            prcp_lag5 = rep(mean(dat$prcp_lag5, na.rm = TRUE), n),
            prcp_lag12 = rep(mean(dat$prcp_lag12, na.rm = TRUE), n),
            tmax_lag6 = seq(min(dat$tmax_lag6, na.rm = TRUE),
                            max(dat$tmax_lag6, na.rm = TRUE), length.out = n),
            tmax_lag12 = rep(mean(dat$tmax_lag12, na.rm = TRUE), n),
            tmin_lag1 = rep(mean(dat$tmin_lag1, na.rm = TRUE), n),
            forest = rep(mean(dat$forest, na.rm = TRUE), n),
            wetlands = rep(mean(dat$wetlands, na.rm = TRUE), n),
            county = NA, site_id = NA, time_step = 180)
  )

# tmin_lag1 (significant effect)
eeev_mod_newdata_tmin_lag1 <- eeev_mod %>%
  predict(se_fit = TRUE, re_form = NA, re_form_iid = NA, type = "response",
          newdata = data.frame(
            prcp_lag1 = rep(mean(dat$prcp_lag1, na.rm = TRUE), n),
            prcp_lag5 = rep(mean(dat$prcp_lag5, na.rm = TRUE), n),
            prcp_lag12 = rep(mean(dat$prcp_lag12, na.rm = TRUE), n),
            tmax_lag6 = rep(mean(dat$tmax_lag6, na.rm = TRUE), n),
            tmax_lag12 = rep(mean(dat$tmax_lag12, na.rm = TRUE), n),
            tmin_lag1 = seq(min(dat$tmin_lag1, na.rm = TRUE),
                            max(dat$tmin_lag1, na.rm = TRUE), length.out = n),
            forest = rep(mean(dat$forest, na.rm = TRUE), n),
            wetlands = rep(mean(dat$wetlands, na.rm = TRUE), n),
            county = NA, site_id = NA, time_step = 180)
  )

# forest (highly significant)
eeev_mod_newdata_forest <- eeev_mod %>%
  predict(se_fit = TRUE, re_form = NA, re_form_iid = NA, type = "response",
          newdata = data.frame(
            prcp_lag1 = rep(mean(dat$prcp_lag1, na.rm = TRUE), n),
            prcp_lag5 = rep(mean(dat$prcp_lag5, na.rm = TRUE), n),
            prcp_lag12 = rep(mean(dat$prcp_lag12, na.rm = TRUE), n),
            tmax_lag6 = rep(mean(dat$tmax_lag6, na.rm = TRUE), n),
            tmax_lag12 = rep(mean(dat$tmax_lag12, na.rm = TRUE), n),
            tmin_lag1 = rep(mean(dat$tmin_lag1, na.rm = TRUE), n),
            forest = seq(min(dat$forest, na.rm = TRUE),
                         max(dat$forest, na.rm = TRUE), length.out = n),
            wetlands = rep(mean(dat$wetlands, na.rm = TRUE), n),
            county = NA, site_id = NA, time_step = 180)
  )

# wetlands (significant effect)
eeev_mod_newdata_wetlands <- eeev_mod %>%
  predict(se_fit = TRUE, re_form = NA, re_form_iid = NA, type = "response",
          newdata = data.frame(
            prcp_lag1 = rep(mean(dat$prcp_lag1, na.rm = TRUE), n),
            prcp_lag5 = rep(mean(dat$prcp_lag5, na.rm = TRUE), n),
            prcp_lag12 = rep(mean(dat$prcp_lag12, na.rm = TRUE), n),
            tmax_lag6 = rep(mean(dat$tmax_lag6, na.rm = TRUE), n),
            tmax_lag12 = rep(mean(dat$tmax_lag12, na.rm = TRUE), n),
            tmin_lag1 = rep(mean(dat$tmin_lag1, na.rm = TRUE), n),
            forest = rep(mean(dat$forest, na.rm = TRUE), n),
            wetlands = seq(min(dat$wetlands, na.rm = TRUE),
                           max(dat$wetlands, na.rm = TRUE), length.out = n),
            county = NA, site_id = NA, time_step = 180)
  )

# Alternative model - significant variables:
# tmin (marginally significant)
eeev_mod_alt_newdata_tmin <- eeev_mod_alt %>%
  predict(se_fit = TRUE, re_form = NA, re_form_iid = NA, type = "response",
          newdata = data.frame(
            tmin = seq(min(dat$tmin, na.rm = TRUE),
                       max(dat$tmin, na.rm = TRUE), length.out = n),
            tmin_lag6 = rep(mean(dat$tmin_lag6, na.rm = TRUE), n),
            tmin_lag12 = rep(mean(dat$tmin_lag12, na.rm = TRUE), n),
            tmax_lag1 = rep(mean(dat$tmax_lag1, na.rm = TRUE), n),
            forest = rep(mean(dat$forest, na.rm = TRUE), n),
            wetlands = rep(mean(dat$wetlands, na.rm = TRUE), n),
            county = NA, site_id = NA, time_step = 180)
  )

# tmin_lag6 (marginally significant)
eeev_mod_alt_newdata_tmin_lag6 <- eeev_mod_alt %>%
  predict(se_fit = TRUE, re_form = NA, re_form_iid = NA, type = "response",
          newdata = data.frame(
            tmin = rep(mean(dat$tmin, na.rm = TRUE), n),
            tmin_lag6 = seq(min(dat$tmin_lag6, na.rm = TRUE),
                            max(dat$tmin_lag6, na.rm = TRUE), length.out = n),
            tmin_lag12 = rep(mean(dat$tmin_lag12, na.rm = TRUE), n),
            tmax_lag1 = rep(mean(dat$tmax_lag1, na.rm = TRUE), n),
            forest = rep(mean(dat$forest, na.rm = TRUE), n),
            wetlands = rep(mean(dat$wetlands, na.rm = TRUE), n),
            county = NA, site_id = NA, time_step = 180)
  )

# tmin_lag12 (marginally significant)
eeev_mod_alt_newdata_tmin_lag12 <- eeev_mod_alt %>%
  predict(se_fit = TRUE, re_form = NA, re_form_iid = NA, type = "response",
          newdata = data.frame(
            tmin = rep(mean(dat$tmin, na.rm = TRUE), n),
            tmin_lag6 = rep(mean(dat$tmin_lag6, na.rm = TRUE), n),
            tmin_lag12 = seq(min(dat$tmin_lag12, na.rm = TRUE),
                             max(dat$tmin_lag12, na.rm = TRUE), length.out = n),
            tmax_lag1 = rep(mean(dat$tmax_lag1, na.rm = TRUE), n),
            forest = rep(mean(dat$forest, na.rm = TRUE), n),
            wetlands = rep(mean(dat$wetlands, na.rm = TRUE), n),
            county = NA, site_id = NA, time_step = 180)
  )

# tmax_lag1 (significant effect)
eeev_mod_alt_newdata_tmax_lag1 <- eeev_mod_alt %>%
  predict(se_fit = TRUE, re_form = NA, re_form_iid = NA, type = "response",
          newdata = data.frame(
            tmin = rep(mean(dat$tmin, na.rm = TRUE), n),
            tmin_lag6 = rep(mean(dat$tmin_lag6, na.rm = TRUE), n),
            tmin_lag12 = rep(mean(dat$tmin_lag12, na.rm = TRUE), n),
            tmax_lag1 = seq(min(dat$tmax_lag1, na.rm = TRUE),
                            max(dat$tmax_lag1, na.rm = TRUE), length.out = n),
            forest = rep(mean(dat$forest, na.rm = TRUE), n),
            wetlands = rep(mean(dat$wetlands, na.rm = TRUE), n),
            county = NA, site_id = NA, time_step = 180)
  )

# forest (highly significant)
eeev_mod_alt_newdata_forest <- eeev_mod_alt %>%
  predict(se_fit = TRUE, re_form = NA, re_form_iid = NA, type = "response",
          newdata = data.frame(
            tmin = rep(mean(dat$tmin, na.rm = TRUE), n),
            tmin_lag6 = rep(mean(dat$tmin_lag6, na.rm = TRUE), n),
            tmin_lag12 = rep(mean(dat$tmin_lag12, na.rm = TRUE), n),
            tmax_lag1 = rep(mean(dat$tmax_lag1, na.rm = TRUE), n),
            forest = seq(min(dat$forest, na.rm = TRUE),
                         max(dat$forest, na.rm = TRUE), length.out = n),
            wetlands = rep(mean(dat$wetlands, na.rm = TRUE), n),
            county = NA, site_id = NA, time_step = 180)
  )

# wetlands (significant effect)
eeev_mod_alt_newdata_wetlands <- eeev_mod_alt %>%
  predict(se_fit = TRUE, re_form = NA, re_form_iid = NA, type = "response",
          newdata = data.frame(
            tmin = rep(mean(dat$tmin, na.rm = TRUE), n),
            tmin_lag6 = rep(mean(dat$tmin_lag6, na.rm = TRUE), n),
            tmin_lag12 = rep(mean(dat$tmin_lag12, na.rm = TRUE), n),
            tmax_lag1 = rep(mean(dat$tmax_lag1, na.rm = TRUE), n),
            forest = rep(mean(dat$forest, na.rm = TRUE), n),
            wetlands = seq(min(dat$wetlands, na.rm = TRUE),
                           max(dat$wetlands, na.rm = TRUE), length.out = n),
            county = NA, site_id = NA, time_step = 180)
  )

# Create list of predictions
eeev_plot_preds <- list(
  eeev_mod_newdata_prcp_lag12,
  eeev_mod_newdata_tmax_lag6,
  eeev_mod_newdata_tmin_lag1,
  eeev_mod_newdata_forest,
  eeev_mod_newdata_wetlands,
  eeev_mod_alt_newdata_tmin,
  eeev_mod_alt_newdata_tmin_lag6,
  eeev_mod_alt_newdata_tmin_lag12,
  eeev_mod_alt_newdata_tmax_lag1,
  eeev_mod_alt_newdata_forest,
  eeev_mod_alt_newdata_wetlands)

# Save data
write_rds(eeev_plot_preds, "output/eeev_response_curves_plot_preds.rds")

# Read in response curves data
eeev_plot_preds <- read_rds("output/eeev_response_curves_plot_preds_ver01.rds")

# Create plots
# Helper function for cloglog inverse transformation
inv_cloglog_scaled <- function(x) {
  (1 - exp(-exp(x)))
}

# Primary model plot
eeev_primary_mod_plot <- ggarrange(
  #prcp_lag12
  ggplot(eeev_plot_preds[[1]], aes(x = prcp_lag12)) +
    geom_ribbon(aes(ymin = inv_cloglog_scaled(est - est_se), ymax = inv_cloglog_scaled(est + est_se)),
                fill = "darkcyan", alpha = 0.4) +
    geom_smooth(aes(y = inv_cloglog_scaled(est)), col = "grey20", se = F, linewidth = 0.5) +
    labs(x = "Precipitation (12 mo. lag)", y = NULL) +
    theme_bw(),
  #tmax_lag6  
  # ggplot(eeev_plot_preds[[2]], aes(x = tmax_lag6)) +
  #   geom_ribbon(aes(ymin = inv_cloglog_scaled(est - est_se), ymax = inv_cloglog_scaled(est + est_se)),
  #               fill = "darkcyan", alpha = 0.4) +
  #   geom_smooth(aes(y = inv_cloglog_scaled(est)), col = "grey20", se = F, linewidth = 0.5) +
  #   labs(x = "Max Temperature (6 mo. lag)", y = NULL) +
  #   theme_bw(),
  #tmin_lag1
  ggplot(eeev_plot_preds[[3]], aes(x = tmin_lag1)) +
    geom_ribbon(aes(ymin = inv_cloglog_scaled(est - est_se), ymax = inv_cloglog_scaled(est + est_se)),
                fill = "darkcyan", alpha = 0.4) +
    geom_smooth(aes(y = inv_cloglog_scaled(est)), col = "grey20", se = F, linewidth = 0.5) +
    labs(x = "Min Temperature (1 mo. lag)", y = NULL) +
    theme_bw(),
  #forest
  ggplot(eeev_plot_preds[[4]], aes(x = forest)) +
    geom_ribbon(aes(ymin = inv_cloglog_scaled(est - est_se), ymax = inv_cloglog_scaled(est + est_se)),
                fill = "darkcyan", alpha = 0.4) +
    geom_smooth(aes(y = inv_cloglog_scaled(est)), col = "grey20", se = F, linewidth = 0.5) +
    labs(x = "Forest Cover", y = NULL) +
    theme_bw(),
  #wetlands
  ggplot(eeev_plot_preds[[5]], aes(x = wetlands)) +
    geom_ribbon(aes(ymin = inv_cloglog_scaled(est - est_se), ymax = inv_cloglog_scaled(est + est_se)),
                fill = "darkcyan", alpha = 0.4) +
    geom_smooth(aes(y = inv_cloglog_scaled(est)), col = "grey20", se = F, linewidth = 0.5) +
    labs(x = "Wetlands Cover", y = NULL) +
    theme_bw(),
  align = "hv",
  nrow = 4, ncol = 1)

eeev_primary_mod_plot_ann <- annotate_figure(
  eeev_primary_mod_plot,
  left = grid::textGrob("EEEV Seroprevalence",
                        rot = 90,
                        hjust = 0.5,
                        vjust = 0.5,
                        gp = grid::gpar(cex = 1, fontface = "bold")))

# Alternate model plot  
eeev_alt_mod_plot <- ggarrange(
  # #tmin
  # ggplot(eeev_plot_preds[[6]], aes(x = tmin)) +
  #   geom_ribbon(aes(ymin = inv_cloglog_scaled(est - est_se), ymax = inv_cloglog_scaled(est + est_se)),
  #               fill = "orange", alpha = 0.4) +
  #   geom_smooth(aes(y = inv_cloglog_scaled(est)), col = "grey20", se = F, linewidth = 0.5) +
  #   labs(x = "Min Temperature", y = NULL) +
  #   theme_bw(),
  # #tmin_lag6
  # ggplot(eeev_plot_preds[[7]], aes(x = tmin_lag6)) +
  #   geom_ribbon(aes(ymin = inv_cloglog_scaled(est - est_se), ymax = inv_cloglog_scaled(est + est_se)),
  #               fill = "orange", alpha = 0.4) +
  #   geom_smooth(aes(y = inv_cloglog_scaled(est)), col = "grey20", se = F, linewidth = 0.5) +
  #   labs(x = "Min Temperature (6 mo. lag)", y = NULL) +
  #   theme_bw(),
  # #tmin_lag12
  # ggplot(eeev_plot_preds[[8]], aes(x = tmin_lag12)) +
  #   geom_ribbon(aes(ymin = inv_cloglog_scaled(est - est_se), ymax = inv_cloglog_scaled(est + est_se)),
  #               fill = "orange", alpha = 0.4) +
  #   geom_smooth(aes(y = inv_cloglog_scaled(est)), col = "grey20", se = F, linewidth = 0.5) +
  #   labs(x = "Min Temperature (12 mo. lag)", y = NULL) +
  #   theme_bw(),
  #tmax_lag1
  ggplot(eeev_plot_preds[[9]], aes(x = tmax_lag1)) +
    geom_ribbon(aes(ymin = inv_cloglog_scaled(est - est_se), ymax = inv_cloglog_scaled(est + est_se)),
                fill = "orange", alpha = 0.4) +
    geom_smooth(aes(y = inv_cloglog_scaled(est)), col = "grey20", se = F, linewidth = 0.5) +
    labs(x = "Max Temperature (1 mo. lag)", y = NULL) +
    theme_bw(),
  #forest
  ggplot(eeev_plot_preds[[10]], aes(x = forest)) +
    geom_ribbon(aes(ymin = inv_cloglog_scaled(est - est_se), ymax = inv_cloglog_scaled(est + est_se)),
                fill = "orange", alpha = 0.4) +
    geom_smooth(aes(y = inv_cloglog_scaled(est)), col = "grey20", se = F, linewidth = 0.5) +
    labs(x = "Forest Cover", y = NULL) +
    theme_bw(),
  #wetlands
  ggplot(eeev_plot_preds[[11]], aes(x = wetlands)) +
    geom_ribbon(aes(ymin = inv_cloglog_scaled(est - est_se), ymax = inv_cloglog_scaled(est + est_se)),
                fill = "orange", alpha = 0.4) +
    geom_smooth(aes(y = inv_cloglog_scaled(est)), col = "grey20", se = F, linewidth = 0.5) +
    labs(x = "Wetlands Cover", y = NULL) +
    theme_bw(),
  align = "hv",
  nrow = 3, ncol = 1)

eeev_alt_mod_plot_ann <- annotate_figure(
  eeev_alt_mod_plot,
  left = grid::textGrob("EEEV Seroprevalence",
                        rot = 90,
                        hjust = 0.5,
                        vjust = 0.5,
                        gp = grid::gpar(cex = 1, fontface = "bold")))

eeev_primary_mod_plot_ann
eeev_alt_mod_plot_ann

# Save plots
ggsave("figs/eeev/primary_model_response_curves_ver02.jpeg", 
       eeev_primary_mod_plot_ann, 
       height = 10, 
       width = 5, 
       units = "in", 
       dpi = 600)
ggsave("figs/eeev/alt_model_response_curves_ver02.jpeg", 
       eeev_alt_mod_plot_ann, 
       height = 7.5, 
       width = 4, 
       units = "in", 
       dpi = 600)