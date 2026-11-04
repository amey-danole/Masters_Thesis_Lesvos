# Wide-angle behavioural analysis: bee activity and microclimate
#
# Questions
#   1. How does overall bee activity (mean number of bees visible at the
#      nesting aggregation) relate to air temperature and to soil temperature
#      at three depths?
#   2. Do these relationships differ between the two sites?
#   3. Are the pooled results robust to how sampling date is handled?
#
# Inputs (working directory)
#   Final_combined_temperature_data.csv   iButton readings from all loggers
#   Wide_behavioural_data.csv             wide-angle camera counts
#   Outdoor_temperature_logger.csv        HOBO air temperature (from 19 April)
#
# Logger key (shared with the focal analysis)
#   Thermochron_1        air temperature at 2 m
#   Thermochron_3 / 6    shallow soil, 2-3 cm  (Site 1 / Site 2)
#   Thermochron_7 / 8    middle soil, 5 cm     (Site 1 / Site 2)
#   Thermochron_9 / 10   deep soil, 10 cm      (Site 1 / Site 2)
#
# Air temperature is taken from the iButton (Thermochron_1), which recorded
# throughout the study. The HOBO logger started on 19 April and is used only
# to check the iButton.
#
# Outputs
#   wide_pooled_gam_summary.csv, wide_site_difference.csv,
#   wide_date_adjustment_comparison.csv, wide_within_day_gam.csv,
#   wide_k_sensitivity.csv, wide_day_structure.csv
#   wide_behavioural_gam_figure.png, wide_behavioural_gam_<depth>_site<n>.png
#   diagnostics/wide_dharma_summary.csv and residual plots


# Packages ----

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(mgcv)
library(lme4)
library(DHARMa)

set.seed(1)


# Shared settings ----

# One colour per thermal layer, used in every figure
layer_cols <- c(
  "Ambient air"           = "#0072B2",
  "Shallow soil (2–3 cm)" = "#E69F00",
  "Middle soil (5 cm)"    = "#009E73",
  "Deep soil (10 cm)"     = "#CC79A7"
)

# Common look for the manuscript figures
theme_ms <- function(base_size = 15) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title   = element_text(face = "bold", size = base_size + 3),
      axis.text    = element_text(face = "bold", size = base_size),
      legend.title = element_text(face = "bold", size = base_size + 1),
      legend.text  = element_text(face = "bold", size = base_size - 1)
    )
}

dir.create("diagnostics", showWarnings = FALSE)


# Load data ----

temperature_raw <- read.csv("Final_combined_temperature_data.csv")
wide_raw        <- read.csv("Wide_behavioural_data.csv")
outdoor_temp    <- read.csv("Outdoor_temperature_logger.csv")


# Build the analysis dataset ----

# One column per iButton, then attach the HOBO air temperature
env_wide <- temperature_raw %>%
  select(Date, Time, Sensors, Avg_sensor_temperature) %>%
  group_by(Date, Time, Sensors) %>%
  summarise(Avg_sensor_temperature = mean(Avg_sensor_temperature, na.rm = TRUE),
            .groups = "drop") %>%
  pivot_wider(names_from = Sensors, values_from = Avg_sensor_temperature) %>%
  left_join(outdoor_temp, by = c("Date", "Time"))

wide_data <- wide_raw %>%
  left_join(env_wide, by = c("Date", "Time")) %>%
  mutate(
    Date   = as.factor(Date),
    Site   = as.factor(Site),
    Camera = as.factor(Camera),
    Time   = as.numeric(Time)
  )


# Recording counts and data availability ----

cat("Wide-angle recordings:", nrow(wide_raw), "\n")

cat("\nRecordings by site and camera:\n")
print(table(Site = wide_raw$Site, Camera = wide_raw$Camera, useNA = "ifany"))

cat("\nRecordings by start hour and site:\n")
print(table(Start_hour = wide_raw$Time, Site = wide_raw$Site, useNA = "ifany"))

cat("\nRecordings by date and site:\n")
print(table(Date = wide_raw$Date, Site = wide_raw$Site, useNA = "ifany"))

cat("\nMissing activity or logger values:\n")
wide_data %>%
  summarise(
    missing_bees          = sum(is.na(Avg_number_of_bees)),
    missing_ibutton_air   = sum(is.na(Thermochron_1)),
    missing_hobo_air      = sum(is.na(Avg_T_ambient)),
    missing_shallow       = sum(is.na(Thermochron_3) & is.na(Thermochron_6)),
    missing_middle        = sum(is.na(Thermochron_7) & is.na(Thermochron_8)),
    missing_deep          = sum(is.na(Thermochron_9) & is.na(Thermochron_10))
  ) %>%
  as.data.frame() %>%
  print()

# Activity is an average of five whole-number counts (2, 4, 6, 8, 10 min),
# so it takes values in steps of 0.2 and is analysed as a continuous index
sort(unique(wide_raw$Avg_number_of_bees))


# Air temperature source ----

# The iButton (Thermochron_1) covers the whole study; the HOBO logger only
# starts on 19 April. Where both exist, how closely do they agree?
ambient_check <- wide_data %>%
  filter(!is.na(Thermochron_1), !is.na(Avg_T_ambient)) %>%
  summarise(
    n_recordings    = n(),
    r               = cor(Thermochron_1, Avg_T_ambient),
    mean_difference = mean(Thermochron_1 - Avg_T_ambient),
    sd_difference   = sd(Thermochron_1 - Avg_T_ambient)
  )
ambient_check


# Pooled GAMs: activity against each thermal layer ----

# Recordings from the same day share weather, so every model includes site
# (fixed) and sampling date (random intercept). Soil layers average the two
# sites' loggers; complete cases only.
gam_data <- wide_data %>%
  mutate(
    ambient = Thermochron_1,
    shallow = rowMeans(cbind(Thermochron_3, Thermochron_6),  na.rm = TRUE),
    middle  = rowMeans(cbind(Thermochron_7, Thermochron_8),  na.rm = TRUE),
    deep    = rowMeans(cbind(Thermochron_9, Thermochron_10), na.rm = TRUE)
  ) %>%
  select(Date, Site, Avg_number_of_bees, ambient, shallow, middle, deep) %>%
  drop_na() %>%
  mutate(Date = droplevels(factor(Date)), Site = droplevels(factor(Site)))

cat("\nRecordings used in the pooled GAMs:", nrow(gam_data),
    "on", n_distinct(gam_data$Date), "dates\n")
print(table(Date = gam_data$Date, Site = gam_data$Site))

layers <- tribble(
  ~variable,  ~layer,
  "ambient",  "Ambient air",
  "shallow",  "Shallow soil (2–3 cm)",
  "middle",   "Middle soil (5 cm)",
  "deep",     "Deep soil (10 cm)"
)

pooled_gams <- lapply(layers$variable, function(v) {
  gam(as.formula(paste0("Avg_number_of_bees ~ Site + s(", v,
                        ", k = 5) + s(Date, bs = 're')")),
      data = gam_data, method = "REML")
})
names(pooled_gams) <- layers$layer

# Statistics for the temperature smooth (first row of the smooth table), one
# row per layer. Deviance explained includes the date effect, so it is not
# the share explained by temperature alone and should not be quoted as such.
pooled_gam_table <- bind_rows(lapply(names(pooled_gams), function(layer) {
  s <- summary(pooled_gams[[layer]])
  data.frame(
    layer         = layer,
    n             = nrow(gam_data),
    edf           = s$s.table[1, "edf"],
    F             = s$s.table[1, "F"],
    p             = s$s.table[1, "p-value"],
    r_sq_adj      = s$r.sq,
    deviance_expl = 100 * s$dev.expl
  )
}))
pooled_gam_table

write.csv(pooled_gam_table, "wide_pooled_gam_summary.csv", row.names = FALSE)


# Robustness 1: with and without the date term ----

# Same 166 recordings, same smooth, fitted by ML so the models can be compared.
# The without-date model treats every recording as independent.
fit_pair <- function(v) {
  m_no   <- gam(as.formula(paste0("Avg_number_of_bees ~ Site + s(", v, ", k = 5)")),
                data = gam_data, method = "ML")
  m_date <- gam(as.formula(paste0("Avg_number_of_bees ~ Site + s(", v,
                                  ", k = 5) + s(Date, bs = 're')")),
                data = gam_data, method = "ML")
  a <- summary(m_no)$s.table[1, ]
  b <- summary(m_date)$s.table[1, ]
  data.frame(layer = v,
             edf_no_date = a["edf"], p_no_date = a["p-value"], AIC_no_date = AIC(m_no),
             edf_with_date = b["edf"], p_with_date = b["p-value"],
             AIC_with_date = AIC(m_date))
}
date_adjustment <- bind_rows(lapply(layers$variable, fit_pair))
rownames(date_adjustment) <- NULL
date_adjustment

write.csv(date_adjustment, "wide_date_adjustment_comparison.csv", row.names = FALSE)


# Robustness 2: how much of each temperature varies between days? ----

# Share of the variance in each thermal layer that lies between days. A high
# share means the layers are mostly describing which day it was.
day_structure <- bind_rows(lapply(layers$variable, function(v) {
  d <- gam_data
  d$x <- d[[v]]
  vc <- as.data.frame(VarCorr(lmer(x ~ 1 + (1 | Date), data = d)))
  data.frame(layer = v,
             between_day_var = vc$vcov[vc$grp == "Date"],
             within_day_var  = vc$vcov[vc$grp == "Residual"]) %>%
    mutate(share_between_days = between_day_var / (between_day_var + within_day_var))
}))
day_structure

write.csv(day_structure, "wide_day_structure.csv", row.names = FALSE)


# Robustness 3: within-day effect of temperature ----

# Each thermal layer is centred within its sampling date, so the smooth
# estimates how activity changes across hours of the same day. It cannot be
# driven by differences between days.
within_gams <- lapply(layers$variable, function(v) {
  d <- gam_data %>%
    group_by(Date) %>%
    mutate(within = .data[[v]] - mean(.data[[v]])) %>%
    ungroup()
  gam(Avg_number_of_bees ~ Site + s(within, k = 5) + s(Date, bs = "re"),
      data = d, method = "REML")
})
names(within_gams) <- layers$variable

within_day_table <- bind_rows(lapply(names(within_gams), function(v) {
  s <- summary(within_gams[[v]])
  data.frame(layer = v, n = nrow(gam_data),
             edf = s$s.table[1, "edf"], F = s$s.table[1, "F"],
             p = s$s.table[1, "p-value"])
}))
within_day_table

write.csv(within_day_table, "wide_within_day_gam.csv", row.names = FALSE)


# Robustness 4: basis dimension ----

# k = 5 caps the wiggliness of each smooth. If the results hold at k = 8, the
# cap is not driving them.
k_sensitivity <- bind_rows(lapply(layers$variable, function(v) {
  bind_rows(lapply(c(5, 8), function(k) {
    m <- gam(as.formula(paste0("Avg_number_of_bees ~ Site + s(", v, ", k = ", k,
                               ") + s(Date, bs = 're')")),
             data = gam_data, method = "REML")
    s  <- summary(m)
    kc <- k.check(m)
    data.frame(layer = v, k = k,
               edf = s$s.table[1, "edf"], F = s$s.table[1, "F"],
               p = s$s.table[1, "p-value"], k_index = kc[1, "k-index"])
  }))
}))
k_sensitivity

write.csv(k_sensitivity, "wide_k_sensitivity.csv", row.names = FALSE)


# Pooled GAM figure ----

# Fitted smooth with 95% confidence band across the observed range, with the
# date effect excluded and site held at its first level
gam_curve <- function(model, data, variable, n = 200) {
  grid <- data.frame(
    seq(
      min(data[[variable]]),
      max(data[[variable]]),
      length.out = n
    )
  )
  names(grid) <- variable
  grid$Site <- factor(levels(data$Site)[1], levels = levels(data$Site))
  grid$Date <- factor(levels(data$Date)[1], levels = levels(data$Date))
  
  p <- predict(
    model,
    newdata = grid,
    se.fit = TRUE,
    exclude = "s(Date)"
  )
  
  data.frame(
    x = grid[[variable]],
    fit = p$fit,
    lower = p$fit - 1.96 * p$se.fit,
    upper = p$fit + 1.96 * p$se.fit
  )
}

pooled_curves <- bind_rows(lapply(seq_len(nrow(layers)), function(i) {
  gam_curve(
    pooled_gams[[layers$layer[i]]],
    gam_data,
    layers$variable[i]
  ) %>%
    mutate(layer = layers$layer[i])
})) %>%
  mutate(layer = factor(layer, levels = layers$layer))

gam_plot <- ggplot(
  pooled_curves,
  aes(x, fit, colour = layer, fill = layer)
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.20,
    colour = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(values = layer_cols) +
  scale_fill_manual(values = layer_cols) +
  labs(
    x = "Temperature (°C)",
    y = "Relative bee activity",
    colour = "Thermal environment",
    fill = "Thermal environment"
  ) +
  theme_ms() +
  theme(
    axis.title.x = element_text(size = 19, face = "bold"),
    axis.title.y = element_text(size = 19, face = "bold"),
    axis.text.x = element_text(size = 17, face = "bold"),
    axis.text.y = element_text(size = 17, face = "bold"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 18, face = "bold")
  )

gam_plot

ggsave("wide_behavioural_gam_figure.png", gam_plot,
       dpi = 600, width = 400, height = 240, units = "mm")


# Do the temperature relationships differ between sites? ----

# Each site is paired with its own soil loggers; air temperature is shared
site_data <- wide_data %>%
  mutate(
    ambient = Thermochron_1,
    shallow = case_when(Site == "1" ~ Thermochron_3,
                        Site == "2" ~ Thermochron_6),
    middle  = case_when(Site == "1" ~ Thermochron_7,
                        Site == "2" ~ Thermochron_8),
    deep    = case_when(Site == "1" ~ Thermochron_9,
                        Site == "2" ~ Thermochron_10)
  ) %>%
  select(Avg_number_of_bees, Site, Date, ambient, shallow, middle, deep)

# Likelihood-ratio comparison: one shared smooth vs a smooth per site,
# with a random intercept for sampling date
test_site_difference <- function(data, variable) {
  d <- data %>%
    select(Avg_number_of_bees, Site, Date, all_of(variable)) %>%
    drop_na()
  
  shared <- gam(
    as.formula(paste0("Avg_number_of_bees ~ Site + s(", variable,
                      ", k = 5) + s(Date, bs = 're')")),
    data = d, method = "ML")
  
  by_site <- gam(
    as.formula(paste0("Avg_number_of_bees ~ Site + s(", variable,
                      ", by = Site, k = 5) + s(Date, bs = 're')")),
    data = d, method = "ML")
  
  list(shared = shared, by_site = by_site, n = nrow(d),
       comparison = anova(shared, by_site, test = "Chisq"))
}

site_predictors <- c("ambient", "shallow", "middle", "deep")

site_tests <- setNames(
  lapply(site_predictors, test_site_difference, data = site_data),
  site_predictors
)

# Chi-square, df and p for each layer. A test with a non-positive chi-square
# or df is invalid (the site-specific model used fewer effective degrees of
# freedom than the shared model) and is set to NA before the Holm correction.
site_difference_table <- bind_rows(lapply(site_predictors, function(v) {
  cmp <- site_tests[[v]]$comparison
  data.frame(
    predictor  = v,
    n          = site_tests[[v]]$n,
    chi_square = cmp$Deviance[2],
    df         = cmp$Df[2],
    p          = cmp[2, "Pr(>Chi)"]
  )
})) %>%
  mutate(
    p      = ifelse(chi_square > 0 & df > 0, p, NA),
    p_holm = p.adjust(p, method = "holm")
  )

site_difference_table

write.csv(site_difference_table, "wide_site_difference.csv", row.names = FALSE)

# Site-specific smooth summaries (edf, F and p for each site's smooth)
lapply(site_tests, function(x) summary(x$by_site))


# Site-specific GAM figures ----

fit_site_gam <- function(data, variable) {
  d <- data %>%
    filter(!is.na(Avg_number_of_bees), !is.na(.data[[variable]]), !is.na(Date)) %>%
    mutate(Date = droplevels(Date))
  
  list(
    model = gam(
      as.formula(paste0("Avg_number_of_bees ~ s(", variable,
                        ", k = 5) + s(Date, bs = 're')")),
      data = d, method = "REML"),
    data = d
  )
}

# Predictions with the date random effect excluded
site_curve <- function(fit, variable, n = 200) {
  d <- fit$data
  grid <- data.frame(seq(min(d[[variable]]), max(d[[variable]]), length.out = n))
  names(grid) <- variable
  grid$Date <- factor(levels(d$Date)[1], levels = levels(d$Date))
  
  p <- predict(fit$model, newdata = grid, se.fit = TRUE, exclude = "s(Date)")
  data.frame(
    temperature = grid[[variable]],
    fit   = p$fit,
    lower = p$fit - 1.96 * p$se.fit,
    upper = p$fit + 1.96 * p$se.fit
  )
}

plot_site_curve <- function(curve, site_label, depth_label) {
  ggplot(curve, aes(temperature, fit)) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      fill = layer_cols[["Shallow soil (2–3 cm)"]],
      alpha = 0.20
    ) +
    geom_line(
      colour = layer_cols[["Shallow soil (2–3 cm)"]],
      linewidth = 1.2
    ) +
    labs(
      title = site_label,
      x = paste0(depth_label, " soil temperature (°C)"),
      y = "Relative bee activity"
    ) +
    theme_ms() +
    theme(
      plot.title = element_text(face = "bold", size = 21, hjust = 0.5),
      axis.title.x = element_text(size = 19, face = "bold"),
      axis.title.y = element_text(size = 19, face = "bold"),
      axis.text.x = element_text(size = 17, face = "bold"),
      axis.text.y = element_text(size = 17, face = "bold")
    )
}

site_figures <- tribble(
  ~site, ~variable, ~depth_label,        ~file,
  "1",   "shallow", "Shallow (2–3 cm)",  "wide_behavioural_gam_shallow_site1.png",
  "2",   "shallow", "Shallow (2–3 cm)",  "wide_behavioural_gam_shallow_site2.png",
  "1",   "middle",  "Middle (5 cm)",     "wide_behavioural_gam_middle_site1.png",
  "2",   "middle",  "Middle (5 cm)",     "wide_behavioural_gam_middle_site2.png"
)

site_plots <- lapply(seq_len(nrow(site_figures)), function(i) {
  row <- site_figures[i, ]
  
  fit <- fit_site_gam(
    filter(site_data, Site == row$site),
    row$variable
  )
  curve <- site_curve(fit, row$variable)
  p <- plot_site_curve(
    curve,
    paste("Site", row$site),
    row$depth_label
  )
  
  ggsave(
    row$file,
    p,
    dpi = 600,
    width = 200,
    height = 160,
    units = "mm"
  )
  p
})

site_plots


# Model diagnostics: basis dimension ----

# k-index close to 1 (and a non-small p) indicates an adequate basis dimension
lapply(pooled_gams, k.check)
lapply(site_tests, function(x) k.check(x$by_site))


# Model diagnostics: residuals ----

# Simulation-based residual checks (DHARMa): KS test of the residual
# distribution, a dispersion test and an outlier test. Plots are saved.
check_residuals <- function(model, name) {
  tryCatch({
    sim <- simulateResiduals(model, n = 1000, plot = FALSE)
    stats <- data.frame(
      analysis     = "wide-angle",
      model        = name,
      ks_p         = testUniformity(sim, plot = FALSE)$p.value,
      dispersion_p = testDispersion(sim, plot = FALSE)$p.value,
      outlier_p    = testOutliers(sim, plot = FALSE)$p.value
    )
    
    png(file.path("diagnostics", paste0("wide_", name, ".png")),
        width = 1800, height = 900, res = 200)
    tryCatch(plot(sim), finally = invisible(dev.off()))
    stats
  }, error = function(e) {
    cat("DHARMa failed for", name, ":", conditionMessage(e), "\n")
    data.frame(analysis = "wide-angle", model = name,
               ks_p = NA, dispersion_p = NA, outlier_p = NA)
  })
}

wide_models <- c(
  setNames(pooled_gams, paste0("pooled_", site_predictors)),
  setNames(lapply(site_tests, `[[`, "by_site"), paste0("by_site_", site_predictors)),
  setNames(within_gams, paste0("within_day_", site_predictors))
)

wide_diagnostics <- bind_rows(lapply(names(wide_models), function(nm) {
  check_residuals(wide_models[[nm]], nm)
}))
wide_diagnostics

write.csv(wide_diagnostics, "diagnostics/wide_dharma_summary.csv", row.names = FALSE)

sessionInfo()