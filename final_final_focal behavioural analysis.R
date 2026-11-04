# Focal behavioural analysis: nest guarding, abundance, nest traffic and timing
#
# Questions
#   1. Does nest guarding in H. scabiosae (HS) depend on soil temperature?
#   2. Do HS and L. malachurum (LM) differ in abundance, and in how abundance
#      responds to temperature and light?
#   3. How do nest entry and exit frequencies relate to the thermal
#      environment in each species?
#   4. Do the species differ in WHEN they are active, and does time of day
#      explain activity beyond temperature?
#
# Inputs (working directory)
#   focal_behavioural_data.csv            focal observations
#   combined_environmental_variables.csv  loggers, matched by Date and Time
#
# Logger key (shared with the wide-angle analysis)
#   Thermochron_1        air temperature at 2 m
#   Thermochron_2 / 4    surface (Site 1 / Site 2)
#   Thermochron_3 / 6    shallow soil, 2-3 cm (Site 1 / Site 2)
#   Thermochron_7 / 8    middle soil, 5 cm    (Site 1 / Site 2)
#   Thermochron_9 / 10   deep soil, 10 cm     (Site 1 / Site 2)
#
# Sampling: the first 10 minutes of every hour. Abundance is the mean of five
# counts per recording (values are multiples of 0.2). Light data cover 95 of
# the 163 observations; analyses using light use those rows only.
#
# Outputs
#   Combined_focal_behavioural_dataset.csv, movement_gam_summary.csv
#   Figures (png) and diagnostics/focal_dharma_summary.csv with residual plots


# Packages ----

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(mgcv)
library(MuMIn)
library(performance)
library(emmeans)
library(DHARMa)


# Shared settings and helpers ----

species_cols <- c("HS" = "#E69F00", "LM" = "#0072B2")
col_HS <- species_cols[["HS"]]
col_LM <- species_cols[["LM"]]

species_labels <- c(
  "HS" = expression(italic("Halictus scabiosae")),
  "LM" = expression(italic("Lasioglossum malachurum"))
)

theme_ms <- function(base_size = 14, axis_face = "bold") {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = axis_face)
    )
}

scale1 <- function(x) as.numeric(scale(x))

# Average the paired site loggers into one value per layer
add_microclimate <- function(d) {
  d %>%
    mutate(
      ambient_temp = Thermochron_1,
      surface_temp = rowMeans(cbind(Thermochron_2, Thermochron_4),  na.rm = TRUE),
      shallow_temp = rowMeans(cbind(Thermochron_3, Thermochron_6),  na.rm = TRUE),
      middle_temp  = rowMeans(cbind(Thermochron_7, Thermochron_8),  na.rm = TRUE),
      deep_temp    = rowMeans(cbind(Thermochron_9, Thermochron_10), na.rm = TRUE)
    )
}

# Standardise each layer (and light) to mean 0, SD 1
add_scaled_env <- function(d) {
  d %>%
    mutate(
      ambient_s = scale1(ambient_temp),
      surface_s = scale1(surface_temp),
      shallow_s = scale1(shallow_temp),
      middle_s  = scale1(middle_temp),
      deep_s    = scale1(deep_temp),
      light_s   = scale1(Avg_lux_light_logger)
    )
}

# Environmental layers, in the order used for the multi-panel figures
env_layers <- tribble(
  ~variable,   ~title,                  ~xlab,                     ~movement_title,
  "ambient_s", "Ambient air",           "Temperature (scaled)",    "Ambient temperature",
  "surface_s", "Surface soil",          "Temperature (scaled)",    "Surface soil temperature",
  "shallow_s", "Shallow soil (2–3 cm)", "Temperature (scaled)",    "Shallow soil temperature",
  "middle_s",  "Middle soil (5 cm)",    "Temperature (scaled)",    "Middle soil temperature",
  "deep_s",    "Deep soil (10 cm)",     "Temperature (scaled)",    "Deep soil temperature",
  "light_s",   "Light intensity",       "Light intensity (scaled)", "Light intensity"
)


# Load and merge data ----

focal_raw    <- read.csv("focal_behavioural_data.csv")
combined_env <- read.csv("combined_environmental_variables.csv")

focal_raw <- focal_raw %>%
  mutate(Date = as.character(Date), Time = as.numeric(Time))

combined_env <- combined_env %>%
  mutate(Date = as.character(Date), Time = as.numeric(Time))

focal <- inner_join(focal_raw, combined_env, by = c("Date", "Time"))

# Merged dataset, saved for reference
write.csv(focal, "Combined_focal_behavioural_dataset.csv", row.names = FALSE)


# Data availability ----

# Each focal row is joined to the loggers by Date and Time; rows without a
# logger match are dropped by the inner join.
cat("Focal rows (raw data):    ", nrow(focal_raw), "\n")
cat("Logger rows:              ", nrow(combined_env), "\n")
cat("Focal rows after the join:", nrow(focal), "\n")

lost <- anti_join(focal_raw, combined_env, by = c("Date", "Time"))
cat("Focal rows without a logger match:", nrow(lost), "\n")
print(as.data.frame(lost))

# Duplicate keys would inflate the join
combined_env %>% count(Date, Time) %>% filter(n > 1) %>% as.data.frame() %>% print()
focal_raw %>%
  count(across(any_of(c("Date", "Site", "Day", "Time")))) %>%
  filter(n > 1) %>%
  as.data.frame() %>%
  print()

# Observations analysed, by hour and site and by day and site
print(table(Hour = focal$Time, Site = focal$Site, useNA = "ifany"))
print(table(Day = focal$Day, Site = focal$Site, useNA = "ifany"))

# Missing logger and response values
colSums(is.na(focal[, c("Thermochron_1", "Thermochron_2", "Thermochron_3",
                        "Thermochron_4", "Thermochron_6", "Thermochron_7",
                        "Thermochron_8", "Thermochron_9", "Thermochron_10",
                        "Avg_lux_light_logger",
                        "Avg_total_number_of_HS", "Avg_total_number_of_LM",
                        "Incoming_frequency_HS", "Outgoing_frequency_HS",
                        "Incoming_frequency_LM", "Outgoing_frequency_LM")]))

# The surface and light loggers have gaps. Analyses that use them are
# restricted to complete rows, so the affected days are listed here.
focal %>%
  group_by(Day, Site) %>%
  summarise(n = n(),
            missing_surface = sum(is.na(Thermochron_2) & is.na(Thermochron_4)),
            missing_light   = sum(is.na(Avg_lux_light_logger)),
            .groups = "drop") %>%
  print(n = Inf)


# Nest guarding in H. scabiosae: occurrence ----

# Binomial GLMM for the presence of guarding against shallow soil
# temperature (2-3 cm); a missing duration means no guarding was recorded
guard_occurrence_data <- focal %>%
  add_microclimate() %>%
  mutate(
    guard_present  = ifelse(is.na(Avg_duration_of_nest_guarding_HS), 0, 1),
    shallow_scaled = scale1(shallow_temp)
  ) %>%
  filter(!is.na(Site), !is.na(Day)) %>%
  mutate(Site = factor(Site), Day = factor(Day))

m_occurrence <- glmer(
  guard_present ~ shallow_scaled + (1 | Site) + (1 | Day),
  family = binomial,
  data = guard_occurrence_data
)
summary(m_occurrence)


# Nest guarding in H. scabiosae: duration ----

# Linear mixed model for log guarding duration against middle soil
# temperature (5 cm), using bouts where guarding occurred
guard_duration_data <- focal %>%
  add_microclimate() %>%
  mutate(middle_scaled = scale1(middle_temp)) %>%
  filter(!is.na(Avg_duration_of_nest_guarding_HS),
         Avg_duration_of_nest_guarding_HS > 0) %>%
  mutate(
    log_guard = log(Avg_duration_of_nest_guarding_HS),
    Site = factor(Site),
    Day  = factor(Day)
  )

m_duration <- lmer(
  log_guard ~ middle_scaled + (1 | Site) + (1 | Day),
  data = guard_duration_data,
  REML = TRUE
)
summary(m_duration)

# Layer-by-layer check. The five temperature layers are strongly correlated
# (r up to about 0.9), so they cannot be entered together; each layer gets its
# own model with the same random effects, and Holm-adjusted p-values correct
# for the five comparisons. n differs between layers because some loggers have
# missing values, so AIC is not compared across these models.
guard_by_layer <- bind_rows(lapply(
  c(ambient = "ambient_temp", surface = "surface_temp", shallow = "shallow_temp",
    middle  = "middle_temp",  deep    = "deep_temp"),
  function(v) {
    d <- focal %>%
      add_microclimate() %>%
      mutate(z = scale1(.data[[v]])) %>%
      filter(!is.na(Avg_duration_of_nest_guarding_HS),
             Avg_duration_of_nest_guarding_HS > 0, !is.na(z)) %>%
      mutate(log_guard = log(Avg_duration_of_nest_guarding_HS),
             Site = factor(Site), Day = factor(Day))
    
    m  <- lmer(log_guard ~ z + (1 | Site) + (1 | Day), data = d)
    cf <- coef(summary(m))
    data.frame(n    = nrow(d),
               beta = cf["z", "Estimate"],
               se   = cf["z", "Std. Error"],
               t    = cf["z", "t value"],
               p    = cf["z", "Pr(>|t|)"])
  }), .id = "layer") %>%
  mutate(p_holm = p.adjust(p, method = "holm"))

guard_by_layer


# Nest guarding figure ----

# Population-level predictions (random effects excluded) over +/- 2 SD.
# predict() with se.fit for mixed models needs a recent version of lme4.
occurrence_curve <- data.frame(shallow_scaled = seq(-2, 2, length.out = 200))
pr <- predict(m_occurrence, newdata = occurrence_curve,
              type = "link", se.fit = TRUE, re.form = NA)
occurrence_curve <- occurrence_curve %>%
  mutate(fit       = plogis(pr$fit),
         conf.low  = plogis(pr$fit - 1.96 * pr$se.fit),
         conf.high = plogis(pr$fit + 1.96 * pr$se.fit))

duration_curve <- data.frame(middle_scaled = seq(-2, 2, length.out = 200))

pr <- predict(m_duration, newdata = duration_curve, se.fit = TRUE, re.form = NA)
duration_curve <- duration_curve %>%
  mutate(fit       = pr$fit,
         conf.low  = pr$fit - 1.96 * pr$se.fit,
         conf.high = pr$fit + 1.96 * pr$se.fit)

plot_guarding_panel <- function(curve, x_var, title, xlab, ylab) {
  ggplot(curve, aes(x = .data[[x_var]])) +
    geom_ribbon(
      aes(ymin = conf.low, ymax = conf.high),
      fill = col_HS,
      alpha = 0.25
    ) +
    geom_line(
      aes(y = fit),
      colour = col_HS,
      linewidth = 1.2
    ) +
    labs(
      title = title,
      x = xlab,
      y = ylab
    ) +
    theme_ms(16) +
    theme(
      plot.title = element_text(size = 22, face = "bold"),
      axis.title.x = element_text(size = 22, face = "bold"),
      axis.title.y = element_text(size = 22, face = "bold"),
      axis.text.x = element_text(size = 18, face = "bold"),
      axis.text.y = element_text(size = 18, face = "bold"),
      legend.title = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 18, face = "bold")
    )
}

guarding_figure <-
  plot_guarding_panel(
    occurrence_curve, "shallow_scaled",
    expression(bold("a) ") * bolditalic("Halictus scabiosae") * bold(" guarding occurrence")),
    xlab = NULL, ylab = "Probability"
  ) /
  plot_guarding_panel(
    duration_curve, "middle_scaled",
    expression(bold("b) ") * bolditalic("Halictus scabiosae") * bold(" guarding duration")),
    xlab = "Soil temperature (scaled)", ylab = "Log(duration)"
  )

guarding_figure

ggsave("HS_guarding_figure_newly_updated.png", guarding_figure,
       width = 200, height = 160, units = "mm", dpi = 600)


# Species abundance and environment ----

# Linear model of log(abundance + 1) with a species-by-environment
# interaction for every predictor. Standardisation is done after reshaping
# to long format, as in the reported analysis.
abundance_predictors <- c("ambient_scaled", "shallow_scaled", "middle_scaled",
                          "deep_scaled", "light_scaled")

abundance_long <- focal %>%
  add_microclimate() %>%
  select(Avg_total_number_of_HS, Avg_total_number_of_LM,
         ambient_temp, shallow_temp, middle_temp, deep_temp,
         Avg_lux_light_logger, Site, Day) %>%
  pivot_longer(cols = c(Avg_total_number_of_HS, Avg_total_number_of_LM),
               names_to = "species", values_to = "abundance") %>%
  mutate(
    species        = factor(ifelse(species == "Avg_total_number_of_HS", "HS", "LM"),
                            levels = c("HS", "LM")),
    abundance_log  = log(abundance + 1),
    ambient_scaled = scale1(ambient_temp),
    shallow_scaled = scale1(shallow_temp),
    middle_scaled  = scale1(middle_temp),
    deep_scaled    = scale1(deep_temp),
    light_scaled   = scale1(Avg_lux_light_logger),
    Site = factor(Site),
    Day  = factor(Day)
  ) %>%
  na.omit()

abundance_model <- lm(
  abundance_log ~ species *
    (ambient_scaled + shallow_scaled + middle_scaled + deep_scaled + light_scaled),
  data = abundance_long
)

summary(abundance_model)

# Marginal F test for each species-by-environment interaction (each term
# dropped from the full model). anova() is not used here because its
# sequential tests depend on the order of the correlated predictors.
drop1(abundance_model, test = "F")


# Abundance: species-specific slopes ----

# Slope of each environmental variable within each species, from the
# interaction model
for (v in abundance_predictors) {
  cat("\nSlope of", v, "by species\n")
  print(as.data.frame(emtrends(abundance_model, ~ species, var = v, infer = TRUE)))
}
summary(abundance_model)$r.squared


# Abundance: sensitivity analysis with random effects ----

# abundance_model is a plain linear model, so the two species rows from one
# observation, and repeated observations within a day, are treated as
# independent. Same fixed effects, plus random intercepts for day and
# observation.
rob_wide <- focal %>%
  add_microclimate() %>%
  mutate(
    obs_id         = factor(row_number()),
    Site           = factor(Site),
    Day            = factor(Day),
    ambient_scaled = scale1(ambient_temp),
    shallow_scaled = scale1(shallow_temp),
    middle_scaled  = scale1(middle_temp),
    deep_scaled    = scale1(deep_temp),
    light_scaled   = scale1(Avg_lux_light_logger)
  ) %>%
  filter(!is.na(Avg_lux_light_logger))

rob_long <- rob_wide %>%
  select(obs_id, Day, Site, ambient_scaled, shallow_scaled, middle_scaled,
         deep_scaled, light_scaled,
         HS = Avg_total_number_of_HS, LM = Avg_total_number_of_LM) %>%
  pivot_longer(c(HS, LM), names_to = "species", values_to = "abundance") %>%
  mutate(species = factor(species, levels = c("HS", "LM")),
         abundance_log = log(abundance + 1)) %>%
  na.omit()

m_abund_re <- lmer(
  abundance_log ~ species *
    (ambient_scaled + shallow_scaled + middle_scaled + deep_scaled + light_scaled) +
    (1 | Day) + (1 | obs_id),
  data = rob_long, REML = FALSE
)
summary(m_abund_re)
drop1(m_abund_re, test = "Chisq")

for (v in abundance_predictors) {
  cat("\nSlope of", v, "by species (random-effects model)\n")
  print(as.data.frame(emtrends(m_abund_re, ~ species, var = v,
                               lmer.df = "satterthwaite", infer = TRUE)))
}
# If signs and significance agree with the linear model, its results are not
# an artefact of ignoring day and observation.


# Species abundance figure ----

abundance_curve <- function(variable) {
  grid <- expand.grid(
    x       = seq(-2, 2, length.out = 200),
    species = factor(c("HS", "LM"), levels = c("HS", "LM"))
  )
  names(grid)[1] <- variable
  
  # Hold every other predictor at its mean (0 on the standardised scale)
  for (v in setdiff(abundance_predictors, variable)) grid[[v]] <- 0
  
  p <- predict(abundance_model, newdata = grid, se.fit = TRUE)
  grid %>%
    mutate(fit       = p$fit,
           conf.low  = p$fit - 1.96 * p$se.fit,
           conf.high = p$fit + 1.96 * p$se.fit)
}

plot_abundance_panel <- function(variable, title, xlab) {
  ggplot(
    abundance_curve(variable),
    aes(
      x = .data[[variable]],
      colour = species,
      fill = species
    )
  ) +
    geom_ribbon(
      aes(ymin = conf.low, ymax = conf.high),
      alpha = 0.2,
      colour = NA
    ) +
    geom_line(aes(y = fit), linewidth = 1) +
    scale_colour_manual(
      name = "Species",
      values = species_cols,
      labels = species_labels
    ) +
    scale_fill_manual(
      name = "Species",
      values = species_cols,
      labels = species_labels
    ) +
    labs(
      title = title,
      x = xlab,
      y = "Log(abundance + 1)"
    ) +
    theme_ms(14, axis_face = "plain") +
    theme(
      plot.title = element_text(size = 20, face = "bold"),
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.text.x = element_text(size = 15, face = "bold"),
      axis.text.y = element_text(size = 15, face = "bold"),
      legend.title = element_text(size = 22, face = "bold"),
      legend.text = element_text(size = 20, face = "bold")
    )
}

abundance_figure <-
  plot_abundance_panel(
    "ambient_scaled", "a) Ambient air",
    "Ambient temperature (scaled)"
  ) +
  plot_abundance_panel(
    "shallow_scaled", "b) Shallow soil",
    "2–3 cm soil temperature (scaled)"
  ) +
  plot_abundance_panel(
    "middle_scaled", "c) Middle soil",
    "5 cm soil temperature (scaled)"
  ) +
  plot_abundance_panel(
    "deep_scaled", "d) Deep soil",
    "10 cm soil temperature (scaled)"
  ) +
  plot_abundance_panel(
    "light_scaled", "e) Light intensity",
    "Light intensity (scaled)"
  ) +
  guide_area() +
  plot_layout(ncol = 3, guides = "collect")

abundance_figure

ggsave("multi_panel_abundance_clean_newly_updated.png", abundance_figure,
       width = 360, height = 240, units = "mm", dpi = 600)


# Summary statistics for the abundance figure caption ----

# Mean and SD of each soil layer across all focal observations. The plotted
# range of -2 to +2 SD corresponds to mean +/- 2 SD (not the observed minimum
# and maximum, which are given here for reference).
focal %>%
  add_microclimate() %>%
  summarise(across(c(shallow_temp, middle_temp, deep_temp),
                   list(mean = mean, sd = sd, min = min, max = max),
                   .names = "{.col}_{.fn}")) %>%
  pivot_longer(everything())


# Nest traffic: GAMs by species and direction ----

# One GAM per species and direction, with a smooth for every layer.
# Layers are standardised on all rows, then rows with complete counts are kept.
movement_data <- focal %>%
  add_microclimate() %>%
  add_scaled_env() %>%
  filter(!is.na(Incoming_frequency_HS), !is.na(Outgoing_frequency_HS),
         !is.na(Incoming_frequency_LM), !is.na(Outgoing_frequency_LM))

movement_gam_data <- movement_data %>%
  select(Incoming_frequency_HS, Outgoing_frequency_HS,
         Incoming_frequency_LM, Outgoing_frequency_LM,
         ambient_s, surface_s, shallow_s, middle_s, deep_s, light_s) %>%
  na.omit()

fit_movement_gam <- function(response) {
  gam(as.formula(paste0(
    "log(", response, " + 1) ~ s(ambient_s) + s(surface_s) + s(shallow_s) + ",
    "s(middle_s) + s(deep_s) + s(light_s)")),
    data = movement_gam_data, method = "REML")
}

movement_models <- list(
  HS_in  = fit_movement_gam("Incoming_frequency_HS"),
  HS_out = fit_movement_gam("Outgoing_frequency_HS"),
  LM_in  = fit_movement_gam("Incoming_frequency_LM"),
  LM_out = fit_movement_gam("Outgoing_frequency_LM")
)

# Smooth-term statistics and deviance explained (basis of the GAM table)
movement_gam_table <- bind_rows(lapply(names(movement_models), function(nm) {
  s <- summary(movement_models[[nm]])
  as.data.frame(s$s.table) %>%
    rownames_to_column("smooth") %>%
    rename(p_value = `p-value`) %>%
    mutate(model = nm, deviance_explained = 100 * s$dev.expl, .before = 1)
}))
movement_gam_table

write.csv(movement_gam_table, "movement_gam_summary.csv", row.names = FALSE)


# Nest traffic figures ----

# Effect of one layer with all other layers held at their mean
predict_smooth <- function(model, variable) {
  grid <- data.frame(
    ambient_s = 0, surface_s = 0, shallow_s = 0,
    middle_s  = 0, deep_s    = 0, light_s   = 0
  )[rep(1, 200), ]
  grid[[variable]] <- seq(-2, 2, length.out = 200)
  
  p <- predict(model, newdata = grid, se.fit = TRUE)
  data.frame(x = grid[[variable]], fit = p$fit,
             lower = p$fit - 1.96 * p$se.fit,
             upper = p$fit + 1.96 * p$se.fit)
}

plot_movement_panel <- function(curve, colour, title, xlab) {
  ggplot(curve, aes(x, fit)) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      fill = colour,
      alpha = 0.25
    ) +
    geom_line(
      colour = colour,
      linewidth = 1.2
    ) +
    labs(
      title = title,
      x = xlab,
      y = "Log(frequency + 1)"
    ) +
    theme_ms(13, axis_face = "plain") +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 17, face = "bold"),
      axis.title.y = element_text(size = 17, face = "bold"),
      axis.text.x = element_text(size = 15, face = "bold"),
      axis.text.y = element_text(size = 15, face = "bold")
    )
}

# Four panels per layer: species (rows) by direction (columns)
movement_figure <- function(variable, title, xlab) {
  panels <- list(
    plot_movement_panel(
      predict_smooth(movement_models$HS_in, variable),
      col_HS, "HS Incoming", xlab
    ),
    plot_movement_panel(
      predict_smooth(movement_models$HS_out, variable),
      col_HS, "HS Outgoing", xlab
    ),
    plot_movement_panel(
      predict_smooth(movement_models$LM_in, variable),
      col_LM, "LM Incoming", xlab
    ),
    plot_movement_panel(
      predict_smooth(movement_models$LM_out, variable),
      col_LM, "LM Outgoing", xlab
    )
  )
  
  wrap_plots(panels, ncol = 2) +
    plot_annotation(
      title = title,
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 23)
      )
    )
}

movement_figures <- lapply(seq_len(nrow(env_layers)), function(i) {
  movement_figure(
    env_layers$variable[i],
    env_layers$movement_title[i],
    env_layers$xlab[i]
  )
})
names(movement_figures) <- env_layers$movement_title

for (nm in names(movement_figures)) {
  ggsave(
    paste0("Fig_", gsub(" ", "_", nm), ".png"),
    movement_figures[[nm]],
    width = 200,
    height = 160,
    units = "mm",
    dpi = 600,
    bg = "white"
  )
}

movement_figures[[1]]


# Species comparison across layers (supplementary) ----

# Single-layer GAMs for each species, drawn on shared axes.
# Layers come from the same averaged loggers as every other chapter.
plot_species_layer <- function(curves, title, xlab, ylab) {
  ggplot(curves, aes(x, fit, colour = species, fill = species)) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      alpha = 0.15,
      colour = NA
    ) +
    geom_line(linewidth = 1.2) +
    scale_colour_manual(
      name = "Species",
      values = species_cols,
      labels = species_labels
    ) +
    scale_fill_manual(
      name = "Species",
      values = species_cols,
      labels = species_labels
    ) +
    labs(
      title = title,
      x = xlab,
      y = ylab
    ) +
    theme_ms(14) +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 19, face = "bold"),
      axis.title.y = element_text(size = 19, face = "bold"),
      axis.text.x = element_text(size = 17, face = "bold"),
      axis.text.y = element_text(size = 17, face = "bold"),
      legend.title = element_text(size = 22, face = "bold"),
      legend.text = element_text(size = 20, face = "bold")
    )
}

species_gam_figure <- function(response_HS, response_LM, y_label,
                               log_response = FALSE, lettered = FALSE,
                               show_legend = FALSE) {
  
  d <- focal %>%
    add_microclimate() %>%
    select(
      all_of(c(HS = response_HS, LM = response_LM)),
      ambient_temp, surface_temp, shallow_temp, middle_temp, deep_temp,
      Avg_lux_light_logger
    ) %>%
    na.omit() %>%
    add_scaled_env()
  
  if (log_response) {
    d <- mutate(d, HS = log(HS + 1), LM = log(LM + 1))
  }
  
  curves_for <- function(layer_var) {
    grid <- data.frame(
      seq(
        min(d[[layer_var]]),
        max(d[[layer_var]]),
        length.out = 100
      )
    )
    names(grid) <- layer_var
    
    bind_rows(lapply(c("HS", "LM"), function(sp) {
      m <- gam(
        as.formula(paste0(sp, " ~ s(", layer_var, ")")),
        data = d,
        method = "REML"
      )
      p <- predict(m, newdata = grid, se.fit = TRUE)
      
      data.frame(
        species = factor(sp, levels = c("HS", "LM")),
        x = grid[[layer_var]],
        fit = p$fit,
        lower = p$fit - 1.96 * p$se.fit,
        upper = p$fit + 1.96 * p$se.fit
      )
    }))
  }
  
  panels <- lapply(seq_len(nrow(env_layers)), function(i) {
    title <- env_layers$title[i]
    if (lettered) {
      title <- paste0(letters[i], ") ", title)
    }
    
    plot_species_layer(
      curves_for(env_layers$variable[i]),
      title,
      env_layers$xlab[i],
      y_label
    )
  })
  
  fig <- wrap_plots(panels, ncol = 3)
  
  if (show_legend) {
    fig + plot_layout(guides = "collect") &
      theme(legend.position = "top")
  } else {
    fig & theme(legend.position = "none")
  }
}

incoming_figure <- species_gam_figure(
  "Incoming_frequency_HS",
  "Incoming_frequency_LM",
  y_label = "Predicted activity frequency"
)
incoming_figure

ggsave("FINAL_GAM_HS_LM_ALL_LAYERS.png", incoming_figure,
       width = 300, height = 200, units = "mm", dpi = 600)

abundance_gam_figure <- species_gam_figure(
  "Avg_total_number_of_HS", "Avg_total_number_of_LM",
  y_label = "Log(abundance + 1)",
  log_response = TRUE, lettered = TRUE, show_legend = TRUE
)
abundance_gam_figure

ggsave("FINAL_GAM_ABUNDANCE_ONLY.png", abundance_gam_figure,
       width = 300, height = 200, units = "mm", dpi = 600)


# Timing of activity: data preparation ----

# Uses the merged focal data frame from the chapters above.
# Logger key: Thermochron_1 ambient; _3/_6 shallow soil (2-3 cm);
# _7/_8 middle soil (5 cm). Site loggers are averaged.
timing <- focal %>%
  mutate(
    obs_id  = factor(row_number()),          # one row = one focal observation
    Site    = factor(Site),
    Day     = factor(Day),
    hour    = as.numeric(Time),
    ambient = Thermochron_1,
    shallow = rowMeans(cbind(Thermochron_3, Thermochron_6), na.rm = TRUE),
    middle  = rowMeans(cbind(Thermochron_7, Thermochron_8), na.rm = TRUE)
  )

# If Time was stored as HHMM (e.g. 900, 1400), convert it to hours
if (max(timing$hour, na.rm = TRUE) > 24) timing$hour <- timing$hour / 100

timing <- timing %>%
  mutate(
    hour_c    = hour - mean(hour, na.rm = TRUE),
    ambient_s = scale1(ambient),
    shallow_s = scale1(shallow),
    middle_s  = scale1(middle)
  )

# One row per species and observation
timing_long <- timing %>%
  transmute(
    obs_id, Day, Site, hour, hour_c, ambient_s, shallow_s, middle_s,
    HS_abund = Avg_total_number_of_HS,
    LM_abund = Avg_total_number_of_LM,
    HS_inc   = Incoming_frequency_HS,
    LM_inc   = Incoming_frequency_LM,
    HS_out   = Outgoing_frequency_HS,
    LM_out   = Outgoing_frequency_LM
  ) %>%
  pivot_longer(
    cols = matches("^(HS|LM)_", ignore.case = FALSE),
    names_to = c("species", ".value"),
    names_sep = "_"
  ) %>%
  mutate(
    species = factor(species, levels = c("HS", "LM")),
    y_abund = log(abund + 1)
  )


# Timing of activity: sampling design ----

# Hours sampled and how many observations fall in each
table(timing$hour, useNA = "ifany")
n_hours <- n_distinct(timing$hour)
n_hours

# Hours covered on each day and site
timing %>%
  group_by(Day, Site) %>%
  summarise(n_obs = n(), first_hour = min(hour), last_hour = max(hour),
            .groups = "drop")

# How entangled are time of day and temperature?
cor(timing[, c("hour", "ambient", "shallow", "middle")],
    use = "pairwise.complete.obs")


# Timing of activity: time versus temperature ----

# Candidate models for abundance:
#   null      species + site only
#   time      species-specific effect of hour
#   temp      species-specific effects of air and shallow soil temperature
#   both      time and temperature together
#   additive  same as both, without species interactions (for collinearity)
# HS and LM come from the same observation, so observation is a random effect;
# site has two levels and is a fixed effect.
fit_timing_set <- function(data, response) {
  d <- data %>%
    filter(!is.na(.data[[response]]), !is.na(hour_c),
           !is.na(ambient_s), !is.na(shallow_s)) %>%
    mutate(y = .data[[response]])
  
  fit <- function(rhs) {
    lmer(as.formula(paste("y ~", rhs, "+ (1 | Day) + (1 | obs_id)")),
         data = d, REML = FALSE)
  }
  
  list(
    data = d,
    models = list(
      null     = fit("species + Site"),
      time     = fit("species * hour_c + Site"),
      temp     = fit("species * (ambient_s + shallow_s) + Site"),
      both     = fit("species * (hour_c + ambient_s + shallow_s) + Site"),
      additive = fit("species + hour_c + ambient_s + shallow_s + Site")
    )
  )
}

report_timing_set <- function(res) {
  cat("\nSingular fits:\n")
  print(sapply(res$models, isSingular))
  
  cat("\nModel comparison (AICc):\n")
  print(model.sel(res$models[c("null", "time", "temp", "both")], rank = "AICc"))
  
  cat("\nCollinearity (additive model):\n")
  print(check_collinearity(res$models$additive))
  
  # Species interactions in the full model:
  #   species:hour_c    species differ in timing beyond temperature
  #   species:ambient_s species differ in their air-temperature response
  #   species:shallow_s species differ in their soil-temperature response
  cat("\nMarginal tests of the species interactions:\n")
  print(drop1(res$models$both, test = "Chisq"))
  
  cat("\nFull model:\n")
  print(summary(res$models$both))
}


# Timing of activity: abundance ----

abund <- fit_timing_set(timing_long, "y_abund")
report_timing_set(abund)


# Timing of activity: species-specific slopes ----

# Slope of each predictor within each species, from one of the fitted models
species_slopes <- function(res, terms = c("hour_c", "ambient_s", "shallow_s"),
                           model = "both") {
  bind_rows(lapply(terms, function(v) {
    as.data.frame(emtrends(res$models[[model]], ~ species, var = v,
                           lmer.df = "satterthwaite", infer = TRUE)) %>%
      rename(trend = 2) %>%
      mutate(term = v, .before = 1)
  }))
}

species_slopes(abund)


# Timing of activity: adding light ----

# The light logger has gaps, so models that include light are fitted on the
# observations where light exists and are compared with each other on those
# same rows. Where are the gaps?
timing %>%
  group_by(Day, Site) %>%
  summarise(n = n(), missing_light = sum(is.na(Avg_lux_light_logger)),
            .groups = "drop") %>%
  print(n = Inf)

# How does light relate to hour and temperature?
timing %>%
  transmute(hour, ambient, shallow, middle, light = Avg_lux_light_logger) %>%
  cor(use = "complete.obs")

light_long <- timing_long %>%
  left_join(timing %>% transmute(obs_id, light_s = scale1(Avg_lux_light_logger)),
            by = "obs_id") %>%
  filter(!is.na(light_s))

fit_light_set <- function(data, response) {
  d <- data %>%
    filter(!is.na(.data[[response]])) %>%
    mutate(y = .data[[response]])
  
  fit <- function(rhs) {
    lmer(as.formula(paste("y ~", rhs, "+ (1 | Day) + (1 | obs_id)")),
         data = d, REML = FALSE)
  }
  
  list(
    data = d,
    models = list(
      null       = fit("species + Site"),
      time       = fit("species * hour_c + Site"),
      light      = fit("species * light_s + Site"),
      temp       = fit("species * (ambient_s + shallow_s) + Site"),
      temp_light = fit("species * (ambient_s + shallow_s + light_s) + Site"),
      time_temp  = fit("species * (hour_c + ambient_s + shallow_s) + Site"),
      all        = fit("species * (hour_c + ambient_s + shallow_s + light_s) + Site"),
      additive   = fit("species + hour_c + ambient_s + shallow_s + light_s + Site")
    )
  )
}

report_light_set <- function(res) {
  cat("\nObservations used:", n_distinct(res$data$obs_id), "\n")
  
  cat("\nModel comparison (AICc):\n")
  print(model.sel(res$models[c("null", "time", "light", "temp", "temp_light",
                               "time_temp", "all")], rank = "AICc"))
  
  cat("\nCollinearity (additive model):\n")
  print(check_collinearity(res$models$additive))
  
  cat("\nMarginal tests of the species interactions (full model):\n")
  print(drop1(res$models$all))
  
  cat("\nSpecies-specific slopes (full model):\n")
  print(species_slopes(res, terms = c("hour_c", "ambient_s", "shallow_s", "light_s"),
                       model = "all"))
}


# Timing of activity: abundance with light ----

light_abund <- fit_light_set(light_long, "y_abund")
report_light_set(light_abund)


# Timing of activity: nest entries and exits as counts ----

# Entries and exits are whole-number counts, so they are modelled with negative
# binomial GLMMs (glmmTMB) using the same candidate models and random effects
# as abundance. Both statements must print TRUE.
for (v in c("inc", "out")) {
  cat(v, "all whole numbers:",
      all(timing_long[[v]] == round(timing_long[[v]]), na.rm = TRUE), "\n")
}
cat("Share of zero counts: entries",
    round(mean(timing_long$inc == 0, na.rm = TRUE), 2), " exits",
    round(mean(timing_long$out == 0, na.rm = TRUE), 2), "\n")

nb_formulas <- list(
  full = c(
    null = "species + Site",
    time = "species * hour_c + Site",
    temp = "species * (ambient_s + shallow_s) + Site",
    both = "species * (hour_c + ambient_s + shallow_s) + Site"
  ),
  light = c(
    null       = "species + Site",
    time       = "species * hour_c + Site",
    light      = "species * light_s + Site",
    temp       = "species * (ambient_s + shallow_s) + Site",
    temp_light = "species * (ambient_s + shallow_s + light_s) + Site",
    time_temp  = "species * (hour_c + ambient_s + shallow_s) + Site",
    all        = "species * (hour_c + ambient_s + shallow_s + light_s) + Site"
  )
)

fit_nb <- function(data, response, formulas) {
  d <- data %>%
    filter(!is.na(.data[[response]])) %>%
    mutate(y = .data[[response]])
  lapply(formulas, function(rhs) {
    glmmTMB(as.formula(paste("y ~", rhs, "+ (1 | Day) + (1 | obs_id)")),
            family = nbinom2, data = d)
  })
}

# Models that did not converge to a positive-definite Hessian are reported and
# excluded from the AICc comparison.
report_nb <- function(models, label, top_model, slope_terms) {
  cat("\n\n==========", label, "==========\n")
  
  converged <- sapply(models, function(m) isTRUE(m$sdr$pdHess))
  cat("\nExcluded (no convergence):",
      if (all(converged)) "none" else paste(names(models)[!converged], collapse = ", "),
      "\n")
  
  cat("\nModel comparison (AICc):\n")
  print(model.sel(models[converged], rank = "AICc"))
  
  m <- models[[top_model]]
  cat("\nFull model (species interactions are the species differences):\n")
  print(summary(m))
  
  cat("\nSpecies-specific slopes (log scale):\n")
  for (v in slope_terms) {
    cat("\n", v, "\n")
    print(as.data.frame(emtrends(m, ~ species, var = v, infer = TRUE)))
  }
}

nb_results <- list()
for (resp in c("inc", "out")) {
  nb_results[[paste0(resp, "_163")]] <- fit_nb(timing_long, resp, nb_formulas$full)
  report_nb(nb_results[[paste0(resp, "_163")]], paste0(resp, "_163 observations"),
            "both", c("hour_c", "ambient_s", "shallow_s"))
  
  nb_results[[paste0(resp, "_95")]] <- fit_nb(light_long, resp, nb_formulas$light)
  report_nb(nb_results[[paste0(resp, "_95")]], paste0(resp, "_95 observations with light"),
            "all", c("hour_c", "ambient_s", "shallow_s", "light_s"))
}
# inc = nest entries, out = nest exits


# Timing of activity: daily activity curve ----

# Species-specific smooth of hour, with random intercepts for day and
# observation. k is capped by the number of distinct hours sampled.
k_h <- max(3, min(5, n_hours - 1))

gam_d <- timing_long %>% filter(!is.na(y_abund))

gam_shared <- gam(
  y_abund ~ species + Site + s(hour, k = k_h) +
    s(Day, bs = "re") + s(obs_id, bs = "re"),
  data = gam_d, method = "ML"
)

gam_by_species <- gam(
  y_abund ~ species + Site + s(hour, by = species, k = k_h) +
    s(Day, bs = "re") + s(obs_id, bs = "re"),
  data = gam_d, method = "ML"
)

# Do the species have different daily curves?
anova(gam_shared, gam_by_species, test = "Chisq")

gam_final <- update(gam_by_species, method = "REML")
summary(gam_final)

# Predictions with the random effects excluded (Site held at its first level)
curve_grid <- expand.grid(
  hour    = seq(min(gam_d$hour), max(gam_d$hour), length.out = 100),
  species = c("HS", "LM")
) %>%
  mutate(
    species = factor(species, levels = c("HS", "LM")),
    Site    = factor(levels(gam_d$Site)[1],   levels = levels(gam_d$Site)),
    Day     = factor(levels(gam_d$Day)[1],    levels = levels(gam_d$Day)),
    obs_id  = factor(levels(gam_d$obs_id)[1], levels = levels(gam_d$obs_id))
  )

pr <- predict(gam_final, newdata = curve_grid, se.fit = TRUE,
              exclude = c("s(Day)", "s(obs_id)"))

curve_grid <- curve_grid %>%
  mutate(fit   = pr$fit,
         lower = pr$fit - 1.96 * pr$se.fit,
         upper = pr$fit + 1.96 * pr$se.fit)

# Hour of peak predicted activity per species (descriptive, no interval)
curve_grid %>%
  group_by(species) %>%
  slice_max(fit, n = 1) %>%
  select(species, hour, fit)


# Timing of activity: centre of activity ----

# One value per day, site and species. The first active hour is limited by
# when observations began each day, so treat it as supporting evidence only.
daily <- timing_long %>%
  filter(!is.na(abund)) %>%
  group_by(Day, Site, species) %>%
  summarise(
    n_hours      = n_distinct(hour),
    total        = sum(abund),
    centre       = ifelse(total > 0, sum(hour * abund) / total, NA_real_),
    first_active = suppressWarnings(min(hour[abund > 0])),
    .groups = "drop"
  ) %>%
  mutate(first_active = ifelse(is.finite(first_active), first_active, NA_real_))

# Days with both species observed in at least three hours
paired_days <- daily %>%
  group_by(Day, Site) %>%
  filter(n() == 2, all(n_hours >= 3)) %>%
  ungroup()

nrow(paired_days) / 2   # number of day-site pairs

m_centre <- lmer(centre ~ species + Site + (1 | Day), data = paired_days)
m_first  <- lmer(first_active ~ species + Site + (1 | Day), data = paired_days)

summary(m_centre)
summary(m_first)


# Timing of activity: figures ----

activity_by_hour <- ggplot() +
  geom_ribbon(
    data = curve_grid,
    aes(hour, ymin = lower, ymax = upper, fill = species),
    alpha = 0.2
  ) +
  geom_line(
    data = curve_grid,
    aes(hour, fit, colour = species),
    linewidth = 1.2
  ) +
  stat_summary(
    data = gam_d,
    aes(hour, y_abund, colour = species),
    fun.data = mean_se,
    geom = "pointrange",
    position = position_dodge(width = 0.15),
    size = 0.5
  ) +
  scale_colour_manual(
    name = "Species",
    values = species_cols,
    labels = species_labels,
    guide = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  scale_fill_manual(
    values = species_cols,
    labels = species_labels,
    guide = "none"
  ) +
  scale_x_continuous(breaks = sort(unique(timing$hour))) +
  labs(
    x = "Hour of day",
    y = "Log(abundance + 1)"
  ) +
  theme_ms(15) +
  theme(
    axis.title.x = element_text(size = 19, face = "bold"),
    axis.title.y = element_text(size = 19, face = "bold"),
    axis.text.x = element_text(size = 17, face = "bold"),
    axis.text.y = element_text(size = 17, face = "bold"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 15, face = "bold"),
    legend.key.width = grid::unit(0.7, "cm"),
    legend.key.height = grid::unit(0.4, "cm"),
    legend.spacing.x = grid::unit(0.2, "cm"),
    legend.margin = margin(0, 0, 0, 0)
  )

activity_by_hour

ggsave(
  "focal_activity_by_hour.png",
  activity_by_hour,
  width = 190,
  height = 125,
  units = "mm",
  dpi = 600
)

temperature_by_hour <- ggplot(timing, aes(hour, ambient)) +
  geom_jitter(
    width = 0.12,
    height = 0,
    size = 2.4,
    alpha = 0.55,
    colour = "#0072B2"
  ) +
  geom_smooth(
    method = "loess",
    formula = y ~ x,
    se = TRUE,
    colour = "#004C75",
    fill = "#0072B2",
    alpha = 0.2,
    linewidth = 1.3
  ) +
  scale_x_continuous(breaks = sort(unique(timing$hour))) +
  labs(
    x = "Hour of day",
    y = "Ambient temperature (°C)"
  ) +
  theme_ms(15) +
  theme(
    axis.title.x = element_text(size = 19, face = "bold"),
    axis.title.y = element_text(size = 19, face = "bold"),
    axis.text.x = element_text(size = 17, face = "bold"),
    axis.text.y = element_text(size = 17, face = "bold")
  )

temperature_by_hour

ggsave("ambient_temperature_by_hour.png", temperature_by_hour,
       width = 160, height = 120, units = "mm", dpi = 600)


# Sample sizes used by each analysis ----

cat("Guarding duration bouts:           ", nrow(m_duration@frame), "\n")
cat("Guarding occurrence rows:          ", nrow(m_occurrence@frame), "\n")
cat("Abundance linear model (obs):      ", nrow(abundance_long) / 2, "\n")
cat("Models with light (obs):           ", n_distinct(light_long$obs_id), "\n")
cat("Nest traffic GAMs (obs):           ", nrow(movement_gam_data), "\n")
cat("Timing models without light (obs): ", n_distinct(timing_long$obs_id), "\n")


# Model diagnostics: residuals ----

# Simulation-based residual checks (DHARMa): KS test of the residual
# distribution, a dispersion test and an outlier test. Plots are saved for the
# supplement. Gaussian GAMs and models that DHARMa cannot simulate are
# reported with NA rather than stopping the script.
check_residuals <- function(model, name) {
  tryCatch({
    sim <- simulateResiduals(model, n = 1000, plot = FALSE)
    stats <- data.frame(
      analysis     = "focal",
      model        = name,
      ks_p         = testUniformity(sim, plot = FALSE)$p.value,
      dispersion_p = testDispersion(sim, plot = FALSE)$p.value,
      outlier_p    = testOutliers(sim, plot = FALSE)$p.value
    )
    
    png(file.path("diagnostics", paste0("focal_", name, ".png")),
        width = 1800, height = 900, res = 200)
    tryCatch(plot(sim), finally = invisible(dev.off()))
    stats
  }, error = function(e) {
    cat("DHARMa failed for", name, ":", conditionMessage(e), "\n")
    data.frame(analysis = "focal", model = name,
               ks_p = NA, dispersion_p = NA, outlier_p = NA)
  })
}

dir.create("diagnostics", showWarnings = FALSE)

diagnostic_models <- list(
  guarding_occurrence       = m_occurrence,
  guarding_duration         = m_duration,
  abundance_lm              = abundance_model,
  abundance_lmer_RE         = m_abund_re,
  timing_abundance          = abund$models$both,
  light_abundance           = light_abund$models$all,
  count_entries_163obs      = nb_results$inc_163$both,
  count_entries_95obs_light = nb_results$inc_95$all,
  count_exits_163obs        = nb_results$out_163$both,
  count_exits_95obs_light   = nb_results$out_95$all,
  gam_HS_in                 = movement_models$HS_in,
  gam_HS_out                = movement_models$HS_out,
  gam_LM_in                 = movement_models$LM_in,
  gam_LM_out                = movement_models$LM_out,
  gam_daily_curve           = gam_final,
  centre_of_activity        = m_centre,
  first_active_hour         = m_first
)

focal_diagnostics <- bind_rows(lapply(names(diagnostic_models), function(nm) {
  check_residuals(diagnostic_models[[nm]], nm)
}))
focal_diagnostics

write.csv(focal_diagnostics, "diagnostics/focal_dharma_summary.csv", row.names = FALSE)

# Zero inflation in the count models
lapply(list(entries_163 = nb_results$inc_163$both, entries_95 = nb_results$inc_95$all,
            exits_163   = nb_results$out_163$both, exits_95   = nb_results$out_95$all),
       function(m) testZeroInflation(simulateResiduals(m, n = 1000, plot = FALSE),
                                     plot = FALSE)$p.value)


# Model diagnostics: GAM basis dimension and concurvity ----

# k-index close to 1 (and a non-small p) indicates an adequate basis dimension;
# concurvity is the GAM analogue of collinearity (values near 1 are problematic).
lapply(movement_models, k.check)
k.check(gam_final)
lapply(movement_models, function(m) concurvity(m, full = TRUE)["worst", ])

