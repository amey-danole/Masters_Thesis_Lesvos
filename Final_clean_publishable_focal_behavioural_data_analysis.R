# ============================================================
# FOCAL BEHAVIOURAL DATA ANALYSIS 
# ============================================================

# ============================================================
# 1. LOAD ESSENTIAL PACKAGES----
# ============================================================

library(dplyr)
library(tidyr)
library(lubridate)
library(readr)
library(lme4)
library(lmerTest)
library(MuMIn)
library(ggplot2)
library(patchwork)
library(cowplot)

# ============================================================
# 2. DATA IMPORT----
# ============================================================

# Raw behavioural and environmental datasets

focal_behavioural_data <- read.csv("focal_behavioural_data.csv")
combined_env <- read.csv("combined_environmental_variables.csv")

# ============================================================
# 3. DATA CLEANING AND MERGING----
# ============================================================

# Align temporal structure and merge datasets by Date and Time

focal_behavioural_data <- focal_behavioural_data %>%
  mutate(Date = as.character(Date),
         Time = as.numeric(Time))

combined_env <- combined_env %>%
  mutate(Date = as.character(Date),
         Time = as.numeric(Time))

Combined_focal_behavioural_dataset <- focal_behavioural_data %>%
  inner_join(combined_env, by = c("Date", "Time"))

# Exploring dataset

str(Combined_focal_behavioural_dataset)

focal_species_summary <- Combined_focal_behavioural_dataset %>%
  summarise(
    
    # --- HS (Halictus scabiosae) ---
    HS_mean_activity   = mean(Avg_total_number_of_HS, na.rm = TRUE),
    HS_guarding_time   = mean(Avg_duration_of_nest_guarding_HS, na.rm = TRUE),
    HS_incoming_freq   = mean(Incoming_frequency_HS, na.rm = TRUE),
    HS_outgoing_freq   = mean(Outgoing_frequency_HS, na.rm = TRUE),
    
    # --- LM (Lasioglossum malachurum) ---
    LM_mean_activity   = mean(Avg_total_number_of_LM, na.rm = TRUE),
    LM_guarding_time   = NA,  # not applicable
    LM_incoming_freq   = mean(Incoming_frequency_LM, na.rm = TRUE),
    LM_outgoing_freq   = mean(Outgoing_frequency_LM, na.rm = TRUE)
  )

focal_species_summary

# ============================================================
# 4. Halictus scabiosae — OCCURRENCE MODEL
# ============================================================

# Binary GLMM: presence/absence of nest guarding behaviour

HS_occ <- Combined_focal_behavioural_dataset %>%
  mutate(
    guard_present = ifelse(is.na(Avg_duration_of_nest_guarding_HS), 0, 1),
    
    # shallow soil temperature (2–3 cm)
    shallow_temp = rowMeans(cbind(Thermochron_3, Thermochron_6), na.rm = TRUE),
    
    shallow_scaled = as.numeric(scale(shallow_temp))
  ) %>%
  filter(!is.na(Site), !is.na(Day)) %>%
  mutate(
    Site = factor(Site),
    Day = factor(Day)
  )

m_HS_occ <- glmer(
  guard_present ~ shallow_scaled + (1 | Site) + (1 | Day),
  family = binomial,
  data = HS_occ
)

summary(m_HS_occ)

# ============================================================
# 5. Halictus scabiosae — DURATION MODEL----
# ============================================================

# LMM: log-transformed nest guarding duration

HS_guard <- Combined_focal_behavioural_dataset %>%
  mutate(
    middle_temp = rowMeans(cbind(Thermochron_7, Thermochron_8), na.rm = TRUE),
    middle_scaled = as.numeric(scale(middle_temp))
  ) %>%
  filter(!is.na(Avg_duration_of_nest_guarding_HS),
         Avg_duration_of_nest_guarding_HS > 0) %>%
  mutate(
    log_guard = log(Avg_duration_of_nest_guarding_HS),
    Site = factor(Site),
    Day = factor(Day)
  )

m_HS_final <- lmer(
  log_guard ~ middle_scaled + (1 | Site) + (1 | Day),
  data = HS_guard,
  REML = TRUE
)

summary(m_HS_final)

# ============================================================
# 6. Halictus scabiosae — COMPOSITE FIGURE
# ============================================================

# -----------------------------
# OCCURRENCE (2–3 cm soil)
# -----------------------------

newdat_occ <- data.frame(
  shallow_scaled = seq(-2, 2, length.out = 200)
)

pred_occ <- predict(
  m_HS_occ,
  newdata = newdat_occ,
  type = "link",
  se.fit = TRUE,
  re.form = NA
)

plot_occ <- newdat_occ %>%
  mutate(
    fit = plogis(pred_occ$fit),
    se = pred_occ$se.fit,
    conf.low = plogis(pred_occ$fit - 1.96 * se),
    conf.high = plogis(pred_occ$fit + 1.96 * se)
  )

p1 <- ggplot(plot_occ, aes(x = shallow_scaled)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              fill = "#E69F00", alpha = 0.25) +
  geom_line(aes(y = fit),
            colour = "#E69F00", linewidth = 1.2) +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text  = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold")
  ) +
  labs(
    title = expression("a) " * italic("Halictus scabiosae") * " occurrence"),
    x = NULL,
    y = "Probability"
  )

# -----------------------------
# DURATION (5 cm soil)
# -----------------------------

newdat_dur <- data.frame(
  middle_scaled = seq(-2, 2, length.out = 200)
)

pred_dur <- predict(
  m_HS_final,
  newdata = newdat_dur,
  se.fit = TRUE,
  re.form = NA
)

plot_dur <- newdat_dur %>%
  mutate(
    fit = pred_dur$fit,
    se = pred_dur$se.fit,
    conf.low = fit - 1.96 * se,
    conf.high = fit + 1.96 * se
  )

p2 <- ggplot(plot_dur, aes(x = middle_scaled)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              fill = "#E69F00", alpha = 0.25) +
  geom_line(aes(y = fit),
            colour = "#E69F00", linewidth = 1.2) +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text  = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold")
  ) +
  labs(
    title = expression("b) " * italic("Halictus scabiosae") * " duration"),
    x = "Soil temperature (scaled)",
    y = "Log(duration)"
  )

# -----------------------------
# COMBINE
# -----------------------------

HS_plot <- p1 / p2

ggsave(
  "HS_guarding_figure_newly_updated.png",
  HS_plot,
  width = 180,
  height = 160,
  units = "mm",
  dpi = 600
)

# ============================================================
# 7. MULTI-SPECIES ABUNDANCE–ENVIRONMENT MODELLING----
# ============================================================

# Joint LMM: HS + LM abundance responses across environmental gradients

bee_long <- Combined_focal_behavioural_dataset %>%
  mutate(
    ambient_temp = Thermochron_1,
    shallow_temp = rowMeans(cbind(Thermochron_3, Thermochron_6), na.rm = TRUE),
    middle_temp  = rowMeans(cbind(Thermochron_7, Thermochron_8), na.rm = TRUE),
    deep_temp    = rowMeans(cbind(Thermochron_9, Thermochron_10), na.rm = TRUE)
  ) %>%
  select(Avg_total_number_of_HS,
         Avg_total_number_of_LM,
         ambient_temp, shallow_temp, middle_temp, deep_temp,
         Avg_lux_light_logger,
         Site, Day) %>%
  pivot_longer(cols = c(Avg_total_number_of_HS, Avg_total_number_of_LM),
               names_to = "species",
               values_to = "abundance") %>%
  mutate(
    species = dplyr::recode(
      species,
      "Avg_total_number_of_HS" = "HS",
      "Avg_total_number_of_LM" = "LM"
    ),
    abundance_log = log(abundance + 1),
    ambient_scaled = as.numeric(scale(ambient_temp)),
    shallow_scaled = as.numeric(scale(shallow_temp)),
    middle_scaled  = as.numeric(scale(middle_temp)),
    deep_scaled    = as.numeric(scale(deep_temp)),
    light_scaled   = as.numeric(scale(Avg_lux_light_logger)),
    Site = factor(Site),
    Day = factor(Day)
  ) %>%
  na.omit()

m_bee_global_fixed <- lm(
  abundance_log ~ species *
    (ambient_scaled + shallow_scaled + middle_scaled + deep_scaled + light_scaled),
  data = bee_long
)

summary(m_bee_global_fixed)

cols <- c("HS" = "#E69F00", "LM" = "#0072B2")

species_labels <- c(
  "HS" = expression(italic("Halictus scabiosae")),
  "LM" = expression(italic("Lasioglossum malachurum"))
)

# ============================================================
# MODEL PREDICTION FUNCTION 
# ============================================================

make_panel <- function(var) {
  
  predictors <- c(
    "ambient_scaled",
    "shallow_scaled",
    "middle_scaled",
    "deep_scaled",
    "light_scaled"
  )
  
  newdat <- expand.grid(
    x = seq(-2, 2, length.out = 200),
    species = c("HS", "LM")
  )
  
  names(newdat)[1] <- var
  
  for (p in predictors) {
    if (!p %in% names(newdat)) {
      newdat[[p]] <- 0
    }
  }
  
  pred <- predict(
    m_bee_global,
    newdata = newdat,
    re.form = NA,
    se.fit = TRUE
  )
  
  newdat %>%
    mutate(
      fit = pred$fit,
      se = pred$se.fit,
      conf.low = fit - 1.96 * se,
      conf.high = fit + 1.96 * se
    )
}

# ============================================================
# PANEL PLOTTING FUNCTION
# ============================================================

plot_panel <- function(df, var, title, xlabel) {
  
  ggplot(df, aes(x = .data[[var]], colour = species, fill = species)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
                alpha = 0.2, colour = NA) +
    geom_line(aes(y = fit), linewidth = 1) +
    scale_colour_manual(values = cols, labels = species_labels) +
    scale_fill_manual(values = cols, labels = species_labels) +
    guides(colour = guide_legend(title = "Species", ncol = 1),
           fill = guide_legend(title = "Species", ncol = 1)) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "top",
      axis.title = element_text(face = "plain"),
      axis.text = element_text(face = "plain"),
      plot.title = element_text(face = "bold")
    ) +
    labs(title = title,
         x = xlabel,
         y = "Log(abundance + 1)")
}

# ============================================================
# FIGURE GENERATION (ENVIRONMENTAL GRADIENTS)
# ============================================================

p1 <- plot_panel(make_panel("ambient_scaled"), "ambient_scaled",
                 "a) Ambient air", "Ambient temperature")

p2 <- plot_panel(make_panel("shallow_scaled"), "shallow_scaled",
                 "b) Shallow soil", "2–3 cm soil temperature (scaled)")

p3 <- plot_panel(make_panel("middle_scaled"), "middle_scaled",
                 "c) Middle soil", "5 cm soil temperature (scaled)")

p4 <- plot_panel(make_panel("deep_scaled"), "deep_scaled",
                 "d) Deep soil", "10 cm soil temperature (scaled)")

p5 <- plot_panel(make_panel("light_scaled"), "light_scaled",
                 "e) Light intensity", "Light intensity (scaled)")

# ============================================================
# FIGURE ASSEMBLY AND EXPORT
# ============================================================

p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")
p4 <- p4 + theme(legend.position = "none")
p5 <- p5 + theme(legend.position = "none")

legend_plot <- ggplot(bee_long, aes(x = ambient_scaled, y = abundance_log, colour = species)) +
  geom_point(size = 3) +
  scale_colour_manual(values = cols, labels = species_labels) +
  guides(colour = guide_legend(ncol = 1,
                               override.aes = list(size = 6))) +
  theme_classic(base_size = 14) +
  theme(legend.position = "right") +
  labs(colour = "Species")

legend <- cowplot::get_legend(legend_plot)
legend_panel <- wrap_elements(full = legend)

final_plot <- (p1 | p2 | p3) /
  (p4 | p5 | legend_panel)

final_plot

ggsave("multi_panel_abundance_clean_newly_updated.png",
       final_plot,
       width = 360, height = 240,
       units = "mm", dpi = 600)

# =====================================================================
# 8. GAM RIDGE PLOTS — HALICTUS SCABIOSAE & LASIOGLOSSUM MALACHURUM----
# Environmental drivers of movement frequency (incoming/outgoing)
# =====================================================================

library(dplyr)
library(ggplot2)
library(mgcv)
library(patchwork)

# ============================================================
# COLOUR PALETTE (species identity)
# ============================================================
col_HS <- "#E69F00"
col_LM <- "#0072B2"

# ============================================================
# DATA PREPARATION
# ============================================================

# Environmental variables are averaged across loggers and scaled
dat <- Combined_focal_behavioural_dataset %>%
  
  mutate(
    ambient_temp = Thermochron_1,
    surface_temp = rowMeans(cbind(Thermochron_2, Thermochron_4), na.rm = TRUE),
    shallow_temp = rowMeans(cbind(Thermochron_3, Thermochron_6), na.rm = TRUE),
    middle_temp  = rowMeans(cbind(Thermochron_7, Thermochron_8), na.rm = TRUE),
    deep_temp    = rowMeans(cbind(Thermochron_9, Thermochron_10), na.rm = TRUE)
  ) %>%
  
  mutate(
    ambient_s = as.numeric(scale(ambient_temp)),
    surface_s = as.numeric(scale(surface_temp)),
    shallow_s = as.numeric(scale(shallow_temp)),
    middle_s  = as.numeric(scale(middle_temp)),
    deep_s    = as.numeric(scale(deep_temp)),
    light_s   = as.numeric(scale(Avg_lux_light_logger))
  ) %>%
  
  filter(
    !is.na(Incoming_frequency_HS),
    !is.na(Outgoing_frequency_HS),
    !is.na(Incoming_frequency_LM),
    !is.na(Outgoing_frequency_LM)
  )

# Reduced dataset for GAM fitting (complete cases only)
dat_gam <- dat %>%
  dplyr::select(
    Incoming_frequency_HS,
    Outgoing_frequency_HS,
    Incoming_frequency_LM,
    Outgoing_frequency_LM,
    ambient_s, surface_s, shallow_s,
    middle_s, deep_s, light_s
  ) %>%
  na.omit()

# ============================================================
# GAM MODELS (species × movement direction)
# ============================================================

gam_HS_in <- gam(log(Incoming_frequency_HS + 1) ~
                   s(ambient_s) + s(surface_s) + s(shallow_s) +
                   s(middle_s) + s(deep_s) + s(light_s),
                 data = dat_gam, method = "REML")

gam_HS_out <- gam(log(Outgoing_frequency_HS + 1) ~
                    s(ambient_s) + s(surface_s) + s(shallow_s) +
                    s(middle_s) + s(deep_s) + s(light_s),
                  data = dat_gam, method = "REML")

gam_LM_in <- gam(log(Incoming_frequency_LM + 1) ~
                   s(ambient_s) + s(surface_s) + s(shallow_s) +
                   s(middle_s) + s(deep_s) + s(light_s),
                 data = dat_gam, method = "REML")

gam_LM_out <- gam(log(Outgoing_frequency_LM + 1) ~
                    s(ambient_s) + s(surface_s) + s(shallow_s) +
                    s(middle_s) + s(deep_s) + s(light_s),
                  data = dat_gam, method = "REML")
summary(gam_HS_in)
summary(gam_HS_out)
summary(gam_LM_in)
summary(gam_LM_out)

# ============================================================
# PREDICTION FUNCTION 
# ============================================================
predict_one <- function(model, varname){
  
  seq_vals <- seq(-2, 2, length.out = 200)
  
  newdata <- data.frame(
    ambient_s = 0,
    surface_s = 0,
    shallow_s = 0,
    middle_s  = 0,
    deep_s    = 0,
    light_s   = 0
  )
  
  newdata <- newdata[rep(1, length(seq_vals)), ]
  newdata[[varname]] <- seq_vals
  
  pr <- predict(model, newdata = newdata, se.fit = TRUE)
  
  data.frame(
    x = seq_vals,
    fit = pr$fit,
    upper = pr$fit + 1.96 * pr$se.fit,
    lower = pr$fit - 1.96 * pr$se.fit
  )
}

# ============================================================
# PLOT TEMPLATE
# ============================================================

plot_panel <- function(df, col, title, xlab = "Temperature (scaled)"){
  
  ggplot(df, aes(x = x, y = fit)) +
    
    geom_ribbon(aes(ymin = lower, ymax = upper),
                fill = col, alpha = 0.25) +
    
    geom_line(color = col, linewidth = 1.2) +
    
    theme_classic(base_size = 13) +
    
    labs(
      title = title,
      x = xlab,
      y = "Log(frequency + 1)"
    ) +
    
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

# ============================================================
# MASTER FUNCTION (4-panel species × direction layout)
# ============================================================

make_env_figure <- function(varname, xlabel){
  
  HS_in  <- predict_one(gam_HS_in,  varname)
  HS_out <- predict_one(gam_HS_out, varname)
  LM_in  <- predict_one(gam_LM_in,  varname)
  LM_out <- predict_one(gam_LM_out, varname)
  
  p1 <- plot_panel(HS_in,  col_HS, "HS Incoming")
  p2 <- plot_panel(HS_out, col_HS, "HS Outgoing")
  p3 <- plot_panel(LM_in,  col_LM, "LM Incoming")
  p4 <- plot_panel(LM_out, col_LM, "LM Outgoing")
  
  (p1 | p2) / (p3 | p4) +
    plot_annotation(
      title = xlabel,
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
      )
    )
}

# ============================================================
# ENVIRONMENTAL VARIABLES (loop)
# ============================================================
vars <- c(
  "ambient_s",
  "surface_s",
  "shallow_s",
  "middle_s",
  "deep_s",
  "light_s"
)

labels <- c(
  "Ambient temperature",
  "Surface soil temperature",
  "Shallow soil temperature",
  "Middle soil temperature",
  "Deep soil temperature",
  "Light intensity"
)

# ============================================================
# GENERATE ALL RIDGE PLOTS
# ============================================================
fig_list <- lapply(seq_along(vars), function(i){
  make_env_figure(vars[i], labels[i])
})

names(fig_list) <- labels

# ============================================================
# SAVE ALL FIGURES
# ============================================================
for(i in seq_along(fig_list)){
  ggsave(
    filename = paste0("Fig_", gsub(" ", "_", names(fig_list)[i]), ".png"),
    plot = fig_list[[i]],
    width = 180,
    height = 140,
    units = "mm",
    dpi = 600
  )
}

# ============================================================
# OPTIONAL DISPLAY
# ============================================================
fig_list[[1]]

# =======================================================================
# LIGHT FIGURE ONLY (Halictus scabiosae + Lasioglossum malachurum IN/OUT)
# =======================================================================

light_HS_in  <- predict_one(gam_HS_in,  "light_s")
light_HS_out <- predict_one(gam_HS_out, "light_s")
light_LM_in  <- predict_one(gam_LM_in,  "light_s")
light_LM_out <- predict_one(gam_LM_out, "light_s")

p1 <- plot_panel(light_HS_in,  col_HS, "HS Incoming", "Light intensity (scaled)")
p2 <- plot_panel(light_HS_out, col_HS, "HS Outgoing", "Light intensity (scaled)")
p3 <- plot_panel(light_LM_in,  col_LM, "LM Incoming", "Light intensity (scaled)")
p4 <- plot_panel(light_LM_out, col_LM, "LM Outgoing", "Light intensity (scaled)")

fig_light <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title = "Light intensity",
    theme = theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 14,
        face = "bold"
      )
    )
  )

fig_light

ggsave(
  filename = "Fig_Light_intensity.png",
  plot = fig_light,
  width = 180,
  height = 140,
  units = "mm",
  dpi = 600,
  bg = "white"
)

# ================================================================
# 9. FINAL GAM RIDGE FIGURE — HS vs LM across environmental layers
# ================================================================

library(mgcv)
library(dplyr)
library(ggplot2)
library(patchwork)

# ============================================================
# COLOUR PALETTE (species identity)
# ============================================================
col_HS <- "#E69F00"
col_LM <- "#0072B2"

# ============================================================
# DATA PREPARATION (complete-case GAM dataset)
# ============================================================
dat <- Combined_focal_behavioural_dataset %>%
  dplyr::select(
    Incoming_frequency_HS,
    Incoming_frequency_LM,
    Avg_T_ambient,
    Thermochron_1,   # surface (0–2 cm)
    Thermochron_2,   # shallow (2–3 cm)
    Thermochron_4,   # middle (5 cm)
    Thermochron_10,  # deep (10 cm)
    Avg_lux_light_logger
  ) %>%
  na.omit()

# ============================================================
# SCALE ENVIRONMENTAL PREDICTORS
# ============================================================
dat <- dat %>%
  mutate(
    amb_s   = as.numeric(scale(Avg_T_ambient)),
    surf_s  = as.numeric(scale(Thermochron_1)),
    shal_s  = as.numeric(scale(Thermochron_2)),
    mid_s   = as.numeric(scale(Thermochron_4)),
    deep_s  = as.numeric(scale(Thermochron_10)),
    light_s = as.numeric(scale(Avg_lux_light_logger))
  )

# ============================================================
# GAM MODELS (HS and LM incoming activity only)
# ============================================================
gam_HS_amb   <- gam(Incoming_frequency_HS ~ s(amb_s),   data = dat, method = "REML")
gam_HS_surf  <- gam(Incoming_frequency_HS ~ s(surf_s),  data = dat, method = "REML")
gam_HS_shal  <- gam(Incoming_frequency_HS ~ s(shal_s),  data = dat, method = "REML")
gam_HS_mid   <- gam(Incoming_frequency_HS ~ s(mid_s),   data = dat, method = "REML")
gam_HS_deep  <- gam(Incoming_frequency_HS ~ s(deep_s),  data = dat, method = "REML")
gam_HS_light <- gam(Incoming_frequency_HS ~ s(light_s), data = dat, method = "REML")

gam_LM_amb   <- gam(Incoming_frequency_LM ~ s(amb_s),   data = dat, method = "REML")
gam_LM_surf  <- gam(Incoming_frequency_LM ~ s(surf_s),  data = dat, method = "REML")
gam_LM_shal  <- gam(Incoming_frequency_LM ~ s(shal_s),  data = dat, method = "REML")
gam_LM_mid   <- gam(Incoming_frequency_LM ~ s(mid_s),   data = dat, method = "REML")
gam_LM_deep  <- gam(Incoming_frequency_LM ~ s(deep_s),  data = dat, method = "REML")
gam_LM_light <- gam(Incoming_frequency_LM ~ s(light_s), data = dat, method = "REML")

# ============================================================
# PREDICTION FUNCTION (95% CI)
# ============================================================
predict_gam <- function(model, newdata, var){
  
  pr <- predict(model, newdata = newdata, se.fit = TRUE)
  
  data.frame(
    x = newdata[[var]],
    fit = pr$fit,
    upper = pr$fit + 1.96 * pr$se.fit,
    lower = pr$fit - 1.96 * pr$se.fit
  )
}

# ============================================================
# GRID GENERATION
# ============================================================
grid_fun <- function(v){
  data.frame(seq(min(dat[[v]]), max(dat[[v]]), length.out = 100))
}

grid_amb   <- grid_fun("amb_s");   names(grid_amb)   <- "amb_s"
grid_surf  <- grid_fun("surf_s");  names(grid_surf)  <- "surf_s"
grid_shal  <- grid_fun("shal_s");  names(grid_shal)  <- "shal_s"
grid_mid   <- grid_fun("mid_s");   names(grid_mid)   <- "mid_s"
grid_deep  <- grid_fun("deep_s");  names(grid_deep)  <- "deep_s"
grid_light <- grid_fun("light_s"); names(grid_light) <- "light_s"

# ============================================================
# PREDICTIONS
# ============================================================
HS_amb   <- predict_gam(gam_HS_amb,   grid_amb,   "amb_s")
HS_surf  <- predict_gam(gam_HS_surf,  grid_surf,  "surf_s")
HS_shal  <- predict_gam(gam_HS_shal,  grid_shal,  "shal_s")
HS_mid   <- predict_gam(gam_HS_mid,   grid_mid,   "mid_s")
HS_deep  <- predict_gam(gam_HS_deep,  grid_deep,  "deep_s")
HS_light <- predict_gam(gam_HS_light, grid_light, "light_s")

LM_amb   <- predict_gam(gam_LM_amb,   grid_amb,   "amb_s")
LM_surf  <- predict_gam(gam_LM_surf,  grid_surf,  "surf_s")
LM_shal  <- predict_gam(gam_LM_shal,  grid_shal,  "shal_s")
LM_mid   <- predict_gam(gam_LM_mid,   grid_mid,   "mid_s")
LM_deep  <- predict_gam(gam_LM_deep,  grid_deep,  "deep_s")
LM_light <- predict_gam(gam_LM_light, grid_light, "light_s")

# ============================================================
# PLOT FUNCTION (with implicit legend support)
# ============================================================
plot_layer <- function(HS, LM, title, xlab){
  
  ggplot() +
    
    geom_ribbon(data = HS,
                aes(x = x, ymin = lower, ymax = upper),
                fill = col_HS, alpha = 0.15) +
    geom_line(data = HS,
              aes(x = x, y = fit),
              color = col_HS, linewidth = 1.2) +
    
    geom_ribbon(data = LM,
                aes(x = x, ymin = lower, ymax = upper),
                fill = col_LM, alpha = 0.15) +
    geom_line(data = LM,
              aes(x = x, y = fit),
              color = col_LM, linewidth = 1.2) +
    
    theme_classic(base_size = 14) +
    
    labs(
      title = title,
      x = xlab,
      y = "Predicted activity frequency"
    ) +
    
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.position = "none"   # no fake legend needed (clean panel layout)
    )
}

# ============================================================
# FINAL 6-PANEL FIGURE
# ============================================================
p1 <- plot_layer(HS_amb,  LM_amb,  "Ambient air", "Temperature (scaled)")
p2 <- plot_layer(HS_surf, LM_surf, "Surface soil", "Temperature (scaled)")
p3 <- plot_layer(HS_shal, LM_shal, "Shallow soil (2–3 cm)", "Temperature (scaled)")
p4 <- plot_layer(HS_mid,  LM_mid,  "Middle soil (5 cm)", "Temperature (scaled)")
p5 <- plot_layer(HS_deep, LM_deep, "Deep soil (10 cm)", "Temperature (scaled)")
p6 <- plot_layer(HS_light, LM_light, "Light intensity", "Light intensity (scaled)")

final_gam_fig <- (p1 | p2 | p3) / (p4 | p5 | p6)

final_gam_fig

# ============================================================
# EXPORT FIGURE
# ============================================================
ggsave(
  "FINAL_GAM_HS_LM_ALL_LAYERS.png",
  final_gam_fig,
  width = 300,
  height = 200,
  units = "mm",
  dpi = 600
)

# =========================================================
# FINAL GAM RIDGE FIGURE FOR HS and LM ABUNDANCE COMPARISON
# =========================================================

# =========================================================
# CLEAN DATA 
# =========================================================
dat <- Combined_focal_behavioural_dataset %>%
  dplyr::select(
    Avg_total_number_of_HS,
    Avg_total_number_of_LM,
    Avg_T_ambient,
    Thermochron_1,
    Thermochron_2,
    Thermochron_4,
    Thermochron_10,
    Avg_lux_light_logger
  ) %>%
  na.omit()

# =========================================================
# SCALE ENVIRONMENTAL VARIABLES
# =========================================================
dat <- dat %>%
  mutate(
    amb_s   = as.numeric(scale(Avg_T_ambient)),
    surf_s  = as.numeric(scale(Thermochron_1)),
    shal_s  = as.numeric(scale(Thermochron_2)),
    mid_s   = as.numeric(scale(Thermochron_4)),
    deep_s  = as.numeric(scale(Thermochron_10)),
    light_s = as.numeric(scale(Avg_lux_light_logger))
  )

# =========================================================
# GAM MODELS (ABUNDANCE ONLY)
# =========================================================
gam_HS_amb  <- gam(log(Avg_total_number_of_HS + 1) ~ s(amb_s),  data = dat, method = "REML")
gam_HS_surf <- gam(log(Avg_total_number_of_HS + 1) ~ s(surf_s), data = dat, method = "REML")
gam_HS_shal <- gam(log(Avg_total_number_of_HS + 1) ~ s(shal_s), data = dat, method = "REML")
gam_HS_mid  <- gam(log(Avg_total_number_of_HS + 1) ~ s(mid_s),  data = dat, method = "REML")
gam_HS_deep <- gam(log(Avg_total_number_of_HS + 1) ~ s(deep_s), data = dat, method = "REML")
gam_HS_light<- gam(log(Avg_total_number_of_HS + 1) ~ s(light_s),data = dat, method = "REML")

gam_LM_amb  <- gam(log(Avg_total_number_of_LM + 1) ~ s(amb_s),  data = dat, method = "REML")
gam_LM_surf <- gam(log(Avg_total_number_of_LM + 1) ~ s(surf_s), data = dat, method = "REML")
gam_LM_shal <- gam(log(Avg_total_number_of_LM + 1) ~ s(shal_s), data = dat, method = "REML")
gam_LM_mid  <- gam(log(Avg_total_number_of_LM + 1) ~ s(mid_s),  data = dat, method = "REML")
gam_LM_deep <- gam(log(Avg_total_number_of_LM + 1) ~ s(deep_s), data = dat, method = "REML")
gam_LM_light<- gam(log(Avg_total_number_of_LM + 1) ~ s(light_s),data = dat, method = "REML")

# =========================================================
# PREDICTION FUNCTION
# =========================================================
predict_gam <- function(model, newdata, var){
  
  pr <- predict(model, newdata = newdata, se.fit = TRUE)
  
  data.frame(
    x = newdata[[var]],
    fit = pr$fit,
    upper = pr$fit + 1.96 * pr$se.fit,
    lower = pr$fit - 1.96 * pr$se.fit
  )
}

# =========================================================
# GRIDS
# =========================================================
grid_fun <- function(v){
  data.frame(seq(min(dat[[v]]), max(dat[[v]]), length.out = 100))
}

grid_amb   <- grid_fun("amb_s");   names(grid_amb)   <- "amb_s"
grid_surf  <- grid_fun("surf_s");  names(grid_surf)  <- "surf_s"
grid_shal  <- grid_fun("shal_s");  names(grid_shal)  <- "shal_s"
grid_mid   <- grid_fun("mid_s");   names(grid_mid)   <- "mid_s"
grid_deep  <- grid_fun("deep_s");  names(grid_deep)  <- "deep_s"
grid_light <- grid_fun("light_s"); names(grid_light) <- "light_s"

# =========================================================
# PREDICTIONS
# =========================================================
HS_amb   <- predict_gam(gam_HS_amb,  grid_amb,   "amb_s")
HS_surf  <- predict_gam(gam_HS_surf, grid_surf,  "surf_s")
HS_shal  <- predict_gam(gam_HS_shal, grid_shal,  "shal_s")
HS_mid   <- predict_gam(gam_HS_mid,  grid_mid,   "mid_s")
HS_deep  <- predict_gam(gam_HS_deep, grid_deep,  "deep_s")
HS_light <- predict_gam(gam_HS_light,grid_light, "light_s")

LM_amb   <- predict_gam(gam_LM_amb,  grid_amb,   "amb_s")
LM_surf  <- predict_gam(gam_LM_surf, grid_surf,  "surf_s")
LM_shal  <- predict_gam(gam_LM_shal, grid_shal,  "shal_s")
LM_mid   <- predict_gam(gam_LM_mid,  grid_mid,   "mid_s")
LM_deep  <- predict_gam(gam_LM_deep, grid_deep,  "deep_s")
LM_light <- predict_gam(gam_LM_light,grid_light, "light_s")

# =========================================================
# PLOT FUNCTION
# =========================================================
plot_layer <- function(HS, LM, title, xlab){
  
  ggplot() +
    
    geom_ribbon(data = HS,
                aes(x = x, ymin = lower, ymax = upper, fill = "HS"),
                alpha = 0.15) +
    geom_line(data = HS,
              aes(x = x, y = fit, colour = "HS"),
              linewidth = 1.2) +
    
    geom_ribbon(data = LM,
                aes(x = x, ymin = lower, ymax = upper, fill = "LM"),
                alpha = 0.15) +
    geom_line(data = LM,
              aes(x = x, y = fit, colour = "LM"),
              linewidth = 1.2) +
    
    scale_colour_manual(
      name = "Species",
      values = c("HS" = col_HS, "LM" = col_LM),
      labels = c(
        "HS" = expression(italic("Halictus scabiosae")),
        "LM" = expression(italic("Lasioglossum malachurum"))
      )
    ) +
    
    scale_fill_manual(
      name = "Species",
      values = c("HS" = col_HS, "LM" = col_LM),
      labels = c(
        "HS" = expression(italic("Halictus scabiosae")),
        "LM" = expression(italic("Lasioglossum malachurum"))
      )
    ) +
    
    theme_classic(base_size = 14) +
    
    labs(
      x = xlab,
      y = "Log(abundance + 1)",
      title = title
    ) +
    
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.position = "none",
      legend.title = element_text(face = "bold")
    )
}

# =========================================================
# FINAL FIGURE (UNCHANGED STRUCTURE)
# =========================================================
p1 <- plot_layer(HS_amb,  LM_amb,  "a) Ambient temperature", "Temperature (scaled)")
p2 <- plot_layer(HS_surf, LM_surf, "b) Surface soil", "Temperature (scaled)")
p3 <- plot_layer(HS_shal, LM_shal, "c) Shallow soil (2–3 cm)", "Temperature (scaled)")
p4 <- plot_layer(HS_mid,  LM_mid,  "d) Middle soil (5 cm)", "Temperature (scaled)")
p5 <- plot_layer(HS_deep, LM_deep, "e) Deep soil (10 cm)", "Temperature (scaled)")
p6 <- plot_layer(HS_light,LM_light,"f) Light intensity", "Light intensity (scaled)")

# extract legend from one plot
legend <- cowplot::get_legend(
  p1 + theme(legend.position = "top")
)

legend_row <- wrap_elements(full = legend)

final_gam_fig <- legend_row / ((p1 | p2 | p3) / (p4 | p5 | p6)) +
  plot_layout(heights = c(0.08, 1))

final_gam_fig

# =========================================================
# SAVE
# =========================================================
ggsave(
  "FINAL_GAM_ABUNDANCE_ONLY.png",
  final_gam_fig,
  width = 300,
  height = 200,
  units = "mm",
  dpi = 600
)

