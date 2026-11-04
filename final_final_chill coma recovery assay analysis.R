# Chill-coma recovery: thermal thresholds, recovery times and body size
# Halictus scabiosae (H.S) and Lasioglossum malachurum (L.M)
#
# Questions
#   1. Do the species differ in body size (intertegular distance, ITD)?
#   2. At what body temperature do bees right themselves and take off after
#      chill coma, and does this depend on ambient temperature, species or size?
#   3. How long does recovery take, and does it depend on the same factors?
#   4. Are the temperature effects explained by recovery duration?
#
# Input (working directory)
#   Chill_coma_recovery.csv
#
# Outputs
#   chill_coma_model_table.csv   fixed effects of every model
#   new_ITD_violin.png, Multipanel_thermal_performance.png,
#   Thermal_traits_vs_ITD_multipanel.png
#   diagnostics/chillcoma_dharma_summary.csv and residual plots


# Packages ----

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(lme4)
library(lmerTest)
library(car)
library(DHARMa)


# Shared settings ----

species_colors <- c("H.S" = "#E69F00", "L.M" = "#0072B2")

species_labels <- c(
  "H.S" = expression(italic("H. scabiosae")),
  "L.M" = expression(italic("L. malachurum"))
)

theme_pub <- theme_classic(base_size = 14) +
  theme(
    axis.title      = element_text(face = "bold"),
    axis.text       = element_text(size = 12),
    legend.position = "none"
  )

dir.create("diagnostics", showWarnings = FALSE)


# Load and prepare data ----

full_data <- read.csv("Chill_coma_recovery.csv")

# Durations are recorded as mm:ss; convert to seconds
mmss_to_seconds <- function(x) {
  sapply(strsplit(as.character(x), ":"),
         function(p) as.numeric(p[1]) * 60 + as.numeric(p[2]))
}

full_data <- full_data %>%
  mutate(
    Bee_species      = factor(Bee_species, levels = c("H.S", "L.M")),
    Date             = factor(Date),
    Righting_seconds = mmss_to_seconds(Righting_Duration),
    Take_off_seconds = mmss_to_seconds(Take_off_Duration),
    log_righting     = log(Righting_seconds + 1),
    log_takeoff      = log(Take_off_seconds + 1)
  )

# Sample sizes by species, and missing values per measurement
table(full_data$Bee_species)
full_data %>%
  summarise(across(c(ITD_Final, Bee_weight, Righting_Temperature,
                     Take_off_Temperature, Righting_seconds, Take_off_seconds,
                     Righting_T_ambient, Take_off_T_ambient),
                   ~ sum(is.na(.x)))) %>%
  print()


# Body size ----

# ITD by species and the size ratio
itd_summary <- full_data %>%
  group_by(Bee_species) %>%
  summarise(n = sum(!is.na(ITD_Final)),
            mean_ITD = mean(ITD_Final, na.rm = TRUE),
            sd_ITD   = sd(ITD_Final, na.rm = TRUE),
            .groups = "drop")
itd_summary

itd_ratio <- itd_summary$mean_ITD[itd_summary$Bee_species == "H.S"] /
  itd_summary$mean_ITD[itd_summary$Bee_species == "L.M"]
itd_ratio   # fold difference, H. scabiosae / L. malachurum

# Test for a species difference in ITD
itd_model <- lm(ITD_Final ~ Bee_species, data = full_data)
summary(itd_model)

# Body mass and ITD are collinear; ITD is retained as the size measure
cor(full_data$ITD_Final, full_data$Bee_weight, use = "complete.obs")


# Thermal thresholds: righting and take-off temperature ----

# Final models contain species and ambient temperature. ITD is not included
# because it is strongly confounded with species (H. scabiosae is 1.65-fold
# larger); models with ITD are fitted for comparison.
righting_temp_final <- lmer(
  Righting_Temperature ~ Bee_species + Righting_T_ambient + (1 | Date),
  data = full_data
)
summary(righting_temp_final)

righting_temp_with_itd <- lmer(
  Righting_Temperature ~ Bee_species + ITD_Final + Righting_T_ambient + (1 | Date),
  data = full_data
)
summary(righting_temp_with_itd)

takeoff_temp_final <- lmer(
  Take_off_Temperature ~ Bee_species + Take_off_T_ambient + (1 | Date),
  data = full_data
)
summary(takeoff_temp_final)

takeoff_temp_with_itd <- lmer(
  Take_off_Temperature ~ Bee_species + ITD_Final + Take_off_T_ambient + (1 | Date),
  data = full_data
)
summary(takeoff_temp_with_itd)

# Collinearity between species, ITD and ambient temperature
vif(righting_temp_with_itd)
vif(takeoff_temp_with_itd)

# Does body size explain take-off temperature within each species?
takeoff_temp_by_species <- lapply(
  split(full_data, full_data$Bee_species),
  function(d) lmer(Take_off_Temperature ~ ITD_Final + Take_off_T_ambient + (1 | Date),
                   data = d)
)
lapply(takeoff_temp_by_species, summary)


# Recovery times: time to right and time to take-off ----

# Durations are log(seconds + 1). The full models contain species, ITD and
# ambient temperature; the reduced models drop species, which is collinear
# with ITD.
righting_duration <- lmer(
  log_righting ~ Bee_species + ITD_Final + Righting_T_ambient + (1 | Date),
  data = full_data
)
summary(righting_duration)
vif(righting_duration)

righting_duration_no_species <- lmer(
  log_righting ~ ITD_Final + Righting_T_ambient + (1 | Date),
  data = full_data
)
summary(righting_duration_no_species)

takeoff_duration <- lmer(
  log_takeoff ~ Bee_species + ITD_Final + Take_off_T_ambient + (1 | Date),
  data = full_data
)
summary(takeoff_duration)

takeoff_duration_no_species <- lmer(
  log_takeoff ~ ITD_Final + Take_off_T_ambient + (1 | Date),
  data = full_data
)
summary(takeoff_duration_no_species)


# Do recovery times explain the temperature effects? ----

# Righting and take-off temperature after accounting for the time taken
righting_temp_with_duration <- lmer(
  Righting_Temperature ~ Bee_species + ITD_Final + Righting_T_ambient +
    log_righting + (1 | Date),
  data = full_data
)
summary(righting_temp_with_duration)

takeoff_temp_with_duration <- lmer(
  Take_off_Temperature ~ Bee_species + ITD_Final + Take_off_T_ambient +
    log_takeoff + (1 | Date),
  data = full_data
)
summary(takeoff_temp_with_duration)


# Does righting time vary with ambient temperature? ----

summary(full_data$Righting_seconds)
sd(full_data$Righting_seconds, na.rm = TRUE)
cor(full_data$Righting_T_ambient, full_data$Righting_seconds, use = "complete.obs")

# Species-specific temperature effect on righting time (likelihood-ratio test;
# fitted by maximum likelihood because the fixed effects differ)
full_data$ambient_c <- as.numeric(scale(full_data$Righting_T_ambient, scale = FALSE))

duration_additive <- lmer(
  log_righting ~ Bee_species + ITD_Final + ambient_c + (1 | Date),
  data = full_data, REML = FALSE
)
duration_interaction <- lmer(
  log_righting ~ Bee_species * ambient_c + ITD_Final + (1 | Date),
  data = full_data, REML = FALSE
)
anova(duration_additive, duration_interaction)
summary(duration_interaction)

ggplot(full_data, aes(Righting_T_ambient, Righting_seconds, colour = Bee_species)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  scale_colour_manual(values = species_colors, labels = species_labels,
                      name = "Species") +
  labs(x = "Ambient temperature (°C)", y = "Righting duration (s)") +
  theme_classic()


# Does body size modify the ambient-temperature effect on take-off? ----

# Nested comparison (maximum likelihood): take-off temperature with and
# without an ITD x ambient interaction
takeoff_temp_additive_ml <- update(takeoff_temp_with_itd, REML = FALSE)
takeoff_temp_interaction <- lmer(
  Take_off_Temperature ~ Bee_species + ITD_Final * Take_off_T_ambient + (1 | Date),
  data = full_data, REML = FALSE
)
anova(takeoff_temp_additive_ml, takeoff_temp_interaction)
summary(takeoff_temp_interaction)


# Influential observations ----

# Cook's distance for the righting-temperature model, and the model refitted
# without observations above the 4/n threshold
influence_check <- tryCatch(influence(righting_temp_final, obs = TRUE),
                            error = function(e) NULL)

if (!is.null(influence_check)) {
  cooks <- cooks.distance(influence_check)
  plot(cooks)
  influential <- which(cooks > 4 / nobs(righting_temp_final))
  print(full_data[influential, ])
  
  righting_temp_no_influential <- update(righting_temp_final,
                                         data = full_data[-influential, ])
  print(summary(righting_temp_no_influential))
} else {
  cat("influence() is not available for this model class; skipped.\n")
}


# Model table ----

chill_models <- list(
  righting_temp_final          = righting_temp_final,
  righting_temp_with_itd       = righting_temp_with_itd,
  takeoff_temp_final           = takeoff_temp_final,
  takeoff_temp_with_itd        = takeoff_temp_with_itd,
  takeoff_temp_HS              = takeoff_temp_by_species[["H.S"]],
  takeoff_temp_LM              = takeoff_temp_by_species[["L.M"]],
  righting_duration            = righting_duration,
  righting_duration_no_species = righting_duration_no_species,
  takeoff_duration             = takeoff_duration,
  takeoff_duration_no_species  = takeoff_duration_no_species,
  righting_temp_with_duration  = righting_temp_with_duration,
  takeoff_temp_with_duration   = takeoff_temp_with_duration
)

# Fixed effects of every model (basis of the model-statistics table)
model_table <- bind_rows(lapply(names(chill_models), function(nm) {
  m  <- chill_models[[nm]]
  cf <- as.data.frame(coef(summary(m)))
  names(cf)[1:5] <- c("estimate", "se", "df", "t", "p")
  cf$term  <- rownames(cf)
  cf$model <- nm
  cf$n     <- nobs(m)
  cf[, c("model", "n", "term", "estimate", "se", "df", "t", "p")]
}))
model_table

write.csv(model_table, "chill_coma_model_table.csv", row.names = FALSE)


# Model diagnostics ----

# Simulation-based residual checks (DHARMa): KS test of the residual
# distribution, a dispersion test and an outlier test. Plots are saved.
check_residuals <- function(model, name) {
  tryCatch({
    sim <- simulateResiduals(model, n = 1000, plot = FALSE)
    stats <- data.frame(
      analysis     = "chill-coma",
      model        = name,
      ks_p         = testUniformity(sim, plot = FALSE)$p.value,
      dispersion_p = testDispersion(sim, plot = FALSE)$p.value,
      outlier_p    = testOutliers(sim, plot = FALSE)$p.value
    )
    
    png(file.path("diagnostics", paste0("chillcoma_", name, ".png")),
        width = 1800, height = 900, res = 200)
    tryCatch(plot(sim), finally = invisible(dev.off()))
    stats
  }, error = function(e) {
    cat("DHARMa failed for", name, ":", conditionMessage(e), "\n")
    data.frame(analysis = "chill-coma", model = name,
               ks_p = NA, dispersion_p = NA, outlier_p = NA)
  })
}

chillcoma_diagnostics <- bind_rows(lapply(names(chill_models), function(nm) {
  check_residuals(chill_models[[nm]], nm)
}))
chillcoma_diagnostics

write.csv(chillcoma_diagnostics, "diagnostics/chillcoma_dharma_summary.csv",
          row.names = FALSE)


# Figure: body size ----

# Violin plot with embedded boxplot (median, IQR, whiskers to 1.5 x IQR),
# individual observations and sample sizes
plot_species_violin <- function(data, y, y_label) {
  counts <- data %>%
    group_by(Bee_species) %>%
    summarise(
      n = sum(!is.na({{ y }})),
      top = max({{ y }}, na.rm = TRUE),
      .groups = "drop"
    )
  
  ggplot(data, aes(Bee_species, {{ y }}, fill = Bee_species)) +
    geom_violin(alpha = 0.55, colour = NA) +
    geom_boxplot(
      width = 0.14,
      fill = "white",
      outlier.shape = NA,
      linewidth = 0.6
    ) +
    geom_jitter(
      aes(colour = Bee_species),
      width = 0.05,
      height = 0,
      size = 1.8,
      alpha = 0.7
    ) +
    geom_text(
      data = counts,
      aes(Bee_species, top, label = paste0("n = ", n)),
      inherit.aes = FALSE,
      vjust = -1.2,
      size = 6,
      fontface = "bold"
    ) +
    scale_fill_manual(values = species_colors) +
    scale_colour_manual(values = species_colors) +
    scale_x_discrete(labels = species_labels) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.14))) +
    labs(x = NULL, y = y_label) +
    theme_pub +
    theme(
      axis.title.y = element_text(size = 19, face = "bold"),
      axis.text.x = element_text(size = 17, face = "bold"),
      axis.text.y = element_text(size = 17, face = "bold"),
      legend.title = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 18, face = "bold")
    )
}

# ITD is recorded in cm (see Methods)
p_itd <- plot_species_violin(full_data, ITD_Final, "Intertegular distance (cm)")
p_itd

ggsave("new_ITD_violin.png", p_itd, dpi = 600, width = 6, height = 5)


# Figure: thermal performance ----

p_righting_temp <- plot_species_violin(full_data, Righting_Temperature,
                                       "Righting temperature (°C)")
p_righting_time <- plot_species_violin(full_data, log_righting,
                                       "Time to right (ln[s + 1])")
p_takeoff_temp  <- plot_species_violin(full_data, Take_off_Temperature,
                                       "Take-off temperature (°C)")
p_takeoff_time  <- plot_species_violin(full_data, log_takeoff,
                                       "Time to take-off (ln[s + 1])")

performance_figure <- (p_righting_temp | p_righting_time) /
  (p_takeoff_temp | p_takeoff_time) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold", size = 18))

performance_figure

ggsave("Multipanel_thermal_performance.png", performance_figure,
       dpi = 600, width = 12, height = 9, units = "in")


# Figure: thermal traits against body size ----

species_labels_bold <- c(
  "H.S" = expression(bolditalic(H.~scabiosae)),
  "L.M" = expression(bolditalic(L.~malachurum))
)

plot_trait_vs_itd <- function(data, y, y_label) {
  ggplot(data, aes(x = ITD_Final, y = {{ y }}, colour = Bee_species)) +
    geom_point(size = 3.5, alpha = 0.85) +
    geom_smooth(aes(fill = Bee_species), method = "lm", formula = y ~ x,
                se = TRUE, alpha = 0.18, linewidth = 1.2) +
    scale_colour_manual(values = species_colors, breaks = c("H.S", "L.M"),
                        labels = species_labels_bold, name = "Bee species") +
    scale_fill_manual(values = species_colors, breaks = c("H.S", "L.M"),
                      labels = species_labels_bold, guide = "none") +
    labs(x = "Intertegular distance (cm)", y = y_label) +
    theme_classic(base_size = 18) +
    theme(
      axis.title.x  = element_text(face = "bold", size = 20),
      axis.title.y  = element_text(face = "bold", size = 20),
      axis.text     = element_text(face = "bold", size = 16),
      legend.title  = element_text(face = "bold", size = 17),
      legend.text   = element_text(face = "bold", size = 16),
      legend.position = "top"
    )
}

traits_figure <-
  (plot_trait_vs_itd(full_data, Righting_Temperature, "Righting temperature (°C)") |
     plot_trait_vs_itd(full_data, log_righting, "Time to right (ln[s + 1])")) /
  (plot_trait_vs_itd(full_data, Take_off_Temperature, "Take-off temperature (°C)") |
     plot_trait_vs_itd(full_data, log_takeoff, "Time to take-off (ln[s + 1])")) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold", size = 20))

traits_figure

ggsave("Thermal_traits_vs_ITD_multipanel.png", traits_figure,
       width = 15, height = 11, units = "in", dpi = 600)
