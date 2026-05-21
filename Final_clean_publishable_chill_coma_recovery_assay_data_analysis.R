# =========================================================
# CHILL COMA RECOVERY & THERMAL TRAIT ANALYSIS
# Halictus scabiosae & Lasioglossum malachurum
# =========================================================

# -------------------------
# 1. LOAD REQUIRED PACKAGES
# -------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(car)
library(DHARMa)
library(ggdag)
library(lubridate)

# -------------------------
# 2. IMPORT DATASETS
# -------------------------

full_data<- read.csv("Chill_coma_recovery_assay.csv")

# Exploring dataset
str(full_data)

# Convert "MM:SS" to total seconds safely
time_to_seconds <- function(x) {
  parts <- strsplit(x, ":")
  sapply(parts, function(z) {
    as.numeric(z[1]) * 60 + as.numeric(z[2])
  })
}

full_data$Righting_Duration_sec <- time_to_seconds(full_data$Righting_Duration)
full_data$Take_off_Duration_sec  <- time_to_seconds(full_data$Take_off_Duration)

# 2. Species-level summary table
species_summary <- full_data %>%
  group_by(Bee_species) %>%
  summarise(
    
    # ITD (rename ITD_Final → ITD)
    ITD = mean(ITD_Final, na.rm = TRUE),
    
    # Righting metrics
    Righting_T_ambient   = mean(Righting_T_ambient, na.rm = TRUE),
    Righting_Temperature = mean(Righting_Temperature, na.rm = TRUE),
    Righting_Duration    = mean(Righting_Duration_sec, na.rm = TRUE),
    
    # Take-off metrics
    Take_off_T_ambient   = mean(Take_off_T_ambient, na.rm = TRUE),
    Take_off_Temperature = mean(Take_off_Temperature, na.rm = TRUE),
    Take_off_Duration    = mean(Take_off_Duration_sec, na.rm = TRUE),
    
    n = n()
  )

species_summary

str(full_data)

full_data %>%
  group_by(Bee_species) %>%
  summarise(
    Mean_ITD = mean(ITD_Final, na.rm = TRUE),
    Mean_Righting_T_ambient = mean(Righting_T_ambient, na.rm = TRUE),
    Mean_Righting_Temperature = mean(Righting_Temperature, na.rm = TRUE),
    Mean_Righting_Duration = mean(Righting_Duration_sec, na.rm = TRUE),
    
    Mean_Takeoff_T_ambient = mean(Take_off_T_ambient, na.rm = TRUE),
    Mean_Takeoff_Temperature = mean(Take_off_Temperature, na.rm = TRUE),
    Mean_Takeoff_Duration = mean(Take_off_Duration_sec, na.rm = TRUE)
  )

print(
  full_data %>%
    group_by(Bee_species) %>%
    summarise(
      Mean_ITD = mean(ITD_Final, na.rm = TRUE),
      Mean_Righting_T_ambient = mean(Righting_T_ambient, na.rm = TRUE),
      Mean_Righting_Temperature = mean(Righting_Temperature, na.rm = TRUE),
      Mean_Righting_Duration = mean(Righting_Duration_sec, na.rm = TRUE),
      
      Mean_Takeoff_T_ambient = mean(Take_off_T_ambient, na.rm = TRUE),
      Mean_Takeoff_Temperature = mean(Take_off_Temperature, na.rm = TRUE),
      Mean_Takeoff_Duration = mean(Take_off_Duration_sec, na.rm = TRUE)
    ),
  width = Inf
)

# -------------------------
# 2.1 BODY SIZE SUMMARY (ITD)
# -------------------------
itd_means <- full_data %>%
  group_by(Bee_species) %>%
  summarise(mean_ITD = mean(ITD_Final, na.rm = TRUE))

itd_means

# Extract species-specific means
H_scab <- itd_means$mean_ITD[itd_means$Bee_species == "H.S"]
L_mal  <- itd_means$mean_ITD[itd_means$Bee_species == "L.M"]

# Compute size ratios
ratio_H_L <- H_scab / L_mal
ratio_L_H <- L_mal / H_scab

ratio_H_L
ratio_L_H

# -------------------------
# 3. DATA PREPARATION & TRANSFORMATIONS
# -------------------------

# Convert behavioural durations (mm:ss → seconds)
full_data$Righting_seconds <- sapply(strsplit(full_data$Righting_Duration, ":"), function(x)
  as.numeric(x[1]) * 60 + as.numeric(x[2]))

full_data$Take_off_seconds <- sapply(strsplit(full_data$Take_off_Duration, ":"), function(x)
  as.numeric(x[1]) * 60 + as.numeric(x[2]))

# Convert grouping variables
full_data$Bee_species <- as.factor(full_data$Bee_species)
full_data$Date <- as.factor(full_data$Date)

# Log-transform behavioural durations
full_data$log_righting <- log(full_data$Righting_seconds + 1)
full_data$log_takeoff  <- log(full_data$Take_off_seconds + 1)

# -------------------------
# 4. EXPLORATORY ANALYSES
# -------------------------

# Check collinearity between body size and mass
cor(full_data$ITD_Final, full_data$Bee_weight, use = "complete.obs")

# Test species differences in body size
lm(ITD_Final ~ Bee_species, data = full_data)

# -------------------------
# 5. THERMAL PERFORMANCE MODELS
# -------------------------

# Righting temperature model
model_righting_temp <- lmer(
  Righting_Temperature ~ Bee_species + ITD_Final + Righting_T_ambient + (1|Date),
  data = full_data
)
summary(model_righting_temp)

model_righting_temp_species <- lmer(
  Righting_Temperature ~ Bee_species + Righting_T_ambient + (1|Date),
  data = full_data
)
summary(model_righting_temp_species)

# Take-off temperature model (full)
model_takeoff_temp <- lmer(
  Take_off_Temperature ~ Bee_species + ITD_Final + Take_off_T_ambient + (1|Date),
  data = full_data
)
summary(model_takeoff_temp)

# Take-off temperature model (reduced collinearity)
model_takeoff_temp_species <- lmer(
  Take_off_Temperature ~ Bee_species + Take_off_T_ambient + (1|Date),
  data = full_data
)
summary(model_takeoff_temp_species)

# -------------------------
# 6. BEHAVIOURAL DURATION MODELS
# -------------------------

# Righting duration model (full)
model_righting_duration <- lmer(
  log_righting ~ Bee_species + ITD_Final + Righting_T_ambient + (1|Date),
  data = full_data
)
summary(model_righting_duration)

# Collinearity check (VIF)
vif(model_righting_duration)

# Righting duration (final reduced models)
model_righting_duration_ITD <- lmer(
  log_righting ~ ITD_Final + Righting_T_ambient + (1|Date),
  data = full_data
)
summary(model_righting_duration_ITD)

model_righting_duration_species <- lmer(
  log_righting ~ Bee_species + Righting_T_ambient + (1|Date),
  data = full_data
)
summary(model_righting_duration_species)

# Take-off duration model
model_takeoff_duration_ITD <- lmer(
  log_takeoff ~ ITD_Final + Take_off_T_ambient + (1|Date),
  data = full_data
)
summary(model_takeoff_duration_ITD)

model_takeoff_duration_species <- lmer(
  log_takeoff ~ Bee_species + Take_off_T_ambient + (1|Date),
  data = full_data
)
summary(model_takeoff_duration_species)

# -------------------------
# 7. MODEL DIAGNOSTICS
# -------------------------

# Righting temperature diagnostics
res_righting_temp <- simulateResiduals(model_righting_temp)
plot(res_righting_temp)
testDispersion(res_righting_temp)
testOutliers(res_righting_temp)
testInfluential(res_righting_temp)

# Identify influential observations
infl_obs <- influence(model_righting_temp, obs = TRUE)
plot(cooks.distance(infl_obs))
which(cooks.distance(infl_obs) > (4 / nrow(full_data)))
full_data[c(20,24), ]

# Sensitivity analysis (excluding influential points)
model_no_influential <- lmer(
  Righting_Temperature ~ Bee_species + ITD_Final + Righting_T_ambient + (1|Date),
  data = full_data[-c(20,24), ]
)
summary(model_no_influential)

# Take-off temperature diagnostics (full)
res_takeoff_temp <- simulateResiduals(model_takeoff_temp)
plot(res_takeoff_temp)
testDispersion(res_takeoff_temp)
testOutliers(res_takeoff_temp)

# Take-off temperature diagnostics (simple)
res_takeoff_temp_simple <- simulateResiduals(model_takeoff_temp_simple)
plot(res_takeoff_temp_simple)
testDispersion(res_takeoff_temp_simple)
testOutliers(res_takeoff_temp_simple)

# Righting duration diagnostics
res_righting_duration <- simulateResiduals(model_righting_duration_ITD)
plot(res_righting_duration)
testDispersion(res_righting_duration)
testOutliers(res_righting_duration)

# Take-off duration diagnostics
res_takeoff_duration <- simulateResiduals(model_takeoff_duration)
plot(res_takeoff_duration)
testDispersion(res_takeoff_duration)
testOutliers(res_takeoff_duration)

# -------------------------
# 8. MODEL INTERACTION TEST
# -------------------------

model_interaction <- lmer(
  Take_off_Temperature ~ ITD_Final * Take_off_T_ambient + (1|Date),
  data = full_data
)

anova(model_takeoff_temp, model_interaction)
summary(model_interaction)

# -------------------------
# 9. CONCEPTUAL DAG
# -------------------------

dag <- dagify(
  ITD ~ Species,
  Take_off_temp ~ ITD + Ambient_temp,
  Righting_temp ~ Ambient_temp,
  Duration ~ ITD + Ambient_temp,
  exposure = c("Species", "ITD", "Ambient_temp")
)

ggdag(dag) + theme_dag()

# =========================================================
# FIGURE GENERATION — PUBLICATION QUALITY PLOTS
# =========================================================

library(ggplot2)
library(dplyr)

# -------------------------
# 1. SPECIES RE-CODING
# -------------------------
full_data$Bee_species <- factor(full_data$Bee_species,
                                levels = c("H.S", "L.M"))

# -------------------------
# 2. COLOUR PALETTE (COLOURBLIND SAFE)
# -------------------------
species_colors <- c(
  "H.S" = "#E69F00",   # Halictus scabiosae
  "L.M" = "#4E4E4E"    # Lasioglossum malachurum
)

# -------------------------
# 3. PUBLICATION THEME
# -------------------------
theme_pub <- theme_classic() +
  theme(
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 11),
    legend.position = "none"
  )

# -------------------------
# 4. ITALIC SPECIES LABELS
# -------------------------
species_labels <- c(
  "H.S" = expression(italic(Halictus~scabiosae)),
  "L.M" = expression(italic(Lasioglossum~malachurum))
)

# -------------------------
# FIGURE 1 — INTERTEGULAR DISTANCE
# -------------------------
p_itd <- ggplot(full_data,
                aes(x = Bee_species, y = ITD_Final, fill = Bee_species)) +
  geom_boxplot(alpha = 0.8, width = 0.6) +
  scale_fill_manual(values = species_colors) +
  scale_x_discrete(labels = species_labels) +
  theme_pub +
  labs(
    x = expression(bold("Species")),
    y = expression(bold("Intertegular distance (mm)"))
  )

# -------------------------
# FIGURE 2 — RIGHTING TEMPERATURE
# -------------------------
p_righting_temp <- ggplot(full_data,
                          aes(x = Bee_species,
                              y = Righting_Temperature,
                              fill = Bee_species)) +
  geom_boxplot(alpha = 0.8) +
  scale_fill_manual(values = species_colors) +
  scale_x_discrete(labels = species_labels) +
  theme_pub +
  labs(
    x = expression(bold("Species")),
    y = expression(bold("Righting temperature (°C)"))
  )

# -------------------------
# FIGURE 3 — TAKE-OFF TEMPERATURE
# -------------------------
p_takeoff_temp <- ggplot(full_data,
                         aes(x = Bee_species,
                             y = Take_off_Temperature,
                             fill = Bee_species)) +
  geom_boxplot(alpha = 0.8) +
  scale_fill_manual(values = species_colors) +
  scale_x_discrete(labels = species_labels) +
  theme_pub +
  labs(
    x = expression(bold("Species")),
    y = expression(bold("Take-off temperature (°C)"))
  )

# -------------------------
# FIGURE 4 — RIGHTING DURATION
# -------------------------
p_righting_dur <- ggplot(full_data,
                         aes(x = Bee_species,
                             y = log_righting,
                             fill = Bee_species)) +
  geom_boxplot(alpha = 0.8) +
  scale_fill_manual(values = species_colors) +
  scale_x_discrete(labels = species_labels) +
  theme_pub +
  labs(
    x = expression(bold("Species")),
    y = expression(bold("Log righting duration (s)"))
  )

# -------------------------
# FIGURE 5 — TAKE-OFF DURATION
# -------------------------
p_takeoff_dur <- ggplot(full_data,
                        aes(x = Bee_species,
                            y = log_takeoff,
                            fill = Bee_species)) +
  geom_boxplot(alpha = 0.8) +
  scale_fill_manual(values = species_colors) +
  scale_x_discrete(labels = species_labels) +
  theme_pub +
  labs(
    x = expression(bold("Species")),
    y = expression(bold("Log take-off duration (s)"))
  )

# -------------------------
# 5. EXPORT FIGURES
# -------------------------
ggsave("Fig1_ITD.png", p_itd, dpi = 600, width = 6, height = 4)
ggsave("Fig2_RightingTemp.png", p_righting_temp, dpi = 600, width = 6, height = 4)
ggsave("Fig3_TakeoffTemp.png", p_takeoff_temp, dpi = 600, width = 6, height = 4)
ggsave("Fig4_RightingDuration.png", p_righting_dur, dpi = 600, width = 6, height = 4)
ggsave("Fig5_TakeoffDuration.png", p_takeoff_dur, dpi = 600, width = 6, height = 4)