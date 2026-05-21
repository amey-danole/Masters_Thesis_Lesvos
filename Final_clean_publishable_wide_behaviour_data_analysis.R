# ============================================================
# WIDE BEHAVIOURAL DATA ANALYSIS
# ============================================================

# -----------------------------
# 1. LOAD PACKAGES----
# -----------------------------
library(dplyr)
library(readr)
library(tidyr)
library(lubridate)
library(lme4)
library(lmerTest)
library(MuMIn)
library(performance)
library(broom.mixed)
library(mgcv)
library(ggplot2)

# ============================================================
# 2. LOAD DATA----
# ============================================================

Final_combined_temperature_data <- read.csv("Final_combined_temperature_data.csv")
Focal_behavioural_data          <- read.csv("Focal_behavioural_data.csv")
Wide_behavioural_data           <- read.csv("Wide_behavioural_data.csv")
Outdoor_temperature_logger      <- read.csv("Outdoor_temperature_logger.csv")
Outdoor_light_logger            <- read.csv("Outdoor_light_logger.csv")

# ============================================================
# 3. CLEAN BEHAVIOURAL DATA----
# ============================================================

Wide_behavioural_data <- Wide_behavioural_data %>%
  rename(Number_of_bees_6min = Number.of.bees_6min)

# ============================================================
# 4. ENVIRONMENTAL DATA STRUCTURE----
# ============================================================

env_wide <- Final_combined_temperature_data %>%
  select(Date, Time, Sensors, Avg_sensor_temperature) %>%
  group_by(Date, Time, Sensors) %>%
  summarise(Avg_sensor_temperature = mean(Avg_sensor_temperature, na.rm = TRUE),
            .groups = "drop") %>%
  pivot_wider(names_from = Sensors,
              values_from = Avg_sensor_temperature)

env_wide <- env_wide %>%
  left_join(Outdoor_temperature_logger, by = c("Date", "Time")) %>%
  left_join(Outdoor_light_logger, by = c("Date", "Time"))

# ============================================================
# 5. MERGE DATASETS----
# ============================================================

final_wide_behavioural_dataset <- Wide_behavioural_data %>%
  left_join(env_wide, by = c("Date", "Time"))

# FIX TYPES
final_wide_behavioural_dataset <- final_wide_behavioural_dataset %>%
  mutate(
    Date = as.factor(Date),
    Site = as.factor(Site),
    Camera = as.factor(Camera),
    Time = as.numeric(Time)
  )

# Exploring dataset
str(final_wide_behavioural_dataset)
names(final_wide_behavioural_dataset)

peak_activity <- final_wide_behavioural_dataset %>%
  group_by(Time) %>%
  summarise(
    mean_activity = mean(Avg_number_of_bees, na.rm = TRUE),
    sd_activity = sd(Avg_number_of_bees, na.rm = TRUE),
    n = n()
  ) %>%
  arrange(desc(mean_activity))

peak_activity

daily_activity <- final_wide_behavioural_dataset %>%
  group_by(Date) %>%
  summarise(
    mean_activity = mean(Avg_number_of_bees, na.rm = TRUE),
    sd_activity = sd(Avg_number_of_bees, na.rm = TRUE),
    n = n()
  ) %>%
  arrange(Date)

daily_activity

# ============================================================
# FULL MIXED MODEL
# ============================================================

model_data <- final_wide_behavioural_dataset %>%
  dplyr::select(
    Avg_number_of_bees,
    Thermochron_1, Thermochron_3, Thermochron_6,
    Thermochron_7, Thermochron_8, Thermochron_9, Thermochron_10,
    Avg_T_ambient,
    Avg_lux_light_logger,
    Site, Date
  ) %>%
  na.omit()

full_model <- lmer(
  Avg_number_of_bees ~
    Thermochron_1 +
    Thermochron_3 +
    Thermochron_6 +
    Thermochron_7 +
    Thermochron_8 +
    Thermochron_9 +
    Thermochron_10 +
    Avg_T_ambient +
    Avg_lux_light_logger +
    (1 | Site) +
    (1 | Date),
  data = model_data,
  REML = FALSE
)

# OPTIONAL MODEL SELECTION
options(na.action = "na.fail")
model_set <- dredge(full_model)

best_model <- get.models(model_set, 1)[[1]]

# ============================================================
# FINAL MODEL
# ============================================================

final_best_wide_model <- lmer(
  Avg_number_of_bees ~
    Thermochron_1 +
    Thermochron_10 +
    Thermochron_3 +
    Thermochron_6 +
    Thermochron_7 +
    Thermochron_8 +
    Thermochron_9 +
    (1 | Date),
  data = final_wide_behavioural_dataset,
  REML = TRUE
)
summary(final_best_wide_model)
# ============================================================
# COEFFICIENT EXTRACTION
# ============================================================

coef_data <- broom.mixed::tidy(final_best_wide_model,
                               effects = "fixed",
                               conf.int = TRUE) %>%
  filter(term != "(Intercept)")

# ============================================================
# LABELS
# ============================================================

coef_data$term <- dplyr::recode(
  coef_data$term,
  "Thermochron_1"  = "Ambient air",
  "Thermochron_3"  = "2–3 cm soil (Site 1)",
  "Thermochron_6"  = "2–3 cm soil (Site 2)",
  "Thermochron_7"  = "5 cm soil (Site 1)",
  "Thermochron_8"  = "5 cm soil (Site 2)",
  "Thermochron_9"  = "10 cm soil (Site 1)",
  "Thermochron_10" = "10 cm soil (Site 2)"
)

# ============================================================
# ORDERING (BIOLOGICAL STRUCTURE)
# ============================================================

coef_data$term <- factor(
  coef_data$term,
  levels = rev(c(
    "Ambient air",
    "2–3 cm soil (Site 1)",
    "2–3 cm soil (Site 2)",
    "5 cm soil (Site 1)",
    "5 cm soil (Site 2)",
    "10 cm soil (Site 1)",
    "10 cm soil (Site 2)"
  ))
)

# ============================================================
# GROUPING FOR COLOURS + PANELS
# ============================================================

coef_data$group <- case_when(
  coef_data$term == "Ambient air" ~ "Ambient air",
  coef_data$term %in% c("2–3 cm soil (Site 1)", "2–3 cm soil (Site 2)") ~ "Shallow soil (2–3 cm)",
  coef_data$term %in% c("5 cm soil (Site 1)", "5 cm soil (Site 2)") ~ "Middle soil (5 cm)",
  coef_data$term %in% c("10 cm soil (Site 1)", "10 cm soil (Site 2)") ~ "Deep soil (10 cm)"
)

coef_data$group <- factor(
  coef_data$group,
  levels = c("Ambient air",
             "Shallow soil (2–3 cm)",
             "Middle soil (5 cm)",
             "Deep soil (10 cm)")
)

# ============================================================
# COLOURS
# ============================================================

cols <- c(
  "Ambient air" = "#0072B2",
  "Shallow soil (2–3 cm)" = "#E69F00",
  "Middle soil (5 cm)" = "#009E73",
  "Deep soil (10 cm)" = "#CC79A7"
)

# ============================================================
# PANEL BLOCK POSITIONS
# ============================================================

y_levels <- levels(coef_data$term)

shallow_y <- which(y_levels %in%
                     c("2–3 cm soil (Site 1)", "2–3 cm soil (Site 2)"))

middle_y <- which(y_levels %in%
                    c("5 cm soil (Site 1)", "5 cm soil (Site 2)"))

deep_y <- which(y_levels %in%
                  c("10 cm soil (Site 1)", "10 cm soil (Site 2)"))

# ============================================================
# PLOT
# ============================================================

effect_plot <- ggplot(coef_data,
                      aes(x = estimate,
                          y = term,
                          colour = group)) +
  
  # PANEL BLOCKS
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = min(shallow_y) - 0.5,
           ymax = max(shallow_y) + 0.5,
           fill = cols["Shallow soil (2–3 cm)"],
           alpha = 0.20) +
  
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = min(middle_y) - 0.5,
           ymax = max(middle_y) + 0.5,
           fill = cols["Middle soil (5 cm)"],
           alpha = 0.20) +
  
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = min(deep_y) - 0.5,
           ymax = max(deep_y) + 0.5,
           fill = cols["Deep soil (10 cm)"],
           alpha = 0.20) +
  
  # EFFECTS
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.2,
                 na.rm = TRUE) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # COLOURS
  scale_colour_manual(values = cols) +
  
  # LABELS
  labs(
    x = "Effect size estimate",
    y = "Thermochron position",
    colour = "Microclimate category"
  ) +
  
  theme_classic(base_size = 15) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text  = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

effect_plot

# ============================================================
# SAVE FIGURE
# ============================================================

ggsave(
  filename = "wide_behavioural_effect_size_figure.png",
  plot = effect_plot,
  dpi = 600,
  width = 300,
  height = 240,
  units = "mm"
)

# ============================================================
# 6. PREPARE DATA FOR GAM ANALYSIS----
# ============================================================

gam_data <- final_wide_behavioural_dataset %>%
  dplyr::select(
    Avg_number_of_bees,
    Avg_T_ambient,
    Thermochron_3,
    Thermochron_6,
    Thermochron_7,
    Thermochron_8,
    Thermochron_9,
    Thermochron_10
  ) %>%
  mutate(
    shallow = rowMeans(cbind(Thermochron_3, Thermochron_6), na.rm = TRUE),
    middle  = rowMeans(cbind(Thermochron_7, Thermochron_8), na.rm = TRUE),
    deep    = rowMeans(cbind(Thermochron_9, Thermochron_10), na.rm = TRUE)
  ) %>%
  drop_na()

# ============================================================
# FIT GAM MODELS
# ============================================================

m_ambient <- gam(Avg_number_of_bees ~ s(Avg_T_ambient, k = 5),
                 data = gam_data, method = "REML")

m_shallow <- gam(Avg_number_of_bees ~ s(shallow, k = 5),
                 data = gam_data, method = "REML")

m_middle <- gam(Avg_number_of_bees ~ s(middle, k = 5),
                data = gam_data, method = "REML")

m_deep <- gam(Avg_number_of_bees ~ s(deep, k = 5),
              data = gam_data, method = "REML")

# ============================================================
# PREDICTION FUNCTION
# ============================================================

pred_fun <- function(model, var_name) {
  
  x_seq <- seq(
    min(gam_data[[var_name]], na.rm = TRUE),
    max(gam_data[[var_name]], na.rm = TRUE),
    length.out = 200
  )
  
  newdata <- data.frame(x_seq)
  names(newdata) <- var_name
  
  preds <- predict(model, newdata = newdata, se.fit = TRUE)
  
  data.frame(
    x = x_seq,
    fit = preds$fit,
    lower = preds$fit - 1.96 * preds$se.fit,
    upper = preds$fit + 1.96 * preds$se.fit
  )
}

# ============================================================
# BUILD PREDICTIONS 
# ============================================================

plot_data <- bind_rows(
  pred_fun(m_ambient, "Avg_T_ambient") %>%
    mutate(group = "Ambient air"),
  
  pred_fun(m_shallow, "shallow") %>%
    mutate(group = "Shallow soil (2–3 cm)"),
  
  pred_fun(m_middle, "middle") %>%
    mutate(group = "Middle soil (5 cm)"),
  
  pred_fun(m_deep, "deep") %>%
    mutate(group = "Deep soil (10 cm)")
)

# enforce order (legend + colours stability)
plot_data$group <- factor(
  plot_data$group,
  levels = c(
    "Ambient air",
    "Shallow soil (2–3 cm)",
    "Middle soil (5 cm)",
    "Deep soil (10 cm)"
  )
)

# ============================================================
# COLOUR PALETTE
# ============================================================

cols <- c(
  "Ambient air" = "#0072B2",
  "Shallow soil (2–3 cm)" = "#E69F00",
  "Middle soil (5 cm)" = "#009E73",
  "Deep soil (10 cm)" = "#CC79A7"
)

# ============================================================
# GAM PLOT 
# ============================================================

gam_plot <- ggplot(
  plot_data,
  aes(x = x, y = fit, colour = group, fill = group, group = group)
) +
  
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.20, colour = NA) +
  
  geom_line(linewidth = 1.2) +
  
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols) +
  
  labs(
    x = "Temperature (°C)",
    y = "Relative bee activity",
    colour = "Thermal environment",
    fill = "Thermal environment"
  ) +
  
  theme_classic(base_size = 15) +
  
  theme(
    axis.title.x = element_text(face = "bold", size = 18),
    axis.title.y = element_text(face = "bold", size = 18),
    
    axis.text.x  = element_text(face = "bold", size = 15),
    axis.text.y  = element_text(face = "bold", size = 15),
    
    legend.title = element_text(face = "bold", size = 16),
    legend.text  = element_text(face = "bold", size = 14)
  )

gam_plot

# ============================================================
# SAVE FIGURE
# ============================================================

ggsave(
  filename = "wide_behavioural_gam_figure.png",
  plot = gam_plot,
  dpi = 600,
  width = 400,
  height = 240,
  units = "mm"
)
