This repository contains R scripts and datasets used in my Islands and Sustainability MSc thesis:

**“The Spring Sting: Activity Patterns and Thermal Performance of Halictid Bees”**

The project investigates how microclimate (air temperature, soil temperature layers, and light intensity) influences activity patterns in:
- *Halictus scabiosae*
- *Lasioglossum malachurum*

and examines species thermal performance in a separate chill coma recovery assay.

---

## Chill Coma Recovery & Thermal Performance Analysis

**Script:** `Chill_coma_recovery_assay_data_analysis.R`  
**Dataset:** `Chill_coma_recovery_assay.csv`

**What it does:**
- Processes chill coma recovery and thermal assay data
- Converts behavioural durations to seconds
- Analyses righting and take-off temperatures and durations
- Tests species differences using linear mixed models (LMMs)
- Examines intertegular distance (body size) effects
- Produces figures

---

## Focal Behavioural Analysis

**Script:** `Final_clean_publishable_focal_behaviorual_data_analysis.R`  
**Datasets:**
- `focal_behavioural_data.csv`
- `combined_environmental_variables.csv`

**What it does:**
- Merges focal nest behaviour with environmental microclimate data
- Models nest guarding presence and duration (GLMM + LMM)
- Analyses incoming and outgoing flight frequencies
- Tests effects of soil temperature layers and light intensity
- Produces predicted response curves and composite figures

---

## Wide Behavioural Analysis

**Script:** `Wide_behavioural_data_analysis.R`  
**Datasets:**
- `Wide_behavioural_data.csv`
- `Final_combined_temperature_data.csv`
- `Outdoor_temperature_logger.csv`
- `Outdoor_light_logger.csv`

**What it does:**
- Merges wide-scale bee activity data with microclimate + light logger data
- Summarises temporal patterns in activity (daily + peak activity)
- Fits mixed models linking bee abundance to temperature across air and soil layers
- Uses model selection (AIC / dredge) to identify key predictors
- Extracts and visualises effect sizes across environmental gradients
- Fits GAMs to capture non-linear relationships between temperature and activity
- Produces:
  - Effect size plots across microclimatic layers
  - Temperature–activity response curves
  - Multi-panel GAM figures across depth gradients

---
