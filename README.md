# *Some like it hot*: Thermal Performance and Activity Patterns of Halictid Bees on Lesvos Island, Greece

This repository contains the R scripts and datasets used in my Islands and Sustainability MSc thesis, which is now being revised for publication.

The project investigates how microclimate—including air temperature, soil temperature at multiple depths, and light intensity—relates to activity in two ground-nesting halictid bee species:

- *Halictus scabiosae*
- *Lasioglossum malachurum*

It also examines thermal performance following chill-coma induction.

## Chill-Coma Recovery and Thermal Performance

**Script:** `final_final_ chill coma recovery assay analysis.R`  
**Dataset:** `Chill_coma_recovery.csv`

**Analyses include:**

- Righting and take-off temperatures and times
- Species differences in intertegular distance (ITD) and thermal performance
- Associations with ambient temperature and ITD
- Whether recovery time accounts for temperature effects
- Model diagnostics and influence checks

**Outputs include:** model summaries, diagnostic results, and figures of ITD and thermal performance.

## Focal Behavioural Analysis

**Script:** `final_final_focal behavioural analysis.R`  
**Datasets:**

- `focal_behavioural_data.csv`
- `combined_environmental_variables.csv`

**Analyses include:**

- Nest-guarding occurrence and duration in *H. scabiosae*
- Species differences in abundance and associations with environmental variables
- Nest-entry and nest-exit frequencies
- Daily activity patterns and the effects of time of day
- Model diagnostics and predicted response curves

**Outputs include:** the merged focal dataset, model summaries, diagnostic results, and figures.

## Wide-Angle Behavioural Analysis

**Script:** `final_final_wide behavioural analysis.R`  
**Datasets:**

- `Wide_behavioural_data.csv`
- `Final_combined_temperature_data.csv`
- `Outdoor_temperature_logger.csv`

**Analyses include:**

- Overall bee activity in relation to ambient air and soil temperatures at three depths
- Site-specific differences in shallow- and middle-soil temperature relationships
- Robustness checks for sampling-date effects, within-day temperature variation, and smooth-term basis dimensions
- Model diagnostics

**Outputs include:** model summaries, diagnostic results, and pooled and site-specific activity–temperature figures.
