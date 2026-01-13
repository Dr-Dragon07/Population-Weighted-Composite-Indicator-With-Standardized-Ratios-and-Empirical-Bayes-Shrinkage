library(readxl)
library(SpatialEpi)
library(tidyverse)
library(reshape2)
library(dplyr)
library(sp)
library(ggplot2)
library(scales)
library(viridis)
library(patchwork)
library(knitr)


# Load data (replace with your actual data loading method)
# df <- read_excel("DM_Related_to_work_2.xlsx")
df <- DM_Related_to_work_2  # Remove this if loading from file

# Load the necessary library
library(dplyr)

# Clean column names to make them easier to work with in R
# R replaces spaces and special characters with dots (.)
names(df) <- make.names(names(df))

# -----------------------------------------------------------------------------
# 1. Calculate Standardized Ratios (SR)
# -----------------------------------------------------------------------------
# We handle cases where the expected count is 0 to avoid division-by-zero errors.
df <- df %>%
  mutate(
    SR_segregated = ifelse(EXP.cases.in.segregated.area == 0, NA, OBS.cases.in.segregated.area / EXP.cases.in.segregated.area),
    SR_complementary = ifelse(EXP.cases.in.complementary.area == 0, NA, OBS.cases.in.complementary.area / EXP.cases.in.complementary.area)
  )

# -----------------------------------------------------------------------------
# 2. Calculate 95% Confidence Intervals (Lower and Upper Limits)
# -----------------------------------------------------------------------------
# This uses a standard log-normal approximation for the CI of a ratio.
# We also handle cases where the observed count is 0.
df <- df %>%
  mutate(
    # Confidence Interval for Segregated Area SR
    lower_ci_seg = ifelse(OBS.cases.in.segregated.area > 0, SR_segregated * exp(-1.96 / sqrt(OBS.cases.in.segregated.area)), 0),
    upper_ci_seg = ifelse(OBS.cases.in.segregated.area > 0, SR_segregated * exp(1.96 / sqrt(OBS.cases.in.segregated.area)), NA),
    
    # Confidence Interval for Complementary Area SR
    lower_ci_comp = ifelse(OBS.cases.in.complementary.area > 0, SR_complementary * exp(-1.96 / sqrt(OBS.cases.in.complementary.area)), 0),
    upper_ci_comp = ifelse(OBS.cases.in.complementary.area > 0, SR_complementary * exp(1.96 / sqrt(OBS.cases.in.complementary.area)), NA)
  )

# -----------------------------------------------------------------------------
# 3. Calculate Indicator Reliability
# -----------------------------------------------------------------------------
# Reliability = (True Variance) / (Total Variance)
# Where True Variance = Total Variance - Noise Variance (random error)
indicator_reliability <- df %>%
  group_by(Indicator) %>%
  summarise(
    # Reliability for Segregated Area
    total_var_seg = var(SR_segregated, na.rm = TRUE),
    noise_var_seg = mean((SR_segregated^2 / OBS.cases.in.segregated.area), na.rm = TRUE),
    true_var_seg = total_var_seg - noise_var_seg,
    reliability_seg = ifelse(total_var_seg > 0, pmax(0, true_var_seg / total_var_seg), NA),
    
    # Reliability for Complementary Area
    total_var_comp = var(SR_complementary, na.rm = TRUE),
    noise_var_comp = mean((SR_complementary^2 / OBS.cases.in.complementary.area), na.rm = TRUE),
    true_var_comp = total_var_comp - noise_var_comp,
    reliability_comp = ifelse(total_var_comp > 0, pmax(0, true_var_comp / total_var_comp), NA),
    .groups = 'drop'
  )

# -----------------------------------------------------------------------------
# 4. Calculate Bayesian Shrinkage Estimator for SR
# -----------------------------------------------------------------------------
# This "shrinks" the provider-level SR towards the indicator's mean SR.
# The amount of shrinkage is based on the reliability calculated above.
df_bayes <- df %>%
  group_by(Indicator) %>%
  mutate(
    mean_sr_seg = mean(SR_segregated, na.rm = TRUE),
    mean_sr_comp = mean(SR_complementary, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  left_join(indicator_reliability, by = "Indicator") %>%
  mutate(
    bayes_sr_seg = (reliability_seg * SR_segregated) + ((1 - reliability_seg) * mean_sr_seg),
    bayes_sr_comp = (reliability_comp * SR_complementary) + ((1 - reliability_comp) * mean_sr_comp)
  )

# -----------------------------------------------------------------------------
# 5. Create a Final Summary Table for Each Indicator
# -----------------------------------------------------------------------------
summary_table <- df_bayes %>%
  group_by(Indicator) %>%
  summarise(
    # --- Segregated Area Summary ---
    mean_sr_seg = mean(SR_segregated, na.rm = TRUE),
    median_sr_seg = median(SR_segregated, na.rm = TRUE),
    mean_lower_ci_seg = mean(lower_ci_seg, na.rm = TRUE),
    mean_upper_ci_seg = mean(upper_ci_seg, na.rm = TRUE),
    reliability_seg = first(reliability_seg),
    mean_bayes_sr_seg = mean(bayes_sr_seg, na.rm = TRUE),
    
    # --- Complementary Area Summary ---
    mean_sr_comp = mean(SR_complementary, na.rm = TRUE),
    median_sr_comp = median(SR_complementary, na.rm = TRUE),
    mean_lower_ci_comp = mean(lower_ci_comp, na.rm = TRUE),
    mean_upper_ci_comp = mean(upper_ci_comp, na.rm = TRUE),
    reliability_comp = first(reliability_comp),
    mean_bayes_sr_comp = mean(bayes_sr_comp, na.rm = TRUE),
    .groups = 'drop'
  )

# Print the final summary table to the console
print(summary_table)

# Optional: Save the detailed provider-level data to a CSV file
write.csv(df_bayes, "provider_level_advanced_stats.csv", row.names = FALSE)

# Optional: Save the summary table to a CSV file
write.csv(summary_table, "indicator_summary_advanced_stats.csv", row.names = FALSE)