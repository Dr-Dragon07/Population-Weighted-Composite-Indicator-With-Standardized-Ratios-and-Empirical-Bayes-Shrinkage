# =============================================================================
# 1. LOAD LIBRARIES
# =============================================================================
# Ensure all necessary packages are installed
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  readr,      # For reading CSV files
  tidyverse,  # For data manipulation and plotting (includes dplyr, ggplot2)
  knitr,      # For creating well-formatted tables
  kableExtra, # For styling kable tables
  scales      # For formatting plot labels
)

# =============================================================================
# 2. DATA LOADING AND PREPARATION
# =============================================================================
# Read the provided CSV data.
# Note: In a real scenario, you would replace 'indicator_summary_advanced_stats.csv'
# with the actual path to your file.
indicator_data <- read_csv("C://Users//waqas//Downloads//ThesisFiles//SR//For each indicator separately//indicator_summary_advanced_stats.csv")

# Reshape the data from wide to long format for easier analysis.
# This creates separate rows for 'Segregated' and 'Complementary' areas for each indicator.
long_data <- indicator_data %>%
  pivot_longer(
    cols = -Indicator,
    names_to = c(".value", "Area"),
    names_sep = "_(?=[^_]+$)" # Splits the column name at the last underscore
  ) %>%
  rename(
    Mean_SR = mean_sr,
    Lower_CI = mean_lower_ci,
    Upper_CI = mean_upper_ci,
    Reliability = reliability
  ) %>%
  mutate(Area = case_when(
    Area == "seg" ~ "Segregated",
    Area == "comp" ~ "Complementary",
    TRUE ~ Area
  )) %>%
  # Calculate Standard Error (SE) from the 95% Confidence Intervals.
  # SE = (Upper CI - Lower CI) / (2 * 1.96)
  mutate(SE = (Upper_CI - Lower_CI) / (2 * 1.96))

# =============================================================================
# 3. VISUALIZATION FUNCTION (UNCHANGED)
# =============================================================================
#' Create a forest plot to visualize and compare the SR for the two areas.
create_forest_plot <- function(data, indicator_name) {
  ggplot(data, aes(y = Area, x = Mean_SR, xmin = Lower_CI, xmax = Upper_CI, color = Area)) +
    geom_point(shape = 15, size = 4) +
    geom_errorbarh(height = 0.15, linewidth = 1) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = "Comparison of Standardized Ratios (SR)",
      subtitle = str_wrap(indicator_name, width = 80),
      x = "Mean Standardized Ratio (SR) with 95% Confidence Interval",
      y = "Area Type"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 11, margin = margin(b = 15)),
      panel.grid.major.y = element_blank(),
      legend.position = "none"
    )
}

# =============================================================================
# 4. EXECUTION LOGIC WITH COMBINED TABLE
# =============================================================================
# Get a unique list of indicators from the data
indicators <- unique(long_data$Indicator)

# Loop through each indicator name
for (indicator in indicators) {
  
  # --- Step 1: Filter data for the current indicator ---
  indicator_subset <- long_data %>% filter(Indicator == indicator)
  
  # --- Step 2: Perform Z-test calculations ---
  seg_data <- indicator_subset %>% filter(Area == "Segregated")
  comp_data <- indicator_subset %>% filter(Area == "Complementary")
  
  # Initialize test result variables to handle cases with missing data
  z_stat <- NA
  p_value <- NA
  
  # Perform calculation only if data for both areas is present
  if (nrow(seg_data) > 0 && nrow(comp_data) > 0) {
    z_stat <- (seg_data$Mean_SR - comp_data$Mean_SR) / 
      sqrt(seg_data$SE^2 + comp_data$SE^2)
    p_value <- 2 * pnorm(-abs(z_stat))
  }
  
  interpretation <- case_when(
    is.na(p_value) ~ "Incomplete data",
    p_value < 0.05 ~ "Significant difference",
    TRUE ~ "No significant difference"
  )
  
  # --- Step 3: Create the Combined Summary Table ---
  combined_table <- indicator_subset %>%
    # Select and format descriptive statistics
    select(Area, Mean_SR, Lower_CI, Upper_CI, Reliability) %>%
    mutate(`95% CI` = paste0("[", round(Lower_CI, 3), ", ", round(Upper_CI, 3), "]")) %>%
    select(Area, Mean_SR, `95% CI`, Reliability) %>%
    # Add the Z-test results as new columns
    mutate(
      `Z-Statistic` = round(z_stat, 3),
      `P-Value` = round(p_value, 4),
      Interpretation = interpretation
    )
  
  # --- Step 4: Create Visualization ---
  forest_plot <- create_forest_plot(indicator_subset, indicator)
  
  # --- Step 5: Print Combined Results ---
  cat("\n", "===================================================================", "\n")
  cat("📊 ANALYSIS FOR INDICATOR:", indicator, "\n")
  cat("===================================================================", "\n\n")
  
  # Print the combined table
  print(
    kable(combined_table, caption = "Combined Summary and Statistical Test Results.") %>%
      kable_styling(bootstrap_options = c("striped", "hover", "responsive"), full_width = FALSE)
  )
  cat("\n") # Add space
  
  # Print the plot
  print(forest_plot)
  
  cat("\n\n") # Add separator
}

# =============================================================================
# 1. LOAD LIBRARIES
# =============================================================================
# Ensure all necessary packages are installed
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  readr,      # For reading CSV files
  dplyr,      # For data manipulation
  tidyr       # For data cleaning
)

# =============================================================================
# 2. DATA LOADING AND ANALYSIS
# ============================================================================

# Perform calculations to create the summary table
summary_table <- indicator_data %>%
  # Calculate Standard Error (SE) for both groups from their 95% CIs
  # The formula is SE = (Upper CI - Lower CI) / (2 * 1.96)
  mutate(
    se_seg = (mean_upper_ci_seg - mean_lower_ci_seg) / (2 * 1.96),
    se_comp = (mean_upper_ci_comp - mean_lower_ci_comp) / (2 * 1.96)
  ) %>%
  # Calculate the Z-statistic to compare the two means
  mutate(
    z_statistic = (mean_sr_seg - mean_sr_comp) / sqrt(se_seg^2 + se_comp^2)
  ) %>%
  # Calculate the two-tailed p-value from the Z-statistic
  mutate(
    p_value = 2 * pnorm(-abs(z_statistic))
  ) %>%
  # Add a clear interpretation of the p-value for the Z-test
  mutate(
    `Z-Test Significance` = case_when(
      is.na(p_value) ~ "Cannot determine",
      p_value < 0.05 ~ "Significant",
      TRUE ~ "Not Significant"
    )
  ) %>%
  # Add placeholder columns for other tests with explanations
  mutate(
    `T-Test` = "N/A (Requires sample sizes)",
    `Wilcoxon-Test` = "N/A (Requires raw data)",
    `Levene-Test` = "N/A (Requires raw data)"
  ) %>%
  # Create formatted 95% CI strings for cleaner output
  mutate(
    `Segregated 95% CI` = paste0("[", round(mean_lower_ci_seg, 3), ", ", round(mean_upper_ci_seg, 3), "]"),
    `Complementary 95% CI` = paste0("[", round(mean_lower_ci_comp, 3), ", ", round(mean_upper_ci_comp, 3), "]")
  ) %>%
  # Select and rename columns for the final table
  select(
    `Indicator` = Indicator,
    `Segregated Mean SR` = mean_sr_seg,
    `Segregated 95% CI`,
    `Complementary Mean SR` = mean_sr_comp,
    `Complementary 95% CI`,
    `Z-Statistic` = z_statistic,
    `P-Value` = p_value,
    `Z-Test Significance`,
    `T-Test`,
    `Wilcoxon-Test`,
    `Levene-Test`
  )

# =============================================================================
# 3. SAVE THE RESULTS
# =============================================================================
# Save the final summary table to a CSV file for easy use in Excel
write.csv(summary_table, "indicator_comparison_all_tests.csv", row.names = FALSE)

# Print a message to confirm the file has been saved
print(summary_table)
View(summary_table)