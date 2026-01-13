# Enhanced R Code for DM2 Analysis - High Quality Tables and Plots
# Load required libraries
library(readxl)
library(tidyverse)
library(knitr)
library(kableExtra)
library(ggplot2)
library(viridis)
library(scales)
library(gridExtra)
library(RColorBrewer)
library(plotly)
library(DT)
library(formattable)
library(openxlsx) # For Excel workbook creation

# Load the data
# Assuming the data is already loaded as DM_Related_to_work_2
if (!exists("DM_Related_to_work_2")) {
  stop("Data object 'DM_Related_to_work_2' not found. Please load your data first.")
}

df <- DM_Related_to_work_2
names(df) <- make.names(names(df))

# Define categories
dm2_categories <- list(
  pharmaceutical = c(
    "FV30a - Prescription - Alimentary tract and metabolism (ATC A)",
    "FV31a - Drug substitution - Alimentary tract and metabolism (ATC A)",
    "FV32a - Drug substitution rate - Alimentary tract and metabolism (ATC A)"
  ),
  general_practitioner_non_nefmi = c(
    "FV13 - Proportion of people aged 40-54 who are switching from medication for diabetes",
    "FV14 - Proportion of people aged 55-69 who are switching from medication for diabetes",
    "FV15 - Proportion of diabetic patients who underwent serum creatinine level determination",
    "FV16a - Proportion of people under 65 years of age receiving influenza vaccination among patients with hypertension, diabetes, ischemic heart disease, or COPD"
  ),
  general_practitioner_nefmi = c(
    "FV06 - Proportion of diabetic and/or hypertensive patients who underwent blood lipid testing",
    "FV08 - Proportion of diabetics who underwent hemoglobin A1c testing",
    "FV09 - Proportion of diabetics who underwent ophthalmological examination"
  ),
  metabolic = c(
    "PK15 - Frequency of limb amputation due to diabetes"
  ),
  professional_care = c(
    "PK16 - Frequency of surgery for diabetic retinopathy"
  )
)

# Create reverse lookup for categorization
indicator_to_category <- list()
for(category in names(dm2_categories)) {
  for(indicator in dm2_categories[[category]]) {
    indicator_to_category[[indicator]] <- category
  }
}

# Categorization function
categorize_indicator <- function(indicator) {
  if (!is.null(indicator_to_category[[indicator]])) {
    return(indicator_to_category[[indicator]])
  }
  
  for (category in names(dm2_categories)) {
    if (any(grepl(indicator, dm2_categories[[category]], fixed = TRUE))) {
      return(category)
    }
  }
  return("Other")
}

# Add category information to the data
df <- df %>%
  rowwise() %>%
  mutate(Category = categorize_indicator(Indicator)) %>%
  ungroup()

# Filter to keep only DM2-related indicators
dm2_data <- df %>%
  filter(Category != "Other")

# Check if required columns exist
required_cols <- c("OBS.cases.in.segregated.area", "EXP.cases.in.segregated.area",  
                   "OBS.cases.in.complementary.area", "EXP.cases.in.complementary.area",
                   "Standardized.Ratio", "ID.GMP", "Indicator")

missing_cols <- required_cols[!required_cols %in% names(dm2_data)]
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
}

# NEW: Calculate Global Parameters for Empirical Bayes Shrinkage
# Global SR (overall mean) to shrink estimates towards
global_sr <- sum(dm2_data$OBS.cases.in.segregated.area, dm2_data$OBS.cases.in.complementary.area, na.rm = TRUE) /
  sum(dm2_data$EXP.cases.in.segregated.area, dm2_data$EXP.cases.in.complementary.area, na.rm = TRUE)

# Shrinkage parameter 'k', estimated as the mean of expected counts.
k_shrinkage <- mean(c(dm2_data$EXP.cases.in.segregated.area, dm2_data$EXP.cases.in.complementary.area), na.rm = TRUE)


# MODIFIED: Function to calculate Standardized Ratios with Bayesian CI and reliability measures
calculate_sr_statistics <- function(data, population_col_obs, population_col_exp, area_type, global_sr, k) {
  
  # Group by Category and calculate statistics
  sr_stats <- data %>%
    group_by(Category) %>%
    summarise(
      Observed = sum(.data[[population_col_obs]], na.rm = TRUE),
      Expected = sum(.data[[population_col_exp]], na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    filter(Expected > 0) %>%
    mutate(
      # Calculate Standardized Ratio
      SR = Observed / Expected,
      Mean = SR, # The aggregated SR is the mean rate for the category
      
      # Bayesian 95% Credible Interval using Gamma-Poisson conjugate
      # Using Jeffrey's prior (Gamma(0.5, 0))
      alpha_posterior = Observed + 0.5,
      beta_posterior = Expected,
      
      # Bayesian Credible Interval
      Bayesian_CI_Lower = qgamma(0.025, alpha_posterior, beta_posterior),
      Bayesian_CI_Upper = qgamma(0.975, alpha_posterior, beta_posterior),
      
      # Frequentist 95% Confidence Interval (Poisson approximation)
      Freq_CI_Lower = ifelse(Observed == 0, 0, qchisq(0.025, 2 * Observed) / (2 * Expected)),
      Freq_CI_Upper = qchisq(0.975, 2 * (Observed + 1)) / (2 * Expected),
      
      # Reliability measures
      Variance = Observed / (Expected^2), # Variance of SR
      Standard_Error = sqrt(Variance),
      
      # Relative precision (inverse of coefficient of variation)
      CV = ifelse(SR > 0, Standard_Error / SR, NA), # Coefficient of Variation
      Reliability = ifelse(is.finite(1/CV), 1/CV, 0),
      
      # NEW: Shrinkage-based Reliability and Empirical Bayes SR (EB-SR)
      Reliability_Shrinkage = Expected / (Expected + k),
      Bayesian_SR = (Reliability_Shrinkage * SR) + ((1 - Reliability_Shrinkage) * global_sr),
      
      # Probability that SR > 1 (Bayesian)
      Prob_SR_greater_than_1 = 1 - pgamma(1, alpha_posterior, beta_posterior),
      
      # Area type
      Area = area_type
    )
  
  return(sr_stats)
}

# Calculate SR statistics for both segregated and complementary populations
sr_segregated <- calculate_sr_statistics(dm2_data, "OBS.cases.in.segregated.area",  
                                         "EXP.cases.in.segregated.area", "Segregated",
                                         global_sr, k_shrinkage)

sr_complementary <- calculate_sr_statistics(dm2_data, "OBS.cases.in.complementary.area",  
                                            "EXP.cases.in.complementary.area", "Complementary",
                                            global_sr, k_shrinkage)

# NEW: Calculate Mean and SD of Indicator-level SRs for each Category
# First, calculate SRs for each individual indicator
dm2_data_with_srs <- dm2_data %>%
  mutate(
    SR_Segregated = OBS.cases.in.segregated.area / EXP.cases.in.segregated.area,
    SR_Complementary = OBS.cases.in.complementary.area / EXP.cases.in.complementary.area
  ) %>%
  # Replace Inf with NA for stable calculations
  mutate(
    SR_Segregated = ifelse(is.infinite(SR_Segregated), NA, SR_Segregated),
    SR_Complementary = ifelse(is.infinite(SR_Complementary), NA, SR_Complementary)
  )

# Calculate Mean and SD for segregated area indicators
indicator_descriptives_seg <- dm2_data_with_srs %>%
  group_by(Category) %>%
  summarise(
    Mean_Indicator_SR = mean(SR_Segregated, na.rm = TRUE),
    SD_Indicator_SR = sd(SR_Segregated, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  rename_with(~ paste0(., "_Segregated"), -Category)

# Calculate Mean and SD for complementary area indicators
indicator_descriptives_comp <- dm2_data_with_srs %>%
  group_by(Category) %>%
  summarise(
    Mean_Indicator_SR = mean(SR_Complementary, na.rm = TRUE),
    SD_Indicator_SR = sd(SR_Complementary, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  rename_with(~ paste0(., "_Complementary"), -Category)

# Join the descriptive stats together
indicator_descriptives <- full_join(indicator_descriptives_seg, indicator_descriptives_comp, by = "Category")


# Combine results for comparison
combined_sr_stats <- bind_rows(sr_segregated, sr_complementary)

# MODIFIED: Create comparison table with both areas side by side, including all new stats
final_comparison_table <- sr_segregated %>%
  select(Category, SR, Mean, Bayesian_SR, CV, Reliability_Shrinkage, Bayesian_CI_Lower, Bayesian_CI_Upper) %>%
  rename_with(~ paste0(., "_Segregated"), -Category) %>%
  full_join(
    sr_complementary %>%
      select(Category, SR, Mean, Bayesian_SR, CV, Reliability_Shrinkage, Bayesian_CI_Lower, Bayesian_CI_Upper) %>%
      rename_with(~ paste0(., "_Complementary"), -Category),
    by = "Category"
  ) %>%
  # NEW: Join the Mean and SD of the individual indicators
  full_join(indicator_descriptives, by = "Category") %>%
  mutate(
    # Calculate difference between areas
    SR_Difference = SR_Segregated - SR_Complementary,
    Relative_Risk = SR_Segregated / SR_Complementary
  ) %>%
  # NEW: Reorder columns for a more logical layout
  select(
    Category,
    SR_Segregated, SR_Complementary,
    Mean_Indicator_SR_Segregated, SD_Indicator_SR_Segregated,
    Mean_Indicator_SR_Complementary, SD_Indicator_SR_Complementary,
    Bayesian_SR_Segregated, Bayesian_SR_Complementary,
    CV_Segregated, CV_Complementary,
    Reliability_Shrinkage_Segregated, Reliability_Shrinkage_Complementary,
    everything()
  )


# Create a styled HTML table using formattable
if (nrow(final_comparison_table) > 0) {
  formatted_comparison_table <- formattable(final_comparison_table,  
                                            list(
                                              area(col = c("SR_Segregated", "SR_Complementary")) ~ color_tile("white", "lightblue"),
                                              area(col = c("Reliability_Shrinkage_Segregated", "Reliability_Shrinkage_Complementary")) ~ color_tile("white", "lightgreen"),
                                              area(col = c("SR_Difference", "Relative_Risk")) ~ color_tile("lightcoral", "lightblue")
                                            )
  )
}

# Comprehensive summary table for all GMPs and categories using existing SR
comprehensive_summary <- dm2_data %>%
  group_by(ID.GMP, Category) %>%
  summarise(
    Mean_Standardized_Ratio = mean(Standardized.Ratio, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  pivot_wider(names_from = Category, values_from = Mean_Standardized_Ratio,  
              values_fill = NA)

# Create a styled HTML table for the comprehensive summary
if (nrow(comprehensive_summary) > 0) {
  formatted_summary_table <- formattable(comprehensive_summary,  
                                         list(area(col = 2:ncol(comprehensive_summary)) ~ color_tile("white", "orange"))
  )
}

# Generate plots

# 1. Comparison plot between segregated and complementary areas
# ==============================================================================

# --- Data Preparation for Plotting ---
# This part remains the same; it correctly pre-sorts the data.
sr_order <- combined_sr_stats %>%
  filter(Area == "Segregated") %>%
  arrange(desc(SR)) %>%
  pull(Category)

combined_sr_stats_ordered <- combined_sr_stats %>%
  mutate(Category = factor(Category, levels = sr_order))


# --- Plot 1: SR Comparison (THE FIX IS HERE) ---
# We add scale_color_identity() to make the white labels visible.

if (nrow(combined_sr_stats_ordered) > 0) {
  p1 <- ggplot(combined_sr_stats_ordered,  
               aes(x = Category, y = SR, fill = Area)) +
    geom_col(position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(ymin = Bayesian_CI_Lower, ymax = Bayesian_CI_Upper),
                  position = position_dodge(width = 0.9),  
                  width = 0.25,  
                  color = "gray20") +
    geom_text(
      aes(label = sprintf("%.2f", SR),
          vjust = ifelse(SR > 0.5, 1.5, -1.0),
          color = ifelse(SR > 0.5, "white", "black")),
      position = position_dodge(width = 0.9),
      size = 3.5
    ) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    # THIS LINE IS THE FIX: It tells ggplot to use the "white" and "black" values directly as colors.
    scale_color_identity(guide = "none") +  
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = "Standardized Ratio Comparison by Area",
      subtitle = "Categories sorted by Segregated Area SR. Error bars are 95% Bayesian Credible Intervals.",
      x = NULL,  
      y = "Standardized Ratio (SR)",
      fill = "Area Type"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", size = 16),
      legend.position = "top"
    )
  
  ggsave("sr_comparison_segregated_vs_complementary_FINAL.png", plot = p1,
         width = 12, height = 8, dpi = 300)
}


# --- Plot 2: Indicator-Level Details ---
# No changes needed here.

if (nrow(dm2_data) > 0) {
  dm2_data_plot <- dm2_data %>%
    mutate(Indicator_Wrapped = str_wrap(Indicator, width = 50))
  
  p2 <- ggplot(dm2_data_plot,  
               aes(x = reorder(Indicator_Wrapped, Standardized.Ratio),  
                   y = Standardized.Ratio,  
                   fill = Category)) +
    geom_col() +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    facet_wrap(~ Category, scales = "free_y", ncol = 2) +
    coord_flip() +
    scale_fill_viridis_d(option = "plasma") +
    labs(
      title = "Detailed View: Standardized Ratio by Indicator",
      subtitle = "Indicators sorted by SR within each category. Red line at SR = 1.0.",
      x = "Indicator",
      y = "Standardized Ratio (SR)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "gray90", color = "black"),
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(size = 9)
    )
  
  ggsave("standardized_ratio_by_category_plot_FINAL.png", plot = p2,
         width = 18, height = 12, dpi = 300)
}


# --- Plot 3: Reliability Comparison ---
# Added smart labels and the scale_color_identity() fix.

if (nrow(combined_sr_stats_ordered) > 0) {
  p3 <- ggplot(combined_sr_stats_ordered,  
               aes(x = Category, y = Reliability, fill = Area)) +
    geom_col(position = position_dodge(width = 0.9)) +
    # Add smart labels to this plot as well
    geom_text(
      aes(label = sprintf("%.1f", Reliability),
          hjust = ifelse(Reliability > 4, 1.2, -0.2), # Position inside or outside the bar
          color = ifelse(Reliability > 4, "white", "black")),
      position = position_dodge(width = 0.9),
      size = 3.5
    ) +
    coord_flip() +
    # Add the fix here too
    scale_color_identity(guide = "none") +  
    scale_fill_brewer(palette = "Accent") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + # Give space for labels
    labs(
      title = "Reliability of Standardized Ratio Estimates",
      subtitle = "Categories are sorted by the Segregated Area's SR. Higher scores are more stable.",
      x = "Category",
      y = "Reliability Score",
      fill = "Area Type"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "top",
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", size = 16)
    )
  
  ggsave("sr_reliability_comparison_FINAL.png", plot = p3,
         width = 12, height = 8, dpi = 300)
}
# Create interactive plots
tryCatch({
  if (exists("p1")) {
    interactive_plot1 <- ggplotly(p1, tooltip = c("x", "y", "fill"))
    print(interactive_plot1)
  }
  if (exists("p2")) {
    interactive_plot2 <- ggplotly(p2, tooltip = c("x", "y", "fill"))
  }
  if (exists("p3")) {
    interactive_plot3 <- ggplotly(p3, tooltip = c("x", "y", "fill"))
  }
}, error = function(e) {
  cat("Error creating interactive plots:", e$message, "\n")
})

# Create and save tables to CSV
write.csv(comprehensive_summary, "comprehensive_summary.csv", row.names = FALSE)
write.csv(final_comparison_table, "sr_comparison_segregated_vs_complementary.csv", row.names = FALSE)
write.csv(combined_sr_stats, "detailed_sr_statistics.csv", row.names = FALSE)

# Create interactive DT tables
dt_comprehensive <- datatable(
  comprehensive_summary,
  caption = 'Comprehensive Summary of Standardized Ratios by GMP and Category',
  extensions = 'Buttons',
  options = list(
    dom = 'Bfrtip',
    buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
    pageLength = 10,
    scrollX = TRUE
  )
)

# MODIFIED: DT Table now includes the new columns
dt_comparison <- datatable(
  final_comparison_table,
  caption = 'Standardized Ratio Comparison: Segregated vs Complementary Areas (with Mean, SD, Bayesian SR, and Reliability)',
  extensions = 'Buttons',
  options = list(
    dom = 'Bfrtip',
    buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
    pageLength = 10,
    scrollX = TRUE
  )
) %>%
  formatRound(columns = which(sapply(final_comparison_table, is.numeric)), digits = 3)


dt_detailed <- datatable(
  combined_sr_stats,
  caption = 'Detailed SR Statistics with Bayesian and Frequentist Analysis',
  extensions = 'Buttons',
  options = list(
    dom = 'Bfrtip',
    buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
    pageLength = 10,
    scrollX = TRUE
  )
) %>%
  formatRound(columns = c("SR", "Bayesian_CI_Lower", "Bayesian_CI_Upper",  
                          "Freq_CI_Lower", "Freq_CI_Upper", "Reliability", "Bayesian_SR", "CV"), digits = 3) %>%
  formatPercentage(columns = "Prob_SR_greater_than_1", digits = 1)

# Print the tables
print(dt_comprehensive)
print(dt_comparison)
print(dt_detailed)

# Save results to Excel workbook
wb <- createWorkbook()

# Add worksheets
addWorksheet(wb, "Comprehensive_Summary")
writeData(wb, "Comprehensive_Summary", comprehensive_summary)

# MODIFIED: Excel sheet now includes the new columns
addWorksheet(wb, "SR_Comparison_All_Stats")
writeData(wb, "SR_Comparison_All_Stats", final_comparison_table)

addWorksheet(wb, "Detailed_Statistics")
writeData(wb, "Detailed_Statistics", combined_sr_stats)

# Apply styling
headerStyle <- createStyle(fgFill = "#4F81BD", fontColour = "white", textDecoration = "bold")
addStyle(wb, "Comprehensive_Summary", headerStyle, rows = 1, cols = 1:ncol(comprehensive_summary))
addStyle(wb, "SR_Comparison_All_Stats", headerStyle, rows = 1, cols = 1:ncol(final_comparison_table))
addStyle(wb, "Detailed_Statistics", headerStyle, rows = 1, cols = 1:ncol(combined_sr_stats))

saveWorkbook(wb, "dm2_sr_analysis_results.xlsx", overwrite = TRUE)

# Print summary of key findings
cat("\n================================================================================")
cat("\n=== STANDARDIZED RATIO (SR) ANALYSIS SUMMARY ===")
cat("\n================================================================================")
cat("\nFiles created:")
cat("\n- Main SR comparison plot: 'sr_comparison_segregated_vs_complementary.png'")
cat("\n- SR difference plot: 'sr_difference_plot.png'")
cat("\n- Reliability comparison plot: 'sr_reliability_comparison.png'")
cat("\n- Original category plot: 'standardized_ratio_by_category_plot.png'")
cat("\n- Three CSV files with detailed results")
cat("\n- Excel workbook: 'dm2_sr_analysis_results.xlsx'")
cat("\n- Interactive HTML tables")

cat("\n\nKey findings from SR comparison:")
for(i in 1:nrow(final_comparison_table)) {
  row <- final_comparison_table[i,]
  cat(sprintf("\n%s:", row$Category))
  cat(sprintf("\n  Segregated SR: %.3f (Bayesian SR: %.3f, CI: %.3f-%.3f)",  
              row$SR_Segregated,
              row$Bayesian_SR_Segregated,
              row$Bayesian_CI_Lower_Segregated,
              row$Bayesian_CI_Upper_Segregated))
  cat(sprintf("\n  Complementary SR: %.3f (Bayesian SR: %.3f, CI: %.3f-%.3f)",  
              row$SR_Complementary,
              row$Bayesian_SR_Complementary,
              row$Bayesian_CI_Lower_Complementary,
              row$Bayesian_CI_Upper_Complementary))
  cat(sprintf("\n  Difference: %.3f, Relative Risk: %.3f",  
              row$SR_Difference, row$Relative_Risk))
}

cat("\n\nAnalysis completed successfully!")
cat("\n================================================================================")

# =============================================================================
# The rest of the script remains unchanged.
# =============================================================================

# =============================================================================
# Unified SR Analysis and Visualization Dashboard
# =============================================================================

library(tidyverse)
library(ggplot2)
library(viridis)
library(patchwork)

# -----------------------------------------------------------------------------
# Example inputs (replace with your own)
# combined_sr_stats <- ...
# dm2_data <- ...
# -----------------------------------------------------------------------------

# 1. Order categories by Segregated SR
sr_order <- combined_sr_stats %>%
  filter(Area == "Segregated") %>%
  arrange(desc(SR)) %>%
  pull(Category)

combined_sr_stats_ordered <- combined_sr_stats %>%
  mutate(Category = factor(Category, levels = sr_order))

# -----------------------------------------------------------------------------
# 2. Plot 1: SR Comparison (with smart labels)
# -----------------------------------------------------------------------------
p1 <- ggplot(combined_sr_stats_ordered, 
             aes(x = Category, y = SR, fill = Area)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = Bayesian_CI_Lower, ymax = Bayesian_CI_Upper),
                position = position_dodge(width = 0.9), 
                width = 0.25, 
                color = "gray20") +
  geom_text(
    aes(label = sprintf("%.2f", SR),
        vjust = ifelse(SR > 0.5, 1.5, -1.0),
        color = ifelse(SR > 0.5, "white", "black")),
    position = position_dodge(width = 0.9),
    size = 3.5
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  scale_color_identity(guide = "none") + 
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Standardized Ratio Comparison by Area",
    subtitle = "Categories sorted by Segregated Area SR. Error bars are 95% Bayesian Credible Intervals.",
    x = NULL, 
    y = "Standardized Ratio (SR)",
    fill = "Area Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title.position = "plot",
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "top"
  )

# -----------------------------------------------------------------------------
# 3. Plot 2: Indicator-Level Details
# -----------------------------------------------------------------------------
dm2_data_plot <- dm2_data %>%
  mutate(Indicator_Wrapped = str_wrap(Indicator, width = 50))

p2 <- ggplot(dm2_data_plot, 
             aes(x = reorder(Indicator_Wrapped, Standardized.Ratio), 
                 y = Standardized.Ratio, 
                 fill = Category)) +
  geom_col() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  facet_wrap(~ Category, scales = "free_y", ncol = 2) +
  coord_flip() +
  scale_fill_viridis_d(option = "plasma") +
  labs(
    title = "Detailed View: Standardized Ratio by Indicator",
    subtitle = "Indicators sorted by SR within each category. Red line at SR = 1.0.",
    x = "Indicator",
    y = "Standardized Ratio (SR)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 9)
  )

# -----------------------------------------------------------------------------
# 4. Plot 3: Reliability Comparison
# -----------------------------------------------------------------------------
required_packages <- c("dplyr", "ggplot2", "stringr", "forcats", "ggalt")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# 2. DATA PREPARATION (Assuming 'dm2_data' is already loaded)
# ------------------------------------------------------------------------------
# Smoothing parameter
alpha <- 10

# Calculate reliability and the DIFFERENCE between the two areas
dm2_data_top5 <- dm2_data %>%
  mutate(
    EXP.cases.in.segregated.area = as.numeric(EXP.cases.in.segregated.area),
    EXP.cases.in.complementary.area = as.numeric(EXP.cases.in.complementary.area),
    Reliability_Segregated = EXP.cases.in.segregated.area / (EXP.cases.in.segregated.area + alpha),
    Reliability_Complementary = EXP.cases.in.complementary.area / (EXP.cases.in.complementary.area + alpha),
    # Calculate the absolute difference, which we'll use for ranking
    Difference = abs(Reliability_Segregated - Reliability_Complementary)
  ) %>%
  # Select the Top 5 indicators with the largest difference within each category
  group_by(Category) %>%
  slice_max(order_by = Difference, n = 5) %>%
  ungroup() %>%
  # Wrap indicator text for plotting
  mutate(Indicator_Wrapped = str_wrap(Indicator, width = 60))

# 3. CREATE THE FOCUSED "TOP 5" PLOT
# ------------------------------------------------------------------------------
p3 <- ggplot(dm2_data_top5, 
             aes(y = fct_reorder2(Indicator_Wrapped, Reliability_Segregated, Difference), # Sort by SR and Difference
                 x = Reliability_Complementary, 
                 xend = Reliability_Segregated)) +
  geom_dumbbell(
    size = 1.2, color = "gray80",
    size_x = 3.5, size_xend = 3.5,
    colour_x = "#1f78b4", # Blue for Complementary
    colour_xend = "#e31a1c"  # Red for Segregated
  ) +
  # Use facets, but now they will only have 5 indicators each
  facet_wrap(~ Category, scales = "free_y", ncol = 1) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray50") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(
    title = "Top 5 Indicators with the Largest Reliability Difference",
    subtitle = "Comparing Segregated vs. Complementary Areas within each category.",
    x = "Reliability Score",
    y = NULL,
    caption = "Blue Point = Complementary Area | Red Point = Segregated Area"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "#e6e6e6"),
    strip.text = element_text(face = "bold", size = 12),
    plot.title.position = "plot",
    plot.title = element_text(face = "bold", size = 18),
    panel.spacing = unit(1.5, "lines")
  )

# Display the new, focused plot
print(p3)

# -----------------------------------------------------------------------------
# 5. Dashboard Layout
# -----------------------------------------------------------------------------
dashboard <- (p1 + p3) / p2 +
  plot_annotation(
    title = "DM2 Healthcare Quality Indicators Dashboard",
    subtitle = "Comparative view of SR, Reliability, and Indicator-level details",
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  )

print(p1)
print(p2)
print(p3)
print(dashboard)