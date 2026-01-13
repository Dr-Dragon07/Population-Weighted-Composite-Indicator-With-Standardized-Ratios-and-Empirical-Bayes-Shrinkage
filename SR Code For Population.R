# =============================================================================
# DM2 Healthcare Quality Indicators Analysis
# =============================================================================

# Load required libraries
library(tidyverse)
library(ggplot2)
library(viridis)
library(patchwork)

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
# Replace with actual data import if needed:
# df <- read_excel("your_file.xlsx")
df <- DM_Related_to_work_2  
# Optional libraries for enhanced functionality
if (!require(kableExtra, quietly = TRUE)) {
  cat("Note: 'kableExtra' not found. Tables will use basic formatting.\n")
  cat("Install with: install.packages('kableExtra')\n")
}

if (!require(plotly, quietly = TRUE)) {
  cat("Note: 'plotly' not found. Interactive plots will not be available.\n")
}

if (!require(RColorBrewer, quietly = TRUE)) {
  cat("Note: 'RColorBrewer' not found. Using default colors.\n")
}

# Custom theme for consistent plotting
theme_custom <- function() {
  theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray60"),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 9),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      strip.text = element_text(size = 10, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

# Data Preparation Functions
# ==========================

prepare_dm2_categories <- function() {
  list(
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
}

create_indicator_lookup <- function(dm2_categories) {
  indicator_to_category <- list()
  for(category in names(dm2_categories)) {
    for(indicator in dm2_categories[[category]]) {
      indicator_to_category[[indicator]] <- category
    }
  }
  return(indicator_to_category)
}

categorize_indicator <- function(indicator, indicator_to_category, dm2_categories) {
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

# Statistical Analysis Functions
# ==============================

calculate_confidence_intervals <- function(observed, expected) {
  ci_results <- map2_dfr(observed, expected, function(obs, exp) {
    if (obs > 0 && exp > 0) {
      test_result <- poisson.test(obs, T = exp)
      return(data.frame(
        CI_lower = test_result$conf.int[1],
        CI_upper = test_result$conf.int[2]
      ))
    } else {
      return(data.frame(CI_lower = NA, CI_upper = NA))
    }
  })
  return(ci_results)
}

empirical_bayes_smooth <- function(observed, expected, alpha = 10) {
  overall_rate <- sum(observed, na.rm = TRUE) / sum(expected, na.rm = TRUE)
  eb_sr <- (observed + alpha * overall_rate) / (expected + alpha)
  return(eb_sr)
}

# Visualization Functions
# =======================

create_sr_forest_plot <- function(composite_data) {
  # Create unique identifier for each row to handle duplicates
  plot_data <- composite_data %>%
    arrange(Category, SR) %>%
    mutate(
      unique_id = paste(`ID GMP`, Category, row_number(), sep = "_"),
      ID_display = `ID GMP`
    )
  
  p <- plot_data %>%
    ggplot(aes(x = SR, y = reorder(unique_id, SR), color = Category)) +
    geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.7, size = 0.8) +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), 
                   height = 0.3, alpha = 0.7, size = 0.5) +
    geom_point(size = 2, alpha = 0.8) +
    facet_wrap(~ Category, scales = "free_y", ncol = 1) +
    scale_color_viridis_d(name = "Category") +
    scale_x_continuous(trans = "log10", 
                       labels = function(x) format(x, digits = 2, nsmall = 1)) +
    scale_y_discrete(labels = function(x) str_extract(x, "^[^_]+")) +  # Extract original ID for display
    labs(
      title = "Standardized Ratios (SR) with 95% Confidence Intervals",
      subtitle = "Forest plot showing SR by healthcare region and category",
      x = "SR (log scale)",
      y = "Healthcare Region (ID GMP)",
      caption = "Dashed line represents SR = 1 (expected rate)"
    ) +
    theme_custom() +
    theme(
      axis.text.y = element_text(size = 7),
      strip.text = element_text(size = 9)
    )
  
  return(p)
}

create_sr_comparison_plot <- function(composite_data) {
  p <- composite_data %>%
    ggplot(aes(x = SR, y = EB_SR, color = Category)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 0.7) +
    geom_point(size = 2.5, alpha = 0.7) +
    scale_color_viridis_d(name = "Category") +
    scale_x_continuous(trans = "log10", 
                       labels = function(x) format(x, digits = 2, nsmall = 1)) +
    scale_y_continuous(trans = "log10", 
                       labels = function(x) format(x, digits = 2, nsmall = 1)) +
    labs(
      title = "Raw SR vs Empirical Bayes Adjusted SR",
      subtitle = "Comparison showing shrinkage effect of EB smoothing",
      x = "Raw SR (log scale)",
      y = "EB-Adjusted SR (log scale)",
      caption = "Points below the diagonal line show shrinkage toward the overall mean"
    ) +
    theme_custom() +
    guides(color = guide_legend(override.aes = list(size = 3)))
  
  return(p)
}

create_category_distribution_plot <- function(composite_data) {
  p <- composite_data %>%
    ggplot(aes(x = reorder(Category, SR, median), y = SR, fill = Category)) +
    geom_violin(alpha = 0.7, trim = FALSE) +
    geom_boxplot(width = 0.3, alpha = 0.8, outlier.alpha = 0.6) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
    scale_fill_viridis_d(name = "Category") +
    scale_y_continuous(trans = "log10", 
                       labels = function(x) format(x, digits = 2, nsmall = 1)) +
    labs(
      title = "Distribution of SR Values by Category",
      subtitle = "Violin plots showing density distribution with box plots and individual points",
      x = "Category (ordered by median SR)",
      y = "SR (log scale)"
    ) +
    theme_custom() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(p)
}

create_reliability_plot <- function(composite_data) {
  composite_data <- composite_data %>%
    mutate(
      SE_SR = sqrt(SR / Total_EXP),
      Reliability = 1 / (1 + SE_SR)
    )
  
  p <- composite_data %>%
    ggplot(aes(x = Total_EXP, y = SE_SR, color = Category, size = SR)) +
    geom_point(alpha = 0.7) +
    scale_color_viridis_d(name = "Category") +
    scale_x_continuous(trans = "log10") +
    scale_size_continuous(name = "SR", range = c(1, 4)) +
    labs(
      title = "Standard Error vs Expected Cases",
      subtitle = "Larger expected cases lead to more reliable SR estimates",
      x = "Total Expected Cases (log scale)",
      y = "Standard Error of SR",
      caption = "Point size represents SR magnitude"
    ) +
    theme_custom()
  
  return(p)
}

# Table Creation Functions
# ========================

create_summary_table <- function(composite_data) {
  summary_stats <- composite_data %>%
    group_by(Category) %>%
    summarise(
      Regions = n(),
      Mean_SR = mean(SR, na.rm = TRUE),
      Median_SR = median(SR, na.rm = TRUE),
      SD_SR = sd(SR, na.rm = TRUE),
      Q25_SR = quantile(SR, 0.25, na.rm = TRUE),
      Q75_SR = quantile(SR, 0.75, na.rm = TRUE),
      Mean_EB_SR = mean(EB_SR, na.rm = TRUE),
      Median_EB_SR = median(EB_SR, na.rm = TRUE),
      Total_OBS = sum(Total_OBS, na.rm = TRUE),
      Total_EXP = sum(Total_EXP, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(
      Overall_SR = Total_OBS / Total_EXP,
      CV = SD_SR / Mean_SR  # Coefficient of variation
    ) %>%
    arrange(desc(Mean_SR))
  
  return(summary_stats)
}

create_formatted_table <- function(summary_stats) {
  formatted_table <- summary_stats %>%
    select(Category, Regions, Mean_SR, Median_SR, SD_SR, 
           Q25_SR, Q75_SR, Overall_SR, CV) %>%
    mutate(
      across(c(Mean_SR:CV), ~ round(.x, 3))
    ) %>%
    kable(
      col.names = c("Category", "N Regions", "Mean", "Median", "SD", 
                    "Q25", "Q75", "Overall", "CV"),
      caption = "Summary Statistics for SR by Category",
      align = c("l", rep("c", 8))
    ) %>%
    kable_styling(
      bootstrap_options = c("striped", "hover", "condensed", "responsive"),
      full_width = FALSE,
      position = "center"
    ) %>%
    add_header_above(c(" " = 2, "SR Distribution" = 5, "Combined" = 1, "Variability" = 1)) %>%
    row_spec(0, bold = TRUE, background = "#f8f9fa")
  
  return(formatted_table)
}

# Main Analysis Pipeline
# ======================

run_dm2_analysis <- function(data) {
  # Data preparation
  dm2_categories <- prepare_dm2_categories()
  indicator_to_category <- create_indicator_lookup(dm2_categories)
  
  # Categorize indicators
  data_categorized <- data %>%
    rowwise() %>%
    mutate(Category = categorize_indicator(Indicator, indicator_to_category, dm2_categories)) %>%
    ungroup() %>%
    filter(Category != "Other")
  
  # Calculate composite indicators
  composite <- data_categorized %>%
    group_by(`ID GMP`, Category) %>%
    summarize(
      Total_OBS = sum(`OBS`, na.rm = TRUE),
      Total_EXP = sum(`Exp`, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    filter(Total_EXP > 0) %>%
    mutate(SR = Total_OBS / Total_EXP)
  
  # Add confidence intervals
  ci_data <- calculate_confidence_intervals(composite$Total_OBS, composite$Total_EXP)
  composite <- bind_cols(composite, ci_data)
  
  # Add Empirical Bayes smoothing
  composite$EB_SR <- empirical_bayes_smooth(composite$Total_OBS, composite$Total_EXP)
  
  # Clean data
  composite <- composite %>%
    filter(!is.na(SR) & !is.na(CI_lower) & !is.na(CI_upper)) %>%
    mutate(Category = as.factor(Category))
  
  return(composite)
}

# Execute Analysis
# ================

# Load and process data
# df <- read_excel("your_file.xlsx")  # Uncomment and modify path as needed
df <- DM_Related_to_work_2  # Replace with your data loading

# Run analysis
composite_results <- run_dm2_analysis(df)

# Create visualizations
plot_forest <- create_sr_forest_plot(composite_results)
plot_comparison <- create_sr_comparison_plot(composite_results)
plot_distribution <- create_category_distribution_plot(composite_results)
plot_reliability <- create_reliability_plot(composite_results)

# Display plots
print(plot_forest)
print(plot_comparison)
print(plot_distribution)
print(plot_reliability)

# Create combined dashboard
dashboard <- (plot_comparison + plot_distribution) / 
  (plot_reliability + plot_forest)

dashboard <- dashboard +
  plot_annotation(
    title = "DM2 Healthcare Quality Indicators Analysis Dashboard",
    subtitle = "Comprehensive analysis of Standardized Ratios across healthcare categories",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

print(dashboard)

# Create and display summary table
summary_stats <- create_summary_table(composite_results)
formatted_table <- create_formatted_table(summary_stats)
print(formatted_table)

# Export results (using base R functions to avoid vroom issues)
write.csv(composite_results, "dm2_composite_indicators_enhanced2.csv", row.names = FALSE)
write.csv(summary_stats, "dm2_summary_statistics2.csv", row.names = FALSE)

# Save plots
ggsave("dm2_forest_plot.png", plot_forest, width = 12, height = 10, dpi = 300)
ggsave("dm2_comparison_plot.png", plot_comparison, width = 10, height = 8, dpi = 300)
ggsave("dm2_distribution_plot.png", plot_distribution, width = 10, height = 8, dpi = 300)
ggsave("dm2_dashboard.png", dashboard, width = 16, height = 12, dpi = 300)