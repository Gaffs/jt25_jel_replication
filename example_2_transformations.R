# Extends Jorda and Taylor (2025) Example 2, Figure 2b
# Author: AG
# Requires .dta files bundled with their Stata replication files:
# - fiscal_consolidation_v032023.dta
# - JSTdatasetR6.dta
# - dcapb.dta
# See https://github.com/ojorda/JEL-Code

# ============================================================================

# Clear workspace
rm(list = ls())

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Source original replication script and LP-IV function
source("example_2_replication.R")
source("lp_iv_panel_fm.R")

# ============================================================================
# EXTENSION: CORRECTING TRANSFORMATIONS OF TREATMENT AND OUTCOME
# ============================================================================

## Prepare data with additional variables for consistent fiscal multiplier calculation

# Create real GDP from population, real GDP per capita and nominal GDP
data_clean <- data %>%
  group_by(iso) %>%
  mutate(rgdpraw = rgdpbarro * pop,
         rgdpraw_2005 = rgdpraw[year == 2005],
         rgdp_idx = (rgdpraw / rgdpraw_2005) * 100,
         gdp_2005 = gdp[year == 2005],
         gdp_idx = (gdp / gdp_2005) * 100,
         pgdp = gdp_idx / rgdp_idx,
         rgdp = gdp / pgdp) %>%
  ungroup() %>%
  select(-rgdpraw_2005, -rgdpraw, -gdp_2005)

# Create real fiscal balance
data_clean <- data_clean %>%
  mutate(fb = revenue - expenditure,
         rfb = fb / pgdp)

# Create additional control variables
data_clean <- data_clean %>%
  group_by(iso) %>%
  mutate(
    L1ygap = dplyr::lag(ygap, 1),
    L1Ddebtgdp = dplyr::lag(debtgdp, 1) - dplyr::lag(debtgdp, 2)
  ) %>%
  ungroup()

# ============================================================================
# Run LP-IV regressions for each transformation using lp_iv_panel_fm()
# ============================================================================

# Set parameters
transformations <- c("HBR", "GK", "CP")
gk_trend <- "hp"  # options: "poly" or "hp"
max_h <- 4
outcome_nlags <- 2
treatment_nlags <- 2
additional_controls <- c("L1ygap", "L1Ddebtgdp")

# Initialize results storage
results_transformed <- data.frame(
  spec = character(),
  h = integer(),
  m = numeric(),
  se = numeric(),
  upper_ci = numeric(),
  lower_ci = numeric(),
  n_obs = integer()
)

# Store all model outputs
all_models <- list()

# Run LP-IV for each transformation
for (trf in transformations) {
  cat("\n======================================\n")
  cat("Running", trf, "transformation\n")
  cat("======================================\n")

  results <- lp_iv_panel_fm(
    df = data_clean,
    outcome = "rgdp",
    treatment = "rfb",
    transformation = trf,
    max_h = max_h,
    gk_trend = gk_trend,
    outcome_nlags = outcome_nlags,
    treatment_nlags = treatment_nlags,
    instrument = "size",
    group_var = "iso",
    time_var = "year",
    additional_controls = additional_controls,
    hp_lambda = 400,
    poly_degree = 3
  )

  # Store models
  all_models[[trf]] <- results

  # Add specification name to results
  results$results$spec <- trf

  # Append to combined results
  results_transformed <- rbind(
    results_transformed,
    results$results[, c("spec", "h", "m", "se", "upper_ci", "lower_ci", "n_obs")]
  )
}

# ============================================================================
# Visualise results
# ============================================================================

# Plot transformed results only
p1 <- ggplot(results_transformed, aes(x = h, y = m, colour = spec, group = spec)) +
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = spec), alpha = 0.2, colour = NA) +
  labs(title = "LP-IV Responses: Different Transformations",
       x = "Horizon",
       y = "Response",
       colour = "Transformation",
       fill = "Transformation") +
  theme_minimal()
print(p1)

# Combine with baseline results from example_2_replication.R for comparison
results_combined <- results_transformed |> 
    select(colnames(results_baseline)) |> 
    rbind(results_baseline)

# Plot combined results with confidence intervals
p2 <- ggplot(results_combined, aes(x = h, y = m, colour = spec, group = spec)) +
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = spec), alpha = 0.2, colour = NA) +
  labs(title = "LP-IV Responses: Baseline vs Different Transformations",
       x = "Horizon",
       y = "Response",
       colour = "Specification",
       fill = "Specification") +
  theme_minimal()
print(p2)

# Plot combined results without ribbons for clarity
p3 <- ggplot(results_combined, aes(x = h, y = m, colour = spec, group = spec)) +
  geom_line() +
  geom_point() +
  labs(title = "LP-IV Responses: Baseline vs Different Transformations",
       x = "Horizon",
       y = "Response",
       colour = "Specification") +
  theme_minimal()
print(p3)

# ============================================================================
# Optional: Access individual transformation results
# ============================================================================

# Example: Access HBR model for horizon 2
# all_models$HBR$models[[3]]$model  # horizon 2 is index 3 (0-indexed + 1)
# all_models$HBR$models[[3]]$coef_test

# Example: Access transformed data for further analysis
# data_hbr <- all_models$HBR$data_transformed
# data_gk <- all_models$GK$data_transformed
# data_cp <- all_models$CP$data_transformed

cat("\n======================================\n")
cat("Analysis complete!\n")
cat("======================================\n")

# END OF SCRIPT
