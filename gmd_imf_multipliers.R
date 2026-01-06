# ============================================================================
# Fiscal Multiplier Estimation using IMF Narrative Shocks and Global Macro Data
# ============================================================================
#
# This script estimates fiscal multipliers using LP-IV (Local Projection-IV)
# methods with three different transformations (HBR, GK, CP) to ensure
# consistent estimation across horizons, following Jordà & Taylor (2025).
#
# Data sources:
# - IMF narrative fiscal consolidation shocks (Alesina et al., 2019)
# - Global Macro Database (globalmacrodata package)
#
# Author: [Your Name]
# Date: 2026-01-07
# ============================================================================

library(globalmacrodata)
library(readr)
library(haven)
library(dplyr)
library(tidyr)
library(ggplot2)

# ============================================================================
# STEP 1: Load and prepare data
# ============================================================================

# Load IMF narrative fiscal consolidation shocks
# Source: Adler, Allen, Ganelli and Leigh (2024) - fiscal_consolidation_v032023.dta from Jorda and Taylor (2025) replication files
agl_data <- read_dta("fiscal_consolidation_v032023.dta")

# Extract country lists from different samples in AGL data
agl_countries <- agl_data |>
    filter(AGL_sample == 1) |>  # Full AGL sample
    pull(iso) |>
    unique()

d_countries <- agl_data |>
    filter(D_sample == 1) |>  # Original 17 OECD economies in IMF sample
    pull(iso) |>
    unique()

# Define variables to extract from Global Macro Database
gmd_vars <- c(
    "rGDP",         # Real GDP (levels)
    "nGDP",         # Nominal GDP
    "deflator",     # GDP deflator
    "govdebt_GDP",  # Government debt to GDP ratio
    "govdef",       # Government deficit (nominal)
    "cbrate",       # Central bank policy rate
    "CA_GDP",       # Current account to GDP ratio
    "USDfx"         # USD exchange rate
)

gmd_version <- "2025_12"  # Version of Global Macro Database to use

# Fetch data from Global Macro Database (commented out - using cached data)
# gmd_data <- gmd(
#     version = gmd_version,
#     variables = gmd_vars,
#     country = agl_countries
# )
# gmd_data <- gmd_data |>
#     rename(iso = ISO3) |>
#     select(-id) |>
#     filter(year >= min(agl_data$year) & year <= max(agl_data$year))
# write_csv(gmd_data, "gmd_data_2025_12_agl_countries.csv")

# Load pre-fetched GMD data from saved .csv file (for reproducibility/speed)
gmd_data <- read_csv("gmd_data_2025_12_agl_countries.csv")

# Keep only the instrumental variables from AGL data
agl_data <- agl_data |>
    select(iso, year, size, tax, spend)  # size = total consolidation shock

# Merge AGL shocks with GMD macro variables
agl_gmd <- agl_data |>
    left_join(gmd_data, by = c("iso", "year"))

# ============================================================================
# STEP 2: Data cleaning
# ============================================================================

# Identify countries with no fiscal consolidation shocks at all
no_shock_isos <- agl_gmd |>
    group_by(iso) |>
    summarize(all_zero = all(size == 0)) |>
    filter(all_zero) |>
    pull(iso)

# Remove countries with no shocks (cannot estimate with IV if no variation)
agl_gmd <- agl_gmd |>
    filter(!(iso %in% no_shock_isos))

# Transform government deficit to real terms using GDP deflator
# This ensures consistency with real GDP in the regression
agl_gmd <- agl_gmd |>
    mutate(rgovdef = govdef / deflator)

# ============================================================================
# STEP 3: Create control variables
# ============================================================================

# Load the LP-IV panel estimation function
source("lp_iv_panel_fm.R")

# Define control variables that will be lagged and first-differenced
# Following Jordà & Taylor (2025), we control for lagged changes in:
control_vars_prep <- c(
    "govdebt_GDP",  # Government debt ratio (controls for initial public indebtedness)
    "cbrate",       # Central bank rate (controls for monetary policy)
    "CA_GDP",       # Current account (controls for external balance)
    "USDfx"         # Exchange rate (controls for currency movements)
)

# Create L1.D(X) = first difference of X lagged once
agl_gmd <- agl_gmd |>
    group_by(iso) |>
    arrange(year) |>
    mutate(across(all_of(control_vars_prep), list(
        L1_D = ~ (. - lag(., 1)) - (lag(., 1) - lag(., 2))  # L1.D(X) = (X_t - X_{t-1}) - (X_{t-1} - X_{t-2})
    ), .names = "{fn}_{col}")) |>
    ungroup()

# ============================================================================
# STEP 4: Define sample
# ============================================================================

# Jordà and Taylor (2025) use 16 OECD countries from 1978-2019
# jt_countries <- c("AUS", "BEL", "CAN", "DEU", "DNK", "ESP", "FIN", "FRA",
#                   "GBR", "IRL", "ITA", "JPN", "NLD", "PRT", "SWE", "USA")

# Define our sample: 10 major OECD economies of interest
sample_countries <- c(
    "USA", "GBR", "FRA", "DEU", "ITA",
    "ESP", "SWE", "CAN", "JPN", "AUS"
)

# Filter data to sample countries
agl_gmd <- agl_gmd |>
    filter(iso %in% sample_countries)

# ============================================================================
# STEP 5: Estimate fiscal multipliers using LP-IV
# ============================================================================
#
# We estimate fiscal multipliers using Local Projection IV (LP-IV) with three
# different transformations to ensure consistent estimation across horizons:
#
# 1. HBR (Hall-Barro-Redlick): Normalizes by lagged outcome level
#    - Y^h = (Y_{t+h} - Y_{t-1}) / Y_{t-1}
#
# 2. GK (Gordon-Krenn): Normalizes by trend outcome
#    - Y^h = Y_{t+h}/Ŷ_{t+h} - Y_{t-1}/Ŷ_{t-1}
#
# 3. CP (Canova-Pappa): Cumulative normalization
#    - Y^h = (Y_{t+h} - (h+1)·Y_{t-1}) / Y_{t-1}
#
# Key specification choices:
# - Outcome: rGDP (real GDP)
# - Treatment: rgovdef (real government deficit)
# - Instrument: size (IMF narrative fiscal consolidation shock)
# - Controls: 2 lags of transformed outcome and treatment
#            + lagged first differences of debt, interest rate, CA, and FX
#
# The LP-IV approach estimates the impulse response at each horizon h by:
# S_h Y_t = β_h · S_h X_t + controls + ε_t
# where S_h denotes cumulative sum from 0 to h, and X is instrumented by size
#
# ============================================================================

# Define transformations to estimate
transformations <- c("HBR", "GK", "CP")

# Parameters for LP-IV estimation
gk_trend <- "hp"           # Use HP filter for GK trend (alternative: "poly")
max_h <- 4                 # Maximum horizon (0 to 4 years)
outcome_nlags <- 2         # Number of lags of transformed outcome
treatment_nlags <- 2       # Number of lags of transformed treatment

# Construct vector of additional control variables
# Format: L1_D_varname for each variable in control_vars_prep
additional_controls <- c()
for (var in control_vars_prep) {
    additional_controls <- c(additional_controls, paste0("L1_D_", var))
}
# Result: c("L1_D_govdebt_GDP", "L1_D_cbrate", "L1_D_CA_GDP", "L1_D_USDfx")


# Initialize data frame to store results across all transformations
results_transformed <- data.frame(
    spec = character(),      # Transformation type (HBR, GK, or CP)
    h = integer(),           # Horizon (0 to max_h)
    m = numeric(),           # Estimated multiplier (coefficient on treatment)
    se = numeric(),          # Standard error (cluster-robust or normal)
    upper_ci = numeric(),    # Upper 95% confidence interval
    lower_ci = numeric(),    # Lower 95% confidence interval
    n_obs = integer()        # Number of observations used in regression
)

# Initialize list to store full model objects for each transformation
all_models <- list()

# Loop through each transformation and estimate LP-IV models
for (trf in transformations) {
    cat("\n======================================\n")
    cat("Running", trf, "transformation\n")
    cat("======================================\n")

    # Call LP-IV panel function with specified parameters
    # This function will:
    # 1. Apply the transformation to outcome and treatment variables
    # 2. Create cumulative sums across horizons
    # 3. Generate lagged transformed variables for controls
    # 4. Run IV regressions with panel fixed effects (within transformation)
    # 5. Compute cluster-robust standard errors (if multiple groups)
    results <- lp_iv_panel_fm(
        df = agl_gmd,                          # Input data
        outcome = "rGDP",                      # Outcome variable
        treatment = "rgovdef",                 # Treatment variable
        transformation = trf,                  # Transformation type (HBR/GK/CP)
        max_h = max_h,                         # Maximum horizon
        gk_trend = gk_trend,                   # Trend method for GK
        outcome_nlags = outcome_nlags,         # Lags of outcome to control for
        treatment_nlags = treatment_nlags,     # Lags of treatment to control for
        instrument = "size",                   # IV: narrative consolidation shock
        group_var = "iso",                     # Panel identifier (country)
        time_var = "year",                     # Time identifier
        additional_controls = additional_controls,  # Extra controls
        hp_lambda = 400,                       # HP filter smoothing (quarterly freq)
        poly_degree = 3                        # Polynomial degree (if using poly trend)
    )

    # Store complete results object (includes models, transformed data, etc.)
    all_models[[trf]] <- results

    # Label results with transformation name
    results$results$spec <- trf

    # Append to combined results data frame
    results_transformed <- rbind(
        results_transformed,
        results$results[, c("spec", "h", "m", "se", "upper_ci", "lower_ci", "n_obs")]
    )
}

# ============================================================================
# STEP 6: Visualize results
# ============================================================================

# Create two plots:
# 1. With confidence intervals (ribbons) to show uncertainty
# 2. Without confidence intervals for cleaner comparison across transformations

# Plot 1: Fiscal multiplier estimates with 95% confidence intervals
p1 <- ggplot(results_transformed, aes(x = h, y = m, colour = spec, group = spec)) +
    geom_line() +
    geom_point() +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = spec), alpha = 0.2, colour = NA) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    labs(title = "LP-IV Fiscal Multipliers: Different Transformations (with CIs)",
         subtitle = "Outcome: rGDP, Treatment: rgovdef, Instrument: size",
         x = "Horizon",
         y = "Multiplier",
         colour = "Transformation",
         fill = "Transformation") +
    theme_minimal()
print(p1)

# Plot 2: Fiscal multiplier estimates without confidence intervals
p2 <- ggplot(results_transformed, aes(x = h, y = m, colour = spec, group = spec)) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    labs(title = "LP-IV Fiscal Multipliers: Different Transformations (without CIs)",
         subtitle = "Outcome: rGDP, Treatment: rgovdef, Instrument: size",
         x = "Horizon",
         y = "Multiplier",
         colour = "Transformation") +
    theme_minimal()
print(p2)

# ============================================================================
# STEP 7: Display and export results
# ============================================================================

cat("\n======================================\n")
cat("Results Summary\n")
cat("======================================\n")
cat("Fiscal multiplier estimates by transformation and horizon\n")
cat("Specification: Outcome = rGDP, Treatment = rgovdef, IV = size\n")
cat("======================================\n\n")
print(results_transformed)

cat("\n======================================\n")
cat("Analysis complete!\n")
cat("======================================\n")

# ============================================================================
# STEP 8: Optional - Access detailed results
# ============================================================================

# Access individual transformation results:
# Example: Access HBR model for horizon 2
# all_models$HBR$models[[3]]$model  # horizon 2 is index 3 (0-indexed + 1)
# all_models$HBR$models[[3]]$coef_test

# Access transformed data for further analysis
data_hbr <- all_models$HBR$data_transformed
data_gk <- all_models$GK$data_transformed
data_cp <- all_models$CP$data_transformed

# ============================================================================
# STEP 9: Additional diagnostics and exploration
# ============================================================================

# Combine transformed variables from all three transformations
# Keep only transformation-specific variables (drop raw and control variables)
data_hbr_clean <- data_hbr |> select(iso, year, starts_with("HBR"))
data_gk_clean <- data_gk |> select(iso, year, starts_with("GK"))
data_cp_clean <- data_cp |> select(iso, year, starts_with("CP"))
data_combined <- data_hbr_clean |>
    left_join(data_gk_clean, by = c("iso", "year")) |>
    left_join(data_cp_clean, by = c("iso", "year"))
# Add instrument (size) back to combined data for diagnostic plots
data_combined <- data_combined |>
    left_join(agl_gmd |> select(iso, year, size), by = c("iso", "year"))

# ============================================================================
# Diagnostic plots
# ============================================================================

# Plot transformed treatment variable over time by country
p3 <- ggplot(data_combined, aes(x = year, y = GK0rgovdef, colour = iso, group = iso)) +
    geom_line() +
    labs(title = "HBR Transformed rGDP over Time by Country",
         x = "Year",
         y = "HBR Transformed rGDP",
         colour = "Country") +
    theme_minimal()

print(p3)

# Demean all transformed variables by country (within transformation)
# This replicates the panel fixed effects transformation used in estimation
data_combined_demeaned <- data_combined |>
    group_by(iso) |>
    mutate(across(starts_with(c("HBR", "GK", "CP")), ~ . - mean(., na.rm = TRUE))) |>
    ungroup()   

# Scatter plot: demeaned outcome vs demeaned treatment at horizon 0
# This visualizes the first-stage relationship used in IV estimation
transformation_to_plot <- "HBR"
horizon_to_plot <- 0
outcome_var <- paste0(transformation_to_plot, horizon_to_plot, "rGDP")
treatment_var <- paste0(transformation_to_plot, horizon_to_plot, "rgovdef")
p4 <- ggplot(
    data_combined_demeaned,
    aes_string(x = treatment_var, y = outcome_var, colour = "iso", group = "iso")
) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(
        title = paste0(transformation_to_plot, " Transformed rGDP vs rgovdef at Horizon ", horizon_to_plot))

print(p4)

