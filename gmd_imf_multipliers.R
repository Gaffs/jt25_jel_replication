library(globalmacrodata)
library(readr)
library(haven)
library(dplyr)
library(tidyr)
library(ggplot2)

# Load IMF narrative shocks and cyclically adjusted primary balance data
agl_data <- read_dta("fiscal_consolidation_v032023.dta")


agl_countries <- agl_data |> 
    filter(AGL_sample == 1) |> 
    pull(iso) |> 
    unique()

d_countries <- agl_data |> 
    filter(D_sample == 1) |> 
    pull(iso) |>
    unique()

gmd_vars <- c(
    "rGDP", "nGDP", "deflator", "govdebt_GDP", "govdef", "cbrate"
)

gmd_version <- "2025_12"

gmd_data <- gmd(
    version = gmd_version,
    variables = gmd_vars,
    country = agl_countries
)

agl_data <- agl_data |>
    select(iso, year, size, tax, spend)

gmd_data <- gmd_data |> rename(iso = ISO3) |> select(-id)

agl_gmd <- agl_data |>
    left_join(gmd_data, by = c("iso", "year"))

agl_gmd |> filter(size == 0)

# identify isos with size == 0 for all observations
no_shock_isos <- agl_gmd |>
    group_by(iso) |>
    summarize(all_zero = all(size == 0)) |>
    filter(all_zero) |>
    pull(iso)

# Remove countries with no shocks
agl_gmd <- agl_gmd |>
    filter(!(iso %in% no_shock_isos))

# Transform govdef to real terms using deflator
agl_gmd <- agl_gmd |>
    mutate(rgovdef = govdef / deflator)

# Load the LP-IV function
source("lp_iv_panel_fm.R")

# Create additional control variable: one lag of first difference of govdebt_GDP
agl_gmd <- agl_gmd |>
    group_by(iso) |>
    mutate(
        D_govdebt_GDP = govdebt_GDP - dplyr::lag(govdebt_GDP, 1),
        L1_D_govdebt_GDP = dplyr::lag(D_govdebt_GDP, 1)
    ) |>
    ungroup()

sample_countries <- c(
    "USA", "GBR", "FRA", "DEU", "ITA", "ESP", "SWE", "CAN", "JPN", "AUS"
)

# Filter agl_gmd to only include sample countries
agl_gmd <- agl_gmd |>
    filter(iso %in% sample_countries)

# ============================================================================
# Run LP-IV regressions for each transformation using lp_iv_panel_fm()
# ============================================================================

# Set parameters
transformations <- c("HBR", "GK", "CP")
gk_trend <- "hp"
max_h <- 4
outcome_nlags <- 2
treatment_nlags <- 2
additional_controls <- c("L1_D_govdebt_GDP")

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
        df = agl_gmd,
        outcome = "rGDP",
        treatment = "rgovdef",
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

# Plot with confidence intervals
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

# Plot without confidence intervals for clarity
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
# Display results summary
# ============================================================================

cat("\n======================================\n")
cat("Results Summary\n")
cat("======================================\n")
print(results_transformed)

cat("\n======================================\n")
cat("Analysis complete!\n")
cat("======================================\n")

# Optional: Access individual transformation results
# Example: Access HBR model for horizon 2
# all_models$HBR$models[[3]]$model  # horizon 2 is index 3 (0-indexed + 1)
# all_models$HBR$models[[3]]$coef_test

# Example: Access transformed data for further analysis
# data_hbr <- all_models$HBR$data_transformed
# data_gk <- all_models$GK$data_transformed
# data_cp <- all_models$CP$data_transformed

