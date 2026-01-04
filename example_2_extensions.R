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

# Source original replication script
source("example_2_replication.R")

# ============================================================================
# EXTENSION 1: CLEANING OUTLIERS
# ============================================================================

# Data wrangling =================================================================

# Demean data by country for fixed effects approximation for treatment and control 
# of baseline specification
demeaned_data <- data %>%
  group_by(iso) %>%
  mutate(across(starts_with("S"), ~ . - mean(., na.rm = TRUE))) %>%
  ungroup()

# Plot horizon 0 outcome vs treatment to identify outliers
# (remember that horizon 0 outliers are cumulated to later horizons)
p <- ggplot(data, aes(x = SdCAPB1, y = S1y, color = iso)) +
  geom_point() +
  labs(title = "Horizon 0 outcome vs treatment",
       x = "Demeaned SdCAPB0",
       y = "Demeaned S0y") +
  theme_minimal()

print(p)

# Outlier countries based on visual inspection of horizon 0
# We also drop Germany because reunification causes large jumps in GDP levels as well as fiscal spending
# (which would affect later use of GDP in levels)
outlier_countries <- c("IRL", "FIN", "DEU")

# Remove outlier countries
data_clean <- data %>%
  filter(!(iso %in% outlier_countries))

# Re-run the LP-IV regressions on cleaned data =======================================

# Initialize results storage
results_clean <- data.frame(
  spec = character(),
  h = integer(),
  m = numeric(),
  upper_ci = numeric(),
  lower_ci = numeric()
)

lp_models <- list()

# Run LP-IV for each horizon (0 to 4)
# In Stata: xi: xtivreg2 S`h'y (SdCAPB`h' = size) _x* , fe cluster(iso)
for (h in 0:4) {
  cat("  Horizon", h, "...\n")

  # Prepare data
  outcome_var <- paste0("S", h, "y")
  treatment_var <- paste0("SdCAPB", h)

  # Select complete cases
  reg_data <- data_clean %>%
    select(all_of(outcome_var), all_of(treatment_var), size,
           all_of(controls), iso) %>%
    na.omit()

  cat("    Observations:", nrow(reg_data), "\n")

  # Demean for fixed effects (i.e. the within transformation)
  reg_data_demeaned <- reg_data %>%
    group_by(iso) %>%
    mutate(across(c(all_of(outcome_var), all_of(treatment_var), size, all_of(controls)),
                  ~ . - mean(., na.rm = TRUE))) %>%
    ungroup()

  # Run IV regression

  # Create string summing control variables for formula
  control_vars_str <- paste(controls, collapse = " + ")

  # Formula: outcome ~ treatment + controls | instrument + controls
  formula_iv <- as.formula(paste0(
    outcome_var, " ~ ", treatment_var, " + ", control_vars_str, " | ",
    "size + ", control_vars_str
  ))

  model <- ivreg(formula_iv, data = reg_data_demeaned)

  # Cluster-robust standard errors by iso
  vcov_cluster <- vcovCL(model, cluster = reg_data$iso, type = "HC1")
  coef_test <- coeftest(model, vcov = vcov_cluster)

  # Store model
  lp_models[[h + 1]] <- list(model = model, coef_test = coef_test)

  # Extract coefficient for treatment variable
  # Try different possible names
  coef_names <- rownames(coef_test)
  treatment_coef_name <- treatment_var  # Default

  if (!(treatment_var %in% coef_names)) {
    # Try with backticks
    if (paste0("`", treatment_var, "`") %in% coef_names) {
      treatment_coef_name <- paste0("`", treatment_var, "`")
    }
  }

  b_coef <- coef_test[treatment_coef_name, "Estimate"]
  se_coef <- coef_test[treatment_coef_name, "Std. Error"]

  cat("    Coefficient:", round(b_coef, 4), ", SE:", round(se_coef, 4), "\n")

  # Store results
  results_clean <- rbind(results_clean, data.frame(
    spec = "no_outliers",
    h = h,
    m = b_coef,
    upper_ci = b_coef + 1.96 * se_coef,
    lower_ci = b_coef - 1.96 * se_coef
  ))
}

# Combine results_baseline and results_clean for plotting
results_combined <- rbind(results_baseline, results_clean)

# Plot combined results ========================================================

# Removing outliers and Germany almost halves the multiplier on impact
# Estimates are also a lot more precise
p <- ggplot(results_combined, aes(x = h, y = m, color = spec, group = spec)) +
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = spec), alpha = 0.2, color = NA) +
  labs(title = "LP-IV Responses: Baseline vs Cleaned Data (No Outliers)",
       x = "Horizon",
       y = "Response") +
  theme_minimal() +
  scale_color_manual(values = c("jt25_baseline" = "blue", "no_outliers" = "red")) +
  scale_fill_manual(values = c("jt25_baseline" = "blue", "no_outliers" = "red"))

print(p)

# ============================================================================
# EXTENSION 2: DROPPING CONTROLS
# ============================================================================

# The output gap and debt controls seem unnecessary given the IV strategy

# Select controls
controls <- c("L1Dy", "L2Dy", "L1dCAPB", "L2dCAPB")

# Re-run the LP-IV regressions without output gap and debt controls =============
# Initialize results storage
results_nocontrols <- data.frame(
  spec = character(),
  h = integer(),
  m = numeric(),
  upper_ci = numeric(),
  lower_ci = numeric()
)

lp_models <- list()

# Run LP-IV for each horizon (0 to 4)
# In Stata: xi: xtivreg2 S`h'y (SdCAPB`h' = size) _x* , fe cluster(iso)
for (h in 0:4) {
  cat("  Horizon", h, "...\n")

  # Prepare data
  outcome_var <- paste0("S", h, "y")
  treatment_var <- paste0("SdCAPB", h)

  # Select complete cases
  reg_data <- data_clean %>%
    select(all_of(outcome_var), all_of(treatment_var), size,
           all_of(controls), iso) %>%
    na.omit()

  cat("    Observations:", nrow(reg_data), "\n")

  # Demean for fixed effects (i.e. the within transformation)
  reg_data_demeaned <- reg_data %>%
    group_by(iso) %>%
    mutate(across(c(all_of(outcome_var), all_of(treatment_var), size, all_of(controls)),
                  ~ . - mean(., na.rm = TRUE))) %>%
    ungroup()

  # Run IV regression

  # Create string summing control variables for formula
  control_vars_str <- paste(controls, collapse = " + ")

  # Formula: outcome ~ treatment + controls | instrument + controls
  formula_iv <- as.formula(paste0(
    outcome_var, " ~ ", treatment_var, " + ", control_vars_str, " | ",
    "size + ", control_vars_str
  ))

  model <- ivreg(formula_iv, data = reg_data_demeaned)

  # Cluster-robust standard errors by iso
  vcov_cluster <- vcovCL(model, cluster = reg_data$iso, type = "HC1")
  coef_test <- coeftest(model, vcov = vcov_cluster)

  # Store model
  lp_models[[h + 1]] <- list(model = model, coef_test = coef_test)

  # Extract coefficient for treatment variable
  # Try different possible names
  coef_names <- rownames(coef_test)
  treatment_coef_name <- treatment_var  # Default

  if (!(treatment_var %in% coef_names)) {
    # Try with backticks
    if (paste0("`", treatment_var, "`") %in% coef_names) {
      treatment_coef_name <- paste0("`", treatment_var, "`")
    }
  }

  b_coef <- coef_test[treatment_coef_name, "Estimate"]
  se_coef <- coef_test[treatment_coef_name, "Std. Error"]

  cat("    Coefficient:", round(b_coef, 4), ", SE:", round(se_coef, 4), "\n")

  # Store results
  results_nocontrols <- rbind(results_nocontrols, data.frame(
    spec = "no_controls",
    h = h,
    m = b_coef,
    upper_ci = b_coef + 1.96 * se_coef,
    lower_ci = b_coef - 1.96 * se_coef
  ))
}
# Combine results_clean and results_nocontrols for plotting
results_combined2 <- rbind(results_clean, results_nocontrols)
# Plot combined results ========================================================

# Dropping controls reduces the estimated multiplier further throughout the horizon
# Estimates are more precise at earlier horizons
p <- ggplot(results_combined2, aes(x = h, y = m, color = spec, group = spec)) +
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = spec), alpha = 0.2, color = NA) +
  labs(title = "LP-IV Responses: Cleaned Data vs No Controls",
       x = "Horizon",
       y = "Response") +
  theme_minimal() +
  scale_color_manual(values = c("no_outliers" = "red", "no_controls" = "green")) +
  scale_fill_manual(values = c("no_outliers" = "red", "no_controls" = "green"))

print(p)

# ============================================================================
# EXTENSION 3: CORRECTING TRANSFORMATIONS OF TREATMENT AND OUTCOME
# ============================================================================

## Create additional variables needed for consistent fiscal multiplier calculation

# Create real GDP from population, real GDP per capita and nominal GDP
data_clean <- data_clean %>%
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

# Fit polynomial to GDP to derive trend for GK transformation below
data_clean$rgdptrend <- NA
countries <- unique(data_clean$iso)
for (c in countries) {
  country_data <- data_clean %>% filter(iso == c) %>% arrange(year)
  y_vals <- log(country_data$rgdp)
  if (sum(!is.na(y_vals)) > 0) {
    x_vals <- country_data$year[!is.na(y_vals)]
    poly_fit <- lm(y_vals[!is.na(y_vals)] ~ poly(x_vals, 5, raw = TRUE))
    full_trend <- predict(poly_fit, newdata = data.frame(x_vals = country_data))
    data_clean$rgdptrend[data_clean$iso == c] <- exp(full_trend)
  }
}

# Create transformations for consistent fiscal multiplier calculation
# vector of variables to transform
vars_to_transform <- c("rgdp", "rfb")
for (var in vars_to_transform) {
  for (h in 0:10) {
    data_clean <- data_clean %>%
      group_by(iso) %>%
      mutate(
        # Hall-Barro-Redlick (HBR) transformation: HBR_h(x) = (x_{t+h} - x_{t-1}) / rgdp_{t-1}
        !!paste0("HBR", h, var) := 
        (dplyr::lead(!!sym(var), h) - dplyr::lag(!!sym(var), 1)) / 
        dplyr::lag(rgdp, 1),
        # Gordon-Krenn (GK) transformation: GK_h(x) = x_{t+h} / rgdptrend_{t+h})
        !!paste0("GK", h, var) := 
        (dplyr::lead(!!sym(var), h)) / dplyr::lead(rgdptrend, h) - 
        (dplyr::lag(!!sym(var), 1)) / dplyr::lag(rgdptrend, 1),
        # Canova-Pappa (CP) transformation: CP_h(x) = (x_{t+h} - (h+1) * x_{t-1}) / rgdp_{t-1}
        !!paste0("CP", h, var) := 
        (dplyr::lead(!!sym(var), h) - (h + 1) * dplyr::lag(!!sym(var), 1)) / 
        dplyr::lag(rgdp, 1)
      ) %>%
      ungroup()
  }
}

# Create horizon-wise cumulative sums for each transformation for each variable
transformations <- c("HBR", "GK", "CP")

for (trf in transformations) {
  for (var in vars_to_transform) {
    for (h in 0:10) {
      cols_to_sum <- paste0(trf, 0:h, var)
      data_clean <- data_clean %>%
        mutate(!!paste0("S", trf, h, var) := rowSums(select(., all_of(cols_to_sum)), na.rm = FALSE))
    }
  }
}

# Create new controls for new transformed variables
data_clean <- data_clean %>%
  group_by(iso) %>%
  mutate(
    rgdpgr = rgdp / dplyr::lag(rgdp, 1) - 1,
    L1rgdpgr = dplyr::lag(rgdpgr, 1),
    L2rgdpgr = dplyr::lag(rgdpgr, 2),
    L1Dfbgdp = dplyr::lag(rfb, 1) / dplyr::lag(rgdp, 1) - dplyr::lag(rfb, 2) / dplyr::lag(rgdp, 2),
    L2Dfbgdp = dplyr::lag(rfb, 2) / dplyr::lag(rgdp, 2) - dplyr::lag(rfb, 3) / dplyr::lag(rgdp, 3),
    L1GK0rgdp = dplyr::lag(GK0rgdp, 1),
    L2GK0rgdp = dplyr::lag(GK0rgdp, 2),
    L1GK0rfb = dplyr::lag(GK0rfb, 1),
    L2GK0rfb = dplyr::lag(GK0rfb, 2)
  ) %>%
  ungroup()

# Re-run the LP-IV regressions for each transformation =================

# Initialize results storage
results_transformed <- data.frame(
  spec = character(),
  h = integer(),
  m = numeric(),
  upper_ci = numeric(),
  lower_ci = numeric()
)

for (trf in transformations) {

lp_models <- list()

# Run LP-IV for each horizon (0 to 4)
# In Stata: xi: xtivreg2 S`h'y (SdCAPB`h' = size) _x* , fe cluster(iso)
for (h in 0:4) {
  cat("  Horizon", h, "...\n")

  # Prepare data
  outcome_var <- paste0("S", trf, h, "rgdp")
  treatment_var <- paste0("S", trf, h, "rfb")

  # Select controls based on transformation
  if (trf == "GK") {
    controls <- c("L1GK0rgdp", "L2GK0rgdp", "L1GK0rfb", "L2GK0rfb")
  } else {
    controls <- c("L1rgdpgr", "L2rgdpgr", "L1Dfbgdp", "L2Dfbgdp")
  }
  # Select complete cases
  reg_data <- data_clean %>%
    select(all_of(outcome_var), all_of(treatment_var), size,
           all_of(controls), iso) %>%
    na.omit()

  cat("    Observations:", nrow(reg_data), "\n")

  # Demean for fixed effects (i.e. the within transformation)
  reg_data_demeaned <- reg_data %>%
    group_by(iso) %>%
    mutate(across(c(all_of(outcome_var), all_of(treatment_var), size, all_of(controls)),
                  ~ . - mean(., na.rm = TRUE))) %>%
    ungroup()

  # Run IV regression

  # Create string summing control variables for formula
  control_vars_str <- paste(controls, collapse = " + ")

  # Formula: outcome ~ treatment + controls | instrument + controls
  formula_iv <- as.formula(paste0(
    outcome_var, " ~ ", treatment_var, " + ", control_vars_str, " | ",
    "size + ", control_vars_str
  ))

  model <- ivreg(formula_iv, data = reg_data_demeaned)

  # Cluster-robust standard errors by iso
  vcov_cluster <- vcovCL(model, cluster = reg_data$iso, type = "HC1")
  coef_test <- coeftest(model, vcov = vcov_cluster)

  # Store model
  lp_models[[h + 1]] <- list(model = model, coef_test = coef_test)

  # Extract coefficient for treatment variable
  # Try different possible names
  coef_names <- rownames(coef_test)
  treatment_coef_name <- treatment_var  # Default

  if (!(treatment_var %in% coef_names)) {
    # Try with backticks
    if (paste0("`", treatment_var, "`") %in% coef_names) {
      treatment_coef_name <- paste0("`", treatment_var, "`")
    }
  }

  b_coef <- coef_test[treatment_coef_name, "Estimate"]
  se_coef <- coef_test[treatment_coef_name, "Std. Error"]

  cat("    Coefficient:", round(b_coef, 4), ", SE:", round(se_coef, 4), "\n")

  # Store results
  results_transformed <- rbind(results_transformed, data.frame(
    spec = trf,
    h = h,
    m = b_coef,
    upper_ci = b_coef + 1.96 * se_coef,
    lower_ci = b_coef - 1.96 * se_coef
  ))
}

}

# Compare new transformed results
p <- ggplot(results_transformed, aes(x = h, y = m, color = spec, group = spec)) +
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = spec), alpha = 0.2, color = NA) +
  labs(title = "LP-IV Responses: Different Transformations",
       x = "Horizon",
       y = "Response") +
  theme_minimal()
print(p)

# Combine results_transformed and results_nocontrols for plotting
results_combined3 <- rbind(results_nocontrols, results_transformed)

# Plot combined results ========================================================
p <- ggplot(results_combined3, aes(x = h, y = m, color = spec, group = spec)) +
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = spec), alpha = 0.2, color = NA) +
  labs(title = "LP-IV Responses: No Controls vs Different Transformations",
       x = "Horizon",
       y = "Response") +
  theme_minimal() 
print(p) 

# Plot combined results without ribbons for clarity
p <- ggplot(results_combined3, aes(x = h, y = m, color = spec, group = spec)) +
  geom_line() +
  geom_point() +
  labs(title = "LP-IV Responses: No Controls vs Different Transformations",
       x = "Horizon",
       y = "Response") +
  theme_minimal() 
print(p)

# END OF SCRIPT