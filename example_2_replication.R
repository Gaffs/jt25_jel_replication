# Replication of Jorda and Taylor (2025) Example 2, Figure 2b
# Author: AG
# Requires .dta files bundled with their Stata replication files:
# - fiscal_consolidation_v032023.dta
# - JSTdatasetR6.dta
# - dcapb.dta
# See https://github.com/ojorda/JEL-Code

# ============================================================================

# Clear workspace
rm(list = ls())

# Load required libraries

library(AER)          # For ivreg (IV regression)
library(lmtest)       # For coeftest
library(sandwich)     # For robust standard errors
library(mFilter)      # For HP filter
library(car)          # For linearHypothesis
library(readr)        # For reading/writing CSV files
library(haven)        # For reading Stata files
library(dplyr)        # For data manipulation
library(tidyr)        # For data reshaping
library(ggplot2)      # For plotting

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ============================================================================
# STEP 1: LOAD AND MERGE DATA
# ============================================================================
cat("Step 1: Loading and merging data...\n")

fiscal_data <- read_dta("fiscal_consolidation_v032023.dta")
jst_data <- read_dta("JSTdatasetR6.dta")
dcapb_data <- read_dta("dcapb.dta")

# Merge datasets (inner join, equivalent to Stata's keep(3) nogen)
data <- fiscal_data %>%
  inner_join(jst_data, by = c("ifs", "year")) %>%
  inner_join(dcapb_data, by = c("ifs", "year"))

# Clean up duplicate columns
data <- data %>%
  select(-any_of(c("iso.y", "wdicode", "country.x", "country.y"))) %>%
  rename(iso = iso.x)

# Create dCAPB variable
data$dCAPB <- data$dnlgxqa

cat("  Data loaded:", nrow(data), "observations\n")
cat("  Countries:", length(unique(data$iso)), "\n\n")

# ============================================================================
# STEP 2: CREATE OUTCOME VARIABLE AND LONG DIFFERENCES
# ============================================================================
cat("Step 2: Creating outcome variables...\n")

# Sort by country and year (equivalent to xtset)
data <- data %>% arrange(iso, year)

# Outcome variable: y = log(rgdpbarro) * 100
data <- data %>%
  group_by(iso) %>%
  mutate(
    y = log(rgdpbarro) * 100,
    Dy = y - dplyr::lag(y, 1)
  ) %>%
  ungroup()

# Create long differences: D_h(y) = y_{t+h} - y_{t-1}
# In Stata: gen D`h'y = f`h'.y - l.y
for (h in 0:10) {
  data <- data %>%
    group_by(iso) %>%
    mutate(!!paste0("D", h, "y") := dplyr::lead(y, h) - dplyr::lag(y, 1)) %>%
    ungroup()
}

# Create sums of long differences: S_h(y) = sum(D_0(y), ..., D_h(y))
# In Stata: egen S`h'y = rowtotal(D0y-D`h'y)
for (h in 0:10) {
  cols_to_sum <- paste0("D", 0:h, "y")
  data <- data %>%
    mutate(!!paste0("S", h, "y") := rowSums(select(., all_of(cols_to_sum)), na.rm = FALSE))
}

cat("  Outcome variables created\n\n")

# ============================================================================
# STEP 3: CREATE TREATMENT VARIABLES
# ============================================================================
cat("Step 3: Creating treatment variables...\n")

# Create cumulative forward dCAPB
# In Stata:
#   if h==0: gen dCAPB`h' = dCAPB
#   if h>0:  gen dCAPB`h' = f`h'.dCAPB + dCAPB`j' where j=h-1

data$dCAPB0 <- data$dCAPB

for (h in 1:10) {
  j <- h - 1
  data <- data %>%
    group_by(iso) %>%
    mutate(!!paste0("dCAPB", h) := dplyr::lead(dCAPB, h) + !!sym(paste0("dCAPB", j))) %>%
    ungroup()
}

# Create sums of dCAPB: SdCAPB_h = sum(dCAPB0, ..., dCAPB_h)
# In Stata: egen SdCAPB`h' = rowtotal(dCAPB0-dCAPB`h')
for (h in 0:10) {
  cols_to_sum <- paste0("dCAPB", 0:h)
  data <- data %>%
    mutate(!!paste0("SdCAPB", h) := rowSums(select(., all_of(cols_to_sum)), na.rm = FALSE))
}

cat("  Treatment variables created\n\n")

# ============================================================================
# STEP 4: HP FILTER FOR OUTPUT GAP
# ============================================================================
cat("Step 4: Applying HP filter for output gap...\n")

data$ygap <- NA

countries <- unique(data$iso)
for (c in countries) {
  country_data <- data %>% filter(iso == c) %>% arrange(year)
  y_vals <- country_data$y
  if (sum(!is.na(y_vals)) > 0) {
    # Adjust for NAs (necessary for Irish data missing final observation)
    valid_indices <- which(!is.na(y_vals))
    y_valid <- y_vals[valid_indices]
    # HP filter with lambda = 400
    hp_result <- hpfilter(y_valid, freq = 400, type = "lambda")
    # Assign back to full data
    full_cycle <- rep(NA, length(y_vals))
    full_cycle[valid_indices] <- hp_result$cycle
    data$ygap[data$iso == c] <- full_cycle
  }
}

cat("  output gap constructed\n\n")

# ============================================================================
# STEP 5: CREATE CONTROL VARIABLES
# ============================================================================
cat("Step 5: Creating control variables...\n")

data <- data %>%
  group_by(iso) %>%
  mutate(
    # Controls: 2 lags of outcome and treatment, l.ygap, ld.debtgdp
    L1Dy = dplyr::lag(Dy, 1),
    L2Dy = dplyr::lag(Dy, 2),
    L1dCAPB = dplyr::lag(dCAPB, 1),
    L2dCAPB = dplyr::lag(dCAPB, 2),
    L1ygap = dplyr::lag(ygap, 1),
    L1Ddebtgdp = dplyr::lag(debtgdp, 1) - dplyr::lag(debtgdp, 2)
  ) %>%
  ungroup()

cat("  Control variables created\n\n")

# ============================================================================
# STEP 6: RUN INDIVIDUAL HORIZON LP-IV REGRESSIONS
# ============================================================================
cat("Step 6: Running individual horizon LP-IV regressions...\n")

# Select controls
controls <- c("L1Dy", "L2Dy", "L1dCAPB", "L2dCAPB", "L1ygap", "L1Ddebtgdp")

# Initialize results storage
results <- data.frame(
  h = integer(),
  b = numeric(),
  d = numeric(),
  u = numeric(),
  Zero = numeric()
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
  reg_data <- data %>%
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
  results <- rbind(results, data.frame(
    h = h,
    b = b_coef,
    d = b_coef + 1.96 * se_coef,
    u = b_coef - 1.96 * se_coef,
    Zero = 0
  ))
}

cat("\nIndividual horizon results:\n")
print(results)
cat("\n")

# ============================================================================
# STEP 7: STACK DATA FOR JOINT REGRESSION
# ============================================================================
cat("Step 7: Stacking data for joint regression...\n")

# In Stata, this saves and appends temp files
# Here we'll create a list and bind

stack_list <- list()

for (h in 0:4) {
  outcome_var <- paste0("S", h, "y")
  treatment_var <- paste0("SdCAPB", h)

  stack_h <- data %>%
    select(all_of(outcome_var), all_of(treatment_var), size,
           all_of(controls), iso, year) %>%
    na.omit() %>%
    mutate(
      Y = !!sym(outcome_var),
      H = h,
      I = iso
    )

  # Create horizon-specific treatment, instruments, and controls
  # Initialize ALL to zero first
  for (j in 0:4) {
    stack_h[[paste0("T_h", j)]] <- 0
    stack_h[[paste0("Z_h", j)]] <- 0
    for (k in controls) {
      stack_h[[paste0(k, "_h", j)]] <- 0
    }
  }

  # Set current horizon to actual values
  stack_h[[paste0("T_h", h)]] <- stack_h[[treatment_var]]
  stack_h[[paste0("Z_h", h)]] <- stack_h$size
  for (k in controls) {
    stack_h[[paste0(k, "_h", h)]] <- stack_h[[paste0(k)]]
  }

  cat("  Horizon", h, ":", nrow(stack_h), "observations\n")
  stack_list[[h + 1]] <- stack_h
}

# Combine all horizons
stacked_data <- bind_rows(stack_list)

# Create fixed effect variable: FE = ctry*1000 + horizon
stacked_data$FE <- as.numeric(factor(stacked_data$I)) * 1000 + stacked_data$H

cat("  Total stacked observations:", nrow(stacked_data), "\n\n")

# ============================================================================
# STEP 8: RUN STACKED LP-IV REGRESSION
# ============================================================================
cat("Step 8: Running stacked LP-IV regression...\n")

# Diagnostic: check FE variable
cat("  Checking FE variable...\n")
cat("  Unique FE values:", length(unique(stacked_data$FE)), "\n")
cat("  Range of FE:", range(stacked_data$FE), "\n")
cat("  Number of factor levels:", nlevels(factor(stacked_data$FE)), "\n")

# Check for any issues with the data
cat("  Checking for missing values in key variables...\n")
cat("    Y:", sum(is.na(stacked_data$Y)), "missing\n")
cat("    T_h0:", sum(is.na(stacked_data$T_h0)), "missing\n")
cat("    Z_h0:", sum(is.na(stacked_data$Z_h0)), "missing\n")
cat("    FE:", sum(is.na(stacked_data$FE)), "missing\n")

# Build formula

# Create string for stacked control variables, e.g., X1_h0 + X1_h1 + ... + X6_h4
control_vars_stacked <- c()
for (k in controls) {
  for (j in 0:4) {
    control_vars_stacked <- c(control_vars_stacked, paste0(k, "_h", j))
  }
}
control_vars_stacked_str <- paste(control_vars_stacked, collapse = " + ")

# All T_h and X variables for each horizon
formula_stacked <- as.formula(paste0(
  "Y ~ T_h0 + T_h1 + T_h2 + T_h3 + T_h4 + ",
  control_vars_stacked_str, " + factor(FE) | ",
  "Z_h0 + Z_h1 + Z_h2 + Z_h3 + Z_h4 + ",
  control_vars_stacked_str, " + factor(FE)"
))

cat("  Estimating stacked model...\n")
stacked_model <- ivreg(formula_stacked, data = stacked_data)

cat("  Model rank:", stacked_model$rank, "\n")
cat("  Number of coefficients:", length(coef(stacked_model)), "\n")

# Cluster-robust standard errors by FE
vcov_stacked <- vcovCL(stacked_model, cluster = stacked_data$FE, type = "HC1")
coef_test_stacked <- coeftest(stacked_model, vcov = vcov_stacked)

# Check which T_h coefficients are in the model
cat("\n  Checking for T_h coefficients...\n")
all_coef_names <- rownames(coef_test_stacked)
T_h_in_model <- grep("^T_h", all_coef_names, value = TRUE)
cat("  T_h coefficients found:", paste(T_h_in_model, collapse = ", "), "\n\n")

# Extract T_h coefficients
T_h_coefs <- character(0)
for (h in 0:4) {
  name <- paste0("T_h", h)
  if (name %in% all_coef_names) {
    T_h_coefs <- c(T_h_coefs, name)
  } else {
    cat("  WARNING: T_h", h, "not found in stacked model!\n")
  }
}

if (length(T_h_coefs) > 0) {
  cat("\nStacked regression results for T_h coefficients:\n")
  print(coef_test_stacked[T_h_coefs, ])

  # Calculate average treatment effect
  avg_coef <- mean(coef_test_stacked[T_h_coefs, "Estimate"])

  # Calculate SE for average (delta method)
  vcov_T <- vcov_stacked[T_h_coefs, T_h_coefs, drop = FALSE]
  n_coefs <- length(T_h_coefs)
  weights <- rep(1/n_coefs, n_coefs)
  avg_se <- as.numeric(sqrt(t(weights) %*% vcov_T %*% weights))

  cat("\nAverage multiplier:", round(avg_coef, 4), "\n")
  cat("Average SE:", round(avg_se, 4), "\n\n")

  # Joint test
  joint_hypotheses <- paste0(T_h_coefs, " = 0")
  joint_test <- linearHypothesis(stacked_model, joint_hypotheses, vcov = vcov_stacked)

  cat("Joint test results:\n")
  print(joint_test)

  df_joint <- joint_test$Df[2]
  chi2_joint <- joint_test$Chisq[2]
  p_joint <- joint_test$`Pr(>Chisq)`[2]

} else {
  stop("No T_h coefficients found in stacked regression!")
}

# ============================================================================
# STEP 9: PREPARE DATA FOR PLOTTING
# ============================================================================
cat("\nStep 9: Preparing plot...\n")

# Extract coefficients from stacked model for plotting
results$B <- NA
results$U <- NA
results$D <- NA

for (h in 0:4) {
  name <- paste0("T_h", h)
  if (name %in% T_h_coefs) {
    results$B[results$h == h] <- coef_test_stacked[name, "Estimate"]
    results$U[results$h == h] <- coef_test_stacked[name, "Estimate"] +
      1.96 * coef_test_stacked[name, "Std. Error"]
    results$D[results$h == h] <- coef_test_stacked[name, "Estimate"] -
      1.96 * coef_test_stacked[name, "Std. Error"]
  } else {
    # Use individual LP estimates if stacked not available
    results$B[results$h == h] <- results$b[results$h == h]
    results$U[results$h == h] <- results$d[results$h == h]
    results$D[results$h == h] <- results$u[results$h == h]
  }
}

# Add average point (h=5)
results <- rbind(results, data.frame(
  h = 5,
  b = avg_coef,
  d = avg_coef + 1.96 * avg_se,
  u = avg_coef - 1.96 * avg_se,
  Zero = 0,
  B = avg_coef,
  U = avg_coef + 1.96 * avg_se,
  D = avg_coef - 1.96 * avg_se
))

# ============================================================================
# STEP 10: CREATE PLOT
# ============================================================================
cat("Step 10: Creating plot...\n")

p <- ggplot(results, aes(x = h)) +
  # Zero line
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  # Confidence band for h < 5
  geom_ribbon(data = filter(results, h < 5),
              aes(ymin = D, ymax = U),
              fill = "blue", alpha = 0.15) +
  # Point estimate line for h < 5
  geom_line(data = filter(results, h < 5),
            aes(y = B),
            color = "blue", linewidth = 1, alpha = 0.5) +
  # Average point with error bar
  geom_errorbar(data = filter(results, h == 5),
                aes(ymin = D, ymax = U),
                width = 0.1, color = "blue", alpha = 0.5, linewidth = 0.8) +
  geom_point(data = filter(results, h == 5),
             aes(y = B),
             color = "blue", size = 3) +
  geom_text(data = filter(results, h == 5),
            aes(y = B, label = sprintf("%.2f", B)),
            color = "blue", vjust = -1) +
  # Joint test annotation
  annotate("text", x = 2, y = 1.9,
           label = sprintf("Joint test, m(h)=0:\n\nChisq(%d)=%.1f (p=%.3f)",
                          df_joint, chi2_joint, p_joint),
           size = 6) +
  # Scales and labels
  scale_x_continuous(breaks = 0:5,
                    labels = c("0", "1", "2", "3", "4", "average"),
                    limits = c(0, 5.4)) +
  scale_y_continuous(breaks = seq(-5, 2, 1),
                    limits = c(-5, 2)) +
  labs(x = "Horizon, years, h",
       y = "Multiplier, m(h)",
       title = NULL) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11),
    plot.margin = margin(10, 10, 10, 10)
  )

print(p)

# Save plot
ggsave("jt25_figure2b.pdf", plot = p, width = 6, height = 6.9, units = "in")

# Save individual horizon results to CSV
results_baseline <- results %>%
  rename(
    m = b,
    lower_ci = u,
    upper_ci = d
  ) |> 
  filter(h < 5) |> 
  mutate(spec = "jt25_baseline") |> 
  select(spec, h, m, lower_ci, upper_ci)

# write_csv(results_baseline, "jt25_figure2b_results.csv")

cat("\n===============================================\n")
cat("Replication Complete!\n")
cat("===============================================\n")
cat("Plot saved as: jt25_figure2b.pdf\n")
cat("This replicates Figure 2b from the paper\n")