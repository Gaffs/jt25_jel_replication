#' Local Projection IV Panel Estimation for Fiscal Multipliers
#'
#' This function performs LP-IV panel estimation with various transformations for
#' consistent fiscal multiplier calculations, inspired by Jorda and Taylor (2025)
#'
#' @param df Data frame containing panel data
#' @param outcome Character string. Name of outcome variable (e.g., "rgdp")
#' @param treatment Character string. Name of treatment variable (e.g., "rfb")
#' @param transformation Character string. Type of transformation: "HBR", "GK", or "CP"
#' @param max_h Integer. Maximum horizon for LP estimation
#' @param gk_trend Character string. Trend type for GK transformation: "hp" or "poly"
#' @param outcome_nlags Integer. Number of lags of transformed outcome to include as controls
#' @param treatment_nlags Integer. Number of lags of transformed treatment to include as controls
#' @param instrument Character string. Name of instrumental variable
#' @param group_var Character string. Name of grouping variable for panel (e.g., "iso")
#' @param time_var Character string. Name of time variable (e.g., "year")
#' @param additional_controls Character vector. Names of additional control variables
#' @param hp_lambda Numeric. Lambda parameter for HP filter (default: 400)
#' @param poly_degree Integer. Degree of polynomial for polynomial trend (default: 3)
#'
#' @return A list containing:
#'   - results: Data frame with coefficients, standard errors, and confidence intervals by horizon
#'   - models: List of ivreg model objects for each horizon
#'   - data_transformed: Data frame with all transformations applied
#'
#' @examples
#' \dontrun{
#' results <- lp_iv_panel_fm(
#'   df = data_clean,
#'   outcome = "rgdp",
#'   treatment = "rfb",
#'   transformation = "HBR",
#'   max_h = 4,
#'   gk_trend = "hp",
#'   outcome_nlags = 2,
#'   treatment_nlags = 2,
#'   instrument = "size",
#'   group_var = "iso",
#'   time_var = "year",
#'   additional_controls = c("L1ygap", "L1Ddebtgdp")
#' )
#' }
lp_iv_panel_fm <- function(df,
                                   outcome,
                                   treatment,
                                   transformation,
                                   max_h,
                                   gk_trend = "hp",
                                   outcome_nlags = 2,
                                   treatment_nlags = 2,
                                   instrument,
                                   group_var,
                                   time_var,
                                   additional_controls = NULL,
                                   hp_lambda = 400,
                                   poly_degree = 3) {

  # Load required packages
  require(dplyr)
  require(AER)
  require(sandwich)
  require(lmtest)
  require(mFilter)

  # Validate inputs
  if (!transformation %in% c("HBR", "GK", "CP")) {
    stop("transformation must be one of 'HBR', 'GK', or 'CP'")
  }

  if (!gk_trend %in% c("hp", "poly")) {
    stop("gk_trend must be one of 'hp' or 'poly'")
  }

  # Create a copy of the data
  data_transformed <- df

  # ============================================================================
  # Step 1: Create trend for GK transformation (if needed)
  # ============================================================================

  if (transformation == "GK") {
    cat("Creating trend for GK transformation using", gk_trend, "method...\n")

    data_transformed[[paste0(outcome, "trend")]] <- NA

    groups <- unique(data_transformed[[group_var]])

    if (gk_trend == "hp") {
      for (g in groups) {
        group_data <- data_transformed %>%
          filter(!!sym(group_var) == g) %>%
          arrange(!!sym(time_var))

        y_vals <- log(group_data[[outcome]])

        if (sum(!is.na(y_vals)) > 0) {
          # Adjust for NAs
          valid_indices <- which(!is.na(y_vals))
          y_valid <- y_vals[valid_indices]

          # HP filter
          hp_result <- hpfilter(y_valid, freq = hp_lambda, type = "lambda")

          # Assign back to full data
          full_trend <- rep(NA, length(y_vals))
          full_trend[valid_indices] <- hp_result$trend

          data_transformed[[paste0(outcome, "trend")]][data_transformed[[group_var]] == g] <- exp(full_trend)
        }
      }
    } else if (gk_trend == "poly") {
      for (g in groups) {
        group_data <- data_transformed %>%
          filter(!!sym(group_var) == g) %>%
          arrange(!!sym(time_var))

        y_vals <- log(group_data[[outcome]])

        if (sum(!is.na(y_vals)) > 0) {
          x_vals <- group_data[[time_var]][!is.na(y_vals)]
          poly_fit <- lm(y_vals[!is.na(y_vals)] ~ poly(x_vals, poly_degree, raw = TRUE))
          full_trend <- predict(poly_fit)

          data_transformed[[paste0(outcome, "trend")]][
            data_transformed[[group_var]] == g &
              data_transformed[[time_var]] %in% x_vals] <- exp(full_trend)
        }
      }
    }
  }

  # ============================================================================
  # Step 2: Create transformations
  # ============================================================================

  cat("Creating", transformation, "transformations...\n")

  vars_to_transform <- c(outcome, treatment)

  # Create transformations based on type
  if (transformation == "HBR") {
    # Hall-Barro-Redlick (HBR) transformation
    for (var in vars_to_transform) {
      for (h in 0:max_h) {
        data_transformed <- data_transformed %>%
          group_by(!!sym(group_var)) %>%
          mutate(
            !!paste0(transformation, h, var) :=
              (dplyr::lead(!!sym(var), h) - dplyr::lag(!!sym(var), 1)) /
              dplyr::lag(!!sym(outcome), 1)
          ) %>%
          ungroup()
      }
    }
  } else if (transformation == "GK") {
    # Gordon-Krenn (GK) transformation
    for (var in vars_to_transform) {
      for (h in 0:max_h) {
        data_transformed <- data_transformed %>%
          group_by(!!sym(group_var)) %>%
          mutate(
            !!paste0(transformation, h, var) :=
              (dplyr::lead(!!sym(var), h)) / dplyr::lead(!!sym(paste0(outcome, "trend")), h) -
              (dplyr::lag(!!sym(var), 1)) / dplyr::lag(!!sym(paste0(outcome, "trend")), 1)
          ) %>%
          ungroup()
      }
    }
  } else if (transformation == "CP") {
    # Canova-Pappa (CP) transformation
    for (var in vars_to_transform) {
      for (h in 0:max_h) {
        data_transformed <- data_transformed %>%
          group_by(!!sym(group_var)) %>%
          mutate(
            !!paste0(transformation, h, var) :=
              (dplyr::lead(!!sym(var), h) - (h + 1) * dplyr::lag(!!sym(var), 1)) /
              dplyr::lag(!!sym(outcome), 1)
          ) %>%
          ungroup()
      }
    }
  } else {
    stop("Invalid transformation type. Must be 'HBR', 'GK', or 'CP'")
  }

  # ============================================================================
  # Step 3: Create cumulative sums across horizons
  # ============================================================================

  cat("Creating cumulative sums...\n")

  for (var in vars_to_transform) {
    for (h in 0:max_h) {
      cols_to_sum <- paste0(transformation, 0:h, var)
      data_transformed <- data_transformed %>%
        mutate(!!paste0("S", transformation, h, var) :=
                 rowSums(select(., all_of(cols_to_sum)), na.rm = FALSE))
    }
  }

  # ============================================================================
  # Step 4: Create lags of transformed variables
  # ============================================================================

  cat("Creating lags of transformed variables...\n")

  # Outcome lags
  outcome_var_base <- paste0("S", transformation, "0", outcome)
  for (lg in 1:outcome_nlags) {
    lagged_var <- paste0("L", lg, outcome_var_base)
    data_transformed <- data_transformed %>%
      group_by(!!sym(group_var)) %>%
      mutate(!!lagged_var := dplyr::lag(!!sym(outcome_var_base), lg)) %>%
      ungroup()
  }

  # Treatment lags
  treatment_var_base <- paste0("S", transformation, "0", treatment)
  for (lg in 1:treatment_nlags) {
    lagged_var <- paste0("L", lg, treatment_var_base)
    data_transformed <- data_transformed %>%
      group_by(!!sym(group_var)) %>%
      mutate(!!lagged_var := dplyr::lag(!!sym(treatment_var_base), lg)) %>%
      ungroup()
  }

  # ============================================================================
  # Step 5: Prepare control variables
  # ============================================================================

  control_vars <- c(
    paste0("L", 1:outcome_nlags, outcome_var_base),
    paste0("L", 1:treatment_nlags, treatment_var_base)
  )

  if (!is.null(additional_controls)) {
    control_vars <- c(control_vars, additional_controls)
  }

  # ============================================================================
  # Step 6: Run LP-IV regressions for each horizon
  # ============================================================================

  cat("Running LP-IV regressions...\n")

  results <- data.frame(
    h = integer(),
    m = numeric(),
    se = numeric(),
    upper_ci = numeric(),
    lower_ci = numeric(),
    n_obs = integer()
  )

  lp_models <- list()

  for (h in 0:max_h) {
    cat("  Horizon", h, "...\n")

    # Define outcome and treatment for this horizon
    outcome_var <- paste0("S", transformation, h, outcome)
    treatment_var <- paste0("S", transformation, h, treatment)

    # Select complete cases
    reg_data <- data_transformed %>%
      select(all_of(outcome_var), all_of(treatment_var),
             all_of(instrument), all_of(control_vars),
             all_of(group_var)) %>%
      na.omit()

    n_obs <- nrow(reg_data)
    cat("    Observations:", n_obs, "\n")

    if (n_obs == 0) {
      warning(paste("No observations for horizon", h))
      next
    }

    # Demean for fixed effects (within transformation)
    reg_data_demeaned <- reg_data %>%
      group_by(!!sym(group_var)) %>%
      mutate(across(c(all_of(outcome_var), all_of(treatment_var),
                      all_of(instrument), all_of(control_vars)),
                    ~ . - mean(., na.rm = TRUE))) %>%
      ungroup()

    # Create formula for IV regression
    control_vars_str <- paste(control_vars, collapse = " + ")

    formula_iv <- as.formula(paste0(
      outcome_var, " ~ ", treatment_var, " + ", control_vars_str, " | ",
      instrument, " + ", control_vars_str
    ))

    # Run IV regression
    model <- ivreg(formula_iv, data = reg_data_demeaned)

    # Cluster-robust standard errors
    vcov_cluster <- vcovCL(model, cluster = reg_data[[group_var]], type = "HC1")
    coef_test <- coeftest(model, vcov = vcov_cluster)

    # Store model
    lp_models[[h + 1]] <- list(
      model = model,
      coef_test = coef_test,
      n_obs = n_obs
    )

    # Extract coefficient for treatment variable
    coef_names <- rownames(coef_test)
    treatment_coef_name <- treatment_var

    if (!(treatment_var %in% coef_names)) {
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
      m = b_coef,
      se = se_coef,
      upper_ci = b_coef + 1.96 * se_coef,
      lower_ci = b_coef - 1.96 * se_coef,
      n_obs = n_obs
    ))
  }

  # ============================================================================
  # Return results
  # ============================================================================

  return(list(
    results = results,
    models = lp_models,
    data_transformed = data_transformed,
    transformation = transformation,
    outcome = outcome,
    treatment = treatment,
    max_h = max_h
  ))
}
