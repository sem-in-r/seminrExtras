# =============================================================================
# helpers.R - Internal helper functions for seminrExtras
# =============================================================================
# This file contains shared utility functions used across the package:
# - Model validation helpers
# - Endogenous construct extraction helpers
# - Loss calculation helpers for CVPAT
# - Bootstrap significance testing helpers
# =============================================================================

#' @importFrom stats pt
#' @importFrom stats t.test
#' @importFrom stats var
#' @importFrom stats sd
#' @importFrom utils head
#' @importFrom seminr predict_DA
#' @importFrom seminr predict_EA
#' @importFrom seminr mean_replacement
#' @importFrom seminr estimate_pls
#' @importFrom seminr predict_pls
#'
NULL

# =============================================================================
# VALIDATION HELPERS
# =============================================================================
# These functions validate that inputs are appropriate SEMinR model objects
# before proceeding with analysis. They provide consistent error handling
# across all exported functions.
# =============================================================================

#' Validate that input is a SEMinR model
#'
#' Checks if the provided object has the "seminr_model" class.
#' Used as the first validation step in all analysis functions.
#'
#' @param model The object to validate
#' @param func_name Name of the calling function for error messages
#' @return TRUE if valid, FALSE otherwise (with warning message)
#' @keywords internal
#' @noRd
validate_seminr_model <- function(model, func_name = "This function") {
  if (!any(class(model) == "seminr_model")) {
    warning(func_name, " only works with SEMinR models.", call. = FALSE)
    return(FALSE)
  }
  return(TRUE)
}

#' Check if model contains higher-order constructs
#'
#' Higher-order constructs (HOCs) require special handling in prediction.
#' Currently, there is no published solution for PLSpredict with HOCs.
#'
#' @param model SEMinR model to check
#' @return TRUE if model has higher-order constructs, FALSE otherwise
#' @keywords internal
#' @noRd
has_higher_order <- function(model) {
  !is.null(model$hoc)
}

#' Validate model for prediction-based analysis
#'
#' Combines seminr_model validation with higher-order construct check.
#' CVPAT functions require models without higher-order constructs because
#' PLSpredict cannot currently handle them.
#'
#' @param model SEMinR model to validate
#' @param func_name Name of the calling function for error messages
#' @return TRUE if valid for prediction, FALSE otherwise (with warning)
#' @keywords internal
#' @noRd
validate_for_prediction <- function(model, func_name = "This function") {
  # First check: is it a seminr_model at all?
  if (!validate_seminr_model(model, func_name)) {
    return(FALSE)
  }
  # Second check: does it have higher-order constructs?
  if (has_higher_order(model)) {
    warning("There is no published solution for applying PLSpredict to higher-order models.",
            call. = FALSE)
    return(FALSE)
  }
  return(TRUE)
}

# =============================================================================
# ENDOGENOUS CONSTRUCT HELPERS
# =============================================================================
# These functions extract information about endogenous (dependent) constructs
# from SEMinR models. Endogenous constructs are those that have incoming paths
# in the structural model - they are the variables we want to predict.
# =============================================================================

#' Get endogenous latent variables from a model
#'
#' Endogenous constructs are those with incoming paths in the structural model.
#' These are the constructs whose prediction accuracy we evaluate in CVPAT.
#'
#' @param model SEMinR model
#' @return Character vector of endogenous construct names
#' @keywords internal
#' @noRd
get_endogenous_constructs <- function(model) {
  seminr:::all_endogenous(model$smMatrix)
}

#' Get measurement items for endogenous constructs
#'
#' Retrieves the indicator/item names that measure the specified constructs.
#' These items are needed to calculate prediction errors at the item level.
#'
#' @param model SEMinR model
#' @param constructs Character vector of construct names (defaults to all endogenous)
#' @return Character vector of item names
#' @keywords internal
#' @noRd
get_endogenous_items <- function(model, constructs = NULL) {
  # Default to all endogenous constructs if none specified
 if (is.null(constructs)) {
    constructs <- get_endogenous_constructs(model)
  }
  # Get items for each construct and flatten into a single vector
  unlist(lapply(constructs, function(x) {
    seminr:::items_of_construct(construct = x, model = model)
  }))
}

# =============================================================================
# LOSS CALCULATION HELPERS
# =============================================================================
# These functions calculate prediction loss (squared error) at the construct
# and overall model level. Loss is the key metric in CVPAT - lower loss
# indicates better predictive performance.
#
# Loss formula: L = mean(error^2) where error = actual - predicted
# =============================================================================

#' Calculate latent variable loss for a construct
#'
#' Computes the mean squared error for a single construct by averaging
#' the squared prediction errors across all items measuring that construct.
#'
#' @param construct Name of the construct
#' @param model SEMinR model (needed to identify which items belong to construct)
#' @param error Error matrix (residuals) with items as columns
#' @return Vector of loss values, one per observation
#' @keywords internal
#' @noRd
lv_loss <- function(construct, model, error) {
  # Get the items that measure this construct
  items <- seminr:::items_of_construct(construct = construct, model = model)

  # Calculate mean squared error across items for each observation
  # Handle both matrix (multiple items) and vector (single item) cases
  if (length(dim(error)) > 1) {
    # Multiple items: average squared errors across items for each row
    loss <- rowMeans(error[, items, drop = FALSE]^2)
  } else {
    # Single item: just square the errors
    loss <- error[, items, drop = FALSE]^2
  }
  return(loss)
}

#' Calculate overall loss across constructs
#'
#' Computes the mean loss across all constructs for each observation.
#' This gives a single overall measure of prediction accuracy.
#'
#' @param error Error matrix with constructs as columns
#' @return Vector of overall loss per observation, or the error if 1D
#' @keywords internal
#' @noRd
overall_loss <- function(error) {
  # Average across columns (constructs) if matrix, otherwise return as-is
  if (length(dim(error)) > 1) {
    return(rowMeans(error))
  }
  return(error)
}

#' Calculate LV losses for multiple constructs
#'
#' Wrapper function that calculates loss for each construct and returns
#' a matrix with one column per construct. This is the main entry point
#' for loss calculation in the CVPAT functions.
#'
#' @param constructs Character vector of construct names
#' @param model SEMinR model
#' @param error_matrix Error matrix (residuals) with items as columns
#' @return Matrix of losses with constructs as columns, observations as rows
#' @keywords internal
#' @noRd
calculate_lv_losses <- function(constructs, model, error_matrix) {
  # Calculate loss for each construct and bind into matrix
  losses <- do.call("cbind", lapply(constructs, function(x) {
    lv_loss(construct = x, model = model, error = error_matrix)
  }))
  colnames(losses) <- constructs
  return(losses)
}

# =============================================================================
# BOOTSTRAP CVPAT HELPERS
# =============================================================================
# These functions implement the bootstrap-based significance testing for CVPAT.
# The bootstrap approach is used because the distribution of loss differences
# is unknown. By resampling, we can estimate the sampling distribution and
# compute p-values without parametric assumptions.
#
# Reference: Liengaard et al. (2021), Decision Sciences, 52(2), 362-392.
# =============================================================================

#' Bootstrap CVPAT significance test
#'
#' Performs bootstrap-based significance testing to compare two loss vectors.
#' Tests whether the difference in predictive loss between two models (or
#' between a model and a benchmark) is statistically significant.
#'
#' The function computes three types of statistics:
#' 1. Standard paired t-test (parametric baseline)
#' 2. Bootstrap t-statistic using resampled variance
#' 3. Percentile-based p-values from bootstrap distribution
#'
#' @param loss_m1 Loss vector for model 1 (or PLS model)
#' @param loss_m2 Loss vector for model 2 (or benchmark)
#' @param testtype "two.sided" tests if losses differ; "greater" tests if m2 > m1
#' @param nboot Number of bootstrap iterations (more = more precise but slower)
#' @return Data frame with t-values and p-values from different methods
#' @keywords internal
#' @noRd
bootstrap_cvpat <- function(loss_m1, loss_m2, testtype = "two.sided", nboot = 2000) {

  n <- length(loss_m1)

  # Step 1: Compute the standard paired t-test statistic
  # This serves as a parametric baseline for comparison
  org_t_test <- t.test(loss_m2, loss_m1,
                       alternative = testtype,
                       paired = TRUE)$statistic

  # Step 2: Compute the original mean difference in losses
  # Negative value means model 1 has lower loss (better prediction)
  org_d_bar <- mean(loss_m2 - loss_m1)

  # Step 3: Center the differences for null hypothesis bootstrap
  # Under H0, the true difference is 0, so we subtract the observed mean
  d_null <- loss_m2 - loss_m1 - org_d_bar

  # Original (uncentered) differences for bootstrap sampling
  d <- loss_m2 - loss_m1

  # Step 4: Allocate storage for bootstrap results
  boot_d_bar <- rep(0, nboot)      # Bootstrap mean differences
  m_losses <- cbind(loss_m1, loss_m2)  # Combined matrix for resampling
  t_stat <- rep(0, nboot)          # Bootstrap t-statistics

  # Step 5: Bootstrap loop - resample and compute statistics
  for (b in 1:nboot) {
    # Resample observations with replacement (paired resampling)
    boot_sample <- m_losses[sample(1:length(d), length(d), replace = TRUE), ]

    # Compute t-statistic for this bootstrap sample
    # mu = mean(d) tests if bootstrap differs from original
    t_stat[b] <- t.test(boot_sample[, 2], boot_sample[, 1],
                        mu = mean(d),
                        alternative = testtype,
                        paired = TRUE)$statistic

    # Compute mean of resampled null-centered differences
    boot_d_bar[b] <- mean(sample(d_null, length(d_null), replace = TRUE))
  }

  # Step 6: Sort bootstrap distributions for percentile calculations
  sorted_t_stat <- sort(t_stat, decreasing = FALSE)
  sorted_boot_d_bar <- sort(boot_d_bar, decreasing = FALSE)

  # Step 7: Compute bootstrap standard error for t-statistic
  std <- sqrt(var(boot_d_bar))

  # Guard against division by zero (can happen with identical predictions)
  if (std < .Machine$double.eps || is.na(std)) {
    t_stat_boot_var <- NA
    warning("Bootstrap variance near zero; t-statistic set to NA", call. = FALSE)
  } else {
    # Bootstrap t-statistic using bootstrap variance estimate
    t_stat_boot_var <- org_d_bar / std
  }

  # Step 8: Calculate p-values based on test type
  if (testtype == "two.sided") {
    # Two-sided: count extreme values in both tails
    p_value_perc_t <- (sum(sorted_t_stat > abs(org_t_test)) +
                       sum(sorted_t_stat <= -abs(org_t_test))) / nboot
    p_value_perc_d <- (sum(sorted_boot_d_bar > abs(org_d_bar)) +
                       sum(sorted_boot_d_bar <= -abs(org_d_bar))) / nboot
    # Parametric p-value from t-distribution
    p_value_var_t <- 2 * pt(-abs(t_stat_boot_var), (n - 1), lower.tail = TRUE)
  }

  if (testtype == "greater") {
    # One-sided: test if model 2 has greater loss than model 1
    idx_t <- which(sorted_t_stat > org_t_test)
    idx_d <- which(sorted_boot_d_bar > org_d_bar)

    # Handle edge case where no bootstrap values exceed observed
    if (length(idx_t) == 0) {
      p_value_perc_t <- 0
      p_value_perc_d <- 0
    } else {
      # Percentile-based p-value: proportion of bootstrap values more extreme
      p_value_perc_t <- 1 - (head(idx_t, 1) - 1) / (nboot + 1)
      p_value_perc_d <- 1 - (head(idx_d, 1) - 1) / (nboot + 1)
    }
    # Parametric p-value for one-sided test
    p_value_var_t <- pt(t_stat_boot_var, (n - 1), lower.tail = FALSE)
  }

  # Step 9: Compile results into data frame
  results <- data.frame(
    "Std. T value" = org_t_test,        # Standard t-test statistic
    "Std. P value" = p_value_perc_t,    # Percentile p-value (t-based)
    "Boot T value" = as.numeric(t_stat_boot_var),  # Bootstrap t-statistic
    "Boot P Value" = as.numeric(p_value_var_t),    # Bootstrap p-value (recommended)
    "Perc. P Value" = as.numeric(p_value_perc_d),  # Percentile p-value (d-based)
    check.names = FALSE
  )

  return(results)
}

#' Apply CVPAT bootstrap to each construct
#'
#' Runs bootstrap_cvpat separately for each endogenous construct,
#' allowing users to see which specific constructs drive overall differences.
#'
#' @param loss_one Loss matrix for model 1, columns = constructs
#' @param loss_two Loss matrix for model 2, columns = constructs
#' @param testtype Either "two.sided" or "greater"
#' @param nboot Number of bootstrap iterations
#' @return Data frame with per-construct CVPAT results
#' @keywords internal
#' @noRd
cvpat_per_construct <- function(loss_one, loss_two, testtype = "two.sided", nboot = 2000) {
  constructs <- colnames(loss_one)
  results <- as.data.frame(matrix(nrow = 0, ncol = 6))

  # Loop through each construct and run bootstrap test
  for (construct in constructs) {
    boot_result <- bootstrap_cvpat(loss_one[, construct],
                                   loss_two[, construct],
                                   testtype = testtype,
                                   nboot = nboot)
    results <- rbind(results, c(construct, unlist(boot_result)))
  }

  colnames(results) <- c("Construct", "Std. T value", "Std. P value",
                         "Boot T value", "Boot P Value", "Perc. P Value")
  return(results)
}
