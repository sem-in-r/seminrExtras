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

# ============================================================================
# Shared Validation Helpers
# ============================================================================

#' Validate that input is a SEMinR model
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
#' @param model SEMinR model to check
#' @return TRUE if model has higher-order constructs, FALSE otherwise
#' @keywords internal
#' @noRd
has_higher_order <- function(model) {
  !is.null(model$hoc)
}

#' Validate model for prediction-based analysis (no higher-order constructs)
#'
#' @param model SEMinR model to validate
#' @param func_name Name of the calling function for error messages
#' @return TRUE if valid for prediction, FALSE otherwise (with warning)
#' @keywords internal
#' @noRd
validate_for_prediction <- function(model, func_name = "This function") {
 if (!validate_seminr_model(model, func_name)) {
    return(FALSE)
  }
  if (has_higher_order(model)) {
    warning("There is no published solution for applying PLSpredict to higher-order models.",
            call. = FALSE)
    return(FALSE)
  }
  return(TRUE)
}

# ============================================================================
# Endogenous Construct Helpers
# ============================================================================

#' Get endogenous latent variables from a model
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
#' @param model SEMinR model
#' @param constructs Character vector of construct names (defaults to all endogenous)
#' @return Character vector of item names
#' @keywords internal
#' @noRd
get_endogenous_items <- function(model, constructs = NULL) {
  if (is.null(constructs)) {
    constructs <- get_endogenous_constructs(model)
  }
  unlist(lapply(constructs, function(x) {
    seminr:::items_of_construct(construct = x, model = model)
  }))
}

# ============================================================================
# Loss Calculation Helpers
# ============================================================================

#' Calculate latent variable loss for a construct
#'
#' @param construct Name of the construct
#' @param model SEMinR model
#' @param error Error matrix (residuals)
#' @return Vector of loss values per observation
#' @keywords internal
#' @noRd
lv_loss <- function(construct, model, error) {
  items <- seminr:::items_of_construct(construct = construct, model = model)
  if (length(dim(error)) > 1) {
    loss <- rowMeans(error[, items, drop = FALSE]^2)
  } else {
    loss <- error[, items, drop = FALSE]^2
  }
  return(loss)
}

#' Calculate overall loss across constructs
#'
#' @param error Error matrix
#' @return Vector of overall loss per observation, or the error if 1D
#' @keywords internal
#' @noRd
overall_loss <- function(error) {
  if (length(dim(error)) > 1) {
    return(rowMeans(error))
  }
  return(error)
}

#' Calculate LV losses for multiple constructs
#'
#' @param constructs Character vector of construct names
#' @param model SEMinR model
#' @param error_matrix Error matrix (residuals)
#' @return Matrix of losses with constructs as columns
#' @keywords internal
#' @noRd
calculate_lv_losses <- function(constructs, model, error_matrix) {
  losses <- do.call("cbind", lapply(constructs, function(x) {
    lv_loss(construct = x, model = model, error = error_matrix)
  }))
  colnames(losses) <- constructs
  return(losses)
}

# ============================================================================
# Bootstrap CVPAT Helpers
# ============================================================================

#' Bootstrap CVPAT significance test
#'
#' Performs bootstrap-based significance testing for comparing two loss vectors.
#'
#' @param loss_m1 Loss vector for model 1
#' @param loss_m2 Loss vector for model 2
#' @param testtype Either "two.sided" or "greater"
#' @param nboot Number of bootstrap iterations
#' @return Data frame with t-values and p-values
#' @keywords internal
#' @noRd
bootstrap_cvpat <- function(loss_m1, loss_m2, testtype = "two.sided", nboot = 2000) {

  n <- length(loss_m1)
  org_t_test <- t.test(loss_m2, loss_m1,
                       alternative = testtype,
                       paired = TRUE)$statistic

  # Original average difference in losses
 org_d_bar <- mean(loss_m2 - loss_m1)

  # Differences in loss functions under the null
  d_null <- loss_m2 - loss_m1 - org_d_bar

  # Differences in loss functions
  d <- loss_m2 - loss_m1

  # Allocate memory for bootstrap
  boot_d_bar <- rep(0, nboot)
  m_losses <- cbind(loss_m1, loss_m2)
  t_stat <- rep(0, nboot)

  for (b in 1:nboot) {
    boot_sample <- m_losses[sample(1:length(d), length(d), replace = TRUE), ]
    t_stat[b] <- t.test(boot_sample[, 2], boot_sample[, 1],
                        mu = mean(d),
                        alternative = testtype,
                        paired = TRUE)$statistic
    boot_d_bar[b] <- mean(sample(d_null, length(d_null), replace = TRUE))
  }

  sorted_t_stat <- sort(t_stat, decreasing = FALSE)
  sorted_boot_d_bar <- sort(boot_d_bar, decreasing = FALSE)

  # Bootstrap variance on d_bar for t-test
  std <- sqrt(var(boot_d_bar))

  # Guard against division by zero
  if (std < .Machine$double.eps || is.na(std)) {
    t_stat_boot_var <- NA
    warning("Bootstrap variance near zero; t-statistic set to NA", call. = FALSE)
  } else {
    t_stat_boot_var <- org_d_bar / std
  }

  # Calculate p-values based on test type
  if (testtype == "two.sided") {
    p_value_perc_t <- (sum(sorted_t_stat > abs(org_t_test)) +
                       sum(sorted_t_stat <= -abs(org_t_test))) / nboot
    p_value_perc_d <- (sum(sorted_boot_d_bar > abs(org_d_bar)) +
                       sum(sorted_boot_d_bar <= -abs(org_d_bar))) / nboot
    p_value_var_t <- 2 * pt(-abs(t_stat_boot_var), (n - 1), lower.tail = TRUE)
  }

  if (testtype == "greater") {
    idx_t <- which(sorted_t_stat > org_t_test)
    idx_d <- which(sorted_boot_d_bar > org_d_bar)

    if (length(idx_t) == 0) {
      p_value_perc_t <- 0
      p_value_perc_d <- 0
    } else {
      p_value_perc_t <- 1 - (head(idx_t, 1) - 1) / (nboot + 1)
      p_value_perc_d <- 1 - (head(idx_d, 1) - 1) / (nboot + 1)
    }
    p_value_var_t <- pt(t_stat_boot_var, (n - 1), lower.tail = FALSE)
  }

  results <- data.frame(
    "Std. T value" = org_t_test,
    "Std. P value" = p_value_perc_t,
    "Boot T value" = as.numeric(t_stat_boot_var),
    "Boot P Value" = as.numeric(p_value_var_t),
    "Perc. P Value" = as.numeric(p_value_perc_d),
    check.names = FALSE
  )

  return(results)
}

#' Apply CVPAT bootstrap to each construct
#'
#' @param loss_one Loss matrix for model 1
#' @param loss_two Loss matrix for model 2
#' @param testtype Either "two.sided" or "greater"
#' @param nboot Number of bootstrap iterations
#' @return Data frame with per-construct CVPAT results
#' @keywords internal
#' @noRd
cvpat_per_construct <- function(loss_one, loss_two, testtype = "two.sided", nboot = 2000) {
  constructs <- colnames(loss_one)
  results <- as.data.frame(matrix(nrow = 0, ncol = 6))

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
