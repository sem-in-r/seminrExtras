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
#' @importFrom stats median
#' @importFrom stats lm
#' @importFrom stats coef
#' @importFrom stats dnorm
#' @importFrom stats cov quantile p.adjust
#' @importFrom utils head combn
#' @importFrom seminr predict_DA
#' @importFrom seminr predict_EA
#' @importFrom seminr mean_replacement
#' @importFrom seminr estimate_pls
#' @importFrom seminr predict_pls
#' @importFrom rpart rpart
#' @importFrom graphics abline axis barplot legend lines matlines matplot
#'   par plot points polygon rect text
#' @importFrom grDevices adjustcolor palette palette.colors
#'
NULL

# =============================================================================
# VALIDATION HELPERS
# =============================================================================

#' @noRd
validate_seminr_model <- function(model, func_name = "This function") {
  if (!inherits(model, "seminr_model")) {
    warning(func_name, " only works with SEMinR models.", call. = FALSE)
    return(FALSE)
  }
  return(TRUE)
}

#' @noRd
has_higher_order <- function(model) {
  !is.null(model$hoc)
}

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

# =============================================================================
# SEMINR INTERNAL REIMPLEMENTATIONS
# =============================================================================
# These replace seminr:::  calls to avoid CRAN NOTEs about unexported imports.

#' Extract unique endogenous (target) constructs from a structural model matrix.
#' @noRd
all_endogenous <- function(smMatrix) {
  unique(smMatrix[, "target"])
}

#' Get the measurement indicator names for a construct.
#' @noRd
items_of_construct <- function(construct, model) {
  model$mmMatrix[model$mmMatrix[, 1] == construct, 2]
}

#' Compute bootstrap percentile confidence interval for a path.
#' @noRd
conf_int <- function(boot_array, from, to, through = NULL, alpha = 0.05) {
  if (is.null(through)) {
    coefficient <- boot_array[from, to, ]
  } else {
    coefficient <- boot_array[from, through, ] * boot_array[through, to, ]
  }
  quantile(coefficient, probs = c(alpha / 2, 1 - alpha / 2))
}

#' Assign the table_output class to a matrix for seminr-compatible printing.
#' @noRd
convert_to_table_output <- function(matrix) {
  class(matrix) <- c("table_output", class(matrix))
  matrix
}

# =============================================================================
# ENDOGENOUS CONSTRUCT HELPERS
# =============================================================================

#' @noRd
get_endogenous_constructs <- function(model) {
  all_endogenous(model$smMatrix)
}

#' @noRd
get_endogenous_items <- function(model, constructs = NULL) {
  if (is.null(constructs)) {
    constructs <- get_endogenous_constructs(model)
  }
  # Interaction constructs have no measurement indicators — exclude them
  constructs <- constructs[!grepl("\\*", constructs)]
  unlist(lapply(constructs, function(x) {
    items_of_construct(construct = x, model = model)
  }))
}

# =============================================================================
# LOSS CALCULATION HELPERS
# =============================================================================

#' @noRd
lv_loss <- function(construct, model, error) {
  items <- items_of_construct(construct = construct, model = model)
  if (length(dim(error)) > 1) {
    loss <- rowMeans(error[, items, drop = FALSE]^2)
  } else {
    loss <- as.vector(error[, items, drop = FALSE]^2)
  }
  return(loss)
}

#' @noRd
overall_loss <- function(error) {
  if (length(dim(error)) > 1) {
    return(rowMeans(error))
  }
  return(error)
}

#' @noRd
calculate_lv_losses <- function(constructs, model, error_matrix) {
  losses <- do.call("cbind", lapply(constructs, function(x) {
    lv_loss(construct = x, model = model, error = error_matrix)
  }))
  colnames(losses) <- constructs
  return(losses)
}

# =============================================================================
# BOOTSTRAP CVPAT HELPERS
# =============================================================================

#' @noRd
bootstrap_cvpat <- function(loss_m1, loss_m2, testtype = "two.sided", nboot = 2000) {

  n <- length(loss_m1)

  org_t_test <- t.test(loss_m2, loss_m1,
                       alternative = testtype,
                       paired = TRUE)$statistic

  org_d_bar <- mean(loss_m2 - loss_m1)

  d_null <- loss_m2 - loss_m1 - org_d_bar
  d <- loss_m2 - loss_m1

  boot_d_bar <- rep(0, nboot)
  m_losses <- cbind(loss_m1, loss_m2)
  t_stat <- rep(0, nboot)

  for (b in seq_len(nboot)) {
    boot_sample <- m_losses[sample(seq_along(d), length(d), replace = TRUE), ]

    t_stat[b] <- t.test(boot_sample[, 2], boot_sample[, 1],
                        mu = mean(d),
                        alternative = testtype,
                        paired = TRUE)$statistic

    boot_d_bar[b] <- mean(sample(d_null, length(d_null), replace = TRUE))
  }

  sorted_t_stat <- sort(t_stat, decreasing = FALSE)
  sorted_boot_d_bar <- sort(boot_d_bar, decreasing = FALSE)

  std <- sqrt(var(boot_d_bar))

  if (std < .Machine$double.eps || is.na(std)) {
    t_stat_boot_var <- NA
    warning("Bootstrap variance near zero; t-statistic set to NA", call. = FALSE)
  } else {
    t_stat_boot_var <- org_d_bar / std
  }

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
    } else {
      p_value_perc_t <- 1 - (head(idx_t, 1) - 1) / (nboot + 1)
    }

    if (length(idx_d) == 0) {
      p_value_perc_d <- 0
    } else {
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
