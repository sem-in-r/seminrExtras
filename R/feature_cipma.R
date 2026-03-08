# =============================================================================
# feature_cipma.R - Combined Importance-Performance Map Analysis (cIPMA)
# =============================================================================
# Implements IPMA and cIPMA (IPMA + NCA) for PLS-SEM models.
#
# References:
# - Ringle, C. M. & Sarstedt, M. (2016). Gain More Insight from Your PLS-SEM
#   Results: The Importance-Performance Map Analysis. Industrial Management &
#   Data Systems, 119(9), 1865-1886.
# - Sarstedt, M., Richter, N. F., Hauff, S. & Ringle, C. M. (2024). Combined
#   Importance-Performance Map Analysis (cIPMA): A SmartPLS 4 Tutorial.
#   Journal of Marketing Analytics, 12, 746-760.
# - Hauff, S. et al. (2024). Importance and Performance in PLS-SEM and NCA:
#   Introducing the cIPMA. Journal of Retailing and Consumer Services, 78,
#   103723.
# =============================================================================

# =============================================================================
# INTERNAL HELPERS
# =============================================================================

#' Detect interaction constructs (names containing "*").
#' @noRd
is_interaction_construct <- function(name) {
  grepl("*", name, fixed = TRUE)
}

#' Check that all outer weights are positive for IPMA constructs.
#' For HOC constructs, also checks LOC-to-indicator weights.
#' Returns character vector of constructs with negative weights.
#' @noRd
check_positive_weights <- function(model, constructs) {
  neg <- character(0)
  construct_cols <- colnames(model$outer_weights)

  for (construct in constructs) {
    items <- seminr:::items_of_construct(construct, model)
    weights <- model$outer_weights[items, construct]
    if (any(weights < 0)) {
      neg <- c(neg, construct)
      next
    }

    # For HOC: also check LOC-to-indicator weights
    loc_names <- items[items %in% construct_cols]
    if (length(loc_names) > 0) {
      loc_neg <- check_positive_weights(model, loc_names)
      if (length(loc_neg) > 0) {
        neg <- c(neg, construct)
      }
    }
  }
  neg
}

#' Compute construct performance scores (0-100 rescaled).
#'
#' Performance is the weighted average of individual indicator performances,
#' using normalized outer weights. Indicator performance:
#' (mean_observed - scale_min) / (scale_max - scale_min) * 100.
#'
#' For higher-order constructs (HOC), performance is computed by chaining
#' through the lower-order constructs (LOCs): first compute each LOC's
#' performance from its actual indicators, then aggregate using HOC-to-LOC
#' weights. This applies to both two-stage and repeated-indicators HOCs.
#'
#' @return Named numeric vector of construct performances (0-100).
#' @noRd
compute_ipma_performance <- function(model, constructs, scale_min, scale_max) {
  scale_range <- scale_max - scale_min
  performance <- numeric(length(constructs))
  names(performance) <- constructs

  construct_cols <- colnames(model$outer_weights)

  for (construct in constructs) {
    items <- seminr:::items_of_construct(construct, model)

    # Detect HOC: items that are themselves construct names in outer_weights
    loc_names <- items[items %in% construct_cols]

    if (length(loc_names) > 0) {
      # HOC: recursively compute LOC performances, then aggregate
      loc_perf <- compute_ipma_performance(model, loc_names, scale_min, scale_max)
      hoc_weights <- model$outer_weights[loc_names, construct]
      performance[construct] <- sum(hoc_weights * loc_perf) / sum(hoc_weights)
    } else {
      # Regular construct: compute from indicators
      available <- items %in% colnames(model$data)
      if (!all(available)) {
        warning("Construct '", construct, "': some items not found in data. ",
                "Performance set to NA.", call. = FALSE)
        performance[construct] <- NA_real_
        next
      }

      weights <- model$outer_weights[items, construct]
      indicator_means <- colMeans(model$data[, items, drop = FALSE])
      indicator_perf <- (indicator_means - scale_min) / scale_range * 100

      # Weighted average using normalized positive weights
      performance[construct] <- sum(weights * indicator_perf) / sum(weights)
    }
  }

  performance
}

#' Compute per-observation performance scores for each construct.
#'
#' Used to obtain SD of rescaled construct scores for unstandardized
#' total effect computation. Returns an N x length(constructs) matrix.
#'
#' For HOC constructs, chains through LOC observation performances
#' and aggregates using HOC-to-LOC weights.
#' @noRd
compute_observation_performance <- function(model, constructs, scale_min, scale_max) {
  scale_range <- scale_max - scale_min
  n <- nrow(model$data)
  perf_scores <- matrix(NA_real_, nrow = n, ncol = length(constructs),
                         dimnames = list(NULL, constructs))

  construct_cols <- colnames(model$outer_weights)

  for (construct in constructs) {
    items <- seminr:::items_of_construct(construct, model)

    # Detect HOC: items that are themselves construct names
    loc_names <- items[items %in% construct_cols]

    if (length(loc_names) > 0) {
      # HOC: recursively compute LOC observation performances, then aggregate
      loc_obs_perf <- compute_observation_performance(model, loc_names, scale_min, scale_max)
      hoc_weights <- model$outer_weights[loc_names, construct]
      perf_scores[, construct] <- loc_obs_perf %*% hoc_weights / sum(hoc_weights)
    } else {
      # Regular construct: compute from indicators
      available <- items %in% colnames(model$data)
      if (!all(available)) next

      weights <- model$outer_weights[items, construct]
      indicator_perf <- (model$data[, items, drop = FALSE] - scale_min) / scale_range * 100
      perf_scores[, construct] <- as.matrix(indicator_perf) %*% weights / sum(weights)
    }
  }

  perf_scores
}

#' Compute total effects matrix from a path coefficient matrix.
#'
#' Total effects T = (I - B)^{-1} - I, where B[from, to] is the
#' direct path coefficient from construct 'from' to construct 'to'.
#' @noRd
compute_total_effects <- function(path_coef_matrix) {
  k <- nrow(path_coef_matrix)
  I_mat <- diag(k)
  dimnames(I_mat) <- dimnames(path_coef_matrix)
  solve(I_mat - path_coef_matrix) - I_mat
}

#' Compute unstandardized total effects for IPMA.
#'
#' Converts standardized path coefficients to unstandardized form using
#' the SDs of rescaled (0-100) construct scores, then computes total effects.
#'
#' Unstandardized path coefficient:
#'   beta_unstd(X -> Y) = beta_std(X -> Y) * SD(Y_rescaled) / SD(X_rescaled)
#'
#' @return Total effects matrix with constructs as row/column names.
#' @noRd
compute_unstd_total_effects <- function(model, constructs, scale_min, scale_max) {
  # Compute SD of rescaled performance scores
  obs_perf <- compute_observation_performance(model, constructs, scale_min, scale_max)
  sd_perf <- apply(obs_perf, 2, sd, na.rm = TRUE)

  # Build unstandardized path coefficient matrix
  B_std <- model$path_coef[constructs, constructs]
  B_unstd <- B_std

  for (i in seq_along(constructs)) {
    for (j in seq_along(constructs)) {
      if (B_std[i, j] != 0) {
        sd_from <- sd_perf[constructs[i]]
        sd_to <- sd_perf[constructs[j]]
        if (is.na(sd_from) || is.na(sd_to) ||
            sd_from < .Machine$double.eps) {
          B_unstd[i, j] <- NA_real_
        } else {
          B_unstd[i, j] <- B_std[i, j] * sd_to / sd_from
        }
      }
    }
  }

  compute_total_effects(B_unstd)
}

#' Classify constructs into cIPMA priority categories.
#' @noRd
classify_cipma_constructs <- function(importance, performance, nca_result) {
  constructs <- names(importance)

  imp_median <- median(importance, na.rm = TRUE)
  high_importance <- importance > imp_median

  necessary <- rep(FALSE, length(constructs))
  names(necessary) <- constructs

  if (!is.null(nca_result) && length(nca_result$necessary_predictors) > 0) {
    necessary[intersect(constructs, nca_result$necessary_predictors)] <- TRUE
  }

  priority <- ifelse(
    high_importance & necessary, "Top priority",
    ifelse(high_importance & !necessary, "Important driver",
      ifelse(!high_importance & necessary, "Bottleneck risk",
        "Low priority")))

  data.frame(
    Construct = constructs,
    Importance = round(importance, 4),
    Performance = round(performance, 2),
    High_Importance = high_importance,
    Necessary = necessary,
    Priority = priority,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# MAIN ENTRY POINTS
# =============================================================================

#' Importance-Performance Map Analysis (IPMA) for PLS-SEM
#'
#' \code{assess_ipma} conducts an Importance-Performance Map Analysis (IPMA)
#' for a given target construct. It computes each predecessor's importance
#' (unstandardized total effect on the target) and performance (rescaled
#' 0--100 mean construct score), identifying constructs with high importance
#' but low performance as priority areas for improvement.
#'
#' This is a convenience wrapper around \code{\link{assess_cipma}} with
#' \code{nca = FALSE}. Use \code{\link{assess_cipma}} instead if you also
#' want Necessary Condition Analysis (NCA) to produce a combined IPMA.
#'
#' @inheritParams assess_cipma
#'
#' @return An object of class \code{cipma_analysis} (see
#'   \code{\link{assess_cipma}} for details). The \code{nca} element will
#'   be \code{NULL}.
#'
#' @references
#' Ringle, C. M. & Sarstedt, M. (2016). Gain More Insight from Your PLS-SEM
#' Results: The Importance-Performance Map Analysis. Industrial Management &
#' Data Systems, 119(9), 1865-1886.
#'
#' @examples
#' library(seminr)
#' library(seminrExtras)
#'
#' mobi_mm <- constructs(
#'   composite("Image",        multi_items("IMAG", 1:5)),
#'   composite("Expectation",  multi_items("CUEX", 1:3)),
#'   composite("Value",        multi_items("PERV", 1:2)),
#'   composite("Satisfaction", multi_items("CUSA", 1:3)),
#'   composite("Loyalty",      multi_items("CUSL", 1:3))
#' )
#'
#' mobi_sm <- relationships(
#'   paths(from = "Image",       to = c("Expectation", "Satisfaction", "Loyalty")),
#'   paths(from = "Expectation", to = c("Value", "Satisfaction")),
#'   paths(from = "Value",       to = "Satisfaction"),
#'   paths(from = "Satisfaction",to = "Loyalty")
#' )
#'
#' mobi_pls <- estimate_pls(data = mobi,
#'                           measurement_model = mobi_mm,
#'                           structural_model  = mobi_sm)
#'
#' \donttest{
#' ipma_result <- assess_ipma(mobi_pls,
#'                             target = "Loyalty",
#'                             scale_min = 1,
#'                             scale_max = 10)
#' print(ipma_result)
#' plot(ipma_result)
#' }
#'
#' @seealso \code{\link{assess_cipma}} for cIPMA (IPMA + NCA)
#'
#' @export
assess_ipma <- function(seminr_model,
                         target,
                         scale_min = 1,
                         scale_max = 7,
                         seed = 123) {
  assess_cipma(seminr_model = seminr_model,
               target = target,
               scale_min = scale_min,
               scale_max = scale_max,
               nca = FALSE,
               seed = seed)
}

#' Combined Importance-Performance Map Analysis (cIPMA) for PLS-SEM
#'
#' \code{assess_cipma} conducts an Importance-Performance Map Analysis (IPMA)
#' and optionally combines it with Necessary Condition Analysis (NCA) to
#' produce a combined IPMA (cIPMA). The analysis identifies constructs that
#' are both important (high total effect) and necessary (without which the
#' outcome cannot reach high levels), following the 7-step cIPMA procedure
#' of Sarstedt et al. (2024).
#'
#' \strong{Importance} is measured by the unstandardized total effect of each
#' construct on the target. Unstandardized effects are obtained by scaling
#' standardized path coefficients using the standard deviations of rescaled
#' (0--100) construct scores.
#'
#' \strong{Performance} is the weighted average of rescaled indicator means,
#' where each indicator is rescaled from the original measurement scale
#' ([\code{scale_min}, \code{scale_max}]) to 0--100. Weights are the PLS
#' outer weights, which must all be positive for valid IPMA rescaling.
#'
#' Interaction constructs (e.g., \code{"X*W"}) are automatically excluded
#' from the IPMA because their performance is not meaningful on a 0--100
#' scale.
#'
#' When \code{nca = TRUE}, NCA is run on construct scores for each
#' predecessor--target pair, and constructs are classified into four
#' categories based on crossing importance (above/below median) with
#' necessity (NCA d >= 0.1 and p < 0.05).
#'
#' @param seminr_model An estimated SEMinR model from \code{estimate_pls()}.
#' @param target Name of the target (endogenous) construct for the IPMA.
#' @param scale_min Minimum of the measurement scale (default 1).
#' @param scale_max Maximum of the measurement scale (default 7).
#' @param nca Logical; if \code{TRUE} (default), run NCA to produce a cIPMA.
#' @param nca_ceilings Character vector of ceiling techniques for NCA
#'   (default \code{c("ce_fdh", "cr_fdh")}).
#' @param nca_test.rep Number of NCA permutation test repetitions (default 0).
#'   Set > 0 to obtain significance for the necessity classification.
#' @param nca_steps Number of bottleneck table steps (default 10).
#' @param seed Random seed for reproducibility (default 123).
#'
#' @return An object of class \code{cipma_analysis} containing:
#'   \item{importance_unstd}{Named numeric vector of unstandardized total effects}
#'   \item{importance_std}{Named numeric vector of standardized total effects}
#'   \item{performance}{Named numeric vector of construct performances (0--100)}
#'   \item{nca}{An \code{nca_analysis} object (NULL if \code{nca = FALSE})}
#'   \item{classification}{Data frame classifying each construct}
#'   \item{target}{Name of the target construct}
#'   \item{constructs}{Character vector of included constructs}
#'   \item{scale_range}{Numeric vector \code{c(scale_min, scale_max)}}
#'   \item{negative_weight_constructs}{Constructs with negative outer weights (if any)}
#'   \item{excluded_interactions}{Interaction constructs excluded from IPMA}
#'   \item{pls_model}{The original estimated seminr model}
#'
#' @references
#' Ringle, C. M. & Sarstedt, M. (2016). Gain More Insight from Your PLS-SEM
#' Results: The Importance-Performance Map Analysis. Industrial Management &
#' Data Systems, 119(9), 1865-1886.
#'
#' Sarstedt, M., Richter, N. F., Hauff, S. & Ringle, C. M. (2024). Combined
#' Importance-Performance Map Analysis (cIPMA): A SmartPLS 4 Tutorial.
#' Journal of Marketing Analytics, 12, 746-760.
#'
#' Hauff, S. et al. (2024). Importance and Performance in PLS-SEM and NCA:
#' Introducing the cIPMA. Journal of Retailing and Consumer Services, 78,
#' 103723.
#'
#' @examples
#' library(seminr)
#' library(seminrExtras)
#'
#' mobi_mm <- constructs(
#'   composite("Image",        multi_items("IMAG", 1:5)),
#'   composite("Expectation",  multi_items("CUEX", 1:3)),
#'   composite("Value",        multi_items("PERV", 1:2)),
#'   composite("Satisfaction", multi_items("CUSA", 1:3)),
#'   composite("Loyalty",      multi_items("CUSL", 1:3))
#' )
#'
#' mobi_sm <- relationships(
#'   paths(from = "Image",       to = c("Expectation", "Satisfaction", "Loyalty")),
#'   paths(from = "Expectation", to = c("Value", "Satisfaction")),
#'   paths(from = "Value",       to = "Satisfaction"),
#'   paths(from = "Satisfaction",to = "Loyalty")
#' )
#'
#' mobi_pls <- estimate_pls(data = mobi,
#'                           measurement_model = mobi_mm,
#'                           structural_model  = mobi_sm)
#'
#' \donttest{
#' cipma_result <- assess_cipma(mobi_pls,
#'                               target = "Loyalty",
#'                               scale_min = 1,
#'                               scale_max = 10)
#' print(cipma_result)
#' summary(cipma_result)
#' plot(cipma_result)
#' }
#'
#' @seealso \code{\link{assess_ipma}} for IPMA without NCA,
#'   \code{\link{assess_nca}} for standalone NCA analysis
#'
#' @export
assess_cipma <- function(seminr_model,
                          target,
                          scale_min = 1,
                          scale_max = 7,
                          nca = TRUE,
                          nca_ceilings = c("ce_fdh", "cr_fdh"),
                          nca_test.rep = 0,
                          nca_steps = 10,
                          seed = 123) {

  # ---------------------------------------------------------------------------
  # Step 1: Validate inputs
  # ---------------------------------------------------------------------------
  if (!validate_seminr_model(seminr_model, "assess_cipma")) return(NULL)

  construct_names <- colnames(seminr_model$construct_scores)

  if (!(target %in% construct_names)) {
    stop("target '", target, "' not found in model constructs: ",
         paste(construct_names, collapse = ", "), call. = FALSE)
  }

  if (!is.numeric(scale_min) || !is.numeric(scale_max) ||
      length(scale_min) != 1 || length(scale_max) != 1) {
    stop("scale_min and scale_max must be single numeric values.", call. = FALSE)
  }

  if (scale_min >= scale_max) {
    stop("scale_min must be less than scale_max.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Step 2: Identify constructs for IPMA
  # ---------------------------------------------------------------------------
  is_interaction <- vapply(construct_names, is_interaction_construct, logical(1))
  excluded_interactions <- construct_names[is_interaction]
  if (length(excluded_interactions) > 0) {
    message("Excluding interaction term(s) from IPMA: ",
            paste(excluded_interactions, collapse = ", "))
  }

  # All non-interaction constructs that appear in path_coef
  pc_names <- colnames(seminr_model$path_coef)
  all_constructs <- intersect(construct_names[!is_interaction], pc_names)

  # IPMA constructs = all except target
  ipma_constructs <- setdiff(all_constructs, target)

  if (length(ipma_constructs) == 0) {
    stop("No constructs available for IPMA (all excluded or only target exists).",
         call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Step 3: Check positive outer weights
  # ---------------------------------------------------------------------------
  neg_weight_constructs <- check_positive_weights(seminr_model, ipma_constructs)
  if (length(neg_weight_constructs) > 0) {
    warning("Negative outer weights detected for: ",
            paste(neg_weight_constructs, collapse = ", "),
            ". IPMA performance rescaling may be unreliable for these constructs. ",
            "Consider reverse-coding indicators or checking the model specification ",
            "(Ringle & Sarstedt, 2016).", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Step 4: Compute performance (0-100 rescaled construct scores)
  # ---------------------------------------------------------------------------
  performance <- compute_ipma_performance(seminr_model, ipma_constructs,
                                           scale_min, scale_max)

  # ---------------------------------------------------------------------------
  # Step 5: Compute importance (total effects on target)
  # ---------------------------------------------------------------------------
  # Standardized total effects
  T_std <- compute_total_effects(seminr_model$path_coef[all_constructs, all_constructs])
  importance_std <- T_std[ipma_constructs, target]

  # Unstandardized total effects
  T_unstd <- compute_unstd_total_effects(seminr_model, all_constructs,
                                          scale_min, scale_max)
  importance_unstd <- T_unstd[ipma_constructs, target]

  # Filter to constructs with non-zero total effect on target
  has_effect <- abs(importance_std) > .Machine$double.eps
  if (!any(has_effect)) {
    stop("No constructs have a non-zero total effect on '", target, "'.",
         call. = FALSE)
  }

  ipma_constructs <- ipma_constructs[has_effect]
  performance <- performance[has_effect]
  importance_std <- importance_std[has_effect]
  importance_unstd <- importance_unstd[has_effect]

  # ---------------------------------------------------------------------------
  # Step 6: Run NCA (for cIPMA)
  # ---------------------------------------------------------------------------
  nca_result <- NULL
  if (nca) {
    nca_result <- assess_nca(seminr_model,
                              target = target,
                              predictors = ipma_constructs,
                              ceilings = nca_ceilings,
                              test.rep = nca_test.rep,
                              steps = nca_steps,
                              seed = seed)
  }

  # ---------------------------------------------------------------------------
  # Step 7: Classify constructs
  # ---------------------------------------------------------------------------
  classification <- classify_cipma_constructs(importance_unstd, performance,
                                               nca_result)

  # ---------------------------------------------------------------------------
  # Step 8: Build and return result
  # ---------------------------------------------------------------------------
  result <- list(
    importance_unstd = importance_unstd,
    importance_std = importance_std,
    performance = performance,
    nca = nca_result,
    classification = classification,
    target = target,
    constructs = ipma_constructs,
    scale_range = c(scale_min, scale_max),
    negative_weight_constructs = neg_weight_constructs,
    excluded_interactions = excluded_interactions,
    pls_model = seminr_model
  )

  class(result) <- c("cipma_analysis", class(result))
  result
}

# =============================================================================
# S3 METHODS
# =============================================================================

#' @export
print.cipma_analysis <- function(x, ...) {
  has_nca <- !is.null(x$nca)

  if (has_nca) {
    cat("Combined Importance-Performance Map Analysis (cIPMA)\n")
    cat("=====================================================\n")
  } else {
    cat("Importance-Performance Map Analysis (IPMA)\n")
    cat("============================================\n")
  }

  cat("Target:", x$target, "\n")
  cat("Constructs:", paste(x$constructs, collapse = ", "), "\n")
  cat("Scale range:", x$scale_range[1], "-", x$scale_range[2], "\n")
  cat("Observations:", nrow(x$pls_model$construct_scores), "\n")

  if (length(x$excluded_interactions) > 0) {
    cat("Excluded (interaction):", paste(x$excluded_interactions, collapse = ", "), "\n")
  }
  if (length(x$negative_weight_constructs) > 0) {
    cat("Warning - negative weights:", paste(x$negative_weight_constructs, collapse = ", "), "\n")
  }

  cat("\nImportance-Performance Results:\n")
  tbl <- data.frame(
    Construct = x$constructs,
    `Unstd. Total Effect` = round(x$importance_unstd, 4),
    `Std. Total Effect` = round(x$importance_std, 4),
    Performance = round(x$performance, 2),
    check.names = FALSE,
    row.names = NULL
  )
  print(tbl, row.names = FALSE, ...)

  if (has_nca) {
    cat("\nNecessary Conditions (NCA):\n")
    nca_tbl <- data.frame(
      Construct = x$constructs,
      check.names = FALSE,
      row.names = NULL
    )
    for (ceil in colnames(x$nca$effect_sizes)) {
      nca_tbl[[paste0("d (", ceil, ")")]] <- round(x$nca$effect_sizes[x$constructs, ceil], 4)
    }
    necessary_flag <- ifelse(x$classification$Necessary, "Yes", "No")
    nca_tbl[["Necessary"]] <- necessary_flag
    print(nca_tbl, row.names = FALSE, ...)
  }

  if (has_nca) {
    cat("\ncIPMA Classification:\n")
  } else {
    cat("\nIPMA Classification:\n")
  }
  class_tbl <- x$classification[, c("Construct", "Priority")]
  print(class_tbl, row.names = FALSE, ...)

  invisible(x)
}

#' @export
summary.cipma_analysis <- function(object, ...) {
  result <- list(
    target = object$target,
    constructs = object$constructs,
    scale_range = object$scale_range,
    n_obs = nrow(object$pls_model$construct_scores),
    importance_unstd = object$importance_unstd,
    importance_std = object$importance_std,
    performance = object$performance,
    classification = object$classification,
    nca = object$nca,
    negative_weight_constructs = object$negative_weight_constructs,
    excluded_interactions = object$excluded_interactions
  )

  class(result) <- c("summary.cipma_analysis", class(result))
  result
}

#' @export
print.summary.cipma_analysis <- function(x, ...) {
  has_nca <- !is.null(x$nca)

  if (has_nca) {
    cat("Combined Importance-Performance Map Analysis (cIPMA) Summary\n")
    cat("=============================================================\n")
  } else {
    cat("Importance-Performance Map Analysis (IPMA) Summary\n")
    cat("====================================================\n")
  }

  cat("Target:", x$target, "\n")
  cat("Scale range:", x$scale_range[1], "-", x$scale_range[2], "\n")
  cat("Observations:", x$n_obs, "\n")

  if (length(x$excluded_interactions) > 0) {
    cat("Excluded interactions:", paste(x$excluded_interactions, collapse = ", "), "\n")
  }
  if (length(x$negative_weight_constructs) > 0) {
    cat("Negative weight constructs:", paste(x$negative_weight_constructs, collapse = ", "), "\n")
  }

  cat("\nImportance (Unstandardized Total Effects on", x$target, "):\n")
  print(round(x$importance_unstd, 4), ...)

  cat("\nImportance (Standardized Total Effects on", x$target, "):\n")
  print(round(x$importance_std, 4), ...)

  cat("\nPerformance (0-100 rescaled):\n")
  print(round(x$performance, 2), ...)

  if (has_nca) {
    cat("\nNCA Effect Sizes:\n")
    print(round(x$nca$effect_sizes[x$constructs, , drop = FALSE], 4), ...)

    if (!is.null(x$nca$significance)) {
      cat("\nNCA Permutation p-values:\n")
      print(round(x$nca$significance[x$constructs, , drop = FALSE], 4), ...)
    }

    if (length(x$nca$necessary_predictors) > 0) {
      cat("\nNecessary conditions:", paste(x$nca$necessary_predictors, collapse = ", "), "\n")
    } else {
      cat("\nNo necessary conditions identified.\n")
    }

    # Print bottleneck for necessary constructs
    if (length(x$nca$necessary_predictors) > 0) {
      for (ceil in names(x$nca$bottleneck)) {
        cat("\nBottleneck table (", ceil, "):\n", sep = "")
        print(x$nca$bottleneck[[ceil]], ...)
      }
    }
  }

  cat("\nConstruct Classification:\n")
  print(x$classification, row.names = FALSE, ...)

  invisible(x)
}

#' Plot cIPMA Results
#'
#' Produces the Importance-Performance Map with optional NCA overlay.
#' When \code{type = "cipma"} (default), necessary conditions are
#' highlighted with filled red circles; other constructs use open blue circles.
#' When \code{type = "ipma"}, all constructs use the same symbol.
#'
#' @param x A \code{cipma_analysis} object from \code{assess_cipma()}.
#' @param type One of \code{"cipma"} (default, with NCA overlay) or
#'   \code{"ipma"} (standard IPMA without NCA distinction).
#' @param importance_metric One of \code{"unstandardized"} (default) or
#'   \code{"standardized"} to choose which total effect is plotted on the
#'   x-axis.
#' @param ... Additional arguments passed to \code{plot()}.
#'
#' @export
plot.cipma_analysis <- function(x,
                                 type = c("cipma", "ipma"),
                                 importance_metric = c("unstandardized", "standardized"),
                                 ...) {
  type <- match.arg(type)
  importance_metric <- match.arg(importance_metric)

  if (type == "cipma" && is.null(x$nca)) {
    message("No NCA results available. Showing standard IPMA plot.")
    type <- "ipma"
  }

  imp <- if (importance_metric == "unstandardized") x$importance_unstd else x$importance_std
  perf <- x$performance
  constructs <- x$constructs

  # Point styles
  if (type == "cipma") {
    necessary <- x$classification$Necessary
    pch <- ifelse(necessary, 19, 1)
    col <- ifelse(necessary, "red", "steelblue")
  } else {
    pch <- rep(19, length(constructs))
    col <- rep("steelblue", length(constructs))
  }

  # Margins
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mar = c(5, 5, 4, 2))

  # Axis limits
  imp_pad <- max(0.05, diff(range(imp, na.rm = TRUE)) * 0.15)
  xlim <- range(imp, na.rm = TRUE) + c(-imp_pad, imp_pad)
  ylim <- c(max(0, min(perf, na.rm = TRUE) - 5),
            min(100, max(perf, na.rm = TRUE) + 5))

  # X-axis label
  xlab <- if (importance_metric == "unstandardized") {
    "Importance (unstandardized total effect)"
  } else {
    "Importance (standardized total effect)"
  }

  # Title
  main_title <- if (type == "cipma") {
    paste("cIPMA:", x$target)
  } else {
    paste("IPMA:", x$target)
  }

  plot(imp, perf, pch = pch, col = col, cex = 1.8,
       xlim = xlim, ylim = ylim,
       xlab = xlab, ylab = "Performance (0-100)",
       main = main_title, bty = "n", ...)

  # Reference lines at means
  abline(v = mean(imp, na.rm = TRUE), col = "gray70", lty = 2)
  abline(h = mean(perf, na.rm = TRUE), col = "gray70", lty = 2)

  # Construct labels
  text(imp, perf, labels = constructs, pos = 3, cex = 0.8)

  # Legend
  if (type == "cipma") {
    legend("bottomleft",
           legend = c("Necessary + sufficient", "Sufficient only"),
           pch = c(19, 1), col = c("red", "steelblue"),
           pt.cex = 1.5, cex = 0.8, bty = "n")
  }
}
