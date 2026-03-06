# =============================================================================
# feature_nca.R - Necessary Condition Analysis (NCA) for PLS-SEM
# =============================================================================
# Implements CE-FDH and CR-FDH ceiling analysis internally. The NCA package
# (Dul, 2016, 2020) is optional -- only required for ceiling techniques
# beyond CE-FDH and CR-FDH, or for NCA's native scatter plots.
#
# References:
# - Dul, J. (2016). "Necessary Condition Analysis (NCA): Logic and
#   Methodology of 'Necessary but Not Sufficient' Causality."
#   Organizational Research Methods, 19(1), 10-52.
# - Richter, N. F., Schubring, S., Hauff, S., Ringle, C. M., &
#   Sarstedt, M. (2020). "When predictors of outcomes are necessary:
#   guidelines for the combined use of PLS-SEM and NCA."
#   Industrial Management & Data Systems, 120(12), 2243-2267.
# =============================================================================

# =============================================================================
# INTERNAL NCA ALGORITHMS
# =============================================================================
# Self-contained CE-FDH and CR-FDH implementations. These cover the two
# most commonly used ceiling techniques in PLS-SEM NCA applications.

INTERNAL_CEILINGS <- c("ce_fdh", "cr_fdh")

#' Compute CE-FDH (Ceiling Envelopment - Free Disposal Hull) ceiling.
#' Returns sorted unique x values and their non-decreasing ceiling y values.
#' The ceiling at any x is the maximum y among all observations with x_obs <= x.
#' @noRd
compute_ce_fdh <- function(x, y) {
  ux <- sort(unique(x))
  max_y_at_x <- vapply(ux, function(xi) max(y[x == xi]), numeric(1))
  list(x = ux, ceiling_y = cummax(max_y_at_x))
}

#' Get CE-FDH peer coordinates (points defining the ceiling envelope).
#' A peer is a unique-x level whose max y equals the running cummax,
#' meaning it actually pushes the ceiling upward.
#' @noRd
get_ce_fdh_peers <- function(x, y) {
  ux <- sort(unique(x))
  max_y_at_x <- vapply(ux, function(xi) max(y[x == xi]), numeric(1))
  cy <- cummax(max_y_at_x)
  is_peer <- max_y_at_x == cy
  data.frame(x = ux[is_peer], y = cy[is_peer])
}

#' Compute CE-FDH effect size d = ceiling_zone / scope.
#' The ceiling zone is the empty area above the step-function ceiling
#' within the scope rectangle [min(x),max(x)] x [min(y),max(y)].
#' @noRd
ce_fdh_effect_size <- function(x, y) {
  x_range <- range(x)
  y_range <- range(y)
  scope <- diff(x_range) * diff(y_range)
  if (scope < .Machine$double.eps) return(0)

  ceil <- compute_ce_fdh(x, y)
  ux <- ceil$x
  cy <- ceil$ceiling_y
  k <- length(ux)
  if (k < 2) return(0)

  # For each interval [ux[i], ux[i+1]), the ceiling is at cy[i].
  # Area above = width * (max_y - cy[i]).
  ceiling_zone <- sum(diff(ux) * (y_range[2] - cy[-k]))
  ceiling_zone / scope
}

#' Compute area above a line y = a + b*x within a rectangle.
#' Handles all cases: positive/negative/zero slope, line partially
#' or fully outside the scope.
#' @noRd
line_ceiling_zone <- function(a, b, x_min, x_max, y_min, y_max) {
  scope_height <- y_max - y_min

  if (abs(b) < .Machine$double.eps) {
    if (a >= y_max) return(0)
    if (a <= y_min) return((x_max - x_min) * scope_height)
    return((x_max - x_min) * (y_max - a))
  }

  # Break points where the line crosses y_min or y_max
  x_cross <- c((y_min - a) / b, (y_max - a) / b)
  inner <- x_cross[x_cross > x_min & x_cross < x_max]
  breaks <- sort(unique(c(x_min, inner, x_max)))

  total <- 0
  for (i in seq_len(length(breaks) - 1)) {
    xl <- breaks[i]
    xr <- breaks[i + 1]
    line_at_mid <- a + b * (xl + xr) / 2

    if (line_at_mid >= y_max) {
      # Line above scope top: no ceiling zone
    } else if (line_at_mid <= y_min) {
      # Line below scope bottom: full ceiling zone
      total <- total + (xr - xl) * scope_height
    } else {
      # Line within scope: integrate y_max - (a + bx)
      total <- total + (y_max - a) * (xr - xl) - b / 2 * (xr^2 - xl^2)
    }
  }

  total
}

#' Compute CR-FDH (Ceiling Regression - Free Disposal Hull) effect size.
#' Fits OLS through CE-FDH peers; falls back to CE-FDH if < 2 peers.
#' @noRd
cr_fdh_effect_size <- function(x, y) {
  x_range <- range(x)
  y_range <- range(y)
  scope <- diff(x_range) * diff(y_range)
  if (scope < .Machine$double.eps) return(0)

  peers <- get_ce_fdh_peers(x, y)
  if (nrow(peers) < 2) return(ce_fdh_effect_size(x, y))

  cf <- coef(lm(y ~ x, data = peers))
  a <- unname(cf[1])
  b <- unname(cf[2])

  ceiling_zone <- line_ceiling_zone(a, b, x_range[1], x_range[2],
                                     y_range[1], y_range[2])
  max(0, min(ceiling_zone / scope, 1))
}

#' Dispatch to the correct effect size function.
#' @noRd
nca_effect_size <- function(x, y, ceiling_type) {
  switch(ceiling_type,
    ce_fdh = ce_fdh_effect_size(x, y),
    cr_fdh = cr_fdh_effect_size(x, y),
    stop("Internal NCA supports: ", paste(INTERNAL_CEILINGS, collapse = ", "),
         ". Install the NCA package for '", ceiling_type, "'.", call. = FALSE)
  )
}

#' Permutation test for NCA effect size significance.
#' Returns p-value as proportion of permuted d >= observed d.
#' @noRd
nca_permutation_test <- function(x, y, ceiling_type, observed_d, n_perm) {
  d_perm <- vapply(seq_len(n_perm), function(i) {
    nca_effect_size(x, sample(y), ceiling_type)
  }, numeric(1))
  mean(d_perm >= observed_d)
}

#' Compute bottleneck column for one predictor: required X percentage
#' at each Y level. NA indicates the condition is not necessary.
#' @noRd
compute_bottleneck_column <- function(x, y, ceiling_type, steps) {
  x_range <- range(x)
  y_range <- range(y)
  y_levels <- seq(0, 100, length.out = steps + 1)
  y_targets <- y_range[1] + y_levels / 100 * diff(y_range)

  if (diff(x_range) < .Machine$double.eps) {
    return(rep(NA_real_, length(y_levels)))
  }

  if (ceiling_type == "ce_fdh") {
    ceil <- compute_ce_fdh(x, y)
    x_needed <- vapply(y_targets, function(yt) {
      idx <- which(ceil$ceiling_y >= yt)
      if (length(idx) == 0) return(NA_real_)
      ceil$x[idx[1]]
    }, numeric(1))
  } else if (ceiling_type == "cr_fdh") {
    peers <- get_ce_fdh_peers(x, y)
    if (nrow(peers) < 2) {
      return(compute_bottleneck_column(x, y, "ce_fdh", steps))
    }
    cf <- coef(lm(y ~ x, data = peers))
    a <- unname(cf[1])
    b <- unname(cf[2])
    if (abs(b) < .Machine$double.eps) {
      return(rep(NA_real_, length(y_levels)))
    }
    x_needed <- (y_targets - a) / b
  } else {
    return(rep(NA_real_, length(y_levels)))
  }

  x_pct <- (x_needed - x_range[1]) / diff(x_range) * 100
  x_pct[x_pct < 0] <- NA_real_
  x_pct[x_pct > 100] <- NA_real_
  round(x_pct, 1)
}

#' Run complete internal NCA analysis for all predictors and ceilings.
#' @noRd
run_nca_internal <- function(data, x_names, y_name, ceilings, test.rep, steps) {
  y <- data[[y_name]]

  effect_sizes <- matrix(NA_real_, nrow = length(x_names), ncol = length(ceilings),
                          dimnames = list(x_names, ceilings))
  p_values <- if (test.rep > 0) {
    matrix(NA_real_, nrow = length(x_names), ncol = length(ceilings),
           dimnames = list(x_names, ceilings))
  } else NULL

  bottlenecks <- list()

  for (ceil in ceilings) {
    y_levels <- seq(0, 100, length.out = steps + 1)
    bn <- data.frame(Y = y_levels)
    colnames(bn)[1] <- y_name

    for (pred in x_names) {
      x_vec <- data[[pred]]
      d <- nca_effect_size(x_vec, y, ceil)
      effect_sizes[pred, ceil] <- d

      if (test.rep > 0) {
        p_values[pred, ceil] <- nca_permutation_test(x_vec, y, ceil, d, test.rep)
      }

      bn[[pred]] <- compute_bottleneck_column(x_vec, y, ceil, steps)
    }

    bottlenecks[[ceil]] <- bn
  }

  list(effect_sizes = effect_sizes, p_values = p_values, bottlenecks = bottlenecks)
}

# =============================================================================
# NCA PACKAGE FALLBACK
# =============================================================================
# Used only when ceiling techniques beyond CE-FDH/CR-FDH are requested.

#' @noRd
check_nca_installed <- function() {
  if (!requireNamespace("NCA", quietly = TRUE)) {
    stop("Package 'NCA' is required for ceiling techniques beyond CE-FDH and CR-FDH. ",
         "Install it with: install.packages('NCA')",
         call. = FALSE)
  }
}

#' @noRd
format_nca_effects <- function(nca_raw, predictors, ceilings) {
  mat <- matrix(NA_real_, nrow = length(predictors), ncol = length(ceilings),
                dimnames = list(predictors, ceilings))
  for (pred in predictors) {
    params <- nca_raw$summaries[[pred]]$params
    if (!is.null(params)) {
      for (ceil in ceilings) {
        if (ceil %in% colnames(params) && "Effect size" %in% rownames(params)) {
          mat[pred, ceil] <- params["Effect size", ceil]
        }
      }
    }
  }
  mat
}

#' @noRd
format_nca_significance <- function(nca_raw, predictors, ceilings) {
  has_pvals <- any(vapply(predictors, function(pred) {
    params <- nca_raw$summaries[[pred]]$params
    !is.null(params) && "p-value" %in% rownames(params) &&
      any(!is.na(params["p-value", ]))
  }, logical(1)))
  if (!has_pvals) return(NULL)

  mat <- matrix(NA_real_, nrow = length(predictors), ncol = length(ceilings),
                dimnames = list(predictors, ceilings))
  for (pred in predictors) {
    params <- nca_raw$summaries[[pred]]$params
    if (!is.null(params)) {
      for (ceil in ceilings) {
        if (ceil %in% colnames(params) && "p-value" %in% rownames(params)) {
          mat[pred, ceil] <- params["p-value", ceil]
        }
      }
    }
  }
  mat
}

#' @noRd
format_nca_bottleneck <- function(nca_raw, ceilings) {
  bn <- list()
  for (ceil in ceilings) {
    if (!is.null(nca_raw$bottlenecks[[ceil]])) {
      bn[[ceil]] <- nca_raw$bottlenecks[[ceil]]
    }
  }
  bn
}

# =============================================================================
# PREDICTOR EXTRACTION
# =============================================================================

#' Extract direct predictor construct names for a given target from the
#' structural model matrix. unname() strips matrix row names that would
#' cause matching issues downstream.
#' @noRd
get_direct_predictors <- function(model, target) {
  sm <- model$smMatrix
  unname(sm[sm[, "target"] == target, "source"])
}

# =============================================================================
# INPUT VALIDATION
# =============================================================================

#' Shared input validation for assess_nca() and assess_nca_esse().
#' Validates test.rep, target, and predictors. Returns the predictor vector
#' (auto-detected or validated). Stops on invalid input.
#' @noRd
validate_nca_inputs <- function(seminr_model, target, predictors, test.rep) {
  if (!is.numeric(test.rep) || length(test.rep) != 1 || test.rep < 0 ||
      test.rep != as.integer(test.rep)) {
    stop("test.rep must be a non-negative integer.", call. = FALSE)
  }

  construct_names <- colnames(seminr_model$construct_scores)
  if (!(target %in% construct_names)) {
    stop("target '", target, "' not found in model constructs: ",
         paste(construct_names, collapse = ", "), call. = FALSE)
  }

  if (is.null(predictors)) {
    predictors <- get_direct_predictors(seminr_model, target)
    if (length(predictors) == 0) {
      stop("No direct predictors found for target '", target,
           "' in the structural model.", call. = FALSE)
    }
  } else {
    invalid <- setdiff(predictors, construct_names)
    if (length(invalid) > 0) {
      stop("Predictor(s) not found in model constructs: ",
           paste(invalid, collapse = ", "), call. = FALSE)
    }
  }

  predictors
}

# =============================================================================
# MAIN ENTRY POINT
# =============================================================================

#' Necessary Condition Analysis (NCA) for PLS-SEM Models
#'
#' `assess_nca` conducts a Necessary Condition Analysis on construct scores
#' from an estimated PLS-SEM model. It identifies whether predictors are
#' necessary conditions for achieving a given level of the target construct.
#'
#' NCA complements PLS-SEM's sufficiency-based regression by testing whether
#' a certain level of a predictor is \emph{necessary} (but not sufficient)
#' for a certain level of the outcome. A predictor is considered a necessary
#' condition when it has an effect size (d) >= 0.1 and a significant p-value
#' (p < 0.05).
#'
#' The function extracts construct scores from the seminr model, auto-detects
#' direct predictors of the target from the structural model, and computes
#' ceiling envelopes using the specified techniques.
#'
#' CE-FDH and CR-FDH ceilings are computed internally. For other ceiling
#' techniques, the \pkg{NCA} package must be installed.
#'
#' @param seminr_model An estimated SEMinR model from \code{estimate_pls()}.
#' @param target Name of the endogenous (outcome) construct (character).
#' @param predictors Optional character vector of predictor construct names.
#'   If \code{NULL} (default), auto-detected from the structural model as
#'   direct predictors of \code{target}.
#' @param ceilings Character vector of ceiling techniques to use
#'   (default \code{c("ce_fdh", "cr_fdh")}). CE-FDH and CR-FDH are computed
#'   internally; other techniques require the \pkg{NCA} package.
#' @param test.rep Number of permutation test repetitions for significance
#'   testing (default 1000). Set to 0 to skip significance testing.
#' @param steps Number of steps in the bottleneck table (default 10).
#' @param seed Random seed for reproducibility (default 123).
#' @param ... Additional arguments passed to \code{NCA::nca_analysis()} when
#'   using ceiling techniques not supported internally.
#'
#' @return An object of class \code{nca_analysis} containing:
#'   \item{nca_raw}{The raw result from \code{NCA::nca_analysis()} (NULL when
#'     using internal ceilings)}
#'   \item{effect_sizes}{Matrix of effect sizes (d) per predictor and ceiling}
#'   \item{significance}{Matrix of p-values per predictor and ceiling (NULL if test.rep = 0)}
#'   \item{bottleneck}{List of bottleneck tables, one per ceiling technique}
#'   \item{pls_model}{The original estimated seminr model}
#'   \item{target}{Name of the target construct}
#'   \item{predictors}{Character vector of predictor names used}
#'   \item{ceilings}{Character vector of ceiling techniques used}
#'
#' @references
#' Dul, J. (2016). Necessary Condition Analysis (NCA): Logic and Methodology
#' of 'Necessary but Not Sufficient' Causality. Organizational Research
#' Methods, 19(1), 10-52.
#'
#' Richter, N. F., Schubring, S., Hauff, S., Ringle, C. M., & Sarstedt, M.
#' (2020). When predictors of outcomes are necessary: guidelines for the
#' combined use of PLS-SEM and NCA. Industrial Management & Data Systems,
#' 120(12), 2243-2267.
#'
#' @examples
#' library(seminr)
#' library(seminrExtras)
#'
#' mobi_mm <- constructs(
#'   composite("Image",        multi_items("IMAG", 1:5)),
#'   composite("Value",        multi_items("PERV", 1:2)),
#'   composite("Satisfaction", multi_items("CUSA", 1:3)),
#'   composite("Loyalty",      multi_items("CUSL", 1:3))
#' )
#'
#' mobi_sm <- relationships(
#'   paths(from = c("Image", "Value"), to = "Satisfaction"),
#'   paths(from = "Satisfaction", to = "Loyalty")
#' )
#'
#' mobi_pls <- estimate_pls(data = mobi,
#'                           measurement_model = mobi_mm,
#'                           structural_model  = mobi_sm)
#'
#' \donttest{
#' nca_result <- assess_nca(mobi_pls,
#'                           target = "Satisfaction",
#'                           test.rep = 100)
#' print(nca_result)
#' summary(nca_result)
#' plot(nca_result, type = "effects")
#' }
#'
#' @export
assess_nca <- function(seminr_model,
                        target,
                        predictors = NULL,
                        ceilings = c("ce_fdh", "cr_fdh"),
                        test.rep = 1000,
                        steps = 10,
                        seed = 123,
                        ...) {

  # ---------------------------------------------------------------------------
  # Step 1: Validate inputs
  # ---------------------------------------------------------------------------
  if (!validate_seminr_model(seminr_model, "assess_nca")) {
    return(NULL)
  }
  predictors <- validate_nca_inputs(seminr_model, target, predictors, test.rep)

  # ---------------------------------------------------------------------------
  # Step 2: Run NCA analysis on construct scores
  # ---------------------------------------------------------------------------
  scores <- as.data.frame(seminr_model$construct_scores)
  use_internal <- all(ceilings %in% INTERNAL_CEILINGS)

  set.seed(seed)

  if (use_internal) {
    nca_result <- run_nca_internal(scores, predictors, target, ceilings,
                                    test.rep, steps)
    effect_sizes <- nca_result$effect_sizes
    significance <- nca_result$p_values
    bottleneck   <- nca_result$bottlenecks
    nca_raw      <- NULL
  } else {
    # Ceiling techniques beyond CE-FDH/CR-FDH require the NCA package
    check_nca_installed()
    nca_raw <- suppressMessages(
      NCA::nca_analysis(
        data = scores, x = predictors, y = target,
        ceilings = ceilings, test.rep = test.rep,
        steps = steps, ...
      )
    )
    effect_sizes <- format_nca_effects(nca_raw, predictors, ceilings)
    significance <- format_nca_significance(nca_raw, predictors, ceilings)
    bottleneck   <- format_nca_bottleneck(nca_raw, ceilings)
  }

  # ---------------------------------------------------------------------------
  # Step 3: Format and identify necessary conditions
  # ---------------------------------------------------------------------------
  comment(effect_sizes) <- "NCA Effect Sizes (d >= 0.1 indicates a necessary condition)"
  class(effect_sizes) <- c("table_output", class(effect_sizes))

  if (!is.null(significance)) {
    comment(significance) <- "NCA Significance (permutation test p-values)"
    class(significance) <- c("table_output", class(significance))
  }

  # Identify necessary conditions: d >= 0.1 AND p < 0.05 (Dul, 2016)
  necessary <- identify_necessary(effect_sizes, significance)
  necessary_predictors <- rownames(effect_sizes)[apply(necessary, 1, any)]

  result <- list(
    nca_raw      = nca_raw,
    effect_sizes = effect_sizes,
    significance = significance,
    bottleneck   = bottleneck,
    necessary_predictors = necessary_predictors,
    pls_model    = seminr_model,
    target       = target,
    predictors   = predictors,
    ceilings     = ceilings
  )

  class(result) <- c("nca_analysis", class(result))
  result
}

# =============================================================================
# NECESSARY CONDITION IDENTIFICATION
# =============================================================================

#' Identify necessary conditions per Dul (2016): a predictor is necessary
#' when effect size d >= 0.1 AND permutation p-value < 0.05.
#' When significance is NULL (no permutation test), uses d threshold only.
#' @noRd
identify_necessary <- function(effect_sizes, significance, d_threshold = 0.1,
                                p_threshold = 0.05) {
  if (is.null(significance)) {
    return(effect_sizes >= d_threshold)
  }
  (effect_sizes >= d_threshold) & (significance < p_threshold)
}

# =============================================================================
# S3 METHODS
# =============================================================================

#' @export
print.nca_analysis <- function(x, ...) {
  cat("Necessary Condition Analysis (NCA)\n")
  cat("==================================\n")
  cat("Target:", x$target, "\n")
  cat("Predictors:", paste(x$predictors, collapse = ", "), "\n")
  cat("Ceilings:", paste(x$ceilings, collapse = ", "), "\n")
  cat("Observations:", nrow(x$pls_model$construct_scores), "\n\n")

  cat("Effect Sizes (d):\n")
  print(round(x$effect_sizes, 4), ...)
  cat("\n")

  if (!is.null(x$significance)) {
    cat("Permutation p-values:\n")
    print(round(x$significance, 4), ...)
    cat("\n")
  }

  if (length(x$necessary_predictors) > 0) {
    cat("Necessary conditions (d >= 0.1, p < 0.05):",
        paste(x$necessary_predictors, collapse = ", "), "\n")
  } else {
    cat("No necessary conditions identified (d >= 0.1, p < 0.05)\n")
  }

  invisible(x)
}

#' @export
summary.nca_analysis <- function(object, ...) {
  result <- list(
    target       = object$target,
    predictors   = object$predictors,
    ceilings     = object$ceilings,
    n_obs        = nrow(object$pls_model$construct_scores),
    effect_sizes = object$effect_sizes,
    significance = object$significance,
    bottleneck   = object$bottleneck,
    necessary_predictors = object$necessary_predictors
  )

  class(result) <- c("summary.nca_analysis", class(result))
  result
}

#' @export
print.summary.nca_analysis <- function(x, ...) {
  cat("Necessary Condition Analysis (NCA) Summary\n")
  cat("============================================\n")
  cat("Target:", x$target, "\n")
  cat("Predictors:", paste(x$predictors, collapse = ", "), "\n")
  cat("Ceilings:", paste(x$ceilings, collapse = ", "), "\n")
  cat("Observations:", x$n_obs, "\n\n")

  cat("Effect Sizes (d):\n")
  print(round(x$effect_sizes, 4), ...)
  cat("\n")

  if (!is.null(x$significance)) {
    cat("Permutation p-values:\n")
    print(round(x$significance, 4), ...)
    cat("\n")
  }

  if (length(x$necessary_predictors) > 0) {
    cat("Necessary conditions (d >= 0.1, p < 0.05):",
        paste(x$necessary_predictors, collapse = ", "), "\n\n")
  } else {
    cat("No necessary conditions identified (d >= 0.1, p < 0.05)\n\n")
  }

  for (ceil in names(x$bottleneck)) {
    cat("Bottleneck table (", ceil, "):\n", sep = "")
    print(x$bottleneck[[ceil]], ...)
    cat("\n")
  }

  invisible(x)
}

#' Plot NCA Results
#'
#' @param x An \code{nca_analysis} object from \code{assess_nca()}.
#' @param type One of \code{"scatter"} (ceiling line scatter plots) or
#'   \code{"effects"} (bar plot of effect sizes).
#' @param ... Additional arguments passed to the underlying plot function.
#'
#' @export
plot.nca_analysis <- function(x, type = c("scatter", "effects"), ...) {
  type <- match.arg(type)

  switch(type,
    scatter = plot_nca_scatter(x, ...),
    effects = plot_nca_effects(x, ...)
  )
}

#' Scatter plot with ceiling lines. Uses NCA package if available and
#' nca_raw exists, otherwise draws scatter with internal ceiling lines.
#' @noRd
plot_nca_scatter <- function(nca_result, ...) {
  if (!is.null(nca_result$nca_raw) && requireNamespace("NCA", quietly = TRUE)) {
    plot(nca_result$nca_raw, ...)
    return(invisible(NULL))
  }

  scores <- nca_result$pls_model$construct_scores
  n_pred <- length(nca_result$predictors)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  if (n_pred > 1) {
    ncol_plot <- min(n_pred, 3)
    nrow_plot <- ceiling(n_pred / ncol_plot)
    par(mfrow = c(nrow_plot, ncol_plot))
  }

  for (pred in nca_result$predictors) {
    x <- scores[, pred]
    y <- scores[, nca_result$target]

    plot(x, y, pch = 19, col = adjustcolor("black", 0.4),
         xlab = pred, ylab = nca_result$target,
         main = paste("NCA:", pred, "->", nca_result$target),
         bty = "n")

    # CE-FDH ceiling: step function
    if ("ce_fdh" %in% nca_result$ceilings) {
      ceil <- compute_ce_fdh(x, y)
      lines(ceil$x, ceil$ceiling_y, type = "s", col = "red", lwd = 2)
    }

    # CR-FDH ceiling: regression line through peers
    if ("cr_fdh" %in% nca_result$ceilings) {
      peers <- get_ce_fdh_peers(x, y)
      if (nrow(peers) >= 2) {
        fit <- lm(y ~ x, data = peers)
        x_seq <- seq(min(x), max(x), length.out = 100)
        y_hat <- coef(fit)[1] + coef(fit)[2] * x_seq
        # Clip to scope
        y_hat <- pmin(pmax(y_hat, min(y)), max(y))
        lines(x_seq, y_hat, col = "blue", lwd = 2, lty = 2)
      }
    }

    # Legend
    ceil_names <- intersect(nca_result$ceilings, INTERNAL_CEILINGS)
    if (length(ceil_names) > 0) {
      cols <- c(ce_fdh = "red", cr_fdh = "blue")[ceil_names]
      ltys <- c(ce_fdh = 1, cr_fdh = 2)[ceil_names]
      legend("bottomright", legend = ceil_names, col = cols, lty = ltys,
             lwd = 2, cex = 0.7, bty = "n")
    }
  }
}

#' @noRd
plot_nca_effects <- function(nca_result, ...) {
  es <- nca_result$effect_sizes
  n_pred <- nrow(es)
  n_ceil <- ncol(es)

  if (n_pred == 0) {
    message("No predictors to plot.")
    return(invisible(NULL))
  }

  bar_width <- 0.8 / n_ceil
  positions <- seq_len(n_pred)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  cols <- if (n_ceil <= 2) {
    c("steelblue", "coral")[seq_len(n_ceil)]
  } else {
    palette.colors(n = n_ceil, palette = "Set1")
  }

  ylim_max <- max(0.5, max(es, na.rm = TRUE) * 1.2)

  par(mar = c(5, 4, 4, 2) + 0.1)
  plot(NULL, NULL,
       xlim = c(0.5, n_pred + 0.5),
       ylim = c(0, ylim_max),
       xaxt = "n", bty = "n",
       xlab = "", ylab = "Effect Size (d)",
       main = paste("NCA Effect Sizes:", nca_result$target),
       ...)

  abline(h = 0.1, col = "darkgray", lty = 2)
  text(n_pred + 0.4, 0.1, "d = 0.1", cex = 0.7, col = "darkgray", pos = 3)

  for (j in seq_len(n_ceil)) {
    offset <- (j - (n_ceil + 1) / 2) * bar_width
    for (i in seq_len(n_pred)) {
      val <- es[i, j]
      if (!is.na(val)) {
        rect(positions[i] + offset - bar_width / 2,
             0,
             positions[i] + offset + bar_width / 2,
             val,
             col = cols[j], border = NA)
      }
    }
  }

  axis(1, at = positions, labels = rownames(es), las = 2)

  if (n_ceil > 1) {
    legend("topright", legend = colnames(es), fill = cols,
           cex = 0.7, bty = "n")
  }
}

# =============================================================================
# NCA-ESSE: EFFECT SIZE SENSITIVITY EXTENSION
# =============================================================================
# Implements the NCA-ESSE method from:
#
# - Becker, J.-M., Richter, N. F., Ringle, C. M., & Sarstedt, M. (2026).
#   "Must-have, or maybe not? A sensitivity-based extension to necessary
#   condition analysis." Journal of Business Research, 206, 115920.
#
# NCA-ESSE assesses the sensitivity of NCA effect sizes to extreme response
# patterns by systematically varying the ECDF ceiling threshold. It compares
# empirical effect size changes against a theoretical benchmark (joint uniform
# distribution) to determine whether a necessary condition is robust.
# =============================================================================

#' Compute the joint empirical CDF for NCA: P(X <= x_i, Y >= y_i).
#'
#' For each observation (x_i, y_i), computes the proportion of observations
#' with X <= x_i AND Y >= y_i. Low values indicate extreme upper-left
#' combinations (low X, high Y) that populate the NCA ceiling zone.
#' @noRd
compute_ecdf_nca <- function(x, y) {
  n <- length(x)
  vapply(seq_len(n), function(i) {
    sum(x <= x[i] & y >= y[i])
  }, numeric(1)) / n
}

#' Compute theoretical CE-FDH benchmark effect size under joint uniform.
#'
#' Under a joint uniform distribution, the CE-FDH effect size at ECDF
#' threshold t is d = t(1 - ln(t)). Derived from the area of the region
#' \{(x,y) : x(1-y) <= t\} on the unit square (Becker et al., 2026).
#' @noRd
benchmark_effect_size <- function(t) {
  ifelse(t == 0, 0, t * (1 - log(t)))
}

#' NCA with Effect Size Sensitivity Extension (NCA-ESSE)
#'
#' \code{assess_nca_esse} applies the NCA-ESSE method (Becker et al., 2026)
#' to assess how sensitive NCA effect sizes are to extreme response patterns.
#' It systematically varies the ECDF ceiling threshold and compares empirical
#' effect size changes against a theoretical benchmark (joint uniform
#' distribution) to determine whether a necessary condition is robust.
#'
#' The method proceeds in three steps:
#' \enumerate{
#'   \item Run standard NCA (threshold 0\%) to identify potential necessary
#'     conditions
#'   \item Apply NCA-ESSE: vary the ECDF threshold from 0\% to a maximum
#'     (default 5\%), computing CE-FDH effect sizes at each level and
#'     comparing against the uniform benchmark
#'   \item (Optional) Select a threshold where the empirical effect size
#'     meaningfully exceeds the benchmark for further evaluation
#' }
#'
#' At each threshold \code{t}, observations whose joint ECDF_NCA value
#' (P(X <= x, Y >= y)) is at most \code{t} are treated as extreme and
#' excluded before computing the CE-FDH ceiling. The uniform benchmark
#' d = t(1 - ln(t)) gives the expected effect size if no necessity exists.
#'
#' CE-FDH and CR-FDH ceilings are computed internally. For other ceiling
#' techniques, the \pkg{NCA} package must be installed.
#'
#' @param seminr_model An estimated SEMinR model from \code{estimate_pls()}.
#' @param target Name of the endogenous (outcome) construct.
#' @param predictors Optional character vector of predictor construct names.
#'   If \code{NULL} (default), auto-detected from the structural model.
#' @param thresholds Numeric vector of ECDF thresholds to evaluate
#'   (default \code{seq(0, 0.05, by = 0.005)}). Values must be in [0, 1].
#' @param ceiling Ceiling technique (default \code{"ce_fdh"}). The analytical
#'   benchmark is derived for CE-FDH; other techniques will trigger a warning.
#' @param test.rep Number of permutation test repetitions at each threshold
#'   (default 0). Set > 0 to obtain p-values, but note this multiplies
#'   computation time by the number of thresholds.
#' @param steps Number of steps in the bottleneck table (default 10).
#' @param seed Random seed for reproducibility (default 123).
#' @param ... Additional arguments passed to \code{NCA::nca_analysis()} when
#'   using ceiling techniques not supported internally.
#'
#' @return An object of class \code{nca_esse} containing:
#'   \item{effect_sizes}{Matrix of empirical effect sizes (thresholds x predictors)}
#'   \item{benchmark}{Matrix of theoretical benchmark effect sizes}
#'   \item{delta}{Matrix of empirical minus benchmark differences}
#'   \item{significance}{Matrix of p-values (NULL if test.rep = 0)}
#'   \item{pls_model}{The original estimated seminr model}
#'   \item{target}{Name of the target construct}
#'   \item{predictors}{Character vector of predictor names}
#'   \item{thresholds}{Numeric vector of ECDF thresholds used}
#'   \item{ceiling}{Ceiling technique used}
#'   \item{n_obs}{Number of observations}
#'
#' @seealso \code{\link{assess_nca}} for standard NCA analysis
#'
#' @references
#' Becker, J.-M., Richter, N. F., Ringle, C. M., & Sarstedt, M. (2026).
#' Must-have, or maybe not? A sensitivity-based extension to necessary
#' condition analysis. Journal of Business Research, 206, 115920.
#'
#' @examples
#' library(seminr)
#' library(seminrExtras)
#'
#' mobi_mm <- constructs(
#'   composite("Image",        multi_items("IMAG", 1:5)),
#'   composite("Value",        multi_items("PERV", 1:2)),
#'   composite("Satisfaction", multi_items("CUSA", 1:3))
#' )
#'
#' mobi_sm <- relationships(
#'   paths(from = c("Image", "Value"), to = "Satisfaction")
#' )
#'
#' mobi_pls <- estimate_pls(data = mobi,
#'                           measurement_model = mobi_mm,
#'                           structural_model  = mobi_sm)
#'
#' \donttest{
#' esse_result <- assess_nca_esse(mobi_pls,
#'                                 target = "Satisfaction",
#'                                 seed = 123)
#' print(esse_result)
#' plot(esse_result, type = "sensitivity")
#' }
#'
#' @export
assess_nca_esse <- function(seminr_model,
                             target,
                             predictors = NULL,
                             thresholds = seq(0, 0.05, by = 0.005),
                             ceiling = "ce_fdh",
                             test.rep = 0,
                             steps = 10,
                             seed = 123,
                             ...) {

  # ---------------------------------------------------------------------------
  # Step 1: Validate inputs
  # ---------------------------------------------------------------------------
  if (!validate_seminr_model(seminr_model, "assess_nca_esse")) return(NULL)

  if (any(thresholds < 0) || any(thresholds > 1)) {
    stop("thresholds must be between 0 and 1.", call. = FALSE)
  }
  predictors <- validate_nca_inputs(seminr_model, target, predictors, test.rep)

  if (ceiling != "ce_fdh") {
    warning("NCA-ESSE benchmark is derived for CE-FDH (Becker et al., 2026). ",
            "Benchmark may not be directly comparable with '", ceiling, "'.",
            call. = FALSE)
  }

  use_internal <- ceiling %in% INTERNAL_CEILINGS
  if (!use_internal) {
    check_nca_installed()
  }

  scores <- as.data.frame(seminr_model$construct_scores)
  n_obs <- nrow(scores)

  set.seed(seed)

  # ---------------------------------------------------------------------------
  # Step 2: Initialize result matrices
  # ---------------------------------------------------------------------------
  threshold_labels <- paste0(thresholds * 100, "%")
  empirical <- matrix(NA_real_, nrow = length(thresholds), ncol = length(predictors),
                       dimnames = list(threshold_labels, predictors))
  significance <- if (test.rep > 0) {
    matrix(NA_real_, nrow = length(thresholds), ncol = length(predictors),
           dimnames = list(threshold_labels, predictors))
  } else NULL

  # ---------------------------------------------------------------------------
  # Step 3: Compute empirical effect sizes at each ECDF threshold
  # ---------------------------------------------------------------------------
  for (p_idx in seq_along(predictors)) {
    pred <- predictors[p_idx]
    x <- scores[[pred]]
    y <- scores[[target]]

    ecdf_nca <- compute_ecdf_nca(x, y)

    for (t_idx in seq_along(thresholds)) {
      t_val <- thresholds[t_idx]

      # At threshold 0, keep all observations (standard NCA)
      if (t_val == 0) {
        filtered_scores <- scores
      } else {
        keep <- ecdf_nca > t_val
        filtered_scores <- scores[keep, , drop = FALSE]
      }

      if (nrow(filtered_scores) < 10) {
        warning("Fewer than 10 observations at threshold ", threshold_labels[t_idx],
                " for ", pred, "; skipping.", call. = FALSE)
        next
      }

      filtered_x <- filtered_scores[[pred]]
      filtered_y <- filtered_scores[[target]]

      if (use_internal) {
        d <- nca_effect_size(filtered_x, filtered_y, ceiling)
        empirical[t_idx, p_idx] <- d

        if (test.rep > 0) {
          significance[t_idx, p_idx] <- nca_permutation_test(
            filtered_x, filtered_y, ceiling, d, test.rep
          )
        }
      } else {
        nca_res <- suppressMessages(
          NCA::nca_analysis(
            data = filtered_scores, x = pred, y = target,
            ceilings = ceiling, test.rep = test.rep,
            steps = steps, ...
          )
        )
        params <- nca_res$summaries[[pred]]$params
        if (!is.null(params) && "Effect size" %in% rownames(params) &&
            ceiling %in% colnames(params)) {
          empirical[t_idx, p_idx] <- params["Effect size", ceiling]
        }
        if (!is.null(significance) && !is.null(params) &&
            "p-value" %in% rownames(params) && ceiling %in% colnames(params)) {
          significance[t_idx, p_idx] <- params["p-value", ceiling]
        }
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Step 4: Compute theoretical benchmark (joint uniform distribution)
  # ---------------------------------------------------------------------------
  benchmark <- matrix(benchmark_effect_size(thresholds),
                       nrow = length(thresholds), ncol = length(predictors),
                       dimnames = list(threshold_labels, predictors))

  # ---------------------------------------------------------------------------
  # Step 5: Compute sensitivity (empirical - benchmark)
  # ---------------------------------------------------------------------------
  delta <- empirical - benchmark

  result <- list(
    effect_sizes = empirical,
    benchmark = benchmark,
    delta = delta,
    significance = significance,
    pls_model = seminr_model,
    target = target,
    predictors = predictors,
    thresholds = thresholds,
    ceiling = ceiling,
    n_obs = n_obs
  )

  class(result) <- c("nca_esse", class(result))
  result
}

# =============================================================================
# NCA-ESSE S3 METHODS
# =============================================================================

#' @export
print.nca_esse <- function(x, ...) {
  cat("NCA-ESSE: Effect Size Sensitivity Extension\n")
  cat("=============================================\n")
  cat("Target:", x$target, "\n")
  cat("Predictors:", paste(x$predictors, collapse = ", "), "\n")
  cat("Ceiling:", x$ceiling, "\n")
  cat("Observations:", x$n_obs, "\n")
  cat("Thresholds:", paste0(min(x$thresholds) * 100, "% to ",
      max(x$thresholds) * 100, "%"), "\n\n")

  cat("Empirical effect sizes by ECDF threshold:\n")
  print(round(x$effect_sizes, 4), ...)
  cat("\n")

  cat("Benchmark (uniform) effect sizes:\n")
  bench_vec <- x$benchmark[, 1]
  names(bench_vec) <- rownames(x$benchmark)
  print(round(bench_vec, 4), ...)
  cat("\n")

  cat("Sensitivity (empirical - benchmark):\n")
  print(round(x$delta, 4), ...)

  if (!is.null(x$significance)) {
    cat("\nPermutation p-values:\n")
    print(round(x$significance, 4), ...)
  }

  invisible(x)
}

#' @export
summary.nca_esse <- function(object, ...) {
  tables <- lapply(object$predictors, function(pred) {
    emp <- object$effect_sizes[, pred]
    bench <- object$benchmark[, pred]
    df <- data.frame(
      ECDF_threshold = object$thresholds,
      Empirical_d = emp,
      Benchmark_d = bench,
      Difference = emp - bench,
      row.names = NULL
    )
    if (!is.null(object$significance)) {
      df$p_value <- object$significance[, pred]
    }
    df
  })
  names(tables) <- object$predictors

  result <- list(
    target = object$target,
    predictors = object$predictors,
    ceiling = object$ceiling,
    n_obs = object$n_obs,
    tables = tables
  )

  class(result) <- c("summary.nca_esse", class(result))
  result
}

#' @export
print.summary.nca_esse <- function(x, ...) {
  cat("NCA-ESSE Summary (Becker et al., 2026)\n")
  cat("=======================================\n")
  cat("Target:", x$target, "\n")
  cat("Ceiling:", x$ceiling, "\n")
  cat("Observations:", x$n_obs, "\n\n")

  for (pred in x$predictors) {
    cat("Predictor:", pred, "\n")
    cat(strrep("-", 60), "\n")
    print(round(x$tables[[pred]], 4), row.names = FALSE, ...)
    cat("\n")
  }

  invisible(x)
}

#' Plot NCA-ESSE Results
#'
#' @param x An \code{nca_esse} object from \code{assess_nca_esse()}.
#' @param type One of \code{"sensitivity"} (effect size vs threshold, Fig. 4
#'   in Becker et al.) or \code{"difference"} (incremental empirical minus
#'   benchmark difference, Fig. 6 in Becker et al.).
#' @param ... Additional arguments passed to the underlying plot function.
#'
#' @export
plot.nca_esse <- function(x, type = c("sensitivity", "difference"), ...) {
  type <- match.arg(type)

  switch(type,
    sensitivity = plot_esse_sensitivity(x, ...),
    difference  = plot_esse_difference(x, ...)
  )
}

#' Sensitivity plot (Fig. 4): empirical and benchmark effect sizes vs threshold.
#' @noRd
plot_esse_sensitivity <- function(esse, ...) {
  n_pred <- length(esse$predictors)
  thresholds <- esse$thresholds

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  if (n_pred > 1) {
    ncol_plot <- min(n_pred, 3)
    nrow_plot <- ceiling(n_pred / ncol_plot)
    par(mfrow = c(nrow_plot, ncol_plot))
  }

  for (pred in esse$predictors) {
    emp <- esse$effect_sizes[, pred]
    bench <- esse$benchmark[, pred]
    ylim <- c(0, max(c(emp, bench, 0.1), na.rm = TRUE) * 1.2)

    plot(thresholds, emp, type = "b", pch = 19,
         ylim = ylim, xlab = "ECDF threshold",
         ylab = "NCA effect size (d)",
         main = paste("NCA-ESSE:", pred, "->", esse$target),
         bty = "n", ...)
    lines(thresholds, bench, type = "b", pch = 17, lty = 2, col = "gray50")
    abline(h = 0.1, col = "darkgray", lty = 3)
    legend("topleft", legend = c("Empirical", "Benchmark (uniform)"),
           pch = c(19, 17), lty = c(1, 2), col = c("black", "gray50"),
           cex = 0.7, bty = "n")
  }
}

#' Difference plot (Fig. 6): incremental empirical - benchmark difference.
#' @noRd
plot_esse_difference <- function(esse, ...) {
  n_pred <- length(esse$predictors)
  thresholds <- esse$thresholds[-1]

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  if (n_pred > 1) {
    ncol_plot <- min(n_pred, 3)
    nrow_plot <- ceiling(n_pred / ncol_plot)
    par(mfrow = c(nrow_plot, ncol_plot))
  }

  for (pred in esse$predictors) {
    emp <- esse$effect_sizes[, pred]
    bench <- esse$benchmark[, pred]

    delta_emp <- diff(emp)
    delta_bench <- diff(bench)
    delta_diff <- delta_emp - delta_bench

    ylim <- range(c(delta_diff, 0), na.rm = TRUE)
    ylim <- ylim + c(-1, 1) * max(0.05, diff(ylim) * 0.1)

    plot(thresholds, delta_diff, type = "b", pch = 19,
         ylim = ylim, xlab = "ECDF threshold",
         ylab = expression(Delta ~ "Empirical" - Delta ~ "Benchmark"),
         main = paste("ESSE Difference:", pred, "->", esse$target),
         bty = "n", ...)
    abline(h = 0, col = "darkgray", lty = 2)
  }
}
