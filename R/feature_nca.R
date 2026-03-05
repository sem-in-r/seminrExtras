# =============================================================================
# feature_nca.R - Necessary Condition Analysis (NCA) for PLS-SEM
# =============================================================================
# This file wraps the NCA R package (Dul, 2016, 2020) for seamless use with
# seminr estimated models. NCA tests whether predictors are necessary
# conditions for an outcome, complementing PLS-SEM's sufficiency logic.
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
# DEPENDENCY CHECK
# =============================================================================

#' @noRd
check_nca_installed <- function() {
  if (!requireNamespace("NCA", quietly = TRUE)) {
    stop("Package 'NCA' is required for assess_nca(). ",
         "Install it with: install.packages('NCA')",
         call. = FALSE)
  }
}

# =============================================================================
# PREDICTOR EXTRACTION
# =============================================================================

#' Extract direct predictor construct names for a given target from the
#' structural model matrix. unname() is required because subsetting smMatrix
#' retains row names, which would cause NCA::nca_analysis() to fail on
#' predictor name matching.
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
#' direct predictors of the target from the structural model, and delegates
#' to \code{NCA::nca_analysis()}.
#'
#' @param seminr_model An estimated SEMinR model from \code{estimate_pls()}.
#' @param target Name of the endogenous (outcome) construct (character).
#' @param predictors Optional character vector of predictor construct names.
#'   If \code{NULL} (default), auto-detected from the structural model as
#'   direct predictors of \code{target}.
#' @param ceilings Character vector of ceiling techniques to use
#'   (default \code{c("ce_fdh", "cr_fdh")}). See \code{NCA::nca_analysis()}
#'   for all available techniques.
#' @param test.rep Number of permutation test repetitions for significance
#'   testing (default 1000). Set to 0 to skip significance testing.
#' @param steps Number of steps in the bottleneck table (default 10).
#' @param seed Random seed for reproducibility (default 123).
#' @param ... Additional arguments passed to \code{NCA::nca_analysis()}.
#'
#' @return An object of class \code{nca_analysis} containing:
#'   \item{nca_raw}{The raw result from \code{NCA::nca_analysis()}}
#'   \item{effect_sizes}{Matrix of effect sizes (d) per predictor and ceiling}
#'   \item{significance}{Matrix of p-values per predictor and ceiling (NULL if test.rep = 0)}
#'   \item{bottleneck}{List of bottleneck tables, one per ceiling technique}
#'   \item{pls_model}{The original estimated seminr model}
#'   \item{target}{Name of the target construct}
#'   \item{predictors}{Character vector of predictor names used}
#'   \item{ceilings}{Character vector of ceiling techniques used}
#'
#' @seealso \code{\link[NCA]{nca_analysis}} for the underlying NCA implementation
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

  check_nca_installed()

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

  # Delegate to NCA package for ceiling envelope fitting and permutation tests
  set.seed(seed)
  nca_raw <- suppressMessages(
    NCA::nca_analysis(
      data = scores,
      x = predictors,
      y = target,
      ceilings = ceilings,
      test.rep = test.rep,
      steps = steps,
      ...
    )
  )

  # ---------------------------------------------------------------------------
  # Step 3: Extract and format results
  # ---------------------------------------------------------------------------
  effect_sizes <- format_nca_effects(nca_raw, predictors, ceilings)
  significance <- format_nca_significance(nca_raw, predictors, ceilings)
  bottleneck <- format_nca_bottleneck(nca_raw, ceilings)

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
# EXTRACTION HELPERS
# =============================================================================

#' Extract effect sizes from NCA result into a predictor x ceiling matrix.
#'
#' NCA stores results per-predictor in nca_raw$summaries[[pred]]$params,
#' a matrix with rows like "Scope", "Effect size", "p-value" and columns
#' for each ceiling technique. We extract the "Effect size" row.
#' @noRd
format_nca_effects <- function(nca_raw, predictors, ceilings) {
  mat <- matrix(NA_real_,
                nrow = length(predictors),
                ncol = length(ceilings),
                dimnames = list(predictors, ceilings))

  for (pred in predictors) {
    if (!is.null(nca_raw$summaries[[pred]]$params)) {
      params <- nca_raw$summaries[[pred]]$params
      for (ceil in ceilings) {
        if (ceil %in% colnames(params) && "Effect size" %in% rownames(params)) {
          mat[pred, ceil] <- params["Effect size", ceil]
        }
      }
    }
  }

  comment(mat) <- "NCA Effect Sizes (d >= 0.1 indicates a necessary condition)"
  class(mat) <- c("table_output", class(mat))
  mat
}

#' Extract p-values from NCA result, same structure as format_nca_effects().
#' Returns NULL if no significance testing was performed (test.rep = 0).
#' @noRd
format_nca_significance <- function(nca_raw, predictors, ceilings) {
  # Check if significance testing was performed
  has_pvals <- any(vapply(predictors, function(pred) {
    params <- nca_raw$summaries[[pred]]$params
    !is.null(params) && "p-value" %in% rownames(params) &&
      any(!is.na(params["p-value", ]))
  }, logical(1)))

  if (!has_pvals) return(NULL)

  mat <- matrix(NA_real_,
                nrow = length(predictors),
                ncol = length(ceilings),
                dimnames = list(predictors, ceilings))

  for (pred in predictors) {
    if (!is.null(nca_raw$summaries[[pred]]$params)) {
      params <- nca_raw$summaries[[pred]]$params
      for (ceil in ceilings) {
        if (ceil %in% colnames(params) && "p-value" %in% rownames(params)) {
          mat[pred, ceil] <- params["p-value", ceil]
        }
      }
    }
  }

  comment(mat) <- "NCA Significance (permutation test p-values)"
  class(mat) <- c("table_output", class(mat))
  mat
}

#' Extract bottleneck tables (one per ceiling) from NCA result.
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
#' @param type One of \code{"scatter"} (ceiling line scatter plots via NCA
#'   package) or \code{"effects"} (bar plot of effect sizes).
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

#' @noRd
plot_nca_scatter <- function(nca_result, ...) {
  check_nca_installed()
  plot(nca_result$nca_raw, ...)
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

  # Bar positions
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

  # Threshold line
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
#' @param ... Additional arguments passed to \code{NCA::nca_analysis()}.
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

  check_nca_installed()

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
  # For each predictor, compute ECDF_NCA(x_i, y_i) = P(X <= x_i, Y >= y_i),
  # then at each threshold t, remove observations with ECDF_NCA <= t (extreme
  # upper-left cases) and run standard NCA on the remaining data.
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

      nca_res <- suppressMessages(
        NCA::nca_analysis(
          data = filtered_scores,
          x = pred,
          y = target,
          ceilings = ceiling,
          test.rep = test.rep,
          steps = steps,
          ...
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

  # ---------------------------------------------------------------------------
  # Step 4: Compute theoretical benchmark (joint uniform distribution)
  # ---------------------------------------------------------------------------
  # Under joint uniform, CE-FDH effect size at threshold t is d = t(1 - ln(t))
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
  # Build Table A2-style output per predictor
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
  # Skip the first threshold (0%) since diff requires pairs
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

    # Incremental changes between consecutive thresholds
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
