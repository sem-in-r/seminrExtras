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

#' @noRd
get_direct_predictors <- function(model, target) {
  sm <- model$smMatrix
  unname(sm[sm[, "target"] == target, "source"])
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

  if (!validate_seminr_model(seminr_model, "assess_nca")) {
    return(NULL)
  }

  # Validate target
  construct_names <- colnames(seminr_model$construct_scores)
  if (!(target %in% construct_names)) {
    stop("target '", target, "' not found in model constructs: ",
         paste(construct_names, collapse = ", "),
         call. = FALSE)
  }

  # Auto-detect or validate predictors
  if (is.null(predictors)) {
    predictors <- get_direct_predictors(seminr_model, target)
    if (length(predictors) == 0) {
      stop("No direct predictors found for target '", target,
           "' in the structural model.",
           call. = FALSE)
    }
  } else {
    invalid <- setdiff(predictors, construct_names)
    if (length(invalid) > 0) {
      stop("Predictor(s) not found in model constructs: ",
           paste(invalid, collapse = ", "),
           call. = FALSE)
    }
  }

  # Build data frame of construct scores
  scores <- as.data.frame(seminr_model$construct_scores)

  # Run NCA analysis
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

  # Extract effect sizes and p-values
  effect_sizes <- format_nca_effects(nca_raw, predictors, ceilings)
  significance <- format_nca_significance(nca_raw, predictors, ceilings)
  bottleneck <- format_nca_bottleneck(nca_raw, ceilings)

  result <- list(
    nca_raw      = nca_raw,
    effect_sizes = effect_sizes,
    significance = significance,
    bottleneck   = bottleneck,
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

  # Identify necessary conditions
  necessary <- identify_necessary(x$effect_sizes, x$significance)
  necessary_preds <- rownames(x$effect_sizes)[apply(necessary, 1, any)]

  if (length(necessary_preds) > 0) {
    cat("Necessary conditions (d >= 0.1, p < 0.05):",
        paste(necessary_preds, collapse = ", "), "\n")
  } else {
    cat("No necessary conditions identified (d >= 0.1, p < 0.05)\n")
  }

  invisible(x)
}

#' @export
summary.nca_analysis <- function(object, ...) {
  necessary <- identify_necessary(object$effect_sizes, object$significance)
  necessary_preds <- rownames(object$effect_sizes)[apply(necessary, 1, any)]

  result <- list(
    target       = object$target,
    predictors   = object$predictors,
    ceilings     = object$ceilings,
    n_obs        = nrow(object$pls_model$construct_scores),
    effect_sizes = object$effect_sizes,
    significance = object$significance,
    bottleneck   = object$bottleneck,
    necessary    = necessary,
    necessary_predictors = necessary_preds
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
    grDevices::palette.colors(n = n_ceil, palette = "Set1")
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
