# =============================================================================
# feature_pcm.R - Predictive Contribution of the Mediator (PCM)
# =============================================================================
# Evaluates the predictive contribution of mediating constructs by comparing
# predictions from Direct Antecedents (DA) vs Earliest Antecedents (EA)
# approaches on isolated sub-models.
#
# References:
# - Danks, N. P. (2021). The Piggy in the Middle: The Role of Mediators in
#   PLS-SEM Prediction. The DATA BASE for Advances in Information Systems,
#   52(SI), 24-42.
# - Shmueli, G., Sarstedt, M., Hair, J. F., Cheah, J.-H., Ting, H.,
#   Vaithilingam, S. & Ringle, C. M. (2019). Predictive Model Assessment in
#   PLS-SEM: Guidelines for Using PLSpredict. European Journal of Marketing,
#   53(11), 2322-2347.
# =============================================================================

# =============================================================================
# INTERNAL HELPERS
# =============================================================================

#' Detect the single final endogenous construct in a model.
#'
#' A "final" endogenous construct is one that appears as a target in the
#' structural model but never as a source.
#' @noRd
detect_final_endogenous <- function(model) {
  sm <- model$smMatrix
  endogenous <- unique(unname(sm[, "target"]))
  sources    <- unique(unname(sm[, "source"]))
  final_endo <- setdiff(endogenous, sources)

  if (length(final_endo) == 1) return(final_endo)
  if (length(final_endo) == 0) {
    stop("No final endogenous construct found (all endogenous constructs ",
         "are also sources). Please specify 'target'.", call. = FALSE)
  }
  stop("Multiple final endogenous constructs found: ",
       paste(final_endo, collapse = ", "),
       ". Please specify 'target'.", call. = FALSE)
}

#' Find all simple mediation paths (X -> M -> Y) to a target.
#'
#' A mediation path exists when a direct predictor of the target (M) is itself
#' endogenous (has its own antecedent X). Interaction constructs and
#' higher-order constructs are excluded (HOC items are LOC construct names,
#' not raw indicators, so sub-model reconstruction is not supported).
#' @noRd
find_mediation_paths <- function(model, target) {
  sm <- model$smMatrix

  # Direct predictors of target, excluding interaction constructs
  direct_preds <- unname(sm[sm[, "target"] == target, "source"])
  direct_preds <- direct_preds[!grepl("\\*", direct_preds)]

  paths <- list()
  for (mediator in direct_preds) {
    antecedents <- unname(sm[sm[, "target"] == mediator, "source"])
    antecedents <- antecedents[!grepl("\\*", antecedents)]

    for (antecedent in antecedents) {
      # Skip paths involving HOC constructs (cannot build sub-model)
      triple <- c(antecedent, mediator, target)
      hoc_in_triple <- sapply(triple, is_hoc_construct, model = model)
      if (any(hoc_in_triple)) {
        hoc_names <- triple[hoc_in_triple]
        warning(sprintf("Skipping %s -> %s -> %s: higher-order construct(s) %s ",
                        antecedent, mediator, target,
                        paste(hoc_names, collapse = ", ")),
                "cannot be used in PCM sub-models.", call. = FALSE)
        next
      }

      paths[[length(paths) + 1]] <- list(
        antecedent = antecedent,
        mediator   = mediator,
        target     = target
      )
    }
  }
  paths
}

#' Build an isolated sub-model for a single mediation triple.
#'
#' Creates and estimates a partial mediation model (X -> M, M -> Y, X -> Y)
#' containing only the three specified constructs.
#' @noRd
build_isolated_sub_model <- function(model, antecedent, mediator, target) {
  constructs_to_keep <- c(antecedent, mediator, target)

  # Reconstruct measurement model from the parent model's mmMatrix
  meas_specs <- list()
  for (cname in constructs_to_keep) {
    indicators <- items_of_construct(cname, model)
    mode <- get_construct_mode(cname, model)

    if (mode == "C") {
      meas_specs[[length(meas_specs) + 1]] <-
        seminr::reflective(cname, indicators)
    } else if (mode == "B") {
      meas_specs[[length(meas_specs) + 1]] <-
        seminr::composite(cname, indicators, weights = seminr::mode_B)
    } else {
      meas_specs[[length(meas_specs) + 1]] <-
        seminr::composite(cname, indicators, weights = seminr::mode_A)
    }
  }
  measurement_model <- do.call(seminr::constructs, meas_specs)

  # Partial mediation structural model: X -> M, M -> Y, X -> Y
  structural_model <- seminr::relationships(
    seminr::paths(from = antecedent, to = mediator),
    seminr::paths(from = c(antecedent, mediator), to = target)
  )

  # Re-estimate on original data
  suppressMessages(
    estimate_pls(
      data = model$rawdata,
      measurement_model = measurement_model,
      structural_model = structural_model
    )
  )
}

#' Compute PCM for a single mediation path.
#'
#' Builds an isolated sub-model, runs DA and EA predictions, and computes the
#' PCM metric per indicator of the target construct.
#' @noRd
compute_pcm_for_path <- function(model, antecedent, mediator, target,
                                 noFolds, reps) {
  # Build isolated sub-model
  sub_model <- build_isolated_sub_model(model, antecedent, mediator, target)

  # Run DA and EA cross-validated predictions
  pred_da <- suppressMessages(
    predict_pls(sub_model, technique = predict_DA,
                noFolds = noFolds, reps = reps)
  )
  pred_ea <- suppressMessages(
    predict_pls(sub_model, technique = predict_EA,
                noFolds = noFolds, reps = reps)
  )

  # Get summary metrics
  sum_da <- summary(pred_da)
  sum_ea <- summary(pred_ea)

  # Extract target indicators
  target_indicators <- items_of_construct(target, sub_model)

  # Get RMSE and MAE for target indicators (metrics = rows, indicators = cols)
  rmse_da <- sum_da$PLS_out_of_sample["RMSE", target_indicators]
  mae_da  <- sum_da$PLS_out_of_sample["MAE",  target_indicators]
  rmse_ea <- sum_ea$PLS_out_of_sample["RMSE", target_indicators]
  mae_ea  <- sum_ea$PLS_out_of_sample["MAE",  target_indicators]

  # PCM = (METRIC_EA - METRIC_DA) / METRIC_EA  (equation 2)
  pcm_rmse <- (rmse_ea - rmse_da) / rmse_ea
  pcm_mae  <- (mae_ea  - mae_da)  / mae_ea

  # Build results matrix
  result_matrix <- cbind(
    RMSE_DA  = rmse_da,
    RMSE_EA  = rmse_ea,
    PCM_RMSE = pcm_rmse,
    MAE_DA   = mae_da,
    MAE_EA   = mae_ea,
    PCM_MAE  = pcm_mae
  )
  rownames(result_matrix) <- target_indicators

  list(
    antecedent = antecedent,
    mediator   = mediator,
    target     = target,
    results    = result_matrix,
    pcm_rmse   = as.vector(pcm_rmse),
    pcm_mae    = as.vector(pcm_mae)
  )
}

#' Check whether a construct is a higher-order construct (HOC).
#'
#' HOC items are LOC construct names (not raw indicators), identifiable by
#' being column names in the outer_weights matrix.
#' @noRd
is_hoc_construct <- function(construct, model) {
  items <- items_of_construct(construct, model)
  construct_names <- colnames(model$outer_weights)
  any(items %in% construct_names)
}

#' Classify a PCM value using rules of thumb from Danks (2021).
#' @noRd
classify_pcm <- function(pcm_value) {
  if (is.na(pcm_value)) return("NA")
  if (pcm_value < 0)    return("Negative")
  if (pcm_value < 0.05) return("Weak")
  if (pcm_value < 0.10) return("Moderate")
  return("Strong")
}

# =============================================================================
# EXPORTED FUNCTIONS
# =============================================================================

#' Predictive Contribution of the Mediator (PCM)
#'
#' Evaluates the predictive contribution of mediating constructs in a PLS-SEM
#' model. For each mediation path (antecedent -> mediator -> target), computes
#' the PCM metric by comparing cross-validated predictions from the Direct
#' Antecedents (DA) and Earliest Antecedents (EA) approaches on an isolated
#' partial mediation sub-model.
#'
#' @param seminr_model A PLS model estimated by \code{seminr::estimate_pls()}.
#' @param target Character string specifying the outcome construct. If
#'   \code{NULL}, auto-detects when a single final endogenous construct exists.
#' @param noFolds Integer, number of folds for cross-validation (default 10).
#' @param reps Integer, number of cross-validation repetitions (default 10).
#'
#' @return An object of class \code{pcm_analysis} containing:
#' \describe{
#'   \item{target}{The outcome construct name.}
#'   \item{mediation_paths}{List of identified mediation triples
#'     (antecedent, mediator, target).}
#'   \item{pcm_results}{List of per-path results, each with RMSE/MAE for DA
#'     and EA approaches plus PCM values per indicator.}
#'   \item{noFolds}{Number of cross-validation folds used.}
#'   \item{reps}{Number of cross-validation repetitions used.}
#' }
#'
#' @details
#' The PCM metric (Danks, 2021) quantifies the predictive improvement due to
#' the mediator. It is computed by isolating each mediation path into a partial
#' mediation sub-model (X -> M -> Y with direct path X -> Y), then comparing
#' out-of-sample predictions from two approaches:
#'
#' \itemize{
#'   \item \strong{DA (Direct Antecedents)}: Predicts Y using its direct
#'     structural predictors (M and X).
#'   \item \strong{EA (Earliest Antecedents)}: Predicts Y using only the
#'     earliest antecedent (X), propagated through the structural model.
#' }
#'
#' \deqn{PCM = \frac{METRIC_{EA} - METRIC_{DA}}{METRIC_{EA}}}
#'
#' Rules of thumb for PCM interpretation:
#' \itemize{
#'   \item \strong{0 to 0.05}: Weak predictive contribution
#'   \item \strong{0.05 to 0.10}: Moderate predictive contribution
#'   \item \strong{Greater than 0.10}: Strong predictive contribution
#'   \item \strong{Negative}: Mediator damages predictive accuracy
#' }
#'
#' @references Danks, N. P. (2021). The Piggy in the Middle: The Role of
#' Mediators in PLS-SEM Prediction. \emph{The DATA BASE for Advances in
#' Information Systems}, 52(SI), 24-42.
#'
#' @seealso \code{\link[seminr]{predict_pls}} for the underlying prediction
#'   methods.
#'
#' @examples
#' \donttest{
#' library(seminr)
#' data(mobi)
#'
#' # Simple mediation: Image -> Satisfaction -> Loyalty
#' mobi_mm <- constructs(
#'   composite("Image",        multi_items("IMAG", 1:5)),
#'   composite("Satisfaction", multi_items("CUSA", 1:3)),
#'   composite("Loyalty",      multi_items("CUSL", 1:3))
#' )
#' mobi_sm <- relationships(
#'   paths(from = "Image", to = "Satisfaction"),
#'   paths(from = "Satisfaction", to = "Loyalty"),
#'   paths(from = "Image", to = "Loyalty")
#' )
#' model <- estimate_pls(mobi, mobi_mm, mobi_sm)
#'
#' # Compute PCM
#' pcm <- assess_pcm(model, target = "Loyalty", noFolds = 10, reps = 10)
#' pcm
#' summary(pcm)
#' }
#'
#' @export
assess_pcm <- function(seminr_model, target = NULL, noFolds = 10, reps = 10) {
  # --- Validation ---
  if (!validate_seminr_model(seminr_model, "assess_pcm")) return(NULL)
  if (!validate_for_prediction(seminr_model, "assess_pcm")) return(NULL)

  if (!is.numeric(noFolds) || noFolds < 2) {
    stop("noFolds must be an integer >= 2.", call. = FALSE)
  }
  if (!is.numeric(reps) || reps < 1) {
    stop("reps must be an integer >= 1.", call. = FALSE)
  }

  # --- Identify target ---
  if (is.null(target)) {
    target <- detect_final_endogenous(seminr_model)
  }
  if (!(target %in% colnames(seminr_model$path_coef))) {
    stop("Target '", target, "' not found in model constructs.", call. = FALSE)
  }

  # --- Find mediation paths ---
  mediation_paths <- find_mediation_paths(seminr_model, target)
  if (length(mediation_paths) == 0) {
    stop("No mediation paths found for target '", target,
         "'. The target must have at least one predictor that is itself ",
         "endogenous.", call. = FALSE)
  }

  # --- Compute PCM per path ---
  pcm_results <- vector("list", length(mediation_paths))
  for (i in seq_along(mediation_paths)) {
    path <- mediation_paths[[i]]
    message(sprintf("Computing PCM for %s -> %s -> %s ...",
                    path$antecedent, path$mediator, path$target))
    pcm_results[[i]] <- compute_pcm_for_path(
      seminr_model, path$antecedent, path$mediator, path$target,
      noFolds, reps
    )
  }

  # --- Return ---
  output <- list(
    target          = target,
    mediation_paths = mediation_paths,
    pcm_results     = pcm_results,
    noFolds         = noFolds,
    reps            = reps
  )
  class(output) <- c("pcm_analysis", class(output))
  output
}

# =============================================================================
# S3 METHODS
# =============================================================================

#' @export
print.pcm_analysis <- function(x, ...) {
  cat("Predictive Contribution of the Mediator (PCM)\n")
  cat("==============================================\n")
  cat("Target:", x$target, "\n")
  cat("Cross-validation:", x$noFolds, "folds,", x$reps, "reps\n")
  cat("Mediation paths:", length(x$mediation_paths), "\n\n")

  for (i in seq_along(x$pcm_results)) {
    res  <- x$pcm_results[[i]]
    path <- x$mediation_paths[[i]]
    cat(sprintf("  %s -> %s -> %s\n",
                path$antecedent, path$mediator, path$target))

    avg_pcm <- mean(res$pcm_rmse)
    cat(sprintf("    Avg PCM (RMSE): %7.4f  [%s]\n",
                avg_pcm, classify_pcm(avg_pcm)))
    avg_pcm_mae <- mean(res$pcm_mae)
    cat(sprintf("    Avg PCM (MAE):  %7.4f  [%s]\n\n",
                avg_pcm_mae, classify_pcm(avg_pcm_mae)))
  }
  invisible(x)
}

#' @export
summary.pcm_analysis <- function(object, ...) {
  out <- list(
    target          = object$target,
    noFolds         = object$noFolds,
    reps            = object$reps,
    mediation_paths = object$mediation_paths,
    pcm_results     = object$pcm_results
  )
  class(out) <- c("summary.pcm_analysis", class(out))
  out
}

#' @export
print.summary.pcm_analysis <- function(x, ...) {
  cat("Predictive Contribution of the Mediator (PCM)\n")
  cat("==============================================\n")
  cat("Target:", x$target, "\n")
  cat("Cross-validation:", x$noFolds, "folds,", x$reps, "reps\n\n")

  for (i in seq_along(x$pcm_results)) {
    res  <- x$pcm_results[[i]]
    path <- x$mediation_paths[[i]]

    cat(sprintf("Mediation: %s -> %s -> %s\n",
                path$antecedent, path$mediator, path$target))
    cat(paste(rep("-", 60), collapse = ""), "\n")

    # Header
    cat(sprintf("  %-12s %8s %8s %9s %8s %8s %9s  %s\n",
                "Indicator", "RMSE_DA", "RMSE_EA", "PCM_RMSE",
                "MAE_DA", "MAE_EA", "PCM_MAE", "Conclusion"))

    display <- res$results
    for (j in seq_len(nrow(display))) {
      cat(sprintf("  %-12s %8.4f %8.4f %9.4f %8.4f %8.4f %9.4f  %s\n",
                  rownames(display)[j],
                  display[j, "RMSE_DA"],  display[j, "RMSE_EA"],
                  display[j, "PCM_RMSE"],
                  display[j, "MAE_DA"],   display[j, "MAE_EA"],
                  display[j, "PCM_MAE"],
                  classify_pcm(display[j, "PCM_RMSE"])))
    }
    cat("\n")
  }

  cat("PCM thresholds: < 0 Negative | 0-0.05 Weak | 0.05-0.10 Moderate | > 0.10 Strong\n")
  cat("Reference: Danks (2021), The DATA BASE for Advances in IS, 52(SI), 24-42.\n")
  invisible(x)
}

#' Plot PCM Results
#'
#' Creates a grouped barplot of PCM values for each mediation path, with one bar
#' per indicator and threshold lines for the weak/moderate/strong classification
#' boundaries.
#'
#' @param x An object of class \code{pcm_analysis} from \code{\link{assess_pcm}}.
#' @param metric Character, either \code{"RMSE"} (default) or \code{"MAE"}.
#' @param legend_pos Character giving the legend position (e.g. \code{"topright"},
#'   \code{"topleft"}, \code{"bottomright"}). Set to \code{NULL} or \code{FALSE}
#'   to suppress the legend entirely. Default \code{"topright"}.
#' @param bar_col Colour for positive-PCM bars (default \code{"steelblue"}).
#' @param neg_col Colour for negative-PCM bars (default \code{"tomato"}).
#' @param cex.labels Numeric expansion factor for the value labels printed on
#'   each bar. Default 0.8.
#' @param cex.legend Numeric expansion factor for the legend text. Default 0.75.
#' @param ... Additional graphical parameters passed to \code{\link[graphics]{barplot}}.
#'
#' @examples
#' \donttest{
#' # Customise legend placement and bar colours
#' # plot(pcm, legend_pos = "topleft", bar_col = "dodgerblue")
#'
#' # Suppress legend entirely
#' # plot(pcm, legend_pos = NULL)
#' }
#'
#' @export
plot.pcm_analysis <- function(x, metric = "RMSE",
                              legend_pos = "topright",
                              bar_col = "steelblue",
                              neg_col = "tomato",
                              cex.labels = 0.8,
                              cex.legend = 0.75, ...) {
  metric <- match.arg(metric, c("RMSE", "MAE"))
  pcm_col <- if (metric == "RMSE") "PCM_RMSE" else "PCM_MAE"

  n_paths <- length(x$pcm_results)

  # Set up multi-panel layout if multiple paths
  if (n_paths > 1) {
    old_par <- par(mfrow = c(1, n_paths), mar = c(5, 4, 4, 1))
    on.exit(par(old_par))
  }

  show_legend <- !is.null(legend_pos) && !identical(legend_pos, FALSE)

  for (i in seq_along(x$pcm_results)) {
    res  <- x$pcm_results[[i]]
    path <- x$mediation_paths[[i]]
    pcm_vals <- res$results[, pcm_col]

    # Bar colours: user-defined positive / negative
    bar_cols <- ifelse(pcm_vals >= 0,
                       adjustcolor(bar_col, 0.7),
                       adjustcolor(neg_col, 0.7))

    y_range <- range(c(pcm_vals, 0, 0.12))
    y_lim   <- c(min(y_range[1], -0.02), max(y_range[2] * 1.15, 0.12))

    bp <- barplot(pcm_vals,
                  names.arg = rownames(res$results),
                  col = bar_cols,
                  ylim = y_lim,
                  ylab = paste0("PCM (", metric, ")"),
                  main = sprintf("%s -> %s -> %s",
                                 path$antecedent, path$mediator, path$target),
                  las = 2, border = NA, ...)

    # Threshold lines
    abline(h = 0,    lty = 1, col = "grey40")
    abline(h = 0.05, lty = 2, col = "orange")
    abline(h = 0.10, lty = 2, col = "darkgreen")

    # Value labels on bars
    text(bp, pcm_vals, labels = sprintf("%.3f", pcm_vals),
         pos = ifelse(pcm_vals >= 0, 3, 1), cex = cex.labels)

    # Legend (only on first panel, if requested)
    if (show_legend && i == 1) {
      legend(legend_pos,
             legend = c("Weak (< 0.05)", "Moderate (0.05-0.10)", "Strong (> 0.10)"),
             lty = c(0, 2, 2),
             col = c(NA, "orange", "darkgreen"),
             pch = c(15, NA, NA),
             pt.cex = 1.5,
             bty = "n", cex = cex.legend)
    }
  }
}
