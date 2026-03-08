# =============================================================================
# feature_pos.R - Prediction-Oriented Segmentation (PLS-POS)
# =============================================================================
# This file implements PLS-POS, a deterministic hill-climbing approach to
# uncovering unobserved heterogeneity in PLS-SEM models, as described in:
#
# Becker, J.-M., Rai, A., Ringle, C. M., & Voelckner, F. (2013). Discovering
# Unobserved Heterogeneity in Structural Equation Models to Avert Validity
# Threats. MIS Quarterly, 37(3), 665-694.
#
# PLS-POS maximizes the sum of R-squared values across all endogenous
# constructs across K segments. Unlike FIMIX-PLS, it makes no distributional
# assumptions, uses hard (deterministic) segment assignment, and can detect
# heterogeneity in both structural and formative measurement models.
#
# Algorithm overview (Appendix B of Becker et al., 2013):
#   1. Create random initial partition of N observations into K groups
#   2. Estimate group-specific PLS models
#   3. Compute distance measure for every observation to every group
#   4. Reassign observations one at a time if it improves the objective
#   5. Repeat until no improving move exists (convergence)
#   6. Use multiple random starts to avoid local optima
# =============================================================================

# =============================================================================
# Internal helpers
# =============================================================================

#' Get endogenous constructs for the PLS-POS objective criterion.
#'
#' Excludes interaction constructs (names containing "*") because their
#' R-squared values are not meaningfully interpretable for segmentation.
#' @noRd
pos_endogenous <- function(model) {
  endo <- get_endogenous_constructs(model)
  endo[!grepl("\\*", endo)]
}

#' Estimate group-specific PLS models for each segment.
#'
#' Re-estimates the full PLS model on each segment's subset of data using
#' \code{seminr::rerun()}. Returns NULL (whole list) if any segment has
#' fewer than 2 observations; returns NULL for individual segments that
#' fail estimation.
#' @noRd
estimate_segment_models <- function(model, assignment, K) {
  segment_models <- vector("list", K)
  for (k in seq_len(K)) {
    idx <- which(assignment == k)
    if (length(idx) < 2) return(NULL)
    segment_data <- model$data[idx, , drop = FALSE]
    segment_models[[k]] <- tryCatch(
      suppressMessages(seminr::rerun(model, data = segment_data)),
      error = function(e) NULL
    )
  }
  segment_models
}

#' Compute the PLS-POS objective criterion.
#'
#' The objective is the sum of R-squared across all endogenous constructs
#' across all K segments: OBJ = sum_{k=1}^{K} sum_{b in endo} R^2_{b,k}.
#' Returns -Inf if any segment model is NULL (failed estimation).
#'
#' Note: seminr's rSquared is a matrix with rows "Rsq"/"AdjRsq" and
#' construct names as columns, so we index as rs["Rsq", endogenous].
#' @noRd
compute_pos_objective <- function(segment_models, endogenous) {
  total <- 0
  for (k in seq_along(segment_models)) {
    m <- segment_models[[k]]
    if (is.null(m)) return(-Inf)
    vals <- m$rSquared["Rsq", endogenous]
    total <- total + sum(vals, na.rm = TRUE)
  }
  total
}

#' Compute squared structural residuals for every observation.
#'
#' For each observation i, endogenous construct b, and segment k, computes
#' the squared structural residual: e^2_{i,b,k} = (y_{i,b} - yhat_{i,b,k})^2
#' where yhat uses global construct scores with segment-specific path
#' coefficients. This follows the PLS-POS distance measure definition
#' (Becker et al., 2013, Appendix B).
#'
#' @return An N x B x K array of squared residuals.
#' @noRd
compute_structural_residuals <- function(model, segment_models, endogenous) {
  N <- nrow(model$construct_scores)
  K <- length(segment_models)
  B <- length(endogenous)
  sm <- model$smMatrix
  scores <- model$construct_scores

  residuals <- array(NA_real_, dim = c(N, B, K))

  for (k in seq_len(K)) {
    seg_model <- segment_models[[k]]
    if (is.null(seg_model)) next

    for (b_idx in seq_len(B)) {
      b <- endogenous[b_idx]
      # Identify direct predictors of this endogenous construct
      predictors <- unname(sm[sm[, "target"] == b, "source"])
      if (length(predictors) == 0) next

      # Predicted value = X %*% beta_k (segment-specific coefficients)
      path_coefs <- seg_model$path_coef[predictors, b]
      yhat <- as.vector(scores[, predictors, drop = FALSE] %*% path_coefs)
      y <- scores[, b]
      residuals[, b_idx, k] <- (y - yhat)^2
    }
  }
  residuals
}

#' Compute the PLS-POS distance matrix D (N x K).
#'
#' For each observation i and segment k, the distance is:
#'   D_{i,k} = sum_b sqrt( e^2_{i,b,k} / sum_i(e^2_{i,b,k}) )
#'
#' This normalizes squared residuals by their column total so that each
#' endogenous construct contributes equally regardless of scale, then
#' aggregates across constructs. Smaller distance = better fit to segment k.
#' @noRd
compute_pos_distances <- function(sq_residuals) {
  N <- dim(sq_residuals)[1]
  B <- dim(sq_residuals)[2]
  K <- dim(sq_residuals)[3]

  distances <- matrix(0, nrow = N, ncol = K)

  for (k in seq_len(K)) {
    for (b in seq_len(B)) {
      e2_i <- sq_residuals[, b, k]
      sum_e2 <- sum(e2_i, na.rm = TRUE)
      # Guard against zero denominator (perfect fit for all observations)
      if (sum_e2 < .Machine$double.eps) next
      distances[, k] <- distances[, k] + sqrt(e2_i / sum_e2)
    }
  }
  distances
}

#' Build sorted candidate list for reassignment.
#'
#' For each observation, computes the distance improvement from moving it
#' from its current segment to its best alternative segment. Only observations
#' with positive improvement potential (d_current > d_best_alt) are included.
#'
#' @return A data.frame sorted by descending improvement potential, or NULL
#'   if no candidates have positive improvement.
#' @noRd
build_candidate_list <- function(distances, assignment, K) {
  N <- nrow(distances)
  obs_vec <- integer(0)
  from_vec <- integer(0)
  to_vec <- integer(0)
  diff_vec <- numeric(0)

  for (i in seq_len(N)) {
    current_k <- assignment[i]
    d_current <- distances[i, current_k]

    # Find the alternative segment with smallest distance
    alt_dists <- distances[i, ]
    alt_dists[current_k] <- Inf
    best_alt <- which.min(alt_dists)
    d_alt <- alt_dists[best_alt]
    diff_val <- d_current - d_alt

    if (diff_val > 0) {
      obs_vec <- c(obs_vec, i)
      from_vec <- c(from_vec, current_k)
      to_vec <- c(to_vec, best_alt)
      diff_vec <- c(diff_vec, diff_val)
    }
  }

  if (length(obs_vec) == 0) return(NULL)

  candidates <- data.frame(
    obs = obs_vec, from = from_vec, to = to_vec, diff = diff_vec,
    stringsAsFactors = FALSE
  )
  candidates[order(candidates$diff, decreasing = TRUE), ]
}

#' Create a random initial partition ensuring minimum segment size.
#'
#' First guarantees min_size observations per segment, then randomly assigns
#' the remaining observations. The final assignment is shuffled to randomize
#' observation order within segments.
#' @noRd
random_partition <- function(N, K, min_size) {
  guaranteed <- rep(seq_len(K), each = min_size)
  remaining <- N - K * min_size
  if (remaining > 0) {
    extra <- sample(seq_len(K), remaining, replace = TRUE)
    assignment <- sample(c(guaranteed, extra))
  } else {
    assignment <- sample(guaranteed)
  }
  assignment
}

#' Get segment-appropriate palette colors.
#' @noRd
pos_palette <- function(K) {
  if (K <= 8) {
    palette.colors(K, palette = "Set2")
  } else {
    palette.colors(n = K)
  }
}

# =============================================================================
# Main exported functions
# =============================================================================

#' PLS-POS: Prediction-Oriented Segmentation
#'
#' Performs prediction-oriented segmentation (PLS-POS) to uncover unobserved
#' heterogeneity in PLS-SEM models. PLS-POS uses a deterministic hill-climbing
#' approach to maximize the sum of R-squared values across all endogenous
#' constructs across K segments (Becker et al., 2013).
#'
#' Unlike FIMIX-PLS, PLS-POS makes no distributional assumptions and can
#' detect heterogeneity in both structural and (formative) measurement models.
#'
#' @param seminr_model An estimated SEMinR model from \code{estimate_pls()}.
#' @param K Integer >= 2. Number of segments to extract.
#' @param nstart Number of random starting partitions (default 10). The best
#'   solution (highest objective) is returned.
#' @param max_iter Maximum hill-climbing iterations per start (default 100).
#' @param search_depth Maximum number of candidates to evaluate per iteration
#'   before accepting the first improvement (default = N, the sample size).
#'   Use a smaller value for faster but potentially less optimal results.
#' @param min_segment_size Minimum observations per segment. Default is
#'   \code{max(10, max_predictors + 2)} where max_predictors is the largest
#'   number of predictors for any endogenous construct.
#' @param seed Random seed for reproducibility.
#'
#' @return An S3 object of class \code{"pos_analysis"} with components:
#'   \describe{
#'     \item{K}{Number of segments}
#'     \item{segment_assignment}{Integer vector of segment assignments (length N)}
#'     \item{segment_sizes}{Named integer vector of segment sizes}
#'     \item{segment_models}{List of K re-estimated seminr model objects}
#'     \item{segment_rsquared}{Matrix of R-squared values (endogenous x K)}
#'     \item{segment_paths}{List of K path coefficient matrices}
#'     \item{objective}{Sum of R-squared (objective criterion)}
#'     \item{converged}{Whether the best start converged}
#'     \item{iterations}{Number of iterations in the best start}
#'     \item{nstart}{Number of random starts attempted}
#'     \item{all_objectives}{Objective values for all starts}
#'   }
#'
#' @seealso [assess_fimix()] for probabilistic (EM-based) segmentation
#'
#' @references
#' Becker, J.-M., Rai, A., Ringle, C. M., & Voelckner, F. (2013). Discovering
#' Unobserved Heterogeneity in Structural Equation Models to Avert Validity
#' Threats. \emph{MIS Quarterly}, 37(3), 665-694.
#'
#' @examples
#' library(seminr)
#'
#' mobi_mm <- constructs(
#'   composite("Image",        multi_items("IMAG", 1:5)),
#'   composite("Expectation",  multi_items("CUEX", 1:3)),
#'   composite("Satisfaction", multi_items("CUSA", 1:3)),
#'   composite("Loyalty",      multi_items("CUSL", 1:3))
#' )
#' mobi_sm <- relationships(
#'   paths(from = "Image",       to = c("Expectation", "Satisfaction", "Loyalty")),
#'   paths(from = "Expectation", to = "Satisfaction"),
#'   paths(from = "Satisfaction", to = "Loyalty")
#' )
#' pls_model <- estimate_pls(data = mobi, measurement_model = mobi_mm,
#'                            structural_model = mobi_sm)
#'
#' pos_result <- assess_pos(pls_model, K = 2, nstart = 3, max_iter = 20,
#'                           seed = 123)
#' print(pos_result)
#' summary(pos_result)
#'
#' @export
assess_pos <- function(seminr_model,
                       K = 2,
                       nstart = 10,
                       max_iter = 100,
                       search_depth = NULL,
                       min_segment_size = NULL,
                       seed = 123) {
  # ---------------------------------------------------------------------------
  # Step 0: Validate inputs
  # ---------------------------------------------------------------------------
  if (!validate_seminr_model(seminr_model, "assess_pos")) {
    return(NULL)
  }

  if (!is.numeric(K) || length(K) != 1 || K < 2 || K != as.integer(K))
    stop("K must be an integer >= 2.")
  K <- as.integer(K)

  N <- nrow(seminr_model$construct_scores)
  endogenous <- pos_endogenous(seminr_model)
  if (length(endogenous) == 0)
    stop("No endogenous constructs found in the model.")

  # Default min_segment_size: enough for OLS estimation in each segment
  sm <- seminr_model$smMatrix
  max_preds <- max(vapply(endogenous, function(b) {
    sum(sm[, "target"] == b)
  }, integer(1)))
  if (is.null(min_segment_size)) {
    min_segment_size <- max(10L, max_preds + 2L)
  }

  if (N < K * min_segment_size) {
    stop("Sample size (", N, ") is too small for K=", K,
         " segments with min_segment_size=", min_segment_size, ".")
  }

  if (is.null(search_depth)) search_depth <- N

  # ---------------------------------------------------------------------------
  # Step 1: Multi-start hill-climbing (Algorithm Steps 1-6, Appendix B)
  # ---------------------------------------------------------------------------
  set.seed(seed)

  best_objective <- -Inf
  best_result <- NULL
  all_objectives <- numeric(nstart)

  for (r in seq_len(nstart)) {

    # --- Step 1.1: Random initial partition ---
    assignment <- random_partition(N, K, min_segment_size)

    # --- Step 1.2: Estimate initial segment models ---
    seg_models <- estimate_segment_models(seminr_model, assignment, K)
    if (is.null(seg_models) || any(vapply(seg_models, is.null, logical(1)))) {
      all_objectives[r] <- -Inf
      next
    }

    # --- Step 1.3: Distance-based initial reassignment ---
    # Assign each observation to the segment with smallest distance
    sq_resid <- compute_structural_residuals(seminr_model, seg_models, endogenous)
    dists <- compute_pos_distances(sq_resid)
    assignment <- apply(dists, 1, which.min)

    # Verify min segment sizes are still met after reassignment
    seg_sizes <- tabulate(assignment, nbins = K)
    if (any(seg_sizes < min_segment_size)) {
      all_objectives[r] <- -Inf
      next
    }

    # --- Step 1.4: Re-estimate after initial assignment ---
    seg_models <- estimate_segment_models(seminr_model, assignment, K)
    if (is.null(seg_models) || any(vapply(seg_models, is.null, logical(1)))) {
      all_objectives[r] <- -Inf
      next
    }

    obj <- compute_pos_objective(seg_models, endogenous)
    converged <- FALSE
    iter <- 0

    # --- Steps 2-6: Hill-climbing loop ---
    # Reassign ONE observation at a time, only if it improves the objective.
    # This ensures monotonic improvement of the objective criterion.
    for (iter_idx in seq_len(max_iter)) {
      iter <- iter_idx

      # Step 2: Compute distances for current partition
      sq_resid <- compute_structural_residuals(seminr_model, seg_models, endogenous)
      dists <- compute_pos_distances(sq_resid)

      # Step 3: Build candidate list sorted by improvement potential
      candidates <- build_candidate_list(dists, assignment, K)
      if (is.null(candidates) || nrow(candidates) == 0) {
        converged <- TRUE
        break
      }

      # Step 4-5: Try reassigning the most promising candidates
      improved <- FALSE
      n_tried <- 0
      for (c_idx in seq_len(nrow(candidates))) {
        if (n_tried >= search_depth) break
        n_tried <- n_tried + 1

        cand <- candidates[c_idx, ]
        from_k <- cand$from
        to_k <- cand$to

        # Enforce minimum segment size constraint
        from_size <- sum(assignment == from_k)
        if (from_size - 1 < min_segment_size) next

        # Tentative reassignment: move this observation
        new_assignment <- assignment
        new_assignment[cand$obs] <- to_k

        # Re-estimate only the two affected segments (efficiency)
        new_seg_models <- seg_models
        for (affected_k in c(from_k, to_k)) {
          idx <- which(new_assignment == affected_k)
          seg_data <- seminr_model$data[idx, , drop = FALSE]
          new_model <- tryCatch(
            suppressMessages(seminr::rerun(seminr_model, data = seg_data)),
            error = function(e) NULL
          )
          if (is.null(new_model)) break
          new_seg_models[[affected_k]] <- new_model
        }

        # Skip if estimation failed for either affected segment
        if (any(vapply(new_seg_models[c(from_k, to_k)], is.null, logical(1)))) next

        # Step 6: Accept move only if objective improves
        new_obj <- compute_pos_objective(new_seg_models, endogenous)
        if (new_obj > obj) {
          assignment <- new_assignment
          seg_models <- new_seg_models
          obj <- new_obj
          improved <- TRUE
          break
        }
      }

      # No improving move found = converged for this start
      if (!improved) {
        converged <- TRUE
        break
      }
    }

    all_objectives[r] <- obj

    # Keep track of the best solution across all starts
    if (obj > best_objective) {
      best_objective <- obj
      best_result <- list(
        assignment = assignment,
        seg_models = seg_models,
        objective = obj,
        converged = converged,
        iterations = iter
      )
    }
  }

  if (is.null(best_result)) {
    stop("PLS-POS failed to find a valid segmentation across all ", nstart,
         " random starts. Consider increasing min_segment_size or decreasing K.")
  }

  # ---------------------------------------------------------------------------
  # Step 7: Build result object from best solution
  # ---------------------------------------------------------------------------
  assignment <- best_result$assignment
  seg_models <- best_result$seg_models
  seg_sizes <- as.integer(tabulate(assignment, nbins = K))
  names(seg_sizes) <- paste("Segment", seq_len(K))

  # R-squared matrix: endogenous constructs (rows) x segments (cols)
  seg_rsq <- matrix(NA_real_, nrow = length(endogenous), ncol = K,
                    dimnames = list(endogenous, paste("Segment", seq_len(K))))
  for (k in seq_len(K)) {
    seg_rsq[, k] <- seg_models[[k]]$rSquared["Rsq", endogenous]
  }

  # Extract path coefficient matrices for convenient access
  seg_paths <- lapply(seg_models, function(m) m$path_coef)

  result <- list(
    K = K,
    segment_assignment = assignment,
    segment_sizes = seg_sizes,
    segment_models = seg_models,
    segment_rsquared = seg_rsq,
    segment_paths = seg_paths,
    objective = best_result$objective,
    converged = best_result$converged,
    iterations = best_result$iterations,
    nstart = nstart,
    all_objectives = all_objectives,
    endogenous = endogenous,
    pls_model = seminr_model,
    n_obs = N
  )
  class(result) <- c("pos_analysis", class(result))
  result
}

#' Compare PLS-POS Solutions Across K Values
#'
#' Runs \code{\link{assess_pos}} for each K in \code{K_range} and returns a
#' comparison table of objective criterion (sum of R-squared) values. This
#' helps identify the optimal number of segments.
#'
#' @inheritParams assess_pos
#' @param K_range Integer vector of K values to compare (default \code{2:5}).
#'
#' @return An S3 object of class \code{"pos_comparison"} with components:
#'   \describe{
#'     \item{solutions}{Named list of \code{pos_analysis} objects}
#'     \item{fit_table}{Data frame comparing K values}
#'     \item{K_range}{The K values compared}
#'   }
#'
#' @seealso [assess_fimix_compare()] for comparing FIMIX solutions across K
#'
#' @examples
#' \donttest{
#' library(seminr)
#'
#' mobi_mm <- constructs(
#'   composite("Image",        multi_items("IMAG", 1:5)),
#'   composite("Expectation",  multi_items("CUEX", 1:3)),
#'   composite("Satisfaction", multi_items("CUSA", 1:3)),
#'   composite("Loyalty",      multi_items("CUSL", 1:3))
#' )
#' mobi_sm <- relationships(
#'   paths(from = "Image",       to = c("Expectation", "Satisfaction", "Loyalty")),
#'   paths(from = "Expectation", to = "Satisfaction"),
#'   paths(from = "Satisfaction", to = "Loyalty")
#' )
#' pls_model <- estimate_pls(data = mobi, measurement_model = mobi_mm,
#'                            structural_model = mobi_sm)
#'
#' pos_comp <- assess_pos_compare(pls_model, K_range = 2:4,
#'                                 nstart = 3, max_iter = 20, seed = 123)
#' print(pos_comp)
#' plot(pos_comp)
#' }
#'
#' @export
assess_pos_compare <- function(seminr_model,
                               K_range = 2:5,
                               nstart = 10,
                               max_iter = 100,
                               search_depth = NULL,
                               min_segment_size = NULL,
                               seed = 123) {
  if (!validate_seminr_model(seminr_model, "assess_pos_compare")) {
    return(NULL)
  }

  solutions <- list()
  fit_rows <- list()

  for (K in K_range) {
    sol <- tryCatch(
      assess_pos(seminr_model, K = K, nstart = nstart, max_iter = max_iter,
                 search_depth = search_depth, min_segment_size = min_segment_size,
                 seed = seed),
      error = function(e) {
        message("K=", K, ": ", e$message)
        NULL
      }
    )

    label <- paste0("K", K)
    solutions[[label]] <- sol

    if (!is.null(sol)) {
      avg_rsq <- mean(colSums(sol$segment_rsquared, na.rm = TRUE))
      fit_rows[[label]] <- data.frame(
        K = K,
        Sum_R2 = sol$objective,
        Avg_R2_per_segment = avg_rsq,
        Converged = sol$converged,
        Iterations = sol$iterations,
        stringsAsFactors = FALSE
      )
    } else {
      fit_rows[[label]] <- data.frame(
        K = K, Sum_R2 = NA_real_, Avg_R2_per_segment = NA_real_,
        Converged = NA, Iterations = NA_integer_,
        stringsAsFactors = FALSE
      )
    }
  }

  fit_table <- do.call(rbind, fit_rows)
  rownames(fit_table) <- NULL

  result <- list(
    solutions = solutions,
    fit_table = fit_table,
    K_range = K_range,
    pls_model = seminr_model
  )
  class(result) <- c("pos_comparison", class(result))
  result
}

#' Extract Segment-Specific PLS Models from PLS-POS Results
#'
#' Returns a list of fully re-estimated PLS models, one per segment, using the
#' hard segment assignment from PLS-POS.
#'
#' @param pos_result A \code{pos_analysis} object from \code{assess_pos()}.
#'
#' @return A named list of K seminr model objects.
#'
#' @seealso [assess_pos()] for running PLS-POS segmentation
#'
#' @export
pos_segments <- function(pos_result) {
  if (!inherits(pos_result, "pos_analysis"))
    stop("pos_result must be a pos_analysis object from assess_pos().")
  pos_result$segment_models
}

# =============================================================================
# S3 methods: print, summary, plot
# =============================================================================

#' @export
print.pos_analysis <- function(x, ...) {
  cat("PLS-POS Analysis\n")
  cat("================\n")
  cat("Segments:", x$K, "\n")
  cat("Observations:", x$n_obs, "\n")
  cat("Converged:", ifelse(x$converged, "Yes", "No"),
      "(", x$iterations, "iterations )\n")
  cat("Random starts:", x$nstart, "\n")
  cat("Objective (Sum R\u00b2):", round(x$objective, 4), "\n\n")

  # Segment sizes table
  cat("Segment Sizes:\n")
  size_df <- data.frame(
    Segment = names(x$segment_sizes),
    Size = as.integer(x$segment_sizes),
    Proportion = round(x$segment_sizes / x$n_obs, 4),
    stringsAsFactors = FALSE
  )
  print(size_df, row.names = FALSE)

  # R-squared comparison: segment vs global
  cat("\nR\u00b2 per Endogenous Construct:\n")
  rsq_df <- as.data.frame(round(x$segment_rsquared, 4))
  rsq_df$Global <- round(x$pls_model$rSquared["Rsq", x$endogenous], 4)
  print(rsq_df)

  # Segment-specific path coefficients
  cat("\nSegment Path Coefficients:\n")
  for (k in seq_len(x$K)) {
    cat("\n  Segment", k, ":\n")
    pm <- x$segment_paths[[k]]
    nonzero <- which(pm != 0, arr.ind = TRUE)
    if (nrow(nonzero) == 0) next
    path_df <- data.frame(
      From = rownames(pm)[nonzero[, 1]],
      To = colnames(pm)[nonzero[, 2]],
      Coefficient = round(pm[nonzero], 4),
      stringsAsFactors = FALSE
    )
    print(path_df, row.names = FALSE)
  }

  invisible(x)
}

#' @export
summary.pos_analysis <- function(object, ...) {
  result <- list(
    K = object$K,
    n_obs = object$n_obs,
    objective = object$objective,
    converged = object$converged,
    iterations = object$iterations,
    nstart = object$nstart,
    all_objectives = object$all_objectives,
    segment_sizes = object$segment_sizes,
    segment_rsquared = object$segment_rsquared,
    segment_paths = object$segment_paths,
    global_rsquared = object$pls_model$rSquared,
    global_paths = object$pls_model$path_coef,
    endogenous = object$endogenous
  )
  class(result) <- c("summary.pos_analysis", class(result))
  result
}

#' @export
print.summary.pos_analysis <- function(x, digits = 4, ...) {
  cat("PLS-POS Analysis \u2014 Detailed Summary\n")
  cat("======================================\n")
  cat("Segments:", x$K, "| Observations:", x$n_obs, "\n")
  cat("Objective (Sum R\u00b2):", round(x$objective, digits), "\n")
  cat("Converged:", ifelse(x$converged, "Yes", "No"),
      "(", x$iterations, "iterations )\n")
  cat("Random starts:", x$nstart, "\n")
  cat("Start objectives:", paste(round(x$all_objectives, digits), collapse = ", "), "\n\n")

  # Segment size breakdown
  cat("Segment Sizes:\n")
  for (k in seq_along(x$segment_sizes)) {
    cat("  Segment ", k, ": ", x$segment_sizes[k],
        " (", round(100 * x$segment_sizes[k] / x$n_obs, 1), "%)\n", sep = "")
  }

  # Side-by-side R-squared: segment vs global
  cat("\nR\u00b2 Comparison (Segment vs Global):\n")
  rsq_df <- as.data.frame(round(x$segment_rsquared, digits))
  rsq_df$Global <- round(x$global_rsquared["Rsq", x$endogenous], digits)
  print(rsq_df)

  # Side-by-side path coefficients: global + each segment
  cat("\nPath Coefficients per Segment:\n")
  global_pm <- x$global_paths
  nonzero <- which(global_pm != 0, arr.ind = TRUE)
  if (nrow(nonzero) > 0) {
    path_labels <- paste(rownames(global_pm)[nonzero[, 1]], "->",
                         colnames(global_pm)[nonzero[, 2]])
    coef_df <- data.frame(Path = path_labels, stringsAsFactors = FALSE)
    coef_df$Global <- round(global_pm[nonzero], digits)
    for (k in seq_along(x$segment_paths)) {
      coef_df[[paste0("Seg.", k)]] <- round(x$segment_paths[[k]][nonzero], digits)
    }
    print(coef_df, row.names = FALSE)
  }

  invisible(x)
}

#' @export
plot.pos_analysis <- function(x,
                              type = c("segments", "rsquared", "paths"),
                              ...) {
  type <- match.arg(type)
  cols <- pos_palette(x$K)

  if (type == "segments") {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    par(mar = c(5, 5, 4, 2))

    barplot(x$segment_sizes / x$n_obs * 100,
            names.arg = paste("Seg.", seq_len(x$K)),
            col = cols, ylim = c(0, 100),
            ylab = "Proportion (%)",
            main = paste("PLS-POS: Segment Proportions (K =", x$K, ")"),
            ...)

  } else if (type == "rsquared") {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    par(mar = c(8, 5, 4, 2))

    barplot(t(x$segment_rsquared), beside = TRUE, col = cols,
            las = 2, ylab = expression(R^2),
            main = paste("PLS-POS: R\u00b2 per Construct (K =", x$K, ")"),
            ...)
    legend("topright", legend = paste("Seg.", seq_len(x$K)),
           fill = cols, bty = "n")

  } else if (type == "paths") {
    pm_ref <- x$segment_paths[[1]]
    nonzero <- which(pm_ref != 0, arr.ind = TRUE)
    if (nrow(nonzero) == 0) {
      message("No non-zero path coefficients to plot.")
      return(invisible(x))
    }

    path_labels <- paste(rownames(pm_ref)[nonzero[, 1]], "->",
                         colnames(pm_ref)[nonzero[, 2]])
    coef_matrix <- matrix(NA_real_, nrow = length(path_labels), ncol = x$K)
    for (k in seq_len(x$K)) {
      coef_matrix[, k] <- x$segment_paths[[k]][nonzero]
    }
    rownames(coef_matrix) <- path_labels
    colnames(coef_matrix) <- paste("Seg.", seq_len(x$K))

    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    par(mar = c(8, 5, 4, 2))

    barplot(t(coef_matrix), beside = TRUE, col = cols,
            las = 2, ylab = "Path Coefficient",
            main = paste("PLS-POS: Path Coefficients (K =", x$K, ")"),
            ...)
    legend("topright", legend = colnames(coef_matrix),
           fill = cols, bty = "n")
  }

  invisible(x)
}

#' @export
print.pos_comparison <- function(x, ...) {
  cat("PLS-POS Comparison\n")
  cat("==================\n")
  cat("K range:", paste(x$K_range, collapse = ", "), "\n\n")
  print(x$fit_table, row.names = FALSE)
  invisible(x)
}

#' @export
summary.pos_comparison <- function(object, ...) {
  object
}

#' @export
plot.pos_comparison <- function(x, ...) {
  ft <- x$fit_table
  valid <- !is.na(ft$Sum_R2)

  if (sum(valid) < 2) {
    message("Not enough valid solutions to plot.")
    return(invisible(x))
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mar = c(5, 5, 4, 2))

  plot(ft$K[valid], ft$Sum_R2[valid],
       type = "b", pch = 19, lwd = 2,
       xlab = "Number of Segments (K)",
       ylab = expression("Objective (Sum " * R^2 * ")"),
       main = "PLS-POS: Objective Criterion vs K",
       xaxt = "n", ...)
  axis(1, at = ft$K[valid])

  invisible(x)
}
