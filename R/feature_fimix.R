# =============================================================================
# feature_fimix.R - FIMIX-PLS (Finite Mixture PLS)
# =============================================================================
# Implements EM-based latent class segmentation for PLS-SEM structural models.
# Observations are probabilistically assigned to K segments, each with
# segment-specific structural path coefficients.
#
# References:
# - Hahn, C., Johnson, M. D., Herrmann, A., & Huber, F. (2002). Capturing
#   Customer Heterogeneity using a Finite Mixture PLS Approach. Schmalenbach
#   Business Review, 54, 243-269.
# - Sarstedt, M., Becker, J.-M., Ringle, C. M. & Schwaiger, M. (2011).
#   Uncovering and Treating Unobserved Heterogeneity with FIMIX-PLS.
#   Schmalenbach Business Review, 63(1), 34-62.
# - Hair, J. F. et al. (2016). Identifying and Treating Unobserved
#   Heterogeneity with FIMIX-PLS: Part I. EBR, 28(1), 63-76.
# - Ringle, C. M., Sarstedt, M. & Mooi, E. A. (2010). Response-Based
#   Segmentation Using Finite Mixture PLS. In Data Mining, Springer, 19-49.
# =============================================================================

# =============================================================================
# INTERNAL HELPERS
# =============================================================================

#' Log-sum-exp trick for numerical stability.
#' @noRd
log_sum_exp <- function(x) {
  m <- max(x)
  if (is.infinite(m)) return(m)
  m + log(sum(exp(x - m)))
}

#' Extract structural equations from a seminr model.
#'
#' Returns a list of equations, each with target name, predictor names,
#' y vector (endogenous scores), and X matrix (intercept + predictor scores).
#' @noRd
extract_structural_equations <- function(model) {
  sm <- model$smMatrix
  scores <- model$construct_scores
  endogenous <- unique(sm[, "target"])

  equations <- list()
  for (target in endogenous) {
    predictors <- unname(sm[sm[, "target"] == target, "source"])
    y <- scores[, target]
    X <- cbind(1, scores[, predictors, drop = FALSE])
    colnames(X)[1] <- "(Intercept)"

    equations[[target]] <- list(
      target = target,
      predictors = predictors,
      y = y,
      X = X
    )
  }

  equations
}

#' Initialize random posterior probabilities via Dirichlet-like sampling.
#' @noRd
init_random_posteriors <- function(N, K) {
  # Sample from Dirichlet(1,...,1) = normalized exponentials of uniform
  raw <- matrix(stats::rexp(N * K, rate = 1), nrow = N, ncol = K)
  raw / rowSums(raw)
}

#' Weighted OLS for one structural equation in one segment.
#'
#' @param y Numeric vector of endogenous construct scores.
#' @param X Matrix with intercept column + predictor scores.
#' @param weights Numeric vector of observation weights (posterior probs).
#' @return List with coefficients and residual variance.
#' @noRd
weighted_ols <- function(y, X, weights) {
  W <- as.vector(weights)
  XtWX <- crossprod(X * W, X)
  XtWy <- crossprod(X * W, y)

  beta <- tryCatch(
    as.vector(solve(XtWX, XtWy)),
    error = function(e) rep(NA_real_, ncol(X))
  )

  if (any(is.na(beta))) {
    return(list(coefficients = beta,
                variance = NA_real_))
  }

  resid <- as.vector(y - X %*% beta)
  sigma2 <- sum(W * resid^2) / sum(W)
  # Floor variance to prevent degenerate densities
  sigma2 <- max(sigma2, .Machine$double.eps * 100)

  list(coefficients = beta, variance = sigma2)
}

#' Count the number of free parameters in a FIMIX model.
#'
#' Q = K * sum_j(p_j + 2) + (K - 1)
#' where p_j = number of predictors for equation j, +1 intercept, +1 variance.
#' @noRd
count_fimix_parameters <- function(equations, K) {
  params_per_eq <- vapply(equations, function(eq) {
    length(eq$predictors) + 2  # coefficients + intercept + variance
  }, numeric(1))

  K * sum(params_per_eq) + (K - 1)
}

#' Compute normed entropy statistic.
#'
#' EN = 1 - (-sum P*log(P)) / (N * log(K))
#' Closer to 1 = better classification quality.
#' @noRd
compute_entropy <- function(posteriors, K) {
  if (K <= 1) return(NA_real_)

  N <- nrow(posteriors)
  # Clamp to avoid log(0)
  P <- pmax(posteriors, .Machine$double.eps)
  raw_entropy <- -sum(P * log(P))
  1 - raw_entropy / (N * log(K))
}

#' Compute all FIMIX information criteria.
#' @noRd
compute_fimix_criteria <- function(lnL, Q, N, K, posteriors) {
  c(
    lnL  = lnL,
    AIC  = -2 * lnL + 2 * Q,
    AIC3 = -2 * lnL + 3 * Q,
    AIC4 = -2 * lnL + 4 * Q,
    BIC  = -2 * lnL + Q * log(N),
    CAIC = -2 * lnL + Q * (log(N) + 1),
    HQ   = -2 * lnL + 2 * Q * log(log(N)),
    MDL5 = -2 * lnL + (Q / 2) * log(N),
    EN   = compute_entropy(posteriors, K)
  )
}

#' Run a single EM iteration set for FIMIX-PLS.
#'
#' @param equations List from extract_structural_equations().
#' @param K Number of segments.
#' @param max_iter Maximum EM iterations.
#' @param stop_criterion Convergence threshold on log-likelihood change.
#' @param posteriors_init N x K matrix of initial posterior probabilities.
#' @return List with final parameters, posteriors, log-likelihood, convergence info.
#' @noRd
run_fimix_em <- function(equations, K, max_iter, stop_criterion, posteriors_init) {
  N <- length(equations[[1]]$y)
  J <- length(equations)
  posteriors <- posteriors_init
  lnL_old <- -Inf
  converged <- FALSE
  iter <- 0L

  # Storage for segment parameters
  # segment_params[[eq_name]][[k]] = list(coefficients, variance)
  segment_params <- list()

  for (it in seq_len(max_iter)) {
    iter <- it

    # --- M-step ---
    N_k <- colSums(posteriors)
    pi_k <- N_k / N

    # Check for degenerate segments
    min_effective <- min(N_k)
    if (min_effective < 1) {
      return(list(converged = FALSE, lnL = -Inf, iter = iter))
    }

    for (eq_name in names(equations)) {
      eq <- equations[[eq_name]]
      segment_params[[eq_name]] <- list()

      for (k in seq_len(K)) {
        result <- weighted_ols(eq$y, eq$X, posteriors[, k])
        if (any(is.na(result$coefficients))) {
          return(list(converged = FALSE, lnL = -Inf, iter = iter))
        }
        segment_params[[eq_name]][[k]] <- result
      }
    }

    # --- E-step ---
    # Compute log component densities: log_dens[i, k]
    log_dens <- matrix(0, nrow = N, ncol = K)

    for (k in seq_len(K)) {
      log_dens[, k] <- log(pi_k[k])
      for (eq_name in names(equations)) {
        eq <- equations[[eq_name]]
        sp <- segment_params[[eq_name]][[k]]
        mu <- as.vector(eq$X %*% sp$coefficients)
        sigma <- sqrt(sp$variance)
        log_dens[, k] <- log_dens[, k] + dnorm(eq$y, mean = mu,
                                                  sd = sigma, log = TRUE)
      }
    }

    # Posteriors via log-sum-exp
    log_denom <- apply(log_dens, 1, log_sum_exp)
    posteriors <- exp(log_dens - log_denom)

    # Clamp to avoid exact 0/1
    posteriors <- pmax(posteriors, .Machine$double.eps)
    posteriors <- posteriors / rowSums(posteriors)

    # --- Convergence check ---
    lnL_new <- sum(log_denom)
    if (abs(lnL_new - lnL_old) < stop_criterion) {
      converged <- TRUE
      lnL_old <- lnL_new
      break
    }
    lnL_old <- lnL_new
  }

  list(
    converged = converged,
    lnL = lnL_old,
    iter = iter,
    posteriors = posteriors,
    pi_k = pi_k,
    segment_params = segment_params
  )
}

#' Build segment-specific path coefficient matrices from EM results.
#'
#' Each segment gets a matrix with the same structure as model$path_coef,
#' but with segment-specific coefficients filled in.
#' @noRd
build_segment_path_matrices <- function(equations, segment_params, model, K) {
  template <- model$path_coef * 0  # Same structure, all zeros
  construct_names <- colnames(template)

  segment_paths <- list()
  segment_intercepts <- list()

  for (k in seq_len(K)) {
    path_mat <- template
    intercepts <- numeric(length(equations))
    names(intercepts) <- names(equations)

    for (eq_name in names(equations)) {
      eq <- equations[[eq_name]]
      sp <- segment_params[[eq_name]][[k]]
      coefs <- sp$coefficients

      intercepts[eq_name] <- coefs[1]  # First coefficient is intercept

      for (p_idx in seq_along(eq$predictors)) {
        pred <- eq$predictors[p_idx]
        if (pred %in% construct_names && eq$target %in% construct_names) {
          path_mat[pred, eq$target] <- coefs[p_idx + 1]
        }
      }
    }

    segment_paths[[k]] <- path_mat
    segment_intercepts[[k]] <- intercepts
  }

  list(paths = segment_paths, intercepts = segment_intercepts)
}

# =============================================================================
# MAIN ENTRY POINTS
# =============================================================================

#' FIMIX-PLS Analysis for PLS-SEM
#'
#' \code{assess_fimix} runs the FIMIX-PLS (Finite Mixture PLS) procedure to
#' identify unobserved heterogeneity in PLS-SEM structural models. The EM
#' algorithm probabilistically assigns observations to \code{K} latent
#' segments, each with segment-specific structural path coefficients.
#'
#' FIMIX-PLS operates on construct scores from the estimated PLS model.
#' It fits a mixture of regressions to each structural equation simultaneously,
#' estimating segment-specific intercepts, path coefficients, and residual
#' variances. Multiple random starts (\code{nstart}) are used to avoid local
#' optima; the solution with the highest log-likelihood is retained.
#'
#' @param seminr_model An estimated SEMinR model from \code{estimate_pls()}.
#' @param K Integer; the number of segments (default 2). Must be >= 2.
#' @param nstart Integer; number of random EM starts (default 10).
#' @param max_iter Integer; maximum EM iterations per start (default 5000).
#' @param stop_criterion Numeric; convergence threshold on log-likelihood
#'   change (default 1e-6).
#' @param seed Random seed for reproducibility (default 123).
#'
#' @return An object of class \code{fimix_analysis} containing:
#'   \item{K}{Number of segments}
#'   \item{segment_proportions}{Named numeric vector of mixing proportions}
#'   \item{posterior}{N x K matrix of posterior segment probabilities}
#'   \item{segment_assignment}{Integer vector of hard assignments (argmax)}
#'   \item{segment_paths}{List of K path coefficient matrices}
#'   \item{segment_intercepts}{List of K intercept vectors}
#'   \item{segment_variances}{J x K matrix of residual variances}
#'   \item{log_likelihood}{Final log-likelihood}
#'   \item{n_parameters}{Number of free parameters}
#'   \item{info_criteria}{Named vector of information criteria}
#'   \item{converged}{Logical; whether EM converged}
#'   \item{iterations}{Number of EM iterations}
#'   \item{n_starts_completed}{Number of random starts completed}
#'   \item{pls_model}{The original estimated seminr model}
#'   \item{n_obs}{Sample size}
#'
#' @references
#' Hahn, C., Johnson, M. D., Herrmann, A., & Huber, F. (2002). Capturing
#' Customer Heterogeneity using a Finite Mixture PLS Approach. Schmalenbach
#' Business Review, 54, 243-269.
#'
#' Sarstedt, M., Becker, J.-M., Ringle, C. M. & Schwaiger, M. (2011).
#' Uncovering and Treating Unobserved Heterogeneity with FIMIX-PLS: Which
#' Model Selection Criterion Provides an Appropriate Number of Segments?
#' Schmalenbach Business Review, 63(1), 34-62.
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
#' fimix_result <- assess_fimix(mobi_pls, K = 2, nstart = 5, seed = 123)
#' print(fimix_result)
#' summary(fimix_result)
#' }
#'
#' @seealso \code{\link{assess_fimix_compare}} for comparing multiple K values
#'
#' @export
assess_fimix <- function(seminr_model,
                          K = 2,
                          nstart = 10,
                          max_iter = 5000,
                          stop_criterion = 1e-6,
                          seed = 123) {

  # ---------------------------------------------------------------------------
  # Step 1: Validate inputs
  # ---------------------------------------------------------------------------
  if (!validate_seminr_model(seminr_model, "assess_fimix")) return(NULL)

  if (!is.numeric(K) || length(K) != 1 || K < 2 || K != as.integer(K)) {
    stop("K must be an integer >= 2.", call. = FALSE)
  }

  if (!is.numeric(nstart) || length(nstart) != 1 || nstart < 1) {
    stop("nstart must be an integer >= 1.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Step 2: Extract structural equations
  # ---------------------------------------------------------------------------
  equations <- extract_structural_equations(seminr_model)

  if (length(equations) == 0) {
    stop("No endogenous constructs found in the structural model.", call. = FALSE)
  }

  N <- nrow(seminr_model$construct_scores)
  Q <- count_fimix_parameters(equations, K)

  # Minimum observations per segment heuristic

  max_predictors <- max(vapply(equations, function(eq) length(eq$predictors), integer(1)))
  min_segment_obs <- max(max_predictors + 2, 10)
  if (N / K < min_segment_obs) {
    warning("Sample size (", N, ") may be too small for K = ", K,
            " segments. Minimum ~", min_segment_obs, " observations per segment ",
            "recommended.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Step 3: Run EM with multiple random starts
  # ---------------------------------------------------------------------------
  set.seed(seed)

  best_result <- NULL
  best_lnL <- -Inf
  n_completed <- 0L

  for (r in seq_len(nstart)) {
    posteriors_init <- init_random_posteriors(N, K)
    result <- run_fimix_em(equations, K, max_iter, stop_criterion,
                            posteriors_init)
    n_completed <- n_completed + 1L

    if (result$converged && result$lnL > best_lnL) {
      best_result <- result
      best_lnL <- result$lnL
    }
  }

  if (is.null(best_result)) {
    warning("FIMIX-PLS did not converge for K = ", K,
            " in any of ", nstart, " random starts. ",
            "Consider reducing K or increasing max_iter.", call. = FALSE)
    return(NULL)
  }

  # ---------------------------------------------------------------------------
  # Step 4: Construct return object
  # ---------------------------------------------------------------------------
  # Build segment path matrices and intercepts
  built <- build_segment_path_matrices(equations, best_result$segment_params,
                                        seminr_model, K)

  # Segment variances: J x K matrix
  J <- length(equations)
  var_matrix <- matrix(NA_real_, nrow = J, ncol = K,
                        dimnames = list(names(equations),
                                        paste0("Segment_", seq_len(K))))
  for (eq_name in names(equations)) {
    for (k in seq_len(K)) {
      var_matrix[eq_name, k] <- best_result$segment_params[[eq_name]][[k]]$variance
    }
  }

  # Segment assignments (hard clustering)
  assignment <- apply(best_result$posteriors, 1, which.max)

  # Segment proportions
  pi_k <- best_result$pi_k
  names(pi_k) <- paste0("Segment_", seq_len(K))

  # Segment sizes (hard assignment counts)
  segment_sizes <- tabulate(assignment, nbins = K)
  names(segment_sizes) <- paste0("Segment_", seq_len(K))

  # Information criteria
  info_criteria <- compute_fimix_criteria(best_result$lnL, Q, N, K,
                                           best_result$posteriors)

  # Column names for posteriors
  colnames(best_result$posteriors) <- paste0("Segment_", seq_len(K))

  result <- list(
    K                   = K,
    segment_proportions = pi_k,
    segment_sizes       = segment_sizes,
    posterior           = best_result$posteriors,
    segment_assignment  = assignment,
    segment_paths       = built$paths,
    segment_intercepts  = built$intercepts,
    segment_variances   = var_matrix,
    log_likelihood      = best_result$lnL,
    n_parameters        = Q,
    info_criteria       = info_criteria,
    converged           = best_result$converged,
    iterations          = best_result$iter,
    n_starts_completed  = n_completed,
    pls_model           = seminr_model,
    n_obs               = N
  )

  class(result) <- c("fimix_analysis", class(result))
  result
}

#' Compare FIMIX-PLS Solutions Across Multiple K Values
#'
#' \code{assess_fimix_compare} runs \code{\link{assess_fimix}} for each value
#' of K in \code{K_range} and produces a comparison table of information
#' criteria to guide segment number selection.
#'
#' The recommended approach is to inspect AIC3 and CAIC jointly: when both
#' agree on K, accuracy reaches approximately 84\% (Sarstedt et al., 2011).
#' AIC4 is the best-performing single criterion. Entropy (EN) indicates
#' classification quality but should not be used alone to select K.
#'
#' @inheritParams assess_fimix
#' @param K_range Integer vector of K values to evaluate (default 2:5).
#'
#' @return An object of class \code{fimix_comparison} containing:
#'   \item{solutions}{Named list of \code{fimix_analysis} objects (one per K)}
#'   \item{fit_table}{Data frame with K as rows and criteria as columns}
#'   \item{K_range}{The evaluated K values}
#'   \item{pls_model}{The original model}
#'
#' @references
#' Sarstedt, M., Becker, J.-M., Ringle, C. M. & Schwaiger, M. (2011).
#' Uncovering and Treating Unobserved Heterogeneity with FIMIX-PLS: Which
#' Model Selection Criterion Provides an Appropriate Number of Segments?
#' Schmalenbach Business Review, 63(1), 34-62.
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
#' fimix_comp <- assess_fimix_compare(mobi_pls, K_range = 2:4,
#'                                      nstart = 5, seed = 123)
#' print(fimix_comp)
#' plot(fimix_comp)
#' }
#'
#' @seealso \code{\link{assess_fimix}} for single-K FIMIX analysis
#'
#' @export
assess_fimix_compare <- function(seminr_model,
                                   K_range = 2:5,
                                   nstart = 10,
                                   max_iter = 5000,
                                   stop_criterion = 1e-6,
                                   seed = 123) {

  if (!validate_seminr_model(seminr_model, "assess_fimix_compare")) return(NULL)

  if (!is.numeric(K_range) || length(K_range) < 1 || any(K_range < 2)) {
    stop("K_range must be a vector of integers >= 2.", call. = FALSE)
  }

  solutions <- list()
  criteria_list <- list()

  for (k in K_range) {
    sol <- assess_fimix(seminr_model, K = k, nstart = nstart,
                          max_iter = max_iter,
                          stop_criterion = stop_criterion,
                          seed = seed)
    sol_name <- paste0("K", k)
    solutions[[sol_name]] <- sol

    if (!is.null(sol)) {
      criteria_list[[sol_name]] <- c(K = k, sol$info_criteria)
    } else {
      criteria_list[[sol_name]] <- c(K = k,
                                      lnL = NA, AIC = NA, AIC3 = NA,
                                      AIC4 = NA, BIC = NA, CAIC = NA,
                                      HQ = NA, MDL5 = NA, EN = NA)
    }
  }

  fit_table <- as.data.frame(do.call(rbind, criteria_list))
  rownames(fit_table) <- NULL

  result <- list(
    solutions = solutions,
    fit_table = fit_table,
    K_range   = K_range,
    pls_model = seminr_model
  )

  class(result) <- c("fimix_comparison", class(result))
  result
}

# =============================================================================
# S3 METHODS — fimix_analysis
# =============================================================================

#' @export
print.fimix_analysis <- function(x, ...) {
  cat("FIMIX-PLS Analysis\n")
  cat("==================\n")
  cat("Segments:", x$K, "\n")
  cat("Observations:", x$n_obs, "\n")
  cat("Converged:", if (x$converged) "Yes" else "No",
      "(", x$iterations, "iterations )\n")
  cat("Random starts:", x$n_starts_completed, "\n")

  cat("\nSegment Proportions:\n")
  seg_tbl <- data.frame(
    Segment = paste0("Segment ", seq_len(x$K)),
    Proportion = round(x$segment_proportions, 4),
    Size = x$segment_sizes,
    check.names = FALSE,
    row.names = NULL
  )
  print(seg_tbl, row.names = FALSE, ...)

  cat("\nFit Criteria:\n")
  ic <- x$info_criteria
  ic_tbl <- data.frame(
    lnL = round(ic["lnL"], 2),
    AIC = round(ic["AIC"], 2),
    AIC3 = round(ic["AIC3"], 2),
    AIC4 = round(ic["AIC4"], 2),
    BIC = round(ic["BIC"], 2),
    CAIC = round(ic["CAIC"], 2),
    EN = round(ic["EN"], 4),
    check.names = FALSE,
    row.names = NULL
  )
  print(ic_tbl, row.names = FALSE, ...)

  cat("\nSegment Path Coefficients:\n")
  for (k in seq_len(x$K)) {
    cat("\n  Segment", k, ":\n")
    pm <- x$segment_paths[[k]]
    # Print only non-zero entries
    nonzero <- which(pm != 0, arr.ind = TRUE)
    if (nrow(nonzero) > 0) {
      path_tbl <- data.frame(
        From = rownames(pm)[nonzero[, 1]],
        To = colnames(pm)[nonzero[, 2]],
        Coefficient = round(pm[nonzero], 4),
        check.names = FALSE,
        row.names = NULL
      )
      print(path_tbl, row.names = FALSE, ...)
    }
  }

  invisible(x)
}

#' @export
summary.fimix_analysis <- function(object, ...) {
  result <- list(
    K                   = object$K,
    n_obs               = object$n_obs,
    converged           = object$converged,
    iterations          = object$iterations,
    n_starts            = object$n_starts_completed,
    segment_proportions = object$segment_proportions,
    segment_sizes       = object$segment_sizes,
    segment_paths       = object$segment_paths,
    segment_intercepts  = object$segment_intercepts,
    segment_variances   = object$segment_variances,
    info_criteria       = object$info_criteria,
    n_parameters        = object$n_parameters
  )

  class(result) <- c("summary.fimix_analysis", class(result))
  result
}

#' @export
print.summary.fimix_analysis <- function(x, ...) {
  cat("FIMIX-PLS Analysis Summary\n")
  cat("==========================\n")
  cat("Segments:", x$K, "\n")
  cat("Observations:", x$n_obs, "\n")
  cat("Free parameters:", x$n_parameters, "\n")
  cat("Converged:", if (x$converged) "Yes" else "No",
      "(", x$iterations, "iterations )\n")
  cat("Random starts:", x$n_starts, "\n")

  cat("\nSegment Proportions:\n")
  print(round(x$segment_proportions, 4), ...)

  cat("\nFit Criteria:\n")
  print(round(x$info_criteria, 4), ...)

  for (k in seq_len(x$K)) {
    cat("\n--- Segment", k, "---\n")

    cat("Intercepts:\n")
    print(round(x$segment_intercepts[[k]], 4), ...)

    cat("Path Coefficients:\n")
    pm <- x$segment_paths[[k]]
    nonzero <- which(pm != 0, arr.ind = TRUE)
    if (nrow(nonzero) > 0) {
      path_tbl <- data.frame(
        From = rownames(pm)[nonzero[, 1]],
        To = colnames(pm)[nonzero[, 2]],
        Coefficient = round(pm[nonzero], 4),
        check.names = FALSE,
        row.names = NULL
      )
      print(path_tbl, row.names = FALSE, ...)
    }

    cat("Residual Variances:\n")
    print(round(x$segment_variances[, k], 4), ...)
  }

  invisible(x)
}

#' Plot FIMIX-PLS Results
#'
#' @param x A \code{fimix_analysis} object from \code{assess_fimix()}.
#' @param type One of \code{"segments"} (segment proportion bar plot) or
#'   \code{"paths"} (grouped bar plot of path coefficients across segments).
#' @param ... Additional arguments passed to \code{plot()}.
#'
#' @export
plot.fimix_analysis <- function(x,
                                 type = c("segments", "paths"),
                                 ...) {
  type <- match.arg(type)

  if (type == "segments") {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    par(mar = c(5, 5, 4, 2))

    cols <- if (x$K <= 8) {
      palette.colors(x$K, palette = "Set2")
    } else {
      palette.colors(n = x$K)
    }

    barplot(x$segment_proportions * 100,
            names.arg = paste("Seg.", seq_len(x$K)),
            col = cols, ylim = c(0, 100),
            ylab = "Proportion (%)",
            main = paste("FIMIX-PLS: Segment Proportions (K =", x$K, ")"),
            ...)

  } else if (type == "paths") {
    # Collect all non-zero path labels and their segment-specific values
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

    cols <- if (x$K <= 8) {
      palette.colors(x$K, palette = "Set2")
    } else {
      palette.colors(n = x$K)
    }

    barplot(t(coef_matrix), beside = TRUE, col = cols,
            las = 2, ylab = "Path Coefficient",
            main = paste("FIMIX-PLS: Path Coefficients (K =", x$K, ")"),
            ...)

    legend("topright", legend = colnames(coef_matrix),
           fill = cols, cex = 0.8, bty = "n")
  }
}

# =============================================================================
# S3 METHODS — fimix_comparison
# =============================================================================

#' @export
print.fimix_comparison <- function(x, ...) {
  cat("FIMIX-PLS Segment Selection\n")
  cat("===========================\n")
  cat("K range:", min(x$K_range), "to", max(x$K_range), "\n\n")

  cat("Fit Criteria:\n")
  ft <- x$fit_table
  ft$K <- as.integer(ft$K)
  for (col_name in setdiff(names(ft), "K")) {
    if (col_name == "EN") {
      ft[[col_name]] <- round(ft[[col_name]], 4)
    } else {
      ft[[col_name]] <- round(ft[[col_name]], 2)
    }
  }
  print(ft, row.names = FALSE, ...)

  # Identify minimum for each criterion (lower = better, except EN: higher = better)
  cat("\nBest K by criterion:\n")
  for (crit in c("AIC", "AIC3", "AIC4", "BIC", "CAIC")) {
    vals <- ft[[crit]]
    if (all(is.na(vals))) next
    best_k <- ft$K[which.min(vals)]
    cat("  ", crit, ":", best_k, "\n")
  }
  en_vals <- ft[["EN"]]
  if (!all(is.na(en_vals))) {
    cat("   EN (best):", ft$K[which.max(en_vals)], "\n")
  }

  invisible(x)
}

#' @export
summary.fimix_comparison <- function(object, ...) {
  result <- list(
    fit_table = object$fit_table,
    K_range   = object$K_range,
    solutions = object$solutions
  )
  class(result) <- c("summary.fimix_comparison", class(result))
  result
}

#' @export
print.summary.fimix_comparison <- function(x, ...) {
  cat("FIMIX-PLS Comparison Summary\n")
  cat("============================\n\n")

  cat("Fit Table:\n")
  ft <- x$fit_table
  for (col_name in setdiff(names(ft), "K")) {
    if (col_name == "EN") {
      ft[[col_name]] <- round(ft[[col_name]], 4)
    } else {
      ft[[col_name]] <- round(ft[[col_name]], 2)
    }
  }
  print(ft, row.names = FALSE, ...)

  # Per-K summaries
  for (sol_name in names(x$solutions)) {
    sol <- x$solutions[[sol_name]]
    if (is.null(sol)) next

    cat("\n--- K =", sol$K, "---\n")
    cat("Converged:", if (sol$converged) "Yes" else "No", "\n")
    cat("Segment sizes:", paste(sol$segment_sizes, collapse = ", "), "\n")
    cat("EN:", round(sol$info_criteria["EN"], 4), "\n")
  }

  invisible(x)
}

#' Plot FIMIX-PLS Comparison
#'
#' @param x A \code{fimix_comparison} object from \code{assess_fimix_compare()}.
#' @param type One of \code{"criteria"} (line plot of IC vs K) or
#'   \code{"entropy"} (EN vs K).
#' @param criteria Character vector of criteria to plot when
#'   \code{type = "criteria"} (default: \code{c("AIC3", "AIC4", "BIC", "CAIC")}).
#' @param ... Additional arguments passed to \code{plot()}.
#'
#' @export
plot.fimix_comparison <- function(x,
                                    type = c("criteria", "entropy"),
                                    criteria = c("AIC3", "AIC4", "BIC", "CAIC"),
                                    ...) {
  type <- match.arg(type)
  ft <- x$fit_table

  if (type == "criteria") {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    par(mar = c(5, 5, 4, 8), xpd = TRUE)

    valid_criteria <- intersect(criteria, names(ft))
    if (length(valid_criteria) == 0) {
      stop("No valid criteria found to plot.", call. = FALSE)
    }

    vals <- ft[, valid_criteria, drop = FALSE]
    K_vals <- ft$K

    cols <- if (length(valid_criteria) <= 8) {
      palette.colors(length(valid_criteria), palette = "Set1")
    } else {
      seq_len(length(valid_criteria))
    }

    ylim <- range(vals, na.rm = TRUE)
    ylim <- ylim + c(-1, 1) * diff(ylim) * 0.05

    plot(K_vals, vals[[1]], type = "b", pch = 19, col = cols[1],
         xlim = range(K_vals), ylim = ylim,
         xlab = "Number of Segments (K)", ylab = "Information Criterion",
         main = "FIMIX-PLS: Model Selection Criteria",
         xaxt = "n", ...)
    axis(1, at = K_vals)

    if (length(valid_criteria) > 1) {
      for (i in 2:length(valid_criteria)) {
        lines(K_vals, vals[[i]], type = "b", pch = 19, col = cols[i])
      }
    }

    legend("topright", inset = c(-0.25, 0),
           legend = valid_criteria, col = cols, lty = 1, pch = 19,
           cex = 0.8, bty = "n")

  } else if (type == "entropy") {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    par(mar = c(5, 5, 4, 2))

    K_vals <- ft$K
    en_vals <- ft$EN

    plot(K_vals, en_vals, type = "b", pch = 19, col = "steelblue",
         xlim = range(K_vals), ylim = c(0, 1),
         xlab = "Number of Segments (K)", ylab = "Normed Entropy (EN)",
         main = "FIMIX-PLS: Classification Quality",
         xaxt = "n", ...)
    axis(1, at = K_vals)
    abline(h = 0.5, col = "gray70", lty = 2)
    text(max(K_vals), 0.5, "EN = 0.50", pos = 1, cex = 0.8, col = "gray50")
  }
}
