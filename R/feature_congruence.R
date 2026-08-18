# =============================================================================
# feature_congruence.R - Congruence Coefficient Testing
# =============================================================================
# This file implements bootstrap-based congruence coefficient testing for
# assessing measure congruence in nomological networks, as described in:
#
# Franke, G. R., Sarstedt, M., & Danks, N. P. (2021). "Assessing measure
# congruence in nomological networks." Journal of Business Research, 130, 318-334.
#
# The congruence coefficient measures the similarity between two constructs
# in terms of their correlational patterns. A coefficient close to 1 indicates
# high congruence (the constructs behave similarly in the nomological network).
# =============================================================================

# Standardised Cronbach's alpha per construct.
#
# Mirrors seminr's internal cronbachs_alpha(): the correlation matrix of a
# construct's items, alpha = (k/(k-1))(1 - sum(diag)/sum(all)), and 1 for any
# single-indicator construct. Reimplemented here rather than taken from
# summary()$reliability because summary() is roughly 400x more expensive and
# this runs once per bootstrap resample -- ~10 minutes of overhead at the
# default nboot = 2000. Uses only exported seminr functions.
#
# @param model A fitted seminr model.
# @param constructs Character vector of construct names.
# @return Named numeric vector of alphas, in the order of `constructs`.
# @noRd
construct_alphas <- function(model, constructs) {
  item_cors <- stats::cor(model$data)
  vapply(constructs, function(cn) {
    items <- seminr::construct_items(model$mmMatrix, cn)
    if (length(items) < 2) {
      return(1)
    }
    cm <- item_cors[items, items, drop = FALSE]
    k <- nrow(cm)
    (k / (k - 1)) * (1 - sum(diag(cm)) / sum(cm))
  }, numeric(1))
}

#' Bootstrap congruence coefficient test
#'
#' `congruence_test` conducts a bootstrapped significance test of congruence
#' coefficients between all pairs of constructs in a PLS-SEM model.
#'
#' The congruence coefficient (rc) measures how similarly two constructs
#' relate to other constructs in the model. Values close to 1 indicate high
#' congruence. The test evaluates H0: rc < threshold (default = 1).
#'
#' @param seminr_model The SEMinR model for congruence analysis
#' @param nboot The number of bootstrap subsamples to execute (defaults to 2000).
#' @param seed The seed for reproducibility (defaults to 123).
#' @param alpha The required level of alpha for statistical testing (defaults
#'   to 0.05). Used to compute confidence intervals.
#' @param threshold The threshold with which to compare significance testing.
#'   H0: rc < threshold (defaults to 1).
#' @param reliability Which reliability estimate to place on the diagonal of the
#'   construct-correlation matrix: `"rhoA"` (default, matches SmartPLS),
#'   `"rhoC"` (composite reliability, the behaviour of seminrExtras <= 1.0.2),
#'   `"cronbach"` (Cronbach's alpha), or `"one"` (unity, as permitted by Franke et
#'   al. (2021) when reliabilities are unknown). Not to be confused with the
#'   `alpha` argument above, which sets the confidence level. The option is named
#'   `"cronbach"` rather than `"alpha"` precisely to keep the two apart.
#'
#'   Franke et al. (2021, Eq. 2) specify "the reliabilities" without fixing an
#'   estimator, so all four are in specification. They agree wherever a
#'   construct has a well-defined internal consistency and diverge where it does
#'   not:
#'
#'   \itemize{
#'     \item **Single-item constructs** get 1 under `"rhoA"`, `"rhoC"` and
#'       `"cronbach"` alike. In a model built entirely from single indicators the
#'       argument has no effect at all.
#'     \item **Mode B (formative) constructs** are where they part company.
#'       `"rhoA"` returns exactly 1, because internal consistency is undefined
#'       for a composite and nothing is estimated. `"rhoC"` and `"cronbach"` both
#'       still compute a value from the indicators. That is arguably the less
#'       honest choice — each presumes a measurement model the construct does
#'       not have — but it does keep one rule for every multi-item construct.
#'   }
#'
#'   Note that this only ever affects pairs that **involve** a Mode B construct.
#'   Each column of the matrix carries exactly one reliability -- its own
#'   construct's -- so a coefficient between two reflective constructs is
#'   invariant to whatever sits on any other construct's diagonal.
#'
#'   `"cronbach"` is offered chiefly for comparison with covariance-based SEM,
#'   where rho_A is not available.
#'
#' @section Model types:
#' **Interaction constructs are excluded** from the construct set — entirely,
#' not merely from the pair list, since Eq. 2 sums over the whole set. An
#' interaction term's measurement is fixed by the product method rather than by
#' theory, so its position in a nomological network is not interpretable. A
#' message names any construct dropped. This matches [assess_cta()],
#' [assess_pos()], [assess_pcm()] and [assess_cipma()].
#'
#' **Higher-order models are not supported** and are refused with a warning. Two-
#' stage estimation replaces the lower-order constructs with a single
#' higher-order composite, and it has not been established what belongs on that
#' composite's diagonal, nor whether a congruence coefficient between a
#' higher-order and a first-order construct is interpretable.
#'
#' @return A list containing a matrix of congruence coefficients and
#'   significance test results for all construct pairs.
#'
#' @seealso [assess_cvpat()] and [assess_cvpat_compare()] for predictive validity testing
#'
#' @references Franke, G. R., Sarstedt, M., & Danks, N. P. (2021). Assessing
#' measure congruence in nomological networks. Journal of Business Research,
#' 130, 318-334.
#'
#' @examples
#' # Load libraries
#' library(seminr)
#' library(seminrExtras)
#'
#' # Create measurement model ----
#' corp_rep_mm <- constructs(
#'   composite("QUAL", multi_items("qual_", 1:8), weights = mode_B),
#'   composite("PERF", multi_items("perf_", 1:5), weights = mode_B),
#'   composite("CSOR", multi_items("csor_", 1:5), weights = mode_B),
#'   composite("ATTR", multi_items("attr_", 1:3), weights = mode_B),
#'   composite("COMP", multi_items("comp_", 1:3)),
#'   composite("LIKE", multi_items("like_", 1:3))
#' )
#'
#' # Create structural model ----
#' corp_rep_sm <- relationships(
#'   paths(from = c("QUAL", "PERF", "CSOR", "ATTR"), to = c("COMP", "LIKE"))
#' )
#'
#' # Estimate the model ----
#' corp_rep_pls_model <- estimate_pls(
#'   data = corp_rep_data,
#'   measurement_model = corp_rep_mm,
#'   structural_model  = corp_rep_sm,
#'   missing = mean_replacement,
#'  missing_value = "-99")
#'
#' # Assess the base model ----
#' congruence_test(seminr_model = corp_rep_pls_model,
#'                 nboot = 20,
#'                 seed = 123,
#'                 alpha = 0.05,
#'                 threshold = 1)
#'
#' @export
congruence_test <- function(seminr_model,
                            nboot = 2000,
                            seed = 123,
                            alpha = 0.05,
                            threshold = 1,
                            reliability = c("rhoA", "rhoC", "cronbach", "one")) {

  reliability <- match.arg(reliability)

  # Set seed for reproducibility of bootstrap resampling
  set.seed(seed)

  # ---------------------------------------------------------------------------
  # Step 1: Validate the model
  # ---------------------------------------------------------------------------
  if (!validate_seminr_model(seminr_model, "congruence_test")) {
    return(NULL)
  }

  # ---------------------------------------------------------------------------
  # Step 1b: Refuse model types that are not yet supported
  # ---------------------------------------------------------------------------
  # Higher-order models are out of scope for now. Two-stage estimation replaces
  # the lower-order constructs with a single higher-order composite, and it has
  # not been established what belongs on the diagonal for that composite or
  # whether a congruence coefficient between a HOC and a first-order construct
  # is interpretable. Refuse rather than return an unvalidated number.
  if (!is.null(seminr_model$hoc)) {
    warning("congruence_test() does not yet support higher-order models.",
            call. = FALSE)
    return(NULL)
  }

  # Get all construct names from the model
  construct_names <- colnames(seminr_model$construct_scores)

  # ---------------------------------------------------------------------------
  # Step 1c: Drop interaction constructs
  # ---------------------------------------------------------------------------
  # An interaction term's measurement is determined by the product method, not
  # by theory, so it is not a construct whose position in the nomological
  # network can be interpreted. It is removed from the construct set entirely,
  # not merely from the pair list: Eq. 2 sums over the whole set, so leaving it
  # in the correlation vectors would still perturb every coefficient. This
  # matches assess_cta(), assess_pos(), assess_pcm() and assess_cipma().
  is_interaction <- grepl("*", construct_names, fixed = TRUE)
  if (any(is_interaction)) {
    message("Excluding interaction constructs (measurement determined by method): ",
            paste(construct_names[is_interaction], collapse = ", "))
    construct_names <- construct_names[!is_interaction]
  }

  if (length(construct_names) < 2) {
    warning("congruence_test() needs at least two non-interaction constructs.",
            call. = FALSE)
    return(NULL)
  }

  # ---------------------------------------------------------------------------
  # Step 2: Define the congruence coefficient calculation
  # ---------------------------------------------------------------------------
  # The congruence coefficient (rc) measures pattern similarity between two
  # vectors. Formula: rc = sum(X*Y) / sqrt(sum(X^2) * sum(Y^2))
  # This is essentially a cosine similarity applied to correlation patterns.
  calc_congruence <- function(mat, X, Y) {
    return(sum(mat[, X] * mat[, Y]) / sqrt(sum(mat[, X]^2) * sum(mat[, Y]^2)))
  }

  # Reliabilities for the diagonal, per the selected estimator. Recomputed from
  # whichever model is passed in, so the bootstrap gets resample-specific values
  # rather than the original fit's.
  diagonal_values <- function(model, constructs) {
    switch(reliability,
      rhoA  = seminr::rho_A(model, constructs)[constructs, 1],
      rhoC  = seminr::rhoC_AVE(x = model)[constructs, 1],
      cronbach = construct_alphas(model, constructs),
      one   = stats::setNames(rep(1, length(constructs)), constructs)
    )
  }

  # ---------------------------------------------------------------------------
  # Step 3: Generate all pairwise construct combinations
  # ---------------------------------------------------------------------------
  # We test congruence for every unique pair of constructs
  combns <- t(utils::combn(construct_names, 2))

  # ---------------------------------------------------------------------------
  # Step 4: Bootstrap loop - resample and compute congruence coefficients
  # ---------------------------------------------------------------------------
  # Create 3D array to store bootstrap results: [constructs x constructs x iterations]
  ret_array <- array(NA,
                     dim = list(length(construct_names), length(construct_names), nboot),
                     dimnames = list(construct_names, construct_names, 1:nboot))

  for (iter in seq_len(nboot)) {
    # Resample data with replacement and re-estimate the model
    resampled_data <- seminr_model$rawdata[
      sample(nrow(seminr_model$rawdata), nrow(seminr_model$rawdata), replace = TRUE),
    ]
    it_model <- suppressMessages(seminr::rerun(seminr_model, data = resampled_data))

    # Compute correlation matrix of construct scores for this bootstrap sample
    # Restricted to construct_names, so any excluded interaction leaves the
    # summation in Eq. 2 rather than merely the pair list.
    ret_mat <- stats::cor(it_model$construct_scores[, construct_names, drop = FALSE])

    # Replace diagonal with the selected reliability estimates
    # This creates a matrix where diagonal = reliability, off-diagonal = correlations
    diag(ret_mat) <- diagonal_values(it_model, colnames(ret_mat))

    # Calculate congruence coefficient for each construct pair.
    # Assign by construct name (not via upper.tri()): apply() returns values in
    # combn() row-major order, whereas `M[upper.tri(M)] <-` fills column-major.
    # Those orderings diverge for >3 constructs, which would attach each
    # coefficient to the wrong construct pair.
    slice <- matrix(NA, length(construct_names), length(construct_names),
                    dimnames = list(construct_names, construct_names))
    for (r in seq_len(nrow(combns))) {
      slice[combns[r, 1], combns[r, 2]] <- calc_congruence(ret_mat, combns[r, 1], combns[r, 2])
    }
    ret_array[, , iter] <- slice
  }

  # ---------------------------------------------------------------------------
  # Step 5: Calculate original (non-bootstrap) congruence coefficients
  # ---------------------------------------------------------------------------
  # Compute correlation matrix from original model
  cor_mat <- stats::cor(seminr_model$construct_scores[, construct_names, drop = FALSE])

  # Replace diagonal with the selected reliability estimates
  diag(cor_mat) <- diagonal_values(seminr_model, colnames(cor_mat))

  # Prepare matrix for original estimates (upper triangle only)
  original_matrix <- cor_mat
  original_matrix[lower.tri(original_matrix)] <- 0
  diag(original_matrix) <- 0

  # Calculate congruence coefficients for original data.
  # Assign by construct name (see note in Step 4) so each coefficient lands in
  # the cell for its actual construct pair rather than the column-major slot.
  for (r in seq_len(nrow(combns))) {
    original_matrix[combns[r, 1], combns[r, 2]] <- calc_congruence(cor_mat, combns[r, 1], combns[r, 2])
  }

  # ---------------------------------------------------------------------------
  # Step 6: Compute bootstrap statistics for each construct pair
  # ---------------------------------------------------------------------------
  boot_array <- ret_array

  # Initialize result vectors
  Path <- c()        # Construct pair labels (e.g., "QUAL -> PERF")
  original <- c()    # Original congruence coefficient
  boot_mean <- c()   # Difference from threshold
  boot_SD <- c()     # Bootstrap standard deviation
  t_stat <- c()      # T-statistic for significance test
  lower <- c()       # Lower confidence interval bound
  upper <- c()       # Upper confidence interval bound

  # Calculate alpha/2 for confidence interval labels
  alpha_text <- alpha / 2 * 100

  # Use upper.tri mask to iterate over construct pairs
  ut_mask <- upper.tri(original_matrix)

  # Loop through upper triangle to extract results for each pair
  for (i in seq_len(nrow(original_matrix))) {
    for (j in seq_len(ncol(original_matrix))) {
      if (ut_mask[i, j]) {
        # Store construct pair label
        Path <- append(Path, paste(rownames(original_matrix)[i], " -> ",
                                   colnames(original_matrix)[j]))

        # Store original coefficient
        original <- append(original, original_matrix[i, j])

        # Compute difference from threshold (for significance test)
        # Diff = threshold - |original| (positive if below threshold)
        boot_mean <- append(boot_mean, (threshold - abs(original_matrix[i, j])))

        # Bootstrap standard deviation of congruence coefficient
        boot_SD <- append(boot_SD, stats::sd(boot_array[i, j, ]))

        # Compute t-statistic: (threshold - |rc|) / SE
        # Guard against division by near-zero SD (indicates perfect stability)
        if (stats::sd(boot_array[i, j, ]) < .Machine$double.eps) {
          t_stat <- append(t_stat, NA)
        } else {
          t_stat <- append(t_stat, (threshold - abs(original_matrix[i, j])) /
                             stats::sd(boot_array[i, j, ]))
        }

        # Compute bootstrap confidence intervals using seminr's internal function
        ci <- conf_int(boot_array,
                                from = rownames(original_matrix)[i],
                                to = colnames(original_matrix)[j],
                                alpha = alpha)
        lower <- append(lower, ci[[1]])
        upper <- append(upper, ci[[2]])
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Step 7: Compile and format results
  # ---------------------------------------------------------------------------
  return_matrix <- cbind(original, boot_mean, boot_SD, t_stat, lower, upper)

  colnames(return_matrix) <- c("Original Est.",        # Congruence coefficient
                               "Diff",                  # threshold - |rc|
                               "Bootstrap SD",          # Standard error
                               "T Stat.",               # t-statistic
                               paste(alpha_text, "% CI", sep = ""),   # Lower CI
                               paste((100 - alpha_text), "% CI", sep = ""))  # Upper CI
  rownames(return_matrix) <- Path

  # Convert to table_output class for consistent printing
  return_matrix <- convert_to_table_output(return_matrix)

  return(list(results = return_matrix))
}
