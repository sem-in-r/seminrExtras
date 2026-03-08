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
                            threshold = 1) {

  # Set seed for reproducibility of bootstrap resampling
  set.seed(seed)

  # ---------------------------------------------------------------------------
  # Step 1: Validate the model
  # ---------------------------------------------------------------------------
  if (!validate_seminr_model(seminr_model, "congruence_test")) {
    return(NULL)
  }

  # Get all construct names from the model
  construct_names <- colnames(seminr_model$construct_scores)

  # ---------------------------------------------------------------------------
  # Step 2: Define the congruence coefficient calculation
  # ---------------------------------------------------------------------------
  # The congruence coefficient (rc) measures pattern similarity between two
  # vectors. Formula: rc = sum(X*Y) / sqrt(sum(X^2) * sum(Y^2))
  # This is essentially a cosine similarity applied to correlation patterns.
  calc_congruence <- function(mat, X, Y) {
    return(sum(mat[, X] * mat[, Y]) / sqrt(sum(mat[, X]^2) * sum(mat[, Y]^2)))
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
    ret_mat <- stats::cor(it_model$construct_scores)

    # Replace diagonal with rhoC (composite reliability) values
    # This creates a matrix where diagonal = reliability, off-diagonal = correlations
    diag(ret_mat) <- seminr::rhoC_AVE(x = it_model)[colnames(ret_mat), 1]

    # Calculate congruence coefficient for each construct pair
    # Store in upper triangle of the result array
    ret_array[, , iter][upper.tri(ret_mat)] <- apply(combns, 1, function(x) {
      calc_congruence(ret_mat, x[1], x[2])
    })
  }

  # ---------------------------------------------------------------------------
  # Step 5: Calculate original (non-bootstrap) congruence coefficients
  # ---------------------------------------------------------------------------
  # Compute correlation matrix from original model
  cor_mat <- stats::cor(seminr_model$construct_scores)

  # Replace diagonal with rhoC values
  diag(cor_mat) <- seminr::rhoC_AVE(x = seminr_model)[colnames(ret_mat), 1]

  # Prepare matrix for original estimates (upper triangle only)
  original_matrix <- cor_mat
  original_matrix[lower.tri(original_matrix)] <- 0
  diag(original_matrix) <- 0

  # Calculate congruence coefficients for original data
  original_matrix[upper.tri(original_matrix)] <- apply(combns, 1, function(x) {
    calc_congruence(cor_mat, x[1], x[2])
  })

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
        ci <- seminr:::conf_int(boot_array,
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
  return_matrix <- seminr:::convert_to_table_output(return_matrix)

  return(list(results = return_matrix))
}
