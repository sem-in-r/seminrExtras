# =============================================================================
# Tests for FIMIX-PLS (Finite Mixture PLS)
# =============================================================================

# --- Setup: pre-compute models and results -----------------------------------
library(seminr)

mobi_mm <- constructs(
  composite("Image",        multi_items("IMAG", 1:5)),
  composite("Expectation",  multi_items("CUEX", 1:3)),
  composite("Value",        multi_items("PERV", 1:2)),
  composite("Satisfaction", multi_items("CUSA", 1:3)),
  composite("Loyalty",      multi_items("CUSL", 1:3))
)
mobi_sm <- relationships(
  paths(from = "Image",       to = c("Expectation", "Satisfaction", "Loyalty")),
  paths(from = "Expectation", to = c("Value", "Satisfaction")),
  paths(from = "Value",       to = "Satisfaction"),
  paths(from = "Satisfaction", to = "Loyalty")
)
pls_model <- estimate_pls(data = mobi, measurement_model = mobi_mm,
                           structural_model = mobi_sm)

# Pre-compute FIMIX for K=2 (reused across tests)
fimix_k2 <- assess_fimix(pls_model, K = 2, nstart = 5, seed = 123)

# =============================================================================
# STRUCTURE TESTS
# =============================================================================

test_that("assess_fimix returns fimix_analysis object", {
  expect_s3_class(fimix_k2, "fimix_analysis")
})

test_that("fimix result contains all expected components", {
  expected <- c("K", "segment_proportions", "segment_sizes", "posterior",
                "segment_assignment", "segment_paths", "segment_intercepts",
                "segment_variances", "log_likelihood", "n_parameters",
                "info_criteria", "converged", "iterations",
                "n_starts_completed", "pls_model", "n_obs")
  expect_true(all(expected %in% names(fimix_k2)))
})

test_that("K is stored correctly", {
  expect_equal(fimix_k2$K, 2)
})

test_that("n_obs matches model data", {
  expect_equal(fimix_k2$n_obs, nrow(pls_model$construct_scores))
})

# =============================================================================
# SEGMENT PROPORTION AND ASSIGNMENT TESTS
# =============================================================================

test_that("segment proportions sum to 1", {
  expect_equal(sum(fimix_k2$segment_proportions), 1, tolerance = 1e-10)
})

test_that("segment proportions are all positive", {
  expect_true(all(fimix_k2$segment_proportions > 0))
})

test_that("segment sizes sum to N", {
  expect_equal(sum(fimix_k2$segment_sizes), fimix_k2$n_obs)
})

test_that("segment assignments are integers from 1 to K", {
  expect_true(all(fimix_k2$segment_assignment %in% 1:2))
})

test_that("segment assignment length matches N", {
  expect_equal(length(fimix_k2$segment_assignment), fimix_k2$n_obs)
})

# =============================================================================
# POSTERIOR PROBABILITY TESTS
# =============================================================================

test_that("posterior matrix has correct dimensions", {
  expect_equal(dim(fimix_k2$posterior), c(fimix_k2$n_obs, 2))
})

test_that("posterior rows sum to 1", {
  row_sums <- rowSums(fimix_k2$posterior)
  expect_equal(row_sums, rep(1, fimix_k2$n_obs), tolerance = 1e-10)
})

test_that("posterior values are in [0, 1]", {
  expect_true(all(fimix_k2$posterior >= 0))
  expect_true(all(fimix_k2$posterior <= 1))
})

test_that("hard assignment matches argmax of posterior", {
  argmax <- apply(fimix_k2$posterior, 1, which.max)
  expect_equal(fimix_k2$segment_assignment, argmax)
})

# =============================================================================
# PATH COEFFICIENT TESTS
# =============================================================================

test_that("segment_paths is a list of length K", {
  expect_true(is.list(fimix_k2$segment_paths))
  expect_equal(length(fimix_k2$segment_paths), 2)
})

test_that("segment path matrices have same dimensions as model path_coef", {
  for (k in 1:2) {
    expect_equal(dim(fimix_k2$segment_paths[[k]]),
                 dim(pls_model$path_coef))
    expect_equal(rownames(fimix_k2$segment_paths[[k]]),
                 rownames(pls_model$path_coef))
    expect_equal(colnames(fimix_k2$segment_paths[[k]]),
                 colnames(pls_model$path_coef))
  }
})

test_that("segment path coefficients are non-zero where model has paths", {
  nonzero_mask <- pls_model$path_coef != 0
  for (k in 1:2) {
    # All structural paths in the model should have coefficients
    seg_coefs <- fimix_k2$segment_paths[[k]][nonzero_mask]
    expect_true(all(seg_coefs != 0))
  }
})

test_that("segment path coefficients are zero where model has no paths", {
  zero_mask <- pls_model$path_coef == 0
  for (k in 1:2) {
    seg_zeros <- fimix_k2$segment_paths[[k]][zero_mask]
    expect_true(all(seg_zeros == 0))
  }
})

test_that("segment intercepts have correct structure", {
  expect_equal(length(fimix_k2$segment_intercepts), 2)
  endogenous <- seminr:::all_endogenous(pls_model$smMatrix)
  for (k in 1:2) {
    expect_equal(length(fimix_k2$segment_intercepts[[k]]), length(endogenous))
    expect_true(all(names(fimix_k2$segment_intercepts[[k]]) %in% endogenous))
  }
})

# =============================================================================
# SEGMENT VARIANCE TESTS
# =============================================================================

test_that("segment_variances has correct dimensions", {
  endogenous <- seminr:::all_endogenous(pls_model$smMatrix)
  expect_equal(nrow(fimix_k2$segment_variances), length(endogenous))
  expect_equal(ncol(fimix_k2$segment_variances), 2)
})

test_that("segment variances are positive", {
  expect_true(all(fimix_k2$segment_variances > 0))
})

# =============================================================================
# INFORMATION CRITERIA TESTS
# =============================================================================

test_that("info_criteria contains expected names", {
  expected <- c("lnL", "AIC", "AIC3", "AIC4", "BIC", "CAIC", "HQ", "MDL5", "EN")
  expect_true(all(expected %in% names(fimix_k2$info_criteria)))
})

test_that("all criteria are finite", {
  expect_true(all(is.finite(fimix_k2$info_criteria)))
})

test_that("AIC < AIC3 < AIC4 (by formula)", {
  ic <- fimix_k2$info_criteria
  expect_true(ic["AIC"] < ic["AIC3"])
  expect_true(ic["AIC3"] < ic["AIC4"])
})

test_that("entropy is in [0, 1]", {
  en <- fimix_k2$info_criteria["EN"]
  expect_true(en >= 0 && en <= 1)
})

test_that("log-likelihood is negative", {
  expect_true(fimix_k2$log_likelihood < 0)
})

test_that("n_parameters is positive integer", {
  expect_true(fimix_k2$n_parameters > 0)
  expect_equal(fimix_k2$n_parameters, as.integer(fimix_k2$n_parameters))
})

# =============================================================================
# CONVERGENCE TESTS
# =============================================================================

test_that("EM converges for K=2", {
  expect_true(fimix_k2$converged)
})

test_that("iterations is positive", {
  expect_true(fimix_k2$iterations > 0)
})

test_that("n_starts_completed equals nstart", {
  expect_equal(fimix_k2$n_starts_completed, 5L)
})

# =============================================================================
# REPRODUCIBILITY TESTS
# =============================================================================

test_that("results are deterministic with same seed", {
  r1 <- assess_fimix(pls_model, K = 2, nstart = 3, seed = 42)
  r2 <- assess_fimix(pls_model, K = 2, nstart = 3, seed = 42)
  expect_equal(r1$log_likelihood, r2$log_likelihood)
  expect_equal(r1$segment_proportions, r2$segment_proportions)
  expect_equal(r1$segment_assignment, r2$segment_assignment)
})

# =============================================================================
# INPUT VALIDATION TESTS
# =============================================================================

test_that("non-seminr model returns NULL with warning", {
  expect_warning(result <- assess_fimix("not a model", K = 2),
                 "SEMinR models")
  expect_null(result)
})

test_that("K < 2 raises error", {
  expect_error(assess_fimix(pls_model, K = 1), "K must be an integer >= 2")
})

test_that("K = 0 raises error", {
  expect_error(assess_fimix(pls_model, K = 0), "K must be an integer >= 2")
})

test_that("non-integer K raises error", {
  expect_error(assess_fimix(pls_model, K = 2.5), "K must be an integer >= 2")
})

test_that("nstart < 1 raises error", {
  expect_error(assess_fimix(pls_model, K = 2, nstart = 0),
               "nstart must be an integer >= 1")
})

# =============================================================================
# S3 METHOD TESTS
# =============================================================================

test_that("print.fimix_analysis returns invisibly", {
  expect_invisible(print(fimix_k2))
})

test_that("print output contains key information", {
  out <- capture.output(print(fimix_k2))
  out_text <- paste(out, collapse = "\n")

  expect_true(grepl("FIMIX-PLS Analysis", out_text))
  expect_true(grepl("Segments: 2", out_text))
  expect_true(grepl("Converged: Yes", out_text))
  expect_true(grepl("Segment Proportions", out_text))
  expect_true(grepl("Fit Criteria", out_text))
  expect_true(grepl("Segment Path Coefficients", out_text))
})

test_that("summary returns summary.fimix_analysis object", {
  s <- summary(fimix_k2)
  expect_s3_class(s, "summary.fimix_analysis")
})

test_that("summary print includes segment details", {
  out <- capture.output(print(summary(fimix_k2)))
  out_text <- paste(out, collapse = "\n")

  expect_true(grepl("Summary", out_text))
  expect_true(grepl("Segment 1", out_text))
  expect_true(grepl("Intercept", out_text))
  expect_true(grepl("Residual Variance", out_text))
})

test_that("plot.fimix_analysis runs without error", {
  expect_no_error(plot(fimix_k2, type = "segments"))
  expect_no_error(plot(fimix_k2, type = "paths"))
})

# =============================================================================
# COMPARISON TESTS
# =============================================================================

test_that("assess_fimix_compare returns fimix_comparison object", {
  comp <- assess_fimix_compare(pls_model, K_range = 2:3,
                                 nstart = 3, seed = 123)
  expect_s3_class(comp, "fimix_comparison")
})

test_that("comparison contains solutions for each K", {
  comp <- assess_fimix_compare(pls_model, K_range = 2:3,
                                 nstart = 3, seed = 123)
  expect_true("K2" %in% names(comp$solutions))
  expect_true("K3" %in% names(comp$solutions))
})

test_that("comparison fit_table has correct structure", {
  comp <- assess_fimix_compare(pls_model, K_range = 2:3,
                                 nstart = 3, seed = 123)
  expect_true("K" %in% names(comp$fit_table))
  expect_true("AIC" %in% names(comp$fit_table))
  expect_true("BIC" %in% names(comp$fit_table))
  expect_true("EN" %in% names(comp$fit_table))
  expect_equal(nrow(comp$fit_table), 2)
})

test_that("comparison print shows best K by criterion", {
  comp <- assess_fimix_compare(pls_model, K_range = 2:3,
                                 nstart = 3, seed = 123)
  out <- capture.output(print(comp))
  out_text <- paste(out, collapse = "\n")

  expect_true(grepl("FIMIX-PLS Segment Selection", out_text))
  expect_true(grepl("Best K by criterion", out_text))
})

test_that("comparison plot runs without error", {
  comp <- assess_fimix_compare(pls_model, K_range = 2:3,
                                 nstart = 3, seed = 123)
  expect_no_error(plot(comp, type = "criteria"))
  expect_no_error(plot(comp, type = "entropy"))
})

test_that("comparison K_range validation works", {
  expect_error(assess_fimix_compare(pls_model, K_range = 1),
               "K_range must be a vector of integers >= 2")
})

# =============================================================================
# MODEL TYPE TESTS
# =============================================================================

test_that("FIMIX works with mediation model", {
  mm_med <- constructs(
    composite("Image",        multi_items("IMAG", 1:5)),
    composite("Expectation",  multi_items("CUEX", 1:3)),
    composite("Satisfaction", multi_items("CUSA", 1:3)),
    composite("Loyalty",      multi_items("CUSL", 1:3))
  )
  sm_med <- relationships(
    paths(from = "Image", to = c("Expectation", "Loyalty")),
    paths(from = "Expectation", to = "Satisfaction"),
    paths(from = "Satisfaction", to = "Loyalty")
  )
  pls_med <- estimate_pls(data = mobi, measurement_model = mm_med,
                            structural_model = sm_med)

  result <- assess_fimix(pls_med, K = 2, nstart = 3, seed = 42)

  expect_s3_class(result, "fimix_analysis")
  expect_true(result$converged)
  expect_equal(sum(result$segment_sizes), nrow(mobi))
})

test_that("FIMIX works with moderation model (two-stage)", {
  mm_mod <- constructs(
    composite("Image",        multi_items("IMAG", 1:5)),
    composite("Value",        multi_items("PERV", 1:2)),
    composite("Satisfaction", multi_items("CUSA", 1:3)),
    composite("Loyalty",      multi_items("CUSL", 1:3)),
    interaction_term(iv = "Image", moderator = "Value", method = two_stage)
  )
  sm_mod <- relationships(
    paths(from = c("Image", "Value", "Image*Value"), to = "Satisfaction"),
    paths(from = "Satisfaction", to = "Loyalty")
  )
  pls_mod <- estimate_pls(data = mobi, measurement_model = mm_mod,
                            structural_model = sm_mod)

  result <- assess_fimix(pls_mod, K = 2, nstart = 3, seed = 42)

  expect_s3_class(result, "fimix_analysis")
  expect_true(result$converged)
  # Interaction construct should be a predictor in segment paths
  expect_true(any(result$segment_paths[[1]]["Image*Value", ] != 0))
})

test_that("FIMIX works with HOC two-stage model", {
  mm_hoc <- constructs(
    composite("Image",        multi_items("IMAG", 1:5)),
    composite("Expectation",  multi_items("CUEX", 1:3)),
    higher_composite("Quality", c("Image", "Expectation"), method = two_stage),
    composite("Satisfaction", multi_items("CUSA", 1:3)),
    composite("Loyalty",      multi_items("CUSL", 1:3))
  )
  sm_hoc <- relationships(
    paths(from = "Quality", to = c("Satisfaction", "Loyalty")),
    paths(from = "Satisfaction", to = "Loyalty")
  )
  pls_hoc <- estimate_pls(data = mobi, measurement_model = mm_hoc,
                            structural_model = sm_hoc)

  result <- assess_fimix(pls_hoc, K = 2, nstart = 3, seed = 42)

  expect_s3_class(result, "fimix_analysis")
  expect_true(result$converged)
  expect_true("Quality" %in% rownames(result$segment_paths[[1]]))
})

# =============================================================================
# SEGMENT DIFFERENTIATION TESTS
# =============================================================================

test_that("segments have different path coefficients", {
  # At least some paths should differ between segments
  diffs <- abs(fimix_k2$segment_paths[[1]] - fimix_k2$segment_paths[[2]])
  nonzero_mask <- pls_model$path_coef != 0
  max_diff <- max(diffs[nonzero_mask])
  expect_true(max_diff > 0.01)
})
