# =============================================================================
# Tests for PLS-POS (Prediction-Oriented Segmentation)
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

# Pre-compute PLS-POS for K=2 (reused across tests)
pos_k2 <- assess_pos(pls_model, K = 2, nstart = 3, max_iter = 20, seed = 123)

# =============================================================================
# STRUCTURE TESTS
# =============================================================================

test_that("assess_pos returns pos_analysis object", {
  expect_s3_class(pos_k2, "pos_analysis")
})

test_that("pos result contains all expected components", {
  expected <- c("K", "segment_assignment", "segment_sizes", "segment_models",
                "segment_rsquared", "segment_paths", "objective",
                "converged", "iterations", "nstart", "all_objectives",
                "endogenous", "pls_model", "n_obs")
  expect_true(all(expected %in% names(pos_k2)))
})

test_that("K is stored correctly", {
  expect_equal(pos_k2$K, 2)
})

test_that("n_obs matches model data", {
  expect_equal(pos_k2$n_obs, nrow(pls_model$construct_scores))
})

# =============================================================================
# SEGMENT ASSIGNMENT AND SIZE TESTS
# =============================================================================

test_that("segment assignments are integers from 1 to K", {
  expect_true(all(pos_k2$segment_assignment %in% 1:2))
})

test_that("segment assignment length matches N", {
  expect_equal(length(pos_k2$segment_assignment), pos_k2$n_obs)
})

test_that("segment sizes sum to N", {
  expect_equal(sum(pos_k2$segment_sizes), pos_k2$n_obs)
})

test_that("segment sizes match assignment counts", {
  counts <- tabulate(pos_k2$segment_assignment, nbins = 2)
  expect_equal(as.integer(pos_k2$segment_sizes), counts)
})

test_that("all segments are non-empty", {
  expect_true(all(pos_k2$segment_sizes > 0))
})

# =============================================================================
# OBJECTIVE CRITERION TESTS
# =============================================================================

test_that("objective is finite and positive", {
  expect_true(is.finite(pos_k2$objective))
  expect_true(pos_k2$objective > 0)
})

test_that("objective equals sum of segment R-squared values", {
  expected_obj <- sum(pos_k2$segment_rsquared, na.rm = TRUE)
  expect_equal(pos_k2$objective, expected_obj, tolerance = 0.01)
})

test_that("all_objectives has length nstart", {
  expect_equal(length(pos_k2$all_objectives), pos_k2$nstart)
})

test_that("best objective is max of all_objectives", {
  valid <- pos_k2$all_objectives[is.finite(pos_k2$all_objectives)]
  expect_equal(pos_k2$objective, max(valid))
})

# =============================================================================
# SEGMENT MODEL TESTS
# =============================================================================

test_that("segment_models is a list of length K", {
  expect_true(is.list(pos_k2$segment_models))
  expect_equal(length(pos_k2$segment_models), 2)
})

test_that("segment models are seminr model objects", {
  for (k in 1:2) {
    expect_true(inherits(pos_k2$segment_models[[k]], "seminr_model"))
  }
})

test_that("segment path matrices have same structure as global model", {
  for (k in 1:2) {
    expect_equal(dim(pos_k2$segment_paths[[k]]),
                 dim(pls_model$path_coef))
    expect_equal(rownames(pos_k2$segment_paths[[k]]),
                 rownames(pls_model$path_coef))
    expect_equal(colnames(pos_k2$segment_paths[[k]]),
                 colnames(pls_model$path_coef))
  }
})

test_that("segment path coefficients are non-zero where model has paths", {
  nonzero_mask <- pls_model$path_coef != 0
  for (k in 1:2) {
    seg_coefs <- pos_k2$segment_paths[[k]][nonzero_mask]
    expect_true(all(seg_coefs != 0))
  }
})

# =============================================================================
# R-SQUARED TESTS
# =============================================================================

test_that("segment_rsquared has correct dimensions", {
  endogenous <- pos_k2$endogenous
  expect_equal(nrow(pos_k2$segment_rsquared), length(endogenous))
  expect_equal(ncol(pos_k2$segment_rsquared), 2)
})

test_that("R-squared values are in [0, 1]", {
  expect_true(all(pos_k2$segment_rsquared >= 0, na.rm = TRUE))
  expect_true(all(pos_k2$segment_rsquared <= 1, na.rm = TRUE))
})

test_that("endogenous constructs exclude interaction terms", {
  expect_false(any(grepl("\\*", pos_k2$endogenous)))
})

# =============================================================================
# CONVERGENCE TESTS
# =============================================================================

test_that("iterations is a positive integer", {
  expect_true(pos_k2$iterations > 0)
})

test_that("nstart is stored correctly", {
  expect_equal(pos_k2$nstart, 3L)
})

# =============================================================================
# REPRODUCIBILITY TESTS
# =============================================================================

test_that("results are deterministic with same seed", {
  r1 <- assess_pos(pls_model, K = 2, nstart = 2, max_iter = 10, seed = 42)
  r2 <- assess_pos(pls_model, K = 2, nstart = 2, max_iter = 10, seed = 42)
  expect_equal(r1$objective, r2$objective)
  expect_equal(r1$segment_assignment, r2$segment_assignment)
  expect_equal(r1$segment_sizes, r2$segment_sizes)
})

# =============================================================================
# SEGMENT DIFFERENTIATION TESTS
# =============================================================================

test_that("segments have different path coefficients", {
  diffs <- abs(pos_k2$segment_paths[[1]] - pos_k2$segment_paths[[2]])
  nonzero_mask <- pls_model$path_coef != 0
  max_diff <- max(diffs[nonzero_mask])
  expect_true(max_diff > 0.01)
})

# =============================================================================
# INPUT VALIDATION TESTS
# =============================================================================

test_that("non-seminr model returns NULL with warning", {
  expect_warning(result <- assess_pos("not a model", K = 2),
                 "SEMinR models")
  expect_null(result)
})

test_that("K < 2 raises error", {
  expect_error(assess_pos(pls_model, K = 1), "K must be an integer >= 2")
})

test_that("K = 0 raises error", {
  expect_error(assess_pos(pls_model, K = 0), "K must be an integer >= 2")
})

test_that("non-integer K raises error", {
  expect_error(assess_pos(pls_model, K = 2.5), "K must be an integer >= 2")
})

test_that("too many segments for sample size raises error", {
  expect_error(assess_pos(pls_model, K = 50, min_segment_size = 20),
               "too small")
})

# =============================================================================
# pos_segments() TESTS
# =============================================================================

test_that("pos_segments returns list of segment models", {
  segs <- pos_segments(pos_k2)
  expect_true(is.list(segs))
  expect_equal(length(segs), 2)
  for (k in 1:2) {
    expect_true(inherits(segs[[k]], "seminr_model"))
  }
})

test_that("pos_segments rejects non-pos_analysis input", {
  expect_error(pos_segments("not_pos"), "pos_analysis object")
})

# =============================================================================
# S3 METHOD TESTS
# =============================================================================

test_that("print.pos_analysis returns invisibly", {
  expect_invisible(print(pos_k2))
})

test_that("print output contains key information", {
  out <- capture.output(print(pos_k2))
  out_text <- paste(out, collapse = "\n")

  expect_true(grepl("PLS-POS Analysis", out_text))
  expect_true(grepl("Segments: 2", out_text))
  expect_true(grepl("Objective", out_text))
  expect_true(grepl("Segment Sizes", out_text))
  expect_true(grepl("Path Coefficients", out_text))
})

test_that("summary returns summary.pos_analysis object", {
  s <- summary(pos_k2)
  expect_s3_class(s, "summary.pos_analysis")
})

test_that("summary print includes detailed output", {
  out <- capture.output(print(summary(pos_k2)))
  out_text <- paste(out, collapse = "\n")

  expect_true(grepl("Detailed Summary", out_text))
  expect_true(grepl("Global", out_text))
  expect_true(grepl("Path Coefficients per Segment", out_text))
})

test_that("plot.pos_analysis runs without error for all types", {
  expect_no_error(plot(pos_k2, type = "segments"))
  expect_no_error(plot(pos_k2, type = "rsquared"))
  expect_no_error(plot(pos_k2, type = "paths"))
})

# =============================================================================
# COMPARISON TESTS
# =============================================================================

test_that("assess_pos_compare returns pos_comparison object", {
  comp <- assess_pos_compare(pls_model, K_range = 2:3,
                              nstart = 2, max_iter = 10, seed = 123)
  expect_s3_class(comp, "pos_comparison")
})

test_that("comparison contains solutions for each K", {
  comp <- assess_pos_compare(pls_model, K_range = 2:3,
                              nstart = 2, max_iter = 10, seed = 123)
  expect_true("K2" %in% names(comp$solutions))
  expect_true("K3" %in% names(comp$solutions))
})

test_that("comparison fit_table has correct structure", {
  comp <- assess_pos_compare(pls_model, K_range = 2:3,
                              nstart = 2, max_iter = 10, seed = 123)
  expect_true("K" %in% names(comp$fit_table))
  expect_true("Sum_R2" %in% names(comp$fit_table))
  expect_true("Converged" %in% names(comp$fit_table))
  expect_equal(nrow(comp$fit_table), 2)
})

test_that("objective increases with more segments", {
  comp <- assess_pos_compare(pls_model, K_range = 2:3,
                              nstart = 2, max_iter = 10, seed = 123)
  r2_vals <- comp$fit_table$Sum_R2
  expect_true(r2_vals[2] >= r2_vals[1])
})

test_that("comparison print runs without error", {
  comp <- assess_pos_compare(pls_model, K_range = 2:3,
                              nstart = 2, max_iter = 10, seed = 123)
  expect_output(print(comp), "PLS-POS Comparison")
})

test_that("comparison plot runs without error", {
  comp <- assess_pos_compare(pls_model, K_range = 2:3,
                              nstart = 2, max_iter = 10, seed = 123)
  expect_no_error(plot(comp))
})

# =============================================================================
# MODEL TYPE TESTS: MEDIATION
# =============================================================================

test_that("PLS-POS works with mediation model", {
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

  result <- assess_pos(pls_med, K = 2, nstart = 2, max_iter = 10, seed = 42)

  expect_s3_class(result, "pos_analysis")
  expect_equal(sum(result$segment_sizes), nrow(mobi))
  expect_true(result$objective > 0)
  # Mediated endogenous constructs should be in the objective
  expect_true("Satisfaction" %in% result$endogenous)
  expect_true("Loyalty" %in% result$endogenous)
})

# =============================================================================
# MODEL TYPE TESTS: MODERATION
# =============================================================================

test_that("PLS-POS works with moderation model (two-stage)", {
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

  result <- assess_pos(pls_mod, K = 2, nstart = 2, max_iter = 10, seed = 42)

  expect_s3_class(result, "pos_analysis")
  expect_equal(sum(result$segment_sizes), nrow(mobi))

  # Interaction construct should NOT be in endogenous list
  expect_false("Image*Value" %in% result$endogenous)

  # Satisfaction and Loyalty should be endogenous
  expect_true("Satisfaction" %in% result$endogenous)
  expect_true("Loyalty" %in% result$endogenous)

  # Interaction construct should still appear as a predictor in segment paths
  expect_true("Image*Value" %in% rownames(result$segment_paths[[1]]))
})

# =============================================================================
# MODEL TYPE TESTS: HOC (Higher-Order Constructs)
# =============================================================================

test_that("PLS-POS works with HOC two-stage model", {
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

  result <- assess_pos(pls_hoc, K = 2, nstart = 2, max_iter = 10, seed = 42)

  expect_s3_class(result, "pos_analysis")
  expect_equal(sum(result$segment_sizes), nrow(mobi))
  expect_true(result$objective > 0)

  # Quality (HOC) should appear as predictor in segment paths
  expect_true("Quality" %in% rownames(result$segment_paths[[1]]))

  # Satisfaction and Loyalty are endogenous
  expect_true("Satisfaction" %in% result$endogenous)
  expect_true("Loyalty" %in% result$endogenous)
})

# =============================================================================
# INTERNAL HELPER TESTS
# =============================================================================

test_that("pos_endogenous excludes interaction constructs", {
  mm_mod <- constructs(
    composite("Image",        multi_items("IMAG", 1:5)),
    composite("Value",        multi_items("PERV", 1:2)),
    composite("Satisfaction", multi_items("CUSA", 1:3)),
    interaction_term(iv = "Image", moderator = "Value", method = two_stage)
  )
  sm_mod <- relationships(
    paths(from = c("Image", "Value", "Image*Value"), to = "Satisfaction")
  )
  pls_mod <- estimate_pls(data = mobi, measurement_model = mm_mod,
                            structural_model = sm_mod)

  endo <- seminrExtras:::pos_endogenous(pls_mod)
  expect_true("Satisfaction" %in% endo)
  expect_false("Image*Value" %in% endo)
})

test_that("compute_pos_objective returns sum of R-squared", {
  seg_models <- pos_k2$segment_models
  endogenous <- pos_k2$endogenous
  obj <- seminrExtras:::compute_pos_objective(seg_models, endogenous)
  expect_true(is.numeric(obj))
  expect_true(obj > 0)
})

test_that("compute_pos_objective returns -Inf for NULL model", {
  obj <- seminrExtras:::compute_pos_objective(list(NULL, NULL), pos_k2$endogenous)
  expect_equal(obj, -Inf)
})

test_that("random_partition produces valid assignment", {
  asgn <- seminrExtras:::random_partition(100, 3, 10)
  expect_equal(length(asgn), 100)
  expect_true(all(asgn %in% 1:3))
  # Each segment has at least min_size
  expect_true(all(tabulate(asgn, nbins = 3) >= 10))
})

test_that("compute_pos_distances returns N x K matrix", {
  sq_resid <- seminrExtras:::compute_structural_residuals(
    pls_model, pos_k2$segment_models, pos_k2$endogenous)
  dists <- seminrExtras:::compute_pos_distances(sq_resid)
  expect_equal(nrow(dists), nrow(pls_model$construct_scores))
  expect_equal(ncol(dists), 2)
  expect_true(all(dists >= 0))
})
