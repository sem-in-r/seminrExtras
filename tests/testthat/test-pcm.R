# =============================================================================
# Tests for PCM (Predictive Contribution of the Mediator)
# =============================================================================
library(testthat)
library(seminr)
library(seminrExtras)

# =============================================================================
# FIXTURES: Pre-compute models once
# =============================================================================

# Simple mediation model: Image -> Satisfaction -> Loyalty
mobi_mm_simple <- constructs(
  composite("Image",        multi_items("IMAG", 1:5)),
  composite("Satisfaction", multi_items("CUSA", 1:3)),
  composite("Loyalty",      multi_items("CUSL", 1:3))
)
mobi_sm_simple <- relationships(
  paths(from = "Image", to = "Satisfaction"),
  paths(from = "Satisfaction", to = "Loyalty"),
  paths(from = "Image", to = "Loyalty")
)
pls_simple <- estimate_pls(mobi, mobi_mm_simple, mobi_sm_simple)

# Multi-mediator model: ECSI paths to Loyalty
mobi_mm_multi <- constructs(
  composite("Image",        multi_items("IMAG", 1:5)),
  composite("Expectation",  multi_items("CUEX", 1:3)),
  composite("Value",        multi_items("PERV", 1:2)),
  composite("Satisfaction", multi_items("CUSA", 1:3)),
  composite("Loyalty",      multi_items("CUSL", 1:3))
)
mobi_sm_multi <- relationships(
  paths(from = "Image",       to = c("Expectation", "Satisfaction", "Loyalty")),
  paths(from = "Expectation", to = c("Value", "Satisfaction")),
  paths(from = "Value",       to = "Satisfaction"),
  paths(from = "Satisfaction", to = "Loyalty")
)
pls_multi <- estimate_pls(mobi, mobi_mm_multi, mobi_sm_multi)

# =============================================================================
# Run PCM on fixtures (low reps for speed)
# =============================================================================
pcm_simple <- assess_pcm(pls_simple, target = "Loyalty", noFolds = 10, reps = 2)
pcm_multi  <- assess_pcm(pls_multi, target = "Loyalty", noFolds = 10, reps = 2)

# =============================================================================
# TESTS: Internal helpers
# =============================================================================

test_that("detect_final_endogenous identifies correct construct", {
  result <- seminrExtras:::detect_final_endogenous(pls_simple)
  expect_equal(result, "Loyalty")

  result2 <- seminrExtras:::detect_final_endogenous(pls_multi)
  expect_equal(result2, "Loyalty")
})

test_that("find_mediation_paths finds correct paths for simple model", {
  paths <- seminrExtras:::find_mediation_paths(pls_simple, "Loyalty")
  expect_length(paths, 1)
  expect_equal(paths[[1]]$antecedent, "Image")
  expect_equal(paths[[1]]$mediator, "Satisfaction")
  expect_equal(paths[[1]]$target, "Loyalty")
})

test_that("find_mediation_paths finds correct paths for multi-mediator model", {
  paths <- seminrExtras:::find_mediation_paths(pls_multi, "Loyalty")
  # Satisfaction has 3 antecedents: Image, Expectation, Value
  expect_length(paths, 3)

  antecedents <- sapply(paths, function(p) p$antecedent)
  mediators   <- sapply(paths, function(p) p$mediator)
  expect_true(all(mediators == "Satisfaction"))
  expect_true("Image" %in% antecedents)
  expect_true("Expectation" %in% antecedents)
  expect_true("Value" %in% antecedents)
})

test_that("find_mediation_paths returns empty for non-mediated construct", {
  # Satisfaction has no mediation paths to Loyalty in simple model
  # because its antecedent (Image) does not have its own antecedent
  paths <- seminrExtras:::find_mediation_paths(pls_simple, "Satisfaction")
  expect_length(paths, 0)
})

test_that("classify_pcm returns correct classifications", {
  expect_equal(seminrExtras:::classify_pcm(-0.05), "Negative")
  expect_equal(seminrExtras:::classify_pcm(0.00),  "Weak")
  expect_equal(seminrExtras:::classify_pcm(0.03),  "Weak")
  expect_equal(seminrExtras:::classify_pcm(0.049), "Weak")
  expect_equal(seminrExtras:::classify_pcm(0.05),  "Moderate")
  expect_equal(seminrExtras:::classify_pcm(0.07),  "Moderate")
  expect_equal(seminrExtras:::classify_pcm(0.10),  "Strong")
  expect_equal(seminrExtras:::classify_pcm(0.20),  "Strong")
  expect_equal(seminrExtras:::classify_pcm(NA),    "NA")
})

# =============================================================================
# TESTS: build_isolated_sub_model
# =============================================================================

test_that("build_isolated_sub_model creates valid sub-model", {
  sub <- seminrExtras:::build_isolated_sub_model(
    pls_simple, "Image", "Satisfaction", "Loyalty"
  )
  expect_s3_class(sub, "seminr_model")
  expect_equal(sort(colnames(sub$construct_scores)),
               sort(c("Image", "Satisfaction", "Loyalty")))

  # Should have 3 structural paths: Image->Satisfaction, Satisfaction->Loyalty, Image->Loyalty
  sm <- sub$smMatrix
  expect_equal(nrow(sm), 3)
})

test_that("build_isolated_sub_model preserves measurement mode", {
  # Use a model with mode_B construct
  mm_mixed <- constructs(
    composite("Image",        multi_items("IMAG", 1:5), weights = mode_B),
    composite("Satisfaction", multi_items("CUSA", 1:3)),
    composite("Loyalty",      multi_items("CUSL", 1:3))
  )
  sm_mixed <- relationships(
    paths(from = "Image", to = "Satisfaction"),
    paths(from = "Satisfaction", to = "Loyalty"),
    paths(from = "Image", to = "Loyalty")
  )
  model_mixed <- estimate_pls(mobi, mm_mixed, sm_mixed)

  sub <- seminrExtras:::build_isolated_sub_model(
    model_mixed, "Image", "Satisfaction", "Loyalty"
  )

  # Image should be mode B in the sub-model
  mode <- sub$mmMatrix[sub$mmMatrix[, "construct"] == "Image", "type"]
  expect_true(all(mode == "B"))
})

# =============================================================================
# TESTS: assess_pcm output structure
# =============================================================================

test_that("assess_pcm returns correct class", {
  expect_s3_class(pcm_simple, "pcm_analysis")
})

test_that("assess_pcm output has expected components", {
  expect_equal(pcm_simple$target, "Loyalty")
  expect_equal(pcm_simple$noFolds, 10)
  expect_equal(pcm_simple$reps, 2)
  expect_length(pcm_simple$mediation_paths, 1)
  expect_length(pcm_simple$pcm_results, 1)
})

test_that("PCM results have correct structure per path", {
  res <- pcm_simple$pcm_results[[1]]
  expect_equal(res$antecedent, "Image")
  expect_equal(res$mediator, "Satisfaction")
  expect_equal(res$target, "Loyalty")

  # Results matrix should have 3 rows (CUSL1-3) and 6 columns
  expect_equal(nrow(res$results), 3)
  expect_equal(ncol(res$results), 6)
  expect_equal(colnames(res$results),
               c("RMSE_DA", "RMSE_EA", "PCM_RMSE", "MAE_DA", "MAE_EA", "PCM_MAE"))
  expect_equal(sort(rownames(res$results)), sort(c("CUSL1", "CUSL2", "CUSL3")))
})

test_that("PCM values are finite and within expected range", {
  for (res in pcm_simple$pcm_results) {
    expect_true(all(is.finite(res$pcm_rmse)))
    expect_true(all(is.finite(res$pcm_mae)))
    # PCM should typically be between -1 and 1
    expect_true(all(res$pcm_rmse > -1 & res$pcm_rmse < 1))
    expect_true(all(res$pcm_mae  > -1 & res$pcm_mae  < 1))
  }
})

test_that("PCM RMSE and MAE are derived correctly from DA and EA metrics", {
  res <- pcm_simple$pcm_results[[1]]
  mat <- res$results
  # PCM = (EA - DA) / EA
  expected_pcm_rmse <- (mat[, "RMSE_EA"] - mat[, "RMSE_DA"]) / mat[, "RMSE_EA"]
  expected_pcm_mae  <- (mat[, "MAE_EA"]  - mat[, "MAE_DA"])  / mat[, "MAE_EA"]
  expect_equal(as.vector(mat[, "PCM_RMSE"]), as.vector(expected_pcm_rmse))
  expect_equal(as.vector(mat[, "PCM_MAE"]),  as.vector(expected_pcm_mae))
})

test_that("DA RMSE and MAE are strictly positive", {
  for (res in pcm_simple$pcm_results) {
    expect_true(all(res$results[, "RMSE_DA"] > 0))
    expect_true(all(res$results[, "MAE_DA"]  > 0))
    expect_true(all(res$results[, "RMSE_EA"] > 0))
    expect_true(all(res$results[, "MAE_EA"]  > 0))
  }
})

# =============================================================================
# TESTS: Multi-mediator model
# =============================================================================

test_that("assess_pcm handles multiple mediation paths", {
  expect_length(pcm_multi$pcm_results, 3)
  mediators <- sapply(pcm_multi$pcm_results, function(r) r$mediator)
  expect_true(all(mediators == "Satisfaction"))
})

test_that("Multi-mediator PCM has results per indicator per path", {
  for (res in pcm_multi$pcm_results) {
    expect_equal(nrow(res$results), 3)  # CUSL1, CUSL2, CUSL3
    expect_equal(ncol(res$results), 6)
  }
})

# =============================================================================
# TESTS: Auto-detection and validation
# =============================================================================

test_that("assess_pcm auto-detects single final endogenous target", {
  pcm_auto <- assess_pcm(pls_simple, noFolds = 10, reps = 2)
  expect_equal(pcm_auto$target, "Loyalty")
})

test_that("assess_pcm rejects non-seminr model", {
  expect_warning(
    result <- assess_pcm(list(not = "a model")),
    "only works with SEMinR"
  )
  expect_null(result)
})

test_that("assess_pcm rejects invalid target", {
  expect_error(
    assess_pcm(pls_simple, target = "Nonexistent"),
    "not found in model"
  )
})

test_that("assess_pcm rejects target with no mediation", {
  # Satisfaction is endogenous but its predictor (Image) has no antecedents
  expect_error(
    assess_pcm(pls_simple, target = "Satisfaction"),
    "No mediation paths"
  )
})

test_that("assess_pcm rejects invalid noFolds", {
  expect_error(assess_pcm(pls_simple, noFolds = 1), "noFolds must be")
  expect_error(assess_pcm(pls_simple, noFolds = "a"), "noFolds must be")
})

test_that("assess_pcm rejects invalid reps", {
  expect_error(assess_pcm(pls_simple, reps = 0), "reps must be")
})

# =============================================================================
# TESTS: S3 methods
# =============================================================================

test_that("print.pcm_analysis produces output", {
  output <- capture.output(print(pcm_simple))
  expect_true(any(grepl("Predictive Contribution", output)))
  expect_true(any(grepl("Loyalty", output)))
  expect_true(any(grepl("Image -> Satisfaction -> Loyalty", output)))
  expect_true(any(grepl("PCM", output)))
})

test_that("summary.pcm_analysis returns correct class", {
  s <- summary(pcm_simple)
  expect_s3_class(s, "summary.pcm_analysis")
})

test_that("print.summary.pcm_analysis shows detailed results", {
  s <- summary(pcm_simple)
  output <- capture.output(print(s))
  expect_true(any(grepl("RMSE_DA", output)))
  expect_true(any(grepl("RMSE_EA", output)))
  expect_true(any(grepl("PCM_RMSE", output)))
  expect_true(any(grepl("CUSL", output)))
  expect_true(any(grepl("Danks", output)))
})

# =============================================================================
# TESTS: Interaction construct exclusion
# =============================================================================

test_that("PCM excludes interaction constructs from mediation paths", {
  mobi_mm_mod <- constructs(
    composite("Image",        multi_items("IMAG", 1:5)),
    composite("Value",        multi_items("PERV", 1:2)),
    composite("Satisfaction", multi_items("CUSA", 1:3)),
    composite("Loyalty",      multi_items("CUSL", 1:3)),
    interaction_term(iv = "Image", moderator = "Value")
  )
  mobi_sm_mod <- relationships(
    paths(from = "Image", to = "Satisfaction"),
    paths(from = "Value", to = "Satisfaction"),
    paths(from = "Image*Value", to = "Satisfaction"),
    paths(from = "Satisfaction", to = "Loyalty"),
    paths(from = "Image", to = "Loyalty")
  )
  model_mod <- estimate_pls(mobi, mobi_mm_mod, mobi_sm_mod)

  paths <- seminrExtras:::find_mediation_paths(model_mod, "Loyalty")
  mediators <- sapply(paths, function(p) p$mediator)
  antecedents <- sapply(paths, function(p) p$antecedent)

  # No interaction construct should appear as mediator or antecedent

  expect_false(any(grepl("\\*", mediators)))
  expect_false(any(grepl("\\*", antecedents)))
})

test_that("PCM works end-to-end with moderation model", {
  mobi_mm_mod <- constructs(
    composite("Image",        multi_items("IMAG", 1:5)),
    composite("Value",        multi_items("PERV", 1:2)),
    composite("Satisfaction", multi_items("CUSA", 1:3)),
    composite("Loyalty",      multi_items("CUSL", 1:3)),
    interaction_term(iv = "Image", moderator = "Value")
  )
  mobi_sm_mod <- relationships(
    paths(from = "Image", to = "Satisfaction"),
    paths(from = "Value", to = "Satisfaction"),
    paths(from = "Image*Value", to = "Satisfaction"),
    paths(from = "Satisfaction", to = "Loyalty"),
    paths(from = "Image", to = "Loyalty")
  )
  model_mod <- estimate_pls(mobi, mobi_mm_mod, mobi_sm_mod)

  pcm_mod <- assess_pcm(model_mod, target = "Loyalty", noFolds = 10, reps = 2)
  expect_s3_class(pcm_mod, "pcm_analysis")
  # Should find mediation paths through Satisfaction, excluding Image*Value
  expect_true(length(pcm_mod$pcm_results) >= 1)
  for (res in pcm_mod$pcm_results) {
    expect_false(grepl("\\*", res$mediator))
    expect_false(grepl("\\*", res$antecedent))
  }
})

# =============================================================================
# TESTS: HOC guard
# =============================================================================

test_that("PCM skips paths involving HOC constructs with warning", {
  # Two-stage HOC model: Image (HOC) with LOCs
  mobi_mm_hoc <- constructs(
    composite("Tangibles",  multi_items("IMAG", 1:3)),
    composite("Intangibles", multi_items("IMAG", 4:5)),
    higher_composite("Image", c("Tangibles", "Intangibles")),
    composite("Satisfaction", multi_items("CUSA", 1:3)),
    composite("Loyalty",      multi_items("CUSL", 1:3))
  )
  mobi_sm_hoc <- relationships(
    paths(from = "Image", to = "Satisfaction"),
    paths(from = "Satisfaction", to = "Loyalty"),
    paths(from = "Image", to = "Loyalty")
  )
  model_hoc <- estimate_pls(mobi, mobi_mm_hoc, mobi_sm_hoc)

  # Should warn about HOC and return no valid paths
  expect_warning(
    paths <- seminrExtras:::find_mediation_paths(model_hoc, "Loyalty"),
    "higher-order"
  )
  expect_length(paths, 0)
})

# =============================================================================
# TESTS: Plot method
# =============================================================================

test_that("plot.pcm_analysis runs without error", {
  expect_no_error(plot(pcm_simple))
})

test_that("plot.pcm_analysis accepts MAE metric", {
  expect_no_error(plot(pcm_simple, metric = "MAE"))
})

test_that("plot.pcm_analysis handles multiple paths", {
  expect_no_error(plot(pcm_multi))
})
