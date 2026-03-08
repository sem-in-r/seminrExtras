# =============================================================================
# Tests for assess_cipma() — Combined Importance-Performance Map Analysis
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

# Pre-compute results for reuse across tests
cipma_result <- assess_cipma(pls_model, target = "Loyalty",
                              scale_min = 1, scale_max = 10,
                              nca_test.rep = 0, seed = 123)

ipma_result <- assess_cipma(pls_model, target = "Loyalty",
                              scale_min = 1, scale_max = 10,
                              nca = FALSE, seed = 123)

# =============================================================================
# STRUCTURE TESTS
# =============================================================================

test_that("assess_cipma returns cipma_analysis object", {
  expect_s3_class(cipma_result, "cipma_analysis")
})

test_that("cipma result contains all expected components", {
  expected_names <- c("importance_unstd", "importance_std", "performance",
                      "nca", "classification", "target", "constructs",
                      "scale_range", "negative_weight_constructs",
                      "excluded_interactions", "pls_model")
  expect_true(all(expected_names %in% names(cipma_result)))
})

test_that("target is stored correctly", {
  expect_equal(cipma_result$target, "Loyalty")
})

test_that("constructs exclude the target", {
  expect_false("Loyalty" %in% cipma_result$constructs)
})

test_that("all included constructs are in the model", {
  model_constructs <- colnames(pls_model$construct_scores)
  expect_true(all(cipma_result$constructs %in% model_constructs))
})

test_that("scale_range is stored correctly", {
  expect_equal(cipma_result$scale_range, c(1, 10))
})

# =============================================================================
# IMPORTANCE TESTS
# =============================================================================

test_that("standardized total effects match seminr", {
  seminr_te <- summary(pls_model)$total_effects
  for (construct in cipma_result$constructs) {
    expect_equal(cipma_result$importance_std[[construct]],
                 seminr_te[construct, "Loyalty"],
                 tolerance = 1e-10)
  }
})

test_that("importance vectors have correct names", {
  expect_equal(names(cipma_result$importance_std), cipma_result$constructs)
  expect_equal(names(cipma_result$importance_unstd), cipma_result$constructs)
})

test_that("unstandardized total effects are finite and non-zero", {
  expect_true(all(is.finite(cipma_result$importance_unstd)))
  expect_true(all(abs(cipma_result$importance_unstd) > .Machine$double.eps))
})

test_that("unstandardized total effects differ from standardized", {
  # They should differ because rescaling changes SDs
  expect_false(all(abs(cipma_result$importance_unstd - cipma_result$importance_std) <
                     .Machine$double.eps))
})

# =============================================================================
# PERFORMANCE TESTS
# =============================================================================

test_that("performance values are within 0-100 range", {
  expect_true(all(cipma_result$performance >= 0))
  expect_true(all(cipma_result$performance <= 100))
})

test_that("performance has correct names", {
  expect_equal(names(cipma_result$performance), cipma_result$constructs)
})

test_that("performance changes with scale range", {
  result_5pt <- assess_cipma(pls_model, target = "Loyalty",
                               scale_min = 1, scale_max = 5,
                               nca = FALSE, seed = 123)
  result_10pt <- assess_cipma(pls_model, target = "Loyalty",
                                scale_min = 1, scale_max = 10,
                                nca = FALSE, seed = 123)

  # Performance differs with different scales
  expect_false(all(abs(result_5pt$performance - result_10pt$performance) <
                     .Machine$double.eps))
})

test_that("performance is correct for known calculation", {
  # Manual calculation for one construct
  items <- seminr:::items_of_construct("Image", pls_model)
  weights <- pls_model$outer_weights[items, "Image"]
  indicator_means <- colMeans(pls_model$data[, items])
  indicator_perf <- (indicator_means - 1) / (10 - 1) * 100
  expected_perf <- sum(weights * indicator_perf) / sum(weights)

  expect_equal(cipma_result$performance[["Image"]], expected_perf, tolerance = 1e-10)
})

# =============================================================================
# NCA INTEGRATION TESTS
# =============================================================================

test_that("cIPMA includes NCA results", {
  expect_s3_class(cipma_result$nca, "nca_analysis")
})

test_that("IPMA-only mode has NULL NCA", {
  expect_null(ipma_result$nca)
})

test_that("NCA target matches cIPMA target", {
  expect_equal(cipma_result$nca$target, cipma_result$target)
})

test_that("NCA predictors match cIPMA constructs", {
  expect_equal(cipma_result$nca$predictors, cipma_result$constructs)
})

# =============================================================================
# CLASSIFICATION TESTS
# =============================================================================

test_that("classification is a data frame with expected columns", {
  expect_s3_class(cipma_result$classification, "data.frame")
  expected_cols <- c("Construct", "Importance", "Performance",
                     "High_Importance", "Necessary", "Priority")
  expect_true(all(expected_cols %in% colnames(cipma_result$classification)))
})

test_that("classification has one row per construct", {
  expect_equal(nrow(cipma_result$classification), length(cipma_result$constructs))
})

test_that("classification priority is a valid category", {
  valid_priorities <- c("Top priority", "Important driver",
                        "Bottleneck risk", "Low priority")
  expect_true(all(cipma_result$classification$Priority %in% valid_priorities))
})

test_that("High_Importance splits at median", {
  imp <- cipma_result$importance_unstd
  high <- cipma_result$classification$High_Importance
  expect_true(all(imp[high] > median(imp)))
  expect_true(all(imp[!high] <= median(imp)))
})

test_that("classification without NCA has all Necessary = FALSE", {
  expect_true(all(!ipma_result$classification$Necessary))
})

test_that("classification without NCA uses only driver/low priority", {
  valid <- c("Important driver", "Low priority")
  expect_true(all(ipma_result$classification$Priority %in% valid))
})

# =============================================================================
# TOTAL EFFECTS COMPUTATION TESTS
# =============================================================================

test_that("total effects include indirect paths", {
  # Image has both direct and indirect effects on Loyalty
  # Direct: Image -> Loyalty

  # Indirect: Image -> Satisfaction -> Loyalty, Image -> Expectation -> ... -> Loyalty
  direct_coef <- pls_model$path_coef["Image", "Loyalty"]
  total_std <- cipma_result$importance_std[["Image"]]

  # Total effect should exceed direct effect when indirect paths exist
  expect_true(total_std > direct_coef)
})

test_that("only constructs with non-zero total effect on target are included", {
  # All included constructs must have a non-zero total effect on the target
  seminr_te <- summary(pls_model)$total_effects
  for (construct in cipma_result$constructs) {
    expect_true(abs(seminr_te[construct, "Loyalty"]) > .Machine$double.eps)
  }
})

# =============================================================================
# INTERACTION TERM EXCLUSION TESTS
# =============================================================================

test_that("interaction constructs are excluded from IPMA", {
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

  expect_message(
    result_mod <- assess_cipma(pls_mod, target = "Satisfaction",
                                 nca = FALSE, scale_min = 1, scale_max = 10),
    "Excluding interaction"
  )

  expect_false("Image*Value" %in% result_mod$constructs)
  expect_true("Image" %in% result_mod$constructs)
  expect_true("Value" %in% result_mod$constructs)
  expect_true("Image*Value" %in% result_mod$excluded_interactions)
})

# =============================================================================
# INPUT VALIDATION TESTS
# =============================================================================

test_that("non-seminr model input returns NULL with warning", {
  expect_warning(result <- assess_cipma("not a model", target = "X"),
                 "SEMinR models")
  expect_null(result)
})

test_that("invalid target construct raises error", {
  expect_error(assess_cipma(pls_model, target = "NonExistent"),
               "not found in model constructs")
})

test_that("invalid scale_min >= scale_max raises error", {
  expect_error(assess_cipma(pls_model, target = "Loyalty",
                              scale_min = 10, scale_max = 1),
               "scale_min must be less than scale_max")
})

test_that("equal scale_min and scale_max raises error", {
  expect_error(assess_cipma(pls_model, target = "Loyalty",
                              scale_min = 5, scale_max = 5),
               "scale_min must be less than scale_max")
})

test_that("non-numeric scale values raise error", {
  expect_error(assess_cipma(pls_model, target = "Loyalty",
                              scale_min = "a", scale_max = 10),
               "single numeric values")
})

# =============================================================================
# S3 METHOD TESTS
# =============================================================================

test_that("print.cipma_analysis returns invisibly", {
  expect_invisible(print(cipma_result))
})

test_that("print output contains key information", {
  out <- capture.output(print(cipma_result))
  out_text <- paste(out, collapse = "\n")

  expect_true(grepl("cIPMA", out_text))
  expect_true(grepl("Loyalty", out_text))
  expect_true(grepl("Image", out_text))
  expect_true(grepl("Top priority", out_text))
})

test_that("print without NCA shows IPMA header", {
  out <- capture.output(print(ipma_result))
  out_text <- paste(out, collapse = "\n")
  expect_true(grepl("Importance-Performance Map Analysis \\(IPMA\\)", out_text))
})

test_that("summary returns summary.cipma_analysis object", {
  s <- summary(cipma_result)
  expect_s3_class(s, "summary.cipma_analysis")
})

test_that("summary print includes bottleneck table", {
  out <- capture.output(print(summary(cipma_result)))
  out_text <- paste(out, collapse = "\n")
  expect_true(grepl("Bottleneck", out_text))
})

test_that("plot.cipma_analysis runs without error", {
  expect_no_error(plot(cipma_result, type = "cipma"))
  expect_no_error(plot(cipma_result, type = "ipma"))
  expect_no_error(plot(cipma_result, importance_metric = "standardized"))
})

test_that("plot with nca=FALSE defaults to ipma type", {
  expect_message(plot(ipma_result, type = "cipma"),
                 "No NCA results available")
})

# =============================================================================
# REPRODUCIBILITY TESTS
# =============================================================================

test_that("results are deterministic with same seed", {
  r1 <- assess_cipma(pls_model, target = "Loyalty",
                       scale_min = 1, scale_max = 10,
                       nca_test.rep = 0, seed = 42)
  r2 <- assess_cipma(pls_model, target = "Loyalty",
                       scale_min = 1, scale_max = 10,
                       nca_test.rep = 0, seed = 42)
  expect_equal(r1$importance_unstd, r2$importance_unstd)
  expect_equal(r1$performance, r2$performance)
})

# =============================================================================
# DIFFERENT TARGET CONSTRUCT TESTS
# =============================================================================

test_that("IPMA works with a different target (Satisfaction)", {
  result_sat <- assess_cipma(pls_model, target = "Satisfaction",
                               scale_min = 1, scale_max = 10,
                               nca = FALSE, seed = 123)
  expect_equal(result_sat$target, "Satisfaction")
  expect_false("Satisfaction" %in% result_sat$constructs)

  # Only constructs with non-zero total effect on Satisfaction should be included
  seminr_te <- summary(pls_model)$total_effects
  constructs_with_effect <- names(which(abs(seminr_te[, "Satisfaction"]) > .Machine$double.eps))
  constructs_with_effect <- setdiff(constructs_with_effect, "Satisfaction")
  expect_true(all(result_sat$constructs %in% constructs_with_effect))
})
