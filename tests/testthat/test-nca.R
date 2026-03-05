library(seminr)

skip_if_not_installed("NCA")

# ============================================================================
# Setup: Estimate model once, reuse across all tests
# ============================================================================

set.seed(123)

mobi_mm <- constructs(
  composite("Image",        multi_items("IMAG", 1:5)),
  composite("Value",        multi_items("PERV", 1:2)),
  composite("Satisfaction", multi_items("CUSA", 1:3)),
  composite("Loyalty",      multi_items("CUSL", 1:3))
)

mobi_sm <- relationships(
  paths(from = c("Image", "Value"), to = "Satisfaction"),
  paths(from = "Satisfaction", to = "Loyalty")
)

pls_model <- estimate_pls(
  data = mobi,
  measurement_model = mobi_mm,
  structural_model  = mobi_sm
)

# Pre-compute NCA result once (minimal reps for speed)
nca_result <- assess_nca(pls_model,
                          target = "Satisfaction",
                          ceilings = c("ce_fdh", "cr_fdh"),
                          test.rep = 0,
                          seed = 123)

# ============================================================================
# Structure tests
# ============================================================================

test_that("assess_nca returns correct S3 class and elements", {
  expect_s3_class(nca_result, "nca_analysis")
  expected_names <- c("nca_raw", "effect_sizes", "significance",
                      "bottleneck", "pls_model", "target",
                      "predictors", "ceilings")
  expect_true(all(expected_names %in% names(nca_result)))
})

test_that("effect_sizes matrix has correct dimensions and names", {
  expect_equal(dim(nca_result$effect_sizes), c(2, 2))
  expect_equal(rownames(nca_result$effect_sizes), c("Image", "Value"))
  expect_equal(colnames(nca_result$effect_sizes), c("ce_fdh", "cr_fdh"))
  expect_true("table_output" %in% class(nca_result$effect_sizes))
})

# ============================================================================
# Auto-detection tests
# ============================================================================

test_that("assess_nca auto-detects direct predictors from smMatrix", {
  expect_equal(nca_result$predictors, c("Image", "Value"))
})

test_that("assess_nca auto-detects single predictor correctly", {
  result <- assess_nca(pls_model, target = "Loyalty", test.rep = 0, seed = 123)
  expect_equal(result$predictors, "Satisfaction")
})

test_that("assess_nca accepts explicit predictors", {
  result <- assess_nca(pls_model, target = "Satisfaction",
                        predictors = "Image", test.rep = 0, seed = 123)
  expect_equal(result$predictors, "Image")
  expect_equal(nrow(result$effect_sizes), 1)
})

# ============================================================================
# Effect sizes tests
# ============================================================================

test_that("effect sizes are non-negative and bounded by 1", {
  expect_true(all(nca_result$effect_sizes >= 0, na.rm = TRUE))
  expect_true(all(nca_result$effect_sizes <= 1, na.rm = TRUE))
})

# ============================================================================
# Significance and bottleneck tests
# ============================================================================

test_that("significance is NULL when test.rep = 0", {
  expect_null(nca_result$significance)
})

test_that("bottleneck is a list keyed by ceiling technique", {
  expect_type(nca_result$bottleneck, "list")
  expect_equal(sort(names(nca_result$bottleneck)),
               sort(c("ce_fdh", "cr_fdh")))
})

# ============================================================================
# Ceiling options tests
# ============================================================================

test_that("assess_nca works with single ceiling technique", {
  result <- assess_nca(pls_model, target = "Satisfaction",
                        ceilings = "ce_fdh", test.rep = 0, seed = 123)
  expect_equal(ncol(result$effect_sizes), 1)
  expect_equal(colnames(result$effect_sizes), "ce_fdh")
})

# ============================================================================
# Validation tests
# ============================================================================

test_that("assess_nca warns on non-seminr model", {
  expect_warning(
    result <- assess_nca(list(not = "a model"), target = "X"),
    "only works with SEMinR models"
  )
  expect_null(result)
})

test_that("assess_nca errors on invalid target", {
  expect_error(assess_nca(pls_model, target = "NONEXISTENT"), "target")
})

test_that("assess_nca errors on invalid predictor", {
  expect_error(
    assess_nca(pls_model, target = "Satisfaction",
                predictors = c("Image", "FAKE")),
    "not found in model constructs"
  )
})

# ============================================================================
# S3 methods tests
# ============================================================================

test_that("print.nca_analysis runs without error", {
  expect_output(print(nca_result), "Necessary Condition Analysis")
})

test_that("summary.nca_analysis returns correct class with expected elements", {
  s <- summary(nca_result)
  expect_s3_class(s, "summary.nca_analysis")
  expect_true(all(c("bottleneck", "necessary_predictors") %in% names(s)))
})

test_that("plot.nca_analysis type='effects' runs without error", {
  expect_no_error(plot(nca_result, type = "effects"))
})

# ============================================================================
# Significance testing (only when test.rep > 0)
# ============================================================================

test_that("significance matrix is correct when test.rep > 0", {
  result_sig <- assess_nca(pls_model, target = "Satisfaction",
                            test.rep = 50, seed = 123)
  expect_false(is.null(result_sig$significance))
  expect_equal(dim(result_sig$significance), dim(result_sig$effect_sizes))
  expect_true(all(result_sig$significance >= 0 &
                    result_sig$significance <= 1, na.rm = TRUE))
})

# ============================================================================
# Edge case: different dataset
# ============================================================================

test_that("assess_nca works with corp_rep_data model", {
  corp_mm <- constructs(
    composite("COMP", multi_items("comp_", 1:3)),
    composite("LIKE", multi_items("like_", 1:3)),
    composite("CUSA", single_item("cusa"))
  )
  corp_sm <- relationships(
    paths(from = c("COMP", "LIKE"), to = "CUSA")
  )
  corp_model <- estimate_pls(
    data = corp_rep_data,
    measurement_model = corp_mm,
    structural_model  = corp_sm,
    missing = mean_replacement,
    missing_value = "-99"
  )

  result <- assess_nca(corp_model, target = "CUSA", test.rep = 0, seed = 123)

  expect_s3_class(result, "nca_analysis")
  expect_equal(result$target, "CUSA")
  expect_equal(sort(result$predictors), sort(c("COMP", "LIKE")))
})
