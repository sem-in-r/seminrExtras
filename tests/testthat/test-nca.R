library(seminr)

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
                      "bottleneck", "necessary_predictors", "pls_model",
                      "target", "predictors", "ceilings")
  expect_true(all(expected_names %in% names(nca_result)))
})

test_that("nca_raw is NULL when using internal ceilings", {
  expect_null(nca_result$nca_raw)
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

test_that("CE-FDH effect size >= CR-FDH for same predictor", {
  # CE-FDH (step function) always has ceiling zone >= CR-FDH (regression line)
  for (pred in nca_result$predictors) {
    expect_true(nca_result$effect_sizes[pred, "ce_fdh"] >=
                nca_result$effect_sizes[pred, "cr_fdh"] - 1e-10)
  }
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

test_that("bottleneck tables have correct structure", {
  for (ceil in names(nca_result$bottleneck)) {
    bn <- nca_result$bottleneck[[ceil]]
    expect_s3_class(bn, "data.frame")
    expect_true("Satisfaction" %in% colnames(bn))
    expect_true(all(nca_result$predictors %in% colnames(bn)))
    expect_equal(nrow(bn), 11)  # 0, 10, 20, ..., 100
  }
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

test_that("necessary_predictors is precomputed in result", {
  expect_type(nca_result$necessary_predictors, "character")
})

test_that("plot.nca_analysis type='effects' runs without error", {
  expect_no_error(plot(nca_result, type = "effects"))
})

test_that("plot.nca_analysis type='scatter' runs without error (internal)", {
  expect_no_error(plot(nca_result, type = "scatter"))
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
# Edge case: test.rep validation
# ============================================================================

test_that("assess_nca errors on negative test.rep", {
  expect_error(
    assess_nca(pls_model, target = "Satisfaction", test.rep = -1),
    "non-negative integer"
  )
})

test_that("assess_nca errors on non-integer test.rep", {
  expect_error(
    assess_nca(pls_model, target = "Satisfaction", test.rep = 1.5),
    "non-negative integer"
  )
})

# ============================================================================
# Edge case: HOC model (should work -- NCA uses construct scores directly)
# ============================================================================

test_that("assess_nca works with higher-order construct model", {
  hoc_mm <- constructs(
    composite("Image", multi_items("IMAG", 1:5)),
    composite("Value", multi_items("PERV", 1:2)),
    higher_composite("Rep", c("Image", "Value"), method = two_stage),
    composite("Satisfaction", multi_items("CUSA", 1:3)),
    composite("Loyalty", multi_items("CUSL", 1:3))
  )
  hoc_sm <- relationships(
    paths(from = "Rep", to = "Satisfaction"),
    paths(from = "Satisfaction", to = "Loyalty")
  )
  hoc_model <- estimate_pls(data = mobi, measurement_model = hoc_mm,
                             structural_model = hoc_sm)

  result <- assess_nca(hoc_model, target = "Satisfaction", test.rep = 0, seed = 123)

  expect_s3_class(result, "nca_analysis")
  expect_equal(result$predictors, "Rep")
})

# ============================================================================
# Internal NCA algorithm tests
# ============================================================================

test_that("CE-FDH effect size is correct for known case", {
  # Perfect diagonal: (0,0), (1,1) -- full ceiling zone
  # Scope = 1, ceiling zone = (1-0)*(1-0) = 1 (step at x=0 has ceiling_y=0)
  x <- c(0, 1)
  y <- c(0, 1)
  d <- seminrExtras:::ce_fdh_effect_size(x, y)
  expect_equal(d, 1.0)
})

test_that("CE-FDH effect size is zero when no empty space exists", {
  # Upper-left point present: (0, 1), (1, 0)
  # ceiling at x=0 is 1, ceiling at x=1 is 1 (cummax) -> no empty space
  x <- c(0, 1)
  y <- c(1, 0)
  d <- seminrExtras:::ce_fdh_effect_size(x, y)
  expect_equal(d, 0.0)
})

test_that("CE-FDH effect size is bounded [0, 1]", {
  set.seed(42)
  x <- rnorm(100)
  y <- rnorm(100)
  d <- seminrExtras:::ce_fdh_effect_size(x, y)
  expect_true(d >= 0 && d <= 1)
})

test_that("CR-FDH effect size is bounded [0, 1]", {
  set.seed(42)
  x <- rnorm(100)
  y <- rnorm(100)
  d <- seminrExtras:::cr_fdh_effect_size(x, y)
  expect_true(d >= 0 && d <= 1)
})

test_that("CE-FDH handles constant x (zero scope) gracefully", {
  d <- seminrExtras:::ce_fdh_effect_size(rep(1, 10), 1:10)
  expect_equal(d, 0)
})

test_that("CE-FDH handles constant y (zero scope) gracefully", {
  d <- seminrExtras:::ce_fdh_effect_size(1:10, rep(5, 10))
  expect_equal(d, 0)
})

test_that("bottleneck has correct length and bounded values", {
  x <- pls_model$construct_scores[, "Image"]
  y <- pls_model$construct_scores[, "Satisfaction"]
  bn <- seminrExtras:::compute_bottleneck_column(x, y, "ce_fdh", 10)
  expect_equal(length(bn), 11)
  # First value (0% of Y) should be 0 or NA
  expect_true(is.na(bn[1]) || bn[1] == 0)
  # Non-NA values should be in [0, 100]
  non_na <- bn[!is.na(bn)]
  expect_true(all(non_na >= 0 & non_na <= 100))
})

# ============================================================================
# NCA-ESSE tests
# ============================================================================

# Pre-compute ESSE result once (no permutation tests for speed)
esse_result <- assess_nca_esse(pls_model,
                                target = "Satisfaction",
                                thresholds = seq(0, 0.05, by = 0.01),
                                seed = 123)

test_that("assess_nca_esse returns correct S3 class and elements", {
  expect_s3_class(esse_result, "nca_esse")
  expected_names <- c("effect_sizes", "benchmark", "delta", "significance",
                      "pls_model", "target", "predictors", "thresholds",
                      "ceiling", "n_obs")
  expect_true(all(expected_names %in% names(esse_result)))
})

test_that("ESSE effect_sizes matrix has correct dimensions", {
  expect_equal(nrow(esse_result$effect_sizes), 6)  # 0, 1, 2, 3, 4, 5%
  expect_equal(ncol(esse_result$effect_sizes), 2)  # Image, Value
  expect_equal(colnames(esse_result$effect_sizes), c("Image", "Value"))
})

test_that("ESSE benchmark matches analytical formula d = t(1 - ln(t))", {
  thresholds <- esse_result$thresholds
  expected <- ifelse(thresholds == 0, 0, thresholds * (1 - log(thresholds)))
  expect_equal(unname(esse_result$benchmark[, 1]), expected)
})

test_that("ESSE delta equals empirical minus benchmark", {
  expect_equal(esse_result$delta,
               esse_result$effect_sizes - esse_result$benchmark)
})

test_that("ESSE effect sizes are non-negative", {
  expect_true(all(esse_result$effect_sizes >= 0, na.rm = TRUE))
})

test_that("ESSE standard NCA (threshold 0%) matches assess_nca()", {
  # Effect sizes at threshold 0% should match standard NCA
  standard_d <- nca_result$effect_sizes[, "ce_fdh"]
  esse_d <- esse_result$effect_sizes["0%", ]
  expect_equal(esse_d, standard_d, tolerance = 1e-10)
})

test_that("ESSE significance is NULL when test.rep = 0", {
  expect_null(esse_result$significance)
})

test_that("ESSE effect sizes increase with higher thresholds", {
  for (pred in esse_result$predictors) {
    d <- esse_result$effect_sizes[, pred]
    expect_true(d[length(d)] >= d[1])
  }
})

test_that("ESSE errors on invalid thresholds", {
  expect_error(
    assess_nca_esse(pls_model, target = "Satisfaction", thresholds = c(-0.1, 0.5)),
    "thresholds"
  )
})

test_that("ESSE print method runs without error", {
  expect_output(print(esse_result), "NCA-ESSE")
})

test_that("ESSE summary returns correct class", {
  s <- summary(esse_result)
  expect_s3_class(s, "summary.nca_esse")
  expect_true(all(c("tables", "target", "predictors") %in% names(s)))
})

test_that("ESSE plot sensitivity runs without error", {
  expect_no_error(plot(esse_result, type = "sensitivity"))
})

test_that("ESSE plot difference runs without error", {
  expect_no_error(plot(esse_result, type = "difference"))
})

test_that("ESSE warns for non-CE-FDH ceiling", {
  expect_warning(
    assess_nca_esse(pls_model, target = "Satisfaction",
                     ceiling = "cr_fdh",
                     thresholds = c(0, 0.01), seed = 123),
    "CE-FDH"
  )
})

# ============================================================================
# Different dataset
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

