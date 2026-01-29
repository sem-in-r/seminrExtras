library(seminr)

# Setup: Create a basic model for testing
set.seed(123)

corp_rep_mm <- constructs(
  composite("COMP", multi_items("comp_", 1:3)),
  composite("LIKE", multi_items("like_", 1:3)),
  composite("CUSA", single_item("cusa")),
  composite("CUSL", multi_items("cusl_", 1:3))
)

corp_rep_sm <- relationships(
  paths(from = c("COMP", "LIKE"), to = c("CUSA")),
  paths(from = c("CUSA"), to = c("CUSL"))
)

test_model <- estimate_pls(
  data = corp_rep_data,
  measurement_model = corp_rep_mm,
  structural_model  = corp_rep_sm,
  missing = mean_replacement,
  missing_value = "-99"
)

# ============================================================================
# Basic Functionality Tests
# ============================================================================

test_that("congruence_test returns proper structure", {
  result <- congruence_test(test_model, nboot = 50, seed = 123)

  # Should return a list

  expect_type(result, "list")


  # Should have 'results' element

  expect_named(result, "results")

  # Results should be a matrix with table_output class
  expect_true(inherits(result$results, "matrix"))
  expect_true("table_output" %in% class(result$results))
})

test_that("congruence_test returns correct dimensions", {
  result <- congruence_test(test_model, nboot = 50, seed = 123)

  # Number of rows should equal number of construct pairs
  n_constructs <- ncol(test_model$construct_scores)
  expected_pairs <- choose(n_constructs, 2)
  expect_equal(nrow(result$results), expected_pairs)

  # Should have 6 columns
  expect_equal(ncol(result$results), 6)
})

test_that("congruence_test returns correct column names", {
  result <- congruence_test(test_model, nboot = 50, seed = 123, alpha = 0.05)

  expected_cols <- c("Original Est.", "Diff", "Bootstrap SD", "T Stat.", "2.5% CI", "97.5% CI")
  expect_equal(colnames(result$results), expected_cols)
})

test_that("congruence_test row names contain construct pairs", {
  result <- congruence_test(test_model, nboot = 50, seed = 123)

  # Row names should contain " -> " pattern
  expect_true(all(grepl(" -> ", rownames(result$results))))
})

# ============================================================================
# Reproducibility Tests
# ============================================================================

test_that("congruence_test is reproducible with same seed", {
  result1 <- congruence_test(test_model, nboot = 50, seed = 42)
  result2 <- congruence_test(test_model, nboot = 50, seed = 42)

  expect_equal(result1$results, result2$results)
})

test_that("congruence_test differs with different seeds", {
  result1 <- congruence_test(test_model, nboot = 50, seed = 42)
  result2 <- congruence_test(test_model, nboot = 50, seed = 99)

  # Bootstrap SD and CI columns should differ
  expect_false(identical(result1$results[, "Bootstrap SD"], result2$results[, "Bootstrap SD"]))
})

# ============================================================================
# Parameter Variation Tests
# ============================================================================

test_that("congruence_test respects alpha parameter for CI columns", {
  result_05 <- congruence_test(test_model, nboot = 50, seed = 123, alpha = 0.05)
  result_10 <- congruence_test(test_model, nboot = 50, seed = 123, alpha = 0.10)

  # Column names should reflect alpha

  expect_true("2.5% CI" %in% colnames(result_05$results))
  expect_true("97.5% CI" %in% colnames(result_05$results))
  expect_true("5% CI" %in% colnames(result_10$results))
  expect_true("95% CI" %in% colnames(result_10$results))
})

test_that("congruence_test works with different nboot values", {
  result_small <- congruence_test(test_model, nboot = 20, seed = 123)
  result_large <- congruence_test(test_model, nboot = 100, seed = 123)

  # Both should return valid results

  expect_true(all(is.finite(result_small$results[, "Original Est."])))
  expect_true(all(is.finite(result_large$results[, "Original Est."])))

  # Original estimates should be identical (same model)
  expect_equal(result_small$results[, "Original Est."], result_large$results[, "Original Est."])
})

test_that("congruence_test works with custom threshold", {
  result <- congruence_test(test_model, nboot = 50, seed = 123, threshold = 0.9)

  # Diff column should reflect threshold - original
  # Diff = threshold - abs(original)
  expect_true(all(is.finite(result$results[, "Diff"])))
})

# ============================================================================
# Input Validation Tests
# ============================================================================

test_that("congruence_test rejects non-seminr model objects", {
  expect_message(
    result <- congruence_test(list(not = "a_model"), nboot = 50),
    "This function only works with SEMinR models"
  )
  expect_null(result)
})

test_that("congruence_test rejects NULL input", {
  expect_message(
    result <- congruence_test(NULL, nboot = 50),
    "This function only works with SEMinR models"
  )
  expect_null(result)
})

test_that("congruence_test rejects data frame input", {
  expect_message(
    result <- congruence_test(data.frame(x = 1:10), nboot = 50),
    "This function only works with SEMinR models"
  )
  expect_null(result)
})

# ============================================================================
# Output Value Tests
# ============================================================================

test_that("congruence_test original estimates are bounded", {
  result <- congruence_test(test_model, nboot = 50, seed = 123)

  # Congruence coefficients should be between -1 and 1
  expect_true(all(result$results[, "Original Est."] >= -1))
  expect_true(all(result$results[, "Original Est."] <= 1))
})

test_that("congruence_test bootstrap SD is positive", {
  result <- congruence_test(test_model, nboot = 50, seed = 123)

  expect_true(all(result$results[, "Bootstrap SD"] >= 0))
})

test_that("congruence_test confidence intervals are ordered correctly", {
  result <- congruence_test(test_model, nboot = 100, seed = 123, alpha = 0.05)

  lower_ci <- result$results[, "2.5% CI"]
  upper_ci <- result$results[, "97.5% CI"]

  # Lower CI should be less than or equal to upper CI
  expect_true(all(lower_ci <= upper_ci))
})
