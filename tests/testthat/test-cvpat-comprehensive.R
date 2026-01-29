library(seminr)

# ============================================================================
# Setup: Create models for testing
# ============================================================================

set.seed(123)

# Basic measurement model
corp_rep_mm <- constructs(
  composite("COMP", multi_items("comp_", 1:3)),
  composite("LIKE", multi_items("like_", 1:3)),
  composite("CUSA", single_item("cusa")),
  composite("CUSL", multi_items("cusl_", 1:3))
)

# Structural model 1
sm_one <- relationships(
  paths(from = c("COMP", "LIKE"), to = c("CUSA")),
  paths(from = c("CUSA"), to = c("CUSL"))
)

# Structural model 2 (alternative paths)
sm_two <- relationships(
  paths(from = c("COMP", "LIKE"), to = c("CUSA", "CUSL")),
  paths(from = c("CUSA"), to = c("CUSL"))
)

model_one <- estimate_pls(
  data = corp_rep_data,
  measurement_model = corp_rep_mm,
  structural_model  = sm_one,
  missing = mean_replacement,
  missing_value = "-99"
)

model_two <- estimate_pls(
  data = corp_rep_data,
  measurement_model = corp_rep_mm,
  structural_model  = sm_two,
  missing = mean_replacement,
  missing_value = "-99"
)

# ============================================================================
# assess_cvpat - Output Structure Tests
# ============================================================================

test_that("assess_cvpat returns list with two elements", {
  result <- assess_cvpat(model_one, nboot = 50, seed = 123, cores = 1)

  expect_type(result, "list")
  expect_named(result, c("CVPAT_compare_LM", "CVPAT_compare_IA"))
})

test_that("assess_cvpat returns matrices with table_output class", {
  result <- assess_cvpat(model_one, nboot = 50, seed = 123, cores = 1)

  expect_true(inherits(result$CVPAT_compare_LM, "matrix"))
  expect_true(inherits(result$CVPAT_compare_IA, "matrix"))
  expect_true("table_output" %in% class(result$CVPAT_compare_LM))
  expect_true("table_output" %in% class(result$CVPAT_compare_IA))
})

test_that("assess_cvpat returns correct column names", {
  result <- assess_cvpat(model_one, nboot = 50, seed = 123, cores = 1)

  lm_cols <- c("PLS Loss", "LM Loss", "Diff", "Boot T value", "Boot P Value")
  ia_cols <- c("PLS Loss", "IA Loss", "Diff", "Boot T value", "Boot P Value")

  expect_equal(colnames(result$CVPAT_compare_LM), lm_cols)
  expect_equal(colnames(result$CVPAT_compare_IA), ia_cols)
})

test_that("assess_cvpat returns correct row names with Overall", {

  result <- assess_cvpat(model_one, nboot = 50, seed = 123, cores = 1)

  # Should have endogenous constructs plus "Overall"
  expect_true("Overall" %in% rownames(result$CVPAT_compare_LM))
  expect_true("Overall" %in% rownames(result$CVPAT_compare_IA))
  expect_true("CUSA" %in% rownames(result$CVPAT_compare_LM))
  expect_true("CUSL" %in% rownames(result$CVPAT_compare_LM))
})

test_that("assess_cvpat has comment attribute", {
  result <- assess_cvpat(model_one, nboot = 50, seed = 123, cores = 1)

  expect_false(is.null(comment(result$CVPAT_compare_LM)))
  expect_false(is.null(comment(result$CVPAT_compare_IA)))
  expect_true(grepl("CVPAT", comment(result$CVPAT_compare_LM)))
})

# ============================================================================
# assess_cvpat - Reproducibility Tests
# ============================================================================

test_that("assess_cvpat is reproducible with same seed", {
  result1 <- assess_cvpat(model_one, nboot = 50, seed = 42, cores = 1)
  result2 <- assess_cvpat(model_one, nboot = 50, seed = 42, cores = 1)

  expect_equal(as.numeric(result1$CVPAT_compare_LM), as.numeric(result2$CVPAT_compare_LM))
  expect_equal(as.numeric(result1$CVPAT_compare_IA), as.numeric(result2$CVPAT_compare_IA))
})

test_that("assess_cvpat differs with different seeds", {
  result1 <- assess_cvpat(model_one, nboot = 50, seed = 42, cores = 1)
  result2 <- assess_cvpat(model_one, nboot = 50, seed = 99, cores = 1)

  # Bootstrap values should differ
  expect_false(identical(
    result1$CVPAT_compare_LM[, "Boot T value"],
    result2$CVPAT_compare_LM[, "Boot T value"]
  ))
})

# ============================================================================
# assess_cvpat - Input Validation Tests
# ============================================================================

test_that("assess_cvpat rejects non-seminr model objects", {
  expect_message(
    result <- assess_cvpat(list(not = "a_model"), nboot = 50, cores = 1),
    "This function only works with SEMinR models"
  )
  expect_null(result)
})

test_that("assess_cvpat rejects data frame input", {
  expect_message(
    result <- assess_cvpat(data.frame(x = 1:10), nboot = 50, cores = 1),
    "This function only works with SEMinR models"
  )
  expect_null(result)
})

# ============================================================================
# assess_cvpat - Output Value Tests
# ============================================================================

test_that("assess_cvpat loss values are positive", {
  result <- assess_cvpat(model_one, nboot = 50, seed = 123, cores = 1)

  expect_true(all(result$CVPAT_compare_LM[, "PLS Loss"] >= 0))
  expect_true(all(result$CVPAT_compare_LM[, "LM Loss"] >= 0))
  expect_true(all(result$CVPAT_compare_IA[, "PLS Loss"] >= 0))
  expect_true(all(result$CVPAT_compare_IA[, "IA Loss"] >= 0))
})

test_that("assess_cvpat p-values are between 0 and 1", {
  result <- assess_cvpat(model_one, nboot = 50, seed = 123, cores = 1)

  expect_true(all(result$CVPAT_compare_LM[, "Boot P Value"] >= 0))
  expect_true(all(result$CVPAT_compare_LM[, "Boot P Value"] <= 1))
  expect_true(all(result$CVPAT_compare_IA[, "Boot P Value"] >= 0))
  expect_true(all(result$CVPAT_compare_IA[, "Boot P Value"] <= 1))
})

test_that("assess_cvpat Diff has correct sign relationship to losses", {
  result <- assess_cvpat(model_one, nboot = 50, seed = 123, cores = 1)

  # Diff should have same sign as (PLS Loss - LM Loss)
  lm_sign <- sign(result$CVPAT_compare_LM[, "PLS Loss"] - result$CVPAT_compare_LM[, "LM Loss"])
  ia_sign <- sign(result$CVPAT_compare_IA[, "PLS Loss"] - result$CVPAT_compare_IA[, "IA Loss"])

  expect_equal(sign(result$CVPAT_compare_LM[, "Diff"]), lm_sign)
  expect_equal(sign(result$CVPAT_compare_IA[, "Diff"]), ia_sign)

  # Diff magnitude should be close to difference of losses
  lm_diff <- result$CVPAT_compare_LM[, "PLS Loss"] - result$CVPAT_compare_LM[, "LM Loss"]
  expect_equal(as.numeric(result$CVPAT_compare_LM[, "Diff"]), as.numeric(lm_diff), tolerance = 0.01)
})

# ============================================================================
# assess_cvpat - Prediction Technique Tests
# ============================================================================

test_that("assess_cvpat works with predict_EA technique", {
  result <- assess_cvpat(model_one, technique = predict_EA, nboot = 50, seed = 123, cores = 1)

  expect_type(result, "list")
  expect_true(all(is.finite(as.numeric(result$CVPAT_compare_LM))))
})

test_that("assess_cvpat works with predict_DA technique", {
  result <- assess_cvpat(model_one, technique = predict_DA, nboot = 50, seed = 123, cores = 1)

  expect_type(result, "list")
  expect_true(all(is.finite(as.numeric(result$CVPAT_compare_LM))))
})

test_that("assess_cvpat gives different results for EA vs DA techniques", {
  result_ea <- assess_cvpat(model_one, technique = predict_EA, nboot = 50, seed = 123, cores = 1)
  result_da <- assess_cvpat(model_one, technique = predict_DA, nboot = 50, seed = 123, cores = 1)

  # Results should differ since techniques compute predictions differently
  expect_false(identical(
    as.numeric(result_ea$CVPAT_compare_LM),
    as.numeric(result_da$CVPAT_compare_LM)
  ))
})

# ============================================================================
# assess_cvpat - Cross-validation Parameter Tests
# ============================================================================

test_that("assess_cvpat works with custom noFolds", {
  result <- assess_cvpat(model_one, nboot = 50, seed = 123, noFolds = 5, reps = 1, cores = 1)

  expect_type(result, "list")
  expect_true(all(is.finite(as.numeric(result$CVPAT_compare_LM))))
})

# ============================================================================
# assess_cvpat_compare - Output Structure Tests
# ============================================================================

test_that("assess_cvpat_compare returns matrix with table_output class", {
  result <- assess_cvpat_compare(
    established_model = model_one,
    alternative_model = model_two,
    nboot = 50, seed = 123, cores = 1
  )

  expect_true(inherits(result, "matrix"))
  expect_true("table_output" %in% class(result))
})

test_that("assess_cvpat_compare returns correct column names", {
  result <- assess_cvpat_compare(
    established_model = model_one,
    alternative_model = model_two,
    nboot = 50, seed = 123, cores = 1
  )

  expected_cols <- c("Base Model Loss", "Alt Model Loss", "Diff", "Boot T value", "Boot P Value")
  expect_equal(colnames(result), expected_cols)
})

test_that("assess_cvpat_compare returns correct row names with Overall", {
  result <- assess_cvpat_compare(
    established_model = model_one,
    alternative_model = model_two,
    nboot = 50, seed = 123, cores = 1
  )

  expect_true("Overall" %in% rownames(result))
})

test_that("assess_cvpat_compare has comment attribute", {
  result <- assess_cvpat_compare(
    established_model = model_one,
    alternative_model = model_two,
    nboot = 50, seed = 123, cores = 1
  )

  expect_false(is.null(comment(result)))
  expect_true(grepl("CVPAT", comment(result)))
})

# ============================================================================
# assess_cvpat_compare - Reproducibility Tests
# ============================================================================

test_that("assess_cvpat_compare is reproducible with same seed", {
  result1 <- assess_cvpat_compare(model_one, model_two, nboot = 50, seed = 42, cores = 1)
  result2 <- assess_cvpat_compare(model_one, model_two, nboot = 50, seed = 42, cores = 1)

  expect_equal(as.numeric(result1), as.numeric(result2))
})

test_that("assess_cvpat_compare differs with different seeds", {
  result1 <- assess_cvpat_compare(model_one, model_two, nboot = 50, seed = 42, cores = 1)
  result2 <- assess_cvpat_compare(model_one, model_two, nboot = 50, seed = 99, cores = 1)

  expect_false(identical(result1[, "Boot T value"], result2[, "Boot T value"]))
})

# ============================================================================
# assess_cvpat_compare - Output Value Tests
# ============================================================================

test_that("assess_cvpat_compare loss values are positive", {
  result <- assess_cvpat_compare(model_one, model_two, nboot = 50, seed = 123, cores = 1)

  expect_true(all(result[, "Base Model Loss"] >= 0))
  expect_true(all(result[, "Alt Model Loss"] >= 0))
})

test_that("assess_cvpat_compare p-values are between 0 and 1", {
  result <- assess_cvpat_compare(model_one, model_two, nboot = 50, seed = 123, cores = 1)

  expect_true(all(result[, "Boot P Value"] >= 0))
  expect_true(all(result[, "Boot P Value"] <= 1))
})

test_that("assess_cvpat_compare Diff equals Base minus Alt Loss", {
  result <- assess_cvpat_compare(model_one, model_two, nboot = 50, seed = 123, cores = 1)

  expected_diff <- result[, "Base Model Loss"] - result[, "Alt Model Loss"]
  expect_equal(result[, "Diff"], expected_diff, tolerance = 1e-10)
})

# ============================================================================
# assess_cvpat_compare - Model Validation Tests
# ============================================================================

test_that("assess_cvpat_compare stops with mismatched endogenous constructs", {
  # Create model with different endogenous construct
  mm_diff <- constructs(
    composite("COMP", multi_items("comp_", 1:3)),
    composite("LIKE", multi_items("like_", 1:3)),
    composite("CUSA", single_item("cusa"))
    # Missing CUSL
  )
  sm_diff <- relationships(
    paths(from = c("COMP", "LIKE"), to = c("CUSA"))
  )
  model_diff <- estimate_pls(
    data = corp_rep_data,
    measurement_model = mm_diff,
    structural_model  = sm_diff,
    missing = mean_replacement,
    missing_value = "-99"
  )

  expect_error(
    assess_cvpat_compare(model_one, model_diff, nboot = 50, cores = 1),
    "identical endogenous"
  )
})

# ============================================================================
# assess_cvpat_compare - Prediction Technique Tests
# ============================================================================

test_that("assess_cvpat_compare works with predict_EA technique", {
  result <- assess_cvpat_compare(
    model_one, model_two,
    technique = predict_EA,
    nboot = 50, seed = 123, cores = 1
  )

  expect_true(inherits(result, "matrix"))
  expect_true(all(is.finite(as.numeric(result))))
})

test_that("assess_cvpat_compare works with predict_DA technique", {
  result <- assess_cvpat_compare(
    model_one, model_two,
    technique = predict_DA,
    nboot = 50, seed = 123, cores = 1
  )

  expect_true(inherits(result, "matrix"))
  expect_true(all(is.finite(as.numeric(result))))
})

# ============================================================================
# assess_cvpat_compare - Testtype Tests
# ============================================================================

test_that("assess_cvpat_compare works with testtype two.sided", {
  result <- assess_cvpat_compare(
    model_one, model_two,
    testtype = "two.sided",
    nboot = 50, seed = 123, cores = 1
  )

  expect_true(inherits(result, "matrix"))
  expect_true(all(is.finite(as.numeric(result))))
})

test_that("assess_cvpat_compare works with testtype greater", {
  result <- assess_cvpat_compare(
    model_one, model_two,
    testtype = "greater",
    nboot = 50, seed = 123, cores = 1
  )

  expect_true(inherits(result, "matrix"))
  expect_true(all(is.finite(as.numeric(result))))
})

test_that("testtype affects p-values appropriately", {
  result_two <- assess_cvpat_compare(
    model_one, model_two,
    testtype = "two.sided",
    nboot = 100, seed = 123, cores = 1
  )
  result_greater <- assess_cvpat_compare(
    model_one, model_two,
    testtype = "greater",
    nboot = 100, seed = 123, cores = 1
  )

  # P-values should generally differ between test types
  # (one-sided vs two-sided tests)
  expect_false(identical(
    result_two[, "Boot P Value"],
    result_greater[, "Boot P Value"]
  ))
})

# ============================================================================
# assess_cvpat - Testtype Tests
# ============================================================================

test_that("assess_cvpat works with testtype greater", {
  result <- assess_cvpat(
    model_one,
    testtype = "greater",
    nboot = 50, seed = 123, cores = 1
  )

  expect_type(result, "list")
  expect_true(all(is.finite(as.numeric(result$CVPAT_compare_LM))))
})

test_that("assess_cvpat testtype affects p-values", {
  result_two <- assess_cvpat(model_one, testtype = "two.sided", nboot = 100, seed = 123, cores = 1)
  result_greater <- assess_cvpat(model_one, testtype = "greater", nboot = 100, seed = 123, cores = 1)

  expect_false(identical(
    result_two$CVPAT_compare_LM[, "Boot P Value"],
    result_greater$CVPAT_compare_LM[, "Boot P Value"]
  ))
})
