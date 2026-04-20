library(seminr)

# ============================================================================
# Setup: Create models and pre-compute results for testing
# ============================================================================

set.seed(123)

corp_rep_mm <- constructs(
  composite("COMP", multi_items("comp_", 1:3)),
  composite("LIKE", multi_items("like_", 1:3)),
  composite("CUSA", single_item("cusa")),
  composite("CUSL", multi_items("cusl_", 1:3))
)

sm_one <- relationships(
  paths(from = c("COMP", "LIKE"), to = c("CUSA")),
  paths(from = c("CUSA"), to = c("CUSL"))
)

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

# Pre-compute assess_cvpat once (expensive: runs k-fold CV + bootstrap)
cvpat_result <- assess_cvpat(model_one, nboot = 50, seed = 123, cores = 1)

# Pre-compute assess_cvpat_compare once
compare_result <- assess_cvpat_compare(
  established_model = model_one,
  alternative_model = model_two,
  nboot = 50, seed = 123, cores = 1
)

# ============================================================================
# assess_cvpat - Output Structure Tests
# ============================================================================

test_that("assess_cvpat returns list with two elements", {
  expect_type(cvpat_result, "list")
  expect_named(cvpat_result, c("CVPAT_compare_LM", "CVPAT_compare_IA"))
})

test_that("assess_cvpat returns matrices with table_output class", {
  expect_true(inherits(cvpat_result$CVPAT_compare_LM, "matrix"))
  expect_true(inherits(cvpat_result$CVPAT_compare_IA, "matrix"))
  expect_true("table_output" %in% class(cvpat_result$CVPAT_compare_LM))
  expect_true("table_output" %in% class(cvpat_result$CVPAT_compare_IA))
})

test_that("assess_cvpat returns correct column names", {
  lm_cols <- c("PLS Loss", "LM Loss", "Diff", "Boot T value", "Boot P Value")
  ia_cols <- c("PLS Loss", "IA Loss", "Diff", "Boot T value", "Boot P Value")

  expect_equal(colnames(cvpat_result$CVPAT_compare_LM), lm_cols)
  expect_equal(colnames(cvpat_result$CVPAT_compare_IA), ia_cols)
})

test_that("assess_cvpat returns correct row names with Overall", {
  expect_true("Overall" %in% rownames(cvpat_result$CVPAT_compare_LM))
  expect_true("Overall" %in% rownames(cvpat_result$CVPAT_compare_IA))
  expect_true("CUSA" %in% rownames(cvpat_result$CVPAT_compare_LM))
  expect_true("CUSL" %in% rownames(cvpat_result$CVPAT_compare_LM))
})

test_that("assess_cvpat has comment attribute", {
  expect_false(is.null(comment(cvpat_result$CVPAT_compare_LM)))
  expect_false(is.null(comment(cvpat_result$CVPAT_compare_IA)))
  expect_true(grepl("CVPAT", comment(cvpat_result$CVPAT_compare_LM)))
})

# ============================================================================
# assess_cvpat - Reproducibility Tests
# ============================================================================

test_that("assess_cvpat is reproducible with same seed", {
  result2 <- assess_cvpat(model_one, nboot = 50, seed = 123, cores = 1)

  expect_equal(as.numeric(cvpat_result$CVPAT_compare_LM),
               as.numeric(result2$CVPAT_compare_LM))
  expect_equal(as.numeric(cvpat_result$CVPAT_compare_IA),
               as.numeric(result2$CVPAT_compare_IA))
})

test_that("assess_cvpat differs with different seeds", {
  result2 <- assess_cvpat(model_one, nboot = 50, seed = 99, cores = 1)

  expect_false(identical(
    cvpat_result$CVPAT_compare_LM[, "Boot T value"],
    result2$CVPAT_compare_LM[, "Boot T value"]
  ))
})

# ============================================================================
# assess_cvpat - Input Validation Tests
# ============================================================================

test_that("assess_cvpat rejects non-seminr model objects", {
  expect_warning(
    result <- assess_cvpat(list(not = "a_model"), nboot = 50, cores = 1),
    "only works with SEMinR models"
  )
  expect_null(result)
})

test_that("assess_cvpat rejects data frame input", {
  expect_warning(
    result <- assess_cvpat(data.frame(x = 1:10), nboot = 50, cores = 1),
    "only works with SEMinR models"
  )
  expect_null(result)
})

# ============================================================================
# assess_cvpat - Output Value Tests
# ============================================================================

test_that("assess_cvpat loss values are positive", {
  expect_true(all(cvpat_result$CVPAT_compare_LM[, "PLS Loss"] >= 0))
  expect_true(all(cvpat_result$CVPAT_compare_LM[, "LM Loss"] >= 0))
  expect_true(all(cvpat_result$CVPAT_compare_IA[, "PLS Loss"] >= 0))
  expect_true(all(cvpat_result$CVPAT_compare_IA[, "IA Loss"] >= 0))
})

test_that("assess_cvpat p-values are between 0 and 1", {
  expect_true(all(cvpat_result$CVPAT_compare_LM[, "Boot P Value"] >= 0))
  expect_true(all(cvpat_result$CVPAT_compare_LM[, "Boot P Value"] <= 1))
  expect_true(all(cvpat_result$CVPAT_compare_IA[, "Boot P Value"] >= 0))
  expect_true(all(cvpat_result$CVPAT_compare_IA[, "Boot P Value"] <= 1))
})

test_that("assess_cvpat Diff has correct sign relationship to losses", {
  lm_diff <- cvpat_result$CVPAT_compare_LM[, "PLS Loss"] -
             cvpat_result$CVPAT_compare_LM[, "LM Loss"]
  expect_equal(as.numeric(cvpat_result$CVPAT_compare_LM[, "Diff"]),
               as.numeric(lm_diff), tolerance = 0.01)
})

# ============================================================================
# assess_cvpat_compare - Output Structure Tests
# ============================================================================

test_that("assess_cvpat_compare returns matrix with table_output class", {
  expect_true(inherits(compare_result, "matrix"))
  expect_true("table_output" %in% class(compare_result))
})

test_that("assess_cvpat_compare returns correct column names", {
  expected_cols <- c("Base Model Loss", "Alt Model Loss", "Diff",
                     "Boot T value", "Boot P Value")
  expect_equal(colnames(compare_result), expected_cols)
})

test_that("assess_cvpat_compare returns correct row names with Overall", {
  expect_true("Overall" %in% rownames(compare_result))
})

test_that("assess_cvpat_compare has comment attribute", {
  expect_false(is.null(comment(compare_result)))
  expect_true(grepl("CVPAT", comment(compare_result)))
})

# ============================================================================
# assess_cvpat_compare - Reproducibility Tests
# ============================================================================

test_that("assess_cvpat_compare is reproducible with same seed", {
  result2 <- assess_cvpat_compare(model_one, model_two,
                                   nboot = 50, seed = 123, cores = 1)
  expect_equal(as.numeric(compare_result), as.numeric(result2))
})

test_that("assess_cvpat_compare differs with different seeds", {
  result2 <- assess_cvpat_compare(model_one, model_two,
                                   nboot = 50, seed = 99, cores = 1)
  expect_false(identical(compare_result[, "Boot T value"],
                         result2[, "Boot T value"]))
})

# ============================================================================
# assess_cvpat_compare - Output Value Tests
# ============================================================================

test_that("assess_cvpat_compare loss values are positive", {
  expect_true(all(compare_result[, "Base Model Loss"] >= 0))
  expect_true(all(compare_result[, "Alt Model Loss"] >= 0))
})

test_that("assess_cvpat_compare p-values are between 0 and 1", {
  expect_true(all(compare_result[, "Boot P Value"] >= 0))
  expect_true(all(compare_result[, "Boot P Value"] <= 1))
})

test_that("assess_cvpat_compare Diff equals Base minus Alt Loss", {
  expected_diff <- compare_result[, "Base Model Loss"] -
                   compare_result[, "Alt Model Loss"]
  expect_equal(compare_result[, "Diff"], expected_diff, tolerance = 1e-10)
})

# ============================================================================
# assess_cvpat_compare - Model Validation Tests
# ============================================================================

test_that("assess_cvpat_compare stops with mismatched endogenous constructs", {
  mm_diff <- constructs(
    composite("COMP", multi_items("comp_", 1:3)),
    composite("LIKE", multi_items("like_", 1:3)),
    composite("CUSA", single_item("cusa"))
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
# Prediction Technique and Testtype Tests (compute fresh, minimal reps)
# ============================================================================

test_that("assess_cvpat works with predict_EA technique", {
  result <- assess_cvpat(model_one, technique = predict_EA,
                          nboot = 20, seed = 123, cores = 1)
  expect_type(result, "list")
  expect_true(all(is.finite(as.numeric(result$CVPAT_compare_LM))))
})

test_that("assess_cvpat_compare works with predict_EA technique", {
  result <- assess_cvpat_compare(model_one, model_two,
                                  technique = predict_EA,
                                  nboot = 20, seed = 123, cores = 1)
  expect_true(inherits(result, "matrix"))
  expect_true(all(is.finite(as.numeric(result))))
})

test_that("assess_cvpat works with testtype greater", {
  result <- assess_cvpat(model_one, testtype = "greater",
                          nboot = 20, seed = 123, cores = 1)
  expect_type(result, "list")
  expect_true(all(is.finite(as.numeric(result$CVPAT_compare_LM))))
})

test_that("assess_cvpat_compare works with testtype greater", {
  result <- assess_cvpat_compare(model_one, model_two,
                                  testtype = "greater",
                                  nboot = 20, seed = 123, cores = 1)
  expect_true(inherits(result, "matrix"))
  expect_true(all(is.finite(as.numeric(result))))
})
