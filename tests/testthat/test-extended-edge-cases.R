# ============================================================================
# Extended edge-case tests
# ============================================================================
# These tests estimate multiple model variants with full cross-validation.
# They are slow (~2-5 min) and skipped during R CMD check / CRAN.
# Run locally with: devtools::test(filter = "extended")

library(seminr)

skip_on_cran()
skip_if(identical(Sys.getenv("NOT_CRAN"), ""), "Skipping extended tests (set NOT_CRAN=true to run)")

# ============================================================================
# COA edge cases
# ============================================================================

test_that("assess_coa works with two-stage moderation", {
  mm_mod <- constructs(
    composite("Image", multi_items("IMAG", 1:5)),
    composite("Value", multi_items("PERV", 1:2)),
    composite("Satisfaction", multi_items("CUSA", 1:3)),
    interaction_term(iv = "Image", moderator = "Value", method = two_stage)
  )
  sm_mod <- relationships(
    paths(from = c("Image", "Value", "Image*Value"), to = "Satisfaction")
  )
  model_mod <- estimate_pls(data = mobi, measurement_model = mm_mod,
                            structural_model = sm_mod)
  result <- assess_coa(model_mod, focal_construct = "Satisfaction",
                        noFolds = 10, cores = 1, seed = 123)

  expect_s3_class(result, "coa_analysis")
  expect_true(result$predictive_deviance$overfit_ratio >= 0)
})

test_that("assess_coa works with mediated moderation", {
  mm_medmod <- constructs(
    composite("Image", multi_items("IMAG", 1:5)),
    composite("Value", multi_items("PERV", 1:2)),
    composite("Satisfaction", multi_items("CUSA", 1:3)),
    composite("Loyalty", multi_items("CUSL", 1:3)),
    interaction_term(iv = "Image", moderator = "Value", method = two_stage)
  )
  sm_medmod <- relationships(
    paths(from = c("Image", "Value", "Image*Value"), to = "Satisfaction"),
    paths(from = "Satisfaction", to = "Loyalty")
  )
  model_medmod <- estimate_pls(data = mobi, measurement_model = mm_medmod,
                               structural_model = sm_medmod)
  result <- assess_coa(model_medmod, focal_construct = "Loyalty",
                        noFolds = 10, cores = 1, seed = 123)

  expect_s3_class(result, "coa_analysis")
  expect_s3_class(result$deviance_tree, "coa_dtree")
})

test_that("assess_coa works with serial mediation", {
  mm_ser <- constructs(
    composite("Image", multi_items("IMAG", 1:5)),
    composite("Quality", multi_items("PERQ", 1:7)),
    composite("Satisfaction", multi_items("CUSA", 1:3)),
    composite("Loyalty", multi_items("CUSL", 1:3))
  )
  sm_ser <- relationships(
    paths(from = "Image", to = "Quality"),
    paths(from = c("Image", "Quality"), to = "Satisfaction"),
    paths(from = "Satisfaction", to = "Loyalty")
  )
  model_ser <- estimate_pls(data = mobi, measurement_model = mm_ser,
                            structural_model = sm_ser)
  result <- assess_coa(model_ser, focal_construct = "Loyalty",
                        noFolds = 10, cores = 1, seed = 123)

  expect_s3_class(result, "coa_analysis")
  expect_equal(length(result$predictive_deviance$PD), nrow(mobi))
})

test_that("assess_coa works with mode_B composites", {
  mm_b <- constructs(
    composite("Image", multi_items("IMAG", 1:5), weights = mode_B),
    composite("Value", multi_items("PERV", 1:2), weights = mode_B),
    composite("Satisfaction", multi_items("CUSA", 1:3)),
    composite("Loyalty", multi_items("CUSL", 1:3))
  )
  sm_b <- relationships(
    paths(from = c("Image", "Value"), to = "Satisfaction"),
    paths(from = "Satisfaction", to = "Loyalty")
  )
  model_b <- estimate_pls(data = mobi, measurement_model = mm_b,
                           structural_model = sm_b)
  result <- assess_coa(model_b, focal_construct = "Loyalty",
                        noFolds = 10, cores = 1, seed = 123)

  expect_s3_class(result, "coa_analysis")
})

test_that("assess_coa works with reflective constructs (PLSc)", {
  mm_ref <- constructs(
    reflective("Image", multi_items("IMAG", 1:5)),
    reflective("Value", multi_items("PERV", 1:2)),
    reflective("Satisfaction", multi_items("CUSA", 1:3)),
    reflective("Loyalty", multi_items("CUSL", 1:3))
  )
  sm_ref <- relationships(
    paths(from = c("Image", "Value"), to = "Satisfaction"),
    paths(from = "Satisfaction", to = "Loyalty")
  )
  model_ref <- estimate_pls(data = mobi, measurement_model = mm_ref,
                             structural_model = sm_ref)
  result <- assess_coa(model_ref, focal_construct = "Loyalty",
                        noFolds = 10, cores = 1, seed = 123)

  expect_s3_class(result, "coa_analysis")
})

test_that("assess_coa works with predict_EA on mediated moderation", {
  mm_ea <- constructs(
    composite("Image", multi_items("IMAG", 1:5)),
    composite("Value", multi_items("PERV", 1:2)),
    composite("Satisfaction", multi_items("CUSA", 1:3)),
    composite("Loyalty", multi_items("CUSL", 1:3)),
    interaction_term(iv = "Image", moderator = "Value", method = two_stage)
  )
  sm_ea <- relationships(
    paths(from = c("Image", "Value", "Image*Value"), to = "Satisfaction"),
    paths(from = "Satisfaction", to = "Loyalty")
  )
  model_ea <- estimate_pls(data = mobi, measurement_model = mm_ea,
                            structural_model = sm_ea)
  result <- assess_coa(model_ea, focal_construct = "Loyalty",
                        technique = predict_EA,
                        noFolds = 10, cores = 1, seed = 123)

  expect_s3_class(result, "coa_analysis")
})

# ============================================================================
# NCA extended tests (with significance testing)
# ============================================================================

test_that("NCA reproducibility with significance testing", {
  skip_if_not_installed("NCA")

  mobi_mm <- constructs(
    composite("Image",        multi_items("IMAG", 1:5)),
    composite("Value",        multi_items("PERV", 1:2)),
    composite("Satisfaction", multi_items("CUSA", 1:3))
  )
  mobi_sm <- relationships(
    paths(from = c("Image", "Value"), to = "Satisfaction")
  )
  pls_model <- estimate_pls(data = mobi, measurement_model = mobi_mm,
                             structural_model = mobi_sm)

  r1 <- assess_nca(pls_model, target = "Satisfaction",
                    test.rep = 200, seed = 42)
  r2 <- assess_nca(pls_model, target = "Satisfaction",
                    test.rep = 200, seed = 42)

  expect_equal(r1$effect_sizes, r2$effect_sizes)
})

# ============================================================================
# CVPAT extended: technique comparison and testtype variants
# ============================================================================

test_that("assess_cvpat EA vs DA give different results", {
  corp_rep_mm <- constructs(
    composite("COMP", multi_items("comp_", 1:3)),
    composite("LIKE", multi_items("like_", 1:3)),
    composite("CUSA", single_item("cusa")),
    composite("CUSL", multi_items("cusl_", 1:3))
  )
  sm <- relationships(
    paths(from = c("COMP", "LIKE"), to = "CUSA"),
    paths(from = "CUSA", to = "CUSL")
  )
  model <- estimate_pls(data = corp_rep_data, measurement_model = corp_rep_mm,
                         structural_model = sm, missing = mean_replacement,
                         missing_value = "-99")

  result_ea <- assess_cvpat(model, technique = predict_EA,
                             nboot = 20, seed = 123, cores = 1)
  result_da <- assess_cvpat(model, technique = predict_DA,
                             nboot = 20, seed = 123, cores = 1)

  expect_false(identical(
    as.numeric(result_ea$CVPAT_compare_LM),
    as.numeric(result_da$CVPAT_compare_LM)
  ))
})
