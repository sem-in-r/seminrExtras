library(seminr)

# Repeated cross-validation in CVPAT: each repetition draws new folds, each
# observation's loss is computed within every repetition, and the losses (not
# the predictions) are averaged over repetitions.

reps_mm <- constructs(
  composite("COMP", multi_items("comp_", 1:3)),
  composite("LIKE", multi_items("like_", 1:3)),
  composite("CUSA", single_item("cusa")),
  composite("CUSL", multi_items("cusl_", 1:3))
)
reps_sm_one <- relationships(
  paths(from = c("COMP", "LIKE"), to = c("CUSA")),
  paths(from = c("CUSA"), to = c("CUSL")))
reps_sm_two <- relationships(
  paths(from = c("COMP", "LIKE"), to = c("CUSA", "CUSL")),
  paths(from = c("CUSA"), to = c("CUSL")))
reps_model_one <- estimate_pls(data = corp_rep_data, measurement_model = reps_mm,
                               structural_model = reps_sm_one,
                               missing = mean_replacement, missing_value = "-99")
reps_model_two <- estimate_pls(data = corp_rep_data, measurement_model = reps_mm,
                               structural_model = reps_sm_two,
                               missing = mean_replacement, missing_value = "-99")

# Per-observation construct losses (mean squared error over the construct's
# items), computed by hand from one set of residuals
by_hand_losses <- function(model, residuals, constructs) {
  mm <- model$mmMatrix
  sapply(constructs, function(lv) {
    items <- mm[mm[, "construct"] == lv, "measurement"]
    rowMeans(as.matrix(residuals)[, items, drop = FALSE]^2)
  })
}

test_that("reps = 1 and reps = NULL give identical results", {
  one  <- assess_cvpat(reps_model_one, nboot = 20, seed = 7, noFolds = 3, reps = 1)
  none <- assess_cvpat(reps_model_one, nboot = 20, seed = 7, noFolds = 3, reps = NULL)
  expect_identical(one, none)

  one  <- assess_cvpat_compare(reps_model_one, reps_model_two, nboot = 20, seed = 7,
                               noFolds = 3, reps = 1)
  none <- assess_cvpat_compare(reps_model_one, reps_model_two, nboot = 20, seed = 7,
                               noFolds = 3, reps = NULL)
  expect_identical(one, none)
})

test_that("assess_cvpat averages each observation's loss over repetitions", {
  seed <- 42
  res <- assess_cvpat(reps_model_one, technique = predict_DA, nboot = 20, seed = seed,
                      noFolds = 3, reps = 2)

  # The same folds, drawn by hand: one predict_pls() run per repetition
  set.seed(seed)
  runs <- lapply(1:2, function(i) {
    suppressWarnings(predict_pls(reps_model_one, technique = predict_DA,
                                 noFolds = 3, reps = NULL))
  })
  endo <- c("CUSA", "CUSL")
  pls <- lapply(runs, function(r) by_hand_losses(reps_model_one, r$items$PLS_out_of_sample_residuals, endo))
  lm  <- lapply(runs, function(r) by_hand_losses(reps_model_one, r$items$lm_out_of_sample_residuals, endo))
  pls_mean <- (pls[[1]] + pls[[2]]) / 2
  lm_mean  <- (lm[[1]] + lm[[2]]) / 2

  # assess_cvpat() returns its table at 7 significant digits
  tab <- res$CVPAT_compare_LM
  expect_equal(unname(tab[endo, "PLS Loss"]), unname(colMeans(pls_mean)), tolerance = 1e-6)
  expect_equal(unname(tab[endo, "LM Loss"]),  unname(colMeans(lm_mean)),  tolerance = 1e-6)
  expect_equal(unname(tab["Overall", "PLS Loss"]), mean(rowMeans(pls_mean)), tolerance = 1e-6)
  expect_equal(unname(tab["Overall", "LM Loss"]),  mean(rowMeans(lm_mean)),  tolerance = 1e-6)

  # Not the loss of the averaged predictions: for squared error that is
  # smaller whenever the repetitions' predictions differ
  avg_resid <- (as.matrix(runs[[1]]$items$lm_out_of_sample_residuals) +
                  as.matrix(runs[[2]]$items$lm_out_of_sample_residuals)) / 2
  lm_of_avg <- by_hand_losses(reps_model_one, avg_resid, endo)
  expect_true(all(colMeans(lm_mean) - colMeans(lm_of_avg) > 1e-4))
  expect_true(all(abs(unname(tab[endo, "LM Loss"]) - unname(colMeans(lm_of_avg))) > 1e-4))
})

test_that("assess_cvpat_compare averages each observation's loss over repetitions", {
  seed <- 42
  res <- assess_cvpat_compare(reps_model_one, reps_model_two, technique = predict_DA,
                              nboot = 20, seed = seed, noFolds = 3, reps = 2)

  # Established model's repetitions first, then the alternative model's
  set.seed(seed)
  runs_one <- lapply(1:2, function(i) predict_pls(reps_model_one, technique = predict_DA,
                                                  noFolds = 3, reps = NULL))
  runs_two <- lapply(1:2, function(i) predict_pls(reps_model_two, technique = predict_DA,
                                                  noFolds = 3, reps = NULL))
  endo <- c("CUSA", "CUSL")
  mean_loss <- function(model, runs) {
    l <- lapply(runs, function(r) by_hand_losses(model, r$items$PLS_out_of_sample_residuals, endo))
    (l[[1]] + l[[2]]) / 2
  }
  one <- mean_loss(reps_model_one, runs_one)
  two <- mean_loss(reps_model_two, runs_two)

  expect_equal(unname(res[endo, "Base Model Loss"]), unname(colMeans(one)), tolerance = 1e-10)
  expect_equal(unname(res[endo, "Alt Model Loss"]),  unname(colMeans(two)), tolerance = 1e-10)
  expect_equal(unname(res["Overall", "Diff"]),
               mean(rowMeans(one)) - mean(rowMeans(two)), tolerance = 1e-10)
})

test_that("invalid reps are rejected", {
  expect_error(assess_cvpat(reps_model_one, nboot = 20, noFolds = 3, reps = 0),
               "reps must be")
  expect_error(assess_cvpat(reps_model_one, nboot = 20, noFolds = 3, reps = 2.5),
               "reps must be")
  expect_error(assess_cvpat_compare(reps_model_one, reps_model_two, nboot = 20,
                                    noFolds = 3, reps = -1),
               "reps must be")
})
