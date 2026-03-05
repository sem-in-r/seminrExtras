library(seminr)

# ============================================================================
# Setup: Estimate models once, reuse across all tests
# ============================================================================

set.seed(123)

# Mediation model: COMP, LIKE -> CUSA -> CUSL
corp_rep_mm <- constructs(
  composite("COMP", multi_items("comp_", 1:3)),
  composite("LIKE", multi_items("like_", 1:3)),
  composite("CUSA", single_item("cusa")),
  composite("CUSL", multi_items("cusl_", 1:3))
)

corp_rep_sm <- relationships(
  paths(from = c("COMP", "LIKE"), to = "CUSA"),
  paths(from = "CUSA", to = "CUSL")
)

pls_model <- estimate_pls(
  data = corp_rep_data,
  measurement_model = corp_rep_mm,
  structural_model  = corp_rep_sm,
  missing = mean_replacement,
  missing_value = "-99"
)

# Pre-compute predictions once (expensive step)
pls_predictions <- predict_pls(
  pls_model,
  technique = predict_DA,
  noFolds = 10,
  reps = 1,
  cores = 1
)

# Pre-compute COA components for reuse
pd_result <- predictive_deviance(
  pls_model,
  focal_construct = "CUSL",
  predict_model = pls_predictions
)

dt_result <- deviance_tree(pd_result)

# ============================================================================
# Unit tests: path_to (tree traversal helper)
# ============================================================================

test_that("path_to computes correct root-to-node path for even node", {
  expect_equal(path_to(468), c(1, 3, 7, 14, 29, 58, 117, 234, 468))
})

test_that("path_to computes correct root-to-node path for odd node", {
  expect_equal(path_to(469), c(1, 3, 7, 14, 29, 58, 117, 234, 469))
})

test_that("path_to returns root for node 1", {
  expect_equal(path_to(1), 1)
})

# ============================================================================
# Unit tests: main_ancestors (deviant group identification)
# ============================================================================

test_that("main_ancestors removes descendant nodes", {
  parent_ids <- c("4", "40", "81", "12", "24", "119", "239", "31")
  result <- main_ancestors(parent_ids)
  expect_equal(result, c("4", "40", "12", "119", "31"))
})

test_that("main_ancestors handles multi-level descendants", {
  parent_ids <- c("2", "12", "24", "204", "26", "29", "117", "469", "15", "31", "62", "124", "63")
  result <- main_ancestors(parent_ids)
  expect_equal(result, c("2", "12", "26", "29", "15"))
})

test_that("main_ancestors returns single node unchanged", {
  expect_equal(main_ancestors("5"), "5")
})

# ============================================================================
# Unit tests: predictive_deviance
# ============================================================================

test_that("predictive_deviance returns correct structure", {
  expect_type(pd_result, "list")
  expect_s3_class(pd_result, "coa_deviance")
  expect_true(all(c("PD", "pd_data", "IS_MSE", "OOS_MSE",
                     "overfit_ratio", "fitted_score",
                     "predicted_score") %in% names(pd_result)))
})

test_that("PD equals fitted minus predicted", {
  expected_pd <- pd_result$fitted_score - pd_result$predicted_score
  expect_equal(pd_result$PD, expected_pd)
})

test_that("PD length matches number of observations", {
  expect_equal(length(pd_result$PD), nrow(pls_model$data))
})

test_that("pd_data contains all construct scores plus PD column", {
  construct_names <- colnames(pls_model$construct_scores)
  expect_true(all(construct_names %in% colnames(pd_result$pd_data)))
  expect_true("PD" %in% colnames(pd_result$pd_data))
})

test_that("overfit_ratio has correct formula", {
  expected_ratio <- (pd_result$OOS_MSE - pd_result$IS_MSE) / pd_result$IS_MSE
  expect_equal(pd_result$overfit_ratio, expected_ratio)
})

test_that("MSE values are non-negative", {
  expect_true(pd_result$IS_MSE >= 0)
  expect_true(pd_result$OOS_MSE >= 0)
})

test_that("predictive_deviance accepts pre-computed predict_model", {
  pd_precomp <- predictive_deviance(
    pls_model,
    focal_construct = "CUSL",
    predict_model = pls_predictions
  )
  pd_fresh <- predictive_deviance(
    pls_model,
    focal_construct = "CUSL",
    technique = predict_DA,
    noFolds = 10,
    reps = 1,
    cores = 1,
    seed = 123
  )

  # Structure should be identical
  expect_s3_class(pd_precomp, "coa_deviance")
  expect_equal(length(pd_precomp$PD), length(pd_fresh$PD))
})

# ============================================================================
# Unit tests: deviance_tree
# ============================================================================

test_that("deviance_tree returns correct structure", {
  expect_type(dt_result, "list")
  expect_s3_class(dt_result, "coa_dtree")
  expect_true(all(c("tree", "deviant_groups", "group_roots",
                     "unique_deviants", "sorted_PD") %in% names(dt_result)))
})

test_that("deviance_tree produces an rpart tree", {
  expect_s3_class(dt_result$tree, "rpart")
})

test_that("deviant_groups are named with uppercase letters", {
  group_names <- names(dt_result$deviant_groups)
  if (length(group_names) > 0) {
    expect_true(all(group_names %in% LETTERS))
  }
})

test_that("group_roots match deviant_groups names", {
  expect_equal(names(dt_result$group_roots), names(dt_result$deviant_groups))
})

test_that("unique_deviants are not in any deviant group", {
  all_grouped <- unlist(dt_result$deviant_groups)
  overlap <- intersect(dt_result$unique_deviants, all_grouped)
  expect_equal(length(overlap), 0)
})

test_that("deviance_tree respects custom bounds", {
  dt_wide <- deviance_tree(pd_result, deviance_bounds = c(0.1, 0.9))
  dt_narrow <- deviance_tree(pd_result, deviance_bounds = c(0.01, 0.99))

  # Wider bounds should produce more deviants
  n_wide <- length(unlist(dt_wide$deviant_groups)) + length(dt_wide$unique_deviants)
  n_narrow <- length(unlist(dt_narrow$deviant_groups)) + length(dt_narrow$unique_deviants)
  expect_true(n_wide >= n_narrow)
})

test_that("sorted_PD is sorted in decreasing order", {
  expect_equal(dt_result$sorted_PD, sort(dt_result$sorted_PD, decreasing = TRUE))
})

# ============================================================================
# Unit tests: unstable_params
# ============================================================================

test_that("unstable_params returns correct structure", {
  skip_if(length(dt_result$deviant_groups) == 0, "No deviant groups found")

  unstable <- unstable_params(
    pls_model,
    deviant_groups = dt_result$deviant_groups,
    params = "path_coef"
  )

  expect_type(unstable, "list")
  expect_s3_class(unstable, "coa_unstable")
  expect_equal(length(unstable), length(dt_result$deviant_groups))
})

test_that("unstable_params diffs have correct dimensions", {
  skip_if(length(dt_result$deviant_groups) == 0, "No deviant groups found")

  unstable <- unstable_params(
    pls_model,
    deviant_groups = dt_result$deviant_groups,
    params = "path_coef"
  )

  first_diff <- unstable[[1]]$param_diffs$path_coef
  expect_equal(dim(first_diff), dim(pls_model$path_coef))
})

test_that("unstable_params works with multiple param types", {
  skip_if(length(dt_result$deviant_groups) == 0, "No deviant groups found")

  unstable <- unstable_params(
    pls_model,
    deviant_groups = dt_result$deviant_groups,
    params = c("path_coef", "outer_weights")
  )

  first_group <- unstable[[1]]$param_diffs
  expect_true("path_coef" %in% names(first_group))
  expect_true("outer_weights" %in% names(first_group))
})

# ============================================================================
# Integration: assess_coa
# ============================================================================

test_that("assess_coa returns correct S3 class and structure", {
  result <- assess_coa(
    pls_model,
    focal_construct = "CUSL",
    noFolds = 10, reps = 1, cores = 1,
    seed = 123
  )

  expect_s3_class(result, "coa_analysis")
  expect_true(all(c("pls_model", "focal_construct", "deviance_bounds",
                     "predictive_deviance", "deviance_tree",
                     "unstable") %in% names(result)))
})

test_that("assess_coa accepts pre-computed predict_model", {
  result <- assess_coa(
    pls_model,
    focal_construct = "CUSL",
    predict_model = pls_predictions,
    seed = 123
  )

  expect_s3_class(result, "coa_analysis")
  expect_s3_class(result$predictive_deviance, "coa_deviance")
  expect_s3_class(result$deviance_tree, "coa_dtree")
})

test_that("assess_coa is reproducible with same seed", {
  r1 <- assess_coa(pls_model, "CUSL",
                    noFolds = 10, reps = 1, cores = 1, seed = 42)
  r2 <- assess_coa(pls_model, "CUSL",
                    noFolds = 10, reps = 1, cores = 1, seed = 42)

  expect_equal(r1$predictive_deviance$PD, r2$predictive_deviance$PD)
})

test_that("assess_coa works with predict_EA technique", {
  result <- assess_coa(
    pls_model,
    focal_construct = "CUSL",
    technique = predict_EA,
    noFolds = 10, reps = 1, cores = 1,
    seed = 123
  )

  expect_s3_class(result, "coa_analysis")
})

# ============================================================================
# Input validation
# ============================================================================

test_that("assess_coa rejects non-seminr model", {
  expect_warning(
    result <- assess_coa(list(not = "a_model"), focal_construct = "X"),
    "only works with SEMinR models"
  )
  expect_null(result)
})

test_that("assess_coa rejects HOC model", {
  hoc_mm <- constructs(
    composite("COMP", multi_items("comp_", 1:3)),
    composite("LIKE", multi_items("like_", 1:3)),
    higher_composite("COKE", c("COMP", "LIKE"), two_stage),
    composite("CUSA", single_item("cusa")),
    composite("CUSL", multi_items("cusl_", 1:3))
  )
  hoc_sm <- relationships(
    paths(from = "COKE", to = "CUSA"),
    paths(from = "CUSA", to = "CUSL")
  )
  hoc_model <- estimate_pls(
    data = corp_rep_data,
    measurement_model = hoc_mm,
    structural_model  = hoc_sm,
    missing = mean_replacement,
    missing_value = "-99"
  )

  expect_warning(
    result <- assess_coa(hoc_model, focal_construct = "CUSL"),
    "higher-order"
  )
  expect_null(result)
})

test_that("assess_coa errors on invalid focal construct", {
  expect_error(
    assess_coa(pls_model, focal_construct = "NONEXISTENT",
               noFolds = 10, reps = 1, cores = 1, seed = 123),
    "focal_construct"
  )
})

# ============================================================================
# Rules extraction, S3 methods, and plots (shared fixture)
# ============================================================================

# Compute once for rules, S3, and plot tests
coa_fixture <- assess_coa(
  pls_model, "CUSL",
  predict_model = pls_predictions, seed = 123
)

test_that("group_rules returns data frame with correct columns", {
  skip_if(length(coa_fixture$deviance_tree$deviant_groups) == 0,
          "No deviant groups found")

  group_name <- names(coa_fixture$deviance_tree$deviant_groups)[1]
  rules <- group_rules(group_name, coa_fixture)

  expect_s3_class(rules, "data.frame")
  expect_true(all(c("construct", "gte", "lt") %in% colnames(rules)))
})

test_that("competes returns data frame with correct columns", {
  skip_if(length(coa_fixture$deviance_tree$deviant_groups) == 0,
          "No deviant groups found")

  group_name <- names(coa_fixture$deviance_tree$deviant_groups)[1]
  node_id <- coa_fixture$deviance_tree$group_roots[[group_name]]
  comp <- competes(node_id, coa_fixture$deviance_tree)

  expect_s3_class(comp, "data.frame")
  expect_true(all(c("criterion", "sign", "value", "improve") %in% colnames(comp)))
})

test_that("competes errors for root node", {
  expect_error(competes(1, coa_fixture$deviance_tree), "root")
})

test_that("print.coa_analysis runs without error", {
  expect_output(print(coa_fixture))
})

test_that("summary.coa_analysis returns summary object", {
  s <- summary(coa_fixture)
  expect_type(s, "list")
  expect_s3_class(s, "summary.coa_analysis")
})

test_that("plot.coa_analysis type='pd' runs without error", {
  expect_no_error(plot(coa_fixture, type = "pd"))
})

test_that("plot.coa_analysis type='groups' runs without error", {
  skip_if(length(coa_fixture$deviance_tree$deviant_groups) == 0,
          "No deviant groups found")
  expect_no_error(plot(coa_fixture, type = "groups"))
})

test_that("plot.coa_analysis type='tree' runs without error", {
  expect_no_error(plot(coa_fixture, type = "tree"))
})

test_that("plot.coa_analysis defaults to 'pd' type", {
  expect_no_error(plot(coa_fixture))
})
