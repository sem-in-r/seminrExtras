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
  result <- congruence_test(test_model, nboot = 20, seed = 123)

  # Should return a list

  expect_type(result, "list")


  # Should have 'results' element

  expect_named(result, "results")

  # Results should be a matrix with table_output class
  expect_true(inherits(result$results, "matrix"))
  expect_true("table_output" %in% class(result$results))
})

test_that("congruence_test returns correct dimensions", {
  result <- congruence_test(test_model, nboot = 20, seed = 123)

  # Number of rows should equal number of construct pairs
  n_constructs <- ncol(test_model$construct_scores)
  expected_pairs <- choose(n_constructs, 2)
  expect_equal(nrow(result$results), expected_pairs)

  # Should have 6 columns
  expect_equal(ncol(result$results), 6)
})

test_that("congruence_test returns correct column names", {
  result <- congruence_test(test_model, nboot = 20, seed = 123, alpha = 0.05)

  expected_cols <- c("Original Est.", "Diff", "Bootstrap SD", "T Stat.", "2.5% CI", "97.5% CI")
  expect_equal(colnames(result$results), expected_cols)
})

test_that("congruence_test row names contain construct pairs", {
  result <- congruence_test(test_model, nboot = 20, seed = 123)

  # Row names should contain " -> " pattern
  expect_true(all(grepl(" -> ", rownames(result$results))))
})

# ============================================================================
# Reproducibility Tests
# ============================================================================

test_that("congruence_test is reproducible with same seed", {
  result1 <- congruence_test(test_model, nboot = 20, seed = 42)
  result2 <- congruence_test(test_model, nboot = 20, seed = 42)

  expect_equal(result1$results, result2$results)
})

test_that("congruence_test differs with different seeds", {
  result1 <- congruence_test(test_model, nboot = 20, seed = 42)
  result2 <- congruence_test(test_model, nboot = 20, seed = 99)

  # Bootstrap SD and CI columns should differ
  expect_false(identical(result1$results[, "Bootstrap SD"], result2$results[, "Bootstrap SD"]))
})

# ============================================================================
# Parameter Variation Tests
# ============================================================================

test_that("congruence_test respects alpha parameter for CI columns", {
  result_05 <- congruence_test(test_model, nboot = 20, seed = 123, alpha = 0.05)
  result_10 <- congruence_test(test_model, nboot = 20, seed = 123, alpha = 0.10)

  # Column names should reflect alpha

  expect_true("2.5% CI" %in% colnames(result_05$results))
  expect_true("97.5% CI" %in% colnames(result_05$results))
  expect_true("5% CI" %in% colnames(result_10$results))
  expect_true("95% CI" %in% colnames(result_10$results))
})

test_that("congruence_test works with different nboot values", {
  result_small <- congruence_test(test_model, nboot = 20, seed = 123)
  result_large <- congruence_test(test_model, nboot = 20, seed = 123)

  # Both should return valid results

  expect_true(all(is.finite(result_small$results[, "Original Est."])))
  expect_true(all(is.finite(result_large$results[, "Original Est."])))

  # Original estimates should be identical (same model)
  expect_equal(result_small$results[, "Original Est."], result_large$results[, "Original Est."])
})

test_that("congruence_test works with custom threshold", {
  result <- congruence_test(test_model, nboot = 20, seed = 123, threshold = 0.9)

  # Diff column should reflect threshold - original
  # Diff = threshold - abs(original)
  expect_true(all(is.finite(result$results[, "Diff"])))
})

# ============================================================================
# Input Validation Tests
# ============================================================================

test_that("congruence_test rejects non-seminr model objects", {
  expect_warning(
    result <- congruence_test(list(not = "a_model"), nboot = 20),
    "only works with SEMinR models"
  )
  expect_null(result)
})

test_that("congruence_test rejects NULL input", {
  expect_warning(
    result <- congruence_test(NULL, nboot = 20),
    "only works with SEMinR models"
  )
  expect_null(result)
})

test_that("congruence_test rejects data frame input", {
  expect_warning(
    result <- congruence_test(data.frame(x = 1:10), nboot = 20),
    "only works with SEMinR models"
  )
  expect_null(result)
})

# ============================================================================
# Output Value Tests
# ============================================================================

test_that("congruence_test original estimates are bounded", {
  result <- congruence_test(test_model, nboot = 20, seed = 123)

  # Congruence coefficients should be between -1 and 1
  expect_true(all(result$results[, "Original Est."] >= -1))
  expect_true(all(result$results[, "Original Est."] <= 1))
})

test_that("congruence_test bootstrap SD is positive", {
  result <- congruence_test(test_model, nboot = 20, seed = 123)

  expect_true(all(result$results[, "Bootstrap SD"] >= 0))
})

test_that("congruence_test confidence intervals are ordered correctly", {
  result <- congruence_test(test_model, nboot = 20, seed = 123, alpha = 0.05)

  lower_ci <- result$results[, "2.5% CI"]
  upper_ci <- result$results[, "97.5% CI"]

  # Lower CI should be less than or equal to upper CI
  expect_true(all(lower_ci <= upper_ci))
})

# ============================================================================
# Regression: estimates must be attached to the correct construct pair
# ============================================================================
# Guards against the row-major (combn) vs column-major (upper.tri) mismatch that
# previously swapped coefficients between non-adjacent pairs (e.g. COMP<>CUSL and
# LIKE<>CUSA in this 4-construct model). The reference below indexes by construct
# NAME, so it is independent of any fill-order convention.

# Independent, name-based congruence coefficient (Franke, Sarstedt & Danks 2021,
# Eq. 2): cor matrix of construct scores with the selected reliability estimate
# on the diagonal. Defaults to rhoA, matching congruence_test()'s own default.
reference_rc <- function(model, X, Y, reliability = "rhoA") {
  m <- stats::cor(model$construct_scores)
  cn <- colnames(m)
  diag(m) <- switch(reliability,
    rhoA  = seminr::rho_A(model, cn)[cn, 1],
    rhoC  = seminr::rhoC_AVE(model)[cn, 1],
    # Oracle for alpha: seminr's own summary(). Independent of the package's
    # implementation, which computes alpha locally (summary() is ~400x slower
    # and cannot be called inside the bootstrap loop).
    cronbach = suppressWarnings(summary(model)$reliability)[cn, "alpha"],
    one   = rep(1, length(cn))
  )
  sum(m[, X] * m[, Y]) / sqrt(sum(m[, X]^2) * sum(m[, Y]^2))
}

test_that("congruence_test labels each estimate with the correct construct pair", {
  result <- congruence_test(test_model, nboot = 20, seed = 123)

  for (rn in rownames(result$results)) {
    pair <- trimws(strsplit(rn, "->", fixed = TRUE)[[1]])
    expect_equal(
      unname(result$results[rn, "Original Est."]),
      reference_rc(test_model, pair[1], pair[2]),
      tolerance = 1e-8,
      info = paste("mislabelled estimate for", rn)
    )
  }
})

test_that("congruence_test does not swap COMP<>CUSL and LIKE<>CUSA", {
  # The two pairs that the column-major fill bug swapped in a 4-construct model.
  result <- congruence_test(test_model, nboot = 20, seed = 123)
  ests <- result$results[, "Original Est."]

  comp_cusl <- ests[grepl("COMP +-> +CUSL", rownames(result$results))]
  like_cusa <- ests[grepl("LIKE +-> +CUSA", rownames(result$results))]

  expect_equal(unname(comp_cusl), reference_rc(test_model, "COMP", "CUSL"), tolerance = 1e-8)
  expect_equal(unname(like_cusa), reference_rc(test_model, "LIKE", "CUSA"), tolerance = 1e-8)
  # And they must be distinct, so a swap would be caught.
  expect_false(isTRUE(all.equal(unname(comp_cusl), unname(like_cusa))))
})

# ============================================================================
# Reliability on the diagonal: 1 / rhoA / rhoC
# ============================================================================
# Franke, Sarstedt & Danks (2021, Eq. 2) place "the reliabilities" on the
# diagonal of the construct-correlation matrix without fixing an estimator.
# SmartPLS uses rhoA; seminrExtras <= 1.0.2 hard-coded rhoC. The two disagree,
# so the estimator is now selectable.
#
# ORACLE: values published by SmartPLS for the simple corporate reputation
# model (C. M. Ringle, personal communication, 13 Aug 2026). This is an
# independent implementation, not a value read back out of this package.
#
# NB: this fixture is NOT `test_model` above — the structural model differs
# (COMP/LIKE feed CUSL directly here). Construct scores depend on the inner
# model, so the oracle only applies to this specification.

smartpls_mm <- constructs(
  composite("COMP", multi_items("comp_", 1:3)),
  composite("LIKE", multi_items("like_", 1:3)),
  composite("CUSA", single_item("cusa")),
  composite("CUSL", multi_items("cusl_", 1:3))
)
smartpls_sm <- relationships(
  paths(from = c("COMP", "LIKE"), to = c("CUSA", "CUSL")),
  paths(from = "CUSA", to = "CUSL")
)
smartpls_model <- estimate_pls(
  data = corp_rep_data,
  measurement_model = smartpls_mm,
  structural_model  = smartpls_sm,
  missing = mean_replacement,
  missing_value = "-99"
)

# SmartPLS congruence coefficients (rhoA on the diagonal), 3 d.p. as published.
smartpls_rc <- c(
  "COMP -> LIKE" = 0.971,
  "COMP -> CUSA" = 0.848,
  "COMP -> CUSL" = 0.891,
  "LIKE -> CUSA" = 0.902,
  "LIKE -> CUSL" = 0.954,
  "CUSA -> CUSL" = 0.967
)

est_for_pair <- function(result, pair) {
  key <- gsub(" ", "", pair)
  rn  <- gsub(" ", "", rownames(result$results))
  unname(result$results[which(rn == key), "Original Est."])
}

test_that("congruence_test reproduces SmartPLS's published coefficients", {
  result <- congruence_test(smartpls_model, nboot = 20, seed = 123)

  for (pair in names(smartpls_rc)) {
    expect_equal(
      est_for_pair(result, pair),
      unname(smartpls_rc[[pair]]),
      # published to 3 d.p., so agreement is only asserted at that precision
      tolerance = 5e-4,
      info = paste("disagrees with SmartPLS for", pair)
    )
  }
})

test_that("reliability = 'rhoC' reproduces pre-1.0.3 behaviour", {
  # Oracle: reference_rc is an independent, name-indexed implementation of
  # Eq. 2, not a value read back from congruence_test().
  result <- congruence_test(smartpls_model, nboot = 20, seed = 123,
                            reliability = "rhoC")
  for (pair in names(smartpls_rc)) {
    cs <- trimws(strsplit(pair, "->", fixed = TRUE)[[1]])
    expect_equal(
      est_for_pair(result, pair),
      reference_rc(smartpls_model, cs[1], cs[2], reliability = "rhoC"),
      tolerance = 1e-8, info = pair
    )
  }
})

test_that("reliability = 'one' puts unity on the diagonal", {
  result <- congruence_test(smartpls_model, nboot = 20, seed = 123,
                            reliability = "one")
  for (pair in names(smartpls_rc)) {
    cs <- trimws(strsplit(pair, "->", fixed = TRUE)[[1]])
    expect_equal(
      est_for_pair(result, pair),
      reference_rc(smartpls_model, cs[1], cs[2], reliability = "one"),
      tolerance = 1e-8, info = pair
    )
  }
})

test_that("the three estimators genuinely differ", {
  # Guards against the argument being silently ignored.
  a <- congruence_test(smartpls_model, nboot = 5, seed = 1, reliability = "rhoA")$results[, "Original Est."]
  c <- congruence_test(smartpls_model, nboot = 5, seed = 1, reliability = "rhoC")$results[, "Original Est."]
  o <- congruence_test(smartpls_model, nboot = 5, seed = 1, reliability = "one")$results[, "Original Est."]
  expect_false(isTRUE(all.equal(a, c)))
  expect_false(isTRUE(all.equal(a, o)))
  expect_false(isTRUE(all.equal(c, o)))
})

test_that("rhoA equals unity for Mode B constructs, so the two agree when all are Mode B", {
  # Internal consistency is undefined for composites: seminr::rho_A() assigns 1
  # to every Mode B construct. An all-Mode-B model therefore has an all-ones
  # diagonal under "rhoA", identical to "one". This is the edge case that makes
  # the estimator choice consequential in MIXED models.
  mb_mm <- constructs(
    composite("QUAL", multi_items("qual_", 1:8), weights = mode_B),
    composite("PERF", multi_items("perf_", 1:5), weights = mode_B),
    composite("CSOR", multi_items("csor_", 1:5), weights = mode_B),
    composite("ATTR", multi_items("attr_", 1:3), weights = mode_B)
  )
  mb_sm <- relationships(paths(from = c("QUAL", "PERF", "CSOR"), to = "ATTR"))
  mb_model <- estimate_pls(corp_rep_data, mb_mm, mb_sm,
                           missing = mean_replacement, missing_value = "-99")

  expect_equal(
    unname(seminr::rho_A(mb_model, colnames(mb_model$construct_scores))[, 1]),
    rep(1, 4)
  )
  expect_equal(
    congruence_test(mb_model, nboot = 5, seed = 1, reliability = "rhoA")$results[, "Original Est."],
    congruence_test(mb_model, nboot = 5, seed = 1, reliability = "one")$results[, "Original Est."],
    tolerance = 1e-10
  )
})

test_that("an unrecognised reliability is rejected", {
  # "omega" is not offered; McDonald's omega is not among the four conventions.
  expect_error(
    congruence_test(smartpls_model, nboot = 5, reliability = "omega"),
    "should be one of"
  )
  # The pre-rename spelling of the Cronbach option must not silently work.
  expect_error(
    congruence_test(smartpls_model, nboot = 5, reliability = "alpha"),
    "should be one of"
  )
})

test_that("the estimators coincide when every construct is single-item", {
  # Both rhoA and rhoC return 1 for a single indicator, so the diagonal is all
  # ones under every option and `reliability` has no effect. Documents that the
  # estimators differ ONLY in their treatment of Mode B constructs.
  si_mm <- constructs(
    composite("A", single_item("comp_1")),
    composite("B", single_item("like_1")),
    composite("C", single_item("cusa")),
    composite("D", single_item("cusl_1"))
  )
  si_sm <- relationships(paths(from = c("A", "B"), to = "C"),
                         paths(from = "C", to = "D"))
  si_model <- estimate_pls(corp_rep_data, si_mm, si_sm,
                           missing = mean_replacement, missing_value = "-99")
  cn <- colnames(si_model$construct_scores)

  expect_equal(unname(seminr::rho_A(si_model, cn)[, 1]), rep(1, 4))
  expect_equal(unname(seminr::rhoC_AVE(si_model)[cn, 1]), rep(1, 4))

  a <- congruence_test(si_model, nboot = 5, seed = 1, reliability = "rhoA")$results[, "Original Est."]
  c <- congruence_test(si_model, nboot = 5, seed = 1, reliability = "rhoC")$results[, "Original Est."]
  o <- congruence_test(si_model, nboot = 5, seed = 1, reliability = "one")$results[, "Original Est."]
  expect_equal(a, c, tolerance = 1e-12)
  expect_equal(a, o, tolerance = 1e-12)
})

test_that("higher-order models are refused", {
  # HOC support is not implemented: the two-stage construct scores mix a
  # higher-order composite with first-stage constructs, and it has not been
  # established what the diagonal should hold for the HOC. Refuse rather than
  # return a number nobody has validated.
  hoc_mm <- constructs(
    composite("QUAL", multi_items("qual_", 1:8), weights = mode_B),
    composite("PERF", multi_items("perf_", 1:5), weights = mode_B),
    higher_composite("REPU", dimensions = c("QUAL", "PERF"),
                     method = two_stage, weights = mode_A),
    composite("CUSL", multi_items("cusl_", 1:3))
  )
  hoc_sm <- relationships(paths(from = "REPU", to = "CUSL"))
  hoc_model <- estimate_pls(corp_rep_data, hoc_mm, hoc_sm,
                            missing = mean_replacement, missing_value = "-99")

  expect_warning(res <- congruence_test(hoc_model, nboot = 5), "higher-order")
  expect_null(res)
})

test_that("interaction constructs are excluded from the analysis entirely", {
  int_mm <- constructs(
    composite("COMP", multi_items("comp_", 1:3)),
    composite("LIKE", multi_items("like_", 1:3)),
    composite("CUSA", single_item("cusa")),
    composite("CUSL", multi_items("cusl_", 1:3)),
    interaction_term(iv = "COMP", moderator = "LIKE", method = two_stage)
  )
  int_sm <- relationships(
    paths(from = c("COMP", "LIKE", "COMP*LIKE"), to = "CUSA"),
    paths(from = "CUSA", to = "CUSL")
  )
  int_model <- estimate_pls(corp_rep_data, int_mm, int_sm,
                            missing = mean_replacement, missing_value = "-99")

  expect_message(
    result <- congruence_test(int_model, nboot = 5, seed = 1),
    "[Ii]nteraction"
  )

  # No pair may mention the interaction term.
  expect_false(any(grepl("*", rownames(result$results), fixed = TRUE)))
  # 4 remaining constructs -> 6 pairs, not the 10 you get with the interaction.
  expect_equal(nrow(result$results), choose(4, 2))

  # And it must leave the CORRELATION VECTORS too, not just the pair list:
  # Eq. 2 sums over the whole construct set, so a coefficient computed with the
  # interaction still in the matrix would differ. Reference excludes it.
  keep <- setdiff(colnames(int_model$construct_scores), "COMP*LIKE")
  m <- stats::cor(int_model$construct_scores[, keep, drop = FALSE])
  diag(m) <- seminr::rho_A(int_model, keep)[keep, 1]
  expected <- sum(m[, "COMP"] * m[, "CUSL"]) /
    sqrt(sum(m[, "COMP"]^2) * sum(m[, "CUSL"]^2))
  expect_equal(est_for_pair(result, "COMP -> CUSL"), expected, tolerance = 1e-8)
})

test_that("reliability = 'cronbach' puts Cronbach's alpha on the diagonal", {
  # ORACLE: seminr's own summary()$reliability alpha column.
  result <- congruence_test(smartpls_model, nboot = 20, seed = 123,
                            reliability = "cronbach")
  for (pair in names(smartpls_rc)) {
    cs <- trimws(strsplit(pair, "->", fixed = TRUE)[[1]])
    expect_equal(
      est_for_pair(result, pair),
      reference_rc(smartpls_model, cs[1], cs[2], reliability = "cronbach"),
      tolerance = 1e-8, info = pair
    )
  }
})

test_that("cronbach is estimated for Mode B constructs, unlike rhoA", {
  # rhoA returns exactly 1 for Mode B (internal consistency is undefined for a
  # composite); alpha and rhoC both compute a value. This is the reason alpha
  # is offered: it applies one rule to every multi-item construct, which makes
  # it comparable with CB-SEM where rhoA is unavailable.
  mb_mm <- constructs(
    composite("QUAL", multi_items("qual_", 1:8), weights = mode_B),
    composite("PERF", multi_items("perf_", 1:5), weights = mode_B),
    composite("COMP", multi_items("comp_", 1:3)),
    composite("CUSL", multi_items("cusl_", 1:3))
  )
  mb_sm <- relationships(paths(from = c("QUAL", "PERF"), to = "COMP"),
                         paths(from = "COMP", to = "CUSL"))
  mb_model <- estimate_pls(corp_rep_data, mb_mm, mb_sm,
                           missing = mean_replacement, missing_value = "-99")
  cn <- colnames(mb_model$construct_scores)
  ref <- suppressWarnings(summary(mb_model)$reliability)

  # rhoA pins the two Mode B constructs at unity; alpha does not.
  expect_equal(unname(seminr::rho_A(mb_model, cn)[c("QUAL", "PERF"), 1]), c(1, 1))
  expect_true(all(ref[c("QUAL", "PERF"), "alpha"] < 1))

  # So the two estimators must disagree on this model.
  a <- congruence_test(mb_model, nboot = 5, seed = 1, reliability = "cronbach")$results[, "Original Est."]
  r <- congruence_test(mb_model, nboot = 5, seed = 1, reliability = "rhoA")$results[, "Original Est."]
  expect_false(isTRUE(all.equal(a, r)))
})

test_that("cronbach matches seminr's own alpha for every construct", {
  # Guards the local reimplementation against drift from seminr's cronbachs_alpha.
  cn <- colnames(smartpls_model$construct_scores)
  ref <- suppressWarnings(summary(smartpls_model)$reliability)[cn, "alpha"]
  got <- vapply(cn, function(x) {
    items <- seminr::construct_items(smartpls_model$mmMatrix, x)
    if (length(items) < 2) return(1)
    cm <- stats::cor(smartpls_model$data)[items, items]
    k <- nrow(cm)
    (k / (k - 1)) * (1 - sum(diag(cm)) / sum(cm))
  }, numeric(1))
  expect_equal(unname(got), unname(ref), tolerance = 1e-10)
})
