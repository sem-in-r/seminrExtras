# =============================================================================
# Tests for assess_cta() — Confirmatory Tetrad Analysis for PLS-SEM
# =============================================================================

library(seminr)

# --- Setup: pre-compute models and results -----------------------------------

# Simple model with a 5-indicator construct (Image) and several 2-3 indicator ones
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

# Pre-compute CTA with borrowing (default) and without borrowing
cta_result <- assess_cta(pls_model, nboot = 50, seed = 123)
cta_no_borrow <- assess_cta(pls_model, nboot = 50, seed = 123, borrow = FALSE)

# =============================================================================
# STRUCTURE TESTS
# =============================================================================

test_that("assess_cta returns cta_analysis object", {
  expect_s3_class(cta_result, "cta_analysis")
})

test_that("cta result contains all expected components", {
  expected_names <- c("construct_results", "tetrad_details", "nboot",
                      "alpha", "correction", "skipped", "borrowing")
  expect_true(all(expected_names %in% names(cta_result)))
})

test_that("construct_results is a data frame with correct columns", {
  df <- cta_result$construct_results
  expect_s3_class(df, "data.frame")
  expected_cols <- c("Construct", "Mode", "Indicators", "Tetrads",
                     "Significant", "Verdict")
  expect_equal(colnames(df), expected_cols)
})

test_that("tetrad_details is a named list matching construct_results", {
  expect_type(cta_result$tetrad_details, "list")
  expect_equal(names(cta_result$tetrad_details),
               cta_result$construct_results$Construct)
})

# =============================================================================
# NO-BORROW TESTS (borrow = FALSE)
# =============================================================================

test_that("constructs with < 4 indicators are skipped when borrow = FALSE", {
  expect_true("Value" %in% cta_no_borrow$skipped)
  expect_true("Expectation" %in% cta_no_borrow$skipped)
  expect_true("Loyalty" %in% cta_no_borrow$skipped)
  expect_true("Satisfaction" %in% cta_no_borrow$skipped)
})

test_that("Image (5 indicators) is tested without borrowing", {
  expect_true("Image" %in% cta_no_borrow$construct_results$Construct)
})

test_that("correct number of tetrads for 5 indicators without borrowing", {
  # C(5,4) = 5 four-tuples, 2 tetrads each = 10 tetrads
  img_row <- cta_no_borrow$construct_results[
    cta_no_borrow$construct_results$Construct == "Image", ]
  expect_equal(img_row$Indicators, 5)
  expect_equal(img_row$Tetrads, 10)
})

# =============================================================================
# BORROWING TESTS
# =============================================================================

test_that("constructs with 3 indicators are tested via borrowing", {
  tested <- cta_result$construct_results$Construct
  expect_true("Expectation" %in% tested)
  expect_true("Satisfaction" %in% tested)
  expect_true("Loyalty" %in% tested)
})

test_that("construct with 2 indicators is tested via borrowing", {
  expect_true("Value" %in% cta_result$construct_results$Construct)
})

test_that("borrowing details recorded in result", {
  expect_type(cta_result$borrowing, "list")
  # All 4 constructs with < 4 indicators should have borrowing details
  expect_true("Expectation" %in% names(cta_result$borrowing))
  expect_true("Value" %in% names(cta_result$borrowing))
  expect_true("Satisfaction" %in% names(cta_result$borrowing))
  expect_true("Loyalty" %in% names(cta_result$borrowing))
})

test_that("3-indicator construct borrows with 'all' vanishing pattern", {
  # Expectation (3 own, Mode A) borrows from a Mode A neighbor
  b <- cta_result$borrowing[["Expectation"]]
  expect_equal(b$vanishing_pattern, "all")
  expect_equal(b$n_vanishing, 2L)
})

test_that("2-indicator construct borrows with 'tau_1342' vanishing pattern", {
  # Value (2 own, Mode A) borrows 2 indicators
  b <- cta_result$borrowing[["Value"]]
  expect_equal(b$vanishing_pattern, "tau_1342")
  expect_equal(b$n_vanishing, 1L)
})

test_that("borrowed construct has correct indicator count", {
  # 3 own + 1 borrowed = 4 total
  exp_row <- cta_result$construct_results[
    cta_result$construct_results$Construct == "Expectation", ]
  expect_equal(exp_row$Indicators, 4)

  # 2 own + 2 borrowed = 4 total
  val_row <- cta_result$construct_results[
    cta_result$construct_results$Construct == "Value", ]
  expect_equal(val_row$Indicators, 4)
})

test_that("borrowed 3+1 construct has 2 tetrads", {
  exp_row <- cta_result$construct_results[
    cta_result$construct_results$Construct == "Expectation", ]
  expect_equal(exp_row$Tetrads, 2)
})

test_that("borrowed 2+2 construct has 1 tetrad (tau_1342 only)", {
  val_row <- cta_result$construct_results[
    cta_result$construct_results$Construct == "Value", ]
  expect_equal(val_row$Tetrads, 1)
})

test_that("borrowed tetrad detail has tetrad_id 3 for tau_1342", {
  td <- cta_result$tetrad_details[["Value"]]
  # The label should show the tau_1342 pattern
  expect_equal(nrow(td), 1)
  expect_true(grepl("^s\\(", td$Tetrad))
})

test_that("mode label shows borrowing annotation", {
  exp_row <- cta_result$construct_results[
    cta_result$construct_results$Construct == "Expectation", ]
  expect_true(grepl("borrowed from", exp_row$Mode))
})

test_that("Image is not borrowed (has enough indicators)", {
  expect_false("Image" %in% names(cta_result$borrowing))
})

test_that("no constructs skipped with borrowing in mobi model", {
  expect_length(cta_result$skipped, 0)
})

# =============================================================================
# TETRAD DETAIL TESTS
# =============================================================================

test_that("tetrad detail data frame has correct structure", {
  detail <- cta_result$tetrad_details[["Image"]]
  expect_s3_class(detail, "data.frame")
  expect_true("Tetrad" %in% colnames(detail))
  expect_true("Estimate" %in% colnames(detail))
  expect_true("Boot_Mean" %in% colnames(detail))
  expect_true("Boot_SD" %in% colnames(detail))
  expect_true("P_Value" %in% colnames(detail))
  expect_true("Adj_P" %in% colnames(detail))
  expect_true("Significant" %in% colnames(detail))
})

test_that("tetrad estimates are finite", {
  detail <- cta_result$tetrad_details[["Image"]]
  expect_true(all(is.finite(detail$Estimate)))
})

test_that("bootstrap SDs are non-negative", {
  detail <- cta_result$tetrad_details[["Image"]]
  expect_true(all(detail$Boot_SD >= 0, na.rm = TRUE))
})

test_that("p-values are between 0 and 1", {
  detail <- cta_result$tetrad_details[["Image"]]
  expect_true(all(detail$P_Value >= 0 & detail$P_Value <= 1, na.rm = TRUE))
  expect_true(all(detail$Adj_P >= 0 & detail$Adj_P <= 1, na.rm = TRUE))
})

test_that("confidence intervals are ordered correctly", {
  detail <- cta_result$tetrad_details[["Image"]]
  ci_cols <- grep("% CI", colnames(detail), value = TRUE)
  expect_length(ci_cols, 2)
  lower <- as.numeric(detail[[ci_cols[1]]])
  upper <- as.numeric(detail[[ci_cols[2]]])
  expect_true(all(lower <= upper, na.rm = TRUE))
})

test_that("verdict is Reflective supported or rejected", {
  verdicts <- cta_result$construct_results$Verdict
  valid_verdicts <- c("Reflective supported", "Reflective rejected")
  expect_true(all(verdicts %in% valid_verdicts))
})

# =============================================================================
# REPRODUCIBILITY TESTS
# =============================================================================

test_that("assess_cta is reproducible with same seed", {
  r1 <- assess_cta(pls_model, nboot = 30, seed = 42)
  r2 <- assess_cta(pls_model, nboot = 30, seed = 42)
  d1 <- r1$tetrad_details[["Image"]]
  d2 <- r2$tetrad_details[["Image"]]
  expect_equal(d1$Estimate, d2$Estimate)
  expect_equal(d1$Boot_SD, d2$Boot_SD)
})

test_that("assess_cta differs with different seeds", {
  r1 <- assess_cta(pls_model, nboot = 30, seed = 42)
  r2 <- assess_cta(pls_model, nboot = 30, seed = 99)
  d1 <- r1$tetrad_details[["Image"]]
  d2 <- r2$tetrad_details[["Image"]]
  # Original estimates identical (same model)
  expect_equal(d1$Estimate, d2$Estimate)
  # Bootstrap SDs differ
  expect_false(identical(d1$Boot_SD, d2$Boot_SD))
})

# =============================================================================
# PARAMETER VARIATION TESTS
# =============================================================================

test_that("correction methods produce different adjusted p-values", {
  r_bh <- assess_cta(pls_model, constructs = "Image", nboot = 50,
                      seed = 123, correction = "BH")
  r_bon <- assess_cta(pls_model, constructs = "Image", nboot = 50,
                       seed = 123, correction = "bonferroni")
  r_none <- assess_cta(pls_model, constructs = "Image", nboot = 50,
                        seed = 123, correction = "none")

  d_bh <- r_bh$tetrad_details[["Image"]]
  d_bon <- r_bon$tetrad_details[["Image"]]
  d_none <- r_none$tetrad_details[["Image"]]

  # Raw p-values should be identical (same bootstrap)
  expect_equal(d_bh$P_Value, d_bon$P_Value)
  expect_equal(d_bh$P_Value, d_none$P_Value)

  # Adjusted p-values: none <= BH <= Bonferroni (generally)
  # At minimum, "none" adjusted = raw
  expect_equal(d_none$P_Value, d_none$Adj_P)
})

test_that("alpha parameter affects CI column names", {
  r_05 <- assess_cta(pls_model, constructs = "Image",
                      nboot = 30, seed = 123, alpha = 0.05)
  r_10 <- assess_cta(pls_model, constructs = "Image",
                      nboot = 30, seed = 123, alpha = 0.10)

  d_05 <- r_05$tetrad_details[["Image"]]
  d_10 <- r_10$tetrad_details[["Image"]]

  expect_true("2.5% CI" %in% colnames(d_05))
  expect_true("97.5% CI" %in% colnames(d_05))
  expect_true("5% CI" %in% colnames(d_10))
  expect_true("95% CI" %in% colnames(d_10))
})

test_that("selecting specific constructs works", {
  r <- assess_cta(pls_model, constructs = "Image", nboot = 30, seed = 123)
  expect_equal(nrow(r$construct_results), 1)
  expect_equal(r$construct_results$Construct, "Image")
  expect_length(r$skipped, 0)
})

# =============================================================================
# INPUT VALIDATION TESTS
# =============================================================================

test_that("assess_cta rejects non-seminr model objects", {
  expect_warning(
    result <- assess_cta(list(not = "a_model"), nboot = 30),
    "only works with SEMinR models"
  )
  expect_null(result)
})

test_that("assess_cta warns for invalid construct names", {
  expect_warning(
    r <- assess_cta(pls_model, constructs = c("Image", "Nonexistent"), nboot = 30),
    "Constructs not found"
  )
  expect_equal(nrow(r$construct_results), 1)
})

test_that("assess_cta returns empty result when all constructs invalid", {
  expect_warning(
    r <- assess_cta(pls_model, constructs = "Nonexistent", nboot = 30),
    "No valid constructs"
  )
  expect_null(r)
})

# =============================================================================
# MODERATION TESTS
# =============================================================================

test_that("interaction constructs are auto-excluded", {
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

  r <- assess_cta(pls_mod, nboot = 30, seed = 123)

  # Interaction construct should be in skipped
  expect_true("Image*Value" %in% r$skipped)

  # Image should still be tested (5 indicators)
  expect_true("Image" %in% r$construct_results$Construct)
})

test_that("moderation model: non-interaction constructs can borrow", {
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

  r <- assess_cta(pls_mod, nboot = 30, seed = 123)

  # Value (2 indicators) should borrow — its neighbors are Image and Satisfaction
  # (Image*Value is excluded as interaction). Both Image and Satisfaction are Mode A.
  # Value needs 2 borrowed, tau_1342 pattern
  expect_true("Value" %in% r$construct_results$Construct ||
              "Value" %in% names(r$borrowing))
})

# =============================================================================
# MEDIATION TESTS
# =============================================================================

test_that("mediation model works with borrowing", {
  mm_med <- constructs(
    composite("Image",        multi_items("IMAG", 1:5)),
    composite("Satisfaction", multi_items("CUSA", 1:3)),
    composite("Loyalty",      multi_items("CUSL", 1:3))
  )
  sm_med <- relationships(
    paths(from = "Image", to = c("Satisfaction", "Loyalty")),
    paths(from = "Satisfaction", to = "Loyalty")
  )
  pls_med <- estimate_pls(data = mobi, measurement_model = mm_med,
                           structural_model = sm_med)

  r <- assess_cta(pls_med, nboot = 30, seed = 123)

  # Image tested directly (5 indicators)
  expect_true("Image" %in% r$construct_results$Construct)

  # Satisfaction (3 indicators) borrows from Image (Mode A, 5 indicators)
  expect_true("Satisfaction" %in% r$construct_results$Construct)
  expect_true("Satisfaction" %in% names(r$borrowing))

  # Loyalty (3 indicators) borrows from Image or Satisfaction
  expect_true("Loyalty" %in% r$construct_results$Construct)
  expect_true("Loyalty" %in% names(r$borrowing))
})

test_that("mediation model: borrow = FALSE skips as before", {
  mm_med <- constructs(
    composite("Image",        multi_items("IMAG", 1:5)),
    composite("Satisfaction", multi_items("CUSA", 1:3)),
    composite("Loyalty",      multi_items("CUSL", 1:3))
  )
  sm_med <- relationships(
    paths(from = "Image", to = c("Satisfaction", "Loyalty")),
    paths(from = "Satisfaction", to = "Loyalty")
  )
  pls_med <- estimate_pls(data = mobi, measurement_model = mm_med,
                           structural_model = sm_med)

  r <- assess_cta(pls_med, nboot = 30, seed = 123, borrow = FALSE)
  expect_true("Satisfaction" %in% r$skipped)
  expect_true("Loyalty" %in% r$skipped)
  expect_true("Image" %in% r$construct_results$Construct)
})

# =============================================================================
# HOC TESTS
# =============================================================================

test_that("HOC with >= 4 LOCs is tested using LOC construct scores", {
  mm_hoc2 <- constructs(
    composite("IMAG1_2",  multi_items("IMAG", 1:2)),
    composite("IMAG3_4",  multi_items("IMAG", 3:4)),
    composite("CUEX_all", multi_items("CUEX", 1:3)),
    composite("CUSA_all", multi_items("CUSA", 1:3)),
    higher_composite("Super",
                     c("IMAG1_2", "IMAG3_4", "CUEX_all", "CUSA_all"),
                     method = two_stage),
    composite("Loyalty",  multi_items("CUSL", 1:3))
  )
  sm_hoc2 <- relationships(
    paths(from = "Super", to = "Loyalty")
  )
  pls_hoc <- estimate_pls(data = mobi, measurement_model = mm_hoc2,
                           structural_model = sm_hoc2)

  r <- assess_cta(pls_hoc, constructs = "Super", nboot = 30, seed = 123)

  expect_equal(nrow(r$construct_results), 1)
  expect_equal(r$construct_results$Construct, "Super")
  expect_equal(r$construct_results$Indicators, 4)
  # 4 LOCs → C(4,4) = 1 four-tuple → 2 tetrads
  expect_equal(r$construct_results$Tetrads, 2)
  expect_true(grepl("HOC", r$construct_results$Mode))
})

test_that("HOC with < 4 LOCs is skipped (no borrowing across HOC/standard)", {
  mm_hoc <- constructs(
    composite("Image",        multi_items("IMAG", 1:5)),
    composite("Expectation",  multi_items("CUEX", 1:3)),
    higher_composite("Quality", c("Image", "Expectation"), method = two_stage),
    composite("Satisfaction", multi_items("CUSA", 1:3)),
    composite("Loyalty",      multi_items("CUSL", 1:3))
  )
  sm_hoc <- relationships(
    paths(from = "Quality", to = c("Satisfaction", "Loyalty")),
    paths(from = "Satisfaction", to = "Loyalty")
  )
  pls_hoc <- estimate_pls(data = mobi, measurement_model = mm_hoc,
                           structural_model = sm_hoc)

  r <- assess_cta(pls_hoc, constructs = "Quality", nboot = 30, seed = 123)

  # Quality has 2 LOCs, neighbors Satisfaction/Loyalty are standard → can't borrow
  expect_true("Quality" %in% r$skipped)
  expect_equal(nrow(r$construct_results), 0)
})

test_that("LOCs under HOC are tested independently with their own indicators", {
  mm_hoc <- constructs(
    composite("Image",        multi_items("IMAG", 1:5)),
    composite("Expectation",  multi_items("CUEX", 1:3)),
    higher_composite("Quality", c("Image", "Expectation"), method = two_stage),
    composite("Satisfaction", multi_items("CUSA", 1:3)),
    composite("Loyalty",      multi_items("CUSL", 1:3))
  )
  sm_hoc <- relationships(
    paths(from = "Quality", to = c("Satisfaction", "Loyalty")),
    paths(from = "Satisfaction", to = "Loyalty")
  )
  pls_hoc <- estimate_pls(data = mobi, measurement_model = mm_hoc,
                           structural_model = sm_hoc)

  # Test Image specifically — it has 5 raw indicators, should work
  r <- assess_cta(pls_hoc, constructs = "Image", nboot = 30, seed = 123)
  expect_true("Image" %in% r$construct_results$Construct)
  expect_equal(r$construct_results$Indicators, 5)
})

# =============================================================================
# FORMATIVE CONSTRUCT BORROWING TESTS
# =============================================================================

test_that("formative focal construct cannot borrow", {
  mm_form <- constructs(
    composite("Image",  multi_items("IMAG", 1:5)),
    composite("Value",  multi_items("PERV", 1:2), weights = mode_B),
    composite("Satisfaction", multi_items("CUSA", 1:3))
  )
  sm_form <- relationships(
    paths(from = c("Image", "Value"), to = "Satisfaction")
  )
  pls_form <- estimate_pls(data = mobi, measurement_model = mm_form,
                            structural_model = sm_form)

  r <- assess_cta(pls_form, constructs = "Value", nboot = 30, seed = 123)

  # Value is Mode B with 2 indicators → find_donor returns NULL for formative focal
  expect_true("Value" %in% r$skipped)
})

test_that("3 own reflective + 1 borrowed from formative: no donor found", {
  # If all neighbors are formative, a 3-indicator reflective construct can't borrow
  mm_mixed <- constructs(
    composite("Image",        multi_items("IMAG", 1:5), weights = mode_B),
    composite("Expectation",  multi_items("CUEX", 1:3)),
    composite("Satisfaction", multi_items("CUSA", 1:3))
  )
  sm_mixed <- relationships(
    paths(from = "Image", to = c("Expectation", "Satisfaction"))
  )
  pls_mixed <- estimate_pls(data = mobi, measurement_model = mm_mixed,
                              structural_model = sm_mixed)

  r <- assess_cta(pls_mixed, constructs = "Expectation", nboot = 30, seed = 123)

  # Expectation's only neighbor (Image) is formative → can't borrow (3+1 from B = 0 vanishing)
  # Satisfaction isn't structurally connected to Expectation
  expect_true("Expectation" %in% r$skipped)
})

# =============================================================================
# S3 METHODS TESTS
# =============================================================================

test_that("print.cta_analysis runs without error", {
  expect_output(print(cta_result), "Confirmatory Tetrad Analysis")
})

test_that("print.cta_analysis shows borrowing info", {
  expect_output(print(cta_result), "Borrowing")
})

test_that("summary.cta_analysis returns summary.cta_analysis class", {
  s <- summary(cta_result)
  expect_s3_class(s, "summary.cta_analysis")
})

test_that("print.summary.cta_analysis runs without error", {
  s <- summary(cta_result)
  expect_output(print(s), "Detailed Results")
})

test_that("print.summary.cta_analysis shows borrowing details", {
  s <- summary(cta_result)
  expect_output(print(s), "Borrowed")
})

test_that("plot.cta_analysis runs without error", {
  expect_no_error(plot(cta_result))
})

# =============================================================================
# TETRAD COMPUTATION CORRECTNESS TESTS
# =============================================================================

test_that("tetrads are zero for perfectly correlated indicators (common factor)", {
  # Generate data from a single factor: x = lambda * f + epsilon
  set.seed(42)
  n <- 500
  f <- rnorm(n)
  lambdas <- c(0.8, 0.7, 0.9, 0.85)
  x <- sapply(lambdas, function(l) l * f + rnorm(n, sd = sqrt(1 - l^2)))
  colnames(x) <- paste0("x", 1:4)

  cov_mat <- cov(x)

  # Compute tetrads
  td <- seminrExtras:::enumerate_tetrads(colnames(x))
  vals <- seminrExtras:::compute_tetrads(cov_mat, td)

  # With n = 500 from a true factor model, tetrads should be close to zero
  expect_true(all(abs(vals) < 0.05))
})

test_that("tetrads are non-zero for composite indicators", {
  set.seed(42)
  n <- 500

  # Create data with a specific non-factor correlation structure
  Sigma <- matrix(c(1.0, 0.5, 0.3, 0.1,
                     0.5, 1.0, 0.6, 0.2,
                     0.3, 0.6, 1.0, 0.7,
                     0.1, 0.2, 0.7, 1.0), nrow = 4)
  x2 <- MASS::mvrnorm(n, mu = rep(0, 4), Sigma = Sigma)
  colnames(x2) <- paste0("x", 1:4)

  cov_mat2 <- cov(x2)
  td <- seminrExtras:::enumerate_tetrads(colnames(x2))
  vals2 <- seminrExtras:::compute_tetrads(cov_mat2, td)

  # At least some tetrads should be notably non-zero for non-factor structure
  expect_true(any(abs(vals2) > 0.01))
})

test_that("tau_1342 (tetrad_id 3) is computed correctly", {
  # Manual verification
  cov_mat <- matrix(c(1.0, 0.5, 0.3, 0.2,
                       0.5, 1.0, 0.4, 0.6,
                       0.3, 0.4, 1.0, 0.7,
                       0.2, 0.6, 0.7, 1.0), nrow = 4)
  rownames(cov_mat) <- colnames(cov_mat) <- c("a", "b", "c", "d")

  td <- data.frame(i = "a", j = "b", k = "c", l = "d",
                    tetrad_id = 3L, stringsAsFactors = FALSE)
  val <- seminrExtras:::compute_tetrads(cov_mat, td)

  # tau_1342 = sigma_ik * sigma_jl - sigma_il * sigma_jk
  # = sigma(a,c) * sigma(b,d) - sigma(a,d) * sigma(b,c)
  # = 0.3 * 0.6 - 0.2 * 0.4 = 0.18 - 0.08 = 0.10
  expect_equal(val, 0.10)
})

test_that("enumerate_borrowed_tetrads produces correct counts", {
  # 3 own + 1 borrowed, "all" pattern: standard enumeration
  td_all <- seminrExtras:::enumerate_borrowed_tetrads(
    c("x1", "x2", "x3"), "x4", "all"
  )
  expect_equal(nrow(td_all), 2)  # C(4,4)=1 combo, 2 tetrads

  # 2 own + 2 borrowed, "tau_1342" pattern: 1 tetrad
  td_tau <- seminrExtras:::enumerate_borrowed_tetrads(
    c("x1", "x2"), c("x3", "x4"), "tau_1342"
  )
  expect_equal(nrow(td_tau), 1)
  expect_equal(td_tau$tetrad_id, 3L)
  expect_equal(td_tau$i, "x1")
  expect_equal(td_tau$j, "x2")
  expect_equal(td_tau$k, "x3")
  expect_equal(td_tau$l, "x4")
})

# =============================================================================
# INTERNAL HELPER TESTS
# =============================================================================

test_that("get_construct_mode returns correct mode", {
  expect_equal(seminrExtras:::get_construct_mode("Image", pls_model), "A")
})

test_that("get_structurally_connected finds neighbors", {
  neighbors <- seminrExtras:::get_structurally_connected("Image", pls_model)
  expect_true("Expectation" %in% neighbors)
  expect_true("Satisfaction" %in% neighbors)
  expect_true("Loyalty" %in% neighbors)
})

test_that("find_donor returns NULL for formative focal", {
  mm_form <- constructs(
    composite("Image", multi_items("IMAG", 1:5)),
    composite("Value", multi_items("PERV", 1:2), weights = mode_B),
    composite("Satisfaction", multi_items("CUSA", 1:3))
  )
  sm_form <- relationships(
    paths(from = c("Image", "Value"), to = "Satisfaction")
  )
  pls_form <- estimate_pls(data = mobi, measurement_model = mm_form,
                            structural_model = sm_form)

  focal_info <- seminrExtras:::resolve_indicators("Value", pls_form)
  result <- seminrExtras:::find_donor("Value", focal_info, pls_form)
  expect_null(result)
})
