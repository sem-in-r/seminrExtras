library(seminr)

mobi_mm <- constructs(
  composite("Image",        multi_items("IMAG", 1:5)),
  composite("Expectation",  multi_items("CUEX", 1:3)),
  composite("Value",        multi_items("PERV", 1:2), weights = mode_B),
  composite("Satisfaction", multi_items("CUSA", 1:3)),
  reflective("Loyalty",     multi_items("CUSL", 1:3))
)
mobi_sm <- relationships(
  paths(from = c("Image", "Expectation"), to = c("Value", "Satisfaction")),
  paths(from = "Value", to = "Satisfaction"),
  paths(from = "Satisfaction", to = "Loyalty")
)
mobi_pls <- suppressMessages(estimate_pls(mobi, mobi_mm, mobi_sm))
mobi_grp <- ifelse(mobi$CUEX1 > stats::median(mobi$CUEX1), "high", "low")

test_that("the internal PLS core reproduces seminr's weights (Mode A, Mode B, reflective; both schemes)", {
  for (iw in list(path_weighting, path_factorial)) {
    m <- suppressMessages(estimate_pls(mobi, mobi_mm, mobi_sm, inner_weights = iw))
    sp <- micom_spec(m, "test")
    W <- micom_pls_core(stats::cor(sp$data), sp$blocks, sp$inner, sp$modes, sp$scheme)
    ow <- m$outer_weights[sp$indicators, sp$constructs]
    for (j in seq_along(sp$constructs)) {
      b <- sp$blocks[[j]]
      w <- W[j, b]
      if (sum(w * ow[b, j]) < 0) w <- -w
      expect_equal(unname(w), unname(ow[b, j]), tolerance = 1e-6)
    }
  }
})

test_that("step-2 c and the step-3 statistics match cSEM::testMICOM", {
  # Reference values from cSEM 0.6.1.9000, testMICOM() on the same model (all
  # composites in Mode A, path weighting) and the same median split of CUEX1.
  mm <- constructs(composite("Image", multi_items("IMAG", 1:5)), composite("Expectation", multi_items("CUEX", 1:3)),
                   composite("Satisfaction", multi_items("CUSA", 1:3)), composite("Loyalty", multi_items("CUSL", 1:3)))
  sm <- relationships(paths(from = c("Image", "Expectation"), to = "Satisfaction"), paths(from = "Satisfaction", to = "Loyalty"))
  m <- suppressMessages(estimate_pls(mobi, mm, sm))
  r <- assess_micom(m, mobi_grp, nperm = 9, seed = 1)
  expect_equal(r$groups, c("high", "low"))
  expect_equal(r$step2$c, c(0.985361, 0.965723, 0.998843, 0.967719), tolerance = 1e-5)
  expect_equal(r$step3$mean_diff, c(0.603450, 1.343900, 0.599691, 0.437276), tolerance = 1e-5)
  expect_equal(r$step3$log_var_ratio, c(0.125523, -0.090378, 0.016500, 0.041257), tolerance = 1e-5)
})

test_that("assess_micom returns one row per construct with valid p-values and a verdict", {
  r <- assess_micom(mobi_pls, mobi_grp, nperm = 19, seed = 1)
  expect_s3_class(r, "micom_analysis")
  expect_equal(r$step2$construct, mobi_pls$constructs)
  expect_equal(nrow(r$step3), length(mobi_pls$constructs))
  ps <- c(r$step2$p_value, r$step3$mean_p, r$step3$var_p)
  expect_true(all(ps >= 0 & ps <= 1))
  expect_true(all(r$step3$invariance %in% c("full", "partial", "none")))
  expect_false("p_studentized" %in% c(names(r$step2), names(r$step3)))
  expect_equal(r$n, c(sum(mobi_grp == "high"), sum(mobi_grp == "low")))
  expect_output(print(r), "Step 2")
})

test_that("results are reproducible with a seed and the group order follows the factor levels", {
  a <- assess_micom(mobi_pls, mobi_grp, nperm = 19, seed = 42)
  b <- assess_micom(mobi_pls, mobi_grp, nperm = 19, seed = 42)
  expect_identical(a$step2, b$step2)
  rev <- assess_micom(mobi_pls, factor(mobi_grp, levels = c("low", "high")), nperm = 19, seed = 42)
  expect_equal(rev$groups, c("low", "high"))
  expect_equal(rev$step3$mean_diff, -a$step3$mean_diff)
  expect_equal(rev$step2$c, a$step2$c)
})

test_that("a grouping column of the model's data can be named", {
  mm <- constructs(composite("COMP", multi_items("comp_", 1:3)), composite("LIKE", multi_items("like_", 1:3)),
                   composite("CUSL", multi_items("cusl_", 1:3)))
  sm <- relationships(paths(from = c("COMP", "LIKE"), to = "CUSL"))
  m <- suppressMessages(estimate_pls(corp_rep_data, mm, sm, missing = mean_replacement, missing_value = "-99"))
  r <- assess_micom(m, "servicetype", nperm = 9, seed = 1)
  expect_equal(r$groups, c("1", "2"))
})

test_that("invalid input is refused", {
  expect_error(assess_micom(mobi_pls, rep(1:3, length.out = nrow(mobi))), "exactly two groups")
  expect_error(assess_micom(mobi_pls, mobi_grp[-1]), "one value per row")
  expect_error(assess_micom(mobi_pls, "no_such_column"), "not in the model's data")
  g <- mobi_grp; g[1] <- NA
  expect_error(assess_micom(mobi_pls, g), "missing values")
  expect_warning(expect_null(assess_micom(list(), mobi_grp)), "SEMinR")
  mm_int <- constructs(composite("Image", multi_items("IMAG", 1:3)), composite("Expectation", multi_items("CUEX", 1:3)),
                       interaction_term(iv = "Image", moderator = "Expectation", method = two_stage),
                       composite("Loyalty", multi_items("CUSL", 1:3)))
  sm_int <- relationships(paths(from = c("Image", "Expectation", "Image*Expectation"), to = "Loyalty"))
  m_int <- suppressMessages(estimate_pls(mobi, mm_int, sm_int))
  expect_warning(expect_null(assess_micom(m_int, mobi_grp)), "interaction")
})

test_that("the research functions are internal, not exported", {
  exports <- getNamespaceExports("seminrExtras")
  expect_true("assess_micom" %in% exports)
  expect_false(any(c("micom_research", "micom_directions", "micom_pls_core") %in% exports))
})

test_that("the internal research test runs and detects a clear difference in composite weights", {
  r0 <- micom_research(mobi_pls, mobi_grp, nperm = 19, jackknife_groups = 5, seed = 1)
  expect_true(is.na(r0$p_studentized[r0$construct == "Value"]))     # Mode B: not tested
  expect_false(any(is.na(r0$p_studentized[r0$mode == "A"])))
  set.seed(7)
  n <- 400
  mk <- function(w) {
    X <- matrix(stats::rnorm(n * 4), n)
    X[, 2] <- X[, 2] + .5 * X[, 1]; X[, 4] <- X[, 4] + .5 * X[, 3]
    eta <- drop(scale(X %*% w))
    data.frame(X, y1 = .6 * eta + stats::rnorm(n), y2 = .6 * eta + stats::rnorm(n))
  }
  d <- rbind(mk(c(1, 1, 1, 1)), mk(c(3, 0, 0, 1)))
  names(d) <- c(paste0("x", 1:4), "y1", "y2")
  mm <- constructs(composite("X", multi_items("x", 1:4)), composite("Y", multi_items("y", 1:2)))
  sm <- relationships(paths(from = "X", to = "Y"))
  m <- suppressMessages(estimate_pls(d, mm, sm))
  r <- micom_research(m, rep(c("a", "b"), each = n), nperm = 49, jackknife_groups = 10, seed = 3)
  expect_lt(r$p_studentized[r$construct == "X"], .05)
  expect_gt(r$angle[r$construct == "X"], 20)
})
