# Do we still produce the numbers PRINTED IN THE BOOK?
#
# This is the test the old harness was missing. `test_02_version_equivalence.R`
# in ~/Documents/Projects/textbook_regression_check asks only "does version A
# equal version B", which is self-referential: two versions can agree with each
# other and both disagree with the page proof. The oracle here is external —
# values transcribed by hand from the proofs (see textbook/README.md).
#
# Motivating failure: seminrExtras 1.0.3 changed the `congruence_test()`
# default from rho_C to rho_A. Chapter 4 calls it WITHOUT `reliability`, so
# every coefficient in Fig. 4.11 moved. Nothing caught it.

# A value agrees with the book if it rounds to what the book prints.
expect_prints_as <- function(actual, expected, label, dp = 3) {
  slack <- 0.5 * 10^(-dp) + 1e-9
  testthat::expect(
    all(abs(actual - expected) <= slack),
    sprintf("%s differs from the printed value.\n  printed: %s\n  computed: %s\n  max|diff|: %.6f",
            label, paste(format(expected, nsmall = dp), collapse = " "),
            paste(format(round(actual, dp), nsmall = dp), collapse = " "),
            max(abs(actual - expected))))
  invisible(actual)
}

norm_pair <- function(x) gsub("\\s+", " ", trimws(x))

test_that("Fig. 4.11 congruence coefficients match the book", {
  skip_unless_textbook()

  g <- golden("fig_4_11_congruence.csv")
  # Exactly the book's chapter 4 call: no `reliability`, so the DEFAULT is
  # under test as much as the arithmetic is.
  res <- congruence_test(textbook_model_simple(), alpha = 0.10)$results

  expect_equal(norm_pair(rownames(res)), g$pair)
  expect_prints_as(res[, "Original Est."], g$original_est, "Fig. 4.11 congruence coefficients")
})

test_that("Fig. 4.11 pair labels are not transposed", {
  skip_unless_textbook()

  # Guards the combn()/upper.tri() ordering bug fixed in 1.0.2, independently
  # of any golden file: rebuild the matrix name-indexed and compare.
  res <- congruence_test(textbook_model_simple(), alpha = 0.10, nboot = 0)$results
  labs <- strsplit(norm_pair(rownames(res)), " -> ", fixed = TRUE)

  expect_true(all(lengths(labs) == 2))
  # COMP -> CUSL and LIKE -> CUSA are the pair that swapped in <= 1.0.1.
  expect_true("COMP -> CUSL" %in% norm_pair(rownames(res)))
  expect_true("LIKE -> CUSA" %in% norm_pair(rownames(res)))
  comp_cusl <- res[norm_pair(rownames(res)) == "COMP -> CUSL", "Original Est."]
  like_cusa <- res[norm_pair(rownames(res)) == "LIKE -> CUSA", "Original Est."]
  # Under the correct labelling LIKE->CUSA exceeds COMP->CUSL; the printed
  # proof has it the other way round, which is how the transposition shows.
  expect_gt(like_cusa, comp_cusl)
})

test_that("Fig. 6.9 CVPAT results match the book", {
  skip_unless_textbook()
  skip_if(Sys.getenv("SEMINR_TEXTBOOK_SLOW", "true") == "false",
          "Slow: 10 reps x 10-fold CV over 2000 bootstrap samples.")

  g <- golden("fig_6_9_cvpat.csv")
  res <- assess_cvpat(seminr_model = textbook_model_extended(),
                      testtype = "greater", nboot = 2000, seed = 123,
                      technique = predict_DA, noFolds = 10, reps = 10)

  for (b in c("LM", "IA")) {
    gb <- g[g$benchmark == b, ]
    tab <- res[[paste0("CVPAT_compare_", b)]]
    expect_false(is.null(tab), label = paste0("CVPAT_compare_", b, " is present"))
    expect_equal(rownames(tab), gb$construct)
    expect_prints_as(tab[, 1], gb$pls_loss,   paste("Fig. 6.9", b, "PLS Loss"))
    expect_prints_as(tab[, 2], gb$bench_loss, paste("Fig. 6.9", b, "benchmark Loss"))
    expect_prints_as(tab[, 3], gb$diff,       paste("Fig. 6.9", b, "Diff"))
  }
})

test_that("the sentences on p. 143 that quote Fig. 6.9 are still true", {
  skip_unless_textbook()

  # The prose makes claims, not just a screenshot. Assert the CLAIMS, which
  # survive small numerical drift that would fail the exact comparison above.
  g <- golden("fig_6_9_cvpat.csv")
  lm <- g[g$benchmark == "LM", ]; ia <- g[g$benchmark == "IA", ]

  # "negative loss differences ... for all endogenous constructs" (vs IA)
  expect_true(all(ia$diff < 0))
  # "all loss differences are significant" (vs IA)
  expect_true(all(ia$boot_p < 0.05))
  # "COMP and CUSA are not better predicted ... than the LM"
  expect_gt(lm$boot_p[lm$construct == "COMP"], 0.05)
  expect_gt(lm$boot_p[lm$construct == "CUSA"], 0.05)
  # "the overall model also yields better predictive power than the LM"
  expect_lt(lm$diff[lm$construct == "Overall"], 0)
  expect_lt(lm$boot_p[lm$construct == "Overall"], 0.001)
})
