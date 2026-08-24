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

# Fig. 4.11 (congruence) was REMOVED FROM THE BOOK on 2026-08-24. Sarstedt ruled
# that an rc significance test is work-in-progress and does not belong in a
# textbook; Ringle had independently said the same. The whole congruence
# treatment comes out of chapter 4 — prose, Table 4.1, the worked example, the
# summary-box sentence, the exercise clause, the glossary and index entries, and
# the demo code line. So there is no longer a printed Fig. 4.11 to test against.
#
# The values were verified correct (all 36 cells, max|diff| 0.00049) and are NOT
# thrown away: they moved to test-congruence.R as a characterisation test of
# congruence_test() itself, which is still an exported function.

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

test_that("Fig. 8.7 specific indirect effects match the book", {
  skip_unless_textbook()

  # Both effects at alpha = 0.05 (corrected 21 Aug 2026; the printed figure has
  # the first at alpha = 0.1 while the surrounding text claims 5%).
  g <- golden("fig_8_7_indirect.csv")
  model <- textbook_model_extended()
  boot <- seminr::bootstrap_model(model, nboot = 1000, seed = 123, cores = 2)

  for (i in seq_len(nrow(g))) {
    parts <- strsplit(g$effect[i], "->", fixed = TRUE)[[1]]
    r <- seminr::specific_effect_significance(boot, from = parts[1],
                                              through = parts[2], to = parts[3],
                                              alpha = 0.05)
    expect_prints_as(unname(r[1, "Original Est."]),  g$original_est[i], paste("Fig. 8.7", g$effect[i], "est"))
    expect_prints_as(unname(r[1, "Bootstrap SD"]),   g$boot_sd[i],      paste("Fig. 8.7", g$effect[i], "SD"))
    expect_prints_as(unname(r[1, "T Stat."]),        g$t_stat[i],       paste("Fig. 8.7", g$effect[i], "t"))
    expect_prints_as(unname(r[1, "2.5% CI"]),        g$ci_low[i],       paste("Fig. 8.7", g$effect[i], "CI low"))
    expect_prints_as(unname(r[1, "97.5% CI"]),       g$ci_high[i],      paste("Fig. 8.7", g$effect[i], "CI high"))
    # The claim the text actually makes: significant at the 5% level.
    expect_gt(unname(r[1, "2.5% CI"]), 0)
  }
})
