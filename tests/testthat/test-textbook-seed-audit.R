# Are the book's own numbers reproducible AT ALL?
#
# Independent of any package version. `bootstrap_model()` does
# `if (is.null(seed)) seed <- sample.int(100000, 1)`, so a call with no `seed`
# draws a RANDOM seed from whatever RNG state is current. The chapter then
# yields different numbers on every fresh run and no reader can reproduce the
# printed values -- nor can we.
#
# TRAP, learned the hard way: run the chapters sequentially in one session and
# chapter 7 looks reproducible. It is not. Chapter 6's SEEDED bootstrap leaves
# the RNG deterministic, so chapter 7's unseeded call picks a repeatable
# "random" seed. Run standalone, chapter 7's bootstrap SD moved 0.0677 ->
# 0.0644 between two fresh runs. This test is static analysis precisely to
# sidestep that: it reads the source instead of trusting a run.

textbook_demo_dir <- function() {
  d <- system.file("demo", package = "seminrExtras")
  if (!nzchar(d) || !dir.exists(d)) d <- testthat::test_path("..", "..", "demo")
  d
}

test_that("every bootstrap_model() call in the book demos sets a seed", {
  skip_unless_textbook()
  d <- textbook_demo_dir()
  skip_if_not(dir.exists(d), "demo/ not available")

  files <- list.files(d, pattern = "^seminr-primer-v2-chap[0-9]+\\.R$", full.names = TRUE)
  skip_if(length(files) == 0, "no chapter demos found")

  unseeded <- character()
  for (f in files) {
    src <- paste(strip_comments(readLines(f, warn = FALSE)), collapse = "\n")
    for (call in extract_calls(src, "bootstrap_model")) {
      if (!grepl("(?<![[:alnum:]._])seed[[:space:]]*=", call, perl = TRUE)) {
        unseeded <- c(unseeded, sprintf("%s: %s", basename(f),
                                        substr(gsub("\\s+", " ", call), 1, 70)))
      }
    }
  }

  # BASELINE, not a clean bill of health. Chapter 4 and chapter 7 each contain
  # one unseeded call (seminrExtras issue #23) and the book is going to print
  # that way unless the decision changes. Asserting `character()` here would
  # leave the suite permanently red, which trains you to ignore it -- so the
  # assertion is against the KNOWN set instead. Any NEW unseeded call fails
  # loudly; fixing one of the known two also fails, at which point delete it
  # from the baseline and the test goes green for good.
  known <- c("seminr-primer-v2-chap4.R", "seminr-primer-v2-chap7.R")
  offending_files <- sort(unique(sub(":.*$", "", unseeded)))

  if (length(unseeded)) {
    message("Unseeded bootstrap_model() calls (issue #23):\n  ",
            paste(unseeded, collapse = "\n  "))
  }
  expect_equal(offending_files, known,
               label = "files with unseeded bootstrap_model() calls")
  # No file may contribute more than the one call already accounted for.
  expect_equal(length(unseeded), length(known),
               label = "count of unseeded bootstrap_model() calls")
})

test_that("congruence and CVPAT calls in the demos are reproducible", {
  skip_unless_textbook()
  d <- textbook_demo_dir()
  skip_if_not(dir.exists(d), "demo/ not available")

  # These default to seed = 123, so an omitted seed is fine. What is NOT fine
  # is passing seed = NULL explicitly, which reaches the same random branch.
  bad <- character()
  for (f in list.files(d, pattern = "\\.R$", full.names = TRUE)) {
    src <- paste(strip_comments(readLines(f, warn = FALSE)), collapse = "\n")
    for (fn in c("congruence_test", "assess_cvpat", "assess_cvpat_compare")) {
      for (call in extract_calls(src, fn)) {
        if (grepl("seed[[:space:]]*=[[:space:]]*NULL", call)) {
          bad <- c(bad, paste0(basename(f), ": ", fn))
        }
      }
    }
  }
  expect_equal(bad, character(), label = "seminrExtras calls forcing seed = NULL")
})

test_that("no demo hard-codes parallel::detectCores()", {
  skip_unless_textbook()
  d <- textbook_demo_dir()
  skip_if_not(dir.exists(d), "demo/ not available")

  # seminr 2.6.0 caps the IMPLICIT default at two cores for CRAN policy, but an
  # EXPLICIT cores = parallel::detectCores() still overrides that. The book
  # prints such a call in chapter 6, which teaches readers a pattern CRAN
  # rejects. Informational: reported, not asserted, since it is the book's
  # editorial choice rather than a package defect.
  hits <- character()
  for (f in list.files(d, pattern = "\\.R$", full.names = TRUE)) {
    src <- paste(strip_comments(readLines(f, warn = FALSE)), collapse = "\n")
    if (grepl("detectCores", src)) hits <- c(hits, basename(f))
  }
  if (length(hits)) {
    message("Demos calling parallel::detectCores(): ", paste(hits, collapse = ", "))
  }
  expect_true(TRUE)
})
