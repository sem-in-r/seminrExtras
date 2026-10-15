# The book PRINTS the argument defaults, so a default is a published number.
#
# Tables 6.2 and 6.3 (pp. 142-143) list the arguments of assess_cvpat() and
# assess_cvpat_compare() with their defaults. Chapter 4 relies on a default it
# does NOT print: it calls congruence_test() with no `reliability`. Changing
# any of these silently falsifies the book. Fast, no model fitting.

test_that("congruence_test() defaults still match what chapter 4 relies on", {
  skip_unless_textbook()

  f <- formals(congruence_test)
  # THE regression of 1.0.3: chapter 4 takes whatever is first here.
  expect_identical(eval(f$reliability)[1], "rhoA")
  expect_identical(eval(f$nboot), 2000)
  expect_identical(eval(f$seed), 123)
  expect_identical(eval(f$threshold), 1)
})

test_that("assess_cvpat() defaults match Table 6.2 as printed", {
  skip_unless_textbook()

  f <- formals(assess_cvpat)
  expect_identical(eval(f$testtype), "two.sided")   # "(default)"
  expect_identical(eval(f$nboot), 2000)             # "defaults to 2000"
  expect_identical(eval(f$seed), 123)               # "defaults to 123"
  expect_identical(deparse(f$technique), "predict_DA")  # "(default)"
})

test_that("assess_cvpat_compare() defaults match Table 6.3 as printed", {
  skip_unless_textbook()

  f <- formals(assess_cvpat_compare)
  expect_identical(eval(f$testtype), "two.sided")
  expect_identical(eval(f$nboot), 2000)
  expect_identical(eval(f$seed), 123)
  expect_identical(deparse(f$technique), "predict_DA")
})

test_that("noFolds has no default, so the book must always pass it", {
  skip_unless_textbook()

  # Tables 6.2/6.3 describe noFolds as "Number of folds for k-fold
  # cross-validation" and quote no default. There is none: it is NULL, which
  # the underlying predict_pls() reads as LEAVE-ONE-OUT -- 344 refits on the
  # corporate reputation data, not 10. That is a performance cliff and a
  # different estimator, so any book code omitting it is a defect.
  expect_null(eval(formals(assess_cvpat)$noFolds))
  expect_null(eval(formals(assess_cvpat_compare)$noFolds))

  demo_dir <- system.file("demo", package = "seminrExtras")
  if (!nzchar(demo_dir)) demo_dir <- testthat::test_path("..", "..", "demo")
  skip_if_not(dir.exists(demo_dir), "demo/ not available")

  offenders <- character()
  for (f in list.files(demo_dir, pattern = "\\.R$", full.names = TRUE)) {
    src <- paste(readLines(f, warn = FALSE), collapse = "\n")
    for (fn in c("assess_cvpat", "assess_cvpat_compare")) {
      for (call in extract_calls(src, fn)) {
        if (!grepl("noFolds", call)) offenders <- c(offenders, paste0(basename(f), ": ", fn))
      }
    }
  }
  expect_equal(offenders, character(),
               label = "demo calls to assess_cvpat*() that omit noFolds (=> LOOCV)")
})
