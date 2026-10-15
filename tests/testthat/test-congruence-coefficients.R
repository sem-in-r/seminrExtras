# =============================================================================
# congruence() -- the descriptive congruence coefficient table
# =============================================================================
# congruence() reports rc and nothing else. Its numbers must agree with the
# established Step 5 calculation in congruence_test(), because both claim to
# implement Franke et al. (2021, Eq. 2) and a silent divergence between them
# would be invisible to a user of either.

library(seminr)

congruence_model <- function() {
  mm <- constructs(
    composite("COMP", multi_items("comp_", 1:3)),
    composite("LIKE", multi_items("like_", 1:3)),
    composite("CUSA", single_item("cusa")),
    composite("CUSL", multi_items("cusl_", 1:3))
  )
  sm <- relationships(
    paths(from = c("COMP", "LIKE"), to = c("CUSA", "CUSL")),
    paths(from = "CUSA", to = "CUSL")
  )
  estimate_pls(corp_rep_data, mm, sm,
               missing = mean_replacement, missing_value = "-99")
}

# ============================================================================
# The formula, against values worked by hand
# ============================================================================

test_that("congruence_rc() is the cosine similarity of two columns", {
  m <- matrix(c(3, 4, 4, 3), nrow = 2,
              dimnames = list(NULL, c("X", "Y")))
  # (3*4 + 4*3) / sqrt((9+16) * (16+9)) = 24 / 25
  expect_equal(congruence_rc(m, "X", "Y"), 0.96)

  identical_cols <- matrix(c(0.4, 0.7, 0.4, 0.7), nrow = 2,
                           dimnames = list(NULL, c("X", "Y")))
  expect_equal(congruence_rc(identical_cols, "X", "Y"), 1)

  orthogonal <- matrix(c(1, 0, 0, 1), nrow = 2,
                       dimnames = list(NULL, c("X", "Y")))
  expect_equal(congruence_rc(orthogonal, "X", "Y"), 0)
})

# ============================================================================
# Structure
# ============================================================================

test_that("congruence() returns a symmetric named matrix with 1 on the diagonal", {
  res <- congruence(congruence_model())

  expect_s3_class(res, "congruence_analysis")
  expect_named(res, c("congruence", "constructs", "reliability"))

  rc <- res$congruence
  expect_equal(dim(rc), c(4, 4))
  expect_equal(rownames(rc), c("COMP", "LIKE", "CUSA", "CUSL"))
  expect_equal(colnames(rc), rownames(rc))

  # rc(X, Y) == rc(Y, X), so either indexing order works
  expect_equal(rc, t(rc))
  expect_equal(rc["COMP", "CUSL"], rc["CUSL", "COMP"])
  expect_equal(unname(diag(rc)), rep(1, 4))
})

test_that("congruence coefficients lie within [-1, 1]", {
  rc <- congruence(congruence_model())$congruence
  expect_true(all(rc >= -1 & rc <= 1))
})

# ============================================================================
# Agreement with the established calculation
# ============================================================================

test_that("congruence() reproduces congruence_test()'s point estimates", {
  model <- congruence_model()

  new <- congruence(model)$congruence
  # nboot only drives the inferential columns; Original Est. is deterministic.
  old <- congruence_test(model, nboot = 5, seed = 123)$results[, "Original Est."]

  for (label in names(old)) {
    pair <- trimws(strsplit(label, "->", fixed = TRUE)[[1]])
    expect_equal(new[pair[1], pair[2]], unname(old[[label]]),
                 tolerance = 1e-12,
                 info = paste("pair", label))
  }
})

test_that("congruence() matches an independent implementation of Eq. 2", {
  model <- congruence_model()
  cs <- colnames(model$construct_scores)

  # Built from seminr only, without touching the package internals under test.
  reference <- stats::cor(model$construct_scores[, cs, drop = FALSE])
  diag(reference) <- seminr::rho_A(model, cs)[cs, 1]
  expected <- function(X, Y) {
    sum(reference[, X] * reference[, Y]) /
      sqrt(sum(reference[, X]^2) * sum(reference[, Y]^2))
  }

  rc <- congruence(model)$congruence
  for (pair in utils::combn(cs, 2, simplify = FALSE)) {
    expect_equal(rc[pair[1], pair[2]], expected(pair[1], pair[2]),
                 tolerance = 1e-12)
  }
})

# ============================================================================
# The reliability argument
# ============================================================================

test_that("reliability selects what sits on the diagonal", {
  model <- congruence_model()

  a <- congruence(model, reliability = "rhoA")
  o <- congruence(model, reliability = "one")

  expect_identical(a$reliability, "rhoA")
  expect_identical(o$reliability, "one")
  expect_false(isTRUE(all.equal(a$congruence, o$congruence)))
})

test_that("reliability = 'one' puts unity on the diagonal of the input matrix", {
  model <- congruence_model()
  cs <- colnames(model$construct_scores)

  reference <- stats::cor(model$construct_scores[, cs, drop = FALSE])
  diag(reference) <- 1
  expected <- sum(reference[, "COMP"] * reference[, "LIKE"]) /
    sqrt(sum(reference[, "COMP"]^2) * sum(reference[, "LIKE"]^2))

  rc <- congruence(model, reliability = "one")$congruence
  expect_equal(rc["COMP", "LIKE"], expected, tolerance = 1e-12)
})

test_that("an unknown reliability is refused", {
  expect_error(congruence(congruence_model(), reliability = "rhoZ"))
})

# ============================================================================
# Model types
# ============================================================================

test_that("congruence() refuses higher-order models", {
  mm <- constructs(
    composite("COMP", multi_items("comp_", 1:3)),
    composite("LIKE", multi_items("like_", 1:3)),
    composite("CUSA", single_item("cusa")),
    higher_composite("HOC", c("COMP", "LIKE"), method = two_stage)
  )
  sm <- relationships(paths(from = "HOC", to = "CUSA"))
  hoc_model <- estimate_pls(corp_rep_data, mm, sm,
                            missing = mean_replacement, missing_value = "-99")

  expect_warning(res <- congruence(hoc_model), "higher-order")
  expect_null(res)
})

test_that("congruence() drops interaction constructs and says so", {
  mm <- constructs(
    composite("COMP", multi_items("comp_", 1:3)),
    composite("LIKE", multi_items("like_", 1:3)),
    composite("CUSA", single_item("cusa")),
    interaction_term(iv = "COMP", moderator = "LIKE", method = two_stage)
  )
  sm <- relationships(
    paths(from = c("COMP", "LIKE", "COMP*LIKE"), to = "CUSA")
  )
  int_model <- estimate_pls(corp_rep_data, mm, sm,
                            missing = mean_replacement, missing_value = "-99")

  expect_message(res <- congruence(int_model), "Excluding interaction constructs")
  expect_false("COMP*LIKE" %in% res$constructs)
  expect_equal(dim(res$congruence), c(3, 3))
})

test_that("a non-seminr object is refused", {
  expect_warning(res <- congruence(list(a = 1)), "only works with SEMinR models")
  expect_null(res)
})

# ============================================================================
# Printing -- this table goes into print, so its shape is load-bearing
# ============================================================================

test_that("print() shows the lower triangle, as seminr prints HTMT", {
  res <- congruence(congruence_model())
  out <- capture.output(print(res))

  expect_true(any(grepl("Congruence Coefficients", out)))

  # This header line is printed in Fig. 4.11 of the PLS-SEM R book, so the
  # wording is load-bearing, not incidental.
  expect_true(any(grepl("Calculation uses rhoA on the diagonal", out, fixed = TRUE)))

  # Lower triangle: the COMP row carries none, the CUSL row carries three.
  comp_row <- out[grepl("^COMP", out)]
  cusl_row <- out[grepl("^CUSL", out)]
  expect_length(regmatches(comp_row, gregexpr("[0-9]\\.[0-9]{3}", comp_row))[[1]], 0)
  expect_length(regmatches(cusl_row, gregexpr("[0-9]\\.[0-9]{3}", cusl_row))[[1]], 3)

  # Empty cells are "." exactly as seminr's HTMT table renders them.
  expect_true(any(grepl("\\.", comp_row)))
})

test_that("summary() reports the range and the interpretation caveat", {
  out <- capture.output(print(summary(congruence(congruence_model()))))

  expect_true(any(grepl("Range:", out)))
  expect_true(any(grepl("Pairs: 6", out)))
  expect_true(any(grepl("No significance test", out)))
})
