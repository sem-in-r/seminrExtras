# Gate for the textbook regression tests.
#
# These reproduce numbers PRINTED IN THE BOOK. They are deliberately excluded
# from every automated context — CRAN, GitHub Actions, cron — because they are
# slow (CVPAT runs 10 x 10-fold cross-validation over 2000 bootstrap samples),
# they assert against page proofs rather than against a specification, and a
# failure needs a human to decide whether the CODE regressed or the FIGURE is
# simply due for regeneration. An automated red build cannot make that call.
#
# To run them:
#   Sys.setenv(SEMINR_TEXTBOOK_TESTS = "true")
#   devtools::test(filter = "textbook")

skip_unless_textbook <- function() {
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  if (!isTRUE(as.logical(Sys.getenv("SEMINR_TEXTBOOK_TESTS", "false")))) {
    testthat::skip("Textbook regression tests: set SEMINR_TEXTBOOK_TESTS=true to run.")
  }
}

# The book's simple corporate reputation model (chapters 3-5).
textbook_model_simple <- function() {
  mm <- seminr::constructs(
    seminr::composite("COMP", seminr::multi_items("comp_", 1:3)),
    seminr::composite("LIKE", seminr::multi_items("like_", 1:3)),
    seminr::composite("CUSA", seminr::single_item("cusa")),
    seminr::composite("CUSL", seminr::multi_items("cusl_", 1:3)))
  sm <- seminr::relationships(
    seminr::paths(from = c("COMP", "LIKE"), to = c("CUSA", "CUSL")),
    seminr::paths(from = "CUSA", to = "CUSL"))
  seminr::estimate_pls(seminr::corp_rep_data, mm, sm,
                       missing = seminr::mean_replacement, missing_value = "-99")
}

# The book's extended corporate reputation model (chapters 6-7).
textbook_model_extended <- function() {
  mm <- seminr::constructs(
    seminr::composite("QUAL", seminr::multi_items("qual_", 1:8), weights = seminr::mode_B),
    seminr::composite("PERF", seminr::multi_items("perf_", 1:5), weights = seminr::mode_B),
    seminr::composite("CSOR", seminr::multi_items("csor_", 1:5), weights = seminr::mode_B),
    seminr::composite("ATTR", seminr::multi_items("attr_", 1:3), weights = seminr::mode_B),
    seminr::composite("COMP", seminr::multi_items("comp_", 1:3)),
    seminr::composite("LIKE", seminr::multi_items("like_", 1:3)),
    seminr::composite("CUSA", seminr::single_item("cusa")),
    seminr::composite("CUSL", seminr::multi_items("cusl_", 1:3)))
  sm <- seminr::relationships(
    seminr::paths(from = c("QUAL", "PERF", "CSOR", "ATTR"), to = c("COMP", "LIKE")),
    seminr::paths(from = c("COMP", "LIKE"), to = c("CUSA", "CUSL")),
    seminr::paths(from = "CUSA", to = "CUSL"))
  seminr::estimate_pls(seminr::corp_rep_data, mm, sm,
                       missing = seminr::mean_replacement, missing_value = "-99")
}

golden <- function(file) {
  utils::read.csv(testthat::test_path("textbook", file), stringsAsFactors = FALSE)
}

# Extract complete calls to `fn` from R source text, matching parentheses
# properly. A regex or an awk paren count gets this WRONG on nested calls such
# as `cores = parallel::detectCores()`; that bug produced a false seed-audit
# table on 21 Aug 2026 before being caught.
extract_calls <- function(src, fn) {
  pat <- paste0("(?<![[:alnum:]._])", fn, "[[:space:]]*\\(")
  starts <- gregexpr(pat, src, perl = TRUE)[[1]]
  if (starts[1] == -1) return(character())
  lens <- attr(starts, "match.length")
  chars <- strsplit(src, "")[[1]]
  out <- character()
  for (k in seq_along(starts)) {
    i <- starts[k] + lens[k] - 1L   # index of the opening "("
    depth <- 0L; j <- i; n <- length(chars); in_str <- ""
    repeat {
      ch <- chars[j]
      if (nzchar(in_str)) {
        if (ch == in_str && chars[j - 1L] != "\\") in_str <- ""
      } else if (ch %in% c("\"", "'")) {
        in_str <- ch
      } else if (ch == "(") depth <- depth + 1L
      else if (ch == ")") {
        depth <- depth - 1L
        if (depth == 0L) break
      }
      j <- j + 1L
      if (j > n) break
    }
    out <- c(out, paste(chars[starts[k]:min(j, n)], collapse = ""))
  }
  out
}

# Drop comment lines before scanning source, so a commented-out example is not
# reported as a live call. Only whole-line comments are removed; a trailing
# comment after real code is harmless for these scans.
strip_comments <- function(lines) lines[!grepl("^\\s*#", lines)]
