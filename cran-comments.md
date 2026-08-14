# CRAN submission: seminrExtras 1.0.2 (resubmission)

## Release summary

Patch release. seminrExtras 1.0.2 fixes a pair-labelling bug in
`congruence_test()` (coefficient values were always correct; only the
construct-pair labels could be transposed for models with four or more
constructs).

This is a resubmission. The previous 1.0.2 submission was rejected by the
incoming auto-check because the test suite ran longer than 10 minutes on
r-devel-windows-x86_64 (Overall checktime NOTE; tests ~11 min).

## Test-timing fix (addressing the rejection)

The test suite has been sped up substantially with no loss of coverage:

* The Cross-Validated Predictive Ability Test (CVPAT) tests dominated the
  runtime. They now use a small number of cross-validation folds
  (`noFolds = 3`) and few bootstrap subsamples (`nboot = 10`-`30`) instead of
  the previous 10-fold CV with up to 2000 bootstrap subsamples. The CVPAT
  fixture was regenerated to match.
* The PLS-POS comparison object is now pre-computed once and reused across
  the comparison tests instead of being re-estimated in each `test_that()`
  block, and the reused fixtures use fewer random starts.

On the maintainer's machine the full `testthat` suite now runs in ~24 s
(previously ~125 s). The slowest extended/edge-case tests were already
gated behind `skip_on_cran()`.

No user-facing API changes. No dependency floor changes.

## R CMD check results

0 errors | 0 warnings | 0 notes

Local check: `devtools::check()` on macOS, R 4.x — clean.
Tested against CRAN seminr 2.4.2.

## Reverse dependencies

None on CRAN.

## Test environments

* local macOS (R release)
* GitHub Actions: macOS and Ubuntu (R release and devel)
