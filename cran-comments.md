# CRAN submission: seminrExtras 1.0.4

## Release summary

Minor release adding one function.

**New function.** `congruence()` reports the congruence coefficient rc for
every pair of constructs in an estimated seminr model, as a lower-triangular
table laid out the way seminr prints HTMT. rc follows Franke, Sarstedt and
Danks (2021, Eq. 2). It is deterministic, a single pass over the estimated
model with no resampling, and reports the coefficients as an effect size only:
no significance test for rc has been validated for PLS-SEM, and `summary()`
prints that caveat.

**No behaviour changes.** `congruence_test()` is unchanged; its formula,
reliability diagonal and construct set are now internals shared with
`congruence()`, and its existing tests pass untouched. No other function
changed.

## Compatibility with seminr

seminrExtras imports seminr (>= 2.4.0). This version was checked against both
the seminr currently on CRAN (2.5.0) and seminr 2.6.0, which is being submitted
alongside it.

## Test environments

- local: macOS (arm64), R 4.5.1, against seminr 2.5.0 and against seminr 2.6.0
- win-builder: R-release and R-devel

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

seminrExtras has no reverse dependencies on CRAN.
