# CRAN submission: seminrExtras 1.0.4

## Release summary

Minor release: one new function, a correction to repeated cross-validation in
CVPAT, and demo fixes. Every change is described in NEWS.md.

**New function.** `congruence()` reports the congruence coefficient rc for
every pair of constructs in an estimated seminr model, as a lower-triangular
table laid out the way seminr prints HTMT. rc follows Franke, Sarstedt and
Danks (2021, Eq. 2). It is deterministic, a single pass over the estimated
model with no resampling, and reports the coefficients as an effect size only:
no significance test for rc has been validated for PLS-SEM, and `summary()`
prints that caveat. `congruence_test()` is unchanged.

**Results that change.** `assess_cvpat()` and `assess_cvpat_compare()` with
`reps` greater than 1 now draw new random folds in every repetition, compute
each observation's loss within each repetition, and average the losses. With
seminr 2.5.0 and earlier every repetition reused the same folds, so `reps` had
no effect. Results with `reps = NULL` or `reps = 1` are unchanged.

**Demos.** Two demos stopped at their first deliberate error (`demo()` stops at
the first error); those calls are now wrapped in `try()`. The chapter 2 demo no
longer calls `install.packages()`.

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
