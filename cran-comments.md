# CRAN submission: seminrExtras 1.0.3

## Release summary

This release supplies the working maintainer address CRAN requested, and adds a
user-selectable reliability estimate to `congruence_test()`.

**Maintainer address.** CRAN reported that mail to the previous maintainer
address bounced and asked for a new version carrying a working one. The address
is now `seminrgroup@gmail.com`. The `Maintainer` field continues to name an
individual, `Nicholas Patrick Danks <seminrgroup@gmail.com>`; only the address
has changed.

**New argument.** `congruence_test()` gains `reliability`, selecting what is
placed on the diagonal of the construct-correlation matrix: `"rhoA"` (new
default), `"rhoC"`, `"cronbach"` or `"one"`. Franke, Sarstedt and Danks (2021,
Eq. 2) specify "the reliabilities" without fixing an estimator, so all four are
in specification. The default moves from `"rhoC"` to `"rhoA"` to agree with
SmartPLS; `reliability = "rhoC"` reproduces results from 1.0.2 and earlier.

**Behaviour changes.** `congruence_test()` now excludes interaction constructs
(their measurement is determined by the product method rather than by theory)
and refuses higher-order models with a warning rather than returning values
whose interpretation has not been established. Both are documented in NEWS.md.

## Test environments

- local: macOS 15.5 (arm64), R 4.6.0
- GitHub Actions: macOS-latest (R release, R devel)
- GitHub Actions: ubuntu-latest (R release, R devel)
- win-builder: R-devel

## R CMD check results

0 errors | 0 warnings | 0 notes

Full test suite: 797 tests passing, 0 failures.

## Downstream dependencies

There are no reverse dependencies on CRAN for this package.

## Notes for the reviewer

The package documentation now uses roxygen2 markdown (`Roxygen: list(markdown =
TRUE)` in DESCRIPTION). This was previously absent, so markdown written in the
roxygen blocks was passing through to the Rd files literally. Enabling it
regenerates the manual pages correctly; the change touches many files under
`man/` but alters rendering only, not content.
