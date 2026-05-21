# CRAN submission: seminrExtras 1.0.1

## Release summary

Patch release. seminrExtras 1.0.1 removes the last `seminr:::` calls from
the test suite by routing five test lines to seminrExtras's own local
helpers of the same name. This keeps the package compatible with the
current CRAN seminr (2.4.2) and unblocks the forthcoming seminr 2.5.0,
which refactors (and renames) those internal helpers.

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
