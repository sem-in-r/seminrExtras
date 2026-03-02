# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

seminrExtras is a supplementary R package that extends SEMinR (Structural Equation Modeling in R) with advanced analysis tools. It is NOT standalone—it requires seminr >= 2.4.0.

**Primary features:**
- Cross-Validated Predictive Ability Test (CVPAT) - `assess_cvpat()` and `assess_cvpat_compare()`
- Congruence testing - `congruence_test()`
- Demo files for the PLS-SEM in R workbook (Hair et al., 2026)

## Build and Development Commands

```r
# Generate documentation (NAMESPACE and man/ files)
roxygen2::roxygenise()

# Run tests
devtools::test()

# Full package check (CRAN-style)
devtools::check()
# or
rcmdcheck::rcmdcheck(args = c("--no-manual", "--as-cran"))

# Install from source
devtools::install()
```

## Code Architecture

### R/ Directory (3 files)

- **feature_cvpat.R** (~530 LOC) - Core CVPAT implementation
  - `assess_cvpat()` - Compare PLS model against Linear Model (LM) and Indicator Average (IA) benchmarks
  - `assess_cvpat_compare()` - Compare two PLS models' predictive performance
  - Uses k-fold cross-validation with bootstrap significance testing
  - Helper functions: `cvpat_per_construct()`, `lv_loss()`, `overall_loss()`, `bootstrap_cvpat()`

- **feature_congruence.R** (~140 LOC) - Congruence coefficient testing
  - `congruence_test()` - Bootstrapped congruence testing for construct validity

- **helpers.R** - Import declarations for roxygen2

### Key Patterns

- Functions accept seminr model objects from `estimate_pls()`
- Results use custom `"table_output"` class with `comment()` for metadata
- Bootstrap resampling via seminr's internal `rerun()` function
- Accesses seminr internals via `seminr:::` (e.g., `seminr:::rhoC_AVE()`)

### Tests

Tests use testthat 3.0+ with fixtures in `tests/fixtures/`. The helper file `tests/testthat/helper-global.R` defines the `test_folder` path.

### Demo Files

Access demos with: `demo("seminr-pls-cvpat", package = "seminrExtras")`

Available: seminr-help-debugging, seminr-pls-cvpat, seminr-primer-v2-chap2 through chap8

## CI/CD

GitHub Actions workflow (`.github/workflows/rcmdcheck.yml`) runs on macOS and Ubuntu (release + devel). Branches ending with `_noci` skip CI.
