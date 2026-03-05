# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

seminrExtras is a supplementary R package that extends SEMinR (Structural Equation Modeling in R) with advanced analysis tools. It is NOT standalone—it requires seminr >= 2.4.0.

**Primary features:**
- Cross-Validated Predictive Ability Test (CVPAT) - `assess_cvpat()` and `assess_cvpat_compare()`
- Composite Overfit Analysis (COA) - `coa()`, `predictive_deviance()`, `deviance_tree()`, `unstable_params()`
- Congruence testing - `congruence_test()`
- Demo files for the PLS-SEM in R workbook (Hair et al., 2026)

## Design Philosophy (follows SEMinR)

seminrExtras follows SEMinR's design patterns and philosophy:

### Three-Stage Pipeline: Specify → Estimate → Evaluate

All analysis functions accept estimated seminr model objects. The workflow is:
1. User specifies and estimates a model in seminr (`estimate_pls()`)
2. seminrExtras functions evaluate/analyze the estimated model
3. Results are returned as S3 objects with `print()` and `summary()` methods

### File Naming Conventions

Follow the SEMinR prefix conventions:
| Prefix | Purpose |
|---|---|
| `feature_*.R` | Feature implementations (CVPAT, COA, congruence) |
| `helpers.R` | Internal utilities and import declarations |

### S3 Class System

- All exported analysis functions return S3 classed objects
- Implement `print()` and `summary()` methods for all result classes
- Results use custom `"table_output"` class with `comment()` for metadata where appropriate
- Use `class(obj) <- c("specific_class", class(obj))` pattern for class assignment

### Leveraging seminr Infrastructure

**Always use seminr's existing functions rather than reimplementing:**
- `seminr::estimate_pls()` — for model estimation and re-estimation
- `seminr::predict_pls()` — for k-fold cross-validated prediction (handles moderation, mediation, HOC)
- `seminr:::rerun()` — for re-estimating models with modified parameters
- `seminr:::all_endogenous()` — for extracting endogenous constructs
- `seminr:::items_of_construct()` — for getting measurement items
- `seminr:::rhoC_AVE()` — for reliability/validity metrics

**Model types that must be supported:**
- Simple path models
- Mediation models
- Moderation models (interaction terms)
- Higher-order construct (HOC) models

### Key Patterns

- Functions accept seminr model objects from `estimate_pls()`
- Bootstrap resampling via seminr's internal `rerun()` function
- Accesses seminr internals via `seminr:::` (e.g., `seminr:::rhoC_AVE()`)
- Avoid reimplementing prediction/CV — use `seminr::predict_pls()` instead

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

### Package Management with usethis

**Prefer `usethis::` and `devtools::` functions over manual file editing or bash commands:**

```r
usethis::use_r("filename")              # Create new R source file
usethis::use_test("testname")           # Create new test file
usethis::use_package("pkgname")         # Add to Imports in DESCRIPTION
usethis::use_package("pkgname", "Suggests")  # Add to Suggests
```

## Code Architecture

### R/ Directory

- **feature_cvpat.R** (~530 LOC) - Core CVPAT implementation
  - `assess_cvpat()` - Compare PLS model against Linear Model (LM) and Indicator Average (IA) benchmarks
  - `assess_cvpat_compare()` - Compare two PLS models' predictive performance
  - Uses k-fold cross-validation with bootstrap significance testing
  - Helper functions: `cvpat_per_construct()`, `lv_loss()`, `overall_loss()`, `bootstrap_cvpat()`

- **feature_coa.R** - Composite Overfit Analysis
  - `coa()` - Main entry point running all COA steps
  - `predictive_deviance()` - Compute predictive deviance via `predict_pls()`
  - `deviance_tree()` - Grow rpart tree to identify deviant case groups
  - `unstable_params()` - Parameter instability analysis by removing deviant groups
  - Depends on `rpart` package for decision tree analysis

- **feature_congruence.R** (~140 LOC) - Congruence coefficient testing
  - `congruence_test()` - Bootstrapped congruence testing for construct validity

- **helpers.R** - Import declarations for roxygen2

### Tests

Tests use testthat 3.0+ with fixtures in `tests/fixtures/`. The helper file `tests/testthat/helper-global.R` defines the `test_folder` path.

### Demo Files

Access demos with: `demo("seminr-pls-cvpat", package = "seminrExtras")`

Available: seminr-help-debugging, seminr-pls-cvpat, seminr-primer-v2-chap2 through chap8

## CI/CD

GitHub Actions workflow (`.github/workflows/rcmdcheck.yml`) runs on macOS and Ubuntu (release + devel). Branches ending with `_noci` skip CI.
