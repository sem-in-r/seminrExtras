# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

seminrExtras is a supplementary R package that extends SEMinR (Structural Equation Modeling in R) with advanced analysis tools. It is NOT standalone—it requires seminr >= 2.4.0.

**Primary features:**
- Cross-Validated Predictive Ability Test (CVPAT) - `assess_cvpat()` and `assess_cvpat_compare()`
- Predictive Contribution of the Mediator (PCM) - `assess_pcm()`
- Combined Importance-Performance Map Analysis (cIPMA) - `assess_cipma()` and `assess_ipma()`
- Composite Overfit Analysis (COA) - `assess_coa()`, `predictive_deviance()`, `deviance_tree()`, `unstable_params()`
- Necessary Condition Analysis (NCA) - `assess_nca()` and `assess_nca_esse()`
- FIMIX-PLS (Finite Mixture PLS) - `assess_fimix()`
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

Use only seminr's **exported** functions (e.g. `seminr::estimate_pls()`,
`seminr::predict_pls()`). Do **not** use `seminr:::` to reach into
internals — if a needed helper isn't exported, add a local copy to
`R/helpers.R`.

**Model types that must be supported:**
- Simple path models
- Mediation models
- Moderation models (interaction terms)
- Higher-order construct (HOC) models

### Key Patterns

- Functions accept seminr model objects from `estimate_pls()`
- Avoid reimplementing prediction/CV — use `seminr::predict_pls()` instead
- No new `seminr:::` calls — use exported seminr functions or local
  helpers in `R/helpers.R`

### Model Type Compatibility

All features must be designed from the start to handle different model types. Key structural details:

**HOC (Higher-Order Constructs):**
- `model$hoc` is `TRUE` when HOC constructs are present
- For **both** two-stage and repeated-indicators approaches, `items_of_construct("HOC", model)` returns LOC names (not raw indicators)
- LOC construct scores are in `model$data` but NOT on the measurement scale — never treat them as raw indicators
- LOC constructs typically do NOT appear in `construct_scores` or `path_coef` — only the HOC does
- To detect HOC items: check if items are also columns in `model$outer_weights` (i.e., they are construct names)
- For indicator-level operations (e.g., IPMA performance), chain recursively: HOC → LOC → actual indicators

**Moderation (interaction terms):**
- Interaction constructs have `*` in their name (e.g., `"Image*Value"`)
- Two-stage: single artificial indicator `"Image*Value_intxn"` with weight = 1.0
- Product indicator: product indicator names like `"IMAG1*PERV1"`
- Interaction constructs ARE in `construct_scores` and `path_coef`
- For IPMA/cIPMA: exclude interaction constructs (performance not meaningful on 0–100 scale)
- For NCA/FIMIX: interaction constructs can be included (they have valid construct scores)

**Mediation:**
- No special handling needed — standard constructs with indirect paths
- Total effects via `(I - B)^{-1} - I` properly captures indirect effects

### Implementing Internal Algorithms vs. External Dependencies

Prefer implementing core algorithms internally when:
- The algorithm is well-defined and not excessively complex
- Full control over PLS-specific aspects is needed
- It avoids adding a runtime dependency

Place external packages in `Suggests` (not `Imports`) when they provide optional enhancements only.

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

- **feature_cvpat.R** (~475 LOC) - Core CVPAT implementation
  - `assess_cvpat()` - Compare PLS model against LM and IA benchmarks
  - `assess_cvpat_compare()` - Compare two PLS models' predictive performance

- **feature_pcm.R** (~410 LOC) - Predictive Contribution of the Mediator
  - `assess_pcm()` - Compare DA vs EA predictions to quantify mediator contribution
  - Isolates each mediation path into sub-models for independent evaluation

- **feature_cipma.R** (~780 LOC) - Combined Importance-Performance Map Analysis
  - `assess_cipma()` - cIPMA (IPMA + NCA integration)
  - `assess_ipma()` - IPMA-only convenience wrapper
  - Handles HOC by chaining LOC → indicator performance recursively

- **feature_nca.R** (~1060 LOC) - Necessary Condition Analysis
  - `assess_nca()` - NCA with internal CE-FDH/CR-FDH algorithms
  - `assess_nca_esse()` - NCA effect size sensitivity extension

- **feature_coa.R** (~870 LOC) - Composite Overfit Analysis
  - `assess_coa()`, `predictive_deviance()`, `deviance_tree()`, `unstable_params()`

- **feature_fimix.R** - FIMIX-PLS (Finite Mixture PLS)
  - `assess_fimix()` - EM-based latent class segmentation for K=1..max_k
  - `fimix_segments()` - Extract segment-specific re-estimated models

- **feature_congruence.R** (~235 LOC) - Congruence coefficient testing
  - `congruence_test()` - Bootstrapped congruence testing

- **helpers.R** (~215 LOC) - Shared validation, loss calculation, and bootstrap helpers

### Tests

Tests use testthat 3.0+ with fixtures in `tests/fixtures/`. The helper file `tests/testthat/helper-global.R` defines the `test_folder` path.

**Testing pattern:** Pre-compute expensive results (model estimation, CV, EM) once at the top of each test file and reuse across tests. Keep bootstrap/permutation reps low (0 or 20-50) for fast CI runs.

### Demo Files

Access demos with: `demo("seminr-pls-cvpat", package = "seminrExtras")`

Available: seminr-help-debugging, seminr-pls-cvpat, seminr-primer-v2-chap2 through chap8

## CI/CD

GitHub Actions workflow (`.github/workflows/rcmdcheck.yml`) runs on macOS and Ubuntu (release + devel). Branches ending with `_noci` skip CI.
