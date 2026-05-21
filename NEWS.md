# seminrExtras 1.0.1

### Fixed

* Tests in `test-cipma-comprehensive.R` and `test-fimix.R` no longer reach
  into seminr's non-exported internals via `seminr:::items_of_construct()`
  and `seminr:::all_endogenous()`. They now use seminrExtras's own local
  helpers of the same name. This avoids breakage against forthcoming
  seminr 2.5.0, which refactors (and renames) those internal helpers.

# seminrExtras 1.0.0

## Major new features

* **Composite Overfit Analysis (COA)**: `assess_coa()`, `predictive_deviance()`,
  `deviance_tree()`, `unstable_params()`, `group_rules()`, `competes()` for
  diagnosing *why* and *for whom* PLS models fail to generalise out-of-sample.

* **Necessary Condition Analysis (NCA)**: `assess_nca()` with fully internal
  CE-FDH and CR-FDH algorithms (no external NCA package dependency).

* **NCA-ESSE**: `assess_nca_esse()` implements the effect size sensitivity
  extension (Becker et al., 2026).

* **Combined Importance-Performance Map Analysis (cIPMA)**: `assess_cipma()`
  integrates IPMA with NCA to classify constructs into actionable priority
  quadrants. `assess_ipma()` provides an IPMA-only convenience wrapper.
  Supports HOC, mediation, and moderation models.

* **FIMIX-PLS**: `assess_fimix()` and `assess_fimix_compare()` for EM-based
  latent class segmentation with multi-start initialisation and information
  criteria comparison.

* **PLS-POS**: `assess_pos()`, `assess_pos_compare()`, and `pos_segments()`
  for prediction-oriented segmentation that maximises the sum of R-squared
  across segments (Becker et al., 2013).

* **CTA-PLS**: `assess_cta()` for confirmatory tetrad analysis with automatic
  indicator borrowing for constructs with fewer than 4 indicators (Gudergan
  et al., 2008).

* **Predictive Contribution of the Mediator (PCM)**: `assess_pcm()` evaluates
  whether a mediating construct improves out-of-sample prediction by comparing
  DA and EA approaches on isolated sub-models (Danks, 2021).

## Improvements

* `assess_cvpat()` and `assess_cvpat_compare()`: fixed bootstrap test
  branches, loss function return types, and reference metadata.

* `congruence_test()`: fixed division guard, upper-triangular masking, and
  bootstrap robustness for `nboot = 0`.

* All features include `print()`, `summary()`, and `plot()` S3 methods.

* Comprehensive test suite (740+ tests).

* Demo files for all features: `demo("seminr-pls-<feature>")`.

* Updated vignette with examples for all features.

## Dependencies

* Requires `seminr` >= 2.4.0.
* `rpart` added to Imports (for COA deviance trees).
* `MASS`, `paran`, `psych`, `learnr` added to Suggests.

# seminrExtras 0.9.0

* Initial CRAN release with CVPAT and congruence testing.
