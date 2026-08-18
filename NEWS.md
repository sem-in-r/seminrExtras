# seminrExtras 1.0.3

### Changed

* `congruence_test()` now **excludes interaction constructs**, naming any it
  drops in a message. An interaction term's measurement is fixed by the product
  method rather than by theory, so its position in a nomological network is not
  interpretable — and it was previously producing coefficients as low as -0.47
  in a framework built around near-redundancy. It is removed from the construct
  set entirely rather than only from the pair list, because Eq. 2 sums over the
  whole set: coefficients for the remaining pairs change accordingly. This
  brings congruence into line with `assess_cta()`, `assess_pos()`,
  `assess_pcm()` and `assess_cipma()`, which already excluded them.

* `congruence_test()` now **refuses higher-order models** with a warning rather
  than returning a number. Two-stage estimation replaces the lower-order
  constructs with a single higher-order composite, and it has not been
  established what belongs on that composite's diagonal, nor whether a
  congruence coefficient between a higher-order and a first-order construct is
  interpretable. Previously such models ran and returned results.

* `congruence_test()` gains a `reliability` argument selecting which estimate
  sits on the diagonal of the construct-correlation matrix: `"rhoA"`, `"rhoC"`,
  `"cronbach"` (Cronbach's alpha) or `"one"`. `"cronbach"` is included chiefly for
  comparison with covariance-based SEM, where rho_A is unavailable. Franke, Sarstedt & Danks (2021, Eq. 2) specify "the reliabilities"
  without fixing an estimator, so all four are in specification.

* **The default is now `"rhoA"`, changed from the previous `"rhoC"`.** This
  aligns seminrExtras with SmartPLS: on the simple corporate reputation model
  the `"rhoA"` diagonal reproduces SmartPLS's published coefficients to three
  decimal places on all six construct pairs, and the estimator choice accounts
  for the whole of the previous disagreement between the two programs. Pass
  `reliability = "rhoC"` to reproduce results generated with 1.0.2 or earlier.

  The estimators differ only where internal consistency is undefined. `"rhoA"`,
  `"rhoC"` and `"cronbach"` all return 1 for single-item constructs, so every
  option coincides in a model built entirely from single indicators. They part
  company on Mode B constructs: `"rhoA"` returns 1, while `"rhoC"` and
  `"cronbach"` compute a value from the indicators. Note this only affects pairs
  that involve a Mode B construct — each matrix column carries exactly one
  reliability, its own construct's, so a coefficient between two reflective
  constructs is invariant to every other construct's diagonal. In a model mixing Mode A and Mode B the `"rhoA"`
  diagonal is therefore set by measurement mode rather than by the data, which
  shifts coefficients systematically — on the extended corporate reputation
  model, reflective-reflective pairs rise while formative-formative pairs fall,
  and the ranking of most-congruent pairs changes. `"one"` removes the
  distinction by treating every construct alike.

### Internal

* Reduced test-suite runtime substantially (CVPAT tests now use fewer
  cross-validation folds and bootstrap subsamples; the PLS-POS comparison is
  pre-computed once and reused) so the package builds well within CRAN's
  check-time limits. No change to user-facing behaviour or coverage.

# seminrExtras 1.0.2

### Fixed

* `congruence_test()` now attaches each congruence coefficient — and its
  bootstrap SD, t-statistic, and confidence interval — to the correct
  construct pair. Coefficients were generated in `combn()` (row-major) order
  but written into the results matrix via `upper.tri()` (column-major) order;
  the two orderings diverge for models with four or more constructs, which
  mislabelled non-adjacent pairs (e.g. `COMP <> CUSL` and `LIKE <> CUSA` were
  swapped in the simple corporate-reputation model, and 10 of 15 pairs were
  mislabelled in the six-construct example). The coefficient values themselves
  were correct; only their pair labels were wrong. Added a regression test
  asserting pair-label correctness against an independent name-indexed
  reference.

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
