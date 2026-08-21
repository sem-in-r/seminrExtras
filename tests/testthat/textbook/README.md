# Golden values from the printed book

Numbers transcribed **by hand from the page proofs** of

> Hair, Hult, Ringle, Sarstedt, Danks & Adler (2026).
> *Partial Least Squares Structural Equation Modeling (PLS-SEM) Using R: A Workbook*
> (2nd ed.). Springer. ISBN 978-3-032-28321-4.

These are an **oracle**, not a snapshot. Never regenerate them by running the
code — that would make the tests self-confirming and defeat the entire point.
Change a value here only when a corrected figure has been read off a proof.

| File | Source | Produced by |
|---|---|---|
| `fig_4_11_congruence.csv` | Fig. 4.11, ch. 4 | `congruence_test(corp_rep_pls_model, alpha = 0.10)` |
| `fig_8_7_indirect.csv` | Fig. 8.7, ch. 8 (p. 185) | `specific_effect_significance(boot_corp_rep_ext, ..., alpha = 0.05)` on a `nboot = 1000, seed = 123` bootstrap |
| `fig_6_9_cvpat.csv` | Fig. 6.9, ch. 6 (p. 143) | `assess_cvpat(corp_rep_pls_model_ext, testtype = "greater", nboot = 2000, seed = 123, technique = predict_DA, noFolds = 10, reps = 10)` |

## Fig. 4.11 — the printed proof is WRONG; these are the replacement values

`Fig.4.11_NEW.png` (the corrected screenshot supplied to Springer) is the
source, **not** the figure currently set on p. 100 of
`978-3-032-28321-4_Book_PrintPDF_CMR2_MSA2_SA.pdf`. That printed figure carries
two independent defects, because it was generated under seminrExtras <= 1.0.1:

1. **Transposed pair labels.** It prints `COMP -> CUSL = 0.893` and
   `LIKE -> CUSA = 0.875`; those two values belong to the opposite pairs.
   Coefficients were generated in `combn()` (row-major) order but written via
   `upper.tri()` (column-major) order — fixed in **1.0.2**.
2. **Superseded reliability diagonal.** It used rho_C. The default became
   **rho_A** in **1.0.3**, to agree with SmartPLS.

Printed (wrong): 0.961, 0.841, 0.893, 0.875, 0.940, 0.961
Correct (rho_A): 0.971, 0.848, 0.891, 0.902, 0.954, 0.967

If Springer fails to swap the figure, the book ships mislabelled coefficients
from a superseded default. **Confirm the replacement was applied.**

The surrounding prose survives either way: every upper CI stays below 1, so
"all constructs have a congruence coefficient significantly smaller than 1.0"
holds under both estimators.

## Fig. 6.9 — the prose quotes these numbers

Unusually, the body text on p. 143 quotes the figure directly: -0.873 overall
and -1.100 for CUSL against IA; -0.060 (p < 0.001) overall against LM; -0.034
(p = 0.124) for COMP and 0.015 (p = 0.729) for CUSA. So a drift here falsifies
sentences, not just a screenshot. `boot_t` is signed opposite to `diff`, as
printed.

The `pls_loss` column is identical across the LM and IA blocks by construction
(same model, same folds) — a free internal consistency check on the transcription.

## Fig. 8.7 — corrected values, superseding the printed figure

Susi requested (21 Aug 2026) that the FIRST indirect effect use `alpha = 0.05`
rather than `alpha = 0.1`, so both effects are tested at the level the text
claims. The printed p. 185 figure has `COMP->CUSA->CUSL` at 5%/95% =
[0.015, 0.136]; these golden values carry the corrected 2.5%/97.5% =
[0.003, 0.146]. Everything else in the figure was verified bit-identical.

This matters beyond formatting: the body text asserts the effects are
"significant at the specified 5% level", which was **not true** of the first
effect as printed — its interval was computed at the 10% level. The corrected
interval still excludes zero, but only just (lower bound 0.003 vs 0.015).

The book's `demo/seminr-primer-v2-chap8.R` and the `ALL CODE FILEs` copy were
updated to `alpha = 0.05`; `Quarto.qmd` already had it. As with the Fig. 8.12
seed discrepancy, the printed figures track the DEMO files, not the Quarto.
