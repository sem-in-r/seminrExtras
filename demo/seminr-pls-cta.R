# Confirmatory Tetrad Analysis (CTA-PLS) with SEMinR
library(seminr)
library(seminrExtras)

# Create measurement model ----
mobi_mm <- constructs(
  composite("Image",        multi_items("IMAG", 1:5)),
  composite("Expectation",  multi_items("CUEX", 1:3)),
  composite("Value",        multi_items("PERV", 1:2)),
  composite("Satisfaction", multi_items("CUSA", 1:3)),
  composite("Loyalty",      multi_items("CUSL", 1:3))
)

# Create structural model ----
mobi_sm <- relationships(
  paths(from = "Image",       to = c("Expectation", "Satisfaction", "Loyalty")),
  paths(from = "Expectation", to = c("Value", "Satisfaction")),
  paths(from = "Value",       to = "Satisfaction"),
  paths(from = "Satisfaction", to = "Loyalty")
)

# Estimate the model ----
mobi_pls <- estimate_pls(data = mobi,
                          measurement_model = mobi_mm,
                          structural_model  = mobi_sm)

# ============================================================================
# CTA-PLS with borrowing (default)
# ============================================================================
# CTA-PLS tests whether a construct's measurement model is consistent with
# a reflective specification. Under the null, all vanishing tetrads = 0.
#
# With borrow = TRUE (default), constructs with 2-3 indicators borrow from
# structurally connected constructs to form testable 4-tuples.

cta_result <- assess_cta(mobi_pls, nboot = 5000, seed = 123)

# Per-construct summary: mode, indicators, tetrads tested, verdict
print(cta_result)

# Detailed results: individual tetrads, t-values, CIs, adjusted p-values
summary(cta_result)

# Visualize adjusted p-values per construct
plot(cta_result)

# ============================================================================
# CTA-PLS without borrowing
# ============================================================================
# Only constructs with >= 4 indicators are tested

cta_no_borrow <- assess_cta(mobi_pls, nboot = 5000, borrow = FALSE, seed = 123)
print(cta_no_borrow)
