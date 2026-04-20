# Combined Importance-Performance Map Analysis (cIPMA) with SEMinR
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
  paths(from = "Image",        to = c("Expectation", "Satisfaction", "Loyalty")),
  paths(from = "Expectation",  to = c("Value", "Satisfaction")),
  paths(from = "Value",        to = "Satisfaction"),
  paths(from = "Satisfaction", to = "Loyalty")
)

# Estimate the model ----
mobi_pls <- estimate_pls(data = mobi,
                          measurement_model = mobi_mm,
                          structural_model  = mobi_sm)

# ============================================================================
# IPMA-only analysis (assess_ipma)
# ============================================================================
# IPMA computes importance (unstandardized total effects) and performance
# (rescaled 0-100 construct scores) for each predictor of a target.

ipma_result <- assess_ipma(mobi_pls,
                            target = "Loyalty",
                            scale_min = 1,
                            scale_max = 10)

print(ipma_result)
summary(ipma_result)

# Standard IPMA map
plot(ipma_result, type = "ipma")

# ============================================================================
# cIPMA: Combined IPMA + NCA (assess_cipma)
# ============================================================================
# cIPMA extends IPMA by integrating NCA to identify which constructs
# are necessary conditions for the target outcome.

cipma_result <- assess_cipma(mobi_pls,
                              target = "Loyalty",
                              scale_min = 1,
                              scale_max = 10,
                              nca_test.rep = 1000,
                              seed = 123)

# Print: importance, performance, necessity classification
print(cipma_result)

# Summary: full IPMA + NCA details
summary(cipma_result)

# cIPMA map (importance vs. performance, with NCA overlay)
plot(cipma_result, type = "cipma")

# Standard IPMA map (without NCA distinction)
plot(cipma_result, type = "ipma")

# Use standardized total effects for importance axis
plot(cipma_result, importance_metric = "standardized")
