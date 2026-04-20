# Congruence Coefficient Testing with SEMinR
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
# Congruence test (congruence_test)
# ============================================================================
# Tests whether the congruence coefficients between PLS composite weights
# and bootstrapped weights are significantly close to 1 (perfect congruence).

cong_result <- congruence_test(mobi_pls,
                                nboot = 2000,
                                seed = 123,
                                alpha = 0.05,
                                threshold = 1)

# Print results
print(cong_result)

# Summary with detailed bootstrap statistics
summary(cong_result)
