# FIMIX-PLS (Finite Mixture PLS) with SEMinR
library(seminr)
library(seminrExtras)

# Create measurement model ----
corp_rep_mm <- constructs(
  composite("QUAL", multi_items("qual_", 1:8), weights = mode_B),
  composite("PERF", multi_items("perf_", 1:5), weights = mode_B),
  composite("CSOR", multi_items("csor_", 1:5), weights = mode_B),
  composite("ATTR", multi_items("attr_", 1:3), weights = mode_B),
  composite("COMP", multi_items("comp_", 1:3)),
  composite("LIKE", multi_items("like_", 1:3)),
  composite("CUSA", single_item("cusa")),
  composite("CUSL", multi_items("cusl_", 1:3))
)

# Create structural model ----
corp_rep_sm <- relationships(
  paths(from = c("QUAL", "PERF", "CSOR", "ATTR"), to = c("COMP", "LIKE")),
  paths(from = c("COMP", "LIKE"), to = c("CUSA", "CUSL")),
  paths(from = c("CUSA"),         to = c("CUSL"))
)

# Estimate the model ----
corp_pls <- estimate_pls(
  data = corp_rep_data,
  measurement_model = corp_rep_mm,
  structural_model  = corp_rep_sm,
  missing = mean_replacement,
  missing_value = "-99")

# ============================================================================
# FIMIX-PLS for a single K (assess_fimix)
# ============================================================================
# FIMIX-PLS uses EM-based latent class segmentation to uncover unobserved
# heterogeneity. Observations are probabilistically assigned to K segments,
# each with segment-specific structural path coefficients.

fimix_k2 <- assess_fimix(corp_pls, K = 2, nstart = 10, seed = 123)

# Print: segment proportions, path coefficients, fit indices
print(fimix_k2)

# Summary: detailed segment-level results
summary(fimix_k2)

# Visualize segment proportions and path coefficients
plot(fimix_k2)

# ============================================================================
# Compare across K values (assess_fimix_compare)
# ============================================================================
# Estimate FIMIX for multiple K values and compare information criteria
# (AIC, BIC, CAIC, etc.) to select the optimal number of segments.

fimix_compare <- assess_fimix_compare(corp_pls,
                                       K_range = 2:4,
                                       nstart = 10,
                                       seed = 123)

# Print: information criteria table across K values
print(fimix_compare)

# Plot: IC comparison across K
plot(fimix_compare)
