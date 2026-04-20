# PLS-POS (Prediction-Oriented Segmentation) with SEMinR
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
# PLS-POS for a single K (assess_pos)
# ============================================================================
# PLS-POS uses a deterministic hill-climbing approach to maximize the sum of
# R-squared across all endogenous constructs across K segments. Unlike
# FIMIX-PLS, PLS-POS makes no distributional assumptions and can detect
# heterogeneity in both structural and (formative) measurement models.

pos_k2 <- assess_pos(corp_pls, K = 2, nstart = 10, max_iter = 100, seed = 123)

# Print: segment sizes, R-squared, path coefficients
print(pos_k2)

# Summary: detailed comparison with global model
summary(pos_k2)

# Visualize segment proportions
plot(pos_k2, type = "segments")

# Visualize R-squared across segments
plot(pos_k2, type = "rsquared")

# Visualize path coefficients across segments
plot(pos_k2, type = "paths")

# ============================================================================
# Extract segment-specific models (pos_segments)
# ============================================================================
# Retrieve the fully re-estimated PLS models for each segment
seg_models <- pos_segments(pos_k2)
summary(seg_models[[1]])

# ============================================================================
# Compare across K values (assess_pos_compare)
# ============================================================================
# Estimate PLS-POS for multiple K values and compare the objective criterion
# (sum of R-squared) to determine the optimal number of segments.

pos_compare <- assess_pos_compare(corp_pls,
                                   K_range = 2:4,
                                   nstart = 10,
                                   max_iter = 100,
                                   seed = 123)

# Print: objective comparison table across K values
print(pos_compare)

# Plot: objective criterion vs K
plot(pos_compare)
