# =============================================================================
# Demonstration of PCM (Predictive Contribution of the Mediator)
# Reference: Danks (2021), The DATA BASE for Advances in IS, 52(SI), 24-42.
# =============================================================================
# PCM evaluates the predictive contribution of a mediating construct by
# comparing predictions from the Direct Antecedents (DA) and Earliest
# Antecedents (EA) approaches on isolated mediation sub-models.
# =============================================================================

library(seminr)
library(seminrExtras)

# Load the MOBI customer satisfaction dataset bundled with seminr
data(mobi)

# --- Step 1: Specify and estimate a mediation model ---
# Image -> Satisfaction -> Loyalty (with direct path Image -> Loyalty)
mobi_mm <- constructs(
  composite("Image",        multi_items("IMAG", 1:5)),
  composite("Satisfaction", multi_items("CUSA", 1:3)),
  composite("Loyalty",      multi_items("CUSL", 1:3))
)
mobi_sm <- relationships(
  paths(from = "Image", to = "Satisfaction"),
  paths(from = "Satisfaction", to = "Loyalty"),
  paths(from = "Image", to = "Loyalty")
)
pls_model <- estimate_pls(mobi, mobi_mm, mobi_sm)
summary(pls_model)

# --- Step 2: Compute PCM ---
# PCM automatically detects the mediation path: Image -> Satisfaction -> Loyalty
pcm_result <- assess_pcm(pls_model,
                         target  = "Loyalty",
                         noFolds = 10,
                         reps    = 10)

# Print a concise overview with average PCM and classification
pcm_result

# Print a detailed summary with per-indicator PCM values
summary(pcm_result)

# Plot PCM barplot with threshold lines
plot(pcm_result)

# --- Step 3: Multi-mediator model ---
# Now use a more complex model with multiple mediation paths to Loyalty
mobi_mm_full <- constructs(
  composite("Image",        multi_items("IMAG", 1:5)),
  composite("Expectation",  multi_items("CUEX", 1:3)),
  composite("Value",        multi_items("PERV", 1:2)),
  composite("Satisfaction", multi_items("CUSA", 1:3)),
  composite("Loyalty",      multi_items("CUSL", 1:3))
)
mobi_sm_full <- relationships(
  paths(from = "Image",       to = c("Expectation", "Satisfaction", "Loyalty")),
  paths(from = "Expectation", to = c("Value", "Satisfaction")),
  paths(from = "Value",       to = "Satisfaction"),
  paths(from = "Satisfaction", to = "Loyalty")
)
pls_full <- estimate_pls(mobi, mobi_mm_full, mobi_sm_full)

# PCM automatically identifies all mediation paths:
# Image -> Satisfaction -> Loyalty
# Expectation -> Satisfaction -> Loyalty
# Value -> Satisfaction -> Loyalty
pcm_full <- assess_pcm(pls_full,
                        target  = "Loyalty",
                        noFolds = 10,
                        reps    = 10)

# Print results for all mediation paths
pcm_full

# Detailed indicator-level results
summary(pcm_full)

# Compare mediation paths visually
plot(pcm_full)
