# NCA (Necessary Condition Analysis) applied to a PLS-SEM model
library(seminr)
library(seminrExtras)

# Create measurement model ----
mobi_mm <- constructs(
  composite("Image",        multi_items("IMAG", 1:5)),
  composite("Value",        multi_items("PERV", 1:2)),
  composite("Satisfaction", multi_items("CUSA", 1:3)),
  composite("Loyalty",      multi_items("CUSL", 1:3))
)

# Create structural model ----
mobi_sm <- relationships(
  paths(from = c("Image", "Value"), to = "Satisfaction"),
  paths(from = "Satisfaction", to = "Loyalty")
)

# Estimate the model ----
mobi_pls <- estimate_pls(
  data = mobi,
  measurement_model = mobi_mm,
  structural_model  = mobi_sm
)

# Run NCA on Satisfaction ----
nca_sat <- assess_nca(mobi_pls,
                       target = "Satisfaction",
                       test.rep = 1000,
                       seed = 123)

# Print effect sizes and significance
print(nca_sat)

# Full summary with bottleneck tables
summary(nca_sat)

# Effect size bar plot
plot(nca_sat, type = "effects")

# Scatter plots with ceiling lines (from NCA package)
plot(nca_sat, type = "scatter")

# Run NCA on Loyalty ----
nca_loy <- assess_nca(mobi_pls,
                       target = "Loyalty",
                       test.rep = 1000,
                       seed = 123)

print(nca_loy)
summary(nca_loy)
