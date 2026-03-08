# Composite Overfit Analysis (COA) with SEMinR
library(seminr)
library(seminrExtras)

# Create measurement model ----
corp_mm <- constructs(
  composite("COMP", multi_items("comp_", 1:3)),
  composite("LIKE", multi_items("like_", 1:3)),
  composite("CUSA", single_item("cusa")),
  composite("CUSL", multi_items("cusl_", 1:3))
)

# Create structural model ----
corp_sm <- relationships(
  paths(from = c("COMP", "LIKE"), to = "CUSA"),
  paths(from = "CUSA", to = "CUSL")
)

# Estimate the model ----
corp_model <- estimate_pls(
  data = corp_rep_data,
  measurement_model = corp_mm,
  structural_model  = corp_sm,
  missing = mean_replacement,
  missing_value = "-99")

# ============================================================================
# Full COA pipeline (assess_coa)
# ============================================================================
# COA detects observation-level overfitting by computing predictive deviance,
# growing a decision tree to identify deviant subgroups, and analyzing
# parameter instability.

coa_result <- assess_coa(corp_model,
                          focal_construct = "CUSL",
                          noFolds = 10, reps = 1, cores = 1,
                          seed = 123)

# Print summary
print(coa_result)
summary(coa_result)

# Predictive deviance distribution
plot(coa_result, type = "pd")

# Decision tree
plot(coa_result, type = "tree")

# Deviant group highlights
plot(coa_result, type = "groups")

# ============================================================================
# Step-by-step COA (individual functions)
# ============================================================================

# Step 1: Compute predictive deviance
pd <- predictive_deviance(corp_model,
                           focal_construct = "CUSL",
                           noFolds = 10, reps = 1, cores = 1,
                           seed = 123)

# Step 2: Identify deviant groups via decision tree
tree <- deviance_tree(pd, bounds = c(0.025, 0.975))

# Step 3: Analyze parameter instability
instab <- unstable_params(corp_model, tree, params = "path_coef")
print(instab)

# Extract decision rules for deviant groups
group_rules(tree)

# Show competing splits at tree nodes
competes(tree)
