# In this demo we seek to illustrate how users can check for potential bugs in their data and models.
# We consider each of potential mistakes and how they might be checked.
library(seminr)
library(seminrExtras)

## Using the `assess_syntax` argument in `estimate_pls()` ----
## Problem 1. Misspelled construct in MM ----
# intentionally misspell the construct `COMP` as `COP`
error_mm <- constructs(
  composite("COP", multi_items("comp_", 1:3)),
  composite("LIKE", multi_items("like_", 1:3)),
  composite("CUSA", single_item("cusa")),
  composite("CUSL", multi_items("cusl_", 1:3)))

simple_sm <- relationships(
  paths(from = c("COMP", "LIKE"), to = c("CUSA", "CUSL")),
  paths(from = c("CUSA"), to = c("CUSL")))

# Note that I have specified the argument `assess_syntax` as TRUE
# Now estimating the model should output an error
estimate_pls(data = corp_rep,
             measurement_model = error_mm,
             structural_model = simple_sm,
             assess_syntax = TRUE)

## Problem 2. Misspelled construct in SM ----
# correct measurement model
simple_mm <- constructs(
  composite("COMP", multi_items("comp_", 1:3)),
  composite("LIKE", multi_items("like_", 1:3)),
  composite("CUSA", single_item("cusa")),
  composite("CUSL", multi_items("cusl_", 1:3)))

# intentionally misspell the construct `COMP` as `COP`
error_sm <- relationships(
  paths(from = c("COP", "LIKE"), to = c("CUSA", "CUSL")),
  paths(from = c("CUSA"), to = c("CUSL")))

# Now estimating the model should output an error
estimate_pls(data = corp_rep,
             measurement_model = simple_mm,
             structural_model = error_sm,
             assess_syntax = TRUE)

## Problem 3. Missing construct in interaction ----
# Specify the measurement model with interaction
interaction_mm <- constructs(
  composite("COMP", multi_items("comp_", 1:3)),
  composite("LIKE", multi_items("like_", 1:3)),
  interaction_term("COMP", "LIKE", method = two_stage),
  composite("CUSL", multi_items("cusl_", 1:3)))

# Specify the structural model including the interaction,
# but omit the moderator (LIKE)
interaction_sm <- relationships(
  paths(from = c("COMP",  "COMP*LIKE"), to = "CUSL"))

# Now estimating the model should output an error
estimate_pls(data = corp_rep_data,
             measurement_model = interaction_mm,
             structural_model = interaction_sm,
             assess_syntax = TRUE)

## Problem 4. Incorrectly specified indicator names vs data
# intentionally misspell the indicators of COMP as co_1 through co_3
# rather than comp_1 through comp_3.
error_mm <- constructs(
  composite("COP", multi_items("co_", 1:3)),
  composite("LIKE", multi_items("like_", 1:3)),
  composite("CUSA", single_item("cusa")),
  composite("CUSL", multi_items("cusl_", 1:3)))

estimate_pls(data = corp_rep_data,
             measurement_model = error_mm,
             structural_model = error_sm,
             assess_syntax = TRUE)

## But what do I do if I still have a bug?
# Plotting the measurement model and structural model pre-estimation can
# be a powerful tool for debugging your model.

# Create measurement model ----
simple_mm <- constructs(
  composite("COMP", multi_items("comp_", 1:3)),
  composite("LIKE", multi_items("like_", 1:3)),
  composite("CUSA", single_item("cusa")),
  composite("CUSL", multi_items("cusl_", 1:3)))

# Plot the conceptual measurement model
# eyeball the construct and indicator names and ensure they match the data
# and structural model.
plot(simple_mm, theme = seminr_theme_academic())

# Create structural model ----
simple_sm <- relationships(
  paths(from = c("COMP", "LIKE"), to = c("CUSA", "CUSL")),
  paths(from = c("CUSA"), to = c("CUSL")))

# Plot the conceptual structural model
# Eyeball the construct names and ensure they match the data and structural
# model, also check that the specified paths are correctly reflected.
plot(simple_sm, theme = seminr_theme_academic())
# Once you are satisfied that the models are correctly specified,
# you can proceed to estimate the model.
