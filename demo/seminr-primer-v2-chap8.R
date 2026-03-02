### Accompanying Code for:
## Partial Least Squares Structural Equation Modeling (PLS-SEM) Using R - A Workbook (2026)
## Hair, J.F. (Jr), Hult, T.M., Ringle, C.M., Sarstedt, M., Danks, N.P., and Adler, S.

## Chapter 8: Mediation analysis

# Load the SEMinR and seminrExtras libraries
library(seminr)
library(seminrExtras)

# Load the data
corp_rep_data <- corp_rep_data

# Create measurement model ----
corp_rep_mm_ext <- constructs(
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
corp_rep_sm_ext <- relationships(
  paths(from = c("QUAL", "PERF", "CSOR", "ATTR"), to = c("COMP", "LIKE")),
  paths(from = c("COMP", "LIKE"), to = c("CUSA", "CUSL")),
  paths(from = c("CUSA"),         to = c("CUSL"))
)

# Estimate the model ----
corp_rep_pls_model_ext <- estimate_pls(
  data = corp_rep_data,
  measurement_model = corp_rep_mm_ext,
  structural_model  = corp_rep_sm_ext,
  missing = mean_replacement,
  missing_value = "-99")

# Extract the summary
summary_corp_rep_ext <- summary(corp_rep_pls_model_ext)

# Bootstrap the model ----
boot_corp_rep_ext <- bootstrap_model(
  seminr_model = corp_rep_pls_model_ext,
  nboot = 1000,
  cores = parallel::detectCores(),
  seed = 123)

# Summarize the results of the bootstrap
summary_boot_corp_rep_ext <- summary(boot_corp_rep_ext,
                                     alpha = 0.05)

# Inspect total indirect effects
summary_corp_rep_ext$total_indirect_effects

# Inspect indirect effects
specific_effect_significance(boot_corp_rep_ext,
                             from = "COMP",
                             through = "CUSA",
                             to = "CUSL",
                             alpha = 0.1)
specific_effect_significance(boot_corp_rep_ext,
                             from = "LIKE",
                             through = "CUSA",
                             to = "CUSL",
                             alpha = 0.05)

# Inspect the direct effects
summary_corp_rep_ext$paths

# Inspect the confidence intervals for direct effects
summary_boot_corp_rep_ext$bootstrapped_paths

# Calculate the sign of p1*p2*p3 for LIKE->CUSA->CUSL
summary_corp_rep_ext$paths["LIKE", "CUSL"] *
  summary_corp_rep_ext$paths["LIKE","CUSA"] *
  summary_corp_rep_ext$paths["CUSA","CUSL"]

# Calculate the sign of p1*p2*p3 for COMP->CUSA->CUSL
summary_corp_rep_ext$paths["COMP", "CUSL"] *
  summary_corp_rep_ext$paths["COMP","CUSA"] *
  summary_corp_rep_ext$paths["CUSA","CUSL"]

# Calculate the effect size v for LIKE->CUSA->CUSL
summary_corp_rep_ext$paths["LIKE","CUSA"]^2 *
  summary_corp_rep_ext$paths["CUSA","CUSL"]^2

# Calculate the effect size v for COMP->CUSA->CUSL
summary_corp_rep_ext$paths["COMP","CUSA"]^2 *
  summary_corp_rep_ext$paths["CUSA","CUSL"]^2

# Conduct a moderated mediation analysis ----
# Create the measurement model ----
corp_rep_mm_mod <- constructs(
  composite("QUAL", multi_items("qual_", 1:8), weights = mode_B),
  composite("PERF", multi_items("perf_", 1:5), weights = mode_B),
  composite("CSOR", multi_items("csor_", 1:5), weights = mode_B),
  composite("ATTR", multi_items("attr_", 1:3), weights = mode_B),
  composite("COMP", multi_items("comp_", 1:3)),
  composite("LIKE", multi_items("like_", 1:3)),
  composite("CUSA", single_item("cusa")),
  composite("SC", multi_items("switch_", 1:4)),
  composite("CUSL", multi_items("cusl_", 1:3)),
  interaction_term(iv = "CUSA", moderator = "SC", method = two_stage))

# Create the structural model ----
corp_rep_sm_mod <- relationships(
  paths(from = c("QUAL", "PERF", "CSOR", "ATTR"), to = c("COMP", "LIKE")),
  paths(from = c("COMP", "LIKE"), to = c("CUSA", "CUSL")),
  paths(from = c("CUSA", "SC", "CUSA*SC"), to = c("CUSL"))
)

# Estimate the new model with moderator ----
corp_rep_pls_model_mod <- estimate_pls(
  data = corp_rep_data,
  measurement_model = corp_rep_mm_mod,
  structural_model = corp_rep_sm_mod,
  missing = mean_replacement,
  missing_value = "-99"
)

# Extract the summary
sum_corp_rep_mod <- summary(corp_rep_pls_model_mod)

# Bootstrap the model ----
boot_corp_rep_mod <- bootstrap_model(
  seminr_model = corp_rep_pls_model_mod,
  seed = 12345,
  nboot = 1000)

# Index of moderated mediation ----
## Get the point estimate for the index of moderated mediation
# Compute p1*p5 for COMP
sum_corp_rep_mod$paths["COMP","CUSA"] *
  sum_corp_rep_mod$paths["CUSA*SC","CUSL"]

# Compute p1*p5 for LIKE
sum_corp_rep_mod$paths["LIKE","CUSA"] *
  sum_corp_rep_mod$paths["CUSA*SC","CUSL"]

## Get confidence intervals and p-values for the index of moderated mediation
# create a helper function to compute p-values
p_val <- function(x){
  2*min(mean(x<=0), mean(x>=0))
}

# Compute p1*p5 and 95% confidence interval for COMP
p1_p5_comp <- boot_corp_rep_mod$boot_path["COMP", "CUSA", ] *
  boot_corp_rep_mod$boot_path["CUSA*SC", "CUSL", ]
quantile(p1_p5_comp, probs = c(0.025, 0.975))

# Compute the p-value of path p1*p5 for COMP
p_val(p1_p5_comp)

# Compute p1*p5 and 95% confidence interval for LIKE
p1_p5_like <- boot_corp_rep_mod$boot_path["LIKE", "CUSA", ] *
  boot_corp_rep_mod$boot_path["CUSA*SC", "CUSL", ]
quantile(p1_p5_like, probs = c(0.025, 0.975))

# Compute the p-value of path p1*p5 for LIKE
p_val(p1_p5_like)
