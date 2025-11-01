### Accompanying Code for:
## Partial Least Squares Structural Equation Modeling (PLS-SEM) Using R - A Workbook (2021)
## Hair, J.F. (Jr), Hult, T.M., Ringle, C.M., Sarstedt, M., Danks, N.P., and Adler, S.

## Chapter 4: Evaluation of reflective measurement models

# Load the SEMinR library
# remember to install.packages("psych") if you have not already done so
library(seminr)
library(seminrExtras)
library(psych)

# Load the data ----
# hint: try changing corp_rep_data to corp_rep_data2 and compare results
corp_rep_data <- corp_rep_data

# Create measurement model ----
corp_rep_mm <- constructs(
  composite("COMP", multi_items("comp_", 1:3)),
  composite("LIKE", multi_items("like_", 1:3)),
  composite("CUSA", single_item("cusa")),
  composite("CUSL", multi_items("cusl_", 1:3)))

# Create structural model ----
corp_rep_sm <- relationships(
  paths(from = c("COMP", "LIKE"), to = c("CUSA", "CUSL")),
  paths(from = c("CUSA"), to = c("CUSL")))

# Estimate the model
corp_rep_pls_model <- estimate_pls(
  data = corp_rep_data,
  measurement_model = corp_rep_mm,
  structural_model  = corp_rep_sm,
  missing = mean_replacement,
  missing_value = "-99")

# Summarize the model results
summary_corp_rep <- summary(corp_rep_pls_model)

# Inspect iterations
summary_corp_rep$iterations

# Evaluate the undimensionality of reflective constructs
## Set up the list of reflective constructs
cons.list <- list(COMP = multi_items("comp_", 1:3),
                  LIKE = multi_items("like_", 1:3),
                  CUSL = multi_items("cusl_", 1:3))

## Set up the data for dimensionality assessment
dim_data <- scale(corp_rep_data[,unlist(cons.list)])

## Set up helper functions
est_dim <- function(construct) {
  ret <- iclust(corp_rep_data[,construct],
                nclusters = 1)$beta
}

## Test for unidimensionality
test_undim <- function(cons.list, sum.model) {
  dims_results <- lapply(cons.list, est_dim)
  dims_results <- unlist(dims_results)
  names(dims_results) <- names(cons.list)
  alphas <- sum.model$reliability[names(cons.list),"alpha"]
  dims_results <- data.frame(
    "Cronbach Alpha" = alphas,
    "Revelle Beta" = dims_results)
  return(dims_results)
}

## Principal Components Analysis Test
PCA_results <- matrix(c(stats::prcomp(dim_data[,cons.list[["COMP"]]])$sdev,
                        stats::prcomp(dim_data[,cons.list[["CUSL"]]])$sdev,
                        stats::prcomp(dim_data[,cons.list[["LIKE"]]])$sdev),
                        nrow = 3, byrow = TRUE,
                        dimnames = list(c("COMP", "CUSL", "LIKE"),
                                        c("PC1 SD", "PC2 SD", "PC3 SD")))

## Display the results of the unidimensionality tests
round(test_undim(cons.list, summary_corp_rep),2)
round(PCA_results,2)

# Inspect the outer loadings
summary_corp_rep$loadings

# Inspect the indicator reliability
summary_corp_rep$loadings^2

# Inspect the internal consistency and reliability
summary_corp_rep$reliability

# Plot the reliabilities of constructs
plot(summary_corp_rep$reliability)

# Table of the FL criteria
summary_corp_rep$validity$fl_criteria

# HTMT Ratio
summary_corp_rep$validity$htmt

# Bootstrap the model
boot_corp_rep <- bootstrap_model(seminr_model = corp_rep_pls_model,
                                 nboot = 1000)

# Store the summary of the bootstrapped model
sum_boot_corp_rep <- summary(boot_corp_rep, alpha = 0.10)

# Extract the bootstrapped HTMT
sum_boot_corp_rep$bootstrapped_HTMT

# Calculate the congruence coefficient rc
# must use the following line of code to load the development version of seminrExtras
# devtools::install_github(repo = "https://github.com/sem-in-r/seminrExtras.git", ref = "textbook")
congruence_test(corp_rep_pls_model, alpha = 0.10)$results
