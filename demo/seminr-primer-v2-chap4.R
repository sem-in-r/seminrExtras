### Accompanying Code for:
## Partial Least Squares Structural Equation Modeling (PLS-SEM) Using R - A Workbook (2026)
## Hair, J.F. (Jr), Hult, T.M., Ringle, C.M., Sarstedt, M., Danks, N.P., and Adler, S.

## Chapter 4: Evaluation of reflective measurement models
# This analysis requires the psych and paran packages (refer Section 2.6)
# Load the SEMinR and seminrExtras libraries
library(seminr)
library(seminrExtras)

# Load the data ----
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
con.list <- list(COMP = multi_items("comp_", 1:3),
                 LIKE = multi_items("like_", 1:3),
                 CUSL = multi_items("cusl_", 1:3))

## Collect the cleaned data from the seminr_model object
## (The missing values have already been treated/imputed)
dim_data <- corp_rep_pls_model$data

## Set up helper functions (using psych package)
est_dim <- function(construct, data) {
  ret <- psych::iclust(data[,construct],
                       nclusters = 1)$beta
}

## Test for unidimensionality
test_unidim <- function(con.list, sum.model, dim_data) {
  dims_results <- lapply(con.list, function(x) est_dim(x, dim_data))
  dims_results <- unlist(dims_results)
  names(dims_results) <- names(con.list)
  alphas <- sum.model$reliability[names(con.list),"alpha"]
  dims_results <- data.frame(
    "Cronbach Alpha" = alphas,
    "Revelle Beta" = dims_results)
  return(dims_results)
}

## Principal Components Analysis Test
PCA_results <- matrix(c(base::eigen(cor(dim_data[,con.list[["COMP"]]]))$values,
                        base::eigen(cor(dim_data[,con.list[["CUSL"]]]))$values,
                        base::eigen(cor(dim_data[,con.list[["LIKE"]]]))$values),
                        nrow = 3, byrow = TRUE,
                        dimnames = list(c("COMP", "CUSL", "LIKE"),
                                        c("PC1 EV", "PC2 EV", "PC3 EV")))

# Parallel analysis (using paran package)
paran_comp <- paran::paran(dim_data[,con.list[["COMP"]]],
                           iterations=1000,
                           centile=95)
paran_cusl <- paran::paran(dim_data[,con.list[["CUSL"]]],
                           iterations=1000,
                           centile=95)
paran_like <- paran::paran(dim_data[,con.list[["LIKE"]]],
                           iterations=1000,
                           centile=95)
para_results <- matrix(c(paran_comp$AdjEv,paran_cusl$AdjEv,paran_like$AdjEv),
                       nrow = 3, byrow = TRUE,
                       dimnames = list(c("COMP", "CUSL", "LIKE"),
                                      c("PC1 EV", "PC2 EV", "PC3 EV")))

## Display the results of the unidimensionality tests
round(test_unidim(con.list, summary_corp_rep, dim_data),2)
round(PCA_results,2)
round(para_results,2)

# Inspect the outer loadings
summary_corp_rep$loadings

# Inspect the indicator reliability
summary_corp_rep$loadings^2

# Inspect the internal consistency and reliability
summary_corp_rep$reliability

# Plot the reliabilities of constructs
plot(summary_corp_rep$reliability)

# HTMT criterion
summary_corp_rep$validity$htmt

# Bootstrap the model
boot_corp_rep <- bootstrap_model(seminr_model = corp_rep_pls_model,
                                 nboot = 1000)

# Store the summary of the bootstrapped model
sum_boot_corp_rep <- summary(boot_corp_rep, alpha = 0.10)

# Extract the bootstrapped HTMT
sum_boot_corp_rep$bootstrapped_HTMT

# Calculate the congruence coefficient rc
congruence_test(corp_rep_pls_model, alpha = 0.10)$results

