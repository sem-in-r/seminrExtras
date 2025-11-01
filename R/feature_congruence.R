## Franke, G. R., Sarstedt, M., & Danks, N. P. (2021).
## "Assessing measure congruence in nomological networks". Journal of
## Business Research, 130, 318-334.
#' SEMinR function to bootstrap calculate the congruence coefficient
#'
#' `congruence_test` conducts a single bootstrapped congruence test.
#'
#' @param seminr_model The SEMinR model for CVPAT analysis
#' @param nboot The number of bootstrap subsamples to execute (defaults to 2000).
#' @param seed The seed for reproducibility (defaults to 123).
#' @param alpha The required level of alpha for statistical testing (defaults
#' to 0.05).
#' @param threshold The threshold with which to compare significance testing
#' H0: rc < 1 (defaults to 1).
#'
#' @return A matrix of the estimated congruence coefficient and results of
#' significance testing.
#'
#' @references Franke, G. R., Sarstedt, M., & Danks, N. P. (2021). Assessing
#' measure congruence in nomological networks. Journal of Business Research,
#' 130, 318-334.
#'
#' @examples
#' # Load libraries
#' library(seminr)
#' library(seminrExtras)
#'
#' # Create measurement model ----
#' corp_rep_mm <- constructs(
#'   composite("QUAL", multi_items("qual_", 1:8), weights = mode_B),
#'   composite("PERF", multi_items("perf_", 1:5), weights = mode_B),
#'   composite("CSOR", multi_items("csor_", 1:5), weights = mode_B),
#'   composite("ATTR", multi_items("attr_", 1:3), weights = mode_B),
#'   composite("COMP", multi_items("comp_", 1:3)),
#'   composite("LIKE", multi_items("like_", 1:3))
#' )
#'
#' # Create structural model ----
#' corp_rep_sm <- relationships(
#'   paths(from = c("QUAL", "PERF", "CSOR", "ATTR"), to = c("COMP", "LIKE"))
#' )
#'
#' # Estimate the model ----
#' corp_rep_pls_model <- estimate_pls(
#'   data = corp_rep_data,
#'   measurement_model = corp_rep_mm,
#'   structural_model  = corp_rep_sm,
#'   missing = mean_replacement,
#'  missing_value = "-99")
#'
#' # Assess the base model ----
#' congruence_test(seminr_model = corp_rep_pls_model,
#'                 nboot = 20,
#'                 seed = 123,
#'                 alpha = 0.05,
#'                 threshold = 1)
#'
#' @export
congruence_test <- function(seminr_model,
                            nboot = 2000,
                            seed = 123,
                            alpha = 0.05,
                            threshold = 1) {

  set.seed(seed)
  # Abort if received a higher-order-model or moderated model
  if (!any(class(seminr_model) == "seminr_model")) {
    message("This function only works with SEMinR models. ")
    return()
  }
  # if (!is.null(seminr_model$hoc)) {
  #   message("There is no published solution for applying PLSpredict to higher-order-models")
  #   return()
  # }
  # if (!is.null(seminr_model$interaction)) {
  #   message("There is no published solution for applying PLSpredict to moderated models")
  #   return()
  # }
  construct_names <- colnames(seminr_model$construct_scores)

  calc_congruence <- function(mat,X,Y) {
    return(sum(mat[,X]*mat[,Y])/sqrt(sum(mat[,X]^2)*sum(mat[,Y]^2)))
  }
  # calc_congruence(mat = mat,
  #                 X = "COMP",
  #                 Y = "CUSA")
  combns <- t(utils::combn(construct_names,2))
  ret_array <- array(,
                     dim = list(length(construct_names),length(construct_names),nboot),
                     dimnames = list(construct_names,construct_names,1:nboot))
  for (iter in 1:nboot) {
    it_model <- suppressMessages(seminr::rerun(seminr_model, data = seminr_model$rawdata[sample(nrow(seminr_model$rawdata),nrow(seminr_model$rawdata), replace = TRUE),]))
    ret_mat <- stats::cor(it_model$construct_scores)
    diag(ret_mat) <- seminr::rhoC_AVE(x = it_model)[colnames(ret_mat),1]
    ret_array[,,iter][upper.tri(ret_mat)] <- apply(combns,1,function(x) calc_congruence(ret_mat,x[1],x[2]))
  }

  cor_mat <- stats::cor(seminr_model$construct_scores)
  diag(cor_mat) <- seminr::rhoC_AVE(x = seminr_model)[colnames(ret_mat),1]

  original_matrix <- cor_mat
  original_matrix[lower.tri(original_matrix)] <- 0
  diag(original_matrix) <- 0
  original_matrix[upper.tri(original_matrix)] <- apply(combns,1,function(x) calc_congruence(cor_mat,x[1],x[2]))

  # diag(original_matrix) <- seminr::rhoC_AVE(x = seminr_model)[colnames(seminr_model$construct_scores),1]
  boot_array <- ret_array
  Path <- c()
  original <- c()
  boot_mean <- c()
  boot_SD <- c()
  t_stat <- c()
  lower <- c()
  upper <- c()
  alpha_text <- alpha/2*100
  original_matrix[is.na(original_matrix)] <- 0
  for (i in 1:nrow(original_matrix)) {
    for (j in 1:ncol(original_matrix)) {
      if (original_matrix[i,j]!=0 ) {
        Path <- append(Path, paste(rownames(original_matrix)[i], " -> ", colnames(original_matrix)[j]))
        original <- append(original, original_matrix[i,j])
        boot_mean <- append(boot_mean, (threshold - abs(original_matrix[i,j])))
        boot_SD <- append(boot_SD, stats::sd(boot_array[i,j,]))
        if (original_matrix[i,j]/ stats::sd(boot_array[i,j,]) > 999999999) {
          t_stat <- append(t_stat, NA)
        } else {
          t_stat <- append(t_stat,  (threshold - abs(original_matrix[i,j]))/ stats::sd(boot_array[i,j,]))
        }
        lower <- append(lower, (seminr:::conf_int(boot_array, from = rownames(original_matrix)[i], to = colnames(original_matrix)[j], alpha = alpha))[[1]])
        upper <- append(upper, (seminr:::conf_int(boot_array, from = rownames(original_matrix)[i], to = colnames(original_matrix)[j], alpha = alpha))[[2]])
      }
    }
  }
  return_matrix <- cbind(original, boot_mean, boot_SD, t_stat, lower, upper)
  colnames(return_matrix) <- c( "Original Est.", "Diff", "Bootstrap SD", "T Stat.",paste(alpha_text, "% CI", sep = ""),paste((100-alpha_text), "% CI", sep = ""))
  rownames(return_matrix) <- Path
  return_matrix <- seminr:::convert_to_table_output(return_matrix)
  return(list(results = return_matrix))

}
