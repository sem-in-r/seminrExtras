#' SEMinR function to compare CVPAT loss of two models
#'
#' `assess_cvpat_compare` conducts a CVPAT significance test of loss between
#' two models.
#'
#' @param established_model The base seminr model for CVPAT comparison.
#' @param alternative_model The alternate seminr model for CVPAT comparison.
#' @param testtype Either "two.sided" (default) or "greater".
#' @param nboot The number of bootstrap subsamples to execute (defaults to 2000).
#' @param seed The seed for reproducibility (defaults to 123).
#' @param technique predict_EA or predict_DA (default).
#' @param noFolds Number of folds for k-fold cross validation.
#' @param reps Number of repetitions for cross validation.
#' @param cores Number of cores for parallelization.
#'
#' @return A matrix of the estimated loss and results of significance testing.
#'
#' @seealso [assess_cvpat()] for single model assessment against benchmarks
#'
#' @references Sharma, P. N., Liengaard, B. D., Hair, J. F., Sarstedt, M., &
#' Ringle, C. M. (2022). Predictive model assessment and selection in
#' composite-based modeling using PLS-SEM: extensions and guidelines for
#' using CVPAT. European journal of marketing, 57(6), 1662-1677.
#'
#' Liengaard, B. D., Sharma, P. N., Hult, G. T. M., Jensen, M. B.,
#' Sarstedt, M., Hair, J. F., & Ringle, C. M. (2021). Prediction: coveted,
#' yet forsaken? Introducing a cross‐validated predictive ability test in
#' partial least squares path modeling. Decision Sciences, 52(2), 362-392.
#'
#' @examples
#' # Load libraries
#' library(seminr)
#' library(seminrExtras)
#'
#' # Create measurement model ----
#' corp_rep_mm_ext <- constructs(
#'  composite("QUAL", multi_items("qual_", 1:8), weights = mode_B),
#'  composite("PERF", multi_items("perf_", 1:5), weights = mode_B),
#'  composite("CSOR", multi_items("csor_", 1:5), weights = mode_B),
#'  composite("ATTR", multi_items("attr_", 1:3), weights = mode_B),
#'  composite("COMP", multi_items("comp_", 1:3)),
#'  composite("LIKE", multi_items("like_", 1:3))
#' )
#'
#' alt_mm <- constructs(
#'  composite("QUAL", multi_items("qual_", 1:8), weights = mode_B),
#'  composite("PERF", multi_items("perf_", 1:5), weights = mode_B),
#'  composite("CSOR", multi_items("csor_", 1:5), weights = mode_B),
#'  composite("COMP", multi_items("comp_", 1:3)),
#'  composite("LIKE", multi_items("like_", 1:3))
#' )
#'
#' # Create structural model ----
#'
#' corp_rep_sm_ext <- relationships(
#'  paths(from = c("QUAL", "PERF", "CSOR", "ATTR"), to = c("COMP", "LIKE"))
#' )
#' alt_sm <- relationships(
#'  paths(from = c("QUAL", "PERF", "CSOR"), to = c("COMP", "LIKE"))
#' )
#'
#' # Estimate the model ----
#' established_model <- estimate_pls(
#'  data = corp_rep_data,
#'  measurement_model = corp_rep_mm_ext,
#'  structural_model  = corp_rep_sm_ext,
#'  missing = mean_replacement,
#'  missing_value = "-99")
#'
#' alternative_model <- estimate_pls(
#'  data = corp_rep_data,
#'  measurement_model = alt_mm,
#'  structural_model  = alt_sm,
#'  missing = mean_replacement,
#'  missing_value = "-99")
#'
#' # Function to compare the Loss of two models
#' assess_cvpat_compare(established_model,
#'                     alternative_model ,
#'                     testtype = "two.sided",
#'                     nboot = 20,
#'                     seed = 123,
#'                     technique = predict_DA,
#'                     noFolds = 5,
#'                     reps = 1,
#'                     cores = 1)
#'
#' @export
assess_cvpat_compare <- function(established_model,
                                 alternative_model,
                                 testtype = "two.sided",
                                 nboot = 2000,
                                 seed = 123,
                                 technique = predict_DA,
                                 noFolds = NULL,
                                 reps = NULL,
                                 cores = NULL) {

  # Validate both models for prediction analysis
  if (has_higher_order(established_model) || has_higher_order(alternative_model)) {
    warning("There is no published solution for applying PLSpredict to higher-order models.",
            call. = FALSE)
    return(NULL)
  }

  set.seed(seed)

  # Retrieve endogenous constructs and items using helpers
  endo_lvs1 <- get_endogenous_constructs(established_model)
  endo_lvs2 <- get_endogenous_constructs(alternative_model)
  endo_mvs1 <- get_endogenous_items(established_model, endo_lvs1)
  endo_mvs2 <- get_endogenous_items(alternative_model, endo_lvs2)

  # Validate that models have identical endogenous constructs and items
  same_mvs <- all(endo_mvs1 %in% endo_mvs2) && all(endo_mvs2 %in% endo_mvs1)
  same_lvs <- all(endo_lvs1 %in% endo_lvs2) && all(endo_lvs2 %in% endo_lvs1)

  if (!(same_mvs && same_lvs)) {
    stop("CVPAT can only be applied to models with identical endogenous constructs and measures",
         call. = FALSE)
  }

  # Calculate PLS predictions for each model
  pls_predict_one <- predict_pls(established_model,
                                 technique = technique,
                                 noFolds = noFolds,
                                 reps = reps,
                                 cores = cores)
  pls_predict_two <- predict_pls(alternative_model,
                                 technique = technique,
                                 noFolds = noFolds,
                                 reps = reps,
                                 cores = cores)

  # Extract and name prediction errors
  pls_error_one <- as.matrix(pls_predict_one$items$PLS_out_of_sample_residuals)
  colnames(pls_error_one) <- endo_mvs1
  pls_error_two <- as.matrix(pls_predict_two$items$PLS_out_of_sample_residuals)
  colnames(pls_error_two) <- endo_mvs2

  pls_error_one <- pls_error_one[, endo_mvs1, drop = FALSE]
  pls_error_two <- pls_error_two[, endo_mvs2, drop = FALSE]

  # Calculate LV losses using helper
  lv_losses_one <- calculate_lv_losses(endo_lvs1, established_model, pls_error_one)
  lv_losses_two <- calculate_lv_losses(endo_lvs2, alternative_model, pls_error_two)

  # Calculate overall loss
  pls_overall_one <- overall_loss(lv_losses_one)
  pls_overall_two <- overall_loss(lv_losses_two)

  # Compare models based on endogenous construct overlap
  if (identical(endo_lvs1, endo_lvs2)) {
    # 100% overlap - direct comparison
    pls_v_pls_overall <- bootstrap_cvpat(pls_overall_one, pls_overall_two,
                                         testtype = testtype, nboot = nboot)
    lv_cvpat <- cvpat_per_construct(loss_one = lv_losses_one, loss_two = lv_losses_two,
                                    testtype = testtype, nboot = nboot)
    mat_one <- cbind(colMeans(lv_losses_one), colMeans(lv_losses_two),
                     colMeans(lv_losses_one) - colMeans(lv_losses_two),
                     lv_cvpat[, -1])
  } else {
    # Partial overlap - compare only overlapping constructs
    overlap <- intersect(endo_lvs1, endo_lvs2)

    if (length(overlap) == 0) {
      return(list(results = "Cannot compare directly"))
    }

    message("Not all endogenous vars co-occur in models 1 and 2. Only comparing overlap.")

    pls_v_pls_overall <- bootstrap_cvpat(pls_overall_one, pls_overall_two,
                                         testtype = testtype, nboot = nboot)
    lv_cvpat <- cvpat_per_construct(loss_one = lv_losses_one[, overlap, drop = FALSE],
                                    loss_two = lv_losses_two[, overlap, drop = FALSE],
                                    testtype = testtype, nboot = nboot)
    mat_one <- cbind(colMeans(lv_losses_one)[overlap], colMeans(lv_losses_two)[overlap],
                     colMeans(lv_losses_one)[overlap] - colMeans(lv_losses_two)[overlap],
                     lv_cvpat[, -1, drop = FALSE])
  }

  # Add overall row
  mat_one <- rbind(mat_one, unlist(c(mean(pls_overall_one),
                                     mean(pls_overall_two),
                                     mean(pls_overall_one) - mean(pls_overall_two),
                                     pls_v_pls_overall)))

  # Format output matrix
  mat_out <- matrix(as.numeric(unlist(mat_one)), nrow = nrow(mat_one))
  rownames(mat_out) <- rownames(mat_one)
  rownames(mat_out)[nrow(mat_one)] <- "Overall"
  mat_out <- mat_out[, c(1, 2, 3, 6, 7)]

  colnames(mat_out) <- c("Base Model Loss", "Alt Model Loss", "Diff",
                         "Boot T value", "Boot P Value")
  comment(mat_out) <- "CVPAT as per Sharma, Liengaard, Hair, Sarstedt, & Ringle, (2023).
  Both models under comparison have identical endogenous constructs with identical measurement models.
  Purely exogenous constructs can differ in regards to their relationships with both nomological
  partners and measurement indicators."
  class(mat_out) <- append(class(mat_out), "table_output")

  return(mat_out)
}

## Sharma, P.N., Liengaard, B.D., Hair, J.F., Sarstedt, M., Ringle, C.M. (2023)
## "Predictive model assessment and selection in composite-based modeling using
## PLS-SEM: extensions and guidelines for using CVPAT", European Journal of
## Marketing, Vol. 57 No. 6, pp. 1662-1677.
## DOI: 10.1108/EJM-08-2020-0636

#' SEMinR function for single model CVPAT assessment
#'
#' `assess_cvpat` conducts a single model CVPAT assessment against item average
#' and linear model prediction benchmarks.
#'
#' @param seminr_model The SEMinR model for CVPAT analysis
#' @param testtype Either "two.sided" (default) or "greater".
#' @param nboot The number of bootstrap subsamples to execute (defaults to 2000).
#' @param seed The seed for reproducibility (defaults to 123).
#' @param technique predict_EA or predict_DA (default).
#' @param noFolds Number of folds for k-fold cross validation.
#' @param reps Number of repetitions for cross validation.
#' @param cores Number of cores for parallelization.
#'
#' @return A list containing two matrices: CVPAT_compare_LM (comparison with
#'   linear model) and CVPAT_compare_IA (comparison with indicator average).
#'
#' @seealso [assess_cvpat_compare()] for comparing two PLS models
#'
#' @references Sharma, P. N., Liengaard, B. D., Hair, J. F., Sarstedt, M., &
#' Ringle, C. M. (2022). Predictive model assessment and selection in
#' composite-based modeling using PLS-SEM: extensions and guidelines for
#' using CVPAT. European journal of marketing, 57(6), 1662-1677.
#'
#' Liengaard, B. D., Sharma, P. N., Hult, G. T. M., Jensen, M. B.,
#' Sarstedt, M., Hair, J. F., & Ringle, C. M. (2021). Prediction: coveted,
#' yet forsaken? Introducing a cross‐validated predictive ability test in
#' partial least squares path modeling. Decision Sciences, 52(2), 362-392.
#'
#' @examples
#' # Load libraries
#' library(seminr)
#' library(seminrExtras)
#'
#' # Create measurement model ----
#' corp_rep_mm_ext <- constructs(
#'   composite("QUAL", multi_items("qual_", 1:8), weights = mode_B),
#'   composite("PERF", multi_items("perf_", 1:5), weights = mode_B),
#'   composite("CSOR", multi_items("csor_", 1:5), weights = mode_B),
#'   composite("ATTR", multi_items("attr_", 1:3), weights = mode_B),
#'   composite("COMP", multi_items("comp_", 1:3)),
#'   composite("LIKE", multi_items("like_", 1:3))
#' )
#'
#' # Create structural model ----
#' corp_rep_sm_ext <- relationships(
#'   paths(from = c("QUAL", "PERF", "CSOR", "ATTR"), to = c("COMP", "LIKE"))
#' )
#'
#' # Estimate the model ----
#' corp_rep_pls_model_ext <- estimate_pls(
#'   data = corp_rep_data,
#'   measurement_model = corp_rep_mm_ext,
#'   structural_model  = corp_rep_sm_ext,
#'   missing = mean_replacement,
#'   missing_value = "-99")
#'
#' # Assess the base model ----
#' assess_cvpat(seminr_model = corp_rep_pls_model_ext,
#'              testtype = "two.sided",
#'              nboot = 20,
#'              seed = 123,
#'              technique = predict_DA,
#'              noFolds = 5,
#'              reps = 1,
#'              cores = 1)
#'
#' @export
assess_cvpat <- function(seminr_model,
                         testtype = "two.sided",
                         nboot = 2000,
                         seed = 123,
                         technique = predict_DA,
                         noFolds = NULL,
                         reps = NULL,
                         cores = NULL) {

  set.seed(seed)

  # Validate model using helper
  if (!validate_for_prediction(seminr_model, "assess_cvpat")) {
    return(NULL)
  }

  # Get endogenous constructs and items using helpers
  endo_lvs <- get_endogenous_constructs(seminr_model)
  endo_mvs <- get_endogenous_items(seminr_model, endo_lvs)

  # Calculate Indicator Average (IA) prediction error
  ia_means <- seminr_model$meanData[endo_mvs, drop = FALSE]
  if (length(endo_mvs) > 1) {
    ia_pred_error <- sweep(seminr_model$data[, endo_mvs, drop = FALSE], 2, ia_means)
  } else {
    ia_pred_error <- seminr_model$data[, endo_mvs, drop = FALSE] - ia_means
  }

  # Calculate PLS and LM predictions
  pls_predict <- suppressWarnings(predict_pls(model = seminr_model,
                                              technique = technique,
                                              noFolds = noFolds,
                                              reps = reps,
                                              cores = cores))

  pls_error <- as.matrix(pls_predict$items$PLS_out_of_sample_residuals)
  colnames(pls_error) <- endo_mvs
  lm_error <- as.matrix(pls_predict$items$lm_out_of_sample_residuals)
  colnames(lm_error) <- endo_mvs

  pls_error <- pls_error[, endo_mvs, drop = FALSE]
  lm_error <- lm_error[, endo_mvs, drop = FALSE]

  # Calculate LV-specific losses using helper
  lv_losses_ia <- calculate_lv_losses(endo_lvs, seminr_model, ia_pred_error)
  lv_losses_lm <- calculate_lv_losses(endo_lvs, seminr_model, lm_error)
  lv_losses_pls <- calculate_lv_losses(endo_lvs, seminr_model, pls_error)

  # Calculate overall loss
  ia_overall <- overall_loss(lv_losses_ia)
  lm_overall <- overall_loss(lv_losses_lm)
  pls_overall <- overall_loss(lv_losses_pls)

  # CVPAT bootstrap tests
  pls_v_ia_overall <- bootstrap_cvpat(pls_overall, ia_overall,
                                      testtype = testtype, nboot = nboot)
  pls_v_lm_overall <- bootstrap_cvpat(pls_overall, lm_overall,
                                      testtype = testtype, nboot = nboot)

  ia_cvpat <- cvpat_per_construct(loss_one = lv_losses_pls, loss_two = lv_losses_ia,
                                  testtype = testtype, nboot = nboot)
  lm_cvpat <- cvpat_per_construct(loss_one = lv_losses_pls, loss_two = lv_losses_lm,
                                  testtype = testtype, nboot = nboot)

  # Build LM comparison matrix
  mat_lm <- cbind(colMeans(lv_losses_pls), colMeans(lv_losses_lm),
                  colMeans(lv_losses_pls) - colMeans(lv_losses_lm),
                  lm_cvpat[, -1, drop = FALSE])
  colnames(mat_lm)[1:3] <- c("PLS Loss", "LM Loss", "Diff")
  mat_lm <- rbind(mat_lm, unlist(c(mean(pls_overall), mean(lm_overall),
                                   mean(pls_overall) - mean(lm_overall),
                                   pls_v_lm_overall)))
  mat_lm <- apply(mat_lm, 2, as.numeric)
  rownames(mat_lm) <- c(endo_lvs, "Overall")
  mat_lm <- mat_lm[, c(1, 2, 3, 6, 7)]

  # Build IA comparison matrix
  mat_ia <- cbind(colMeans(lv_losses_pls), colMeans(lv_losses_ia),
                  colMeans(lv_losses_pls) - colMeans(lv_losses_ia),
                  ia_cvpat[, -1, drop = FALSE])
  colnames(mat_ia)[1:3] <- c("PLS Loss", "IA Loss", "Diff")
  mat_ia <- rbind(mat_ia, unlist(c(mean(pls_overall), mean(ia_overall),
                                   mean(pls_overall) - mean(ia_overall),
                                   pls_v_ia_overall)))
  mat_ia <- apply(mat_ia, 2, as.numeric)
  rownames(mat_ia) <- c(endo_lvs, "Overall")
  mat_ia <- mat_ia[, c(1, 2, 3, 6, 7)]

  # Add metadata
  comment(mat_lm) <- comment(mat_ia) <- "CVPAT as per Sharma et al. (2023)."
  class(mat_lm) <- class(mat_ia) <- append(class(mat_lm), "table_output")

  return(list(CVPAT_compare_LM = mat_lm,
              CVPAT_compare_IA = mat_ia))
}
