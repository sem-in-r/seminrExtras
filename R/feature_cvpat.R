# =============================================================================
# feature_cvpat.R - Cross-Validated Predictive Ability Test (CVPAT)
# =============================================================================
# This file implements CVPAT for PLS-SEM models, as described in:
#
# - Liengaard et al. (2021). "Prediction: coveted, yet forsaken? Introducing a
#   cross-validated predictive ability test in partial least squares path
#   modeling." Decision Sciences, 52(2), 362-392.
#
# - Sharma et al. (2022). "Predictive model assessment and selection in
#   composite-based modeling using PLS-SEM: extensions and guidelines for
#   using CVPAT." European Journal of Marketing, 57(6), 1662-1677.
#
# CVPAT compares the predictive loss of PLS models against benchmarks (LM, IA)
# or against each other using bootstrap-based significance testing.
# =============================================================================

#' Compare CVPAT loss between two PLS models
#'
#' `assess_cvpat_compare` conducts a CVPAT significance test of loss between
#' two models. Use this to determine which of two competing models has
#' superior predictive performance.
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

  # ---------------------------------------------------------------------------
  # Step 1: Validate inputs
  # ---------------------------------------------------------------------------
  # Higher-order constructs cannot be used with PLSpredict (no published method)
  if (has_higher_order(established_model) || has_higher_order(alternative_model)) {
    warning("There is no published solution for applying PLSpredict to higher-order models.",
            call. = FALSE)
    return(NULL)
  }

  # Set seed for reproducibility of cross-validation and bootstrap
  set.seed(seed)

  # ---------------------------------------------------------------------------
  # Step 2: Extract endogenous constructs and their measurement items
  # ---------------------------------------------------------------------------
  # Endogenous constructs are those we predict (have incoming paths)
  endo_lvs1 <- get_endogenous_constructs(established_model)
  endo_lvs2 <- get_endogenous_constructs(alternative_model)
  endo_mvs1 <- get_endogenous_items(established_model, endo_lvs1)
  endo_mvs2 <- get_endogenous_items(alternative_model, endo_lvs2)

  # ---------------------------------------------------------------------------
  # Step 3: Validate model compatibility
  # ---------------------------------------------------------------------------
  # CVPAT comparison requires identical endogenous constructs and items
  # (only the structural paths to exogenous constructs can differ)
  same_mvs <- all(endo_mvs1 %in% endo_mvs2) && all(endo_mvs2 %in% endo_mvs1)
  same_lvs <- all(endo_lvs1 %in% endo_lvs2) && all(endo_lvs2 %in% endo_lvs1)

  if (!(same_mvs && same_lvs)) {
    stop("CVPAT can only be applied to models with identical endogenous constructs and measures",
         call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Step 4: Generate out-of-sample predictions via k-fold cross-validation
  # ---------------------------------------------------------------------------
  # predict_pls uses k-fold CV to generate unbiased prediction errors
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

  # ---------------------------------------------------------------------------
  # Step 5: Extract prediction errors (residuals)
  # ---------------------------------------------------------------------------
  # Residuals = actual - predicted; used to calculate loss
  pls_error_one <- as.matrix(pls_predict_one$items$PLS_out_of_sample_residuals)
  colnames(pls_error_one) <- endo_mvs1
  pls_error_two <- as.matrix(pls_predict_two$items$PLS_out_of_sample_residuals)
  colnames(pls_error_two) <- endo_mvs2

  # Ensure correct column ordering
  pls_error_one <- pls_error_one[, endo_mvs1, drop = FALSE]
  pls_error_two <- pls_error_two[, endo_mvs2, drop = FALSE]

  # ---------------------------------------------------------------------------
  # Step 6: Calculate construct-level losses (mean squared error)
  # ---------------------------------------------------------------------------
  # Loss aggregates item-level squared errors to construct level
  lv_losses_one <- calculate_lv_losses(endo_lvs1, established_model, pls_error_one)
  lv_losses_two <- calculate_lv_losses(endo_lvs2, alternative_model, pls_error_two)

  # ---------------------------------------------------------------------------
  # Step 7: Calculate overall model loss
  # ---------------------------------------------------------------------------
  # Overall loss = mean of construct losses (equal weighting)
  pls_overall_one <- overall_loss(lv_losses_one)
  pls_overall_two <- overall_loss(lv_losses_two)

  # ---------------------------------------------------------------------------
  # Step 8: Perform CVPAT comparison based on construct overlap
  # ---------------------------------------------------------------------------
  if (identical(endo_lvs1, endo_lvs2)) {
    # Case A: 100% overlap - direct comparison of all constructs
    pls_v_pls_overall <- bootstrap_cvpat(pls_overall_one, pls_overall_two,
                                         testtype = testtype, nboot = nboot)
    lv_cvpat <- cvpat_per_construct(loss_one = lv_losses_one, loss_two = lv_losses_two,
                                    testtype = testtype, nboot = nboot)
    # Build results matrix: columns = [Loss1, Loss2, Diff, bootstrap stats]
    mat_one <- cbind(colMeans(lv_losses_one), colMeans(lv_losses_two),
                     colMeans(lv_losses_one) - colMeans(lv_losses_two),
                     lv_cvpat[, -1])
  } else {
    # Case B: Partial overlap - compare only shared constructs
    overlap <- intersect(endo_lvs1, endo_lvs2)

    if (length(overlap) == 0) {
      # No shared constructs - cannot compare
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

  # ---------------------------------------------------------------------------
  # Step 9: Add overall row and format output
  # ---------------------------------------------------------------------------
  mat_one <- rbind(mat_one, unlist(c(mean(pls_overall_one),
                                     mean(pls_overall_two),
                                     mean(pls_overall_one) - mean(pls_overall_two),
                                     pls_v_pls_overall)))

  # Convert to numeric matrix with proper row/column names
  mat_out <- matrix(as.numeric(unlist(mat_one)), nrow = nrow(mat_one))
  rownames(mat_out) <- rownames(mat_one)
  rownames(mat_out)[nrow(mat_one)] <- "Overall"

  # Select final columns: Loss1, Loss2, Diff, Boot T, Boot P
  # (columns 4-5 are Std T and Std P which we exclude)
  mat_out <- mat_out[, c(1, 2, 3, 6, 7)]

  colnames(mat_out) <- c("Base Model Loss", "Alt Model Loss", "Diff",
                         "Boot T value", "Boot P Value")

  # Add metadata for print method
  comment(mat_out) <- "CVPAT as per Sharma, Liengaard, Hair, Sarstedt, & Ringle, (2023).
  Both models under comparison have identical endogenous constructs with identical measurement models.
  Purely exogenous constructs can differ in regards to their relationships with both nomological
  partners and measurement indicators."
  class(mat_out) <- append(class(mat_out), "table_output")

  return(mat_out)
}

# =============================================================================
# Single Model CVPAT Assessment
# =============================================================================

#' Assess single model CVPAT against benchmarks
#'
#' `assess_cvpat` conducts a single model CVPAT assessment against item average
#' and linear model prediction benchmarks. Use this to determine whether your
#' PLS model has meaningful predictive power.
#'
#' Two benchmarks are used:
#' - **Linear Model (LM)**: Multiple regression predicting each item
#' - **Indicator Average (IA)**: Naive benchmark using training set means
#'
#' If PLS loss < LM loss (significantly), the structural model adds predictive value.
#' If PLS loss < IA loss (significantly), the model has basic predictive relevance.
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

  # Set seed for reproducibility
  set.seed(seed)

  # ---------------------------------------------------------------------------
  # Step 1: Validate the model
  # ---------------------------------------------------------------------------
  if (!validate_for_prediction(seminr_model, "assess_cvpat")) {
    return(NULL)
  }

  # ---------------------------------------------------------------------------
  # Step 2: Extract endogenous constructs and items
  # ---------------------------------------------------------------------------
  endo_lvs <- get_endogenous_constructs(seminr_model)
  endo_mvs <- get_endogenous_items(seminr_model, endo_lvs)

  # ---------------------------------------------------------------------------
  # Step 3: Calculate Indicator Average (IA) benchmark errors
  # ---------------------------------------------------------------------------
  # IA benchmark: predict each item using its training set mean
  # This is the simplest possible prediction - no model structure used
  ia_means <- seminr_model$meanData[endo_mvs, drop = FALSE]

  # IA error = actual - mean (for each item)
  if (length(endo_mvs) > 1) {
    # sweep subtracts ia_means from each row across columns
    ia_pred_error <- sweep(seminr_model$data[, endo_mvs, drop = FALSE], 2, ia_means)
  } else {
    # Single item case
    ia_pred_error <- seminr_model$data[, endo_mvs, drop = FALSE] - ia_means
  }

  # ---------------------------------------------------------------------------
  # Step 4: Generate PLS and LM predictions via cross-validation
  # ---------------------------------------------------------------------------
  # predict_pls returns both PLS predictions and LM (linear model) predictions
  # LM benchmark: multiple regression predicting each item from other items
  pls_predict <- suppressWarnings(predict_pls(model = seminr_model,
                                              technique = technique,
                                              noFolds = noFolds,
                                              reps = reps,
                                              cores = cores))

  # Extract prediction errors (residuals)
  pls_error <- as.matrix(pls_predict$items$PLS_out_of_sample_residuals)
  colnames(pls_error) <- endo_mvs
  lm_error <- as.matrix(pls_predict$items$lm_out_of_sample_residuals)
  colnames(lm_error) <- endo_mvs

  # Ensure correct column ordering
  pls_error <- pls_error[, endo_mvs, drop = FALSE]
  lm_error <- lm_error[, endo_mvs, drop = FALSE]

  # ---------------------------------------------------------------------------
  # Step 5: Calculate construct-level losses for all three methods
  # ---------------------------------------------------------------------------
  lv_losses_ia <- calculate_lv_losses(endo_lvs, seminr_model, ia_pred_error)
  lv_losses_lm <- calculate_lv_losses(endo_lvs, seminr_model, lm_error)
  lv_losses_pls <- calculate_lv_losses(endo_lvs, seminr_model, pls_error)

  # ---------------------------------------------------------------------------
  # Step 6: Calculate overall model losses
  # ---------------------------------------------------------------------------
  ia_overall <- overall_loss(lv_losses_ia)
  lm_overall <- overall_loss(lv_losses_lm)
  pls_overall <- overall_loss(lv_losses_pls)

  # ---------------------------------------------------------------------------
  # Step 7: Bootstrap significance tests (PLS vs each benchmark)
  # ---------------------------------------------------------------------------
  # Test 1: PLS vs IA (does model beat naive mean prediction?)
  pls_v_ia_overall <- bootstrap_cvpat(pls_overall, ia_overall,
                                      testtype = testtype, nboot = nboot)
  # Test 2: PLS vs LM (does structural model add value over regression?)
  pls_v_lm_overall <- bootstrap_cvpat(pls_overall, lm_overall,
                                      testtype = testtype, nboot = nboot)

  # Per-construct bootstrap tests
  ia_cvpat <- cvpat_per_construct(loss_one = lv_losses_pls, loss_two = lv_losses_ia,
                                  testtype = testtype, nboot = nboot)
  lm_cvpat <- cvpat_per_construct(loss_one = lv_losses_pls, loss_two = lv_losses_lm,
                                  testtype = testtype, nboot = nboot)

  # ---------------------------------------------------------------------------
  # Step 8: Build output matrix for LM comparison
  # ---------------------------------------------------------------------------
  mat_lm <- cbind(colMeans(lv_losses_pls), colMeans(lv_losses_lm),
                  colMeans(lv_losses_pls) - colMeans(lv_losses_lm),
                  lm_cvpat[, -1, drop = FALSE])
  colnames(mat_lm)[1:3] <- c("PLS Loss", "LM Loss", "Diff")

  # Add overall row
  mat_lm <- rbind(mat_lm, unlist(c(mean(pls_overall), mean(lm_overall),
                                   mean(pls_overall) - mean(lm_overall),
                                   pls_v_lm_overall)))
  mat_lm <- apply(mat_lm, 2, as.numeric)
  rownames(mat_lm) <- c(endo_lvs, "Overall")

  # Select columns: Loss1, Loss2, Diff, Boot T, Boot P
  mat_lm <- mat_lm[, c(1, 2, 3, 6, 7)]

  # ---------------------------------------------------------------------------
  # Step 9: Build output matrix for IA comparison
  # ---------------------------------------------------------------------------
  mat_ia <- cbind(colMeans(lv_losses_pls), colMeans(lv_losses_ia),
                  colMeans(lv_losses_pls) - colMeans(lv_losses_ia),
                  ia_cvpat[, -1, drop = FALSE])
  colnames(mat_ia)[1:3] <- c("PLS Loss", "IA Loss", "Diff")

  # Add overall row
  mat_ia <- rbind(mat_ia, unlist(c(mean(pls_overall), mean(ia_overall),
                                   mean(pls_overall) - mean(ia_overall),
                                   pls_v_ia_overall)))
  mat_ia <- apply(mat_ia, 2, as.numeric)
  rownames(mat_ia) <- c(endo_lvs, "Overall")

  # Select columns: Loss1, Loss2, Diff, Boot T, Boot P
  mat_ia <- mat_ia[, c(1, 2, 3, 6, 7)]

  # ---------------------------------------------------------------------------
  # Step 10: Add metadata and return
  # ---------------------------------------------------------------------------
  comment(mat_lm) <- comment(mat_ia) <- "CVPAT as per Sharma et al. (2023)."
  class(mat_lm) <- class(mat_ia) <- append(class(mat_lm), "table_output")

  return(list(CVPAT_compare_LM = mat_lm,
              CVPAT_compare_IA = mat_ia))
}
