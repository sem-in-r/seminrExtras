#' SEMinR function to compare CV-PAT loss of two models
#'
#' `assess_cvpat_compare` conducts a CV-PAT significance test of loss between
#' two models.
#'
#' @param base_model The base model for CV-PAT comparison.
#' @param alt_sm The alternate structural model for CV-PAT comparison.
#' @param testtype Either "two.sided" (default) or "greater".
#' @param BootSamp The number of bootstrap subsamples to execute (defaults to 2000).
#' @param seed The seed for reproducibility (defaults to 123).
#' @param technique predict_EA or predict_DA (default).
#' @param noFolds Mumber of folds for k-fold cross validation.
#' @param reps Number of repetitions for cross validation.
#' @param cores Number of cores for parallelization.
#'
#' @return A matrix of the estimated loss and results of significance testing.
#'
#' @examples
#' # Load libraries
#'library(seminr)
#'
#'# Create measurement model ----
#'corp_rep_mm_ext <- constructs(
#'  composite("QUAL", multi_items("qual_", 1:8), weights = mode_B),
#'  composite("PERF", multi_items("perf_", 1:5), weights = mode_B),
#'  composite("CSOR", multi_items("csor_", 1:5), weights = mode_B),
#'  composite("ATTR", multi_items("attr_", 1:3), weights = mode_B),
#'  composite("COMP", multi_items("comp_", 1:3)),
#'  composite("LIKE", multi_items("like_", 1:3))
#')
#'
#'# Create structural model ----
#'
#'corp_rep_sm_ext <- relationships(
#'  paths(from = c("QUAL", "PERF", "CSOR", "ATTR"), to = c("COMP"))
#')
#'alt_sm <- relationships(
#'  paths(from = c("QUAL", "PERF", "CSOR", "ATTR"), to = c("COMP", "LIKE"))
#')
#'
#'# Estimate the model ----
#'corp_rep_pls_model_ext <- estimate_pls(
#'  data = corp_rep_data,
#'  measurement_model = corp_rep_mm_ext,
#'  structural_model  = corp_rep_sm_ext,
#'  missing = mean_replacement,
#'  missing_value = "-99")
#'
#'# Function to compare the Loss of two models
#'assess_cvpat_compare(base_model = corp_rep_pls_model_ext,
#'                     alt_sm = alt_sm,
#'                     testtype = "two.sided",
#'                     BootSamp = 100,
#'                     technique = predict_DA,
#'                     seed = 123,
#'                     cores = 1)
#'
#' @export
assess_cvpat_compare <- function(base_model,
                                 alt_sm ,
                                 testtype = "two.sided",
                                 BootSamp = 2000,
                                 seed = 123,
                                 technique = predict_DA,
                                 noFolds = NULL,
                                 reps = NULL,
                                 cores = NULL) {
  # Abort if received a higher-order-model or moderated model
  if (!is.null(base_model$hoc)) {
    message("There is no published solution for applying PLSpredict to higher-order-models")
    return()
  }
  if (!is.null(base_model$interaction)) {
    message("There is no published solution for applying PLSpredict to moderated models")
    return()
  }
  set.seed(seed)
  # Estimate model two
  pls_two <- estimate_pls(
    data = base_model$data,
    measurement_model = base_model$measurement_model,
    structural_model  = alt_sm,
    missing = mean_replacement,
    missing_value = base_model$settings$missing_value)


  endo_lvs1 <- seminr:::all_endogenous(base_model$smMatrix)
  endo_lvs2 <- seminr:::all_endogenous(alt_sm)

  endo_mvs1 <- unlist(lapply(endo_lvs1,
                             function(x) seminr:::items_of_construct(construct = x,
                                                                     model = base_model)))
  endo_mvs2 <- unlist(lapply(endo_lvs2,
                             function(x) seminr:::items_of_construct(construct = x,
                                                                     model = pls_two)))

  # Calculate PLS predictions for each model
  pls_predict_model_one <- predict_pls(base_model,
                                       technique = technique,
                                       noFolds = noFolds,
                                       reps = reps,
                                       cores = cores)
  pls_predict_model_two <- predict_pls(pls_two,
                                       technique = technique,
                                       noFolds = noFolds,
                                       reps = reps,
                                       cores = cores)

  pls_predict_error_one_item <- as.matrix(pls_predict_model_one$PLS_out_of_sample_residuals)
  colnames(pls_predict_error_one_item) <- endo_mvs1
  pls_predict_error_two_item <- as.matrix(pls_predict_model_two$PLS_out_of_sample_residuals)
  colnames(pls_predict_error_two_item) <- endo_mvs2

  PLS_predict_error_one <- pls_predict_error_one_item[,endo_mvs1,drop = F]
  PLS_predict_error_two <- pls_predict_error_two_item[,endo_mvs2,drop = F]

  ## Calculate LV losses for each PLS model
  ## model one
  LV_losses_PLS_one <- do.call("cbind", lapply(endo_lvs1,
                                               function(x) lv_loss(construct = x,
                                                                   model = base_model,
                                                                   error = PLS_predict_error_one)))
  ## model two
  LV_losses_PLS_two <- do.call("cbind", lapply(endo_lvs2,
                                               function(x) lv_loss(construct = x,
                                                                   model = pls_two,
                                                                   error = PLS_predict_error_two)))

  # Name LVs
  colnames(LV_losses_PLS_one) <- endo_lvs1
  colnames(LV_losses_PLS_two) <- endo_lvs2

  # Calculate overall loss
  # for PLS model one (base)
  PLS_overall_one <- overall_loss(LV_losses_PLS_one)

  # for PLS model two (alt)
  PLS_overall_two <- overall_loss(LV_losses_PLS_two)

  # If there is 100% overlap in endogenous, then we direct compare
  if (identical(endo_lvs1, endo_lvs2)) {
    # CVPAT: PLS1 vs PLS2 overall
    PLS_v_PLS_overall <- bootstrap_cvpat(PLS_overall_one,
                                        PLS_overall_two,
                                        testtype = testtype,
                                        BootSamp = BootSamp)
    LV_cvpat <- cvpat_per_construct(loss_one = LV_losses_PLS_one,
                                    loss_two = LV_losses_PLS_two,
                                    testtype = testtype,
                                    BootSamp = BootSamp)
    mat_one <- cbind(colMeans(LV_losses_PLS_one),colMeans(LV_losses_PLS_two),
                     colMeans(LV_losses_PLS_one) - colMeans(LV_losses_PLS_two),
                     LV_cvpat[,-1])
  }
  # if there is less than 100% overlap in endogneous, we compare only the
  # relevant endogenous
  if (!identical(endo_lvs1, endo_lvs2)) {
    # CVPAT: PLS1 vs PLS2 overall
    overlap <- intersect(endo_lvs1, endo_lvs2)
    PLS_v_PLS_overall <- bootstrap_cvpat(PLS_overall_one,
                                         PLS_overall_two,
                                         testtype = testtype,
                                         BootSamp = BootSamp)

    LV_cvpat <- cvpat_per_construct(loss_one = LV_losses_PLS_one[,overlap,drop = F],
                                    loss_two = LV_losses_PLS_two[,overlap, drop = F],
                                    testtype = testtype,
                                    BootSamp = BootSamp)
    message("Not all endogenous vars co-occur in models 1 and 2. Only comparing overlap. ")
    mat_one <- cbind(colMeans(LV_losses_PLS_one)[overlap],colMeans(LV_losses_PLS_two)[overlap],
                     colMeans(LV_losses_PLS_one)[overlap] - colMeans(LV_losses_PLS_two)[overlap],
                     LV_cvpat[,-1,drop = F])

  }
  if (length(intersect(endo_lvs1, endo_lvs2) ) == 0) {
    return(list(results = "Cannot compare directly"))
  }

  # mat_one <- cbind(colMeans(LV_losses_PLS_one),colMeans(LV_losses_PLS_two),
  #                  colMeans(LV_losses_PLS_one) - colMeans(LV_losses_PLS_two),
  #                  LV_cvpat[,-1])
  colnames(mat_one)[1:3] <- c("Base Model Loss", "Alt Model Loss", "Diff")
  mat_one <- rbind(mat_one, c(mean(PLS_overall_one),
                              mean(PLS_overall_two),
                              mean(PLS_overall_one) -  mean(PLS_overall_two),
                              PLS_v_PLS_overall))
  rownames(mat_one)[nrow(mat_one)] <- "Overall"
  mat_one <- mat_one[,c(1,2,3,6,7)]
  return(mat_one)
}

## Sharma, P.N., Liengaard, B.D., Hair, J.F., Sarstedt, M., Ringle, C.M. (2023)
## "Predictive model assessment and selection in composite-based modeling using
## PLS-SEM: extensions and guidelines for using CVPAT", European Journal of
## Marketing, Vol. 57 No. 6, pp. 1662-1677.
## DOI: 10.1108/EJM-08-2020-0636
#' SEMinR function to compare CV-PAT loss of two models
#'
#' `assess_cvpat` conducts a single model CV-PAT assessment against item average
#' and linear model benchmarks.
#'
#' @param model The SEMinR model for CV-PAT comparison.
#' @param testtype Either "two.sided" (default) or "greater".
#' @param BootSamp The number of bootstrap subsamples to execute (defaults to 2000).
#' @param seed The seed for reproducibility (defaults to 123).
#' @param technique predict_EA or predict_DA (default).
#' @param noFolds Mumber of folds for k-fold cross validation.
#' @param reps Number of repetitions for cross validation.
#' @param cores Number of cores for parallelization.
#'
#' @return A matrix of the estimated loss and results of significance testing.
#'
#' @examples
#' # Load libraries
#' library(seminr)
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
#' assess_cvpat(corp_rep_pls_model_ext,
#'              testtype = "two.sided",
#'              BootSamp = 100,
#'              seed = 123,
#'              technique = predict_DA,
#'              noFolds = 10,
#'              reps = 10,
#'              cores = 1)
#'
#' @export
assess_cvpat <- function(model,
                         testtype = "two.sided",
                         BootSamp = 2000,
                         seed = 123,
                         technique = predict_DA,
                         noFolds = NULL,
                         reps = NULL,
                         cores = NULL) {

  set.seed(seed)
  # Abort if received a higher-order-model or moderated model
  if (!is.null(model$hoc)) {
    message("There is no published solution for applying PLSpredict to higher-order-models")
    return()
  }
  if (!is.null(model$interaction)) {
    message("There is no published solution for applying PLSpredict to moderated models")
    return()
  }
  # First we must calculate a IA model which is the "indicator average" model ----
  # we must identify endogenous latents and measures
  endo_lvs <- seminr:::all_endogenous(model$smMatrix)
  endo_mvs <- unlist(lapply(endo_lvs,
                            function(x) seminr:::items_of_construct(construct = x,
                                                                    model = model)))
  # Indicator average (IA) from the training model
  IA <- model$meanData[endo_mvs,drop = F]

  # Calculate IA predictive error
  if (length(endo_mvs) > 1) {
    IA_pred_error <- sweep(model$data[,endo_mvs,drop = F],2 ,IA)
  }
  if (length(endo_mvs) == 1) {
    IA_pred_error <- model$data[,endo_mvs,drop = F] - IA
  }

  # Calculate PLS and LM predictions
  pls_predict_model <- predict_pls(model = model,
                                   technique = technique,
                                   noFolds = noFolds,
                                   reps = reps,
                                   cores = cores)
#####
  pls_predict_error <- as.matrix(pls_predict_model$PLS_out_of_sample_residuals)
  colnames(pls_predict_error) <- endo_mvs
  LM_predict_error <- as.matrix(pls_predict_model$lm_out_of_sample_residuals)
  colnames(LM_predict_error) <- endo_mvs
#####
  PLS_predict_error <- pls_predict_error[,endo_mvs,drop = F]
  LM_predict_error <- LM_predict_error[,endo_mvs,drop = F]

  # Calculate LV-specific losses
  ## for IA model
  LV_losses_IA <- do.call("cbind", lapply(endo_lvs,
                                          function(x) lv_loss(construct = x,
                                                              model = model,
                                                              error = IA_pred_error)))
  ## for LM model
  LV_losses_LM <- do.call("cbind", lapply(endo_lvs,
                                          function(x) lv_loss(construct = x,
                                                              model = model,
                                                              error = LM_predict_error)))
  ## for PLS model
  LV_losses_PLS <- do.call("cbind", lapply(endo_lvs,
                                           function(x) lv_loss(construct = x,
                                                               model = model,
                                                               error = PLS_predict_error)))
  # Name LVs
  colnames(LV_losses_IA) <-  colnames(LV_losses_LM) <- colnames(LV_losses_PLS) <- endo_lvs

  # Calculate overall loss
  ## for IA model
  IA_overall <- overall_loss(LV_losses_IA)
  ## for LM model
  LM_overall <- overall_loss(LV_losses_LM)
  # for PLS model
  PLS_overall <- overall_loss(LV_losses_PLS)

  # CVPAT: PLS vs IA overall
  PLS_v_IA_overall <- bootstrap_cvpat(PLS_overall,
                                      IA_overall,
                                      testtype = testtype,
                                      BootSamp = BootSamp)

  # CVPAT: PLS vs LM overall
  PLS_v_LM_overall <- bootstrap_cvpat(LossM1 = PLS_overall,
                                      LossM2 = LM_overall,
                                      testtype = testtype,
                                      BootSamp = BootSamp)
  ia_cvpat <- cvpat_per_construct(loss_one = LV_losses_PLS,
                                  loss_two = LV_losses_IA,
                                  testtype = testtype,
                                  BootSamp = BootSamp)
  lm_cvpat <- cvpat_per_construct(loss_one = LV_losses_PLS,
                                  loss_two = LV_losses_LM,
                                  testtype = testtype,
                                  BootSamp = BootSamp)

  mat_one <- cbind(colMeans(LV_losses_PLS),colMeans(LV_losses_LM),
                   colMeans(LV_losses_PLS) - colMeans(LV_losses_LM),
                   lm_cvpat[,-1,drop = F])


  colnames(mat_one)[1:3] <- c("PLS Loss", "LM Loss", "Diff")
  mat_one
  mat_one <- rbind(mat_one, c(mean(PLS_overall),
                              mean(LM_overall),
                              mean(PLS_overall) -  mean(LM_overall),
                              PLS_v_LM_overall))
  rownames(mat_one)[nrow(mat_one)] <- "Overall"
  mat_one <- mat_one[,c(1,2,3,6,7)]

  mat_two <- cbind(colMeans(LV_losses_PLS),colMeans(LV_losses_IA),
                   colMeans(LV_losses_PLS) - colMeans(LV_losses_IA),
                   ia_cvpat[,-1,drop = F])
  colnames(mat_two)[1:3] <- c("PLS Loss", "IA Loss", "Diff")
  PLS_v_IA_overall

  mat_two <- rbind(mat_two, c(mean(PLS_overall),
                              mean(IA_overall),
                              mean(PLS_overall) -  mean(IA_overall),
                              PLS_v_IA_overall))

  rownames(mat_two)[nrow(mat_two)] <- "Overall"
  mat_two <- mat_two[,c(1,2,3,6,7)]
  return(list(CVPAT_compare_LM = mat_one,
              CVPAT_compare_IA = mat_two))

}

#function to apply bootstrap to every LV
cvpat_per_construct <- function(loss_one,
                                loss_two,
                                testtype = "two.sided",
                                BootSamp = 2000) {

  index <- colnames(loss_one)
  results <- matrix(nrow = 0, ncol = 6)
  colnames(results) <- c("Construct","Std. T value", "Std. P value", "Boot T value", "Boot P Value", "Perc. P Value")
  for (iter in index) {
    results <- rbind(results,c(iter,bootstrap_cvpat(loss_one[,iter],
                                                    loss_two[,iter],
                                                    testtype = testtype,
                                                    BootSamp = BootSamp)))
  }
  return(results)
}

lv_loss <- function(construct, model, error) {

  if(length(dim((error))) > 1) {
    loss <- rowMeans(error[,seminr:::items_of_construct(construct = construct,
                                                        model = model),drop = F]^2)
  }
  if (length(dim((error))) == 1) {
    loss <- error[,seminr:::items_of_construct(construct = construct,
                                               model = model),drop = F]^2
  }
  return(loss)
}

overall_loss <- function(error) {
  if(length(dim((error))) > 1) {return(rowMeans(error))}
  return(error)
}

bootstrap_cvpat <- function(LossM1,
                            LossM2,
                            testtype = "two.sided",
                            BootSamp = 2000) {

  N <- length(LossM1)
  OrgTtest <- t.test(LossM2,
                     LossM1,
                     alternative = testtype,
                     paired=TRUE)$statistic

  # Originial average difference in losses
  OrgDbar <- mean(LossM2-LossM1)

  # Differences in loss functions under the null
  D_0 <- LossM2-LossM1-OrgDbar

  # Differences in loss functions
  D <- LossM2-LossM1

  #Allocating memory to bootrap
  BootSample <- matrix(0,ncol=2,nrow=length(D))
  BootDbar <- rep(0,BootSamp)
  m_losses<-cbind(LossM1,LossM2)
  tStat <- rep(0,BootSamp)
  for (b in 1:BootSamp) {
    BootSample <- m_losses[sample((1:length(D)), length(D), replace=TRUE),]
    tStat[b]<-t.test(BootSample[,2],BootSample[,1],mu=mean(D),alternative=testtype, paired=TRUE)$statistic
    BootDbar[b] <- mean(sample(D_0, length(D_0), replace=TRUE))
  }
  SorttStat<-sort(tStat, decreasing = FALSE)
  SortBootDbar<-sort(BootDbar, decreasing = FALSE)
  # Bootstrap variance on Dbar for t-test
  std<-sqrt(var(BootDbar))
  tstat_boot_Var<-OrgDbar/std
  # Calculating p-values
  if (testtype=="two.sided") {
    p.value_perc_Ttest<-(sum(SorttStat>abs(OrgTtest))+sum(SorttStat<=(-abs(OrgTtest))))/BootSamp
    p.value_perc_D<-(sum(SortBootDbar>abs(OrgDbar))+sum(SortBootDbar<=(-abs(OrgDbar))))/BootSamp
    p.value_var_ttest<-2*pt(-abs(tstat_boot_Var),(N-1), lower.tail = TRUE)
  }
  if (testtype=="greater") {
    p.value_perc_Ttest<-1-(head(which(SorttStat>OrgTtest),1)-1)/(BootSamp+1)
    p.value_perc_D<-1-(head(which(SortBootDbar>OrgDbar),1)-1)/(BootSamp+1)
    p.value_var_ttest<-pt(tstat_boot_Var,(N-1), lower.tail = FALSE)
    if (length(which(SorttStat>OrgTtest))==0){
      p.value_perc_Ttest=0
      p.value_perc_D=0
    }
  }
  # Load outputs 1
  Results1 <- c("Std. T value"=OrgTtest, "Std. P value"=p.value_perc_Ttest,
                "Boot T value" = tstat_boot_Var, "Boot P Value" = p.value_var_ttest,
                "Perc. P Value" = p.value_perc_D)
  return(Results1)
}

