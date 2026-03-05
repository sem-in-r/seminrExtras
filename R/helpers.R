#' @importFrom stats pt
#' @importFrom stats t.test
#' @importFrom stats var
#' @importFrom utils head
#' @importFrom seminr predict_DA
#' @importFrom seminr predict_EA
#' @importFrom seminr mean_replacement
#' @importFrom seminr estimate_pls
#' @importFrom seminr predict_pls
#' @importFrom rpart rpart
#' @importFrom graphics abline axis legend matlines matplot
#'   par plot points polygon text
#' @importFrom grDevices palette
#'
NULL

#' @noRd
validate_seminr_model <- function(model, func_name = "This function") {
  if (!any(class(model) == "seminr_model")) {
    warning(func_name, " only works with SEMinR models.", call. = FALSE)
    return(FALSE)
  }
  return(TRUE)
}

#' @noRd
has_higher_order <- function(model) {
  !is.null(model$hoc)
}

#' @noRd
validate_for_prediction <- function(model, func_name = "This function") {
  if (!validate_seminr_model(model, func_name)) {
    return(FALSE)
  }
  if (has_higher_order(model)) {
    warning("There is no published solution for applying PLSpredict to higher-order models.",
            call. = FALSE)
    return(FALSE)
  }
  return(TRUE)
}
