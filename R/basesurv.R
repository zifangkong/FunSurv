#' A function to obtain the baseline survival function
#'
#' @param A AR1 FRAILTY object
#'
#' @returns a data frame including time and baseline survival
#' @export

basesurv <- function(object){
   cbind(object$time, object$basesurv)
}
