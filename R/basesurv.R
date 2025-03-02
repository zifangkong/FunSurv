#' A function to obtain the baseline survival function
#'
#' @param object An AR1 FRAILTY object
#'
#' @returns A data frame including time and baseline survival
#' @export

basesurv <- function(object){
  if (!is.funsurv(object)) stop("Must be a funsurv object")
   data.frame(time = object$time, basesurv = object$basesurv)
}
