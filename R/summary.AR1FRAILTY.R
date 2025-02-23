#' A summary function for the AR1 FRAILTY object
#'
#' @param object A AR1 FRAILTY object
#' @export
summary.AR1FRAILTY <- function(object) {
   print(object$ebeta)
   cat("-------------------------------\n")
   print(object$eAR_var)
}
