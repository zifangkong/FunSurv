#' construct the AR(1) correlation structure
#'
#' @param n number of events for each subject
#' @param rho autoregressive correlation
#'
#' @returns A n by n matrix
#' @export
ar1_cor <- function(n, rho) {
  exponent = abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - (1:n - 1))
  return(rho^exponent/(1-rho^2))
}
