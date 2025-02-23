#' Construct an AR(1) correlation matrix
#'
#' @param n number of events for each subject
#' @param rho autoregressive correlation
#'
#' @returns A n by n matrix
#' 
#' @export
ar1_cor <- function(n, rho) {
  d <- c(n, n)
  rho^abs(.col(d) - .row(d)) / (1 - rho^2)
}


## old version from https://www.r-bloggers.com/2020/02/generating-correlation-matrix-for-ar1-mode/l
## why? you don't even need the "-1" at these 2 places in exponent
## ar1_cor <- function(n, rho) {
##   exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - (1:n - 1))
##   rho^exponent / (1 - rho^2)
## }

