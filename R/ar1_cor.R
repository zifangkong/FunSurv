#' Construct an AR(1) correlation matrix
#'
#' @param n number of events for each subject
#' @param rho autoregressive correlation
#'
#' @returns A n by n matrix
#' 
#' @export
#' @example inst/examples/ex_ar1_cor.R
ar1_cor <- function(n, rho) {
  ## Check if rho is between 0 and 1
  if (rho < 0 || rho > 1) stop("Error: rho must be between 0 and 1.")
  d <- c(n, n)
  rho^abs(.col(d) - .row(d)) / (1 - rho^2)
}


#' First derivative of AR(1) correlation matrix with respect to the auto-regressive coefficient
#'
#' @param n number of events for each subject
#' @param rho autoregressive correlation
#'
#' @returns A n by n inverse matrix
#' @export
dar1_cor.drho <- function(n, rho) {
  ## Check if rho is between 0 and 1
  if (rho < 0 || rho > 1) stop("Error: rho must be between 0 and 1.")
  d <- c(n, n)
  idx <- abs(.col(d) - .row(d))
  numerator <- idx * rho^(idx - 1) * (1 - rho^2) + 2 * rho * rho^(idx)
  return(numerator / (1 - rho^2)^2)
}
