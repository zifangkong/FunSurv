#' @noRd
get_G_inv <- function(rho20, I, J, K) {
  IR <- (1 + rho20 ^ 2) * I - rho20 * J - rho20 ^ 2 * K # IR is G^{-1}
  difIR <- 2 * rho20 * I - J - 2 * rho20 * K # first derivative of IR wrt rho
  list(IR = IR, difIR = difIR)
}
