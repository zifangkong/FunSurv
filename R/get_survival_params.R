#' Function to get the random parameters in the survival part
#' @noRd
get_survival_params <- function(tau2, ni, rho20, J, K) {
  ni <- as.vector(ni) # ni
  M <- length(ni) # number of subjects
  N <- sum(ni) # total number of observations
  ## calculate B1,B2,B3
  B1 <- sum(diag(tau2))
  B2 <- sum(diag(J %*% tau2)) / 2
  B3 <- sum(diag(K %*% tau2))
  ## calculate the phi and theta
  A1 <- (N - M) * (B1 - B3)
  A2 <- (2 * M - N) * B2
  A3 <- N * B3 - (N + M) * B1
  A4 <- N * B2
  rh2 <- (rho20 - (A1 * rho20 ^ 3 + A2 * rho20 ^ 2 + A3 * rho20 + A4) / (3 * A1 * rho20 ^ 2 + 2 * A2 * rho20 + A3))
  th2 <- (1 / (N)) * ((1 + rh2 ^ 2) * B1 - 2 * rh2 * B2 - (rh2 ^ 2) * B3)
  list(rho2 = rh2, theta2 = th2)
}

