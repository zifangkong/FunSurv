#' Function to construct I J K matrices for the estimation of AR(1) correlation structure
#'
#' @noRd
IJK <- function(ni) {
  dp <- max(ni)
  N <- sum(ni) # total number of observations
  I <- K <- diag(N)
  J <- diag(0, N) # all 0 matrix
  ln <- length(ni) # number of subjects
  counter <- c(0, cumsum(ni))
  # K
  if (ln == 1 & dp == 1) {
    K[1, 1] <- 2
  } else{
    for (i in 1:ln) {
      if (ni[i] == 1) {
        K[counter[i] + 1, counter[i] + 1] <- 2
      }
      else if (ni[i] > 2) {
        K[counter[i] + 1 + 1:(ni[i] - 2), counter[i] + 1 + 1:(ni[i] - 2)] <- 0
      }
    }
  }
  # J
  if (dp > 1) {
    Jp <- diag(0, dp)
    for (j in 2:dp) Jp[j - 1, j] <- 1
    for (i in 2:dp) Jp[i, i - 1] <- 1
    for (i in 1:ln) {
      if (ni[i] > 1)
        J[counter[i] + 1:ni[i], counter[i] + 1:ni[i]] <- Jp[1:ni[i], 1:ni[i]]
    }
  }
  list(I = I, K = K, J = J) # all have dimension N by N
}
