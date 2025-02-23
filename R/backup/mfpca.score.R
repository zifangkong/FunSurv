#' mfpc score calculation
mfpca.score <- function(predXi, Cms) {
  rho <- matrix(NA, nrow = nrow(predXi), ncol = dim(Cms)[2])
  for (i in 1:nrow(predXi)) {
    for (m in 1:dim(Cms)[2]) {
      rho[i, m] <- predXi[i, ] %*% Cms[, m]
    }
  }
  return(rho)
}
