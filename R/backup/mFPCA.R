#' multivariate FPCA based on results from uPACE
#'
#' @param
#' 
mFPCA <- function(PACE1, PACE2, I) {
  ## urineList <- c('upH', 'uCit')
  scores <- cbind(PACE1$scores, PACE2$scores)
  npc <- c(PACE1$npc, PACE2$npc)  
  phi <- NULL
  phi[[1]] <- t(PACE1$functions@X)
  phi[[2]] <- t(PACE2$functions@X)  
  Xi <- scores
  phi <- phi
  p <- length(urineList)
  L <- npc  
  ## eigenanalysis on matrix M
  M <- t(Xi) %*% Xi / (I - 1)
  # M <- (Xi)%*%t(Xi)/(I-1)
  eigen.M <- eigen(M)
  values <- eigen.M$values # eigenvalues
  pve <- cumsum(values) / sum(values)
  Cms <- eigen.M$vectors # eigenvectors
  index <- unlist(lapply(1:length(L), function(x) rep(x, L[x])))
  ## MFPCA score
  rho <- mfpca.score(predXi = Xi, Cms)  
  ## MFPCA eigenfunction
  psis <- NULL
  for (j in 1:p) {
    psi <- NULL
    for (m in 1:ncol(Cms)) {
      psi <- cbind(psi, phi[[j]] %*% Cms[which(index == j), m])
    }
    psis[[j]] <- psi
  }
  ## out <- list(eigenvalue <- values, Cms <- Cms, pve <- pve, index=index, rho <- rho, psis=psis)
  out <- list(
    values = values,
    functions = psis,
    scores = rho,
    npc = npc
  )
  return(out)
}
