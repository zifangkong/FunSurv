#' Fit the Functional Regression with AutoregressIve fraiLTY (FRAILTY) model
#' 
#' @param formula A Event formula, with the response on the left of a ~ operator, 
#' and scalar covariates on the right. The response must be a Event object using Event function. 
#' @param data A data frame, including time-to-event outcomes, and scalar coefficients
#' @param fpca_obj A FPCA object returned by the AR1_PACE function
#' @param para0 Initial values for \eqn{\theta^2} and auto-regressive coefficient \eqn{\rho}.
#' @param iter.max Maximum number of iterations for both inner iteration and outer iteration. Defaults to \code{50}.
#' @param eps Tolerance criteria for a possible infinite coefficient value. Defaults to \code{1e-6}.
#'
#' @importFrom MASS ginv
#' @importFrom Matrix bdiag
#' @returns A AR1_FRAILTY object
#' @export
#'
#' @examples inst/examples/ex_AR1_FRAILTY.R
AR1_FRAILTY <- function(formula, 
                        data, 
                        fpca_obj, 
                        para0, 
                        iter.max = 50, 
                        eps = 1e-6) {  
  if (missing(formula)) stop("A Event formula is required.")
  if (missing(data)) stop("A dataset of time-to-event outcomes is required.")
  if (missing(fpca_obj)) stop("A FPC object, returned by the AR1_PACE function is required.")  
  if (missing(para0)) stop("Initial values for theta^2 and rho are required") 
  terms_obj <- terms(formula)
  response <- as.character(attr(terms_obj, "variables")[[2]])[-1]
  covariates <- attr(terms_obj, "term.labels")
  ## extract data
  DF <- model.frame(reformulate(c(response, covariates)), data)
  scores <- fpca_obj$window_scores_fpca
  if(nrow(DF) != nrow(scores)) stop("Dimension of data does not match with the dimension of FPC scores")
  DF <- cbind(DF, scores)
  colnames(DF) <- c("t_start", "t_stop", "id", "status", covariates, paste0("score", seq(ncol(scores))))
  DF$gap_time = DF$t_stop - DF$t_start  
  p <- length(covariates) + ncol(scores)  # mumber of predictors in the survival model
  N <- nrow(DF)  # total number of event observations
  M <- length(unique(DF$id))  # number of unique subjects
  R <- diag(N)  # Identity matrix
  DF <- cbind(DF, R)  
  ni <- table(DF$id)  # Observations per subject
  ijk <- IJK(ni)
  DF_sorted <- DF[order(DF$gap_time), ]  # Sort by gap time
  indi <- as.vector(DF_sorted$status)  # Event indicator
  X_surv <- as.matrix(DF_sorted[, c(covariates, paste0("score", seq(ncol(scores))))]) 
  R <- as.matrix(DF_sorted[, (ncol(DF_sorted)-ncol(R)+1):ncol(DF_sorted)])
  XX <- cbind(X_surv, R) # design matrix
  W <- matrix(1, N, N)
  W[upper.tri(W)] <- 0
  H22 <- diag(0, (p + N)) ## initial of Hessian matrix  
  ## Initial values
  beta0 <- rep(0, p)
  V0 <- rep(0, N)
  par0 <- c(beta0, V0)
  rho0 <- para0[1]
  theta20 <- para0[2]
  convergence <- 0
  for (outer.iter in 1:iter.max) {
    ## Optimize beta and V
    AR_inv <- (1+rho0^2) * ijk$I - rho0 * ijk$J - rho0^2 * ijk$K
    H22[(p + 1):(p + N), (p + 1):(p + N)] <- AR_inv / theta20    
    for (inner.iter in 1:iter.max) {
      eta <- as.vector(crossprod(t(XX), par0))
      ## ####################################################
      w <- diag(as.vector(exp(eta)))
      A <- diag(as.vector(indi / crossprod(W, exp(eta)) ))
      B <- diag(as.vector(W %*% A %*% rep(1, N)))
      ## dll.eta <- w %*% B - w %*% W %*% A %*% A %*% t(W) %*% w # second derivative wrt eta
      ## dl.eta <- as.vector(indi - w %*% W %*% A %*% rep(1, N)) # first derivative wrt eta
      tmp <- w %*% W %*% A
      dll.eta <- crossprod(w, B) - tcrossprod(tmp)
      dl.eta <- as.vector(indi - rowSums(tmp))
      dl.dbeta <- crossprod(X_surv, dl.eta)
      dl.dV <- crossprod(R, dl.eta) - (1 / theta20) * tcrossprod(AR_inv, t(V0))
      H2 <- solve(crossprod(XX, dll.eta) %*% XX + H22)
      Svec <- as.vector(c(dl.dbeta, dl.dV))
      par <- par0 + tcrossprod(H2, t(Svec))
      if (max(abs(par - par0)) < eps) {
        break
      }      
      ## Update beta and V
      par0 <- par
      beta0 <- beta <- par[1:p]
      V0 <- V <- par[(p + 1):(p + N)]
    }
    ## UPDATE rho and theta2
    tau <- tcrossprod(V) + H2[(p + 1):(p + N), (p + 1):(p + N)]
    ar_var <- calculate_variance_components(ni, tau, rho0, ijk$J, ijk$K)
    rho <- ar_var$rho
    theta2 <- ar_var$theta2    
    if (pmax(abs(theta2 - theta20), abs(rho - rho0)) < eps) {
      convergence <- 1
      break
    }
    theta20 <- theta2
    rho0 <- rho
  }
  ## baseline survival function
  basesurv <- exp(-W %*% (indi / crossprod(W, exp(eta)) ))
  ## Calculate variance-covariance matrix using numerical Hessian
  se.beta <- sqrt(abs(diag(H2)[1:p]))
  ebeta <- cbind(beta, se.beta, 2*(1-pnorm(abs(beta/se.beta))))
  dimnames(ebeta) <- list(c(covariates, colnames(scores)), c("Estimate", "SE", "p-value"))
  ## SE for the variance components
  block_list <- lapply(ni, function(ni) ar1_cor(ni, rho))
  U <- bdiag(block_list)  # Sparse block diagonal matrix
  U <- as.matrix(U)
  Q1 <- ginv(dll.eta) + theta2 * R %*% U %*% t(R)
  U.inv <- bdiag(lapply(ni, function(ni) dar1_cor.drho(ni, rho))) 
  U.inv <- as.matrix(U.inv)
  dQ1.theta2 <- R %*% U %*% t(R)
  dQ1.rho <- theta2 * R %*% U.inv %*% t(R)
  Q2 <- solve(Q1)  - solve(Q1) %*% X_surv %*% solve(t(X_surv) %*% solve(Q1) %*% X_surv) %*% t(X_surv) %*% solve(Q1)
  b11 <- sum(diag(Q2 %*% dQ1.theta2 %*% Q2 %*% dQ1.theta2))
  b12 <- sum(diag(Q2 %*% dQ1.theta2 %*% Q2 %*% dQ1.rho))
  b22 <- sum(diag(Q2 %*% dQ1.rho %*% Q2 %*% dQ1.rho))
  
  varmat <- 2 * solve(matrix(c(b11, b12, b12, b22), ncol = 2))
  eAR_var <- cbind(c(theta2, rho), sqrt(diag(varmat)), 2 * (1 - pnorm(abs(c(theta2, rho) / sqrt(diag(varmat))))))
  dimnames(eAR_var) <- list(c("theta2", "rho2"), c("estimate", "SE", "p-value"))
  
  print(ebeta)
  cat("-------------------------------\n")
  print(round(eAR_var, 3))
  out <- list(ebeta = ebeta, V = V, eAR_var = eAR_var,
              basesurv=basesurv, time=sort(DF$t_stop), PACE=fpca_obj$uni.PACE)
  class(out) <- "funsurv"
  return(out)
}


                                        # 
#' Estimate the variance component \eqn{\theta^2} and auto-regressive coefficient \eqn{\rho}
#'
#' @param ni A vector representing the number of events for each subject.
#' @param tau 
#' @param rho \eqn{\rho} value at previous iteration
#' @param J J matrix, return by IJK function
#' @param K K matrix, returned by IJK function
#'
#' @returns estimation for \eqn{\theta^2} and \eqn{\rho}
#' @noRd
calculate_variance_components <- function(ni, tau, rho, J, K) {
  N <- sum(ni)
  M <- length(ni)
  ## calculate A1, A2, A3
  A1 <- sum(diag(tau))
  A2 <- sum(diag(J %*% tau)) / 2
  A3 <- sum(diag(K %*% tau))
  ## calculate B1, B2, B3, B4
  B1 <- (N - M) * (A1 - A3)
  B2 <- (2 * M - N) * A2
  B3 <- N * A3 - (N + M) * A1
  B4 <- N * A2
  ## estimate rho and theta
  rho0 <- (rho - (B1*rho^3 + B2*rho^2 + B3*rho + B4) / (3*B1*rho^2 + 2*B2*rho + B3))
  theta20 <- (1/N) * ((1 + rho0^2)*A1 - 2*rho0*A2 - (rho0^2)*A3)
  list(rho = rho0, theta2 = theta20)
}


#' Estimate I J K matrices are used for the estimation of AR(1) correlation structure
#'
#' @param ni A vector representing the number of events for each subject.
#'
#' @returns I, J, K matrices
#' @noRd
IJK <- function(ni){
  N <- sum(ni)
  M <- length(ni)
  counter <- c(0, cumsum(ni))
  ## I
  I <- diag(N) 
  ## J
  J <- diag(0,N)
  if(max(ni)>1){
    maxJ <- diag(0,max(ni))
    diag(maxJ[,-1]) <- 1
    diag(maxJ[-1,]) <- 1
    for(i in 1:M){
      if(ni[i]>1) J[counter[i]+1:ni[i],counter[i]+1:ni[i]] <- maxJ[1:ni[i],1:ni[i]]
    }
  }  
  ## K
  K <- diag(N)
  counter <- c(0,cumsum(ni))
  for(k in 1:length(ni)){
    if(ni[k]==1){ K[counter[k]+1,counter[k]+1] = 2 }
    else if(ni[k]>2) { K[counter[k]+1+1:(ni[k]-2),counter[k]+1+1:(ni[k]-2)] = 0 }
  }  
  list(I=I,J=J,K=K)
}

is.funsurv <- function(x) inherits(x, "funsurv")
