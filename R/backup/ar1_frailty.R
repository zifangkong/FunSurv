#' This looks like one of the main functions
#' Need to modify the function so that predictors can be included in function arguments
#'
#' @param dat a dataset with columns c("id", "gap_time","status", scalar covariate, window-specific FPC scores)
#' @param theta20 initial value for variance component in the AR(1) structure theta
#' @param rho20 initial value for autoregressive coefficient rho
#' @param itmax maximum number of iteration
#'
#' @returns estimated beta, alpha, rho, theta, and their corresponding standard deviations
#'
AR1_frailty <- function(dat, theta20, rho20, itmax) {
  p <- ncol(dat) - 3 # 3 for c("id", "gap_time","status")  # number of predictors in the survival model
  N <- nrow(dat) # total number of event observations
  patient <- dat[, 1]
  M <- length(unique(patient)) #  number of unique subjects
  R <- diag(N) # N by N identity matrix
  ni <- table(dat[, 1]) # number of observations per subject
  ijk <- IJK(ni)
  RR <- cbind(dat, R)
  RR_sorted <- RR[sort.list(RR[, 2]), ] # sort by gap time (to get distinct gap time)
  indi <- as.vector(RR_sorted[, 3]) # corresponding event indicator (1 for event, 0 for censored)
  R <- as.matrix(RR_sorted[, -(1:(3 + p))]) # for the hazard model (random effect design matrix) (sorted based on gap time)
  X_surv <- matrix(as.numeric(unlist(RR_sorted[, 4:(3 + p)])), ncol = p) # design matrix
  M1 <- matrix(1, N, N)
  M1[upper.tri(M1)] <- 0
  
  for (ii in 1:N) {
    if (indi[N - ii + 1] == 1) {
      largest <- N - ii + 1
      break
    } # largest index of event outcome
  }
  
  ## initial values
  beta0 <- as.vector(rep(0, p)) # beta: the coefficients in the hazard model without intercept term
  V <- as.vector(rep(0, N)) # frailty terms
  par0 <- as.vector(c(beta0, V)) # parameters in the hazard model
  eta <- as.vector(X_surv %*% beta0 + R %*% V)
  flag.var <- 0
  eps.reg <- eps.var <- 0.00001
  
  
  for (outer.iter in 1:itmax) {
    ## itmax is user-defined max number of iterations
    flag.reg <- 0
    ob3 <- get_G_inv(rho20, ijk$I, ijk$J, ijk$K)
    UG2 <- diag(0, (p + N))
    UG2[(p + 1):(p + N), (p + 1):(p + N)] <- ob3$IR / theta20
    
    
    for (inner.iter in 1:itmax) {
      eta <- as.vector(X_surv %*% beta0 + R %*% V)
      
      ######################################################
      w <- diag(as.vector(exp(eta)))
      A <- diag(as.vector(indi / (t(M1) %*% (exp(eta)))))
      B <- diag(as.vector(M1 %*% A %*% rep(1, N)))
      f2.eta <- w %*% B - w %*% M1 %*% A %*% A %*% t(M1) %*% w # second derivative
      f1.eta <- as.vector(indi - w %*% M1 %*% A %*% rep(1, N)) # first derivative
      dl.dbeta <- t(X_surv) %*% f1.eta
      dl.dV <- t(R) %*% f1.eta - (1 / theta20) * (ob3$IR %*% V)
      XX2 <- cbind(X_surv, R)
      H2 <- solve(t(XX2) %*% f2.eta %*% XX2 + UG2)
      Svec <- as.vector(c(dl.dbeta, dl.dV))
      par <- par0 + H2 %*% Svec
      ## update initial values
      if (max(abs(par - par0)) < eps.reg) {
        flag.reg <- 1
        break
      }
      par0 <- par
      beta0 <- par[1:p]
      V <- par[(p + 1):(p + N)]
    }
    
    
    
    if (flag.reg == 0)
      stop("not reach convergence")
    flag.reg = 0
    
    ####### Variance parameter##################
    tau2 <- (V %*% t(V) + H2[(p + 1):(p + N), (p + 1):(p + N)])
    objQ <- get_survival_params(tau2, ni, rho20, ijk$J, ijk$K)
    rho2 <- objQ$rho2
    theta2 <- objQ$theta2
    if (max(abs(c((theta2 - theta20), (rho2 - rho20)))) < eps.var) {
      flag.var <- 1
      break
    }
    theta20 <- theta2
    rho20 <- rho2
    bres.etail <- breslow(M1, indi, eta)
    survbase <- bres.etail$survbase
  }
  
  
  
  
  bres.etail <- breslow(M1, indi, eta)
  survbase <- bres.etail$survbase
  eta0 <- as.vector(X_surv[, 1] * beta0[1] + R %*% V)
  LR <- 2 * (sum(indi * eta) - log(indi %*% M1 %*% exp(eta))) - 2 * (sum(indi * eta0) - log(indi %*% M1 %*% exp(eta0)))
  
  ###############################################################
  ## SE for the variance components
  IR <- ob3$IR
  difIR <- ob3$difIR
  K1.lat <- (H2[(p + 1):(p + N), (p + 1):(p + N)] %*% IR) / theta2
  K2.lat <- (H2[(p + 1):(p + N), (p + 1):(p + N)] %*% difIR) / theta2
  K3.lat <- solve(IR) %*% difIR
  b11 <- sum(diag((diag(N) - K1.lat) %*% (diag(N) - K1.lat))) / theta2^2
  b12 <- -sum(diag((diag(N) - K1.lat) %*% (diag(N) - K1.lat) %*% K3.lat)) / theta2
  b21 <- b12
  b22 <- sum(diag((K2.lat - K3.lat) %*% (K2.lat - K3.lat)))
  varmat <- 2 * solve(matrix(c(b11, b12, b21, b22), ncol = 2))
  se.var <- sqrt(diag(varmat))
  stdvar <- cbind(c(theta2, rho2), sqrt(diag(varmat)), 2 * (1 - pnorm(abs(c(theta2, rho2) / sqrt(diag(varmat))))))
  dimnames(stdvar) <- list(c("theta2", "rho2"), c("estimate", "s.e.", "p-value"))
  stdvar <- round(stdvar, 3)
  ## #########################################
  beta <- par0[1:p]
  se.beta <- sqrt(abs(diag(H2)[1:p]))
  ## #################
  ebeta <- cbind(beta, se.beta, 2 * (1 - pnorm(abs(beta / se.beta))))
  names(beta) <- colnames(dat)[4:ncol(dat)]
  dimnames(ebeta) <- list(names(beta), c("Estimate", "SE", "p-value"))
  score_stat <- t(dl.dbeta[2:p]) %*% (t(XX2) %*% f2.eta %*% XX2 + UG2)[2:p, 2:p] %*% (dl.dbeta[2:p])
  return(list(ebeta = ebeta,
              stdvar = stdvar,
              df = p - 1,
              survbase = survbase[rank(RR[, 2])]
  ))
}