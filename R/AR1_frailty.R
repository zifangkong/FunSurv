#' This looks like one of the main functions
#'
#' @param 
#' 
AR1_frailty <- function(dat, theta20, rho20, itmax, q) {
  p <- ncol(dat) - 3 # 3 (x1, x2, y) # number of predictors in the survival model
  M <- nrow(dat) # M<-2093 observations
  patient <- dat[, 1]
  N <- length(unique(patient)) # N<-500 unique subjects  
  R2 <- diag(M) # M by M identity matrix
  ni <- table(dat[, 1]) # number of observations per subject
  obj3 <- IJK(ni)
  XZ2 <- cbind(dat, R2)
  XZ2.r <- XZ2[sort.list(XZ2[, 2]), ] # sort by gap time (to get distinct gap time)
  time <- as.vector(XZ2.r[, 2]) #sorted gap times
  indi <- as.vector(XZ2.r[, 3]) #corresponding event indicator (1 for event, 0 for censored)
  R2 <- as.matrix(XZ2.r[, -(1:(3 + p))]) # for the hazard model (random effect design matrix) (sorted based on gap time)
  X_surv <- matrix(as.numeric(unlist(XZ2.r[, 4:(3 + p)])), ncol = p) # design matrix
  M1 <- matrix(1, M, M)
  M1[upper.tri(M1)] <- 0  
  for (ii in 1:M) {
    if (indi[M - ii + 1] == 1) {
      largest <- M - ii + 1
      break
    } # largest index of event outcome
  }  
  Constrain <- diag(c(rep(1, largest), rep(0, (M - largest))))
  Constrain0 <- diag(c(rep(0, largest), rep(1, (M - largest))))
  ## initial values
  beta0 <- as.vector(rep(0, p)) # beta is the coefficients in the hazard model without intercept term
  V <- as.vector(rep(0, M)) # M is the number of observations
  par0 <- as.vector(c(beta0, V)) # parameters in the hazard model
  eta <- as.vector(X_surv %*% beta0 + R2 %*% V)
  flag.var <- 0
  eps.reg <- eps.var <- 0.0001
  for (outer.iter in 1:itmax) {
    ## itmax is user-defined max number of iterations
    flag.reg <- 0
    ob3 <- get_G_inv(rho20, obj3$I, obj3$J, obj3$K)
    UG2 <- diag(0, (p + M)) # what are UG1 and UG2??
    UG2[(p + 1):(p + M), (p + 1):(p + M)] <- ob3$IR / theta20
    for (inner.iter in 1:itmax) {
      eta <- as.vector(X_surv %*% beta0 + R2 %*% V)      
      ## ############ Latency part ########################################
      w <- diag(as.vector(exp(eta)))
      A <- diag(as.vector(indi / (t(M1) %*% (exp(eta)))))
      B <- diag(as.vector(M1 %*% A %*% rep(1, M)))
      f2.eta <- w %*% B - w %*% M1 %*% A %*% A %*% t(M1) %*% w # second derivative??
      f1.eta <- as.vector(indi - w %*% M1 %*% A %*% rep(1, M)) # first derivative?
      dl.dbeta <- t(X_surv) %*% f1.eta
      dl.dV <- t(R2) %*% f1.eta - (1 / theta20) * (ob3$IR %*% V)
      XX2 <- cbind(X_surv, R2)
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
      V <- par[(p + 1):(p + M)]
    }    
    if (flag.reg == 0)
      stop("not reach convergence")
    flag.reg = 0    
    ## ##### Variance parameter of the latency part##################
    tau2 <- (V %*% t(V) + H2[(p + 1):(p + M), (p + 1):(p + M)])
    objQ <- get_survival_params(tau2, ni, rho20, obj3$J, obj3$K)
    rho2 <- objQ$rho2
    theta2 <- objQ$theta2
    if (max(abs(c((theta2 - theta20), (rho2 - rho20)))) < eps.var) {
      flag.var <- 1
      break
    }
    theta20 <- theta2
    rho20 <- rho2
    bres.etail <- breslow(Constrain, Constrain0, largest, time, M1, indi, eta)
    survbase <- bres.etail$survbase
  }
  bres.etail <- breslow(Constrain, Constrain0, largest, time, M1, indi, eta)
  survbase <- bres.etail$survbase
  survprob <- survbase ^ (exp(eta))
  eta0 <- as.vector(X_surv[, 1] * beta0[1] + R2 %*% V)
  LR <- 2 * (sum(indi * eta) - log(indi %*% M1 %*% exp(eta))) - 2 * (sum(indi * eta0) - log(indi %*% M1 %*% exp(eta0)))
  ## #############################################################
  ## SE for the variance components
  IR <- ob3$IR
  difIR <- ob3$difIR
  K1.lat <- (H2[(p + 1):(p + M), (p + 1):(p + M)] %*% IR) / theta2
  K2.lat <- (H2[(p + 1):(p + M), (p + 1):(p + M)] %*% difIR) / theta2
  K3.lat <- solve(IR) %*% difIR
  b11 <- sum(diag((diag(M) - K1.lat) %*% (diag(M) - K1.lat))) / theta2^2
  b12 <- -sum(diag((diag(M) - K1.lat) %*% (diag(M) - K1.lat) %*% K3.lat)) / theta2
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
  return(
    list(ebeta = ebeta,
         stdvar = stdvar,
         bres.etail = bres.etail,
         H2 = H2[1:p, 1:p],
         score_stat = score_stat,
         LR = LR,
         df = p - 1,
         survprob = survprob[rank(XZ2[, 2])],
         survbase = survbase[rank(XZ2[, 2])],
         eta = eta[rank(XZ2[, 2])],
         frailty = V)
  )
}
