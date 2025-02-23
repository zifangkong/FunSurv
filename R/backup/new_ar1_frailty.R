surv_data = read.table("../data/surv_data.dat", header=TRUE)
fun_data = read.table("../data/fun_data.dat", header=TRUE)
scores = AR1_PACE(fun_data, surv_data, nbasis=10, pve=0.9)

# estimate theta2 and rho (random parameters in the survival part)
calculate_variance_components <- function(ni, tau, rho, J, K) {
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


# get_G_inv <- function(rho0, I, J, K) {
#    IR <- (1 + rho0 ^ 2) * I - rho0 * J - rho0 ^ 2 * K # IR is G^{-1}
#    difIR <- 2 * rho0 * I - J - 2 * rho0 * K # first derivative of IR wrt rho
#    list(IR = IR, difIR = difIR)
# }

# I J K matrices are used for the estimation of AR(1) correlation structure
IJK <- function(ni) {
   N <- sum(ni) # total number of observations
   
   ## I
   I <- diag(N)
   
   ## J
   J <- diag(0, N)
   diag(J[,-1]) <- 1 
   diag(J[-1,]) <- 1
   
   ## K
   K <- diag(N)
   counter <- c(0,cumsum(ni))
   for(k in 1:length(ni)){
      if(ni[k]==1){ K[counter[k]+1,counter[k]+1] = 2 }
      else if(ni[k]>2) { K[counter[k]+1+1:(ni[k]-2),counter[k]+1+1:(ni[k]-2)] = 0 }
   }
   list(I = I, J = J, K = K) # all have dimension N by N
}


# breslow <- function(W, indi, eta) {
#    survbase <- exp(-W %*% (indi / (t(W) %*% (exp(eta)))))
#    list(survbase = survbase)
# }


AR1_frailty <- function(formula, data, scores, para0, itmax, eps.reg) {
   
   if (missing(formula)) stop("A Event formula is required.")
   
   terms_obj <- terms(formula)
   response <- as.character(attr(terms_obj, "variables")[[2]])[-1]
   covariates <- attr(terms_obj, "term.labels")
   
   # extract data
   DF <- model.frame(reformulate(c(response, covariates)), data)
   
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
   X_surv <- as.matrix(DF_sorted[, c(covariates, paste0("score", seq(ncol(scores))))])  # Design matrix
   R <- as.matrix(DF_sorted[, (ncol(DF_sorted)-ncol(R)+1):ncol(DF_sorted)])
   XX <- cbind(X_surv, R)
   W <- matrix(1, N, N)
   W[upper.tri(W)] <- 0
   H22 <- diag(0, (p + N)) ## initial of Hessian matrix
   
   ## Initial values
   beta0 <- rep(0, p)
   V0 <- rep(0, N)
   par0 <- c(beta0, V0)
   rho0 <- params[1]
   theta20 <- params[2]
   eps <- 1e-6
   convergence <- 0

   
   for (outer.iter in 1:itmax) {
      # Optimize beta and V
      AR_inv <- (1+rho0^2) * ijk$I - rho0 * ijk$J - rho0^2 * ijk$K
      H22[(p + 1):(p + N), (p + 1):(p + N)] <- AR_inv / theta20
      
      for (inner.iter in 1:itmax) {
         eta <- as.vector(X_surv %*% beta0 + R %*% V0)
         
         ######################################################
         w <- diag(as.vector(exp(eta)))
         A <- diag(as.vector(indi / (t(W) %*% (exp(eta)))))
         B <- diag(as.vector(W %*% A %*% rep(1, N)))
         dll.eta <- w %*% B - w %*% W %*% A %*% A %*% t(W) %*% w # second derivative wrt eta
         dl.eta <- as.vector(indi - w %*% W %*% A %*% rep(1, N)) # first derivative wrt eta
         dl.dbeta <- t(X_surv) %*% dl.eta
         dl.dV <- crossprod(R, dl.eta) - (1 / theta20) * (AR_inv %*% V0)
         H2 <- solve(crossprod(XX, dll.eta) %*% XX + H22)
         Svec <- as.vector(c(dl.dbeta, dl.dV))
         par <- par0 + H2 %*% Svec
         if (max(abs(par - par0)) < eps) {
            convergence <- 1
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
   basesurv <- exp(-W %*% (indi / (t(W) %*% (exp(eta)))))
   
   ## Calculate variance-covariance matrix using numerical Hessian
   se.beta <- sqrt(abs(diag(H2)[1:p]))
   ebeta <- cbind(beta, se.beta, 2*(1-pnorm(abs(beta/se.beta))))
   dimnames(ebeta) <- list(covariates, colnames(scores)), c("Estimate", "SE", "p-value"))
   
   ## SE for the variance components


   --------
   difIR <- ob3$difIR
   K1.lat <- (H2[(p + 1):(p + N), (p + 1):(p + N)] %*% AR_inv) / theta2
   K2.lat <- (H2[(p + 1):(p + N), (p + 1):(p + N)] %*% difIR) / theta2
   K3.lat <- solve(AR_inv) %*% difIR
   
   b11 <- sum(diag((diag(N) - K1.lat) %*% (diag(N) - K1.lat))) / theta2^2
   b12 <- -sum(diag((diag(N) - K1.lat) %*% (diag(N) - K1.lat) %*% K3.lat)) / theta2
   b21 <- b12
   b22 <- sum(diag((K2.lat - K3.lat) %*% (K2.lat - K3.lat)))
   varmat <- 2 * solve(matrix(c(b11, b12, b21, b22), ncol = 2))
   stdvar <- cbind(c(theta2, rho2), sqrt(diag(varmat)), 2 * (1 - pnorm(abs(c(theta2, rho2) / sqrt(diag(varmat))))))
   dimnames(stdvar) <- list(c("theta2", "rho2"), c("estimate", "SE", "p-value"))

   
   return(list(
      beta = beta,
      V = V,
      theta2 = theta2,
      rho = rho,
      logLik = -fit$value
   ))
}


summary.ar1frailty <- function(object, ...) {
   se <- sqrt(diag(object$vcov))
   z <- object$coefficients/se[1:length(object$coefficients)]
   p <- 2 * (1 - pnorm(abs(z)))
   
   coef_table <- cbind(
      Estimate = object$coefficients,
      "Std. Error" = se[1:length(object$coefficients)],
      "z value" = z,
      "Pr(>|z|)" = p
   )
   
   results <- list(
      coefficients = coef_table,
      phi = c(object$phi, se[length(se)-1]),
      theta = c(object$theta, se[length(se)]),
      loglik = object$loglik,
      convergence = object$convergence
   )
   
   class(results) <- "summary.ar1frailty"
   return(results)
}
