#' Fit a Functional Regression with AutoregressIve fraiLTY (FRAILTY) model for Recurrnet Event Data
#' 
#' Jointly model longitudinal measurements and recurrent events, 
#' accommodating both scalar and functional predictors while capturing time-dependent correlations among events. 
#' The FRAILTY method employs a two-step estimation procedure. 
#' First, functional principal component analysis through conditional expectation (PACE) is applied to extract 
#' key temporal features from sparse, irregular longitudinal data. 
#' Second, the obtained scores are incorporated into a dynamic recurrent frailty model 
#' with an autoregressive structure to account for within-subject correlations across recurrent events.
#' This function works only for univariate functional data.
#' 
#' @details
#' \bold{Model specification: }
#' 
#' Let \eqn{T_{ij}} denote the time of the \eqn{j}th event for subject \eqn{i}, 
#' and let \eqn{C_i} represent the censoring time. 
#' The observed event time, accounting for right censoring, is \eqn{\widetilde{T}_{ij}=\min(T_{ij}, C_i)}, 
#' and \eqn{\delta_{ij}=I(T_{ij}\leq C_i)} serves as an indicator of whether the \eqn{j}th event for subject \eqn{i} is observed.
#' The hazard function is specified as
#' \deqn{h(t; \boldsymbol{Z}_i, {X}_i(\cdot))=h_{0}(t-t_{i,j-1}) \exp \left(\eta_{ij}\right),}
#' where \eqn{h_0(\cdot)} is the baseline hazard function, 
#' and \eqn{\eta_{ij} = \bm{\alpha}^{\top}\boldsymbol{Z}_i +\int_{t_{i, j-1}}^{t}{X}_{i}(s)\beta(s)ds + v_{ij}}. 
#' Here, \eqn{t_{i, j-1}} is the previous event time with \eqn{t_{i0} = 0}. 
#' \eqn{\bm{\alpha}} is the fixed effect parameter associated with the time-invariant covariates \eqn{\boldsymbol{Z}_i}, 
#' and \eqn{\beta(t)} is a time-varying coefficient that captures the effect of functional predictor \eqn{X_{i}(t)} on the hazard rate of recurrent events.
#' 
#' 
#' @param formula A formula, with the response on the left of a ~ operator being a \code{Recur} object as returned by function \code{Recur} in \strong{reda}, 
#' and scalar covariates on the right.
#' @param sdat A data frame containing subject IDs,  time-to-event outcomes (starting time, end point, censoring time and event status), and scalar covariates
#' @param fdat A data frame containing subject IDs, longitudinal measurements, and the corresponding time points for each measurement.
#' @param para0 A vector of initial values for \eqn{\theta^2} and auto-regressive coefficient \eqn{\rho}. Both default to 0.5.
#' @param nbasis An integer, representing the number of B-spline basis functions 
#'    used for estimation of the mean function and bivariate smoothing of the covariance surface. 
#'    Defaults to \code{10} (cf. \code{fpca.sc} in \strong{refund}).
#' @param pve A numeric value between 0 and 1, the proportion of variance 
#'   explained: used to choose the number of principal components. Defaults to 
#'   \code{0.9} (cf. \code{fpca.sc} in \strong{refund}).
#' @param npc An integer, giving a prespecified value for the number of 
#'   principal components. Defaults to \code{NULL}. If given, this overrides 
#'   \code{pve} (cf. \code{fpca.sc} in \strong{refund}).
#' @param makePD Logical: should positive definiteness be enforced for the 
#'   covariance surface estimate? Defaults to \code{FALSE} (cf. 
#'   \code{fpca.sc} in \strong{refund}).
#' @param cov.weight.type The type of weighting used for the smooth covariance 
#'   estimate. Defaults to \code{"none"}, i.e. no weighting. Alternatively, 
#'   \code{"counts"} (corresponds to \code{fpca.sc} in \strong{refund}) weights the
#'   pointwise estimates of the covariance function by the number of observation
#'   points.
#' @param iter.max Maximum number of iterations for both inner iteration and outer iteration. Defaults to \code{50}.
#' @param eps Tolerance criteria for a possible infinite coefficient value. Defaults to \code{1e-6}.
#'   
#' @returns A funsurv object containing the following components:
#'  \item{beta}{Estimation of coefficients of scalar covariates and FPC scores. 
#' Including estimated values, standard errors, and p-values}
#'  \item{beta_vcov}{Estimated variance-covariance of the estimates of beta}
#'  \item{eAR}{Estimation of variance components (\eqn{\theta^2} and \eqn{\rho})}
#'  \item{eAR_vcov}{Estimated variance of estimates of \eqn{\theta^2} and \eqn{\rho}}
#'  \item{frailties}{Estimated frailty terms (random effects)}
#'  \item{basesurv}{Estimated baseline survival probability}
#'  \item{time}{Time points associated with baseline survival probability}
#'  \item{FPC}{Functional principal components}
#'  
#' @importFrom MASS ginv
#' @importFrom stats model.matrix model.extract model.frame
#' @importFrom Matrix bdiag
#'
#' @seealso \code{\link{Recur}}
#' @seealso \code{\link{PACE}}
#' @export
#'
#' @example inst/examples/ex_AR1_FRAILTY.R
AR1_FRAILTY <- function(formula, 
                        sdat,
                        fdat,
                        para0 = c(0.5,0.5), 
                        nbasis = 10, 
                        pve = 0.90,
                        npc = NULL, 
                        makePD = FALSE, 
                        cov.weight.type = c("none", "counts"), 
                        iter.max = 50, 
                        eps = 1e-6) {  
  if (missing(formula)) stop("A Event formula is required.")
  if (missing(sdat)) stop("A dataset of time-to-event outcomes is required.")
  if (missing(fdat)) stop("A dataset of longitudinal measurements is required.")
  cov.weight.type <- match.arg(cov.weight.type)
  data <- sdat
  Call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$data <- data
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(model.frame)
  mf <- eval(mf, parent.frame())
  mm <- model.matrix(formula, data = mf)
  obj <- model.extract(mf, "response")
  DF <- cbind(obj, mm)
  DF <- DF[,colnames(DF) != "(Intercept)"]
  DF <- data.frame(DF)
  DF$gap_time <- DF$time2 - DF$time1  
  ## AR1_PACE
  fpca_obj <- AR1_PACE(sdat, fdat, nbasis = nbasis, pve = pve, npc = npc,
                       makePD = makePD, cov.weight.type = cov.weight.type)
  scores <- fpca_obj$window_scores_fpca
  if(nrow(DF) != nrow(scores)) stop("Dimension of sdat does not match with the dimension of FPC scores")
  DF <- cbind(DF, scores)
  colnames(DF)[which(colnames(DF)=="event")] <- "status"
  covariates <- colnames(mm)[-1]
  
  p <- length(covariates) + ncol(scores)  # number of predictors in the survival model
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
  beta_vcov <- H2[1:p, 1:p]
  names(beta) <- colnames(H2)[1:p]
  ## SE for the variance components
  block_list <- lapply(ni, function(ni) ar1_cor(ni, rho))
  U <- as.matrix(bdiag(block_list))
  Q1 <- ginv(dll.eta) + theta2 * R %*% U %*% t(R)
  U.inv <- bdiag(lapply(ni, function(ni) dar1_cor.drho(ni, rho))) 
  U.inv <- as.matrix(U.inv)
  dQ1.theta2 <- R %*% U %*% t(R)
  dQ1.rho <- theta2 * R %*% U.inv %*% t(R)
  Q1.inv <- solve(Q1)
  Q2 <- Q1.inv - Q1.inv %*% X_surv %*% solve(t(X_surv) %*% Q1.inv %*% X_surv) %*% t(X_surv) %*% Q1.inv
  b11 <- sum(diag(Q2 %*% dQ1.theta2 %*% Q2 %*% dQ1.theta2))
  b12 <- sum(diag(Q2 %*% dQ1.theta2 %*% Q2 %*% dQ1.rho))
  b22 <- sum(diag(Q2 %*% dQ1.rho %*% Q2 %*% dQ1.rho))
  eAR <- c(theta2, rho)
  names(eAR) <- c("theta2", "rho")
  eAR_vcov <- 2 * solve(matrix(c(b11, b12, b12, b22), ncol = 2))
  out <- list(beta = beta, beta_vcov = beta_vcov, 
              eAR = eAR, eAR_vcov = eAR_vcov, 
              frailties = as.vector(V),
              basesurv = as.vector(basesurv), 
              time = as.vector(sort(DF$time2)), 
              FPC_argvals = fpca_obj$FPC_argvals,
              FPC_X = fpca_obj$FPC_X,
              call = Call)
  class(out) <- "funsurv"
  return(out)
}

#' Estimate the variance component \eqn{\theta^2} and auto-regressive coefficient \eqn{\rho}
#'
#' @param ni A vector representing the number of events for each subject.
#' @param tau 
#' @param rho \eqn{\rho} value at previous iteration
#' @param J J matrix, returned by IJK function
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
  theta20 <- ((1 + rho0^2)*A1 - 2*rho0*A2 - (rho0^2)*A3) / N
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
  ## I
  I <- diag(N)    
  ## J
  J <- diag(0,N)
  if(max(ni)>1){
    maxJ <- diag(0,max(ni))
    diag(maxJ[,-1]) <- diag(maxJ[-1,]) <- 1
    J_block_list <- lapply(ni, function(x) maxJ[1:x,1:x])
    J <- as.matrix(bdiag(J_block_list))
  } 
  ## K
  K_block_list <- lapply(ni, function(x) {
    if(x==1) k <- as.matrix(2)
    else{
      k <- diag(0, x)
      k[1,1] <- k[x, x] <- 1
    } 
    return(k)})
  K <- as.matrix(bdiag(K_block_list))
  list(I = I, J = J, K = K)
}

## #' @importFrom methods getClass
is.funsurv <- function(x) inherits(x, "funsurv")
