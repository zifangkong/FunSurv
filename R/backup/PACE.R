## the PACE function is adopted from the MFPCA package

#' Calculate univariate functional PCA
#'
#' This function is a slightly adapted version of the
#' \code{fpca.sc} function in the \strong{refund} package for
#' calculating univariate functional principal components based on a smoothed
#' covariance function. The smoothing basis functions are penalized splines.
#'
#' @param X A vector of xValues.
#' @param Y A matrix of observed functions (by row).
#' @param Y.pred A matrix of functions (by row) to be approximated using the
#'   functional principal components. Defaults to \code{NULL}, i.e. the
#'   prediction is made for the functions in \code{Y}.
#' @param nbasis An integer, giving the number of B-spline basis to use.
#'   Defaults to \code{10}.
#' @param pve A value between 0 and 1, giving the percentage of variance
#'   explained in the data by the functional principal components. This value is
#'   used to choose the number of principal components. Defaults to \code{0.99}
#' @param npc The number of principal components to be estimated. Defaults to
#'   \code{NULL}. If given, this overrides \code{pve}.
#' @param makePD Logical, should positive definiteness be enforced for the
#'   covariance estimate? Defaults to \code{FALSE}.
#' @param cov.weight.type The type of weighting used for the smooth covariance
#'   estimate. Defaults to \code{"none"}, i.e. no weighting. Alternatively,
#'   \code{"counts"} (corresponds to \code{fpca.sc} in \strong{refund}) weights the pointwise estimates of the covariance function
#'   by the number of observation points.
#'
#' @return \item{fit}{The approximation of \code{Y.pred} (if \code{NULL}, the
#'   approximation of \code{Y}) based on the functional principal components.}
#'   \item{scores}{A matrix containing the estimated scores (observations by
#'   row).} \item{mu}{The estimated mean function.} \item{efunctions}{A matrix
#'   containing the estimated eigenfunctions (by row).} \item{evalues}{The
#'   estimated eigenvalues.} \item{npc}{The number of principal comopnents that
#'   were calculated.} \item{sigma2}{The estimated variance of the measurement
#'   error.}  \item{estVar}{The estimated smooth variance function of the data.}
#'
#' @seealso \code{\link{PACE}}
#'
#' @references Di, C., Crainiceanu, C., Caffo, B., and Punjabi, N. (2009).
#'   Multilevel functional principal component analysis. Annals of Applied
#'   Statistics, 3, 458--488. Yao, F., Mueller, H.-G., and Wang, J.-L. (2005).
#'   Functional data analysis for sparse longitudinal data. Journal of the
#'   American Statistical Association, 100, 577--590.
#'
#' @importFrom mgcv gam predict.gam s te
#'
#' @keywords internal
.PACE <- function(X,
                  Y,
                  Y.pred = NULL,
                  nbasis = 10,
                  pve = 0.99,
                  npc = NULL,
                  makePD = FALSE,
                  cov.weight.type = "none") {
    ## X is time points
    ## Y is the functional data
    if (is.null(Y.pred)) Y.pred <- Y
    D <- NCOL(Y)
    if (D != length(X))
      # check if number of observation points in X & Y are identical
      stop("different number of (potential) observation points differs in X and Y!")
    I <- NROW(Y)
    I.pred <- NROW(Y.pred)
    d.vec <- rep(X, each = I) # use given X-values for estimation of mu
    gam0 <- mgcv::gam(as.vector(Y) ~ s(d.vec, k = nbasis))
    mu <- mgcv::predict.gam(gam0, newdata = data.frame(d.vec = X))
    Y.tilde <- Y - matrix(mu, I, D, byrow = TRUE)
    cov.sum <- cov.count <- cov.mean <- matrix(0, D, D)
    for (i in seq_len(I)) {
      obs.points <- which(!is.na(Y[i,]))
      cov.count[obs.points, obs.points] <- cov.count[obs.points, obs.points] + 1
      cov.sum[obs.points, obs.points] <- cov.sum[obs.points, obs.points] + tcrossprod(Y.tilde[i, obs.points])
    }
    G.0 <- ifelse(cov.count == 0, NA, cov.sum / cov.count)
    G0 <- G.0
    diag.G0 <- diag(G.0)
    diag(G.0) <- NA
    row.vec <- rep(X, each = D) # use given X-values
    col.vec <- rep(X, D) # use given X-values
    cov.weights <- switch(cov.weight.type,
                          none = rep(1, D ^ 2),
                          counts = as.vector(cov.count),
                          stop(
                            "cov.weight.type ",
                            cov.weight.type,
                            " unknown in smooth covariance estimation"
                          ))    
    npc.0 <- matrix(predict.gam(gam(as.vector(G.0) ~ te(row.vec, col.vec, k = nbasis), weights = cov.weights),
                                newdata = data.frame(row.vec = row.vec, col.vec = col.vec)
    ), D, D)
    npc.0 = (npc.0 + t(npc.0)) / 2
    # no extra-option (useSymm) as in fpca.sc-method
    if (makePD) {
      # see fpca.sc
      npc.0 <- {
        tmp <- Matrix::nearPD(
          npc.0,
          corr = FALSE,
          keepDiag = FALSE,
          do2eigen = TRUE,
          trace = options()$verbose
        )
        as.matrix(tmp$mat)
      }
    }
    ## numerical integration for calculation of eigenvalues (see Ramsay & Silverman, Chapter 8)
    w <- funData::.intWeights(X, method = "trapezoidal")
    Wsqrt <- diag(sqrt(w))
    Winvsqrt <- diag(1 / (sqrt(w)))
    V <- Wsqrt %*% npc.0 %*% Wsqrt
    evalues <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
    evalues <- replace(evalues, which(evalues <= 0), 0)
    npc <- ifelse(is.null(npc), min(which(cumsum(evalues) / sum(evalues) > pve)), npc)
    efunctions <- matrix(Winvsqrt %*% eigen(V, symmetric = TRUE)$vectors[, seq(len = npc)],
                         nrow = D, ncol = npc)
    evalues <- eigen(V, symmetric = TRUE, only.values = TRUE)$values[seq_len(npc)]  # use correct matrix for eigenvalue problem
    cov.hat <- efunctions %*% tcrossprod(diag(evalues, nrow = npc, ncol = npc), efunctions)
    ## numerical integration for estimation of sigma2
    T.len <- X[D] - X[1] # total interval length
    T1.min <- min(which(X >= X[1] + 0.25 * T.len)) # left bound of narrower interval T1
    T1.max <- max(which(X <= X[D] - 0.25 * T.len)) # right bound of narrower interval T1
    DIAG <- (diag.G0 - diag(cov.hat))[T1.min:T1.max] # function values
                                        # weights
    w <- funData::.intWeights(X[T1.min:T1.max], method = "trapezoidal")
    sigma2 <- max(1 / (X[T1.max] - X[T1.min]) * sum(DIAG * w, na.rm = TRUE), 0) #max(1/T.len * sum(DIAG*w), 0)
    D.inv <- diag(1 / evalues, nrow = npc, ncol = npc)
    Z <- efunctions
    Y.tilde <- Y.pred - matrix(mu, I.pred, D, byrow = TRUE)
    fit <- matrix(0, nrow = I.pred, ncol = D)
    scores <- matrix(NA, nrow = I.pred, ncol = npc)
    ## no calculation of confidence bands, no variance matrix
    for (i.subj in seq_len(I.pred)) {
      obs.points <- which(!is.na(Y.pred[i.subj,]))
      if (sigma2 == 0 & length(obs.points) < npc) {
        stop(
          "Measurement error estimated to be zero and there are fewer observed points than PCs; scores cannot be estimated."
        )
      }
      Zcur <- matrix(Z[obs.points,], nrow = length(obs.points), ncol = dim(Z)[2])
      ZtZ_sD.inv <- solve(crossprod(Zcur) + sigma2 * D.inv)
      scores[i.subj,] <- ZtZ_sD.inv %*% crossprod(Zcur, Y.tilde[i.subj, obs.points])
      fit[i.subj,] <- t(as.matrix(mu)) + tcrossprod(scores[i.subj,], efunctions)
    }
    ret.objects <- c("fit",
                     "scores",
                     "mu",
                     "efunctions",
                     "evalues",
                     "npc",
                     "sigma2") # add sigma2 to output
    ret <- lapply(seq_len(length(ret.objects)), function(u) get(ret.objects[u]))
    names(ret) <- ret.objects
    ret$estVar <- diag(cov.hat)
    ret$G0 <- G0    
    return(ret)
  }


#' A main function?
#'
#' 
PACE <- function(funDataObject,
                 predData = NULL,
                 nbasis = 10,
                 pve = 0.99,
                 npc = NULL,
                 makePD = FALSE,
                 cov.weight.type = "none") {
  ## check inputs
  if (!class(funDataObject) %in% c("funData", "irregFunData"))
      stop("Parameter 'funDataObject' must be a funData or irregFunData object.")
    if (dimSupp(funDataObject) != 1)
      stop("PACE: Implemented only for funData objects with one-dimensional support.")
    if (methods::is(funDataObject, "irregFunData"))
      # for irregular functional data, use funData representation
      funDataObject <- as.funData(funDataObject)
  if (is.null(predData))
      Y.pred <- NULL # use only funDataObject
    else
    {
      if (!isTRUE(all.equal(funDataObject@argvals, predData@argvals)))
        stop("PACE: funDataObject and predData must be defined on the same domains!")
      Y.pred <- predData@X
    }
    if (!all(is.numeric(nbasis), length(nbasis) == 1, nbasis > 0))
      stop("Parameter 'nbasis' must be passed as a number > 0.")    
    if (!all(is.numeric(pve), length(pve) == 1, 0 <= pve, pve <= 1))
      stop("Parameter 'pve' must be passed as a number between 0 and 1.")    
    if (!is.null(npc) &
        !all(is.numeric(npc), length(npc) == 1, npc > 0))
      stop("Parameter 'npc' must be either NULL or passed as a number > 0.")    
    if (!is.logical(makePD))
      stop("Parameter 'makePD' must be passed as a logical.")    
    if (!is.character(cov.weight.type))
      stop("Parameter 'cov.weight.type' must be passed as a character.")
    res <-
      .PACE(
        X = funDataObject@argvals[[1]],
        funDataObject@X,
        Y.pred = Y.pred,
        nbasis = nbasis,
        pve = pve,
        npc = npc,
        makePD = makePD,
        cov.weight.type = cov.weight.type
      )
  return(
      list(
        mu = funData(funDataObject@argvals, matrix(res$mu, nrow = 1)),
        values = res$evalues,
        functions = funData(funDataObject@argvals, t(res$efunctions)),
        scores = res$scores,
        fit = funData(funDataObject@argvals, res$fit),
        npc = res$npc,
        sigma2 = res$sigma2,
        G0 = res$G0,
        estVar = funData(funDataObject@argvals, matrix(res$estVar, nrow = 1))
      )
    )
  }
