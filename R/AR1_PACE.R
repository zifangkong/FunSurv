#' Calculate event-specific FPC scores
#' 
#' This function works only for univariate functional data.
#' This function calculates the FPC scores associated with each event gap time by 
#' 1) calculating the overall FPC scores
#' 2) estimating the event-specific FPC scores 
#' \code{fpca.sc} in package \strong{refund}.
#' \code{PACE} in package \strong{MFPCA}.
#'   
#' @param sdat A data frame containing subject IDs,  time-to-event outcomes (starting time, end point, censoring time and event status), and scalar coefficients
#' @param fdat A data frame containing subject IDs, longitudinal measurements, and the corresponding time points for each measurement.
#' @param nbasis An integer, representing the number of  B-spline basis 
#'   functions used for estimation of the mean function and bivariate smoothing 
#'   of the covariance surface. Defaults to \code{10} (cf. 
#'   \code{fpca.sc} in \strong{refund}).
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
#'   
#' @return A list containing:
#'         \item{scores}{An matrix of estimated event-specific scores for the 
#'   observations in \code{funDataObject}. Each row corresponds to the scores of
#'   one observation.}
#'         \item{FPC}{Functional principal components}
#' @seealso \code{\link[funData]{funData}}, \code{\link{fpcaBasis}}, \code{\link{univDecomp}}
#'
#' @importFrom MFPCA PACE
#' @importFrom funData .intWeights irregFunData
#'   
#' @noRd
AR1_PACE <- function(sdat, fdat, nbasis = 10, pve = 0.90,
                     npc = NULL, makePD = FALSE,
                     cov.weight.type = c("none", "counts")) {
  cov.weight.type <- match.arg(cov.weight.type)
  if(length(unique(sdat$id)) > length(unique(fdat$id))){
    removed_subj = unique(sdat$id)[! unique(sdat$id) %in% unique(fdat$id)]
    warning("Subjects" , paste0(removed_subj, collapse = ","), "were removed from analysis.")
  }
  sdat <- sdat[which(sdat$id %in% unique(fdat$id)), ]  
  ## transform fdat to a functional data object
  data_split <- split(fdat, fdat$id)
  argvals <- lapply(data_split, function(df) df$time)
  xList <- lapply(data_split, function(df) df$x)
  x_FunObject <- irregFunData(argvals = argvals, X = xList)  
  ## apply functional principal component analysis conditional expectation to the functional object
  uni.PACE <- PACE(x_FunObject,
                   nbasis = nbasis, pve=pve, npc = npc,
                   makePD = makePD, cov.weight.type = cov.weight.type)  
  sigma2 <- uni.PACE$sigma2
  argvals_irregular <- uni.PACE$mu@argvals[[1]]
  w0 <- .intWeights(argvals_irregular, method = "trapezoidal")
  npc <- uni.PACE$npc  
  X_matrix <- uni.PACE$functions@X  
  # function to calculate scores for one event gap time
  calculate_event_score <- function(obs_points) {
    if (length(obs_points) == 2 && obs_points[1] > obs_points[2]) {
      obs_points <- obs_points[2]
    }    
    if (sigma2 == 0 && length(obs_points) < npc) {
      stop("Measurement error estimated to be zero and there are fewer observed points than PCs; scores cannot be estimated.")
    }    
    ## Extract relevant portions of X matrix
    X_obs <- X_matrix[, obs_points, drop=FALSE]
    W_obs <- diag(w0)[obs_points, obs_points, drop=FALSE]    
    ## Calculate J_fpca_basis using tcrossprod
    J_fpca_basis <- X_obs %*% tcrossprod(W_obs, X_obs)
    ## Return the scores
    tcrossprod(uni.PACE$scores[id,], J_fpca_basis)
  }  
  window_scores_fpca <- matrix(NA, nrow = nrow(sdat), ncol = npc)
  counter <- 1
  for (id in seq_along(unique(sdat$id))) {
    temp_df <- sdat[sdat$id == unique(sdat$id)[id], ]   
    idx <- c(0, vapply(temp_df$t_stop, function(x) {
      max(which(argvals_irregular <= x))
    }, FUN.VALUE = numeric(1)))    
    ## get observation points for each event time window
    obs_points_list <- lapply(2:length(idx), function(i) { (idx[i-1] + 1):idx[i] })
    ## calculate scores for each event time
    scores <- lapply(obs_points_list, calculate_event_score)
    ## assign scores to result matrix
    window_scores_fpca[counter:(counter + length(scores) - 1), ] <- do.call(rbind, scores)
    counter <- counter + length(scores)
  }
  colnames(window_scores_fpca) <- paste0("score", seq(npc))
  return(list(window_scores_fpca = window_scores_fpca,
              FPC_argvals = uni.PACE$functions@argvals[[1]],
              FPC_X = uni.PACE$functions@X))
}
