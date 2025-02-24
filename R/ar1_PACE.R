#' Calculate event-specific FPC scores
#' 
#' This function works only for univariate functional data.
#' This function calculates the FPC scores associated with each event gap time by 
#' 1) calculating the overall FPC scores
#' 2) estimating the event-specific FPC scores 
#' \code{fpca.sc} in package \strong{refund}.
#' \code{PACE} in package \strong{MFPCA}.
#'   
#' @param funDataObject An object of class \code{\link[funData]{funData}} or 
#'   \code{\link[funData]{irregFunData}} containing the functional data 
#'   observed, for which the functional principal component analysis is 
#'   calculated. If the data is sampled irregularly (i.e. of class 
#'   \code{\link[funData]{irregFunData}}), \code{funDataObject} is transformed 
#'   to a \code{\link[funData]{funData}} object first.
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
#' @return \item{scores}{An matrix of estimated event-specific scores for the 
#'   observations in \code{funDataObject}. Each row corresponds to the scores of
#'   one observation.}
#' @return \item{uni.PACE}{A univariate functional principal components analysis}
#' @seealso \code{\link[funData]{funData}}, \code{\link{fpcaBasis}}, \code{\link{univDecomp}}
#' @import MFPCA
#' @export AR1_PACE
#'   

AR1_PACE <- function(fun_data, surv_data, nbasis = 10, pve = 0.90,
                     npc = NULL, makePD = FALSE, cov.weight.type = "none"){
  if(length(unique(surv_data$id)) > length(unique(fun_data$id))){
    removed_subj = unique(surv_data$id)[! unique(surv_data$id) %in% unique(fun_data$id)]
    warning("Subjects" , paste0(removed_subj, collapse = ","), "were removed from analysis.")
  }
  surv_data <- surv_data[which(surv_data$id %in% unique(fun_data$id)), ]
  
  ## transform fun_data to a functional data object
  data_split <- split(fun_data, fun_data$id)
  argvals <- lapply(data_split, function(df) df$time)
  xList <- lapply(data_split, function(df) df$x)
  x_FunObject <- irregFunData(argvals = argvals, X = xList)
  
  ## apply functional principal component analyssis conditional expectation to the functional object
  uni.PACE <- MFPCA::.PACE(x_FunObject, nbasis=nbasis, pve=pve)
  
  sigma2 <- uni.PACE$sigma2
  argvals_irregular <- uni.PACE$mu@argvals[[1]]
  w0 <- funData::.intWeights(argvals_irregular, method = "trapezoidal")
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
    
    # Extract relevant portions of X matrix
    X_obs <- X_matrix[, obs_points, drop=FALSE]
    W_obs <- diag(w0)[obs_points, obs_points, drop=FALSE]
    
    # Calculate J_fpca_basis using tcrossprod
    J_fpca_basis <- X_obs %*% tcrossprod(W_obs, X_obs)
    
    # Return the scores
    tcrossprod(uni.PACE$scores[id,], J_fpca_basis)
  }
  
  window_scores_fpca = matrix(NA, nrow=nrow(surv_data), ncol=npc)
  counter <- 1
  for (id in seq_along(unique(surv_data$id))) {
    temp_df <- surv_data[surv_data$id == unique(surv_data$id)[id], ]
    
    # Vectorize the index calculation
    idx <- c(0, vapply(temp_df$t_stop, function(x) {
      max(which(argvals_irregular <= x))
    }, FUN.VALUE = numeric(1)))
    
    # get observation points for each event time window
    obs_points_list <- lapply(2:length(idx), function(i) { (idx[i-1] + 1):idx[i] })
    
    # calculate scores for each event time
    scores <- lapply(obs_points_list, calculate_event_score)
    
    # assign scores to result matrix
    window_scores_fpca[counter:(counter + length(scores) - 1), ] <- do.call(rbind, scores)
    counter <- counter + length(scores)
  }
  colnames(window_scores_fpca) <- paste0("score", seq(npc))
  return(list(window_scores_fpca=window_scores_fpca, uni.PACE=uni.PACE))
}




