#' Doesn't look like a main function?!
#'
#' @noRd
#' 
ar1_get_window_scores <- function(uni.PACE,
                                  which.pace,
                                  surv_data,
                                  bspline_basis) {
  current_PACE <- uni.PACE[[which.pace]]
  sigma2 <- current_PACE$sigma2
  argvals_irregular <- current_PACE$mu@argvals[[1]]
  w0 <- current_PACE$w0
  ## Z <- current_PACE$Z
  ## D.inv <- current_PACE$D.inv
  ## mu <- current_PACE$mu@X
  ## Y.tilde.fitted <- current_PACE$fit@X - matrix(mu, I.pred, D, byrow <- TRUE)
  ZtZ_sD.inv <- current_PACE$ZtZ_sD.inv
  Zcur <- current_PACE$Zcur
  npc <- current_PACE$npc
  window_scores_fpca <- matrix(NA, nrow = nrow(surv_data), ncol = npc)
  window_scores_bspline <- matrix(NA,
                                  nrow = nrow(surv_data),
                                  ncol = ncol(bspline_basis))  
  counter <- 1
  for (i.subj in 1:length(unique(surv_data$id))) {
    temp_df <- surv_data[surv_data$id == unique(surv_data$id)[i.subj],]
    idx <- sapply(temp_df$cum_time, function(x) max(which(argvals_irregular <= x)))
    idx <- c(0, idx)
    for (i in 2:length(idx)) {
      obs.points <- (idx[i - 1] + 1):idx[i]
      if (length(obs.points) == 2 & (obs.points[1] > obs.points[2])) {
        obs.points <- obs.points[2]
      }
      # obs.points = (1):idx[i]
      if (sigma2 == 0 & length(obs.points) < npc) {
        stop(
          "Measurement error estimated to be zero and there are fewer observed points than PCs; scores cannot be estimated."
        )
      }
      if (length(obs.points) == 1) {
        J_fpca_basis <- (current_PACE$functions@X[, obs.points]) %*% t(current_PACE$functions@X[, obs.points]) * diag(w0)[obs.points, obs.points]
        J_bspline_basis <- (current_PACE$functions@X[, obs.points]) %*% t(bspline_basis[obs.points, ]) * diag(w0)[obs.points, obs.points]
      } else if (current_PACE$npc == 1) {
        J_fpca_basis <- t(as.vector(current_PACE$functions@X[, obs.points])) %*% diag(w0)[obs.points, obs.points] %*% (as.vector(current_PACE$functions@X[, obs.points]))
        J_bspline_basis <- (current_PACE$functions@X[, obs.points]) %*% diag(w0)[obs.points, obs.points] %*% (bspline_basis[obs.points, ])
      }
      else{
        J_fpca_basis <- (current_PACE$functions@X[, obs.points]) %*% diag(w0)[obs.points, obs.points] %*% t(current_PACE$functions@X[, obs.points])
        J_bspline_basis <- (current_PACE$functions@X[, obs.points]) %*% diag(w0)[obs.points, obs.points] %*% (bspline_basis[obs.points, ])
      }
      ss_fpca <- (t(current_PACE$scores[i.subj, ])) %*% J_fpca_basis
      ss_bspline <- (t(current_PACE$scores[i.subj, ])) %*% J_bspline_basis
      window_scores_fpca[counter,] <- ss_fpca
      window_scores_bspline[counter,] <- ss_bspline
      counter <- counter + 1
    }
  }
  colnames(window_scores_fpca) <- paste0("ws", seq(npc))
  colnames(window_scores_bspline) <- paste0("ws", seq(ncol(bspline_basis)))
  return(
    list(
      window_scores_fpca = window_scores_fpca,
      window_scores_bspline = window_scores_bspline
    )
  )
}
