#' Not a main function? 
#'
#' @noRd
ar1_get_window_scores_mfpca <- function(mFPCA.mod,
                                        current_PACE,
                                        which_urine,
                                        surv_data,
                                        bspline_basis) {
  sigma2 <- current_PACE$sigma2
  argvals_irregular <- current_PACE$mu@argvals[[1]]
  w0 <- current_PACE$w0
  ZtZ_sD.inv <- current_PACE$ZtZ_sD.inv
  Zcur <- current_PACE$Zcur
  npc <- sum(mFPCA.mod$npc)
  mfpca_functions <- t(mFPCA.mod$functions[[which_urine]])
  window_scores_fpca <- matrix(NA, nrow = nrow(surv_data), ncol = npc)
  counter <- 1
  for (i.subj in 1:length(unique(surv_data$id))) {
    temp_df <- surv_data[surv_data$id == unique(surv_data$id)[i.subj],]
    idx <- sapply(temp_df$cum_time, function(x)
      max(which(argvals_irregular <= x)))
    idx <- (c(0, idx))
    for (i in 2:length(idx)) {
      obs.points <- (idx[i - 1] + 1):idx[i]
      if (length(obs.points) == 2 &
          (obs.points[1] > obs.points[2])) {
        obs.points <- obs.points[2]
      }      
      if (sigma2 == 0 & length(obs.points) < npc) {
        stop(
          "Measurement error estimated to be zero and there are fewer observed points than PCs; scores cannot be estimated."
        )
      }
      if (length(obs.points) == 1) {
        J_fpca_basis <- (mfpca_functions[, obs.points]) %*% t(mfpca_functions[, obs.points]) * diag(w0)[obs.points, obs.points]
      } else{
        J_fpca_basis <- (mfpca_functions[, obs.points]) %*% diag(w0)[obs.points, obs.points] %*% t(mfpca_functions[, obs.points])
      }
      ss_fpca <- (t(mFPCA.mod$scores[i.subj, ])) %*% J_fpca_basis
      window_scores_fpca[counter,] <- ss_fpca
      counter <- counter + 1
    }
  }
  colnames(window_scores_fpca) <- paste0("ws", seq(npc))
  return(list(window_scores_fpca = window_scores_fpca))
}
