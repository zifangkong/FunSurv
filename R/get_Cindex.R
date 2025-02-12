#' Calculate C-index given time
#' 
get_Cindex <- function(tt, surv_data, one_simulation, betaType = c("curve", "constant", "sincos")) {
  uni_PACE <- one_simulation$uni.PACE
  betaType <- match.arg(betaType)
  if (betaType == "sincos") {
    surv_data$frailty_fpca_1 <- one_simulation$mm_fpca_1$frailty
    surv_data$frailty_fpca_2 <- one_simulation$mm_fpca_2$frailty
    surv_data$frailty_bspline <- one_simulation$mm_bspline$frailty
  } else{
    surv_data$frailty_fpca_1 <- one_simulation$mm_fpca_pve0.90$frailty
    surv_data$frailty_fpca_2 <- one_simulation$mm_fpca_pve0.95$frailty
    surv_data$frailty_bspline <- one_simulation$mm_bspline$frailty
  }
  ## at time tt, for each subject i, tt, censorTime, numEvent, riskScore,
  argvals_irregular <- uni_PACE[[1]]$mu@argvals[[1]]
  if (betaType == "constant") {
    bspline_df <- 4
    bspline_intercept <- TRUE
  } else{
    bspline_df <- 6
    bspline_intercept <- FALSE
  }
  bspline_basis <- bSpline(
    argvals_irregular,
    df <- bspline_df,
    degree <- 3,
    intercept <- bspline_intercept
  )  
  dat.risk <- c()
  for (idindex in 1:length(unique(surv_data$id))) {
    id <- unique(surv_data$id)[idindex]
    temp_df <- surv_data[surv_data$id == id, ]
    temp_df <- temp_df[temp_df$cum_time <= tt,]
    if (nrow(temp_df) > 0) {
      ## number of events before time tt
      temp_df$numEvent <- sum(temp_df$status)
      ## to find the int(x(t)beta(t)) from initial_time to tt
      initial_time <- temp_df$cum_time[nrow(temp_df)]
      ## initial_time <- 0
      initial_idx <- max(which(argvals_irregular <= initial_time))
      end_idx <- max(which(argvals_irregular <= tt))
      obs.points <- initial_idx:end_idx
      if (length(obs.points) == 1) {
        J_fpca_1 <- (uni_PACE[[1]]$functions@X[, obs.points]) %*% t(uni_PACE[[1]]$functions@X[, obs.points]) * diag(uni_PACE[[1]]$w0)[obs.points, obs.points]
        J_fpca_2 <- (uni_PACE[[2]]$functions@X[, obs.points]) %*% t(uni_PACE[[2]]$functions@X[, obs.points]) * diag(uni_PACE[[2]]$w0)[obs.points, obs.points]
        ## J_fpca_3 <- (uni_PACE[[3]]$functions@X[,obs.points]) %*% t(uni_PACE[[3]]$functions@X[,obs.points]) * diag(uni_PACE[[3]]$w0)[obs.points,obs.points]
        J_bspline <- (uni_PACE[[1]]$functions@X[, obs.points]) %*% t(bspline_basis[obs.points, ]) * diag(uni_PACE[[1]]$w0)[obs.points, obs.points]
      } else{
        J_fpca_1 <- (uni_PACE[[1]]$functions@X[, obs.points]) %*% diag(uni_PACE[[1]]$w0)[obs.points, obs.points] %*% t(uni_PACE[[1]]$functions@X[, obs.points])
        J_fpca_2 <- (uni_PACE[[2]]$functions@X[, obs.points]) %*% diag(uni_PACE[[2]]$w0)[obs.points, obs.points] %*% t(uni_PACE[[2]]$functions@X[, obs.points])
        ## J_fpca_3 <- (uni_PACE[[3]]$functions@X[, obs.points]) %*% diag(uni_PACE[[3]]$w0)[obs.points, obs.points] %*% t(uni_PACE[[3]]$functions@X[, obs.points])
        J_bspline <- (uni_PACE[[1]]$functions@X[, obs.points]) %*% diag(uni_PACE[[1]]$w0)[obs.points, obs.points] %*% (bspline_basis[obs.points, ])
      }
      if (betaType == "sincos") {
        temp_df$risk_fpca_1 <- sum((t(uni_PACE[[1]]$scores[surv_data$id[idindex], ])) %*% J_fpca_1 * one_simulation$mm_fpca_1$ebeta[-1, 1]) + temp_df$z1[1] *
          one_simulation$mm_fpca_1$ebeta[1, 1] + temp_df$frailty_fpca_1[nrow(temp_df)]
        temp_df$risk_fpca_2 <- sum((t(uni_PACE[[2]]$scores[surv_data$id[idindex], ])) %*% J_fpca_2 * one_simulation$mm_fpca_2$ebeta[-1, 1]) + temp_df$z1[1] *
          one_simulation$mm_fpca_2$ebeta[1, 1] + temp_df$frailty_fpca_2[nrow(temp_df)]
        ## temp_df$risk_fpca_3 <- sum((t(uni_PACE[[3]]$scores[surv_data$id[idindex],])) %*% J_fpca_3 * one_simulation$mm_fpca_3$ebeta[-1,1]) + temp_df$z1[1]*one_simulation$mm_fpca_3$ebeta[1,1] + temp_df$frailty_fpca_3[nrow(temp_df)]
        temp_df$risk_bspline <- sum(t(uni_PACE[[1]]$scores[surv_data$id[idindex], ]) %*% J_bspline * one_simulation$mm_bspline$ebeta[-1, 1]) + temp_df$z1[1] *
          one_simulation$mm_bspline$ebeta[1, 1] + temp_df$frailty_bspline[nrow(temp_df)]
      } else{
        temp_df$risk_fpca_1 <- sum((t(uni_PACE[[1]]$scores[surv_data$id[idindex], ])) %*% J_fpca_1 * one_simulation$mm_fpca_pve0.90$ebeta[-1, 1]) + temp_df$z1[1] *
          one_simulation$mm_fpca_pve0.90$ebeta[1, 1] + temp_df$frailty_fpca_1[nrow(temp_df)]
        temp_df$risk_fpca_2 <- sum((t(uni_PACE[[2]]$scores[surv_data$id[idindex], ])) %*% J_fpca_2 * one_simulation$mm_fpca_pve0.95$ebeta[-1, 1]) + temp_df$z1[1] *
          one_simulation$mm_fpca_pve0.95$ebeta[1, 1] + temp_df$frailty_fpca_2[nrow(temp_df)]
        ## temp_df$risk_fpca_3 <- sum((t(uni_PACE[[3]]$scores[surv_data$id[idindex],])) %*% J_fpca_3 * one_simulation$mm_fpca_3$ebeta[-1,1]) + temp_df$z1[1]*one_simulation$mm_fpca_3$ebeta[1,1] + temp_df$frailty_fpca_3[nrow(temp_df)]
        temp_df$risk_bspline <- sum(t(uni_PACE[[1]]$scores[surv_data$id[idindex], ]) %*% J_bspline * one_simulation$mm_bspline$ebeta[-1, 1]) + temp_df$z1[1] *
          one_simulation$mm_bspline$ebeta[1, 1] + temp_df$frailty_bspline[nrow(temp_df)]
      }
      temp_df <- temp_df[1, ]
    }
    dat.risk <- rbind(dat.risk, temp_df)
  }
  nComp <- nConc_fpca_1 <- nConc_fpca_2 <- nConc_fpca_3 <- nConc_bspline <- 0
  for (pair_i in 1:(nrow(dat.risk) - 1)) {
    for (pair_j in (pair_i + 1):nrow(dat.risk)) {
      if ((dat.risk$numEvent[pair_i] > dat.risk$numEvent[pair_j]) &
          (dat.risk$censoring_time[pair_i] >= tt) &
          (dat.risk$censoring_time[pair_j] >= tt)) {
        nComp <- nComp + 1
        nConc_fpca_1 <- nConc_fpca_1 + as.numeric((dat.risk$risk_fpca_1[pair_i] > dat.risk$risk_fpca_1[pair_j]))
        nConc_fpca_2 <- nConc_fpca_2 + as.numeric((dat.risk$risk_fpca_2[pair_i] > dat.risk$risk_fpca_2[pair_j]))
        ## nConc_fpca_3 <- nConc_fpca_3 + as.numeric((dat.risk$risk_fpca_3[pair_i] > dat.risk$risk_fpca_3[pair_j]))
        nConc_bspline <- nConc_bspline + as.numeric((dat.risk$risk_bspline[pair_i] > dat.risk$risk_bspline[pair_j]))
      } else if ((dat.risk$numEvent[pair_i] < dat.risk$numEvent[pair_j]) &
                 (dat.risk$censoring_time[pair_i] >= tt) &
                 (dat.risk$censoring_time[pair_j] >= tt)) {
        nComp <- nComp + 1
        nConc_fpca_1 <- nConc_fpca_1 + as.numeric((dat.risk$risk_fpca_1[pair_i] < dat.risk$risk_fpca_1[pair_j]))
        nConc_fpca_2 <- nConc_fpca_2 + as.numeric((dat.risk$risk_fpca_2[pair_i] < dat.risk$risk_fpca_2[pair_j]))
        ## nConc_fpca_3 <- nConc_fpca_3 + as.numeric((dat.risk$risk_fpca_3[pair_i] < dat.risk$risk_fpca_3[pair_j]))
        nConc_bspline <- nConc_bspline + as.numeric((dat.risk$risk_bspline[pair_i] < dat.risk$risk_bspline[pair_j]))
      }
    }
  }
  c(nConc_fpca_1 / nComp,
    nConc_fpca_2 / nComp,
    nConc_bspline / nComp)
}
