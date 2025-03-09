#' Plot method for 'funsurv' objects
#'
#' @param x A funsurv object
#' @param which A character string specifying what to be plotted. 
#' Use \code{which = "beta"} to plot the estimated \eqn{\beta(t)}.
#' Use \code{which = "fpc"} to plot the functional principal components associated with the the longitudinal measurements.
#' Use \code{which = "basesurv"} to plot the baseline survival probabilities.
#' 
#' @return A ggplot object ... 
#' 
#' @importFrom ggplot2 ggplot geom_ribbon aes geom_line xlab ylab labs
#' @importFrom dplyr %>% 
#' @importFrom tidyr pivot_longer
#' @exportS3Method plot funsurv
#' 
#' @example inst/examples/ex_plot.R
plot.funsurv <- function(x, which = c("beta", "fpc", "basesurv")) {
  if (!is.funsurv(x)) stop("Must be a funsurv object")
  which <- match.arg(which) 
  if(which == "beta"){
    npc <- nrow(x$FPC@X)
    beta <- crossprod(x$FPC@X, x$beta[-(1:(length(x$beta) - npc))])
    moe <- 1.96*sqrt(diag(t(x$FPC@X) %*% x$beta_vcov[-(1:(length(x$beta) - npc)) , -(1:(length(x$beta) - npc))] %*% x$FPC@X))
    lower <- beta - moe
    upper <- beta + moe
    ggplot(NULL, aes(x = x$FPC@argvals[[1]], y = beta)) +
      geom_line() + xlab("Follow-up time (year)") + ylab(expression("Estimate of " * beta(t))) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "blue")
  }else if(which == "fpc"){
    df <- data.frame(time = x$FPC@argvals[[1]])
    for (i in 1:nrow(x$FPC@X)) {
      df[[paste0("FPC", i)]] <- x$FPC@X[i, ]
    }
    df_long <- df %>%
      pivot_longer(cols = -time, names_to = "FPC", values_to = "value")
    ggplot(df_long, aes(x = time, y = value, color = FPC)) +
      geom_line() +
      labs( x = "Follow-up time (year)", y= "", color = "FPC")
  }else if(which == "basesurv"){
    ggplot(NULL, aes(x = x$time, y = x$basesurv)) +
      geom_line() + xlab("Follow-up time (year)") + ylab("Baseline survival probability")
  }else{
    stop("Invalid input: 'which' must be 'beta', 'fpc', or basesurv'.")
  }
}
