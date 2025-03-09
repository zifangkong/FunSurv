#' Plot method for 'funsurv' objects
#'
#' @param x A funsurv object
#' @param what A character string specifying what to be plotted. 
#' Use \code{what = "beta"} to plot the estimated \eqn{\beta(t)}.
#' Use \code{what = "fpc"} to plot the functional principal components associated with the the longitudinal measurements.
#' Use \code{what = "basesurv"} to plot the baseline survival probabilities.
#' @param ... additional graphical parameters to be passed to methods.
#' 
#' @return A ggplot object ... 
#' 
#' @importFrom ggplot2 ggplot geom_ribbon aes geom_line geom_step xlab ylab labs
#' @exportS3Method plot funsurv
#' 
#' @example inst/examples/ex_plot.R
plot.funsurv <- function(x, what = c("beta", "fpc", "basesurv"), ...) {
  if (!is.funsurv(x)) stop("Must be a funsurv object")
  what <- match.arg(what)
  npc <- nrow(x$FPC_X)
  if(what == "beta"){
    beta <- crossprod(x$FPC_X, x$beta[-(1:(length(x$beta) - npc))])
    moe <- 1.96*sqrt(diag(t(x$FPC_X) %*%
                          x$beta_vcov[-(1:(length(x$beta) - npc)) , -(1:(length(x$beta) - npc))] %*%
                          x$FPC_X))
    lower <- beta - moe
    upper <- beta + moe
    ggplot(NULL, aes(x = x$FPC_argvals[[1]], y = beta)) +
      geom_line() + xlab("Follow-up time (year)") + ylab(expression("Estimate of " * beta(t))) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "blue")
  }
  if(what == "fpc"){
    df <- data.frame(time = rep(x$FPC_argvals, each = npc),
                     FPC = c(x$FPC_X), group = 1:npc)
    ggplot(NULL, aes(x = df$time, y = df$FPC, color = factor(df$group))) +
      geom_line() +
      labs( x = "Follow-up time (year)", y= "", color = "FPC")
  } 
  if(what == "basesurv"){
    ggplot(NULL, aes(x = x$time, y = x$basesurv)) +
      geom_step() + xlab("Follow-up time (year)") + ylab("Baseline survival probability")
  } 
}
