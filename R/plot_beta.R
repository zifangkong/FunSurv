<<<<<<< HEAD
#' Plot functional coefficient
#'
#' @param object An AR1_FRAILTY object
#'
#' @returns a plot of functional coefficient
#' @export
#' @example inst/examples/ex_plot_beta.R
plot_beta <- function(object){
  npc <- object$PACE$npc
  beta <- crossprod(object$PACE$functions@X, object$ebeta[-(1:(nrow(object$ebeta) - npc)), 1])
  plot(object$PACE$functions@argvals[[1]], beta, type="l", xlab="Follow-up time",
       ylab=expression("Estimate of " * beta(t)))
}

=======
>>>>>>> 678d4ad (improve print)
#' Plot method for 'funsurv' objects
#'
#' Plots beta(t) ...
#'
#' @return A ggplot object ... 
#' 
#' @importFrom ggplot2 ggplot aes geom_line xlab ylab
#'
#' @exportS3Method plot funsurv
plot.funsurv <- function(x) {
  if (!is.funsurv(x)) stop("Must be a funsurv object")
  npc <- x$PACE$npc
  beta <- crossprod(x$PACE$functions@X, x$beta[-(1:(length(x$beta) - npc))])
  ggplot(NULL, aes(x = x$PACE$functions@argvals[[1]], y = beta)) +
    geom_line() + xlab("Follow-up time") + ylab(expression("Estimate of " * beta(t)))
}
