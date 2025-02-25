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
   plot(object$PACE$functions@argvals[[1]], beta, type="l", xlab="Follow-up time", ylab=expression("Estimate of " * beta(t)))
}
