#' Plot functional coefficient
#'
#' @param object An AR1_FRAILTY object
#'
#' @returns a plot of functional coefficient
#' @export
plot_beta <- function(object){
   npc <- object$PACE$npc
   beta <- object$ebeta[-(1:(nrow(object$ebeta)-npc)), 1] %*% object$PACE$functions@X
   plot(object$PACE$functions@argvals[[1]], beta, type="l", xlab="Follow-up time", ylab=expression("Estimate of " * beta(t)))
}
