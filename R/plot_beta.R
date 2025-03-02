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
