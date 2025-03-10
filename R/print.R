#' @exportS3Method print funsurv
print.funsurv <- function(x, ...) {
  if (!is.funsurv(x)) stop("Must be a funsurv object")
  cat("Call: \n")
  dput(x$call) 
  cat("\n Coefficients: \n")
  print(x$beta)
  cat("\n Variance component and auto-regressive coefficient: \n")
  print(x$eAR)
  cat("\n")
}


#' @importFrom stats pnorm
#' @exportS3Method summary funsurv
summary.funsurv <- function(object, ...) {
  if (!is.funsurv(object)) stop("Must be a funsurv object")
  se.beta <- sqrt(abs(diag(object$beta_vcov)))
  tabbeta <- cbind(object$beta, se.beta, 2 * (1 - pnorm(abs(object$beta / se.beta))))
  se.eAR <- sqrt(diag(object$eAR_vcov))
  tabAR <- cbind(object$eAR, se.eAR, 2 * (1 - pnorm(abs(object$eAR / se.eAR))))
  colnames(tabbeta) <- colnames(tabAR) <- c("Estimate", "SE", "p-value")
  out <- list(call = object$call, tabbeta = tabbeta, tabAR = tabAR)
  class(out) <- "summary.funsurv"
  return(out)
}

#' @importFrom stats printCoefmat
#' @exportS3Method print summary.funsurv
print.summary.funsurv <- function(x, ...) {
  cat("Call: \n")
  dput(x$call) 
  cat("\n Coefficients: \n")
  printCoefmat(x$tabbeta, P.values = TRUE, has.Pvalue = TRUE, signif.legend = FALSE)
  cat("\n Variance component and auto-regressive coefficient: \n")
    printCoefmat(x$tabAR, P.values = TRUE, has.Pvalue = TRUE, signif.legend = TRUE)
  cat("\n")
}
 
