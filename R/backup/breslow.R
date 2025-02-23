#' Function to ...
#' @noRd
breslow <- function(Constrain, Constrain0, largest, time, M1, indi, eta) {
  survbase <- exp(-M1 %*% (indi / (t(M1) %*% (exp(eta)))))
  list(survbase = survbase)
}
