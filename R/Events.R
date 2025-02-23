#' Define Recurrent Event Times
#'
#' @param t_start Entry time
#' @param t_stop Exit time
#' @param id Subject id
#' @param status Event indicator (1 = event, 0 = censored)
#'
#' @return A structured list with event attributes
#' @export
Events <- function(t_start, t_stop, id, status) {
   structure(
      list(t_start = t_start, t_stop = t_stop, id = id, status = status),
      class = "Events"
   )
}
