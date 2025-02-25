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
   dat <- data.frame(id=id, t_start=t_start, t_stop=t_stop, status=status)
   sorted_dat <- dat[order(dat$id),]
   char_rec <- lapply( unique(sort(id)),
      function(one_id) {
         idx <- which(dat$id==one_id)
         time1 <- sorted_dat[idx, "t_start"]
         time2 <- sorted_dat[idx, "t_stop"]
         is_censored <- sorted_dat[idx, "status"] == 0
         sign <- ifelse(is_censored, "", "+")
         intervals <- paste0("(", sprintf("%.2f", time1), ", ", sprintf("%.2f", time2), sign, "]")
         if (length(intervals) > 4) {
            intervals <- c(intervals[1], intervals[2], "...", intervals[length(intervals)])
         }
         paste(one_id, ": ", paste(intervals, collapse = ", "), sep = "")
      })
   cat(paste(char_rec, collapse = "\n"), "\n")
   structure(
      list(t_start = t_start, t_stop = t_stop, id = id, status = status),
      class = "Events"
   )
}
