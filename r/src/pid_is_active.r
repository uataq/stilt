#' pid_is_active checks pid against ps -A table
#' @author Ben Fasoli
#'
#' @param pid process identifier to match against ps -A pid column
#'
#' @export

pid_is_active <- function(pid) {

  pid <- as.character(pid)

  cmd  <- paste('ps -A -o pid')
  any(grepl(pid, system(cmd, intern = T), fixed = T))
}
