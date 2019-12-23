#' Count lines in a file
#'
#' \code{count_lines} uses UNIX system command \code{wc} to determine the number
#'   of lines in the given file
#'
#' @param file path to file in question
#'
#' @export

count_lines <- function(file) {
  as.numeric(system(paste('wc -l <', file), intern = T))
}