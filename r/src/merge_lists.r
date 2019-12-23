#' Merge two lists by named keys
#'
#' \code{merge_lists} joins two lists by key, giving preference to values in 
#'   \code{y} when duplicate keys are found
#' 
#' @param x list
#' @param y list
#'
#' @export

merge_lists <- function(x, y = list()) {
  z <- x
  keys <- names(y)
  for (k in keys) {
    z[[k]] <- y[[k]]
  }
  z
}