#' Find nearest neighbors
#'
#' \code{find_neighbor} finds the closest value index in y for each x,
#'   using findInterval. Returns the index to be applied to y that most
#'   closely matches x.
#'
#' @param x desired values
#' @param y table to search values
#'
#' @export

find_neighbor <- function(x, y){
  x <- as.numeric(x)
  y <- as.numeric(y)
  small <- findInterval(x, y, all.inside=TRUE)
  large <- small + 1
  is.larger <- 2*x > y[small] + y[large]
  return(small + is.larger)
}