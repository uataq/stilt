#' Linearly interpolate NA values
#'
#' \code{na_interp} linearly interpolates NA values found in vector y with
#'   respect to index x (e.g. timestamp).
#'
#' @param y numeric vector in which to fill bracketed NA values
#' @param x vector giving index of y. NULL if y is equally spaced
#'
#' @importFrom stats approx
#' @importFrom utils head tail
#' 
#' @export

na_interp <- function (y, x = NULL) {
  if (is.null(x)) 
    x <- 1:length(y)
  nona <- which(!is.na(y))
  start <- head(nona, 1)
  end <- tail(nona, 1)
  if (length(end - start) < 1) 
    return(y)
  xsub <- x[start:end]
  ysub <- y[start:end]
  idx <- which(is.na(ysub))
  ysub[idx] <- approx(xsub, ysub, xout = xsub[idx], method = "linear")$y
  y[start:end] <- ysub
  return(y)
}