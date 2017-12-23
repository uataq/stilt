#' write_winderr writes a WINDERR file to control the model behavior
#' @author Ben Fasoli
#'
#' Controls wind error covariance options for horizontal transport error.
#'
#' @param siguverr standard deviation of horizontal wind errors [m/s]
#' @param tluverr standard deviation of horizontal wind error timescale [min]
#' @param zcoruverr vertical correlation lengthscale [m]
#' @param horcoruverr horizontal correlation lengthscale [km]
#' @param file path and name for output file
#'
#' @export

write_winderr <- function(siguverr = NULL, tluverr = NULL, zcoruverr = NULL,
                          horcoruverr = NULL, file = 'WINDERR') {

  txt <- c(siguverr,
           tluverr,
           zcoruverr,
           horcoruverr)

  if (any(is.null(txt))) {
    if (all(is.null(txt)))
      return()
    stop('write_winderr(): Incorrect specification for ',
         paste(input[is.null(input)], sep = ','))
  }

  write(as.character(txt), file)
  file
}
