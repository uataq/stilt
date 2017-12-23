#' write_zierr writes a ZIERR file to control the model behavior
#' @author Ben Fasoli
#'
#' Controls wind error covariance options for vertical transport error.
#'
#' @param sigzierr standard deviation of mixed layer height errors [%]
#' @param tlzierr standard deviation of mixed layer height timescale [min]
#' @param horcorzierr horizontal correlation lengthscale of mixed layer height
#'   errors [km]
#' @param file path and name for output file
#'
#' @export

write_zierr <- function(sigzierr = NULL, tlzierr = NULL, horcorzierr = NULL,
                        file = 'ZIERR') {

  txt <- c(sigzierr,
           tlzierr,
           horcorzierr)

  if (any(is.null(txt))) {
    if (all(is.null(txt)))
      return()
    stop('write_zierr(): Incorrect specification for ',
         paste(input[is.null(input)], sep = ','))
  }

  write(as.character(txt), file)
  file
}
