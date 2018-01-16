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

write_zierr <- function(sigzierr = NA, tlzierr = NA, horcorzierr = NA,
                        file = 'ZIERR') {

  txt <- c(sigzierr,
           tlzierr,
           horcorzierr)

  if (any(is.na(txt))) {
    if (all(is.na(txt)))
      return()
    stop('write_zierr(): Incorrect specification for ',
         paste(input[is.na(input)], sep = ','))
  }

  write(as.character(txt), file)
  file
}
