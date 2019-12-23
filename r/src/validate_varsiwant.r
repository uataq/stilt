#' validate_varsiwant checks for valid 4-letter variable codes
#' @author Ben Fasoli
#'
#' @param varsiwant character vector of 4-letter hysplit variables
#'
#' @export

validate_varsiwant <- function(varsiwant) {

  valid <- c('time', 'indx', 'long', 'lati', 'zagl', 'sigw', 'tlgr', 'zsfc',
             'icdx', 'temp', 'samt', 'foot', 'shtf', 'tcld', 'dmas', 'dens',
             'rhfr', 'sphu', 'solw', 'lcld', 'zloc', 'dswf', 'wout', 'mlht', 
             'rain', 'crai', 'pres')
  invalid <- setdiff(varsiwant, valid)

  if (length(invalid) > 0)
    stop('validate_varsiwant(): varsiwant argument can only include:\n',
         paste(collapse = '\n', valid))
  return(invisible(T))
}
