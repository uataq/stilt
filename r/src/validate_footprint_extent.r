#' validate_footprint_extent checks that user has set the footprint domain
#' @author Ben Fasoli
#'
#' @param xmn minimum x (longitude) position
#' @param xmx maximum x (longitude) position
#' @param ymn minimum y (latitude) position
#' @param ymx maximum y (latitude) position
#'
#' @export

validate_footprint_extent <- function(xmn, xmx, ymn, ymx) {
  
  if (any(is.na(c(xmn, xmx, ymn, ymx)))) {
      stop('validate_footprint_extent(): Must specify xmn, xmx, ymn, ymx ',
           'for your domain')
  }
  
  return(invisible(T))
}
