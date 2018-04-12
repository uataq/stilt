#' validate_projection checks for valid PROJ.4 crs strings
#' @author Ben Fasoli
#'
#' @param projection character string containing PROJ.4 compatible CRS string
#'
#' @export

validate_projection <- function(projection) {
  
  is_longlat <- grepl('+proj=longlat', projection, fixed = T)
  if (is_longlat)
    return(invisible(T))
  
  if (!'proj4' %in% installed.packages())
    warning('validate_projection(): Manual installation of PROJ4 and the proj4',
            ' R package may be required for alternative map projections.')
  
  load_libs('proj4')
  if (class(try(crs(projection))) != 'CRS')
    stop('validate_projection(): Invalid PROJ.4 projection string. See ',
         'http://www.remotesensing.org/geotiff/proj_list/')
  
  return(invisible(T))
}
