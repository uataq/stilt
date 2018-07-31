#' write_zicontrol writes a ZICONTROL file to control the model behavior
#' @author Ben Fasoli
#'
#' Controls wind error covariance options for vertical transport error.
#'
#' @param ziscale vector with which to scale the mixed-layer height, with each 
#'   element specifying a scaling factor for each simulation hour (ziscale can
#'   be of length that is smaller than abs(nhrs)
#' @param file path and name for output file
#'
#' @export

write_zicontrol <- function(ziscale, file = 'ZICONTROL') {
  txt <- c(length(ziscale), ziscale)
  write(as.character(txt), file)
  file
}
