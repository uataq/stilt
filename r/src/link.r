#' Symbolically link a file to another location
#' @author Ben Fasoli
#'
#' @export

link <- function(file, destination) {
  if (file.exists(destination))
    system(paste('unlink', destination))
  system(paste('ln -s', file, destination))
}