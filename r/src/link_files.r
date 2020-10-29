#' link_files symlinks files in exe/* to the target path
#' @author Ben Fasoli
#'
#' @param from location of files
#' @param to location to create links
#'
#' @export

link_files <- function(from, to) {
  files <- dir(from, full.names = T)
  suppressWarnings(file.symlink(files, to))
}
