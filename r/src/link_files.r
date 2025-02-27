#' link_files symlinks files to the target path using relative paths
#' @author Ben Fasoli & James Mineau
#'
#' @param from location of files
#' @param to location to create links
#'
#' @import R.utils
#' @export

link_files <- function(from, to) {
  require(R.utils)

  # # Get all files in `from`
  x <- lapply(from, list.files, full.names = T)
  is_file <- sapply(x, function(x) length(x) == 0)
  files <- c(from[is_file], unlist(x[!is_file]))

  # Convert to paths relative to stilt_wd
  relative_files <- R.utils::getRelativePath(files, to)

  # Create the links
  suppressWarnings(file.symlink(relative_files, to))
}
