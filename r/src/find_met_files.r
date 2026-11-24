#' list_met_files lists the files in a meteorological data directory
#' @author Ben Fasoli & James Mineau
#'
#' @param met_path directory to find meteorological data
#'
#' @return paths of the files under \code{met_path}, each physical file once,
#'   without lock files
#'
#' @export
list_met_files <- function(met_path) {
  available <- dir(met_path, full.names = T, recursive = T)
  available <- available[!grepl('.lock', available)]

  # Some met archives expose the same file via multiple paths -- e.g. a
  # top-level symlink pointing into a YYYY/MM/ subdirectory.  recursive=T
  # finds both the symlink and the real file, and downstream grep yields
  # duplicate entries that produce a malformed HYSPLIT CONTROL file.
  # Dedupe by resolved path so each physical file is listed once.
  available[!duplicated(normalizePath(available, mustWork = FALSE))]
}
