#' find_met_files searches for meteorological data files
#' @author Ben Fasoli & James Mineau
#'
#' Searches for available meteorological files matching the given strftime
#' compatible file naming convention
#'
#' @param t_start time of simulation start
#' @param n_hours number of hours to run each simulation; negative indicates
#'   backward in time
#' @param n_hours_per_met_file number of hours of meteorological data in each
#'   met file
#' @param met_file_format grep compatible file naming convention to identify
#'   meteorological data files necessary for the timing of the simulation
#'   indicated by \code{t_start} and \code{n_hours}
#' @param met_path directory to find meteorological data
#'
#' @import dplyr
#' @export

find_met_files <- function(t_start, n_hours, n_hours_per_met_file,
                           met_file_format, met_path) {
  require(dplyr)

  ts <- as.POSIXct(t_start, tz = 'UTC')
  is_backward <- as.numeric(n_hours < 0)
  ib <- ifelse(is_backward, 1, -1)

  # Generate the list of files to search for
  request <- seq(
    ts - (ib * as.difftime(abs(n_hours) + ib, units = 'hours')),
    ts + (ib * as.difftime(n_hours_per_met_file - is_backward, units = 'hours')),
    by = ib * 3600
  ) %>%
    strftime(tz = 'UTC', format = met_file_format) %>%
    unique()

  # Find the available files
  available <- dir(met_path, full.names = T, recursive = T)
  available <- available[!grepl('.lock', available)]

  # Find the files that match the request
  idx <- do.call(c, lapply(request, function(pattern) {
    grep(pattern = pattern, x = available)
  }))

  if (any(idx < 1))
    return()

  available[idx]  # Available files that match the request
}
