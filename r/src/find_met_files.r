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
  is_backward <- n_hours < 0
  met_bracket <- n_hours_per_met_file - 1  # ts can be in the middle of a met file

  # Generate the hours to search for
  if (is_backward) {
    met_hours <- seq(
      ts - as.difftime(abs(n_hours) + met_bracket, units = 'hours'),
      ts,
      by = 3600
    )
  } else {
    met_hours <- seq(
      ts - as.difftime(met_bracket, units = 'hours'),
      ts + as.difftime(n_hours, units = 'hours'),
      by = 3600
    )
  }

  # Format the request and remove duplicates
  request <- met_hours %>%
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
