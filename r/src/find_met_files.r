#' find_met_files searches for meteorological data files
#' @author Ben Fasoli & James Mineau
#'
#' Searches for available meteorological files matching the given strftime
#' compatible file naming convention
#'
#' @param t_start time of simulation start
#' @param n_hours number of hours to run each simulation; negative indicates
#'   backward in time
#' @param met_path directory to find meteorological data
#' @param met_file_format grep compatible file naming convention to identify
#'   meteorological data files necessary for the timing of the simulation
#'   indicated by \code{t_start} and \code{n_hours}
#' @param met_file_tres time resolution of meteorological data files
#'
#' @import dplyr
#' @import lubridate
#' @export

find_met_files <- function(t_start, n_hours, met_path,
                           met_file_format, met_file_tres) {
  require(dplyr)
  require(lubridate)

  # Simulation timing
  sim_start <- as.POSIXct(t_start, tz = 'UTC')
  sim_end <- sim_start + as.difftime(n_hours, units = 'hours')

  # Generate the hours to search for
  met_start <- floor_date(min(sim_start, sim_end), unit = met_file_tres)
  met_end <- max(sim_start, sim_end)
  met_end_ceil <- ceiling_date(met_end, unit = met_file_tres)
  if (n_hours < 0 
      && (hour(met_end) + 1) %% 24 == hour(met_end_ceil)
      && minute(met_end) > 0) {
    # If the end time is at the end of a met file,
    # add an hour to the end time to include the next file.
    # This is necessary to interpolate the last hour of data.
    met_end <- met_end_ceil
  }
  met_times <- seq(met_start, met_end, by = met_file_tres)

  # Format the request and remove duplicates
  request <- met_times %>%
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

  unique(available[idx])  # Available files that match the request
}
