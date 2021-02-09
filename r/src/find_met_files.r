#' find_met_files searches for meteorological data files
#' @author Ben Fasoli
#'
#' Searches for available meteorological files matching the given strftime
#' compatible file naming convention
#'
#' @param t_start time of simulation start
#' @param met_file_format grep compatible file naming convention to identify
#'   meteorological data files necessary for the timing of the simulation
#'   indicated by \code{t_start} and \code{n_hours}
#' @param n_hours number of hours to run each simulation; negative indicates
#'   backward in time
#' @param met_path directory to find meteorological data
#'
#' @import dplyr
#' @export

find_met_files <- function(t_start, met_file_format, n_hours, met_path) {
  require(dplyr)
  
  is_backward <- n_hours < 0
  
  # TODO: implement n_hours_per_met_file to better determine file names at
  # varying time resolutions
  request <- as.POSIXct(t_start, tz='UTC') %>%
    c(. + c(1, -1, n_hours, is_backward * (n_hours - 5)) * 3600) %>%
    range() %>%
    (function(x) seq(x[1], x[2], by = 'hour')) %>%
    strftime(tz = 'UTC', format = met_file_format)
  
  available <- dir(met_path, full.names = T)
  available <- available[!grepl('\\.lock$', available)]
  
  idx <- do.call(c, lapply(request, function(pattern) {
    grep(pattern = pattern, x = available)
  }))
  
  if (any(idx < 1))
    return()
  
  unique(available[idx])
}
