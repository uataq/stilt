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
  needed <- seq(met_start, met_end, by = met_file_tres) %>%
    strftime(tz = 'UTC', format = met_file_format) %>%
    unique()

  # A run that ends partway through a file may need that file's next step, which
  # is the first step of the next file when the end is past the file's last step
  # (after 21:00 in a daily file of 3-hourly steps). Which file holds that step
  # depends on the met's time step, so include the next file in either direction
  # when it exists, without requiring it.
  next_file <- ceiling_date(met_end, unit = met_file_tres) %>%
    strftime(tz = 'UTC', format = met_file_format)
  request <- unique(c(needed, next_file))

  # Find the available files
  available <- dir(met_path, full.names = T, recursive = T)
  available <- available[!grepl('.lock', available)]

  # Some met archives expose the same file via multiple paths -- e.g. a
  # top-level symlink pointing into a YYYY/MM/ subdirectory.  recursive=T
  # finds both the symlink and the real file, and downstream grep yields
  # duplicate entries that produce a malformed HYSPLIT CONTROL file.
  # Dedupe by resolved path so each physical file is listed once.
  available <- available[!duplicated(normalizePath(available, mustWork = FALSE))]

  # Find the files that match the request
  matches <- lapply(request, function(pattern) {
    grep(pattern = pattern, x = available)
  })

  # HYSPLIT stops its particles where the met runs out, so a run with a file
  # missing from its span would come out short with no error
  missing <- request[lengths(matches) == 0 & request %in% needed]
  if (length(missing) > 0) {
    warning('No met file in ', met_path, ' matching ',
            paste(missing, collapse = ', '))
    return()
  }

  idx <- do.call(c, matches)
  unique(available[idx])  # Available files that match the request
}
