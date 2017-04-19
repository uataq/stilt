#' write_control writes out a namelist file to control the model behavior
#' @author Ben Fasoli
#'
#' @param receptor data frame containing columns for \code{run_times} as a
#'   POSIXct formatted timestamp, \code{lati} (degrees), \code{long} (degrees),
#'   and \code{zagl} (meters above ground level) as numeric specifications of
#'   the receptor location
#' @param n_hour number of hours to run each simulation; negative indicates
#'   backward in time
#' @param w_option vertical motion calculation method. 0: use vertical velocity
#'   from data, 1: isob, 2: isen, 3: dens, 4: sigma; defaults to 0
#' @param z_top top of model domain, in meters above ground level; defaults to
#'   25000.0
#' @param met_files full paths to meteorological data files in .arl format
#' @param file path and name for output file
#'
#' @import dplyr
#' @export

write_control <- function(receptor, n_hour, w_option = 0,
                          z_top = 25000, met_files, file = 'exe/CONTROL') {
  require(dplyr)

  if (!'CONTROL' %in% basename(file))
    stop('write_control(): file argument must end with CONTROL')

  receptor <- mutate(receptor,
                     print = paste(lati, long, zagl))
  n_loc <- nrow(receptor)

  txt = c(
    strftime(receptor$run_time, tz = 'UTC', format = '%y %m %d %H'),
    n_loc,
    receptor$print,
    n_hour,
    w_option,
    format(z_top, nsmall = 1),
    length(met_files),
    paste0(dirname(met_files), '/ \n', basename(met_files)),
    '1',
    'test',
    '1',
    '0.01',
    '00 00 00 00 00',
    '1',
    '0.0 0.0',
    '0.5 0.5',
    '30.0 30.0',
    './',
    'cdump',
    '1',
    '100',
    '00 00 00 00 00',
    '00 00 00 00 00',
    '00 2 00',
    '1',
    '0.0 0.0 0.0',
    '0.0 0.0 0.0 0.0 0.0',
    '0.0 0.0 0.0',
    '0.0',
    '0.0'
  )

  write(txt, file)
  return(file)
}
