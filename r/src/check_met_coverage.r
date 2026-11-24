#' read_met_times reads the time steps an ARL meteorological file holds
#' @author James Mineau
#'
#' Reads the index record at the start of each time step of an ARL packed file.
#' Each time step is one index record followed by one data record per variable
#' per level, every record 50 + nx * ny bytes long, so the first index record
#' gives the length of a time step and the file size gives their number.
#'
#' @param path path to the ARL file
#'
#' @return list with \code{times} (POSIXct, UTC, one per time step) and
#'   \code{whole} (FALSE when the file size is not a whole number of time
#'   steps, as in a file cut off partway). \code{NULL} when the file does not
#'   start with an index record.
#'
#' @export
read_met_times <- function(path) {
  size <- file.size(path)
  con <- file(path, 'rb')
  on.exit(close(con))

  first <- read_index_record(con, 0)
  if (is.null(first)) return(NULL)

  record_length <- 50 + first$nx * first$ny
  step_length <- (1 + first$n_records) * record_length
  n_steps <- size %/% step_length

  times <- first$time
  for (i in seq_len(n_steps - 1)) {
    index <- read_index_record(con, i * step_length)
    if (is.null(index)) return(list(times = times, whole = FALSE))
    times <- c(times, index$time)
  }
  list(times = times, whole = size %% step_length == 0)
}


#' read_index_record reads one ARL index record
#' @author James Mineau
#'
#' @param con binary connection to an ARL file
#' @param position byte offset of the record
#'
#' @return list with \code{time}, \code{nx}, \code{ny} and \code{n_records}
#'   (the data records in its time step), or \code{NULL} when the bytes there
#'   are not an index record.
read_index_record <- function(con, position) {
  seek(con, where = position)
  # 50 byte label, 108 byte fixed part, then each level's variables
  raw <- readBin(con, 'raw', n = 158)
  if (length(raw) < 158) return(NULL)
  text <- rawToChar(replace(raw, raw == 0, charToRaw(' ')))
  if (substr(text, 15, 18) != 'INDX') return(NULL)

  field <- function(first, last) suppressWarnings(as.integer(substr(text, first, last)))
  year <- field(1, 2)
  year <- year + ifelse(year < 40, 2000, 1900)
  time <- ISOdatetime(year, field(3, 4), field(5, 6), field(7, 8), 0, 0,
                      tz = 'UTC') + 60 * field(58, 59)

  # The label's grid letters carry the thousands of nx and ny (A = 1000);
  # '@' or a digit means none
  thousands <- function(char) {
    k <- match(char, LETTERS)
    if (is.na(k)) 0 else 1000 * k
  }
  nx <- thousands(substr(text, 13, 13)) + field(144, 146)
  ny <- thousands(substr(text, 14, 14)) + field(147, 149)
  nz <- field(150, 152)
  if (anyNA(c(time, nx, ny, nz))) return(NULL)

  # Each level: height (6), number of variables (2), then 8 per variable
  n_records <- 0
  for (k in seq_len(nz)) {
    level <- rawToChar(readBin(con, 'raw', n = 8))
    n_vars <- suppressWarnings(as.integer(substr(level, 7, 8)))
    if (is.na(n_vars)) return(NULL)
    n_records <- n_records + n_vars
    if (n_vars > 0) readBin(con, 'raw', n = 8 * n_vars)
  }

  list(time = time, nx = nx, ny = ny, n_records = n_records)
}


#' check_met_coverage checks that met files hold every time step a run needs
#' @author James Mineau
#'
#' HYSPLIT needs the met step at or before the run's first time, the one at or
#' after its last time, and every step between. When one is missing it stops
#' its particles where the met runs out, often without an error. The file
#' names cannot tell which time steps a file holds, so this reads them from
#' the files' index records.
#'
#' @param met_files paths of the files the run will read
#' @param t_start time of simulation start
#' @param n_hours number of hours to run each simulation; negative indicates
#'   backward in time
#'
#' @return \code{NULL} when the files cover the run, else a message saying
#'   what is missing.
#'
#' @export
check_met_coverage <- function(met_files, t_start, n_hours) {
  sim_start <- as.POSIXct(t_start, tz = 'UTC')
  sim_end <- sim_start + as.difftime(n_hours, units = 'hours')
  first <- min(sim_start, sim_end)
  last <- max(sim_start, sim_end)
  fmt <- function(t) strftime(t, '%Y-%m-%d %H:%M', tz = 'UTC')

  times <- list()
  for (f in met_files) {
    met <- read_met_times(f)
    if (is.null(met))
      return(paste('Not an ARL meteorological file:', f))
    if (!met$whole)
      return(paste('Meteorological file cut off partway through a time step:', f))
    times[[f]] <- met$times
  }
  times <- sort(unique(do.call(c, unname(times))))

  if (length(times) < 2)
    return(paste('The meteorological files hold fewer than two time steps:',
                 paste(met_files, collapse = ', ')))
  if (times[1] > first)
    return(paste('No meteorological time step at or before', fmt(first),
                 '- the files start at', fmt(times[1])))
  if (times[length(times)] < last)
    return(paste('No meteorological time step at or after', fmt(last),
                 '- the files end at', fmt(times[length(times)])))

  # A gap between the steps the run reads, wider than the met's time step
  used <- times[max(which(times <= first)):min(which(times >= last))]
  gaps <- diff(as.numeric(used))
  step <- min(diff(as.numeric(times)))
  wide <- which(gaps > step)
  if (length(wide) > 0)
    return(paste('No meteorological time steps between',
                 fmt(used[wide[1]]), 'and', fmt(used[wide[1] + 1])))

  NULL
}
