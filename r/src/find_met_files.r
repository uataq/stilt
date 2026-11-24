#' find_met_files finds the meteorological files a simulation needs
#' @author Ben Fasoli & James Mineau
#'
#' HYSPLIT needs the met step at or before the run's first time, the one at or
#' after its last time, and every step between. These are the files of
#' \code{catalog} whose time steps reach into that span, in time order.
#'
#' @param t_start time of simulation start
#' @param n_hours number of hours to run each simulation; negative indicates
#'   backward in time
#' @param catalog table from \code{met_catalog}
#'
#' @return paths of the files, or \code{character(0)} with attribute
#'   \code{reason} saying what is missing
#'
#' @export
find_met_files <- function(t_start, n_hours, catalog) {
  sim_start <- as.POSIXct(t_start, tz = 'UTC')
  sim_end <- sim_start + as.difftime(n_hours, units = 'hours')
  first <- min(sim_start, sim_end)
  last <- max(sim_start, sim_end)
  fmt <- function(t) strftime(t, '%Y-%m-%d %H:%M', tz = 'UTC')
  none <- function(...) structure(character(0), reason = paste0(...))

  usable <- catalog[catalog$error == '' & catalog$n_steps > 0, ]
  if (nrow(usable) == 0)
    return(none('No ARL meteorological files found'))
  step <- ifelse(is.na(usable$step), 0, usable$step)

  # The latest step at or before first, and the earliest at or after last,
  # in any file
  f0 <- as.numeric(usable$first)
  f1 <- as.numeric(usable$last)
  before <- ifelse(f0 > as.numeric(first), NA,
                   pmin(f1, f0 + ifelse(step > 0, floor((as.numeric(first) - f0) / step) * step, 0)))
  after <- ifelse(f1 < as.numeric(last), NA,
                  pmax(f0, f0 + ifelse(step > 0, ceiling((as.numeric(last) - f0) / step) * step, 0)))
  if (all(is.na(before)))
    return(none('No meteorological file holds a time step at or before ', fmt(first)))
  if (all(is.na(after)))
    return(none('No meteorological file holds a time step at or after ', fmt(last)))
  from <- max(before, na.rm = TRUE)
  to <- min(after, na.rm = TRUE)

  hit <- usable[f1 >= from & f0 <= to, ]
  hit <- hit[order(hit$first), ]

  # Two files holding the same times are two products, or two copies, in one
  # directory, and HYSPLIT would read whichever it was given first. One shared
  # step at a boundary (a file ending at 06:00 and the next starting there) is
  # one product.
  if (nrow(hit) > 1) {
    overlap <- which(as.numeric(hit$first[-1]) < as.numeric(hit$last[-nrow(hit)]))
    if (length(overlap) > 0) {
      k <- overlap[1]
      return(none('Two meteorological files hold the same times: ',
                  hit$path[k], ', ', hit$path[k + 1]))
    }
  }
  hit$path
}


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
