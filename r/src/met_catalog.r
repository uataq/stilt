#' met_catalog lists what each meteorological file in a directory holds
#' @author James Mineau
#'
#' Reads the first and last time step and the time step of every ARL file under
#' \code{met_path} from its index records (\code{read_met_times}), so files can
#' be found by the times they hold whatever they are named. Reading a large
#' archive takes minutes, so the table is saved and a later call reads only the
#' files that are new or whose size or modification time changed.
#'
#' The table is saved in \code{$XDG_CACHE_HOME/stilt/met_catalogs/}
#' (\code{~/.cache} by default), never in \code{met_path}, and rewritten only
#' when it changed. A \code{.stilt-met-catalog.csv} in \code{met_path}, put
#' there on purpose by whoever keeps the archive, is read as a starting point.
#'
#' @param met_path directory to find meteorological data
#' @param save save the table for later calls. \code{FALSE} reads every file
#'   and saves nothing.
#'
#' @return data frame with one row per file: \code{path}, \code{size},
#'   \code{mtime}, \code{first} and \code{last} (POSIXct, UTC), \code{step}
#'   (seconds), \code{n_steps}, and \code{error} (\code{''} for a whole ARL file,
#'   else the reason it cannot be used). Attribute \code{file} is where it was
#'   saved, or \code{NULL}.
#'
#' @export
met_catalog <- function(met_path, save = TRUE) {
  met_path <- normalizePath(met_path, mustWork = TRUE)
  shared <- file.path(met_path, '.stilt-met-catalog.csv')
  cache <- file.path(
    Sys.getenv('XDG_CACHE_HOME', file.path(Sys.getenv('HOME'), '.cache')),
    'stilt', 'met_catalogs', paste0(gsub('[^A-Za-z0-9]+', '_', met_path), '.csv'))

  # Rows already read: the user's own table, else the shared one
  known <- NULL
  if (save) {
    for (f in unique(c(cache, shared))) {
      if (file.exists(f)) {
        known <- read_met_catalog(f)
        break
      }
    }
  }

  paths <- list_met_files(met_path)
  info <- file.info(paths)
  size <- info$size
  # Whole seconds, as the saved table keeps them
  mtime <- floor(as.numeric(info$mtime))

  # Keep the rows of files unchanged since they were read; read the rest
  j <- if (is.null(known)) rep(NA_integer_, length(paths)) else match(paths, known$path)
  same <- !is.na(j)
  same[same] <- known$size[j[same]] == size[same] & known$mtime[j[same]] == mtime[same]
  fresh <- lapply(which(!same), function(i) read_catalog_row(paths[i], size[i], mtime[i]))
  catalog <- do.call(rbind, c(list(empty_met_catalog()),
                              if (any(same)) list(known[j[same], ]),
                              fresh))
  catalog <- catalog[order(catalog$path), ]
  rownames(catalog) <- NULL

  # Rewrite only when something changed, so a batch of CLI calls doesn't
  # rewrite the same table every call
  changed <- is.null(known) || !all(same) || nrow(known) != sum(same) ||
    !file.exists(cache)
  file <- NULL
  if (save && changed) {
    dir.create(dirname(cache), recursive = TRUE, showWarnings = FALSE)
    write_met_catalog(catalog, cache)
  }
  if (save) file <- cache
  attr(catalog, 'file') <- file
  catalog
}


#' read_met_catalog reads a table saved by met_catalog
#'
#' @param file path of the saved table
#'
#' @return data frame as from \code{met_catalog}
#'
#' @export
read_met_catalog <- function(file) {
  catalog <- read.csv(file, colClasses = c(path = 'character', error = 'character'),
                      na.strings = 'NA', stringsAsFactors = FALSE)
  catalog$error[is.na(catalog$error)] <- ''
  catalog$first <- as.POSIXct(catalog$first, tz = 'UTC', format = '%Y-%m-%d %H:%M:%S')
  catalog$last <- as.POSIXct(catalog$last, tz = 'UTC', format = '%Y-%m-%d %H:%M:%S')
  catalog
}


write_met_catalog <- function(catalog, file) {
  out <- catalog
  out$first <- strftime(out$first, '%Y-%m-%d %H:%M:%S', tz = 'UTC')
  out$last <- strftime(out$last, '%Y-%m-%d %H:%M:%S', tz = 'UTC')
  # Written aside and renamed, so a reader never sees half a table
  tmp <- paste0(file, '.', Sys.getpid(), '.tmp')
  write.csv(out, tmp, row.names = FALSE)
  file.rename(tmp, file)
}


empty_met_catalog <- function() {
  data.frame(path = character(0), size = numeric(0), mtime = numeric(0),
             first = as.POSIXct(character(0), tz = 'UTC'),
             last = as.POSIXct(character(0), tz = 'UTC'),
             step = numeric(0), n_steps = integer(0), error = character(0),
             stringsAsFactors = FALSE)
}


read_catalog_row <- function(path, size, mtime) {
  row <- data.frame(path = path, size = size, mtime = mtime,
                    first = as.POSIXct(NA, tz = 'UTC'),
                    last = as.POSIXct(NA, tz = 'UTC'),
                    step = NA_real_, n_steps = 0L, error = '',
                    stringsAsFactors = FALSE)
  met <- tryCatch(read_met_times(path), error = function(e) NULL)
  if (is.null(met)) {
    row$error <- 'not an ARL file'
    return(row)
  }
  times <- met$times
  row$first <- times[1]
  row$last <- times[length(times)]
  row$n_steps <- length(times)
  if (length(times) > 1) row$step <- min(diff(as.numeric(times)))
  if (!met$whole) row$error <- 'cut off partway through a time step'
  row
}


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
