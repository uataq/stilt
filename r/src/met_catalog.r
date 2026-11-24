#' met_catalog lists what each meteorological file in a directory holds
#' @author James Mineau
#'
#' Reads the first and last time step and the time step of every ARL file under
#' \code{met_path} from its index records (\code{read_met_times}), so files can
#' be found by the times they hold whatever they are named. Reading a large
#' archive takes minutes, so the table is saved and a later call reads only the
#' files that are new or whose size or modification time changed.
#'
#' The table is saved as \code{.stilt-met-catalog.csv} in \code{met_path} when
#' that is writable, so everyone who reads the directory can use it, and in
#' \code{$XDG_CACHE_HOME/stilt/met_catalogs/} (\code{~/.cache} by default)
#' otherwise, and rewritten only when it changed. A table in \code{met_path}
#' that the user cannot write is read as a starting point.
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
  writable <- file.access(met_path, 2) == 0
  cache <- if (writable) shared else file.path(
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
  catalog <- do.call('rbind', c(list(empty_met_catalog()),
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
