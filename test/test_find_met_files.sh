#!/bin/bash
# Checks that find_met_files picks, from a catalog of the files' headers, every
# file a run needs and no other, in either direction, for met laid out as NOAA
# ARL's products are: 6-hour files of hourly steps (HRRR), daily files of
# 3-hourly steps (NAM, GFS), and weekly files of 3-hourly steps with names no
# date format can match (gdas1). HYSPLIT needs the met steps from the one at or
# before the run's first time to the one at or after its last time. Then checks
# that the catalog is saved, reused, updated and shared.
set -e

Rscript - <<'REOF'
for (f in c('r/src/check_met_coverage.r', 'r/src/find_met_files.r',
            'r/src/met_catalog.r', 'test/arl_files.r')) source(f)
Sys.setenv(XDG_CACHE_HOME = tempfile('cache'))

utc <- function(x) as.POSIXct(x, tz = 'UTC')
steps <- function(from, n, by) utc(from) + by * 3600 * (seq_len(n) - 1)

# Each layout: file names and the time steps each holds
layouts <- list(
  hrrr = lapply(0:15, function(i) {
    t <- utc('2016-02-27 00:00:00') + i * 6 * 3600
    list(name = paste0(strftime(t, '%Y%m%d_%H', tz = 'UTC'), '-xx_hrrr'),
         times = steps(t, 6, 1))
  }),
  daily = lapply(0:4, function(i) {
    t <- utc('2016-02-27 00:00:00') + i * 86400
    list(name = paste0(strftime(t, '%Y%m%d', tz = 'UTC'), '_nam12'),
         times = steps(t, 8, 3))
  }),
  weekly = lapply(1:3, function(w) {
    t <- utc('2016-10-01 00:00:00') + (w - 1) * 7 * 86400
    list(name = paste0('gdas1.oct16.w', w), times = steps(t, 56, 3))
  })
)

check <- function(name, files) {
  met <- tempfile(name)
  dir.create(met)
  for (f in files) write_arl(file.path(met, f$name), f$times)
  catalog <- met_catalog(met)
  all_times <- sort(do.call(c, lapply(files, `[[`, 'times')))

  for (t in seq(all_times[1] + 86400, all_times[length(all_times)] - 86400, by = 3600)) {
    for (n_hours in c(-24, -6, -1, 1, 6, 24)) {
      t0 <- utc('1970-01-01') + t
      ends <- sort(c(t0, t0 + n_hours * 3600))
      if (ends[1] < all_times[1] || ends[2] > all_times[length(all_times)]) next
      # Brute force: the steps the run needs, and the files that hold them
      from <- max(all_times[all_times <= ends[1]])
      to <- min(all_times[all_times >= ends[2]])
      needed <- all_times[all_times >= from & all_times <= to]
      expected <- sort(vapply(Filter(function(f) any(f$times %in% needed), files),
                              `[[`, '', 'name'))
      found <- find_met_files(t0, n_hours, catalog)
      if (!identical(sort(basename(found)), expected) ||
          !is.null(check_met_coverage(found, t0, n_hours))) {
        stop(name, ' at ', format(t0, tz = 'UTC'), ', n_hours ', n_hours,
             ': found ', paste(basename(found), collapse = ' '),
             '; expected ', paste(expected, collapse = ' '))
      }
    }
  }
  message(name, ': ok')
}
for (name in names(layouts)) check(name, layouts[[name]])

# Outside the files, and a gap
met <- tempfile('gap')
dir.create(met)
write_arl(file.path(met, 'a'), steps('2016-10-08 00:00:00', 6, 1))
write_arl(file.path(met, 'c'), steps('2016-10-08 12:00:00', 6, 1))
catalog <- met_catalog(met)
says <- function(x, pattern) isTRUE(grepl(pattern, attr(x, 'reason')))
stopifnot(
  says(find_met_files('2016-10-08 03:00:00', -6, catalog), 'at or before 2016-10-07 21:00'),
  says(find_met_files('2016-10-08 15:00:00', 6, catalog), 'at or after 2016-10-08 21:00'),
  isTRUE(grepl('between 2016-10-08 05:00 and 2016-10-08 12:00',
               check_met_coverage(find_met_files('2016-10-08 14:00:00', -10, catalog),
                                  '2016-10-08 14:00:00', -10)))
)
# One step shared at a boundary is one product
write_arl(file.path(met, 'b'), steps('2016-10-08 05:00:00', 8, 1))
stopifnot(length(find_met_files('2016-10-08 07:00:00', -4, met_catalog(met))) == 2)
# Two files holding the same times
write_arl(file.path(met, 'a_copy'), steps('2016-10-08 00:00:00', 6, 1))
stopifnot(says(find_met_files('2016-10-08 04:00:00', -2, met_catalog(met)),
               'Two meteorological files hold the same times'))
message('missing and doubled: ok')

# Saved in the user's cache, never in met_path, reused, and updated when a
# file changes
Sys.setenv(XDG_CACHE_HOME = tempfile('cache'))
met <- tempfile('saved')
dir.create(met)
write_arl(file.path(met, 'a'), steps('2016-10-08 00:00:00', 6, 1))
catalog <- met_catalog(met)
saved <- attr(catalog, 'file')
stopifnot(startsWith(saved, Sys.getenv('XDG_CACHE_HOME')),
          !file.exists(file.path(met, '.stilt-met-catalog.csv')),
          identical(read_met_catalog(saved)$first, utc('2016-10-08 00:00:00')))
# Not rewritten when nothing changed
t0 <- file.mtime(saved)
Sys.sleep(1.1)
invisible(met_catalog(met))
stopifnot(file.mtime(saved) == t0)
# A row is reused while the file's size and time are unchanged
table <- read.csv(saved, stringsAsFactors = FALSE)
table$n_steps <- 99
write.csv(table, saved, row.names = FALSE)
stopifnot(met_catalog(met)$n_steps == 99)
# and read again when the file changes
write_arl(file.path(met, 'a'), steps('2016-10-09 00:00:00', 3, 1))
catalog <- met_catalog(met)
stopifnot(catalog$n_steps == 3, catalog$first == utc('2016-10-09 00:00:00'))
message('saved: ok')

# A catalog put in met_path by its keeper is the starting point for a user
# with no cached copy
table <- read.csv(saved, stringsAsFactors = FALSE)
table$n_steps <- 77
write.csv(table, file.path(met, '.stilt-met-catalog.csv'), row.names = FALSE)
Sys.setenv(XDG_CACHE_HOME = tempfile('cache'))
catalog <- met_catalog(met)
stopifnot(catalog$n_steps == 77,
          startsWith(attr(catalog, 'file'), Sys.getenv('XDG_CACHE_HOME')))
message('shared: ok')
REOF

echo "met file search test successful"
