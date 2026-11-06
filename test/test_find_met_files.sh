#!/bin/bash
# Checks that find_met_files lists every file a run needs, in either direction,
# for each product in the NOAA ARL archive. HYSPLIT needs the met steps from the
# one at or before the run's first time to the one at or after its last time.
set -e

Rscript - <<'REOF'
source('r/src/find_met_files.r')

six_hour_block <- function(t) floor_date(t, '6 hours')

# File names, time resolution, and hours between met steps of each product
archives <- list(
  hrrr = list(
    format = '%Y%m%d_%H', tres = '6 hours', step = 1,
    file = function(t) {
      b <- six_hour_block(t)
      paste0(strftime(b, '%Y%m%d_%H', tz = 'UTC'), '-',
             sprintf('%02d', hour(b) + 5), '_hrrr')
    }),
  hrrr.v1 = list(
    format = 'hysplit.%Y%m%d.%Hz.hrrra', tres = '6 hours', step = 1,
    file = function(t) strftime(six_hour_block(t), 'hysplit.%Y%m%d.%Hz.hrrra',
                                tz = 'UTC')),
  nam12 = list(
    format = '%Y%m%d_nam12', tres = '24 hours', step = 3,
    file = function(t) strftime(t, '%Y%m%d_nam12', tz = 'UTC')),
  nams = list(
    format = '%Y%m%d_hysplit.t00z.namsa$', tres = '24 hours', step = 1,
    file = function(t) strftime(t, '%Y%m%d_hysplit.t00z.namsa', tz = 'UTC')),
  nams.ak = list(
    format = '%Y%m%d_hysplit.t00z.namsa.AK', tres = '24 hours', step = 1,
    file = function(t) strftime(t, '%Y%m%d_hysplit.t00z.namsa.AK', tz = 'UTC')),
  nams.hi = list(
    format = '%Y%m%d_hysplit.t00z.namsa.HI', tres = '24 hours', step = 1,
    file = function(t) strftime(t, '%Y%m%d_hysplit.t00z.namsa.HI', tz = 'UTC')),
  gfs0p25 = list(
    format = '%Y%m%d_gfs0p25', tres = '24 hours', step = 3,
    file = function(t) strftime(t, '%Y%m%d_gfs0p25', tz = 'UTC')),
  gdas0p5 = list(
    format = '%Y%m%d_gdas0p5', tres = '24 hours', step = 3,
    file = function(t) strftime(t, '%Y%m%d_gdas0p5', tz = 'UTC')),
  narr = list(
    format = 'NARR%Y%m', tres = '1 month', step = 3,
    file = function(t) strftime(t, 'NARR%Y%m', tz = 'UTC')),
  reanalysis = list(
    format = 'RP%Y%m.gbl', tres = '1 month', step = 6,
    file = function(t) strftime(t, 'RP%Y%m.gbl', tz = 'UTC'))
)

utc <- function(x) as.POSIXct(x, tz = 'UTC')

# The end of a file, of a leap February, and of a year, every 30 minutes
releases <- c(
  seq(utc('2016-02-27 00:00:00'), utc('2016-03-02 00:00:00'), by = 1800),
  seq(utc('2016-12-30 00:00:00'), utc('2017-01-02 00:00:00'), by = 1800)
)
hours <- seq(utc('2016-01-25 00:00:00'), utc('2017-02-05 00:00:00'), by = 3600)

check <- function(name, a, n_hours) {
  met <- tempfile(name)
  dir.create(met)
  file.create(file.path(met, unique(a$file(hours))))

  for (i in seq_along(releases)) {
    t <- releases[i]
    span <- range(t, t + n_hours * 3600)
    step <- paste(a$step, 'hours')
    steps <- seq(floor_date(span[1], step), ceiling_date(span[2], step),
                 by = a$step * 3600)
    needed <- unique(a$file(steps))
    listed <- basename(suppressMessages(
      find_met_files(t, n_hours, met, a$format, a$tres)))
    missing <- setdiff(needed, listed)
    if (length(missing) > 0) {
      stop(sprintf('%s: %s %+d h needs %s, listed %s', name, format(t, tz = 'UTC'),
                   n_hours, paste(missing, collapse = ' '),
                   paste(listed, collapse = ' ')))
    }
  }
  unlink(met, recursive = TRUE)
}

suppressMessages(library(lubridate))
for (name in names(archives)) {
  for (n_hours in c(-24, -6, -1, 1, 6, 24)) check(name, archives[[name]], n_hours)
  message(name, ': ok')
}

# gdas1 is weekly (gdas1.oct16.w2), and strftime has no week of the month and no
# lowercase month, so a met_file_format can only name one month's files
met <- tempfile('gdas1')
dir.create(met)
weeks <- paste0('gdas1.oct16.w', 1:5)
file.create(file.path(met, weeks))
listed <- basename(find_met_files('2016-10-07 22:00:00', -6, met,
                                  'gdas1.oct16.w', '7 days'))
stopifnot(all(weeks[1:2] %in% listed))
message('gdas1: ok')

# A run that ends on a file boundary needs nothing past that step
met <- tempfile('boundary')
dir.create(met)
days <- paste0('201610', c('07', '08', '09', '10'))
file.create(file.path(met, days))
stopifnot(
  identical(basename(find_met_files('2016-10-09 00:00:00', -6, met, '%Y%m%d',
                                    '24 hours')), days[2:3]),
  identical(basename(find_met_files('2016-10-08 18:00:00', 6, met, '%Y%m%d',
                                    '24 hours')), days[2:3])
)
message('file boundary: ok')
REOF

echo "met file search test successful"
