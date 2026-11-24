#!/bin/bash
# Checks that check_met_coverage reads the time steps of ARL files from their
# index records and refuses a run the files do not cover: a missing step at
# either end, a gap between steps, a file cut off partway, and a file that is
# not ARL. The files are small ARL files written here, with one level of two
# variables on a 20 x 20 grid (or 1020 x 3, to read the grid's thousands).
set -e

Rscript - <<'REOF'
source('r/src/check_met_coverage.r')

utc <- function(x) as.POSIXct(x, tz = 'UTC')
# A message that says pattern; a NULL message (covered) does not
says <- function(msg, pattern) isTRUE(grepl(pattern, msg))

# Write an ARL file with one time step per time in times. Each step is an
# index record and one data record per variable, all 50 + nx * ny bytes.
write_arl <- function(path, times, nx = 20, ny = 20, vars = c('PRSS', 'T02M')) {
  record_length <- 50 + nx * ny
  thousands <- function(n) if (n >= 1000) LETTERS[n %/% 1000] else '@'
  grid <- if (nx < 1000 && ny < 1000) '99' else paste0(thousands(nx), thousands(ny))
  pad <- function(s) {
    stopifnot(nchar(s) <= record_length)
    charToRaw(formatC(s, width = -record_length))
  }
  label <- function(t, var) {
    sprintf('%2d%2d%2d%2d%2d%2d%2s%4s%4d%14.7E%14.7E',
            as.integer(format(t, '%y', tz = 'UTC')),
            as.integer(format(t, '%m', tz = 'UTC')),
            as.integer(format(t, '%d', tz = 'UTC')),
            as.integer(format(t, '%H', tz = 'UTC')),
            0, 0, grid, var, 0, 0, 0)
  }
  con <- file(path, 'wb')
  on.exit(close(con))
  for (i in seq_along(times)) {
    t <- times[i]
    levels <- paste0(sprintf('%6.1f%2d', 0, length(vars)),
                     paste(sprintf('%-4s%3d ', vars, 0), collapse = ''))
    fixed <- paste0(sprintf('%-4s%3d%2d', 'TEST', 0, 0),
                    paste(rep(sprintf('%7.2f', 0), 12), collapse = ''),
                    sprintf('%3d%3d%3d%2d%4d', nx %% 1000, ny %% 1000, 1, 1,
                            108 + nchar(levels)))
    writeBin(pad(paste0(label(t, 'INDX'), fixed, levels)), con)
    for (v in vars) writeBin(pad(label(t, v)), con)
  }
}

met <- tempfile('met')
dir.create(met)
f1 <- file.path(met, 'a')
f2 <- file.path(met, 'b')
hours <- function(from, n) utc(from) + 3600 * (seq_len(n) - 1)
write_arl(f1, hours('2016-10-08 00:00:00', 6))
write_arl(f2, hours('2016-10-08 06:00:00', 6))

# The times read back, and a grid of more than 1000 points
stopifnot(identical(read_met_times(f1)$times, hours('2016-10-08 00:00:00', 6)))
big <- file.path(met, 'big')
write_arl(big, hours('2016-10-08 00:00:00', 3), nx = 1020, ny = 3)
stopifnot(identical(read_met_times(big)$times, hours('2016-10-08 00:00:00', 3)))
message('read times: ok')

# Covered, backward and forward, ending on a step and between steps
stopifnot(
  is.null(check_met_coverage(c(f1, f2), '2016-10-08 11:00:00', -11)),
  is.null(check_met_coverage(c(f1, f2), '2016-10-08 10:30:00', -10)),
  is.null(check_met_coverage(c(f1, f2), '2016-10-08 00:00:00', 11)),
  is.null(check_met_coverage(c(f1, f2), '2016-10-08 00:30:00', 10))
)
message('covered: ok')

# A forward run that ends after the last step of its last file
msg <- check_met_coverage(f1, '2016-10-08 01:00:00', 4.5)
stopifnot(says(msg, 'at or after 2016-10-08 05:30'))
# A backward run that starts before the first step
msg <- check_met_coverage(f2, '2016-10-08 08:00:00', -3)
stopifnot(says(msg, 'at or before 2016-10-08 05:00'))
message('ends: ok')

# A time step missing inside a file
gap <- file.path(met, 'gap')
write_arl(gap, hours('2016-10-08 00:00:00', 6)[-3])
msg <- check_met_coverage(gap, '2016-10-08 05:00:00', -5)
stopifnot(says(msg, 'between 2016-10-08 01:00 and 2016-10-08 03:00'))
# The same gap outside the run's span does not matter
stopifnot(is.null(check_met_coverage(gap, '2016-10-08 05:00:00', -2)))
message('gaps: ok')

# A file cut off partway through a time step
cut <- file.path(met, 'cut')
invisible(file.copy(f1, cut))
size <- file.size(cut)
writeBin(readBin(f1, 'raw', size)[1:(size - 100)], cut)
stopifnot(says(check_met_coverage(cut, '2016-10-08 03:00:00', -2), 'cut off partway'))

# A file that is not ARL
junk <- file.path(met, 'junk')
writeLines(rep('not met', 100), junk)
stopifnot(says(check_met_coverage(junk, '2016-10-08 03:00:00', -2), 'Not an ARL'))
message('damaged files: ok')
REOF

echo "met coverage check test successful"
