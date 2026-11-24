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

source('test/arl_files.r')

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
# ... and one cut off inside its first time step
writeBin(readBin(f1, 'raw', 1000), cut)
stopifnot(says(check_met_coverage(cut, '2016-10-08 03:00:00', -2), 'cut off partway'))

# A file that is not ARL
junk <- file.path(met, 'junk')
writeLines(rep('not met', 100), junk)
stopifnot(says(check_met_coverage(junk, '2016-10-08 03:00:00', -2), 'Not an ARL'))
message('damaged files: ok')
REOF

echo "met coverage check test successful"
