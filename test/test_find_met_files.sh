#!/bin/bash
# Checks that find_met_files lists the file holding the met step after a run's
# end, in either direction, for met with steps longer than an hour
set -e

Rscript - <<'REOF'
source('r/src/find_met_files.r')

met <- tempfile()
dir.create(met)
days <- paste0('201610', c('07', '08', '09'))
file.create(file.path(met, days))

cases <- list(
  # Daily files of 3-hourly steps: 20161008 ends at 21:00
  list('2016-10-08 22:00:00', -6, days[2:3]),
  list('2016-10-08 21:30:00', -6, days[2:3]),
  list('2016-10-08 23:00:00', -6, days[2:3]),
  list('2016-10-08 22:00:00', -24, days[1:3]),
  list('2016-10-08 16:00:00', 6, days[2:3]),
  # A run that ends on a file boundary needs nothing past that step
  list('2016-10-09 00:00:00', -6, days[2:3]),
  list('2016-10-08 18:00:00', 6, days[2:3])
)

for (case in cases) {
  files <- basename(find_met_files(case[[1]], case[[2]], met, '%Y%m%d',
                                   '24 hours'))
  message(sprintf('%s %+d h: %s', case[[1]], case[[2]],
                  paste(files, collapse = ' ')))
  stopifnot(identical(files, case[[3]]))
}
REOF

echo "met file search test successful"
