#!/bin/bash
# Checks that calc_plume_dilution() starts the plume at the release height above
# ground when release heights are above sea level (kmsl = 1)
set -e

Rscript - <<'REOF'
suppressMessages(library(dplyr))
source('r/src/calc_plume_dilution.r')

# Two backward particles over rising terrain, oldest row first. The first output
# row of particle 1 is over 1471 m and that of particle 2 over 1473 m
p <- expand.grid(indx = 1:2, time = -c(10, 5, 2, 1)) %>%
  mutate(foot = 0.001, samt = 1, sigw = 0.2, tlgr = 20, dens = 1.1, mlht = 200,
         zsfc = 1469 + 2 * indx + 2 * (abs(time) - 1))
f <- function(...) calc_plume_dilution(p, 2, ..., veght = 0.5)$foot

# A point receptor gives the foot of the same heights above ground
twin <- calc_plume_dilution(mutate(p, xhgt = c(5, 3)[indx]), 2, NA, 0.5)
stopifnot(isTRUE(all.equal(f(1476, kmsl = 1), twin$foot)))

# The plume is below the mixing height, so the correction applies
stopifnot(max(f(1476, kmsl = 1)) > 0.001)

# So does a column, keeping its release heights above sea level
col_msl <- mutate(p, xhgt = c(1571, 1673)[indx])
col_agl <- mutate(p, xhgt = c(100, 200)[indx])
out <- calc_plume_dilution(col_msl, 2, c(1500, 1700), 0.5, kmsl = 1)
stopifnot(isTRUE(all.equal(out$foot,
                           calc_plume_dilution(col_agl, 2, c(0, 200), 0.5)$foot)),
          all(out$xhgt %in% c(1571, 1673)))

# A release height below the terrain starts at the ground
below <- f(1400, kmsl = 1)
stopifnot(all(is.finite(below)), all(below > 0),
          isTRUE(all.equal(below, f(0))))

# That start has no depth while the particle is above the counting height, which
# keeps HYSPLIT's foot rather than dividing by zero
above <- p %>% mutate(samt = ifelse(abs(time) == 1, 0, samt),
                      foot = ifelse(abs(time) == 1, 0, foot))
stopifnot(all(is.finite(calc_plume_dilution(above, 2, 1400, 0.5, kmsl = 1)$foot)))

# Above ground nothing changes, and zsfc is not read
stopifnot(isTRUE(all.equal(f(5), f(5, kmsl = 0))))
f_nozsfc <- calc_plume_dilution(select(p, -zsfc), 2, 5, 0.5)$foot
stopifnot(isTRUE(all.equal(f(5), f_nozsfc)))

# Above sea level zsfc is required
err <- tryCatch(calc_plume_dilution(select(p, -zsfc), 2, 1476, 0.5, kmsl = 1),
                error = conditionMessage)
stopifnot(is.character(err), grepl('varsiwant must include', err), grepl('zsfc', err))
REOF

echo "plume dilution test successful"
