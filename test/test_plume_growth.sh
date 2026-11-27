#!/bin/bash
# Checks that calc_plume_dilution() grows the near-field plume by Taylor's
# spread, sigma_z^2 = 2 sigw^2 (tlgr t + tlgr^2 exp(-t / tlgr) - tlgr^2), which
# is never wider than sigw t
set -e

Rscript - <<'REOF'
suppressMessages(library(dplyr))
source('r/src/calc_plume_dilution.r')

# One backward particle, one minute old, well below the mixing height
p <- data.frame(indx = 1, time = -1, foot = 0.001, samt = 1, sigw = 0.5,
                tlgr = 100, dens = 1.2, mlht = 1000)
out <- calc_plume_dilution(p, 1, 5, veght = 0.5)

# Recover the plume depth from the corrected foot
plume <- 0.02897 / (out$foot * p$dens) * p$samt * 60
t <- 60
taylor <- sqrt(2) * p$sigw * sqrt(p$tlgr * t + p$tlgr^2 * exp(-t / p$tlgr) -
                                    p$tlgr^2)
stopifnot(isTRUE(all.equal(plume - 5, taylor)),
          plume - 5 <= p$sigw * t)

# At the release the plume has no depth to add
p0 <- mutate(p, time = 0)
stopifnot(isTRUE(all.equal(
  0.02897 / (calc_plume_dilution(p0, 1, 5, veght = 0.5)$foot * p$dens) * 60, 5)))
REOF

echo "plume growth test successful"
