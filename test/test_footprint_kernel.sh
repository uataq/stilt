#!/bin/bash
# Checks that permute.f90 spreads a particle over the smoothing kernel's own
# cells, longitude along the first index of the footprint and latitude along
# the second, including kernels that are wider in one direction (xres != yres)
set -e

Rscript - <<'REOF'
source('r/src/calc_footprint.r')
dyn.load('r/src/permute.so')

for (res in list(c(0.01, 0.01), c(0.01, 0.05), c(0.05, 0.01), c(0.02, 0.025))) {
  k <- make_gauss_kernel(res, 0.1, '+proj=longlat')
  nlon <- 201L
  nlat <- 151L
  loi <- 101L
  lai <- 76L

  # Called as in calc_footprint(), where nax is the latitude count
  out <- .Fortran('permute', ans = matrix(0, nrow = nlon, ncol = nlat),
                  nax = nlat, nay = nlon, k = k, nkx = ncol(k), nky = nrow(k),
                  len = 1L, lai = lai, loi = loi, foot = 1)

  xh <- (ncol(k) - 1L) / 2L
  yh <- (nrow(k) - 1L) / 2L
  expected <- matrix(0, nrow = nlon, ncol = nlat)
  expected[(loi - xh):(loi + xh), (lai - yh):(lai + yh)] <- t(k)

  message(sprintf('xres %g yres %g: kernel %d x %d, total %g',
                  res[1], res[2], ncol(k), nrow(k), sum(out$ans)))
  stopifnot(isTRUE(all.equal(out$ans, expected)))
}
REOF

echo "footprint kernel test successful"
