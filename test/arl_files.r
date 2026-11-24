# Helpers for the tests: small ARL files written in R.

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
