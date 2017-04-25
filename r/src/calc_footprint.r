#' calc_footprint generates upstream influence footprint
#' @author Ben Fasoli
#'
#' Aggregates the upstream particle trajectories into a time integrated
#' footprint, expanding particle influence using variable 2d gaussian kernels
#' with bandwidths proportional to the mean pairwise distance between all
#' particles at each time step.
#'
#' @param p data frame containing particle trajectories, typically obtained
#'   from PARTICLE.rds but can be derived from PARTICLE.dat with column names
#'   equivalent to \code{varsiwant}. Must contain colums specifying long, lati,
#'   indx, foot, and time.
#' @param output filename argument passed to \code{raster::writeRaster}
#' @param n_cores_grid number of cores to use in parallel while calculating
#'   gaussian kernel matrices
#' @param xmn sets grid start longitude
#' @param xmx sets grid end longitude
#' @param xres resolution for longitude grid
#' @param ymn sets grid start latitude
#' @param ymx sets grid end latitude
#' @param yres resolution for latitude grid
#'
#' @import dplyr, parallel, raster, uataq
#' @export

calc_footprint <- function(p, output = NULL, n_cores_grid = 1,
                           xmn = -180, xmx = 180, xres = 0.1,
                           ymn = -90, ymx = 90, yres = xres) {
  require(dplyr)
  require(parallel)
  require(raster)
  require(uataq)

  glong <- seq(xmn, xmx, by = xres)
  glati <- seq(ymn, ymx, by = yres)

  p <- p %>%
    filter(long >= min(glong), long <= max(glong),
           lati >= min(glati), lati <= max(glati))

  # Interpolate particle locations
  times <- c(seq(0, -10, by = -0.1),
             seq(-10.2, -20, by = -0.2),
             seq(-20.5, -100, by = -0.5))
  i <- p %>%
    dplyr::select(indx, time, long, lati, foot) %>%
    full_join(expand.grid(time = times,
                          indx = unique(p$indx)), by = c('indx', 'time')) %>%
    arrange(indx, -time) %>%
    group_by(indx) %>%
    mutate(long = na_interp(long, x = time),
           lati = na_interp(lati, x = time),
           foot = na_interp(foot, x = time)) %>%
    ungroup() %>%
    na.omit() %>%
    mutate(time = round(time, 1))

  # Scale interpolated values to retain total field
  mi <- i$time >= -10
  mp <- p$time >= -10
  i$foot[mi] <- i$foot[mi] / (sum(i$foot[mi], na.rm = T) / sum(p$foot[mp], na.rm = T))
  mi <- i$time < -10 & i$time >= -20
  mp <- p$time < -10 & p$time >= -20
  i$foot[mi] <- i$foot[mi] / (sum(i$foot[mi], na.rm = T) / sum(p$foot[mp], na.rm = T))
  mi <- i$time < -20 & i$time >= -100
  mp <- p$time < -20 & p$time >= -100
  i$foot[mi] <- i$foot[mi] / (sum(i$foot[mi], na.rm = T) / sum(p$foot[mp], na.rm = T))

  particle <- i
  np <- max(particle$indx, na.rm = T)

  # Calculate pairwise distances
  pd <- particle %>%
    group_by(time) %>%
    summarize(dist = matrix(c(long, lati), ncol = 2) %>%
                dist() %>%
                mean(na.rm = T),
              lati = mean(lati, na.rm = T))

  # Generate gaussian kernels
  make_gauss_kernel <- function (rs, sigma) {
    # Modified from raster:::.Gauss.weight()
    d <- 3 * sigma
    nx <- 1 + 2 * floor(d/rs[1])
    ny <- 1 + 2 * floor(d/rs[2])
    m <- matrix(ncol = nx, nrow = ny)
    xr <- (nx * rs[1])/2
    yr <- (ny * rs[2])/2
    r <- raster(m, xmn = -xr[1], xmx = xr[1], ymn = -yr[1], ymx = yr[1],
                crs = "+proj=utm +zone=1 +datum=WGS84")
    p <- xyFromCell(r, 1:ncell(r))^2
    m <- 1/(2 * pi * sigma^2) * exp(-(p[, 1] + p[, 2])/(2 * sigma^2))
    m <- matrix(m, ncol = nx, nrow = ny, byrow = TRUE)
    m/sum(m)
  }

  message('Gridding to ', xres, 'x', yres, ' with ', n_cores_grid, ' threads...')
  cl <- makeForkCluster(n_cores_grid, outfile = '')

  message('Generating gaussian particles...')
  gk <- parLapply(cl, pd$time, function(x) {
    traj <- particle %>%
      filter(time == x)

    if (nrow(traj) < 2) return()

    traj <- traj %>%
      dplyr::select(long, lati, foot) %>%
      group_by(long, lati) %>%
      summarize(foot = sum(foot, na.rm = T)) %>%
      ungroup()

    res <- c(xres, yres)

    # d gives gaussian bandwidth sensitivity scalar
    # lower numbers yield larger gaussian matrices
    # higher numbers allow faster computation
    # pd$dist is mean pairwise distance within particle cloud at timestep
    idx <- pd$time == x
    d <- pd$dist[idx] / (40 * cos(pd$lati[idx] * pi/180)) + max(res) / 8
    k <- make_gauss_kernel(res, d)
    nk <- nrow(k)

    if (nk < 2)
      traj %>%
      transmute(long = glong[uataq::find_neighbor(long, glong + xres / 2)],
                lati = glati[uataq::find_neighbor(lati, glati + yres / 2)],
                foot) %>%
      group_by(long, lati) %>%
      summarize(foot = sum(foot, na.rm = T)) %>%
      ungroup() %>%
      return()

    lo <- 1:nk - ((nk + 1) / 2)

    df <- lapply(1:nk, function(i) {
      lapply(1:nk, function(j) {
        xi <- lo[i]
        yi <- lo[j]
        ki <- k[i, j]
        traj %>%
          transmute(long = long + xi * xres,
                    lati = lati + yi * yres,
                    foot = foot * ki)

      })
    })

    lapply(df, bind_rows) %>%
      bind_rows() %>%
      mutate(long = glong[uataq::find_neighbor(long, glong + xres / 2)],
             lati = glati[uataq::find_neighbor(lati, glati + yres / 2)]) %>%
      group_by(long, lati) %>%
      summarize(foot = sum(foot, na.rm = T)) %>%
      ungroup() %>%
      return()
  })

  stopCluster(cl)

  message('Aggregating kernels...')
  gk <- bind_rows(gk) %>%
    group_by(long, lati) %>%
    summarize(foot = sum(foot, na.rm = T) / np) %>%
    ungroup() %>%
    rasterFromXYZ(crs = '+proj=longlat +ellps=WGS84')

  if (!is.null(output))
    raster::writeRaster(gk, output, overwrite = T)

  return(gk)
}
