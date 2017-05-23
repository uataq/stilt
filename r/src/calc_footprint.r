#' calc_footprint generates upstream influence footprint
#' @author Ben Fasoli
#'
#' Aggregates the upstream particle trajectories into a time integrated
#' footprint, expanding particle influence using variable 2d gaussian kernels
#' with bandwidths proportional to the mean pairwise distance between all
#' particles at each time step. Requires compiled permute.so to build the
#' gaussian kernels with fortran.
#'
#' @param p data frame containing particle trajectories, typically obtained
#'   from PARTICLE.rds but can be derived from PARTICLE.dat with column names
#'   equivalent to \code{varsiwant}. Must contain colums specifying long, lati,
#'   indx, foot, and time.
#' @param output filename argument. Must end with .nc (for ncdf output,
#'   preferred) or .rds (for serialized R data output). .nc files are saved
#'   in a format natively compatible with raster::brick() and raster::raster()
#'   while .rds files do not require any additional libraries and have better
#'   compression
#' @param xmn sets grid start longitude
#' @param xmx sets grid end longitude
#' @param xres resolution for longitude grid
#' @param ymn sets grid start latitude
#' @param ymx sets grid end latitude
#' @param yres resolution for latitude grid
#'
#' @import dplyr, raster, uataq
#' @export

calc_footprint <- function(p, output = NULL, r_run_time = NULL,
                           time_integrate = F, dist_factor = 1, time_factor = 1,
                           xmn = -180, xmx = 180, xres = 0.1,
                           ymn = -90, ymx = 90, yres = xres) {
  
  require(dplyr)
  require(raster)
  require(uataq)
  
  np <- max(p$indx, na.rm = T)
  
  glong <- head(seq(xmn, xmx, by = xres), -1)
  glati <- head(seq(ymn, ymx, by = yres), -1)
  
  # Interpolate particle locations during initial time steps
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
  
  # Bootstrap pairwise distance calculation
  calc_dist <- function(x, y) {
    df <- data_frame(x, y)
    foo <- function(df) {
      matrix(c(df$x, df$y), ncol = 2) %>%
        dist() %>%
        mean(na.rm = T)
    }
    mean(bootstrap(df, foo, size = 50, iter = 5), na.rm = T)
  }
  
  pd <- i %>%
    group_by(time) %>%
    summarize(dist = calc_dist(long, lati),
              lati = mean(lati, na.rm = T))
  
  # Generate gaussian kernels
  make_gauss_kernel <- function (rs, sigma) {
    # Modified from raster:::.Gauss.weight()
    require(raster)
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
  
  # Gaussian kernel bandwidth scaling
  calc_bandwidth <- function(dist, dist_factor, time, time_factor, lati, xyres) {
    (time_factor * log10(1 + abs(time))/60 + dist_factor * dist / 10) /
      (25 * cos(lati * pi/180)) + max(xyres) / 4
  }

  # Determine maximum kernel size
  xyres <- c(xres, yres)
  sigma <- max(calc_bandwidth(pd$dist, dist_factor, pd$time, time_factor,
                              pd$lati, xyres))
  max_k <- make_gauss_kernel(xyres, sigma)
  xbuf <- (ncol(max_k) - 1) / 2
  ybuf <- (nrow(max_k) - 1) / 2
  
  max_glong <- seq(xmn - (xbuf*xres), xmx + ((xbuf - 1)*xres), by = xres)
  max_glati <- seq(ymn - (ybuf*yres), ymx + ((ybuf - 1)*yres), by = yres)
  
  # Remove zero influence particles and positions outside of domain
  xyzt <- i %>%
    filter(foot > 0, long >= xmn, long <= xmx,
           lati >= ymn, lati <= ymx)
  pd <- pd[is.element(pd$time, xyzt$time), ]
  
  # Pre grid particle locations
  xyzt <- xyzt %>%
    transmute(loi = as.integer(findInterval(long, max_glong)),
              lai = as.integer(findInterval(lati, max_glati)),
              foot = foot,
              time = time) %>%
    group_by(loi, lai, time) %>%
    summarize(foot = sum(foot, na.rm = T)) %>%
    ungroup()
  
  grd <- matrix(0, ncol = length(max_glong), nrow = length(max_glati))
  
  # Build gaussian kernels by time step
  foot <- sapply(pd$time, simplify = 'array', function(x) {
    print(x)
    step <- xyzt %>%
      filter(time == x)
    
    if (nrow(step) < 2) {
      return(grd)
    }
    
    # Dispersion kernel
    idx <- pd$time == x
    sigma <- calc_bandwidth(pd$dist[idx], dist_factor, x, time_factor,
                            pd$lati[idx], xyres)
    k <- make_gauss_kernel(xyres, sigma)
    
    # Array dimensions
    len <- nrow(step)
    nkx <- ncol(k)
    nky <- nrow(k)
    nax <- ncol(grd)
    nay <- nrow(grd)
    
    # Call permute fortran subroutine to build and aggregate kernels
    out <- .Fortran('permute', ans = grd, nax = nax, nay = nay,
                    k = k, nkx = nkx, nky = nky,
                    len = len, lai = step$lai, loi = step$loi, foot = step$foot)
    
    foot <- out$ans
    return(foot)
  })
  
  # Remove added buffer around domain
  size <- dim(foot)
  foot <- array(foot, dim = size[c(2, 1, 3)])
  foot <- foot[xbuf:(size[2]-xbuf-1), ybuf:(size[1]-ybuf-1), ] / np
  
  # Time integrate footprint by aggregating across 3rd dimension
  if (time_integrate) {
    foot <- apply(foot, c(1, 2), sum)
    time_out <- as.numeric(r_run_time)
  } else {
    hid <- floor(pd$time / 60)
    foot <- sapply(unique(hid), simplify = 'array', function(hour) {
      apply(foot[ , ,hid == hour], c(1, 2), sum)
    })
    time_out <- as.numeric(r_run_time + unique(hid) * 3600)
  }
  
  # Format footprints to CF-1.4 convention and save to file
  if (is.null(output)) return()
  if (file.exists(output)) system(paste('rm', output))
  
  if (grepl('\\.nc$', output, ignore.case = T) &&
      'ncdf4' %in% names(sessionInfo()$otherPkgs)) {
    xdim <- ncdim_def('longitude', 'degrees E', glong + xres/2)
    ydim <- ncdim_def('latitude', 'degrees N', glati + yres/2)
    tdim <- ncdim_def('time', 'seconds since 1970-01-01',
                      as.numeric(time_out))
    fvar <- ncvar_def('Footprint', 'umol CO2 m-2 s-1',
                      list(xdim, ydim, tdim), -1)
    
    nc <- nc_create(output, fvar)
    ncvar_put(nc, fvar, foot)
    ncatt_put(nc, 'longitude', 'position', 'grid center')
    ncatt_put(nc, 'latitude', 'position', 'grid center')
    ncatt_put(nc, 'time', 'timezone', 'UTC')
    ncatt_put(nc, 0, 'crs', '+proj=longlat +ellpsWGS84')
    ncatt_put(nc, 0, 'crs_format', 'PROJ.4')
    ncatt_put(nc, 0, 'Conventions', 'CF-1.4')
    ncatt_put(nc, 0, 'Compatibility', 'raster::raster() and raster::brick()')
    ncatt_put(nc, 0, 'Documentation', 'benfasoli.github.io/stilt')
    return(output)
  }
  
  if (grepl('\\.rds$', output, ignore.case = T)) {
    out_custom <- list(
      longitude = list(unit = 'degrees E',
                       position = 'grid center',
                       values = glong + xres/2),
      latitude = list(unit = 'degrees N',
                      position = 'grid center',
                      values = glati + yres/2),
      time = list(unit = 'seconds since 1970-01-01',
                  position = 'beginning of hour',
                  values = as.numeric(time_out)),
      Footprint = list(unit = 'umol CO2 m-2 s-1',
                       dimensions = c('longitude', 'latitude', 'time'),
                       values = foot),
      attributes = list(crs = '+proj=longlat +ellpsWGS84',
                        crs_format = 'PROJ.4',
                        Conventions = 'CF-1.4',
                        Documentation = 'benfasoli.github.io/stilt')
    )
    saveRDS(out_custom, output)
    return(output)
  }
  
  warning('Output format must end in .nc (preferred) or .rds')
  return()
}
