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
#'   preferred), .rds (for serialized R data output), .csv (for comma separated
#'   value table), or NULL to return the footprint object for programatic use.
#'   .nc files are saved in the CF-1.4 (Climate and Forcast Metadata) convention
#'   for native use with raster::brick() and raster::raster(). rds files do not
#'   require any additional libraries and have better compression
#' @param r_run_time receptor run time as a POSIXct object. Can be NULL
#'   resulting in NULL timestamp outputs
#' @param projection proj4 string defining the map projection of the footprint
#'   netCDF output
#' @param time_integrate logical indicating whether to integrate footprint over
#'   time or retain discrete time steps
#' @param smooth_factor factor by which to linearly scale footprint smoothing;
#'   0 to disable, defaults to 1
#' @param xmn sets grid start longitude
#' @param xmx sets grid end longitude
#' @param xres resolution for longitude grid
#' @param ymn sets grid start latitude
#' @param ymx sets grid end latitude
#' @param yres resolution for latitude grid
#'
#' @import dplyr, proj4, raster
#' @export

calc_footprint <- function(p, output = NULL, r_run_time,
                           projection = '+proj=longlat',
                           smooth_factor = 1, time_integrate = F,
                           xmn, xmx, xres, ymn, ymx, yres = xres) {
  
  require(dplyr)
  require(raster)
  
  np <- length(unique(p$indx))
  
  # Interpolate particle locations during initial time steps
  times <- c(seq(0, 10, by = 0.1),
             seq(10.2, 20, by = 0.2),
             seq(20.5, 100, by = 0.5))
  time_sign <- sign(median(p$time))
  times <- times * time_sign
  
  # Preserve total field prior to split-interpolating particle positions
  aptime <- abs(p$time)
  foot_0_10_sum <- sum(p$foot[aptime <= 10], na.rm = T)
  foot_10_20_sum <- sum(p$foot[aptime > 10 & aptime <= 20], na.rm = T)
  foot_20_100_sum <- sum(p$foot[aptime > 20 & aptime <= 100], na.rm = T)

  # Split particle incluence along linear trajectory to sub-minute timescales
  p <- p %>%
    dplyr::select(indx, time, long, lati, foot) %>%
    full_join(expand.grid(time = times,
                          indx = unique(p$indx)),
              by = c('indx', 'time')) %>%
    arrange(indx, -time) %>%
    group_by(indx) %>%
    mutate(long = na_interp(long, x = time),
           lati = na_interp(lati, x = time),
           foot = na_interp(foot, x = time)) %>%
    ungroup() %>%
    na.omit() %>%
    mutate(time = round(time, 1))

  # Scale interpolated values to retain total field
  aptime <- abs(p$time)
  mi <- aptime <= 10
  p$foot[mi] <- p$foot[mi] / (sum(p$foot[mi], na.rm = T) / foot_0_10_sum)
  mi <- aptime > 10 & aptime <= 20
  p$foot[mi] <- p$foot[mi] / (sum(p$foot[mi], na.rm = T) / foot_10_20_sum)
  mi <- aptime > 20 & aptime <= 100
  p$foot[mi] <- p$foot[mi] / (sum(p$foot[mi], na.rm = T) / foot_20_100_sum)
  
  # Preserve time relative to individual particle release as rtime
  is_backward <- median(p$time) < 0
  p <- p %>%
    group_by(indx) %>%
    mutate(rtime = time - (time_sign) * min(abs(time))) %>%
    ungroup()

  # Translate x, y coordinates into desired map projection
  is_longlat <- grepl('+proj=longlat', projection, fixed = T)
  if (!is_longlat) {
    require(proj4)
    p[, c('long', 'lati')] <- project(p[, c('long', 'lati')], projection)
    grid_lims <- project(list(c(xmn, xmx), c(ymn, ymx)), projection)
    xmn <- min(grid_lims$x)
    xmx <- max(grid_lims$x)
    ymn <- min(grid_lims$y)
    ymx <- max(grid_lims$y)
  }
  
  # Set footprint grid
  glong <- head(seq(xmn, xmx, by = xres), -1)
  glati <- head(seq(ymn, ymx, by = yres), -1)
  
  # Gaussian kernel bandwidth scaling
  kernel <- p %>%
    group_by(rtime) %>%
    dplyr::summarize(varsum = var(long, na.rm = T) + var(lati, na.rm = T),
                     lati = mean(lati, na.rm = T)) %>%
    na.omit()
  di <- kernel$varsum^(1/4)
  ti <- abs(kernel$rtime/1440)^(1/2)
  grid_conv <- if (is_longlat) cos(kernel$lati * pi/180) else 1
  w <- smooth_factor * 0.06 * di * ti / grid_conv
  
  # Gaussian kernel weighting calculation
  make_gauss_kernel <- function (rs, sigma, projection) {
    # Modified from raster:::.Gauss.weight()
    require(raster)
    d <- 3 * sigma
    nx <- 1 + 2 * floor(d/rs[1])
    ny <- 1 + 2 * floor(d/rs[2])
    m <- matrix(ncol = nx, nrow = ny)
    xr <- (nx * rs[1])/2
    yr <- (ny * rs[2])/2
    r <- raster(m, xmn = -xr[1], xmx = xr[1], ymn = -yr[1], ymx = yr[1],
                crs = projection)
    p <- xyFromCell(r, 1:ncell(r))^2
    m <- 1/(2 * pi * sigma^2) * exp(-(p[, 1] + p[, 2])/(2 * sigma^2))
    m <- matrix(m, ncol = nx, nrow = ny, byrow = TRUE)
    w <- m/sum(m)
    w[is.na(w)] <- 1
    w
  }
  
  # Determine maximum kernel size
  xyres <- c(xres, yres)
  max_k <- make_gauss_kernel(xyres, max(w), projection)
  
  # Expand grid extent using maximum kernel size
  xbuf <- ncol(max_k)
  xbufh <- (xbuf - 1) / 2
  ybuf <- nrow(max_k)
  ybufh <- (ybuf - 1) / 2
  max_glong <- seq(xmn - (xbuf*xres), xmx + ((xbuf - 1)*xres), by = xres)
  max_glati <- seq(ymn - (ybuf*yres), ymx + ((ybuf - 1)*yres), by = yres)
  
  # Remove zero influence particles and positions outside of domain
  p <- p %>%
    dplyr::filter(foot > 0,
                  long >= (xmn - xbufh*xres), long < (xmx + xbufh*xres),
                  lati >= (ymn - ybufh*yres), lati < (ymx + ybufh*yres))
  if (nrow(p) == 0) return(NULL)
  
  # Pre grid particle locations
  p <- p %>%
    transmute(loi = as.integer(findInterval(long, max_glong)),
              lai = as.integer(findInterval(lati, max_glati)),
              foot = foot,
              time = time,
              rtime) %>%
    group_by(loi, lai, time, rtime) %>%
    dplyr::summarize(foot = sum(foot, na.rm = T)) %>%
    ungroup()
  
  # Dimensions in accordance with CF convention (x, y, t)
  nx <- length(max_glati)
  ny <- length(max_glong)
  grd <- matrix(0, nrow = ny, ncol = nx)
  
  # Split particle data by footprint layer
  interval <- 3600
  interval_mins <- interval / 60
  p$layer <- if (time_integrate) 0 else floor(p$time / interval_mins)
  
  layers <- sort(unique(p$layer))
  nlayers <- length(layers)
  
  # Preallocate footprint output array
  foot <- array(grd, dim = c(dim(grd), nlayers))
  for (i in 1:nlayers) {
    layer_subset <- dplyr::filter(p, layer == layers[i])
    
    rtimes <- unique(layer_subset$rtime)
    for (j in 1:length(rtimes)) {
      step <- dplyr::filter(layer_subset, rtime == rtimes[j])
      
      # Create dispersion kernel based using nearest-in-time kernel bandwidth w
      step_w <- w[find_neighbor(rtimes[j], kernel$rtime)]
      k <- make_gauss_kernel(xyres, step_w, projection)
      
      # Array dimensions
      len <- nrow(step)
      nkx <- ncol(k)
      nky <- nrow(k)
      
      # Call permute fortran subroutine to build and aggregate kernels
      out <- .Fortran('permute', ans = grd, nax = nx, nay = ny, k = k, 
                      nkx = nkx, nky = nky, len = len, lai = step$lai, 
                      loi = step$loi, foot = step$foot)
      foot[ , , i] <- foot[ , , i] + out$ans
    }
  }
  
  # Remove spatial buffer around domain used in kernel aggregation
  size <- dim(foot)
  foot <- foot[(xbuf+1):(size[1]-xbuf), (ybuf+1):(size[2]-ybuf), ] / np
  
  # Determine timestamp to use in output files
  if (time_integrate) {
    time_out <- as.numeric(r_run_time) 
  } else {
    time_out <- as.numeric(r_run_time + layers * interval)
  }
  
  # Set footprint metadata and write to file
  write_footprint(foot, output = output, glong = glong, glati = glati,
                  projection = projection, xres = xres, yres = yres,
                  time_out = time_out)
}
