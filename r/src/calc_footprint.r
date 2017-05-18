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
#' @param output filename argument passed to \code{saveRDS}
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
                           time_integrate = F,
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
    mean(bootstrap(df, foo, size = 50, iter = 4), na.rm = T)
  }
  
  pd <- xyzt %>%
    group_by(time) %>%
    summarize(dist = calc_dist(long, lati),
              lati = mean(lati, na.rm = T))
  
  # Remove zero influence particles and positions outside of domain
  xyzt <- i %>%
    filter(foot > 0, long >= xmn, long <= xmx, lati >= ymn, lati <= ymx)
  
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
  calc_bandwidth <- function(dist, lati, xyres) {
    dist / (40 * cos(lati * pi/180)) + max(xyres) / 8
  }
  
  xyres <- c(xres, yres)
  
  # Determine maximum kernel size
  max_k <- make_gauss_kernel(xyres, calc_bandwidth(max(pd$dist), min(pd$lati), xyres))
  xbuf <- (ncol(max_k) - 1) / 2
  ybuf <- (nrow(max_k) - 1) / 2
  
  max_glong <- seq(xmn - (xbuf*xres), xmx + ((xbuf - 1)*xres), by = xres)
  max_glati <- seq(ymn - (ybuf*yres), ymx + ((ybuf - 1)*yres), by = yres)
  
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
    step <- xyzt %>%
      filter(time == x)
    
    # Dispersion kernel
    idx <- pd$time == x
    d <- calc_bandwidth(pd$dist[idx], pd$lati[idx], xyres)
    k <- make_gauss_kernel(xyres, d)
    
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
  
  size <- dim(foot)
  
  # Flip matrix vertically to align with lati domain and fill by row for
  # compatibility with raster
  for (i in 1:size[3]) {
    foot[size[1]:1, ,i] <- matrix(foot[ , ,i], byrow = T,
                                  nrow = size[1], ncol = size[2])
  }
  
  # Subset domain to extent specified with xmn/xmx/ymn/ymx
  # Use tolerance range to eliminate issues with floating point comparisons
  tol <- 1e-8
  foot <- foot[rev(max_glati < ymx-tol & max_glati >= ymn-tol),
               max_glong < xmx-tol & max_glong >= xmn-tol, ] / np
  
  # Time integrate footprint by summing across 3rd dimension
  if (time_integrate) {
    foot <- apply(foot, c(1, 2), sum)
    time_out <- list(unit = 'Time Integrated',
                     data = NULL)
    foot_dims <- c('lati', 'long')
  } else {
    size <- dim(foot)
    hid <- floor(pd$time / 60)
    foot_hour <- array(0, dim = c(size[1], size[2], length(unique(hid))))
    for (i in 1:length(unique(hid))) {
      foot_hour[ , ,i] <- apply(foot[ , ,hid == hid[i]], c(1, 2), sum)
    }
    foot <- foot_hour
    time_out <- list(unit = 'POSIXct Time',
                     data = r_run_time + unique(hid) * 3600)
    foot_dims <- c('lati', 'long', 'time')
  }
  
  # Format output object
  foot <- list(
    lati = list(unit = 'Degrees',
                res = yres,
                extent = c(ymn, ymx),
                data = rev(glati)),
    long = list(unit = 'Degrees',
                res = xres,
                extent = c(xmn, xmx),
                data = glong),
    time = time_out,
    foot = list(unit = 'ppm umol-1 m2 s',
                dims = foot_dims, data = foot)
  )
  
  if (!is.null(output))
    saveRDS(foot, output)
  
  return(foot)
}
