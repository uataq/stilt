#' calc_met_subgrids generates necessary meteorological subdomains
#'
#' \code{na_interp} linearly interpolates NA values found in vector y with
#'   respect to index x (e.g. timestamp).
#'
#' @param met_files vector of files returned by find_met_files
#' @param output_met directory to contain subgrids
#' @param exe path to exe directory containing xtrct_grid binary
#' @param projection proj4 string defining the map projection of the grid extent
#' @param xmn minimum x coordinate of subdomain
#' @param xmx maximum x coordinate of subdomain
#' @param ymn minimum y coordinate of subdomain
#' @param ymx maximum y coordinate of subdomain
#' @param levels number of vertical levels to include, defaults to NA which
#'   returns all vertical levels available in the input file
#' 
#' @export

calc_met_subgrids <- function (met_files,
                               output_met,
                               exe,
                               projection,
                               xmn,
                               xmx,
                               ymn,
                               ymx,
                               levels = NA,
                               met_subgrid_buffer = 0.1) {
  dir.create(output_met, showWarnings = F)
  
  timeout <- 60
  len <- length(met_files)
  eval_start <- Sys.time()
  i <- 1
  while (i <= len) {
    f <- met_files[i]
    
    # skip file if already exists
    subgrid <- file.path(output_met, basename(f))
    if (file.exists(subgrid)) {
      i <- i + 1
      eval_start <- Sys.time()
      next
    }
    
    # throw error if xtrct_grid duration exceeds timeout
    elapsed <- as.double.difftime(Sys.time() - eval_start, units = 'secs')
    if (elapsed > timeout) {
      msg <- paste('Failed to extract subgrid in', lock)
      warning(msg)
      return()
    }
    
    # skip file if another process already acquired lock
    lock <- paste0(subgrid, '.lock')
    if (dir.exists(lock)) {
      if (i == len) {
        Sys.sleep(0.1)
      }
      next
    }
    
    # attempt to acquire lock
    if (!dir.create(lock, showWarnings = F)) {
      next
    }
    
    # calculate subgrid and remove lock
    link_files(exe, lock)
    
    is_longlat <- grepl('+proj=longlat', projection, fixed = T)
    if (!is_longlat) {
      require(proj4)
      grid_lims <- project(list(c(xmn, xmx), c(ymn, ymx)), projection)
      xmn <- min(grid_lims$x)
      xmx <- max(grid_lims$x)
      ymn <- min(grid_lims$y)
      ymx <- max(grid_lims$y)
    }
    xrng <- diff(range(xmn, xmx))
    yrng <- diff(range(ymn, ymx))
    xbuf <- xrng * met_subgrid_buffer
    ybuf <- yrng * met_subgrid_buffer
    
    success <- xtrct_grid(f,
                          subgrid, 
                          lock,
                          xmn = max(-180, xmn - xbuf),
                          xmx = min(180, xmx + xbuf),
                          ymn = max(-90, ymn - ybuf),
                          ymx = min(90, ymx + ybuf),
                          levels = NA)
    
    if (!success) {
      msg <- paste('Failed to extract subgrid in', lock)
      warning(msg)
      return()
    }
    
    unlink(lock, recursive = T)
  }
}