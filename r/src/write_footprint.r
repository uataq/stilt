#' write_control writes out a namelist file to control the model behavior
#' @author Ben Fasoli
#'
#' @param foot footprint array calculated by calc_footprint()
#' @param output filename argument. Must end with .nc (for ncdf output,
#'   preferred), .rds (for serialized R data output), .csv (for comma separated
#'   value table), or NULL to return the footprint object for programatic use.
#'   .nc files are saved in the CF-1.4 (Climate and Forcast Metadata) convention
#'   for native use with raster::brick() and raster::raster(). rds files do not
#'   require any additional libraries and have better compression
#' @param glong longitude array corresponding to first dimension of foot
#' @param glati latitude array corresponding to second dimension of foot
#' @param xres resolution of longitude grid
#' @param yres resolution of latitude grid
#' @param time_out time array corresponding to third dimension of foot
#'
#' @export

write_footprint <- function(foot, output, glong, glati, xres, yres) {
  
  # Save footprint fo file
  if (!is.null(output) && file.exists(output))
    system(paste('rm', output))
  
  if (!is.null(output) && grepl('\\.nc$', output, ignore.case = T) &&
      'ncdf4' %in% names(sessionInfo()$otherPkgs)) {
    xdim <- ncdim_def('lon', 'degrees_east', glong + xres/2)
    ydim <- ncdim_def('lat', 'degrees_north', glati + yres/2)
    tdim <- ncdim_def('time', 'seconds since 1970-01-01 00:00:00',
                      as.numeric(time_out))
    fvar <- ncvar_def('foot', 'ppm (umol-1 m2 s)',
                      list(xdim, ydim, tdim), -1)
    
    nc <- nc_create(output, fvar)
    ncvar_put(nc, fvar, foot)
    
    ncatt_put(nc, 'lon', 'standard_name', 'longitude')
    ncatt_put(nc, 'lon', 'long_name', 'cell center longitude')
    ncatt_put(nc, 'lat', 'standard_name', 'latitude')
    ncatt_put(nc, 'lat', 'long_name', 'cell center latitude')
    ncatt_put(nc, 'time', 'standard_name', 'time')
    ncatt_put(nc, 'time', 'long_name', 'time')
    ncatt_put(nc, 'time', 'calendar', 'standard')
    ncatt_put(nc, 'time', 'timezone', 'UTC')
    ncatt_put(nc, 'foot', 'standard_name', 'footprint')
    ncatt_put(nc, 'foot', 'long_name', 'footprint')
    ncatt_put(nc, 0, 'crs', '+proj=longlat +ellpsWGS84')
    ncatt_put(nc, 0, 'crs_format', 'PROJ.4')
    ncatt_put(nc, 0, 'conventions', 'CF-1.4')
    ncatt_put(nc, 0, 'documentation', 'github.com/uataq/stilt')
    ncatt_put(nc, 0, 'title', 'STILT Footprint Output')
    return(output)
  }
  
  out_custom <- list(
    lon = list(unit = 'degrees_east',
               values = glong),
    lat = list(unit = 'degrees_north',
               values = glati),
    time = list(unit = 'seconds since 1970-01-01 00:00:00',
                position = 'beginning of hour',
                values = as.numeric(time_out)),
    footprint = list(unit = 'ppm (umol-1 m2 s)',
                     dimensions = c('longitude', 'latitude', 'time'),
                     values = foot),
    attributes = list(crs = '+proj=longlat +ellpsWGS84',
                      crs_format = 'PROJ.4',
                      conventions = 'CF-1.4',
                      documentation = 'uataq.github.io/stilt',
                      title = 'STILT Footprint Output')
  )
  
  if (!is.null(output) && grepl('\\.csv$', output, ignore.case = T)) {
    csv <- data_frame(expand.grid(longitude = glong, latitude  = glati),
                      c(foot)) %>%
      filter(foot > 0)
    write('STILT Footprint. For documentation, see uataq.github.io/stilt',
          file = output)
    write('latitude and longitude positions indicate cell center',
          file = output, append = T)
    write.table(csv, append = T, quote = F, sep = ',', row.names = F)
    return(output)
  }
  
  if (!is.null(output) && grepl('\\.rds$', output, ignore.case = T)) {
    saveRDS(out_custom, output)
    return(output)
  }
  
  invisible(out_custom)
}
