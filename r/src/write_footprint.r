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
#' @param projection proj4 string defining the map projection of the footprint
#'   netCDF output
#' @param time_out time array corresponding to third dimension of foot
#' @param xres resolution of longitude grid
#' @param yres resolution of latitude grid
#'
#' @export

write_footprint <- function(foot, output, glong, glati, projection, time_out, 
                            xres, yres) {
  
  is_longlat <- grepl('+proj=longlat', projection, fixed = T)
  
  # CF projection name lookup table
  cf_proj_names <- c('aea'     = 'albers_conical_equal_area',
                     'aeqd'    = 'azimuthal_equidistant',
                     'geos'    = 'vertical_perspective',
                     'laea'    = 'labert_azimuthal_equal_area',
                     'latlon'  = 'latitude_longitude',
                     'latlong' = 'latitude_longitude', 
                     'lcc'     = 'lambert_conformal_conic',
                     'lonlat'  = 'latitude_longitude',
                     'longlat' = 'latitude_longitude',
                     'merc'    = 'mercator',
                     'ortho'   = 'orthographic',
                     'ups'     = 'polar_stereographic',
                     'stere'   = 'stereographic',
                     'tmerc'   = 'transverse_mercator')
  
  # CF projection attributes lookup table
  cf_proj_attrs <- c('h'      = 'perspective_point_height',
                     'k0'     = 'scale_factor_at_central_meridian',
                     'lat_0'  = 'latitude_of_projection_origin',
                     'lat_ts' = 'latitude_of_true_scale',
                     'lon_0'  = 'longitude_of_projection_origin',
                     'x_0'    = 'false_easting',
                     'y_0'    = 'false_northing',
                     'lat_1'  = 'standard_parallel',
                     'lat_2'  = 'standard_parallel')
  
  # Define projection name
  proj <- regmatches(projection, regexpr('\\+proj=[^ ]+', projection))
  ncdf_proj_name <- cf_proj_names[strsplit(proj, '=')[[1]][2]]
  
  # Define projection attributes
  projection_attrs <- gsub('\\+proj=[^ ]+ |\\+', '', projection) %>%
    (function(x) strsplit(x, '\\s+')[[1]]) %>%
    strsplit('=')
  ncdf_proj_attr <- list()
  for (i in 1:length(projection_attrs)) {
    key <- cf_proj_attrs[projection_attrs[[i]][1]]
    if (is.na(key)) next
    val <- projection_attrs[[i]][2]
    ncdf_proj_attr[[key]] <- ifelse(key %in% names(ncdf_proj_attr), 
                                    append(ncdf_proj_attr[[key]], val),
                                    val)
  }
  
  
  # Save footprint fo file
  if (!is.null(output) && file.exists(output))
    system(paste('rm', output))
  
  # netCDF output
  if (!is.null(output) && grepl('\\.nc$', output, ignore.case = T) &&
      'ncdf4' %in% names(sessionInfo()$otherPkgs)) {
    
    # xy dimensions in lat/lon or alternative projection
    if (is_longlat) {
      xdim <- ncdim_def('lon', 'degrees_east', glong + xres/2)
      ydim <- ncdim_def('lat', 'degrees_north', glati + yres/2)
    } else {
      xdim <- ncdim_def('x', 'm', glong + xres/2)
      ydim <- ncdim_def('y', 'm', glati + yres/2)
      pvar <- ncvar_def('projection', '', list())
    }
    tdim <- ncdim_def('time', 'seconds since 1970-01-01 00:00:00Z',
                      as.numeric(time_out))
    fvar <- ncvar_def('foot', 'ppm (umol-1 m2 s)',
                      list(xdim, ydim, tdim), -1)
    
    
    # Projection specific xy definitions
    if (is_longlat) {
      nc <- nc_create(output, list(fvar), force_v4 = T)
      ncatt_put(nc, 'lon', 'standard_name', 'longitude')
      ncatt_put(nc, 'lon', 'long_name', 'longitude at cell center')
      ncatt_put(nc, 'lat', 'standard_name', 'latitude')
      ncatt_put(nc, 'lat', 'long_name', 'latitude at cell center')
    } else {
      nc <- nc_create(output, list(fvar, pvar), force_v4 = T)
      ncatt_put(nc, 'x', 'standard_name', 'projection_x_coordinate')
      ncatt_put(nc, 'x', 'long_name', 'x coordinate of projection')
      ncatt_put(nc, 'y', 'standard_name', 'projection_y_coordinate')
      ncatt_put(nc, 'y', 'long_name', 'y coordinate of projection')
      ncatt_put(nc, 'projection', 'grid_mapping_name', ncdf_proj_name)
      ncatt_put(nc, 'foot', 'grid_mapping', 'projection')
      for (i in names(ncdf_proj_attr)) {
        ncatt_put(nc, 'projection', i, ncdf_proj_attr[[i]])
      }
    }
    
    # Insert footprint data
    ncvar_put(nc, fvar, foot)
    
    ncatt_put(nc, 'time', 'standard_name', 'time')
    ncatt_put(nc, 'time', 'long_name', 'utc time')
    ncatt_put(nc, 'time', 'calendar', 'standard')
    
    ncatt_put(nc, 'foot', 'standard_name', 'footprint')
    ncatt_put(nc, 'foot', 'long_name', 'stilt surface influence footprint')
    
    ncatt_put(nc, 0, 'crs', projection)
    ncatt_put(nc, 0, 'crs_format', 'PROJ.4')
    ncatt_put(nc, 0, 'documentation', 'github.com/uataq/stilt')
    ncatt_put(nc, 0, 'title', 'STILT Footprint')
    ncatt_put(nc, 0, 'time_created', format(Sys.time(), tz = 'UTC'))
    return(output)
  }
  
  # List output opbject
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
  
  # Alternative .csv output
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
  
  # Alternative .rds output
  if (!is.null(output) && grepl('\\.rds$', output, ignore.case = T)) {
    saveRDS(out_custom, output)
    return(output)
  }
  
  # No file output - return list object
  invisible(out_custom)
}
