#' simulation_step runs STILT for the given timestep
#' @author Ben Fasoli
#'
#' Executes trajectory calculations with (and optionally without) transport
#' error and calculates kernel density derived footprint grids.
#'
#' For documentation, see https://uataq.github.io/stilt/
#'
#' @export

simulation_step <- function(X, rm_dat = T, stilt_wd = getwd(), lib.loc = NULL,
                            conage = 48, cpack = 1, delt = 0, dxf = 1, dyf = 1,
                            dzf = 0.01, emisshrs = 0.01, frhmax = 3, frhs = 1,
                            frme = 0.1, frmr = 0, frts = 0.1, frvs = 0.1, hnf_plume = T,
                            hscale = 10800, horcoruverr = NA, horcorzierr = NA,
                            ichem = 0, iconvect = 0, initd = 0, isot = 0,
                            kbls = 1, kblt = 1, kdef = 1, khmax = 9999,
                            kmix0 = 250, kmixd = 3, kmsl = 0, kpuff = 0,
                            krnd = 6, kspl = 1, kzmix = 1, maxdim = 1,
                            maxpar = 10000, met_file_format, met_loc,
                            mgmin = 2000, n_hours = -24, n_met_min = 1,
                            ncycl = 0, ndump = 0, ninit = 1, nturb = 0,
                            numpar = 200, outdt = 0, outfrac = 0.9,
                            output_wd = file.path(stilt_wd, 'out'), p10f = 1,
                            projection = '+proj=longlat', qcycle = 0, r_run_time,
                            r_lati, r_long, r_zagl, random = 1, run_trajec = T,
                            siguverr = NA, sigzierr = NA, smooth_factor = 1,
                            splitf = 1, time_integrate = F, timeout = 3600,
                            tkerd = 0.18, tkern = 0.18, tlfrac = 0.1,
                            tluverr = NA, tlzierr = NA, tratio = 0.9, tvmix = 1,
                            varsiwant = NULL, veght = 0.5, vscale = 200,
                            w_option = 0, xmn = -180, xmx = 180, xres = 0.1,
                            ymn = -90, ymx = 90, yres = xres, zicontroltf = 0,
                            z_top = 25000, zcoruverr = NA) {
  try({
    # If using lapply or parLapply, receptors are passed as vectors and need to
    # be subsetted for the specific simulation index
    if (length(r_run_time) > 1) {
      r_run_time <- r_run_time[X]
      r_lati <- r_lati[X]
      r_long <- r_long[X]
      r_zagl <- r_zagl[X]
    }
    
    # Column trajectories use r_zagl passed as a list of values for lapply and 
    # parLapply but a vector in slurm_apply
    r_zagl <- unlist(r_zagl)
    
    # Ensure dependencies are loaded for current node/process
    source(file.path(stilt_wd, 'r/dependencies.r'), local = T)
    
    if (is.null(varsiwant)) {
      varsiwant <- c('time', 'indx', 'long', 'lati', 'zagl', 'sigw', 'tlgr',
                     'zsfc', 'icdx', 'temp', 'samt', 'foot', 'shtf', 'tcld',
                     'dmas', 'dens', 'rhfr', 'sphu', 'solw', 'lcld', 'zloc',
                     'dswf', 'wout', 'mlht', 'rain', 'crai')
    } else if (any(grepl('/', varsiwant))) {
      varsiwant <- unlist(strsplit(varsiwant, '/', fixed = T))
    }
    
    # Creates subdirectories in out for each model run time. Each of these
    # subdirectories is populated with symbolic links to the shared datasets
    # below and a run-specific SETUP.CFG and CONTROL
    rundir_format <- paste0('%Y%m%d%H_', r_long, '_', r_lati, '_', 
                            ifelse(length(r_zagl) > 1, 'X', r_zagl))
    rundir  <- file.path(output_wd, 'by-id',
                         strftime(r_run_time, rundir_format, 'UTC'))
    uataq::br()
    message(paste('Running simulation ID:  ', basename(rundir)))
    
    # Calculate particle trajectories ------------------------------------------
    # run_trajec determines whether to try using existing trajectory files or to
    # recycle existing files
    output <- list()
    output$file <- file.path(rundir, paste0(basename(rundir), '_traj.rds'))
    
    if (run_trajec) {
      # Ensure necessary files and directory structure are established in the
      # current rundir
      dir.create(rundir)
      link_files <- dir(file.path(stilt_wd, 'exe'))
      file.symlink(file.path(stilt_wd, 'exe', link_files), rundir)
      
      # Find necessary met files
      met_files <- find_met_files(r_run_time, met_file_format, n_hours, met_loc)
      if (length(met_files) < n_met_min) {
        warning('Insufficient amount of meteorological data found...')
        cat('Insufficient amount of meteorological data found. Check ',
            'specifications in run_stilt.r\n',
            file = file.path(rundir, 'ERROR'))
        return()
      }
      
      # Execute particle trajectory simulation, and read results into data frame
      output$receptor <- list(run_time = r_run_time,
                              lati = r_lati,
                              long = r_long,
                              zagl = r_zagl)
      particle <- calc_trajectory(varsiwant, conage, cpack, delt, dxf, dyf, dzf,
                                  emisshrs, frhmax, frhs, frme, frmr, frts, frvs,
                                  hscale, ichem, iconvect, initd, isot, ivmax,
                                  kbls, kblt, kdef, khmax, kmix0, kmixd, kmsl,
                                  kpuff, krnd, kspl, kzmix, maxdim, maxpar,
                                  met_files, mgmin, ncycl, ndump, ninit, numpar,
                                  nturb, n_hours, outdt, outfrac, output, p10f,
                                  qcycle, random, splitf, tkerd, tkern, rm_dat,
                                  timeout, tlfrac, tratio, tvmix, veght, vscale,
                                  0, w_option, zicontroltf, z_top, rundir)
      if (is.null(particle)) return()
      
      # Bundle trajectory configuration metadata with trajectory informtation
      output$particle <- particle
      output$params <- read_config(file = file.path(rundir, 'CONC.CFG'))
      
      # Optionally execute second trajectory simulations to quantify transport
      # error using parameterized correlation length and time scales
      xyerr <- write_winderr(siguverr, tluverr, zcoruverr, horcoruverr,
                             file.path(rundir, 'WINDERR'))
      zerr <- write_zierr(sigzierr, tlzierr, horcorzierr,
                          file = file.path(rundir, 'ZIERR'))
      winderrtf <- (!is.null(xyerr)) + 2 * !is.null(zerr)
      if (winderrtf > 0) {
        particle_error <- calc_trajectory(varsiwant, conage, cpack, delt, dxf,
                                          dyf, dzf, emisshrs, frhmax, frhs, frme,
                                          frmr, frts, frvs, hscale, ichem,
                                          iconvect, initd, isot, ivmax, kbls,
                                          kblt, kdef, khmax, kmix0, kmixd, kmsl,
                                          kpuff, krnd, kspl, kzmix, maxdim,
                                          maxpar, met_files, mgmin, ncycl, ndump,
                                          ninit, numpar, nturb, n_hours, outdt,
                                          outfrac, output, p10f, qcycle, random,
                                          splitf, tkerd, tkern, rm_dat, timeout,
                                          tlfrac, tratio, tvmix, veght, vscale,
                                          winderrtf, w_option, zicontroltf,
                                          z_top, rundir)
        if (is.null(particle_error)) return()
        output$particle_error <- particle_error
        output$particle_error_params <- list(siguverr = siguverr,
                                             tluverr = tluverr,
                                             zcoruverr = zcoruverr,
                                             horcoruverr = horcoruverr,
                                             sigzierr = sigzierr,
                                             tlzierr = tlzierr,
                                             horcorzierr = horcorzierr)
      }
      
      # Save output object to compressed rds file and symlink to out/particles
      # directory for convenience
      saveRDS(output, output$file)
      file.symlink(output$file, file.path(output_wd, 'particles',
                                          basename(output$file))) %>%
        invisible()
      
    } else {
      # If user opted to recycle existing trajectory files, read in the recycled
      # file to a data_frame with an adjusted timestamp and index for the
      # simulation step. If none exists, report an error and proceed
      if (!file.exists(output$file)) {
        warning('simulation_step(): No _traj.rds file found in ', rundir,
                '\n    skipping this timestep and trying the next...')
        return()
      }
      particle <- readRDS(output$file)$particle
    }
    
    # Calculate near-field dilution height based on gaussian plume width
    # approximation and recalculate footprint sensitivity for cases when the
    # plume height is less than the PBL height scaled by veght
    if (hnf_plume)
      particle <- calc_plume_dilution(particle, numpar, r_zagl, veght)
    
    # Produce footprint --------------------------------------------------------
    # Aggregate the particle trajectory into surface influence footprints. This
    # outputs a .rds file, which can be read with readRDS() containing the
    # resultant footprint and various attributes
    foot_file <- file.path(rundir, paste0(basename(rundir), '_foot.nc'))
    foot <- calc_footprint(particle, output = foot_file,
                           r_run_time = r_run_time,
                           smooth_factor = smooth_factor,
                           time_integrate = time_integrate,
                           xmn = xmn, xmx = xmx, xres = xres,
                           ymn = ymn, ymx = ymx, yres = yres)
    if (is.null(foot)) {
      warning('No non-zero footprint values found within the footprint domain.')
      cat('No non-zero footprint values found within the footprint domain.\n',
          file = file.path(rundir, 'ERROR'))
      return()
    } else {
      file.symlink(foot_file, file.path(output_wd, 'footprints',
                                        basename(foot_file))) %>%
        invisible()
    }
    invisible(gc())
    return(foot)
  })
}
