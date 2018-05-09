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
                            delt = 0, emisshrs = 0.01, frmr = 0, hnf_plume = T,
                            horcoruverr = NA, horcorzierr = NA, iconvect = 0,
                            initd = 0, isot = 0, khmax = 9999, kmix0 = 250,
                            kmixd = 3, krnd = 6, met_file_format, met_loc,
                            mgmin = 2000, n_hours = -24, n_met_min = 1,
                            ndump = 0, nturb = 0, numpar = 200, outdt = 0,
                            outfrac = 0.9, projection = '+proj=longlat',
                            random = 1, run_trajec = T, r_run_time, r_lati,
                            r_long, r_zagl, siguverr = NA, sigzierr = NA,
                            smooth_factor = 1, time_integrate = F,
                            timeout = 3600, tlfrac = 0.1, tluverr = NA,
                            tlzierr = NA, tratio = 0.9, varsiwant = NULL,
                            veght = 0.5, w_option = 0, xmn = -180, xmx = 180,
                            xres = 0.1, ymn = -90, ymx = 90, yres = xres,
                            zicontroltf = 0, z_top = 25000,  zcoruverr = NA) {
  try({
    # If using lapply or parLapply, receptors are passed as vectors and need to
    # be subsetted for the specific simulation index
    if (length(r_run_time) > 1) {
      r_run_time <- r_run_time[X]
      r_lati <- r_lati[X]
      r_long <- r_long[X]
      r_zagl <- r_zagl[X]
    }
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
    rundir  <- file.path(file.path(stilt_wd, 'out', 'by-id'),
                         strftime(r_run_time, tz = 'UTC',
                                  format = paste0('%Y%m%d%H_', r_long, '_',
                                                  r_lati, '_', r_zagl)))
    uataq::br()
    message(paste('Running simulation ID:  ', basename(rundir)))

    # Calculate particle trajectories ------------------------------------------
    # run_trajec determines whether to try using existing trajectory files or to
    # recycle existing files
    output <- list()
    output$runtime <- r_run_time
    output$file <- file.path(rundir, paste0(basename(rundir), '_traj.rds'))
    output$receptor <- list(run_time = r_run_time,
                            lati = r_lati,
                            long = r_long,
                            zagl = r_zagl)

    if (run_trajec) {
      # Ensure necessary files and directory structure are established in the
      # current rundir
      dir.create(rundir)
      link_files <- c('ASCDATA.CFG', 'CONC.CFG', 'hymodelc',
                      'LANDUSE.ASC', 'ROUGLEN.ASC')
      file.symlink(file.path(file.path(stilt_wd, 'exe'), link_files), rundir)

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
      particle <- calc_trajectory(varsiwant, delt, emisshrs, frmr, iconvect,
                                  initd, isot, khmax, kmix0, kmixd, krnd,
                                  met_files, mgmin, ndump, numpar, nturb,
                                  n_hours, outdt, outfrac, output, random,
                                  rm_dat, timeout, tlfrac, tratio, veght, 0,
                                  w_option, zicontroltf, z_top, rundir)
      if (is.null(particle)) return()
      output$particle <- particle

      # Optionally execute second trajectory simulations to quantify transport
      # error using parameterized correlation length and time scales
      xyerr <- write_winderr(siguverr, tluverr, zcoruverr, horcoruverr,
                             file.path(rundir, 'WINDERR'))
      zerr <- write_zierr(sigzierr, tlzierr, horcorzierr,
                          file = file.path(rundir, 'ZIERR'))
      winderrtf <- (!is.null(xyerr)) + 2 * !is.null(zerr)
      if (winderrtf > 0) {
        particle_error <- calc_trajectory(varsiwant, delt, emisshrs, iconvect,
                                          isot, khmax, kmix0, kmixd, krnd,
                                          met_files, mgmin, ndump, numpar,
                                          nturb, n_hours, outdt, outfrac,
                                          output, random, rm_dat, tlfrac,
                                          tratio, veght, winderrtf, w_option,
                                          zicontroltf, z_top, rundir)
        if (is.null(particle_error)) return()
        output$particle_error_params <- list(siguverr = siguverr,
                                             tluverr = tluverr,
                                             zcoruverr = zcoruverr,
                                             horcoruverr = horcoruverr,
                                             sigzierr = sigzierr,
                                             tlzierr = tlzierr,
                                             horcorzierr = horcorzierr)
        output$particle_error <- particle_error
      }

      # Save output object to compressed rds file and symlink to out/particles
      # directory for convenience
      saveRDS(output, output$file)
      file.symlink(output$file, file.path(stilt_wd, 'out', 'particles',
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
      file.symlink(foot_file, file.path(stilt_wd, 'out', 'footprints',
                                        basename(foot_file))) %>%
        invisible()
    }
    invisible(gc())
    return(foot)
  })
}
