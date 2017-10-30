#' simulation_step runs STILT for the given timestep
#' @author Ben Fasoli
#'
#' @export

simulation_step <- function(X, rm_dat = T, stilt_wd = getwd(), lib.loc = NULL,
                            met_file_format, met_loc, delt = 0,
                            iconvect = 0, isot = 0, khmax = 9999, kmix0 = 250,
                            kmixd = 3, krnd = 6, mgmin = 2000, n_hours = -72,
                            n_met_min = 1, ndump = 0, nturb = 0, numpar = 100,
                            outdt = 0, outfrac = 0.9, run_trajec = T, random = 1,
                            r_run_time, r_lati, r_long, r_zagl,
                            smooth_factor = 1, time_integrate = F, timeout = 3600,
                            tlfrac = 0.1, tratio = 0.9, varsiwant = NULL,
                            veght = 0.5, w_option = 0, winderrtf = 0,
                            zicontroltf = 0, z_top = 25000,
                            xmn = -180, xmx = 180, xres = 0.1,
                            ymn = -90, ymx = 90, yres = xres) {
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
    # subdirectories is populated with symbolic links to the shared datasets below
    # and a run-specific SETUP.CFG and CONTROL
    rundir  <- file.path(file.path(stilt_wd, 'out', 'by-id'),
                         strftime(r_run_time, tz = 'UTC',
                                  format = paste0('%Y%m%d%H_', r_long, '_',
                                                  r_lati, '_', r_zagl)))
    uataq::br()
    message(paste('Running simulation ID:  ', basename(rundir)))

    # Calculate particle trajectories --------------------------------------------
    # run_trajec determines whether to try using existing trajectory files or to
    # recycle existing files
    output <- list()
    output$runtime <- r_run_time
    output$file <- file.path(rundir, paste0(basename(rundir), '_traj.rds'))
    output$receptor <- data_frame(run_time = r_run_time,
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

      # Save workspace to recreate simulation in settings.RData
      save.image(file = file.path(rundir, 'settings.RData'))

      # Find necessary met files
      met_files <- find_met_files(r_run_time, met_file_format, n_hours, met_loc)
      if (length(met_files) < n_met_min) {
        warning('Insufficient amount of meteorological data found...')
        cat('Insufficient amount of meteorological data found. Check ',
            'specifications in run_stilt.r\n', file = file.path(rundir, 'ERROR'))
        return()
      }

      # Write SETUP.CFG, CONTROL, and runhymodelc.sh files to control model
      write_setup(varsiwant, delt, iconvect, isot, khmax, kmix0, kmixd, krnd,
                  mgmin, ndump, numpar, nturb, outdt, outfrac, random, tlfrac,
                  tratio, veght, winderrtf, zicontroltf,
                  file.path(rundir, 'SETUP.CFG'))
      write_control(output$receptor, n_hours, w_option, z_top, met_files,
                    file.path(rundir, 'CONTROL'))
      sh <- write_runhymodelc(file.path(rundir, 'runhymodelc.sh'))

      of <- file.path(rundir, 'hymodelc.out')

      # Simulation timeout -------------------------------------------------------
      # Monitors time elapsed running hymodelc. If elapsed time exceeds timeout
      # specified in run_stilt.r, kills hymodelc and moves on to next simulation
      elapsed <- 0
      eval_start <- Sys.time()
      pid <- system(paste('bash', sh), intern = T)
      on.exit(tools::pskill(pid))
      repeat {
        elapsed <- as.double.difftime(Sys.time() - eval_start, units = 'secs')
        if (file.exists(of) && length(readLines(of)) > 0) {
          on.exit()
          break
        } else if (elapsed > timeout) {
          warning(basename(rundir), ' timeout. Killing hymodelc pid ', pid, '\n')
          cat('hymodelc timeout after ', elapsed, ' seconds\n',
              file = file.path(rundir, 'ERROR'))
          return()
        }
        Sys.sleep(1)
      }

      # Error check hymodelc output
      pf <- file.path(rundir, 'PARTICLE.DAT')
      if (!file.exists(pf)) {
        warning('Failed to output PARTICLE.DAT in ', basename(rundir))
        cat('No PARTICLE.DAT found. Check for errors in hymodelc.out\n',
            file = file.path(rundir, 'ERROR'))
        return()
      }

      n_lines <- uataq::count_lines(pf)
      if (n_lines < 2) {
        warning('No trajectory data found in ', pf)
        cat('PARTICLE.DAT does not contain any trajectory data. Check for ',
            'errors in hymodelc.out\n', file = file.path(rundir, 'ERROR'))
        return()
      }

      # Convert PARTICLE.DAT data into compressed .rds data_frame
      particle <- read_particle(file = pf, varsiwant = varsiwant)
      output$particle <- particle
      saveRDS(output, output$file)

      link(output$file,
           file.path(stilt_wd, 'out', 'particles', basename(output$file)))

      if (rm_dat) system(paste('rm', pf))
    } else {
      # If user opted to recycle existing trajectory files, read in the recycled
      # file to a data_frame with an adjusted timestamp and index for the
      # simulation step. If none exists, report an error and try the next timestep
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
    particle <- calc_plume_dilution(particle, numpar, r_zagl, veght)

    # Produce footprint ----------------------------------------------------------
    # Aggregate the particle trajectory into surface influence footprints. This
    # outputs a .rds file, which can be read with readRDS() containing the
    # resultant footprint and various attributes
    foot_file <- file.path(rundir, paste0(basename(rundir), '_foot.nc'))
    foot <- calc_footprint(particle, output = foot_file,
                           r_run_time = r_run_time, smooth_factor = smooth_factor,
                           time_integrate = time_integrate,
                           xmn = xmn, xmx = xmx, xres = xres,
                           ymn = ymn, ymx = ymx, yres = yres)
    if (is.null(foot)) {
      warning('No non-zero footprint values found within the footprint domain.')
      cat('No non-zero footprint values found within the footprint domain.\n',
          file = file.path(rundir, 'ERROR'))
      return()
    } else {
      link(foot_file,
           file.path(stilt_wd, 'out', 'footprints', basename(foot_file)))
    }
    invisible(gc())
    return(foot)
  })
}
