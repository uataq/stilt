#' simulation_step runs STILT for the given receptor
#' @author Ben Fasoli
#'
#' Executes trajectory calculations with (and optionally without) transport
#' error and calculates kernel density derived footprint grids.
#'
#' For documentation, see https://uataq.github.io/stilt/
#'
#' @export

simulation_step <- function(before_footprint = list(function() {output}),
                            before_trajec = list(function() {output}),
                            conage = 48,
                            cpack = 1,
                            delt = 0,
                            dxf = 1,
                            dyf = 1,
                            dzf = 0.01,
                            emisshrs = 0.01,
                            frhmax = 3,
                            frhs = 1,
                            frme = 0.1,
                            frmr = 0,
                            frts = 0.1,
                            frvs = 0.1, 
                            hnf_plume = T,
                            hscale = 10800,
                            horcoruverr = NA, 
                            horcorzierr = NA,
                            ichem = 0,
                            iconvect = 0, 
                            initd = 0,
                            isot = 0,
                            kbls = 1,
                            kblt = 1,
                            kdef = 1,
                            khmax = 9999,
                            kmix0 = 250,
                            kmixd = 3,
                            kmsl = 0, 
                            kpuff = 0,
                            krnd = 6,
                            kspl = 1,
                            kzmix = 1,
                            lib.loc = NULL,
                            maxdim = 1,
                            maxpar = max(10000, numpar),
                            met_file_format, 
                            met_loc,
                            mgmin = 2000,
                            n_hours = -24,
                            n_met_min = 1,
                            ncycl = 0,
                            ndump = 0,
                            ninit = 1,
                            nturb = 0,
                            numpar = 200,
                            outdt = 0,
                            outfrac = 0.9,
                            output_wd = file.path(stilt_wd, 'out'),
                            p10f = 1,
                            projection = '+proj=longlat',
                            qcycle = 0, 
                            r_run_time,
                            r_lati,
                            r_long,
                            r_zagl,
                            random = 1, 
                            rm_dat = T,
                            run_foot = T,
                            run_trajec = T,
                            siguverr = NA,
                            sigzierr = NA, 
                            smooth_factor = 1,
                            splitf = 1,
                            stilt_wd = getwd(),
                            time_integrate = F,
                            timeout = 3600,
                            tkerd = 0.18,
                            tkern = 0.18,
                            tlfrac = 0.1,
                            tluverr = NA,
                            tlzierr = NA, 
                            tratio = 0.9,
                            tvmix = 1,
                            varsiwant = c('time', 'indx', 'long', 'lati', 'zagl', 
                                          'sigw', 'tlgr', 'zsfc', 'icdx', 'temp',
                                          'samt', 'foot', 'shtf', 'tcld', 'dmas', 
                                          'dens', 'rhfr', 'sphu', 'solw', 'lcld', 
                                          'zloc', 'dswf', 'wout', 'mlht', 'rain', 
                                          'crai', 'pres'), 
                            veght = 0.5,
                            vscale = 200,
                            w_option = 0,
                            xmn,
                            xmx,
                            xres,
                            ymn,
                            ymx, 
                            yres = xres,
                            zicontroltf = 0,
                            ziscale = 0,
                            z_top = 25000,
                            zcoruverr = NA,
                            ...) {
  try({
    setwd(stilt_wd)
    
    args <- list(...)
    
    # Ensure user specified functions reference the simulation_step environment
    if (is.list(before_footprint)) {
      before_footprint <- before_footprint[[1]]
    }
    if (is.list(before_trajec)) {
      before_trajec <- before_trajec[[1]]
    }
    environment(before_footprint) <- environment()
    environment(before_trajec) <- environment()

    # Vector style arguments passed as a list
    r_zagl <- unlist(r_zagl)
    varsiwant <- unlist(varsiwant)
    ziscale <- unlist(ziscale)

    # Validate arguments
    if (!run_trajec && !run_foot)
      stop('simulation_step(): Nothing to do, set run_trajec or run_foot to T')
    
    # Ensure dependencies are loaded for current node/process
    source(file.path(stilt_wd, 'r/dependencies.r'), local = T)
    
    # Creates subdirectories in out for each model run time. Each of these
    # subdirectories is populated with symbolic links to the shared datasets
    # below and a run-specific SETUP.CFG and CONTROL
    rundir_format <- paste0('%Y%m%d%H_', r_long, '_', r_lati, '_', 
                            ifelse(length(r_zagl) > 1, 'X', r_zagl))
    rundir  <- file.path(output_wd, 'by-id',
                         strftime(r_run_time, rundir_format, 'UTC'))
    dir.create(rundir, showWarnings = F, recursive = T)
    dir.create(file.path(output_wd, 'particles'), showWarnings = F, recursive = T)
    dir.create(file.path(output_wd, 'footprints'), showWarnings = F, recursive = T)
    message(paste('Running simulation ID:  ', basename(rundir)))
    
    # Calculate particle trajectories ------------------------------------------
    # run_trajec determines whether to try using existing trajectory files or to
    # recycle existing files
    output <- list()
    output$file <- file.path(rundir, paste0(basename(rundir), '_traj.rds'))
    
    if (run_trajec) {
      # Ensure necessary files and directory structure are established in the
      # current rundir
      if (!dir.exists(rundir)) dir.create(rundir)
      link_files <- dir(file.path(stilt_wd, 'exe'))
      suppressWarnings(
        file.symlink(file.path(stilt_wd, 'exe', link_files), rundir)
      )
      
      # Find necessary met files
      met_files <- find_met_files(r_run_time, met_file_format, n_hours, met_loc)
      if (length(met_files) < n_met_min) {
        msg <- paste('Insufficient number of meteorological files found. Check',
                     'specifications in run_stilt.r')
        warning(msg)
        cat(msg, '\n', file = file.path(rundir, 'ERROR'))
        return()
      }
      
      # Execute particle trajectory simulation, and read results into data frame
      output$receptor <- list(run_time = r_run_time,
                              lati = r_lati,
                              long = r_long,
                              zagl = r_zagl)

      # User defined function to mutate the output object
      output <- before_trajec()

      particle <- calc_trajectory(varsiwant, conage, cpack, delt, dxf, dyf, dzf,
                                  emisshrs, frhmax, frhs, frme, frmr, frts, frvs,
                                  hnf_plume, hscale, ichem, iconvect, initd, isot,
                                  ivmax, kbls, kblt, kdef, khmax, kmix0, kmixd,
                                  kmsl, kpuff, krnd, kspl, kzmix, maxdim, maxpar,
                                  met_files, mgmin, ncycl, ndump, ninit, numpar,
                                  nturb, n_hours, outdt, outfrac, output, p10f,
                                  qcycle, random, splitf, tkerd, tkern, rm_dat,
                                  timeout, tlfrac, tratio, tvmix, veght, vscale,
                                  0, w_option, zicontroltf, ziscale, z_top, 
                                  rundir)
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
                                          frmr, frts, frvs, hnf_plume, hscale,
                                          ichem, iconvect, initd, isot, ivmax,
                                          kbls, kblt, kdef, khmax, kmix0, kmixd,
                                          kmsl, kpuff, krnd, kspl, kzmix, maxdim,
                                          maxpar, met_files, mgmin, ncycl, ndump,
                                          ninit, numpar, nturb, n_hours, outdt,
                                          outfrac, output, p10f, qcycle, random,
                                          splitf, tkerd, tkern, rm_dat, timeout,
                                          tlfrac, tratio, tvmix, veght, vscale,
                                          winderrtf, w_option, zicontroltf,
                                          ziscale, z_top, rundir)
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
      saveRDS(output, output$file)
      invisible(file.symlink(output$file, file.path(output_wd, 'particles',
                                                    basename(output$file))))
      # Exit if not performing footprint calculations
      if (!run_foot) return(invisible(output$file))

    } else {
      # If user opted to recycle existing trajectory files, read in the recycled
      # file to a data frame with an adjusted timestamp and index for the
      # simulation step. If none exists, report an error and proceed
      if (!file.exists(output$file)) {
        warning('simulation_step(): No _traj.rds file found in ', rundir,
                '\n    skipping this receptor and trying the next...')
        return()
      }
      output <- readRDS(output$file)
    }
    
    # User defined function to mutate the output object
    output <- before_footprint()

    # Produce footprint --------------------------------------------------------
    # Aggregate the particle trajectory into surface influence footprints. This
    # outputs a .rds file, which can be read with readRDS() containing the
    # resultant footprint and various attributes
    foot_file <- file.path(rundir, paste0(basename(rundir), '_foot.nc'))
    foot <- calc_footprint(output$particle, output = foot_file,
                           r_run_time = r_run_time,
                           smooth_factor = smooth_factor,
                           time_integrate = time_integrate,
                           xmn = xmn, xmx = xmx, xres = xres,
                           ymn = ymn, ymx = ymx, yres = yres)
    if (is.null(foot)) {
      msg <- 'No non-zero footprint values found within the footprint domain.'
      warning(msg)
      cat(msg, '\n', file = file.path(rundir, 'ERROR'))
      return()
    }
    
    # Symlink footprint to out/footprints
    invisible(file.symlink(foot_file, file.path(output_wd, 'footprints',
                                                basename(foot_file))))

    invisible(gc())
    return(foot)
  })
}
