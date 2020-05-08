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
                            capemin = -1,
                            cmass = 0,
                            conage = 48,
                            cpack = 1,
                            delt = 1,
                            dxf = 1,
                            dyf = 1,
                            dzf = 0.01,
                            efile = '',
                            emisshrs = 0.01,
                            frhmax = 3,
                            frhs = 1,
                            frme = 0.1,
                            frmr = 0,
                            frts = 0.1,
                            frvs = 0.1, 
                            hnf_plume = T,
                            horcoruverr = NA, 
                            horcorzierr = NA,
                            hscale = 10800,
                            ichem = 8,
                            idsp = 2,
                            initd = 0,
                            isot = -99,
                            k10m = 1,
                            kagl = 1,
                            kbls = 1,
                            kblt = 5,
                            kdef = 1,
                            khinp = 0,
                            khmax = 9999,
                            kmix0 = 150,
                            kmixd = 3,
                            kmsl = 0, 
                            kpuff = 0,
                            krand = 2,
                            krnd = 6,
                            kspl = 1,
                            kwet = 1,
                            kzmix = 0,
                            lib.loc = NULL,
                            maxdim = 1,
                            maxpar = numpar,
                            met_file_format, 
                            met_loc,
                            mgmin = 10,
                            mhrs = 9999,
                            n_hours = -24,
                            n_met_min = 1,
                            ncycl = 0,
                            ndump = 0,
                            ninit = 1,
                            nstr = 0,
                            nturb = 0,
                            nver = 0,
                            numpar = 1000,
                            outdt = 0,
                            output_wd = file.path(stilt_wd, 'out'),
                            p10f = 1,
                            pinbc = '',
                            pinpf = '',
                            poutf = '',
                            projection = '+proj=longlat',
                            qcycle = 0, 
                            r_run_time,
                            r_lati,
                            r_long,
                            r_zagl,
                            rhb = 80,
                            rht = 60, 
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
                            tout = 0, 
                            tratio = 0.75,
                            tvmix = 1,
                            varsiwant = c('time', 'indx', 'long', 'lati',
                                          'zagl', 'foot', 'mlht', 'dens',
                                          'samt', 'sigw', 'tlgr'),
                            veght = 0.5,
                            vscale = 200,
                            vscaleu = 200,
                            vscales = -1,
                            w_option = 0,
                            wbbh = 0,
                            wbwf = 0,
                            wbwr = 0,
                            wvert = FALSE,
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
    
    # Aggregate STILT/HYSPLIT namelist
    namelist <- list(
      capemin = capemin,
      cmass = cmass,
      conage = conage,
      cpack = cpack,
      delt = delt,
      dxf = dxf,
      dyf = dyf,
      dzf = dzf,
      efile = efile,
      frhmax = frhmax,
      frhs = frhs,
      frme = frme,
      frmr = frmr,
      frts = frts,
      frvs = frvs,
      hnf_plume = hnf_plume,
      hscale = hscale,
      ichem = ichem,
      idsp = idsp,
      initd = initd,
      isot = isot,
      k10m = k10m,
      kagl = kagl,
      kbls = kbls,
      kblt = kblt,
      kdef = kdef,
      khinp = khinp,
      khmax = khmax,
      kmix0 = kmix0,
      kmixd = kmixd,
      kmsl = kmsl,
      kpuff = kpuff,
      krand = krand,
      krnd = krnd,
      kspl = kspl,
      kwet = kwet,
      kzmix = kzmix,
      maxdim = maxdim,
      maxpar = maxpar,
      mgmin = mgmin,
      ncycl = ncycl,
      ndump = ndump,
      ninit = ninit,
      nstr = nstr,
      nturb = nturb,
      numpar = numpar,
      nver = nver,
      outdt = outdt,
      p10f = p10f,
      pinbc = pinbc,
      pinpf = pinpf,
      poutf = poutf,
      qcycle = qcycle,
      rhb = rhb,
      rht = rht,
      splitf = splitf,
      tkerd = tkerd,
      tkern = tkern,
      tlfrac = tlfrac,
      tout = tout,
      tratio = tratio,
      tvmix = tvmix,
      varsiwant = varsiwant,
      veght = veght,
      vscale = vscale,
      vscaleu = vscaleu,
      vscales = vscales,
      wbbh = wbbh,
      wbwf = wbwf,
      wbwr = wbwr,
      winderrtf = 0,
      wvert = wvert,
      zicontroltf = zicontroltf,
      ziscale = ziscale
    )
    
    # Creates subdirectories in out for each model run time. Each of these
    # subdirectories is populated with symbolic links to the shared datasets
    # below and a run-specific SETUP.CFG and CONTROL
    rundir_format <- paste0('%Y%m%d%H%M_', r_long, '_', r_lati, '_', 
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
      particle <- calc_trajectory(namelist, rundir, emisshrs, hnf_plume, 
                                  met_files, n_hours, output, rm_dat, timeout,
                                  w_option, z_top)
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
        particle_error <- calc_trajectory(
          merge_lists(namelist, list(winderrtf = winderrtf)), rundir, emisshrs,
          hnf_plume, met_files, n_hours, output, rm_dat, timeout, w_option, 
          z_top)
        if (is.null(particle_error)) return()
        output$particle_error <- particle_error
        output$particle_error_params <- list(siguverr = siguverr,
                                             tluverr = tluverr,
                                             zcoruverr = zcoruverr,
                                             horcoruverr = horcoruverr,
                                             sigzierr = sigzierr,
                                             tlzierr = tlzierr,
                                             horcorzierr = horcorzierr,
                                             winderrtf = winderrtf)
      }
      
      # Save output object to compressed rds file and symlink to out/particles
      saveRDS(output, output$file)
      
      link <- file.path(output_wd, 'particles', basename(output$file))
      suppressWarnings(file.symlink(output$file, link))

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
    link <- file.path(output_wd, 'footprints', basename(foot_file))
    suppressWarnings(file.symlink(foot_file, link))
    
    invisible(gc())
    return(foot)
  })
}
