#' simulation_step runs STILT for the given timestep
#' @author Ben Fasoli
#'
#' @export

simulation_step <- function(X, rm_dat = T, stilt_wd = getwd(), lib.loc = NULL,
                            slurm = F, met_file_format, met_loc,
                            delt = 0, iconvect = 0, isot = 0,
                            mgmin = 2000, n_cores_grid = 1, n_hours = -72,
                            ndump = 0, nturb = 0, numpar = 100, outdt = 0,
                            outfrac = 0.9, run_trajec = T, random = 1,
                            r_run_time, r_lati, r_long, r_zagl,
                            tlfrac = 0.1, tratio = 0.9,
                            varsiwant = NULL, veght = 0.5, w_option = 0,
                            winderrtf = 0, zicontroltf = 0, z_top = 25000,
                            xmn = -180, xmx = 180, xres = 0.1,
                            ymn = -90, ymx = 90, yres = xres) {

  # If using lapply or parLapply, receptors are passed as vectors and need to
  # be subsetted for the specific simulation index
  if (!slurm) {
    r_run_time <- r_run_time[X]
    r_lati <- r_lati[X]
    r_long <- r_long[X]
    r_zagl <- r_zagl[X]
  } else {
    # If job submission via sbatch, dependencies need to be loaded for each
    # unique process. Receptor information is passed as single values for each
    # simulation step
    source(file.path(stilt_wd,'r/dependencies.r'), local = T)
  }

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
  rundir  <- file.path(file.path(stilt_wd, 'out'),
                       strftime(r_run_time, tz = 'UTC',
                                format = paste0(X, '_%Y%m%d%H_', r_long, '_',
                                                r_lati, '_', r_zagl)))

  # Generate PARTICLE.DAT ------------------------------------------------------
  # run_trajec determines whether to try using existing PARTICLE.DAT files or to
  # recycle existing files
  output <- list()
  output$runtime <- r_run_time
  output$file <- file.path(rundir, 'STILT_OUTPUT.rds')
  output$receptor <- data_frame(run_time = r_run_time,
    lati = r_lati,
    long = r_long,
    zagl = r_zagl)

  if (run_trajec) {
    dir.create(rundir)
    link_files <- c('ASCDATA.CFG', 'CONC.CFG', 'hymodelc',
                    'LANDUSE.ASC', 'ROUGLEN.ASC')
    file.symlink(file.path(file.path(stilt_wd, 'exe'), link_files), rundir)

    # Find met files necessary for simulation ------------------------------------
    met_files <- find_met_files(r_run_time, met_file_format, n_hours, met_loc)
    if (length(met_files) < 1) return(NULL)

    write_setup(numpar, delt, tratio, isot, tlfrac, ndump, random, outdt, nturb,
                veght, outfrac, iconvect, winderrtf, zicontroltf, mgmin,
                varsiwant, file.path(rundir, 'SETUP.CFG'))
    write_control(output$receptor, n_hours, w_option, z_top, met_files,
                  file.path(rundir, 'CONTROL'))
    write_runhymodelc(file.path(rundir, 'runhymodelc.bat')) %>%
      paste('bash', .) %>%
      system()

    pf <- file.path(rundir, 'PARTICLE.DAT')

    n_lines <- uataq::count_lines(pf)
    if (n_lines < 2) return()

    particle <- read_particle(file = pf,
                              varsiwant = varsiwant)
    if (rm_dat) system(paste('rm', pf))
  } else {
    # If user opted to recycle existing PARTICLE.DAT files, read in the recycled
    # file to a data_frame with an adjusted timestamp and index for the
    # simulation step. If none exists, report an error and try the next timestep
    if (!file.exists(output$file)) {
      warning(paste(
        'simulation_step(): No STILT_OUTPUT.rds file found in', rundir, '\n',
        'skipping this timestep and trying the next...'))
      return()
    }
    particle <- readRDS(output$file)$particle
  }


  # Produce FOOTPRINT.rds ------------------------------------------------------
  # Aggregate the particle trajectory into surface influence footprints. This
  # outputs a FOOTPRINT.rds file, which can be read with readRDS() containing
  # the resultant footprint and various attributes
  footprint <- list()
  footprint$data <- calc_footprint(particle, output = NULL,
                 n_cores_grid = n_cores_grid,
                 xmn = xmn, xmx = xmx, xres = xres,
                 ymn = ymn, ymx = ymx, yres = yres)

  output$footprint <- footprint
  output$particle <- particle
  saveRDS(output, output$file)
  return(T)
}
