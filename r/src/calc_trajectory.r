#' calc_footprint generates upstream influence footprint
#' @author Ben Fasoli
#'
#' Aggregates the upstream particle trajectories into a time integrated
#' footprint, expanding particle influence using variable 2d gaussian kernels
#' with bandwidths proportional to the mean pairwise distance between all
#' particles at each time step. Requires compiled permute.so to build the
#' gaussian kernels with fortran.
#'
#' For documentation, see https://uataq.github.io/stilt/
#'
#' @export

calc_trajectory <- function(namelist,
                            rundir,
                            emisshrs,
                            hnf_plume,
                            met_files,
                            n_hours,
                            output,
                            rm_dat,
                            timeout,
                            w_option,
                            z_top) {
  
  # Enable manual rescaling of mixed layer height
  if (as.logical(namelist[['zicontroltf']])) {
    write_zicontrol(namelist[['ziscale']], file.path(rundir, 'ZICONTROL'))
  }
  
  # Write SETUP.CFG and CONTROL files to control model
  do.call(write_setup, 
          merge_lists(namelist, list(file = file.path(rundir, 'SETUP.CFG'))))
  write_control(output$receptor, emisshrs, n_hours, w_option, z_top, met_files,
                file.path(rundir, 'CONTROL'))
  
  # Simulation timeout ---------------------------------------------------------
  # Monitors time elapsed running hysplit. If elapsed time exceeds timeout
  # specified in run_stilt.r, kills hysplit and moves on to next simulation
  #
  # TODO: as of R 3.5, system() and system2() have introduced a timout arg that
  # may enable this to be depracated in the future. For now, most linux package
  # libraries are not up to date so waiting to implement edge requirements
  eval_start <- Sys.time()
  cmd <- paste('(cd', rundir, '&& (./hycs_std &> log.txt & echo $!))')
  pid <- system(cmd, intern = T)
  on.exit(tools::pskill(pid))
  repeat {
    elapsed <- as.double.difftime(Sys.time() - eval_start, units = 'secs')
    if (!pid_is_active(pid)) {
      on.exit()
      break
    } else if (elapsed > timeout) {
      msg <- paste('timeout after', elapsed, ' seconds\n')
      warning(msg)
      cat(msg, '\n', file = file.path(rundir, 'ERROR'))
      return()
    }
    Sys.sleep(0.1)
  }
  
  # Exit if running in HYSPLIT mode
  if (namelist[['ichem']] != 8) return()

  pf <- file.path(rundir, 'PARTICLE_STILT.DAT')
  if (!file.exists(pf)) {
    msg <- paste('Failed to output ', pf, ' check for errors in log.txt')
    warning(msg)
    cat(msg, '\n', file = file.path(rundir, 'ERROR'))
    return()
  }
  
  n_lines <- count_lines(pf)
  if (n_lines < 2) {
    msg <- paste(pf, 'does not contain any trajectory data.',
                 'Check for errors in log.txt')
    warning(msg)
    cat(msg, '\n', file = file.path(rundir, 'ERROR'))
    return()
  }
  
  # Read particle file, optionally remove PARTICLE.DAT in favor of compressed
  # .rds file, and return particle data frame
  p <- read_particle(file = pf, varsiwant = namelist[['varsiwant']])
  if (rm_dat) {
    system(paste('rm', pf))
    system(paste('rm', file.path(rundir, 'PARTICLE.DAT')))
  }
  
  numpar <- max(p$indx)
  
  # For column trajectories, preserve release height as xhgt
  if (length(output$receptor$zagl) > 1) {
    x_heights <- output$receptor$zagl
    px <- data.frame(indx = 1:numpar)
    px$xhgt <- rep(x_heights, each = length(px$indx) / length(x_heights))
    p <- merge(p, px, by = 'indx', sort = F)
  }
  
  # Calculate near-field dilution height based on gaussian plume width
  # approximation and recalculate footprint sensitivity for cases when the
  # plume height is less than the PBL height scaled by veght
  if (hnf_plume) 
    p <- calc_plume_dilution(p, numpar, output$receptor$zagl, namelist[['veght']])
  p
}
