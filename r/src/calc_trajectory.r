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
  # Monitors time elapsed running hycs_std If elapsed time exceeds timeout
  # specified in run_stilt.r, kills hycs_std and moves on to next simulation
  cmd <- paste0('cd ', rundir, ' && ./hycs_std >> stilt.log 2>&1')
  system(cmd, timeout = timeout)

  # Exit if running in HYSPLIT mode
  if (namelist[['ichem']] != 8) return()

  pf <- file.path(rundir, 'PARTICLE_STILT.DAT')
  if (!file.exists(pf)) {
    msg <- paste('Failed to output ', pf)
    warning(msg)
    cat(msg, '\n', file = file.path(rundir, 'stilt.log'), append = T)
    return()
  }
  
  n_lines <- count_lines(pf)
  if (n_lines < 2) {
    msg <- paste(pf, 'does not contain any trajectory data.')
    warning(msg)
    cat(msg, '\n', file = file.path(rundir, 'stilt.log'), append = T)
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
    xhgt_min <- min(output$receptor$zagl)
    xhgt_max <- max(output$receptor$zagl)
    xhgt_rng <- xhgt_max - xhgt_min
    xhgt_step <- xhgt_rng / numpar

    px <- data.frame(indx = 1:numpar)    
    px$xhgt <- (px$indx - 0.5) * xhgt_step + xhgt_min
    p <- merge(p, px, by = 'indx', sort = F)
  }
  
  # Calculate near-field dilution height based on gaussian plume width
  # approximation and recalculate footprint sensitivity for cases when the
  # plume height is less than the PBL height scaled by veght
  if (hnf_plume) 
    p <- calc_plume_dilution(p, numpar, output$receptor$zagl, namelist[['veght']])
  p
}
