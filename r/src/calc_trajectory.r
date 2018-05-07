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
#' @import uataq
#' @export

calc_trajectory <- function(varsiwant, delt, emisshrs, frmr, iconvect, initd,
                            isot, khmax, kmix0, kmixd, krnd, met_files, mgmin,
                            ndump, numpar, nturb, n_hours, outdt, outfrac,
                            output, random, rm_dat, timeout, tlfrac, tratio,
                            veght, winderrtf, w_option, zicontroltf, z_top,
                            rundir) {

  require(uataq)

  # Write SETUP.CFG, CONTROL, and runhymodelc.sh files to control model
  write_setup(varsiwant, delt, frmr, iconvect, initd, isot, khmax, kmix0, kmixd,
              krnd, mgmin, ndump, numpar, nturb, outdt, outfrac, random, tlfrac,
              tratio, veght, winderrtf, zicontroltf,
              file.path(rundir, 'SETUP.CFG'))
  write_control(output$receptor, emisshrs, n_hours, w_option, z_top, met_files,
                file.path(rundir, 'CONTROL'))
  sh <- write_runhymodelc(file.path(rundir, 'runhymodelc.sh'))

  # Simulation timeout ---------------------------------------------------------
  # Monitors time elapsed running hymodelc. If elapsed time exceeds timeout
  # specified in run_stilt.r, kills hymodelc and moves on to next simulation
  of <- file.path(rundir, 'hymodelc.out')
  eval_start <- Sys.time()
  pid <- system(paste('bash', sh), intern = T)
  on.exit(tools::pskill(pid))
  repeat {
    elapsed <- as.double.difftime(Sys.time() - eval_start, units = 'secs')
    if (!pid_is_active(pid)) {
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
    cat('PARTICLE.DAT does not contain any trajectory data. Check for errors ',
        'in hymodelc.out\n', file = file.path(rundir, 'ERROR'))
    return()
  }

  # Read particle file, optionally remove PARTICLE.DAT in favor of compressed
  # .rds file, and return particle data frame
  p <- read_particle(file = pf, varsiwant = varsiwant)
  if (rm_dat) system(paste('rm', pf))
  p
}
