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

calc_trajectory <- function(varsiwant, conage, cpack, delt, dxf, dyf, dzf,
                            emisshrs, frhmax, frhs, frme, frmr, frts, frvs,
                            hscale, ichem, iconvect, initd, isot, ivmax, kbls,
                            kblt, kdef, khmax, kmix0, kmixd, kmsl, kpuff, krnd,
                            kspl, kzmix, maxdim, maxpar, met_files, mgmin, ncycl,
                            ndump, ninit, numpar, nturb, n_hours, outdt, outfrac,
                            output, p10f, qcycle, random, splitf, tkerd, tkern,
                            rm_dat, timeout, tlfrac, tratio, tvmix, veght,
                            vscale, winderrtf, w_option, zicontroltf, ziscale, 
                            z_top, rundir) {

  require(uataq)

  # Enable manual rescaling of mixed layer height
  if (as.logical(zicontroltf)) {
    write_zicontrol(ziscale, file.path(rundir, 'ZICONTROL'))
  }

  # Write SETUP.CFG, CONTROL, and runhymodelc.sh files to control model
  write_setup(varsiwant, conage, cpack, delt, dxf, dyf, dzf, frhmax, frhs, frme,
              frmr, frts, frvs, hscale, ichem, iconvect, initd, isot, kbls, kblt, kdef,
              khmax, kmix0, kmixd, kmsl, kpuff, krnd, kspl, kzmix, maxdim,
              maxpar, mgmin, ncycl, ndump, ninit, numpar, nturb, outdt, outfrac,
              p10f, qcycle, random, splitf, tkerd, tkern, tlfrac, tratio, tvmix,
              veght, vscale, winderrtf, zicontroltf,
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
      msg <- paste('hymodelc timeout after', elapsed, ' seconds\n')
      warning(msg)
      cat(msg, '\n', file = file.path(rundir, 'ERROR'))
      return()
    }
    Sys.sleep(1)
  }

  # Error check hymodelc output
  pf <- file.path(rundir, 'PARTICLE.DAT')
  if (!file.exists(pf)) {
    msg <- paste('Failed to output PARTICLE.DAT in', basename(rundir),
                 'Check for errors in hymodelc.out')
    warning(msg)
    cat(msg, '\n', file = file.path(rundir, 'ERROR'))
    return()
  }

  n_lines <- uataq::count_lines(pf)
  if (n_lines < 2) {
    msg <- paste(pf, 'does not contain any trajectory data.',
                 'Check for errors in hymodelc.out')
    warning(msg)
    cat(msg, '\n', file = file.path(rundir, 'ERROR'))
    return()
  }

  # Read particle file, optionally remove PARTICLE.DAT in favor of compressed
  # .rds file, and return particle data frame
  p <- read_particle(file = pf, varsiwant = varsiwant)
  if (rm_dat) system(paste('rm', pf))
  
  # For column trajectories, preserve release height as xhgt
  if (length(output$receptor$zagl) > 1) {
    x_heights <- output$receptor$zagl
    px <- data.frame(indx = 1:numpar)
    px$xhgt <- rep(x_heights, each = length(px$indx) / length(x_heights))
    p <- merge(p, px, by = 'indx', sort = F)
  }
  p
}
