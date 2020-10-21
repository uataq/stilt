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

calc_trajectory <- function(varsiwant,
                            conage,
                            cpack,
                            delt,
                            dxf,
                            dyf,
                            dzf,
                            emisshrs,
                            frhmax,
                            frhs,
                            frme,
                            frmr,
                            frts,
                            frvs,
                            hnf_plume,
                            hscale,
                            ichem,
                            iconvect,
                            initd,
                            isot,
                            ivmax,
                            kbls,
                            kblt,
                            kdef,
                            khmax,
                            kmix0,
                            kmixd,
                            kmsl,
                            kpuff,
                            krnd,
                            kspl,
                            kzmix,
                            maxdim,
                            maxpar,
                            met_files,
                            mgmin,
                            ncycl,
                            ndump,
                            ninit,
                            numpar,
                            nturb,
                            n_hours,
                            outdt,
                            outfrac,
                            output,
                            p10f,
                            qcycle,
                            random,
                            splitf,
                            tkerd,
                            tkern,
                            rm_dat,
                            timeout,
                            tlfrac,
                            tratio,
                            tvmix,
                            veght,
                            vscale,
                            winderrtf,
                            w_option,
                            zicontroltf,
                            ziscale, 
                            z_top,
                            rundir) {
  
  # Enable manual rescaling of mixed layer height
  if (as.logical(zicontroltf)) {
    write_zicontrol(ziscale, file.path(rundir, 'ZICONTROL'))
  }

  # Write SETUP.CFG and CONTROL files to control model
  write_setup(varsiwant, conage, cpack, delt, dxf, dyf, dzf, frhmax, frhs, frme,
              frmr, frts, frvs, hscale, ichem, iconvect, initd, isot, kbls, kblt, 
              kdef, khmax, kmix0, kmixd, kmsl, kpuff, krnd, kspl, kzmix, maxdim,
              maxpar, mgmin, ncycl, ndump, ninit, numpar, nturb, outdt, outfrac,
              p10f, qcycle, random, splitf, tkerd, tkern, tlfrac, tratio, tvmix,
              veght, vscale, winderrtf, zicontroltf,
              file.path(rundir, 'SETUP.CFG'))
  write_control(output$receptor, emisshrs, n_hours, w_option, z_top, met_files,
                file.path(rundir, 'CONTROL'))
  
  # Simulation timeout ---------------------------------------------------------
  # Monitors time elapsed running hymodelc. If elapsed time exceeds timeout
  # specified in run_stilt.r, kills hymodelc and moves on to next simulation
  cmd <- paste0('cd ', rundir, ' && ./hymodelc > hymodelc.out')
  system(cmd, timeout=timeout)

  # Error check hymodelc output
  pf <- file.path(rundir, 'PARTICLE.DAT')
  if (!file.exists(pf)) {
    msg <- paste('Failed to output PARTICLE.DAT in', basename(rundir),
                 'Check for errors in hymodelc.out')
    warning(msg)
    cat(msg, '\n', file = file.path(rundir, 'ERROR'))
    return()
  }

  n_lines <- count_lines(pf)
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

  # Calculate near-field dilution height based on gaussian plume width
  # approximation and recalculate footprint sensitivity for cases when the
  # plume height is less than the PBL height scaled by veght
  if (hnf_plume) 
    p <- calc_plume_dilution(p, numpar, output$receptor$zagl, veght)
  p
}
