#' write_setup writes a SETUP.CFG namelist file to control the model behavior
#' @author Ben Fasoli
#'
#' Variable descriptions from "STILT How To v8", John C. Lin
#' @param numpar number of particles to be run; defaults to 100
#' @param delt integration timestep [min]; if set to 0.0, then timestep is
#'   dynamically determined
#' @param tratio maximum fraction of gridcell to be travelled by a particle in a
#'   single integration timestep. This determines the timestep if DELT is set to
#'   be dynamic
#' @param isot flag used to set the isotropic turbulence option; defaults to 0 to
#'   compute horizontal turbulence from wind field deformation. Setting to 1
#'   results in the horizontal turbulence to be the same in both the u and v
#'   directions
#' @param tlfrac the fraction of the lagrangian timescale TL to set as timestep
#'   in dispersion subroutine. The smaller this fraction is, the more finely the
#'   turbulence is resolved; defaults to 0.1
#' @param ndump flag to dump all particle/puff points at the end of a simulation
#'   to a file called PARDUMP. This can be read at the start of a new simulation
#'   to continue the previous calculation. Valid settings include 0 (no i/o),
#'   1 (read/write), 2 (read only), 3 (write only); defaults to 0
#' @param random flag that tells the random number generator whether to have a
#'   different random sequence for each model run (0 - false, 1 - true); defaults
#'   to 1
#' @param outdt interval [min] to output data to PARTICLE.DAT; defaults to 0.0,
#'   which outputs at every timestep
#' @param nturb no turbulence flag; defaults to 0, which includes turbulence
#'   rather than simulating mean trajectories
#' @param veght height below which a particle's time spent is tallied; defaults
#'   to 0.5, which specifies half of the PBL. Setting <=1.0 specifies a fraction
#'   of the boundary layer height, and setting >1.0 specifies a height above
#'   ground in meters
#' @param outfrac the fraction of the particles that are allowed to leave the
#'   model domain (given by met data); defaults to 0.9. If exceeded, the model
#'   stops
#' @param iconvect flag for convection. If set to 1, then runs excessive
#'   convection as described in Gerbig et al., 2003. For specialized RAMS output,
#'   the particles will be vertically redistributed according to the output
#'   convective mass fluxes; defaults to 0
#' @param winderrtf flag that specifies whether to have particle motions be
#'   affected by horizontal wind errors; defaults to 0. If set to 1, then STILT
#'   looks for a file called "WINDERR" that has four lines: 1. Standard deviation
#'   of errors [m/s], 2. Correlation timescale of errors [min], 3. Vertical
#'   correlation lengthscale of errors [m], 4. Horizontal correlation lengthscale
#'   of errors [km]. All statistical properties are applied equally in the u and
#'   v wind components
#' @param zicontroltf flag that specifies whether to scale the PBL heights in
#'   STILT uniformly in the entire model domain; defaults to 0. If set to 1, then
#'   STILT looks for a file called "ZICONTROL" that specifies the scaling for the
#'   PBL height. The first line indicates the number of hours that the PBL height
#'   will be changed, and each subsequent line indicates the scaling factor for
#'   that hour
#' @param varsiwant character vector of 4-letter hymodelc variables. Options:
#'   'time', 'indx', 'long', 'lati', 'zagl', 'sigw', 'tlgr', 'zsfc',
#'   'icdx', 'temp', 'samt', 'foot', 'shtf', 'tcld', 'dmas',
#'   'dens', 'rhfr', 'sphu', 'solw', 'lcld', 'zloc', 'dswf', 'wout',
#'   'mlht', 'rain', 'crai'
#'   defaults to 'time', 'indx', 'lati', 'long', 'zagl', 'foot'
#' @param file path and name for output file
#'
#' @export

write_setup <- function(numpar = 100, delt = 0, tratio = 0.9, isot = 0,
                        tlfrac = 0.1, ndump = 0, random = 1, outdt = 0,
                        nturb = 0, veght = 0.5, outfrac = 0.9, iconvect = 0,
                        winderrtf = 0, zicontroltf = 0, mgmin = 2000,
                        varsiwant = c('time', 'indx', 'lati', 'long', 'zagl'),
                        file = 'exe/SETUP.CFG') {

  if (basename(file) != 'SETUP.CFG')
    stop('write_setup(): file argument must end with SETUP.CFG')

  ivmax <- length(varsiwant)
  varsiwant <- paste0('\'', paste(varsiwant, collapse = '\', \''), '\'')

  eq <- function(lhs, rhs)
    paste0(lhs, '=', rhs, ',')

  txt <- c('$SETUP',
           eq('DELT', delt),
           'FRMR=0.0',
           eq('ICONVECT', iconvect),
           'INITD=0',
           eq('ISOT', isot),
           eq('IVMAX', ivmax),
           'KHMAX=9999',
           'KMIX0=250',
           'KMIXD=3',
           'KRND=6',
           eq('MGMIN', mgmin),
           eq('NDUMP', ndump),
           eq('NTURB', nturb),
           eq('NUMPAR', numpar ),
           eq('OUTDT', outdt),
           eq('OUTFRAC', outfrac),
           'QCYCLE=0',
           eq('RANDOM', random),
           eq('TLFRAC', tlfrac),
           eq('TRATIO', tratio),
           eq('VARSIWANT', varsiwant),
           eq('VEGHT', veght),
           eq('WINDERRTF', winderrtf),
           eq('ZICONTROLTF', zicontroltf),
           '$END')

  write(txt, file)
  return(file)
}
