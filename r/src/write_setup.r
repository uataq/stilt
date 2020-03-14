#' write_setup writes a SETUP.CFG namelist file to control the model behavior
#' @author Ben Fasoli
#'
#' SETUP.CFG contains namelist for controlling transport and dispersion behavior
#'
#' @param conage particle/puff conversions at conage (hours); defaults to 48
#' @param cpack binary concentration packing. 0 - all grid points written to
#'   file. 1 - only nonzero points written to file. 2 - special non-regular grid.
#'   Defaults to 1
#' @param delt integration timestep [min]; if set to 0.0, then timestep is
#'   dynamically determined
#' @param dxf horizontal x grid adjustment factor for ensemble; defaults to 1
#' @param dyf horizontal y grid adjustment factor for ensemble; defaults to 1
#' @param dzf vertical factor for ensemble; defaults to 0.01 ~ 250m
#' @param frhmax maximum value for horizontal rounding parameter; defaults to 3
#' @param frhs horizontal puff rounding fraction for merge; defaults to 1
#' @param frhs mass rounding fraction for enhanced merging; defaults to 0.1
#' @param frmr the fraction of the mass that is permitted to be removed at krnd
#'   intervals. For certain simulations, such as when a pollutant has a high
#'   ambient background relative, a small removal rate will significantly reduce
#'   the number of puffs on the grid at no loss in accuracy; defaults to 0
#' @param frts temporal puff rounding fraction for merge; defaults to 0.1
#' @param frvs vertical puff rounding fraction for merge; defaults to 0.1
#' @param hscale horizontal Lagrangian time scale (sec); defaults to 10800
#' @param ichem special chemistry or conversion modules. 1 - concentration grid
#'   treated as source-receptor matrix format. 2 - convert pollutant from species
#'   #1 to species #2. 3 - enable pm10 dust storm emission module. 4 - configure
#'   concentration grid similar to meteorology grid. 5 - treat 3D particle
#'   deposition using probability function. 7 - enable water surface transport
#'    of particle deposition. Defaults to 0
#' @param iconvect flag for convection. If set to 1, then runs excessive
#'   convection as described in Gerbig et al., 2003. For specialized RAMS output,
#'   the particles will be vertically redistributed according to the output
#'   convective mass fluxes; defaults to 0
#' @param initd determines model configuration as a particle or puff model. 0 -
#'   horizontal and vertical particle. 1 - horizontal gaussian puff, vertical top
#'   hat puff. 2 - horizontal and vertical top hat puff. 3 - horizontal gaussian
#'   puff, vertical particle. 4 - horizontal top-hat puff, vertical particle
#'   (HYSPLIT default). Defaults to 0 (horizontal and vertical particle)
#' @param isot flag used to set the isotropic turbulence option; defaults to 0 to
#'   compute horizontal turbulence from wind field deformation. Setting to 1
#'   results in the horizontal turbulence to be the same in both the u and v
#'   directions
#' @param kbls boundary layer stability derived from 1 - heat and momentum fluxes
#'   or 2 - wind and temperature profiles. Defaults to 1
#' @param kblt boundary layer turbulence parameterization. 1 - Beljaars/Holtslag
#'   and Betchov/Yaglom. 2 - Kanthar/Clayson. 3 - TKE field from input
#'   meteorology data file. 4 - velocity variances from input meteorology.
#'   Defaults to 1
#' @param kdef horizontal turbulence computation. 0 - in proportion to vertical
#'   turbulence. 1 - computed from velocity deformation. Defaults to 1
#' @param khmax maximum duration in hours of any particle/puff; defaults to 9999
#' @param kmix0 minimum mixing depth (abs(kmix0) used as minimum depth) in
#'   meters. Negative values force mixing heights coincident with model levels.
#'   Defaults to 250
#' @param kmixd mixing depth formulation. 0 - use input data MIXD if available,
#'   otherwise compute. 1 - compute from temperature profile (used for kmixd = 0
#'   or 2 if data is missing). 2 - compute from TKE profile. 3 - compute from
#'   bulk Ri profile. > = 10 use this value as a constant. Defaults to 3
#' @param kmsl starting heights default to AGL = 0 or MSL = 1; defaults to 0
#' @param kpuff horizontal puff dispersion. 0 - linear. 1 - square root. Defaults
#'   to 0
#' @param krnd enhanced merge interval (hours); defaults to 6
#' @param kspl standard splitting interval (hours); defaults to 1
#' @param kzmix vertical mixing adjustments. 0 - none, vertical diffusivity in
#'   PBL varies with height. 1 - vertical diffusivity in PBL single average
#'   value. 2 - scale PBL values by tvmix. 3 - scale free-troposphere values by
#'   tvmix. Defaults to 1
#' @param maxdim maximum number of pollutants to carry on one mass particle;
#'   defaults to 1
#' @param maxpar maximum number of particles carried in simulation; defaults to
#'   10,000
#' @param mgmin minimum meteorological subgrid size; defaults to 2000
#' @param ncycl pardump output cycle time; defaults to 0
#' @param ndump flag to dump all particle/puff points at the end of a simulation
#'   to a file called PARDUMP. This can be read at the start of a new simulation
#'   to continue the previous calculation. Valid settings include 0 (no i/o),
#'   1 (read/write), 2 (read only), 3 (write only); defaults to 0
#' @param ninit particle initialization. 0 - none. 1 - once. 2 - add. 3 -
#'   replace. Defaults to 1
#' @param nturb no turbulence flag; defaults to 0, which includes turbulence
#'   rather than simulating mean trajectories
#' @param numpar number of particles to be run; defaults to 100
#' @param outdt interval [min] to output data to PARTICLE.DAT; defaults to 0.0,
#'   which outputs at every timestep
#' @param outfrac the fraction of the particles that are allowed to leave the
#'   model domain (given by met data); defaults to 0.9. If exceeded, the model
#'   stops
#' @param p10f dust threshold velocity sensitivity factor; defaults to 1
#' @param qcycle optional cycling of emissions (hours); defaults to 0
#' @param random flag that tells the random number generator whether to have a
#'   different random sequence for each model run (0 - false, 1 - true); defaults
#'   to 1
#' @param splitf horizontal flow splitting. 0 - disable. 1 - automatic size
#'   adjustment. >1 - constant value. Defaults to 1
#' @param tkerd day (unstable) turbulent kinetic energy ratio; defaults to 0.18
#' @param tkern night (stable) turbulent kinetic energy ratio; defaults to 0.18
#' @param tlfrac the fraction of the lagrangian timescale TL to set as timestep
#'   in dispersion subroutine. The smaller this fraction is, the more finely the
#'   turbulence is resolved; defaults to 0.1
#' @param tratio maximum fraction of gridcell to be travelled by a particle in a
#'   single integration timestep. This determines the timestep if DELT is set to
#'   be dynamic
#' @param tvmix vertical mixing scale factor; defaults to 1
#' @param varsiwant character vector of 4-letter hymodelc variables. Options:
#'   'time', 'indx', 'long', 'lati', 'zagl', 'sigw', 'tlgr', 'zsfc',
#'   'icdx', 'temp', 'samt', 'foot', 'shtf', 'tcld', 'dmas',
#'   'dens', 'rhfr', 'sphu', 'solw', 'lcld', 'zloc', 'dswf', 'wout',
#'   'mlht', 'rain', 'crai'
#' @param veght height below which a particle's time spent is tallied; defaults
#'   to 0.5, which specifies half of the PBL. Setting <=1.0 specifies a fraction
#'   of the boundary layer height, and setting >1.0 specifies a height above
#'   ground in meters
#' @param vscale vertical Lagrangian time scale (sec); defaults to 200
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
#' @param file path and name for output file
#'
#' @export

write_setup <- function(varsiwant, conage = 48, cpack = 1, delt = 0, dxf = 1, 
                        dyf = 1, dzf = 0.01, frhmax = 3, frhs = 1, frme = 0.1,
                        frmr = 0, frts = 0.1, frvs = 0.1, hscale = 10800,
                        ichem = 0, iconvect = 0, initd = 0, isot = 0, kbls = 1,
                        kblt = 1, kdef = 1, khmax = 9999, kmix0 = 250, 
                        kmixd = 3, kmsl = 0, kpuff = 0, krnd = 6, kspl = 1, 
                        kzmix = 1, maxdim = 1, maxpar = 10000, mgmin = 2000, 
                        ncycl = 0, ndump = 0, ninit = 1, numpar = 100, 
                        nturb = 0, outdt = 0, outfrac = 0.9, p10f = 1, 
                        qcycle = 0, random = 1, splitf = 1, tkerd = 0.18, 
                        tkern = 0.18, tlfrac = 0.1, tratio = 0.9, tvmix = 1, 
                        veght = 0.5, vscale = 200, winderrtf = 0,
                        zicontroltf = 0, file = 'SETUP.CFG') {

  if (basename(file) != 'SETUP.CFG')
    stop('write_setup(): file argument must end with SETUP.CFG')

  ivmax <- length(varsiwant)
  varsiwant <- paste0('\'', paste(varsiwant, collapse = '\', \''), '\'')

  eq <- function(lhs, rhs) {
    if (is.logical(rhs))
      rhs <- as.numeric(rhs)
    paste0(lhs, '=', format(rhs, scientific = F), ',')
  }

  txt <- c('$SETUP',
           eq('CONAGE', conage),
           eq('CPACK', cpack),
           eq('DELT', delt),
           eq('DXF', dxf),
           eq('DYF', dyf),
           eq('DZF', dzf),
           eq('FRHMAX', frhmax),
           eq('FRHS', frhs),
           eq('FRME', frme),
           eq('FRMR', frmr),
           eq('FRTS', frts),
           eq('FRVS', frvs),
           eq('HSCALE', hscale),
           eq('ICHEM', ichem),
           eq('ICONVECT', iconvect),
           eq('INITD', initd),
           eq('ISOT', isot),
           eq('IVMAX', ivmax),
           eq('KBLS', kbls),
           eq('KBLT', kblt),
           eq('KDEF', kdef),
           eq('KHMAX', khmax),
           eq('KMIX0', kmix0),
           eq('KMIXD', kmixd),
           eq('KMSL', kmsl),
           eq('KPUFF', kpuff),
           eq('KRND', krnd),
           eq('KSPL', kspl),
           eq('KZMIX', kzmix),
           eq('MAXDIM', maxdim),
           eq('MAXPAR', maxpar),
           eq('MGMIN', mgmin),
           eq('NCYCL', ncycl),
           eq('NDUMP', ndump),
           eq('NINIT', ninit),
           eq('NTURB', nturb),
           eq('NUMPAR', numpar ),
           eq('OUTDT', outdt),
           eq('OUTFRAC', outfrac),
           eq('P10F', p10f),
           eq('QCYCLE', qcycle),
           eq('RANDOM', random),
           eq('SPLITF', splitf),
           eq('TKERD', tkerd),
           eq('TKERN', tkern),
           eq('TLFRAC', tlfrac),
           eq('TRATIO', tratio),
           eq('TVMIX', tvmix),
           eq('VARSIWANT', varsiwant),
           eq('VEGHT', veght),
           eq('VSCALE', vscale),
           eq('WINDERRTF', winderrtf),
           eq('ZICONTROLTF', zicontroltf),
           '$END')

  write(txt, file)
  file
}
