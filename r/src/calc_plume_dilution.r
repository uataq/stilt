#' Recalculate influence in HNF by modifying dilution depth
#' @author Ben Fasoli
#'
#' Modifies foot influence by recalculating dilution depth using estimated
#' mixing depth based on vertical velocity standard deviations (vertical
#' turbulence) and the Lagrangian decorrelation timescale.
#'
#' @param p particle data frame including named dens, samt, sigw, tlgr,
#'   and foot columns
#' @param numpar number of particles in simulation
#' @param r_zagl receptor height, in meters above ground, or above sea level
#'   when kmsl is 1
#' @param veght height below which a particle is counted, as HYSPLIT reads
#'   VEGHT: a fraction of the PBL height when 1 or less, meters above ground
#'   otherwise
#' @param kmsl 1 if release heights are above sea level, as HYSPLIT reads KMSL,
#'   which also requires zsfc in p; defaults to 0, above ground
#'
#' @import dplyr
#' @export

calc_plume_dilution <- function(p, numpar, r_zagl, veght, kmsl = 0) {

  require(dplyr)

  varsineed <- c('dens', 'samt', 'sigw', 'tlgr', 'foot', if (kmsl == 1) 'zsfc')
  if (!all(varsineed %in% names(p)))
    stop('calc_plume_dilution(): varsiwant must include: ',
         paste(collapse = ', ', varsineed))

  # HYSPLIT reads VEGHT as a fraction of the mixed-layer height when it is 1 or
  # less and as meters above ground otherwise. Compare the plume with the same
  # height HYSPLIT counted samt below
  p %>%
    mutate(foot_no_hnf_dilution = foot,
           sigma = samt * sqrt(2) * sigw *
             sqrt(tlgr*abs(time*60) + tlgr^2 * exp(-abs(time*60)/tlgr)-1),
           pbl_mixing = if (veght <= 1) veght * mlht else veght) %>%
    # Accumulate the spread in order of particle age. time is negative on a
    # backward run and positive on a forward run, so sort on its magnitude
    arrange(abs(time)) %>%
    group_by(indx) %>%
    # pbl_mixing is above ground, so the plume is too. Above sea level, it starts
    # at the release height minus the terrain at each particle's first output
    # row, or at the ground if the release height is below that terrain
    mutate(release = ifelse('xhgt' %in% names(p), xhgt, r_zagl),
           release = if (kmsl == 1) pmax(release - zsfc[1], 0) else release,
           plume = release + cumsum(sigma),
           # A plume starting at the ground has no depth until sigma is above 0
           foot = ifelse(plume > 0 & plume < pbl_mixing,
                         0.02897 / (plume * dens) * samt*60,
                         foot)) %>%
    dplyr::select(-sigma, -pbl_mixing, -release, -plume) %>%
    ungroup()
}
