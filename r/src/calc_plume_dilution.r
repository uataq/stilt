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
#' @param r_zagl receptor height above ground, in meters
#' @param veght fraction of pbl height below which particle is counted
#'
#' @import dplyr
#' @export

calc_plume_dilution <- function(p, numpar, r_zagl, veght) {

  require(dplyr)

  varsineed <- c('dens', 'samt', 'sigw', 'tlgr', 'foot')
  if (!all(varsineed %in% names(p)))
    stop('calc_plume_dilution(): varsiwant must include: ',
         paste(collapse = ', ', varsineed))

  p %>%
    mutate(foot_no_hnf_dilution = foot,
           sigma = samt * sqrt(2) * sigw *
             sqrt(tlgr*abs(time*60) + tlgr^2 * exp(-abs(time*60)/tlgr)-1),
           pbl_mixing = veght * mlht) %>%
    arrange(-time) %>%
    group_by(indx) %>%
    mutate(plume = ifelse('xhgt' %in% names(p), xhgt, r_zagl) + cumsum(sigma),
           foot = ifelse(plume < pbl_mixing,
                         0.02897 / (plume * dens) * samt*60,
                         foot)) %>%
    dplyr::select(-sigma, -pbl_mixing, -plume) %>%
    ungroup()
}
