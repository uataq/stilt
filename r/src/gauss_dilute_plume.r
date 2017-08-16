#' gauss_dilute_plume modifies foot influence by recalculating dilution depth
#'   using estimated gaussian plume height based on average vertical velocity
#'   standard deviations at each time step
#' @author Ben Fasoli
#' 
#' @param p particle data frame including named dens, tlgr, sigw, and foot
#'   columns
#' @param numpar number of particles in simulation
#' @param r_zagl receptor height above ground, in meters
#' @param veght fraction of pbl height below which particle is counted
#'
#' @import dplyr
#' @export

gauss_dilute_plume <- function(p, numpar, r_zagl, veght) {
  
  require(dplyr)
  
  p %>%
    mutate(sigma = samt * sqrt(2) * sigw *
             sqrt(tlgr*abs(time*60) + tlgr^2 * exp(-abs(time*60)/tlgr)-1),
           pbl_mixing = veght * mlht) %>%
    arrange(-time) %>%
    group_by(indx) %>%
    mutate(gauss_plume = r_zagl + cumsum(sigma),
           foot = ifelse(gauss_plume < pbl_mixing,
                         0.02884 / (gauss_plume * dens) * samt*60,
                         foot)) %>%
    ungroup() %>%
    return()
}
