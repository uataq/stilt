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
  
  # Calculate gaussian plume height using vertical velocity standard deviations
  p <- p %>%
    mutate(sigma = samt * sqrt(2) * sigw *
             sqrt(tlgr*abs(time*60) + tlgr^2 * exp(-abs(time*60)/tlgr)-1),
           pbl_mixing = veght * mlht)
  
  # Compare average gaussian plume height with dilution height calculated at
  # each time stepusing veght * mlht, or the fraction of the mixed layer height
  # (defaults to veght = 0.5) to apply the influence of a particle
  dilute <- p %>%
    group_by(time) %>%
    summarize_at(c('pbl_mixing', 'sigma'),
                 funs(mean(., na.rm = T))) %>%
    arrange(-time) %>%
    mutate(gauss_plume = r_zagl + cumsum(sigma),
           mask = gauss_plume < pbl_mixing)
  
  # Recalculate near-field foot values until the average gaussian plume height
  # grows to veght * mlht
  if (dilute$mask[1]) {
    dilute_recalc_end_time <- dilute$time[which(!dilute$mask)[1] - 1]
    mask <- p$time >= dilute_recalc_end_time
    p$foot[mask] <- (0.02884) / (p$gauss_plume[mask] * p$dens[mask]) * p$samt[mask]*60
  }
  
  return(p)
}
