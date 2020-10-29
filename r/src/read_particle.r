#' read_particle parses the hysplit PARTICLE.DAT output file
#' @author Ben Fasoli
#'
#' Reads PARTICLE.DAT hysplit output file into data frame
#'
#' @param file location of PARTICLE.DAT file
#' @param varsiwant character vector of 4-letter hysplit variables
#'
#' @import dplyr
#' @export

read_particle <- function(file, varsiwant) {

  require(dplyr)

  n_lines <- count_lines(file)

  if (n_lines < 2) {
    warning(paste('read_particle(): only 1 line found in', file))
    return(NULL)
  }

  read.table(file, header = F, skip = 1,
             colClasses = 'numeric', col.names = varsiwant,
             stringsAsFactors = F)
}
