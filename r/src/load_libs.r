#' load_libs multiple package loader/installer
#' @author Ben Fasoli
#'
#' Load and attach add-on packages. If the packages are not installed, load_libs
#' will attempt to install the package from CRAN.
#'
#' @param ... names of packages to load and attach
#' @param lib.loc a character vector describing the location of R library trees
#'   to search or install packages. The default value of NULL corresponds to
#'   libraries currently known to .libPaths()
#'
#' @examples
#' load_libs('dplyr', 'parallel')
#'
#' @import devtools
#' @export


load_libs <- function(..., lib.loc = NULL) {
  args <- list(...)

  repo <- 'https://cran.rstudio.com/'

  load_check <- function(pkg, lib.loc) {
    suppressWarnings(
      require(pkg, character.only = T, quiet = T,
              warn.conflicts = F, lib.loc = lib.loc)
    )
  }

  invisible({
    lapply(args, function(pkg) {
      # If package is not installed
      if (!load_check(pkg, lib.loc)) {
        # Try installing from CRAN
        try(
          suppressWarnings(
            install.packages(pkg, repo = repo, lib = lib.loc)
          )
        )
        # If the package is not found on CRAN
        if (!load_check(pkg, lib.loc)) {
          # Try Ben Fasoli's github
          if (!load_check('devtools', lib.loc))
            install.packages('devtools', repo = repo, lib = lib.loc)
          try(devtools::install_github(paste0('benfasoli/', pkg)))
          # Throw error if the package can't be found on CRAN or github
          if (!load_check(pkg, lib.loc)) {
            stop(paste('load_libs(): Failed to load', pkg))
          }
        }
      }
    })
    return(paste(args))
  })
}
