# Ben Fasoli

load_libs <- function(..., lib.loc = NULL) {
  args <- list(...)
  
  load_check <- function()
    require(pkg, character.only = T, quiet = T, lib.loc = lib.loc)
  
  lapply(args, load_check = load_check, function(pkg, load_check) {
    # If package is not installed
    if (!load_check()) {
      # Try installing from CRAN
      try( install.packages(pkg) )
      # If the package is not found on CRAN
      if (!load_check()) {
        # Try Ben Fasoli's github
        try( devtools::install_github(pkg, username = 'benfasoli') )
        # Stop if the package can't be found on CRAN or github
        if (!load_check()) {
          stop(paste('load_libs(): Failed to load', pkg))
        }
      }
    }
  })
  paste(args)
}
