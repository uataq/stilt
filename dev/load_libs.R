load_libs <- function(...) {
  args <- list(...)
  lapply(args, function(pkg) {
    if (!pkg %in% installed.packages()) {
      try(install.packages(pkg))
      if (!require(pkg, character.only = T, quiet = T)) {
        try(devtools::install_github(pkg, username = 'benfasoli'))
      }
    }
    if (!require(pkg, character.only = T, quiet = T)) {
      stop(paste('load_libs(): Failed to load', pkg))
    }
  })
  paste(args)
}

load_libs('dplyr', 'uataq')