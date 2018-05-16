#' read_config parses the CONC.CFG output file
#' @author Ben Fasoli
#'
#' Reads CONC.CFG hymodelc output file into named list
#'
#' @param file location of CONC.CFG file
#'
#' @export

read_config <- function(file) {

  n_lines <- uataq::count_lines(file)

  if (n_lines < 2) {
    warning(paste('read_config(): only 1 line found in', file))
    return(NULL)
  }

  config <- grep('=', readLines(file), fixed = T, value = T)
  config_min <- strsplit(gsub('\\s+|\'|,', '', config), '=')

  config_list <- list()
  for (i in 1:length(config_min)) {
    key <- config_min[[i]][1]
    val <- config_min[[i]][2]
    config_list[[key]] <- type.convert(val, as.is = T)
  }

  config_list
}
