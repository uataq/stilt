#' xtrct_grid extract spatial subdomain from ARL packed file.
#' @author Ben Fasoli
#'
#' @param input arl formatted meteorological file path
#' @param output subsetted arl formatted meteorological file path
#' @param workdir path to working directory containing xtrct_grid binary
#' @param xmn minimum x coordinate of subdomain
#' @param xmx maximum x coordinate of subdomain
#' @param ymn minimum y coordinate of subdomain
#' @param ymx maximum y coordinate of subdomain
#' @param levels number of vertical levels to include, defaults to NA which
#'   returns all vertical levels available in the input file
#'
#' @export

xtrct_grid <- function(input, output, workdir, xmn, xmx, ymn, ymx, levels = NA) {
  
  input_dirname <- file.path(dirname(input), '')  # ensures trailing slash
  input_basename <- basename(input)
  
  levels <- ifelse(is.na(levels), read_met_header(input)$nz, levels)
  
  stdin <- paste(c(
    input_dirname,
    input_basename,
    paste(ymn, xmn),
    paste(ymx, xmx),
    levels
  ), collapse = '\n')
  
  write(stdin, file.path(workdir, 'input'))
  cmd <- paste('cd', workdir, '&& ./xtrct_grid < input')
  proc <- system(cmd, input = stdin, ignore.stdout = T, ignore.stderr = T, intern = T)
  
  success <- is.null(attributes(proc)$status)
  if (!success) {
    return(F)
  }
  
  dir.create(dirname(output), showWarnings = F, recursive = T)
  file.rename(file.path(workdir, 'extract.bin'), output)
  
  T
}
