#' read_traj reads in the output assosciated with each simulation step
#' @author James Mineau
#'
#' @param file filename argument. Must end with .rds (for serialized R data
#'   output, preferred), .h5 (for Hierarchical Data Format output),
#'   or NULL to return the footprint object for programatic use.
#'   rds files do not require any additional libraries and have
#'   better compression. The \code{rhdf5} package is required for .h5 output.
#'   HDF5 output is slightly less efficient than .rds output, but is more
#'   portable and can be read by other languages.
#'
#' @import rhdf5
#' @export


read_traj <- function(file) {

  # .rds output (preferred)
  if (!is.null(file) && grepl('\\.rds$', file, ignore.case = T)) {
    return(readRDS(file))
  }

  output <- list()  # Define output list to rebuild

  # .h5 output
  if (!is.null(file) && grepl('\\.h5$', file, ignore.case = T)) {

    # Assign current file path
    output$file <- file

    # Get receptor attributes
    output$receptor <- rhdf5::h5readAttributes(file, 'receptor')

    # Get particle data
    output$particle <- rhdf5::h5read(file, 'particle')

    # Get params attributes
    output$params <- rhdf5::h5readAttributes(file, 'params')
    output$params[output$params == "NA"] <- NA

    # Get particle error info if it exists
    if ('particle_error' %in% rhdf5::h5ls(file)$name) {
      output$particle_error <- rhdf5::h5read(file, 'particle_error')
      output$particle_error_params <- rhdf5::h5readAttributes(file, 'particle_error_params')
      output$particle_error_params[output$particle_error_params == "NA"] <- NA
    }

    return(output)
  }
}
