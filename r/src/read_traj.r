#' read_traj reads in the output assosciated with each simulation step
#' @author James Mineau
#'
#' @param file filename argument. Must end with .rds (for serialized R data
#'   output, preferred), .h5 (for Hierarchical Data Format output),
#'   or .parquet (for Apache Parquet output).
#'   rds files do not require any additional libraries and have
#'   better compression. Alternative output formats require additional
#'   libraries.
#'
#' @import arrow
#' @import rhdf5
#' @export


read_traj <- function(file) {

  # .rds output (preferred)
  if (!is.null(file) && grepl('\\.rds$', file, ignore.case = T)) {
    return(readRDS(file))
  }

  output <- list()  # Define output list to rebuild

  # Alternative .h5 output
  if (!is.null(file) && grepl('\\.h5$', file, ignore.case = T)) {

    output$file <- file  # Assign current file path

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

  # Alternative .parquet output
  if (!is.null(file) && grepl('\\.parquet$', file, ignore.case = T)) {
    extract_meta <- function(meta, prefix) {
      prefix <- paste0('^', prefix, ':')
      keys <- grep(prefix, names(meta), value = T)
      data <- meta[keys]
      names(data) <- gsub(prefix, '', names(data))
      return(data)
    }

    convert_to_numeric <- function(x) {
      if (x == "NA") {
        return(NA)  # Handle the string "NA" case
      }
      num_val <- suppressWarnings(as.numeric(x))  # Attempt conversion
      if (is.na(num_val) && !is.na(x)) {
        return(x)  # Return original string if conversion fails
      }
      return(num_val)  # Otherwise, return the numeric value
    }

    output$file <- file  # Assign current file path

    # Read in parquet file
    particle <- arrow::read_parquet(file, as_data_frame = F)
    metadata <- particle$metadata

    # Get receptor metadata
    receptor <- extract_meta(metadata, 'receptor')
    receptor$run_time <- as.POSIXct(receptor$run_time, tz='UTC')
    receptor$lati <- as.numeric(receptor$lati)
    receptor$long <- as.numeric(receptor$long)
    receptor$zagl <- as.numeric(receptor$zagl)

    # Get param metadata
    params <- extract_meta(metadata, 'param')
    params <- lapply(params, convert_to_numeric)

    particle <- as.data.frame(particle)  # Convert to data frame

    # Handle particle error data
    err_cols <- grepl('_err$', names(particle))
    if (any(err_cols)) {
      # Split particle data into particle and particle_error
      particle <- particle[!err_cols]
      particle_error <- particle[err_cols]

      # Strip '_err' from particle_error names
      names(particle_error) <- gsub('_err$', '', names(particle_error))

      # Get particle error params
      particle_error_params <- extract_meta(metadata, 'error_param')
      particle_error_params <- lapply(particle_error_params, convert_to_numeric)

      output$particle_error <- particle_error
      output$particle_error_params <- particle_error_params
    }

    output$receptor <- receptor
    output$particle <- particle
    output$params <- params

    return(output)
  }
}
