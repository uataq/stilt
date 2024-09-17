#' write_traj writes out the output assosciated with each simulation step
#' @author James Mineau
#'
#' @param traj output created during simulation_step
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

write_traj <- function(traj, file) {

  # Remove existing trajectory file
  if (!is.null(file) && file.exists(file))
    system(paste('rm', file))

  # .rds output (preferred)
  if (!is.null(file) && grepl('\\.rds$', file, ignore.case = T)) {
    saveRDS(traj, file)
    return(file)
  }

  # .h5 output
  if (!is.null(file) && grepl('\\.h5$', file, ignore.case = T)) {
    # Check if rhdf5 package is installed
    if (!requireNamespace("rhdf5", quietly = TRUE)) {
      stop("The 'rhdf5' package is required for hdf5 output. Please install it.")
    }

    # Create HDF5 file
    fid <- rhdf5::H5Fcreate(file)

    # Write file info
    rhdf5::h5writeAttribute(file, fid, 'file')

    # Create groups
    rhdf5::h5createGroup(fid, 'receptor')
    rhdf5::h5createGroup(fid, 'params')

    # Write receptor data
    receptor <- rhdf5::H5Gopen(fid, 'receptor')
    rhdf5::h5writeAttribute(as.numeric(traj$receptor$run_time),
                            receptor, 'run_time')
    rhdf5::h5writeAttribute(traj$receptor$lati, receptor, 'lati')
    rhdf5::h5writeAttribute(traj$receptor$long, receptor, 'long')
    rhdf5::h5writeAttribute(traj$receptor$zagl, receptor, 'zagl')
    rhdf5::H5Gclose(receptor)

    # Write particle dataset
    rhdf5::h5write(traj$particle, fid, 'particle')

    # Create params group
    params <- rhdf5::H5Gopen(fid, 'params')
    for (p in names(traj$params)) {
      val <- traj$params[[p]]
      if (is.na(val)) val <- 'NA'
      rhdf5::h5writeAttribute(val, params, p)
    }
    rhdf5::H5Gclose(params)

    # Create particle error group if it exists
    if ('particle_error' %in% names(traj)) {

      # Write particle_error dataset
      rhdf5::h5write(traj$particle_error, fid, 'particle_error')

      # Create particle_error_param group
      rhdf5::h5createGroup(fid, 'particle_error_params')
      particle_error_params <- rhdf5::H5Gopen(fid, 'particle_error_params')
      for (p in names(traj$particle_error_params)) {
        val <- traj$particle_error_params[[p]]
        if (is.na(val)) val <- 'NA'
        rhdf5::h5writeAttribute(val, particle_error_params, p)
      }
      rhdf5::H5Gclose(particle_error_params)
    }

    # Close file
    rhdf5::H5Fclose(fid)

    return(file)
  }

  return(traj)
}
