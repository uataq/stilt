#' stilt_apply parallel apply function selection
#' @author Ben Fasoli
#'
#' Chooses apply function based on user parameters and available parallelization
#' options. Uses lapply if n_cores is 1, slurm_apply if SLURM job management is
#' available, and mclapply if SLURM is not available.
#'
#' @param X a vector to apply function FUN over
#' @param FUN the function to be applied to each element of X
#' @param slurm logical that forces job submission via SLURM
#' @param slurm_options a named list of options recognized by \code{sbatch};
#'   passed to rslurm::slurm_apply()
#' @param n_nodes number of nodes to submit SLURM jobs to using \code{sbatch}
#' @param n_cores number of CPUs to utilize per node
#' @param ... arguments to FUN
#'
#' @return if using slurm, returns sjob information. Otherwise, will return a
#'   TRUE for every successful model completion
#'
#' @export

stilt_apply <- function(X, FUN, slurm = F, slurm_options = list(),
                        n_nodes = 1, n_cores = 1, ...) {

  if (slurm) {
    sbatch_avail <- system('which sbatch', intern = T)
    if (length(sbatch_avail) == 0 || nchar(sbatch_avail[1]) == 0)
      stop('Problem identifying sbatch executable for slurm...')

    print('Parallelization using slurm. Dispatching jobs...')
    load_libs('rslurm')
    Y <- data.frame(X = X, slurm = T, ..., stringsAsFactors = F)
    sjob <- rslurm::slurm_apply(FUN, Y,
                                jobname = basename(getwd()), pkgs = 'base',
                                nodes = n_nodes, cpus_per_node = n_cores,
                                slurm_options = slurm_options)
    uataq::br(2)
    return(sjob)
  } else if (n_cores > 1) {
    print('Parallelization using multiple R jobs. Dispatching processes...')
    load_libs('parallel')
    cl <- makeForkCluster(n_cores, outfile = '')
    out <- parallel::parLapply(cl, X, FUN, ...)
    stopCluster(cl)
    return(out)
  }

  print('Parallelization disabled...')
  return(lapply(X, FUN, ...))
}
