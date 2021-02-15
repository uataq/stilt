#' stilt_apply parallel apply function selection
#' @author Ben Fasoli
#'
#' Chooses apply function based on user parameters and available parallelization
#' options. Uses lapply if n_cores is 1, slurm_apply if SLURM job management is
#' available, and mclapply if SLURM is not available.
#'
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

stilt_apply <- function(FUN, slurm = F, slurm_options = list(),
                        n_nodes = 1, n_cores = 1, ...) {
  
  if (!slurm && n_nodes > 1) {
    stop('n_nodes > 1 but but slurm is disabled. ',
         'Did you mean to set slurm = T in run_stilt.r?')
  }
  
  # Expand arguments to form a data frame where rows serve as iterations of FUN
  # using named columns as arguments to FUN. Base R's data.frame class expands
  # list arguments to multiple rows when constructing a new dataframe (not when
  # adding a new column) but tibble retains list classed objects
  Y <- tibble(...)

  message('Initializing STILT')
  message('Commit ID: ', find_git_commit_id())
  message('Number of receptors: ', nrow(Y))
  message('Number of parallel workers: ', n_nodes * n_cores)

  if (slurm) {
    # Confirm availability of sbatch executable and dispatch simulation
    # configurations to SLURM
    sbatch_avail <- system('which sbatch', intern = T)
    if (length(sbatch_avail) == 0 || nchar(sbatch_avail[1]) == 0)
      stop('Problem identifying sbatch executable for slurm...')
    
    # Shuffle receptor order for quasi-load balancing
    if (n_nodes > 1 || n_cores > 1) Y <- Y[sample.int(nrow(Y), nrow(Y)), ]
    
    message('Multi node parallelization using slurm. Dispatching jobs...')
    load_libs('rslurm')
    sjob <- rslurm::slurm_apply(FUN, Y,
                                jobname = basename(getwd()), pkgs = 'base',
                                nodes = n_nodes,
                                cpus_per_node = n_cores,
                                preschedule_cores = F,
                                slurm_options = slurm_options)
    return(invisible(sjob))
  }
  
  if (n_cores > 1) {
    # Load parallel backend and dispatch simulations to worker processes using
    # dynamic load balancing
    message('Single node parallelization. Dispatching worker processes...')
    load_libs('parallel')
    out <- do.call(parallel::mcmapply, c(
      FUN = FUN, 
      Y,
      mc.cores = n_cores,
      mc.preschedule = F,
      SIMPLIFY = F))
    return(invisible(out))
  }
  
  # Call FUN for each row of Y
  message('Parallelization disabled. Executing simulations sequentially...')
  for (i in 1:nrow(Y)) out <- do.call(FUN, Y[i, ])
  return(invisible(out))
}
