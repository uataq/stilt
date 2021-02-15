#' Find git commit id
#' @author Ben Fasoli
#'
#' @export

find_git_commit_id <- function() {
  system('git rev-parse HEAD', intern = T)
}
