#' write_runhymodelc writes a runhymodelc.sh batch file to execute model
#' @author Ben Fasoli
#'
#' @param file full path to desired runhymodelc.sh file that contains a
#'   hymodelc executable
#'
#' @export

write_runhymodelc <- function(file) {

  if (basename(file) != 'runhymodelc.sh')
    stop('write_runhymodelc(): file argument must end with runhymodelc.sh')

  txt  <- c(paste('cd', dirname(file)),
            './hymodelc &> ./hymodelc.out 2>&1 & echo $!')

  write(txt, file)
  return(file)
}
