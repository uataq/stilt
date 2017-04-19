#' write_runhymodelc writes a runhymodelc.bat batch file to execute model
#' @author Ben Fasoli
#'
#' @param file full path to desired runhymodelc.bat file that contains a
#'   hymodelc executable
#'
#' @export

write_runhymodelc <- function(file) {

  if (basename(file) != 'runhymodelc.bat')
    stop('write_runhymodelc(): file argument must end with runhymodelc.bat')

  txt  <- c(paste('cd', dirname(file)),
            './hymodelc &> hymodelc.out')

  write(txt, file)
  return(file)
}
