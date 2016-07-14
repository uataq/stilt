#
#  $Id: my.file.exists.r,v 1.2 2007-06-27 11:54:03 skoerner Exp $
#---------------------------------------------------------------------------------------------------

my.file.exists <- function(fname) {
  out <- file.exists(fname)
  if (!out) {
    # check again using ls command
    test.out <- system(paste("/bin/ls",fname,sep=" "),ignore.stderr=T)
    out <- test.out == 0
  }
  out
}
