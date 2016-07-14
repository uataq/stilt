existsr <- function(xname, path="", checkit=FALSE) {
# see also assignr() and getr()
# similar to exists(), but
# checks for object on stored location (file path/.Rdataxname)
# checkit: flag if want function to attach and find r object
# 2/27/04 by CHG
#
#  $Id: existsr.r,v 1.7 2008-08-12 08:49:58 skoerner Exp $
#---------------------------------------------------------------------------------------------------

     zip <- TRUE
     xxname <- paste(path, ".RData", xname, sep="")
     if (!any(file.exists(xxname, paste(xxname, ".gz", sep="")))) return(FALSE) # not there
     if (file.exists(xxname)) zip <- FALSE # check if zipped
     if (!zip&checkit) {
       attach(xxname, pos=2)
       isthere <- exists(xname, where=2)
       detach(2)
       return(isthere)
     }
     if (zip) {print("compressed file; not tested"); return(TRUE)}
     if (!checkit) return(TRUE)
}
