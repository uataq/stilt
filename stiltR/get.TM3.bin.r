#***************************************************************************************************
#  Reads initial values from a TM3 file and attaches them to the "result" matrix.
#***************************************************************************************************

get.TM3.bin <- function(yr4=NULL, mon=NULL, day=NULL, hr=NULL, co2inifile=NULL,
                   result=NULL, result.sel=NULL) {

# --------------------------------------------------------------------------------------------------
#  Interface
#  =========
#
#  yr4          YYYY of the receptor time
#  mon          month "
#  day          day   "
#  co2inifile   file name (absolute path) of the TM3 binary data
#  result       matrix with columns "btime" [hours backtime from yr4/mon/day], "lat", "lon" [degree
#               north, east] and "pres" [mb]
#               gets the existing (!) column "co2ini" overwritten on output
#  result.sel   selector of same length as the columns of "result": which rows should be considered
#
#
# --------------------------------------------------------------------------------------------------
#  $Id: get.TM3.bin.r,v 1.2 2007-11-29 16:23:15 skoerner Exp $
# --------------------------------------------------------------------------------------------------

if (!is.loaded(tolower("getTM3bin"))) dyn.load(paste(shlibpath,'getTM3bin.so',sep=""))


#  length of Fortran string argument co2inifile255 must be 255
tmp <- "                                                                                          "
tmp <- paste(tmp,tmp,tmp)
co2inifile255 <- substring(paste(co2inifile,tmp),1,255)

btime <- result[result.sel,"btime"]
n     <- length(btime)
lat   <- result[result.sel,"lat"]
lon   <- result[result.sel,"lon"]
p     <- result[result.sel,"pres"]
xco2  <- p                                                  # dummy

x <- .Fortran("getTM3bin", as.integer(yr4), as.integer(mon), as.integer(day), as.integer(hr),
         co2inifile255, as.integer(n), as.single(as.vector(btime)), as.single(as.vector(lat)),
         as.single(as.vector(lon)), as.single(as.vector(p)), xco2=as.single(as.vector(xco2)) )$xco2

result[,"co2ini"] <- rep(0,nrow(result))
x[x < -998.] <- NA
result[result.sel,"co2ini"] <- x


return(result)

}
