#***************************************************************************************************
#  Reads initial CO2 concentration values from a TM3 file and attaches them to the "result" matrix.
#***************************************************************************************************

get.TM3CH4.bin <- function(yr4=NULL, mon=NULL, day=NULL, hr=NULL, ch4inifile=NULL,
                   result=NULL, result.sel=NULL) {
#
# --------------------------------------------------------------------------------------------------
#  Interface
#  =========
#
#  yr4          YYYY of the receptor time
#  mon          month "
#  day          day   "
#  co2inifile   File name (absolute path) of the TM3 binary data file, e.g.
#               "/Net/Groups/BSY/STILT/fluxes_input/TM3/mu1.0_070_mix_1995_fg.b". To auto-select
#               the year, specify a "YYYY" for it, it's getting substituted by the receptor year.
#  result       matrix with columns "btime" [hours backtime from yr4/mon/day], "lat", "lon" [degree
#               north, east] and "pres" [mb]
#               gets the existing (!) column "co2ini" overwritten on output
#  result.sel   selector of same length as the columns of "result": which rows should be considered
#
#
#  This is only a wrapper to the Fortran workhorse getTM3CH4.f90. See there for additional
#  documentation.
#
# --------------------------------------------------------------------------------------------------
#  $Id: get.TM3CH4.bin.r,v 1.3 2008-01-17 17:25:21 gerbig Exp $
# --------------------------------------------------------------------------------------------------

if (!is.loaded(tolower("getTM3CH4bin"))) dyn.load(paste(shlibpath,'getTM3CH4bin.so',sep=""))


#  length of Fortran string argument co2inifile255 must be 255
tmp <- "                                                                                          "
tmp <- paste(tmp,tmp,tmp)
ch4inifile255 <- substring(paste(ch4inifile,tmp),1,255)

btime <- result[result.sel,"btime"]
n     <- length(btime)
lat   <- result[result.sel,"lat"]
lon   <- result[result.sel,"lon"]
p     <- result[result.sel,"pres"]
xch4  <- p                                                  # dummy

x <- .Fortran("getTM3CH4bin", as.integer(yr4), as.integer(mon), as.integer(day), as.integer(hr),
         ch4inifile255, as.integer(n), as.single(as.vector(btime)), as.single(as.vector(lat)),
         as.single(as.vector(lon)), as.single(as.vector(p)), xch4=as.single(as.vector(xch4)) )$xch4

result[,"ch4ini"] <- rep(0,nrow(result))
x[x < -998.] <- NA
result[result.sel,"ch4ini"] <- x
return(result)

}
