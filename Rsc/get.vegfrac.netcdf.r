get.vegfrac.netcdf <- function(date=NULL, i=NULL, j=NULL, ires=NULL, jres=NULL, dpath=NULL) {
#---------------------------------------------------------------------------------------------------
#
# Function to extract vegetation fraction with spatial aggregation from a map stored as netCDF file.
# This is only a wrapper for the Fortran subroutine vegfracread in "vegfrac_read.f90". See there for
# additional documentation.
#
# Interface
# ---------
# date      : a vector of reference dates, each date of format c(YYYY, MM, DD, HH, mm),
#             the total vector concatenates the individual dates.
# i, j      : particle location indices, referencing to the respective grids at resolution ires, jres
# ires, jres: aggregation parameters for the surface indices
# dpath     : path to the netCDF files
#
#
# On return, VEG_FRAC[n,1:8], e.g., contains the vegetation fractions of the 8 classes corresponding
# to date n and location (i(n),j(n)). Always TRUE is sum(VEG_FRA[2,]) == 1.
#
#---------------------------------------------------------------------------------------------------
# $Id: get.vegfrac.netcdf.r,v 1.1 2007-11-23 13:57:31 gerbig Exp $
#---------------------------------------------------------------------------------------------------

###################
### for testing
#shlibpath <- "/Net/Groups/BSY/tpobw/STILT/VPRM/"
#n     <- 2
#date  <- as.vector(c( c(2005,5,17,15,9), c(2005,4,30,23,59) ))
#i     <- as.vector(c(1,18))
#j     <- as.vector(c(39,145))
#ires  <- as.vector(c(1,4))
#jres  <- as.vector(c(1,2))
#dpath <- "/Net/Groups/BSY/STILT/fluxes_input/veg-related/"
### end for testing
###################

if (!any(names(getLoadedDLLs())=="getvegfracNetCDF")) dyn.load(paste(shlibpath,"getvegfracNetCDF.so",sep=""))

n<-length(i)

#  length of dpath MUST be 255 when calling a Fortran subroutine
tmp <- "                                                                                          "
dpath <- substring(paste(dpath,tmp,tmp,tmp),1,255)

res <- .Fortran('getvegfracNetCDF',
                as.character(dpath), as.integer(n), as.integer(date),
                as.integer(i), as.integer(j), as.integer(ires), as.integer(jres),
		VEG_FRA=single(length=8*n))

VEG_FRA <- matrix(res$VEG_FRA,ncol=8,nrow=n, byrow=T)
return(VEG_FRA)

}
