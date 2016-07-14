get.modis.netcdf <- function(date=NULL, i=NULL, j=NULL, ires=NULL, jres=NULL, dpath=NULL, 
                             subill.loc=subill, subjll.loc=subjll, subiur.loc=subiur, subjur.loc=subjur,nv=n.vegclass) {
#---------------------------------------------------------------------------------------------------
#
# Function to extract surface indices with spatial aggregation from maps stored in netCDF files.
# This is only a wrapper for the Fortran subroutine VPRMread in "VPRM_read.f90". See there for
# additional documentation.
#
# Interface
# ---------
# date      : a vector of reference dates, each date of format c(YYYY, MM, DD, HH, mm),
#             the total vector concatenates the individual dates.
# i, j      : particle location indices, referencing to the respective grids at resolution ires, jres
# ires, jres: aggregation parameters for the surface indices
# dpath     : path to the netCDF files
# subill.loc etc.: nested area corner points (grid indices, not lat-lon)
#
# On return, the list component $EVI[n,1:8], e.g., contains the EVI index of the 8 vegetation classes
# corresponding to date n. Always TRUE is sum(SfIndices$VEG_FRA[2,]) == 1.
#
#---------------------------------------------------------------------------------------------------
# $Id: get.modis.netcdf.r,v 1.4 2009-07-03 14:53:43 gerbig Exp $
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

if (!any(names(getLoadedDLLs())=="getMODISnetCDF")) dyn.load(paste(shlibpath,"getMODISnetCDF.so",sep=""))

n<-length(i)

#  length of dpath MUST be 255 when calling a Fortran subroutine
tmp <- "                                                                                          "
dpath <- substring(paste(dpath,tmp,tmp,tmp),1,255)

res <- .Fortran('getMODISnetCDF',
                as.character(dpath), as.integer(n), as.integer(date),
                as.integer(i), as.integer(j), as.integer(ires), as.integer(jres),
                as.integer(subill.loc), as.integer(subjll.loc), as.integer(subiur.loc), as.integer(subjur.loc),
                as.integer(nv), 
                EVI=single(length=nv*n), EVI_amax=single(length=nv*n), EVI_amin=single(length=nv*n),
                LSWI=single(length=nv*n), LSWI_amax=single(length=nv*n), LSWI_amin=single(length=nv*n),
		VEG_FRA=single(length=nv*n))


SfIndices <- list(matrix(res$EVI,ncol=nv,nrow=n, byrow=T), matrix(res$EVI_amax,ncol=nv,nrow=n, byrow=T), matrix(res$EVI_amin,ncol=nv,nrow=n, byrow=T),
                  matrix(res$LSWI,ncol=nv,nrow=n, byrow=T), matrix(res$LSWI_amax,ncol=nv,nrow=n, byrow=T), matrix(res$LSWI_amin,ncol=nv,nrow=n, byrow=T),
		  matrix(res$VEG_FRA,ncol=nv,nrow=n, byrow=T))
names(SfIndices) <- c("EVI", "EVI_amax", "EVI_amin", "LSWI", "LSWI_amax", "LSWI_amin", "VEG_FRA")
return(SfIndices)

}
