#***************************************************************************************************
#  read dimension of VPRM fields from netcdf files
#***************************************************************************************************

get.vprm.dim <- function(evilswipath) {
#evilswipath: input path for VPRM EVI and LSWI files

      if (!is.element("package:ncdf", search())) library("ncdf") # Load ncdf library

      #check directory, get EVI_MAX file
      fnames<-dir(evilswipath)
      fname<-fnames[substring(fnames,1,nchar("VPRM_input_VEG_FRA_"))=="VPRM_input_VEG_FRA_"][1]
      evifile <- open.ncdf(paste(evilswipath,fname,sep=""), write=F)

      # get dimensions
      if(regexpr("d01",fname)>0){ #old format (before June 2009)
        nx<-length(get.var.ncdf(evifile,"west_east"))
        ny<-length(get.var.ncdf(evifile,"south_north"))
        nv<-length(get.var.ncdf(evifile,"vprm_classes"))
      } else {#new format
        nx<-length(get.var.ncdf(evifile,"lon"))
        ny<-length(get.var.ncdf(evifile,"lat"))
        nv<-length(get.var.ncdf(evifile,"vprm_classes"))
      }
      close.ncdf(evifile)
      result<-c(nx,ny,nv)
      names(result)<-c("nx","ny","nv")
      return(result)
}
