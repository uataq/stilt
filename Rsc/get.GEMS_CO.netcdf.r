#***************************************************************************************************
#  read initial values from CarbonTracker file and attach them to the "result" matrix
#***************************************************************************************************

get.GEMS_CO.netcdf <- function(yr4=NULL, mon=NULL, day=NULL, hr=NULL, co2inifile=NULL,
                             result=NULL, result.sel=NULL, spec="coini") {
#ls -hl /Net/Groups/BSY/people/vbeck/WRF_runs/BC_WRF/GEMS_reanalysis
#co2inifile="/Net/Groups/BSY/people/vbeck/WRF_runs/BC_WRF/GEMS_reanalysis/co_nov_dec_08.nc"
# --------------------------------------------------------------------------------------------------
#  Interface
#  =========
#
#  yr4          YYYY of the receptor time
#  mon          month "
#  day          day   "
#  hr           hour  "
#  co2inifile   netCDF file name (absolute path) of the GEMS fields
#  result       matrix with columns "btime", "lat", "lon" and "pres"
#               gets the column "co2ini" attached/overwritten on output
#  result.sel   selector of same length as the columns of "result": which rows should be considered
#
# --------------------------------------------------------------------------------------------------
#  $Id: get.TM3.netcdf.r,v 1.1 2009-02-24 14:27:54 gerbig Exp $
# --------------------------------------------------------------------------------------------------


   if (!is.element("package:ncdf", search())) library("ncdf") # Load ncdf library
   refdate <- month.day.year(floor(julian(mon, day, yr4) + hr/24 - result[result.sel, "btime"]/24))

#

   tm3file <- open.ncdf(co2inifile, write=F)
   centerlats <- get.var.ncdf(tm3file, varid="latitude")           # Center (assumed here, not sure!)
   centerlons <- get.var.ncdf(tm3file, varid="longitude")           # Center
   times <- get.var.ncdf(tm3file, varid="time")           # Center (instantaneous)
   levs <- get.var.ncdf(tm3file, varid="levelist")           # levels
   nlevs<-length(levs)
   ini <- att.get.ncdf(tm3file, varid="co",attname="add_offset")$value           # global CO2 offset
   scl <- att.get.ncdf(tm3file, varid="co",attname="scale_factor")$value           # global CO2 offset

   # Get Ending positions
   p4ct <- result[result.sel, "pres"]
   lat4ct <- result[result.sel, "lat"]
   lon4ct <- result[result.sel, "lon"]
   lon4ct[lon4ct<0] <- lon4ct[lon4ct<0]+360
   ps4ct <- 1013*exp(-result[result.sel,"grdht"]/8000)#get surface pressure using simple scale height

   #-- lat lon pointer
   dlat<-unique(round(diff(centerlats),4))
   dlon<-unique(round(diff(centerlons),4))
   latpt<-round((lat4ct-centerlats[1])/dlat)+1
   lonpt<-round((lon4ct-centerlons[1])/dlon)+1
   point<-cbind(lonpt,latpt) #spatial pointer (2D)

   #-- time pointer
   delt.h<-(times[2]-times[1])
   timept<-round((julian(mon, day, yr4,c(1, 1, 1900)) + hr/24 - result[result.sel, "btime"]/24)*24/delt.h)*delt.h #hours since 1900-01-01 00:00:0.0
   timept<-(timept-times[1])/delt.h+1
   timept[timept<1]<-1 #default to first
   timept[timept>length(times)]<-length(times) #default to last

   tm3boundary <- rep(NA, length(timept))

   #prepare for height pointer
   # read heights from ECM website for 60 levels (http://www.ecmwf.int/products/data/technical/model_levels/model_def_60.html)
   plevs<-c(  0.000,    0.200,    0.384,    0.636,    0.956,    1.345,    1.806,
       2.348,    2.985,    3.740,    4.646,    5.757,    7.132,    8.837,
      10.948,   13.565,   16.806,   20.823,   25.799,   31.964,   39.603,
      49.067,   60.180,   73.066,   87.727,  104.229,  122.614,  142.902,
     165.089,  189.147,  215.025,  242.652,  272.059,  303.217,  336.044,
     370.407,  406.133,  443.009,  480.791,  519.209,  557.973,  596.777,
     635.306,  673.240,  710.263,  746.064,  780.346,  812.830,  843.263,
     871.420,  897.112,  920.189,  940.551,  958.148,  972.987,  985.140,
     994.747, 1002.024, 1007.264, 1010.849, 1013.250)

   # loop over unique end times
   for (time.i in unique(timept)) {
      ctsel <- which(timept == time.i)
      if(length(ctsel)==1){
        #get vertical pointer: 
        #1) use heights from ECM website for 60 levels (in plevs)
        #2) get surface pressure using simple scale height
        psurf<-ps4ct[ctsel]
        #3) apply compression factor (simplified ...)
        pres<-p4ct[ctsel]*(1013.25/psurf);pres[pres>max(plevs)]<-max(plevs);pres[pres<=min(plevs)]<-min(plevs)+0.001
        #now get vertical pointer
        hpt<-cut(pres,breaks=plevs,labels=FALSE)
        hpt[hpt<min(levs)]<-min(levs)
        hpt[hpt>max(levs)]<-max(levs)
        tm3boundary[ctsel] <- get.var.ncdf(tm3file, varid="co", start=c(point[ctsel,],hpt, time.i),count=c(1,1,1,1))
        #Note that scale factor and ini (global offset) has been applied in get.var.ncdf
      } else { #more than one element to extract
        start.rd<-apply(point[ctsel,],2,min)
        count.rd<-apply(point[ctsel,],2,max)-start.rd + c(1,1)
        count.rd[count.rd<2]<-2 #read at least 2 items per dimension to have right overall dimension in object
        #get vertical pointer: 
        #1) use heights from ECM website for 60 levels (in plevs)
        #2) get surface pressure using simple scale height
        psurf<-ps4ct[ctsel]
        #3) apply compression factor (simplified ...)
        pres<-p4ct[ctsel]*(1013.25/psurf);pres[pres>max(plevs)]<-max(plevs);pres[pres<=min(plevs)]<-min(plevs)+0.001
        #now get vertical pointer
        hpt<-cut(pres,breaks=plevs,labels=FALSE)
        hpt[hpt<min(levs)]<-min(levs)
        hpt[hpt>max(levs)]<-max(levs)
        pt3d<-cbind(point[ctsel,],hpt) #turn into 3 D pointer
        start.rd<-apply(pt3d,2,min)
        count.rd<-apply(pt3d,2,max)-start.rd + c(1,1,1)
        count.rd[count.rd<2]<-2 #read at least 2 items per dimension to have right overall dimension in object
        pt3d<-pt3d+matrix(-start.rd+c(1,1,1),nrow=nrow(pt3d),ncol=ncol(pt3d),byrow=TRUE) #offset correction, to extract from small array to be read from ncdf
        tm3boundary[ctsel] <- get.var.ncdf(tm3file, varid="co", start=c(start.rd, time.i),count=c(count.rd,1))[pt3d]
        #Note that scale factor and ini (global offset) has been applied in get.var.ncdf
      }
   }

   close.ncdf(tm3file)
   result[,spec] <- rep(0,nrow(result))
   #unit conversion to ppb
   kgkg2ppb<-28.97/28.01*1E9
   result[result.sel,spec] <- tm3boundary*kgkg2ppb

   return(result)

}
