#***************************************************************************************************
#  read initial values from CarbonTracker file and attach them to the "result" matrix
#***************************************************************************************************

get.LMDZ.netcdf <- function(yr4=NULL, mon=NULL, day=NULL, hr=NULL, co2inifile=NULL,
                             result=NULL, result.sel=NULL) {

# --------------------------------------------------------------------------------------------------
#  Interface
#  =========
#
#  yr4          YYYY of the receptor time
#  mon          month "
#  day          day   "
#  hr           hour  "
#  co2inifile   netCDF file name (absolute path) of the CarbonTracker data
#  result       matrix with columns "btime", "lat", "lon" and "pres"
#               gets the column "co2ini" attached/overwritten on output
#  result.sel   selector of same length as the columns of "result": which rows should be considered
#
# --------------------------------------------------------------------------------------------------
#  $Id: get.LMDZ.netcdf.r,v 1.2 2008-07-03 11:32:20 gerbig Exp $
# --------------------------------------------------------------------------------------------------


   if (!is.element("package:ncdf", search())) library("ncdf") # Load ncdf library
   refdate <- month.day.year(floor(julian(mon, day, yr4) + hr/24 - result[result.sel, "btime"]/24))

#kludge to limit to may-june only
   unimon<-unique(refdate$month)
   unimon[unimon<5]<-5
   unimon[unimon>6]<-6
   monis<-refdate$month; monis[monis<5]<-5; monis[monis>6]<-6
   for(moni in unimon){
      selm<-monis==moni
      resultm<-result[selm,]
      result.selm<-result.sel[selm]
      co2inifilem<-paste(substring(co2inifile,1,nchar(co2inifile)-8),tolower(month.abb[moni]),"05.nc",sep="")
      lmdzfile <- open.ncdf(co2inifilem, write=F)
      centerlats <- get.var.ncdf(lmdzfile, varid="nav_lat")           # Center
      centerlons <- get.var.ncdf(lmdzfile, varid="nav_lon")           # Center

      plevs <- lmdzfile$dim$sig_s$vals

      # Get Ending positions
      p4ct <- resultm[result.selm, "pres"]*100                   # mb to Pa (stilt output in mb; carbontracker uses Pa)
      lat4ct <- resultm[result.selm, "lat"]
      lon4ct <- resultm[result.selm, "lon"]

      #-- lat lon pointer
      getllpt<-function(latlon,lat.mat,lon.mat){
        return(which(
            abs(lat.mat-latlon[1])==min(abs(lat.mat-latlon[1]))&
            abs(lon.mat-latlon[2])==min(abs(lon.mat-latlon[2])),arr.ind=T)[1:2])
      }
      lonpt<-getllpt(c(lat4ct[1],lon4ct[1]),centerlats,centerlons)
      llpt<-apply(cbind(lat4ct,lon4ct),1,getllpt,centerlats,centerlons) #row index is in 1st row, col index is in 2nd row
      #-- pressure pointer
      getpt<-function(pres,pres.prof){
        return(which(abs(pres-pres.prof)==min(abs(pres-pres.prof)))[1])
      }
      ppt<-lapply(p4ct,getpt,plevs) #pressure pointer

      #-- time pointer
      timept<-(julian(mon, day, yr4) + hr/24 - resultm[result.selm, "btime"]/24-julian(moni, 1, 2005))*24+1
      timept[timept<1]<-1 #default to first
      timept[timept>lmdzfile$dim$time_counter$len]<-lmdzfile$dim$time_counter$len #default to last

      point<-cbind(as.numeric(llpt[1,]),as.numeric(llpt[2,]),as.numeric(ppt)) #spatial pointer (3D)

      lmdzboundary <- rep(NA, length(timept))

      # loop over unique end days
      for (time.i in unique(timept)) {
         ctsel <- which(timept == time.i)
         if(length(ctsel)==1){
           lmdzboundary[ctsel] <- get.var.ncdf(lmdzfile, varid="fos", start=c(point[ctsel,], time.i),count=c(1,1,1,1)) +
                                  get.var.ncdf(lmdzfile, varid="sib_hr", start=c(point[ctsel,], time.i),count=c(1,1,1,1)) +
                                  get.var.ncdf(lmdzfile, varid="taka", start=c(point[ctsel,], time.i),count=c(1,1,1,1))
         } else { #more than one element to extract
           start.rd<-apply(point[ctsel,],2,min)
           count.rd<-apply(point[ctsel,],2,max)-start.rd + c(1,1,1)
           count.rd[count.rd<2]<-2 #read at least 2 items per dimension to have right overall dimension in object
           point[ctsel,1]<-point[ctsel,1]-start.rd[1]+1
           point[ctsel,2]<-point[ctsel,2]-start.rd[2]+1
           point[ctsel,3]<-point[ctsel,3]-start.rd[3]+1
           lmdzboundary[ctsel] <- get.var.ncdf(lmdzfile, varid="fos", start=c(start.rd, time.i),count=c(count.rd,1))[point[ctsel,]] +
                                  get.var.ncdf(lmdzfile, varid="sib_hr", start=c(start.rd, time.i),count=c(count.rd,1))[point[ctsel,]] +
                                  get.var.ncdf(lmdzfile, varid="taka", start=c(start.rd, time.i),count=c(count.rd,1))[point[ctsel,]]
           test <- get.var.ncdf(lmdzfile, varid="fos", start=c(start.rd, time.i),count=c(count.rd,1))[point[ctsel,]] +
                                  get.var.ncdf(lmdzfile, varid="sib_hr", start=c(start.rd, time.i),count=c(count.rd,1))[point[ctsel,]] +
                                  get.var.ncdf(lmdzfile, varid="taka", start=c(start.rd, time.i),count=c(count.rd,1))[point[ctsel,]]
         }
      }

      resultm[,"co2ini"] <- rep(0,nrow(resultm))
      resultm[result.selm,"co2ini"] <- lmdzboundary*2.416E6 + 376 #add scaling and offset

      close.ncdf(lmdzfile)
      result[selm,"co2ini"] <- resultm[,"co2ini"]
   } #loop over months

   return(result)

}
