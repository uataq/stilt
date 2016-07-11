#***************************************************************************************************
#  read initial values from CarbonTracker file and attach them to the "result" matrix
#***************************************************************************************************

get.CarbonTracker.netcdf <- function(yr4=NULL, mon=NULL, day=NULL, hr=NULL, co2inifile=NULL,
                             result, result.sel) {

# --------------------------------------------------------------------------------------------------
#  Interface
#  =========
#
#  yr4          YYYY of the receptor time
#  mon          month "
#  day          day   "
#  hr           hour  "
#  co2inifile   netCDF file name (absolute path) of the CarbonTracker data
#               However, for 2009 version ("CT2009"), give only the FIRST PART of file (e.g., "CT2009.molefrac_glb3x2_")--since files are no longer monthly.
#               The code would then decide smartly exactly which daily file to use
#  result       matrix with columns "btime", "lat", "lon" and "pres"
#               gets the column "co2ini" attached/overwritten on output
#  result.sel   selector of same length as the columns of "result": which rows should be considered
#
# --------------------------------------------------------------------------------------------------
#  $Id: get.CarbonTracker.netcdf.r,v 1.5 2009/03/04 22:13:09 johnlin Exp $
# --------------------------------------------------------------------------------------------------


#JCL(100221): check if these are newer versions of CarbonTracker (CT 2007B or 2008 or 2009)
newTF<-length(grep("CT2008",co2inifile))==1|length(grep("CT2009",co2inifile))==1
#if(length(grep("CT2007B",co2inifile))==1|length(grep("CT2008",co2inifile))==1|length(grep("CT2009",co2inifile))==1){ 
if(length(grep("CT2007B",co2inifile))==1|newTF){
   fjday<-julian(mon, day, yr4) + hr/24 - result[result.sel, "btime"]/24  
   mons<-month.day.year(floor(fjday))$month;days<-month.day.year(floor(fjday))$day
   datecodes.all<-paste(sprintf('%2.2d',mons),sprintf('%2.2d',days))
   datecodes<-unique(datecodes.all)
   for(i in 1:length(datecodes)){
      datecode<-datecodes[i]
      moni<-as.numeric(substring(datecode,1,2))
      #JCL(100222):  deal with multiple trajectory end dates
      #selm<-month.day.year(floor(fjday))$month==moni
      selm<-datecodes.all==datecode
      resultm<-result[result.sel,][selm,]
      if(sum(selm)==1)resultm<-t(as.matrix(resultm))

      co2inifilei<-paste(co2inifile,yr4,'-',sprintf('%2.2d',moni),'.nc',sep="")
      #newest version:  netCDF files are no longer monthly, but daily
      if(length(grep("CT2009",co2inifile))==1){
        #dayi<-unique(month.day.year(floor(fjday))$day)
        dayi<-as.numeric(substring(datecode,4,5))
        co2inifilei<-paste(co2inifile,yr4,'-',sprintf('%2.2d',moni),'-',sprintf('%2.2d',dayi),'.nc',sep="")
      } #if(length(grep("CT2009",co2inifile))==1){

      # CarbonTracker notes
      if (!is.element("package:ncdf", search())) library("ncdf") # Load ncdf library
      print(paste("Opening CarbonTracker file:",co2inifilei))
      trackerfile <- open.ncdf(co2inifilei, write=F)
      ctcenterlats <- trackerfile$dim$lat$vals           # Center
      ctcenterlons <- trackerfile$dim$lon$vals           # Center
      ctlatres <- unique(diff(ctcenterlats))
      ctlonres <- unique(diff(ctcenterlons))
      ctnumlat <- length(ctcenterlats)
      ctnumlon <- length(ctcenterlons)
      ctlowerleftlats <- as.numeric(ctcenterlats)-0.5*ctlatres
      ctlowerleftlons <- as.numeric(ctcenterlons)-0.5*ctlonres
      level <- get.var.ncdf(trackerfile, varid="level")   
      ctnumlevs <- length(level)                     # Coefficients denote mid-points

      # Get Ending positions
      p4ct <- resultm[, "pres"]*100                   # mb to Pa (stilt output in mb; carbontracker uses Pa)
      lat4ct <- resultm[, "lat"]
      lon4ct <- resultm[, "lon"]

      #-- check reference dates for CT access and reset year to CT range, if necessary

      tfjday<-trackerfile$dim$time$vals+julian(1, 1, 2000)
      #JCL(100221):  also allow for 2009 updates; these new files use "date" rather than "time"
      #if(length(grep("CT2008",co2inifile))==1){tfjday<-trackerfile$dim$date$vals+julian(1, 1, 2000)}
      if(newTF){tfjday<-trackerfile$dim$date$vals+julian(1, 1, 2000)}
      tpointer <- round((fjday[selm]-tfjday[1])/unique(diff(tfjday))+1+1E-9) #1E-9 to make proper round ...
      ctboundary <- rep(NA, length(tpointer))


      # loop over unique end days
      for (jjct in unique(tpointer)) {

         ctsel <- which(tpointer == jjct)
         sublat4ct <- lat4ct[ctsel]
         sublon4ct <- lon4ct[ctsel]
         subp4ct <- p4ct[ctsel]

         trackerpress <- get.var.ncdf(trackerfile, varid="press", start=c(1,1,1,jjct),
                                      count=c(ctnumlon, ctnumlat,ctnumlevs,1))
         # Make Pointer
         ctlatlimlow <- min(ctlowerleftlats)
         ctlatlimhi <- max(ctlowerleftlats)+unique(diff(ctlowerleftlats))
         ctlonlimlow <- min(ctlowerleftlons)
         ctlonlimhi <- max(ctlowerleftlons)+unique(diff(ctlowerleftlons))

         ctlatpointer <- findInterval(sublat4ct, c(ctlowerleftlats, ctlatlimhi), rightmost.closed=TRUE,all.inside=TRUE)
         ctlonpointer <- findInterval(sublon4ct, c(ctlowerleftlons, ctlonlimhi), rightmost.closed=TRUE,all.inside=TRUE)

         # Vertical pressure pointer will vary with surface pressure due to hybrid-sigma coordinates, loop through
         ctprespointer <- rep(0, length(ctlatpointer))
         for (ppct in 1:length(subp4ct)) {
            trackerpresbreaks <- trackerpress[ctlonpointer[ppct], ctlatpointer[ppct],]
            if (subp4ct[ppct]>max(trackerpresbreaks))
               ctprespointer[ppct] <- 1
            else
               ctprespointer[ppct] <- findInterval(-subp4ct[ppct],-trackerpresbreaks,
                                                   rightmost.closed=TRUE,all.inside=TRUE)
         } # for ppct

         ctpointer <- cbind(ctlonpointer, ctlatpointer, ctprespointer)

         # Get CO2
         trackerco2.bg    <- get.var.ncdf(trackerfile, varid="bg", start=c(1,1,1, jjct),
                                    count=c(ctnumlon, ctnumlat, ctnumlevs,1))
         trackerco2.bio   <- get.var.ncdf(trackerfile, varid="bio", start=c(1,1,1, jjct),
                                    count=c(ctnumlon, ctnumlat, ctnumlevs,1))
         trackerco2.ff    <- get.var.ncdf(trackerfile, varid="ff", start=c(1,1,1, jjct),
                                    count=c(ctnumlon, ctnumlat, ctnumlevs,1))
         trackerco2.fires <- get.var.ncdf(trackerfile, varid="fires", start=c(1,1,1, jjct),
                                    count=c(ctnumlon, ctnumlat, ctnumlevs,1))
         trackerco2.ocean <- get.var.ncdf(trackerfile, varid="ocean", start=c(1,1,1, jjct),
                                    count=c(ctnumlon, ctnumlat, ctnumlevs,1))
         trackerco2 <- trackerco2.bg + trackerco2.bio + trackerco2.ff + trackerco2.fires + trackerco2.ocean

         # Select CO2 from array
         ctboundary[ctsel] <- trackerco2[ctpointer]
      }

      resultm[,"co2ini"] <- ctboundary

      close.ncdf(trackerfile)
      result[result.sel,][selm,"co2ini"] <- resultm[,"co2ini"]
   } #loop over months

   return(result)

} else { #use carbon tracker 2007 a version

   # CarbonTracker notes
   if (!is.element("package:ncdf", search())) library("ncdf") # Load ncdf library
   trackerfile <- open.ncdf(co2inifile, write=F)
   ctcenterlats <- trackerfile$dim$latitude$vals            # Center
   ctcenterlons <- trackerfile$dim$longitude$vals           # Center
   ctlatres <- unique(diff(ctcenterlats))
   ctlonres <- unique(diff(ctcenterlons))
   ctnumlat <- length(ctcenterlats)
   ctnumlon <- length(ctcenterlons)
   ctlowerleftlats <- as.numeric(ctcenterlats)-0.5*ctlatres
   ctlowerleftlons <- as.numeric(ctcenterlons)-0.5*ctlonres
   hysig.coeffa <- get.var.ncdf(trackerfile, varid="at")    # hybrid sigma coordinate coefficients
   hysig.coeffb <- get.var.ncdf(trackerfile, varid="bt")
   ctnumlevs <- length(hysig.coeffa)-1                      # Coefficients denote breaks
   ctStartDate  <- as.integer(att.get.ncdf(trackerfile,varid=0,attname="StartDate")$value) # expect YYYYMMDD
   ctStartYear  <- ctStartDate %/% 10000
   ctStartMonth <- (ctStartDate %/% 100) %% 100
   ctStartDay   <- ctStartDate %% 100
   ctEndDate    <- as.integer(att.get.ncdf(trackerfile,varid=0,attname="EndDate")$value)
   ctEndYear    <- ctEndDate %/% 10000

   # Get Ending positions
   p4ct <- result[result.sel, "pres"]*100                   # mb to Pa (stilt output in mb; carbontracker uses Pa)
   lat4ct <- result[result.sel, "lat"]
   lon4ct <- result[result.sel, "lon"]


   #-- check reference dates for CT access and reset year to CT range, if necessary

   refdate <- month.day.year(julian(mon, day, yr4) + hr/24 - result[result.sel, "btime"]/24)
   refdate$year <- pmin(pmax(refdate$year,ctStartYear), ctEndYear)
   if (trackerfile$dim$date$units != "days since 2000-01-01")
      stop("getCarbonTracker: wrong date unit in CarbonTracker file. Stop.") # cttimeindices formula wrong
   cttimeindices <- floor(julian(refdate$month, refdate$day, refdate$year)
                          - julian(ctStartMonth,ctStartDay,ctStartYear)) + 1
   ctboundary <- rep(NA, length(cttimeindices))


   # loop over unique end days
   for (jjct in unique(cttimeindices)) {

      ctsel <- which(cttimeindices == jjct)
      sublat4ct <- lat4ct[ctsel]
      sublon4ct <- lon4ct[ctsel]
      subp4ct <- p4ct[ctsel]

      trackersurfp <- get.var.ncdf(trackerfile, varid="psurf", start=c(1,1, jjct),
                                   count=c(ctnumlon, ctnumlat,1))

      # Make Pointer
      ctlatlimlow <- min(ctlowerleftlats)
      ctlatlimhi <- max(ctlowerleftlats)+unique(diff(ctlowerleftlats))
      ctlonlimlow <- min(ctlowerleftlons)
      ctlonlimhi <- max(ctlowerleftlons)+unique(diff(ctlowerleftlons))

      ctlatpointer <- findInterval(sublat4ct, c(ctlowerleftlats, ctlatlimhi), rightmost.closed=TRUE,all.inside=TRUE)
      ctlonpointer <- findInterval(sublon4ct, c(ctlowerleftlons, ctlonlimhi), rightmost.closed=TRUE,all.inside=TRUE)

      # Vertical pressure pointer will vary with surface pressure due to hybrid-sigma coordinates, loop through
      ctprespointer <- rep(0, length(ctlatpointer))
      for (ppct in 1:length(subp4ct)) {
         trackerpresbreaks <- hysig.coeffa+hysig.coeffb*trackersurfp[ctlonpointer[ppct], ctlatpointer[ppct]]
         if (subp4ct[ppct]>max(trackerpresbreaks))
            ctprespointer[ppct] <- 1
         else
            ctprespointer[ppct] <- findInterval(-subp4ct[ppct],-trackerpresbreaks,
                                                rightmost.closed=TRUE,all.inside=TRUE)
      } # for ppct

      ctpointer <- cbind(ctlonpointer, ctlatpointer, ctprespointer)

      # Get CO2
      trackerco2 <- get.var.ncdf(trackerfile, varid="co2mf", start=c(1,1,1, jjct),
                                 count=c(ctnumlon, ctnumlat, ctnumlevs,1))

      # Select CO2 from array
      ctboundary[ctsel] <- trackerco2[ctpointer]
   }

   result[,"co2ini"] <- rep(0,nrow(result))
   result[result.sel,"co2ini"] <- ctboundary

   close.ncdf(trackerfile)

   return(result)


}
}

