#***************************************************************************************************
# Function that loops over all starting times
#***************************************************************************************************

Trajecmod <- function(partarg=NULL, totpartarg=NULL, nodeoffset=NULL) {

#---------------------------------------------------------------------------------------------------
# Calls 'Trajec' for each starting time
# arguments assigned from call to setStiltparam.r
#
# $Id: Trajecmod.r,v 1.22 2010/02/23 19:24:19 trn Exp $
#---------------------------------------------------------------------------------------------------


# need to assign parameters; also save parameter setting in archive file with date in name
source("setStiltparam.r")
savename <- gsub(" ", ".", date())
savename <- substring(savename,4)
runs.done.dir <- NULL
if (file.exists('./Runs.done')) runs.done.dir <- './Runs.done/'
if (is.null(runs.done.dir) && file.exists(paste(sourcepath,'Runs.done',sep='')))
  runs.done.dir <- paste(sourcepath,'/Runs.done/',sep='')
if (is.null(runs.done.dir) &&
    substring(path, nchar(path)-nchar("Runs.done"), nchar(path)) == "Runs.done/")
  runs.done.dir <- sourcepath
if (!is.null(runs.done.dir)) {
   file.copy("setStiltparam.r",
             paste(runs.done.dir, "setStiltparam", savename, ".r", sep=""),
             overwrite=T)
   cat("Saving copy of setStiltparam.r in ",
       paste(runs.done.dir, "setStiltparam", savename, ".r", sep=""),
       "\n",sep="")
} else {
   cat("No Runs.done directory; parameter information not saved\n")
}



totpart <- 1
if (!is.null(totpartarg)) {
  cat('resetting totpart=', totpart, ' to totpartarg=', totpartarg, '\n', sep='')
  totpart <- totpartarg
}
part <- 1
if (!is.null(partarg)) {
  cat('Using totpart=', totpart, ' resetting part=', part, ' to partarg=', partarg, '\n', sep='')
  part <- partarg
}
if (!is.null(nodeoffset)) {
  nummodel <- part+nodeoffset
  cat('Using nodeoffset= ', nodeoffset, ' ,results in nummmodel= ', nummodel, '\n', sep='')
} else {
  nummodel <- part
}



# get Starting info
if (!existsr(Timesname, path)) stop(paste("cannot find object ", Timesname, " in directory ", path, sep=""))
StartInfo <- getr(paste(Timesname, sep=""), path) # object containing fractional julian day, lat, lon, agl for starting position and time
# SELECTION OF A FEW Receptors for testing!
if (Times.startrow > 0) StartInfo <- StartInfo[Times.startrow:Times.endrow,, drop=FALSE] # can be just one (Times.startrow=Times.endrow)

# divide job into "totpart" parts to speed up
if (dim(StartInfo)[1] < totpart) {
  cat ('Warning: resetting totpart=', totpart, ' to dim(StartInfo)[1]=', dim(StartInfo)[1], '\n', sep='')
  totpart <- dim(StartInfo)[1]
}
if (part > totpart) {
  stop.message <- paste('Specified part=', part, ' > totpart=', totpart, ', stopping\n')
  cat(stop.message)
  stop(stop.message)
}
start.rows <- c(1 + (0:(totpart-1))*floor(dim(StartInfo)[1]/totpart), dim(StartInfo)[1]+1)
StartInfo <- StartInfo[start.rows[part]:(start.rows[part+1]-1),, drop=FALSE]
dimnames(StartInfo)[[2]] <- toupper(dimnames(StartInfo)[[2]])

if (biomassburnTF) {
   biomassburnoutmat <- matrix(nrow=dim(StartInfo)[[1]], ncol=2)
   dimnames(biomassburnoutmat) <- list(NULL, c("ident", "CO"))
}

# OVERWRITE WARNING
if(existsr(paste("stiltresult",part,sep=""),path=path)) {
   warning("You are attempting to overwrite an existing stiltresult object")
   warning("Notice: If you have changed parameters and Trajecmod fails, first try to move or remove the existing stiltresult object")
}

nrows <- length(StartInfo[,1]) # 1 row for each trajectory
rownum <- 1
firsttraj <- T
firstflux <- T
l.remove.Trajfile <- FALSE
if (exists('remove.Trajfile')) l.remove.Trajfile <- remove.Trajfile
for (j in 1:nrows) {

  ###############################################
  ##### run trajectories and save output ########
  ###############################################
  lat <- StartInfo[j, "LAT"]; lon <- StartInfo[j, "LON"]; agl <- StartInfo[j, "AGL"]
  identname <- pos2id(StartInfo[j,1], lat, lon, agl)
  cat("Trajecmod(): ", identname, " running at ", date(), "\n", sep="")
  dat <- month.day.year(floor(StartInfo[j,1])) # from julian to mmddyy
  yr4 <- dat$year # 4 digit year
  yr <- yr4%%100 # 2 digit year (or 1 digit...)
  mon <- dat$month
  day <- dat$day
  hr <- round((StartInfo[j,1]-floor(StartInfo[j,1]))*24)
  l.ziscale <- NULL
  if (exists('ziscale')) l.ziscale <- ziscale
  l.zsg.name <- NULL
  if (exists('zsg.name')) l.zsg.name <- zsg.name
  l.create.X0 <- FALSE
  if (exists('create.X0')) l.create.X0 <- create.X0
  l.use.multi <- TRUE
  if (exists('use.multi')) l.use.multi <- use.multi
  l.hymodelc.exe <- NULL
  if (exists('hymodelc.exe')) l.hymodelc.exe <- hymodelc.exe
  if (l.use.multi) {
    l.setup.list <- list()
    if (exists('setup.list')) l.setup.list <- setup.list
    info <- Trajecmulti(yr=yr, mon=mon, day=day, hr=hr, lat=lat, lon=lon, agl=agl, nhrs=nhrs,
                   numpar=nparstilt, doublefiles=T, metd=metsource, metlib=metpath,
                   conv=convect, overwrite=overwrite, outpath=path, varsout=varstrajec, rundir=rundir,
                   nummodel=nummodel, sourcepath=sourcepath, ziscale=l.ziscale, zsg.name=l.zsg.name,
                   create.X0=l.create.X0,setup.list=l.setup.list,hymodelc.exe=l.hymodelc.exe)
  } else {
    info <- Trajec(yr=yr, mon=mon, day=day, hr=hr, lat=lat, lon=lon, agl=agl, nhrs=nhrs,
                   maxdist=stepsize, numpar=nparstilt, doublefiles=T, metd=metsource, metlib=metpath,
                   conv=convect, overwrite=overwrite, outpath=path, varsout=varstrajec, rundir=rundir,
                   nummodel=nummodel, sourcepath=sourcepath, ziscale=l.ziscale, zsg.name=l.zsg.name,
                   create.X0=l.create.X0)
  }
  if (firsttraj) { # set up array for run info
    run.info <- matrix(NA, nrow=nrows, ncol=length(info))
    dimnames(run.info) <- list(NULL, names(info))
    firsttraj <- F
  } else {
    havenames <- dimnames(run.info)[[2]]
    for (nm in names(info)) {
      ndx <- which(havenames == nm)[1] # check if there are new column names
      if (is.na(ndx)) { # new column name, need to add column
        run.info <- cbind(run.info, rep(NA, dim(run.info)[1]))
        dimnames(run.info)[[2]][dim(run.info)[2]] <- nm
      } else { havenames[ndx] <- paste(havenames[ndx], "done", sep="")} # dummy
    }
  }
  run.info[j, names(info)] <- info

  #########################################################################
  ###### TM3-STILT ################
  if(writeBinary == T){

    ####### variables for writing to binary files ######
    dlat   <-as.numeric(lat,digits=10)
    dlon   <-as.numeric(lon,digits=10)
    dagl   <-as.numeric(agl,digits=10)
    dlatres<-as.numeric(lat.res,digits=10)
    dlonres<-as.numeric(lon.res,digits=10)

    ####### construct filename for binary file #########
    cyr<-as.character(2000+yr)
    cmon<-as.character(mon)
    cday<-as.character(day)
    chr<-as.character(hr)
    x1<-""; x2<-""; x3<-""
    if(mon<10) x1<-paste(x1,"0",sep="")
    if(day<10) x2<-paste(x2,"0",sep="")
    if(hr<10)  x3<-paste(x3,"0",sep="")

    pathBinFootprintstation<-paste(pathBinFootprint,station,"/",sep="")
    if (file.access(pathBinFootprintstation,0)!=0) {
     system(paste("mkdir ",pathBinFootprintstation,sep=""))
     }

    filename<-paste(pathBinFootprintstation,station,"_",as.character(-nhrs),"h0",as.character(ftintr),
                    "h_",cyr,x1,cmon,x2,cday,x3,chr,"_",gridtag,cendian,".d",sep="")
# Jan's version with height in filename
#    filename<-paste(pathBinFootprintstation,station,"_",as.character(-nhrs),"h0",as.character(ftintr),
#                    "h_",cyr,x1,cmon,x2,cday,x3,chr,"_",sprintf("%5.5d",as.integer(agl)),"_",gridtag,cendian,".d",sep="")

    if (file.exists(filename)) {
    print(paste("Binary footprint file ",filename," already exists"))
    print(paste("not replaced !!"))
    }else{

    ident<-info["outname"]
    print(paste(path,".RData",ident,sep=""))
    if (file.exists(paste(path,".RData",ident,sep=""))) {

    #for longer than hourly intervals for footprints, first make sure to match time intervals of flux fields
    #assume those are e.g. 0-3, 3-6 etc. UTC, or 0-24 UTC 
    #NOTE: only for hourly intervals variables "foottimes", "nfoottimes" and "nftpix" are computed in setStiltparam.r 
    if(ftintr>1){	
      nfoottimes <- -nhrs/ftintr+2               #number of footprints computed
      foottimes<-rep(c(0),nfoottimes)            #vector of times (backtimes) in hours between which footprint is computed
      nftpix<-rep(c(0),nfoottimes)               #vector of numbers of pixels in each footprint
      for(ft in 2:nfoottimes){ 
        foottimes[ft]<- hr+(ft-2)*ftintr 
      }
      foottimes[nfoottimes]<- -nhrs
      if(hr==0){				 #special case when starting at midnight
        foottimes<-foottimes[2:nfoottimes]
        nftpix<-nftpix[2:nfoottimes]
        nfoottimes<-nfoottimes-1
      }
    }

    ####### call Trajecfoot ############################ 
    ident<-info["outname"]
     foot <- Trajecfoot(ident, pathname=path, foottimes=foottimes, zlim=c(zbot, ztop),fluxweighting=NULL, coarse=1, vegpath=vegpath,
                   numpix.x=numpix.x, numpix.y=numpix.y,lon.ll=lon.ll, lat.ll=lat.ll, lon.res=lon.res, lat.res=lat.res)
#    foot<-Trajecfoot(ident=ident,pathname=path,foottimes=foottimes,zlim=c(zbot,ztop),fluxweighting=NULL,coarse=1,vegpath=vegpath,
#                  numpix.x=numpix.x,numpix.y=numpix.y,lon.ll=lon.ll,lat.ll=lat.ll,lon.res=lon.res,lat.res=lat.res)
    nameslat<-rownames(foot)
    nameslon<-colnames(foot)
#    print(paste("nameslat, nameslon: ",nameslat, nameslon))
    if(is.null(foot)){
      print(paste("is.null(foot): ",is.null(foot),foot))
      print(paste("No binary footprint file for TM3 written!!!!"))

    }else{
 
#    #### write the output file ####
#    cyr<-as.character(2000+yr)
#    cmon<-as.character(mon)
#    cday<-as.character(day)
#    chr<-as.character(hr)
#    x1<-""; x2<-""; x3<-""
#    if(mon<10) x1<-paste(x1,"0",sep="")
#    if(day<10) x2<-paste(x2,"0",sep="")
#    if(hr<10)  x3<-paste(x3,"0",sep="")
#
#    pathBinFootprintstation<-paste(pathBinFootprint,station,"/",sep="")
#    if (file.access(pathBinFootprintstation,0)!=0) {
#     system(paste("mkdir ",pathBinFootprintstation,sep=""))
#     }
#
##    filename<-paste(pathBinFootprintstation,station,"_",as.character(-nhrs),"h0",as.character(ftintr),
##                    "h_",cyr,x1,cmon,x2,cday,x3,chr,"_",gridtag,cendian,".d",sep="")
## Jan's version with height in filename
#    filename<-paste(pathBinFootprintstation,station,"_",as.character(-nhrs),"h0",as.character(ftintr),
#                    "h_",cyr,x1,cmon,x2,cday,x3,chr,"_",sprintf("%5.5d",as.integer(agl)),"_",gridtag,cendian,".d",sep="")

    print(paste("Writing binary footprint file: ",filename))
    con<-file(filename, "wb")
# write first: hours backward, hours interval, lat, lon, altitude, lat resolution, and lon resolution
    writeBin(as.numeric(c(nhrs,ftintr,lat,lon,agl,lat.res,lon.res)),con,size=4,endian=endian) 

    nftpix<-apply(foot,c(3),function(x)return(sum(x>0)))
    ftmax<-which(nftpix>0)[sum(nftpix>0)]
    for(ft in 1:(nfoottimes-1)){
      writeBin(as.numeric(c(ft,nftpix[ft])),con,size=4,endian=endian)          # index of footprint and # grids in footprint
      #### pixel by pixel #####
      if(nftpix[ft]>0){
        id<-which(t(foot[,,ft])>0,arr.ind=T)
        foot_data<-cbind(as.numeric(nameslat[id[,2]]),as.numeric(nameslon[id[,1]]),foot[,,ft][foot[,,ft]>0])# save coordinates of pixels of each footprint and footprint value
        writeBin(as.vector(t(foot_data)),con,size=4,endian=endian)
        if(any(foot_data[,3]<0))browser()
      } #if(nftpix[ft]>0)
    } #for(ft   
     print(paste("ftmax: ",ftmax))
 
   close(con)
   } #if is.null(foot) 

   }else{
   print(paste(path,outname," does exist -> no new STILT run !!",sep=""))
   } #if (file.exists(paste(path,outname,sep=""))) 

   } #if (file.exists(filename))


  }  #end if(writeBinary == T)

  ###### end of TM3-STILT output ################
  #########################################################################

  #########################################################################
  ##### map trajectories to flux grids and vegetation maps ################
  ##### calculate mixing ratios at receptor points, save in result ########
  if (fluxTF) {
     print(paste("Trajecmod(): rownumber j:", j))


     traj <- Trajecvprm(ident=identname, pathname=path, tracers=fluxtracers, coarse=aggregation,
                dmassTF=T, nhrs=nhrs, vegpath=vegpath, evilswipath=evilswipath,
                vprmconstantspath=vprmconstantspath, vprmconstantsname=vprmconstantsname, nldaspath=nldaspath,
                nldasrad=usenldasrad, nldastemp=usenldastemp, pre2004=pre2004,
                keepevimaps=keepevimaps, detailsTF=detailsTF, bios=fluxmod, landcov=landcov,
                numpix.x=numpix.x, numpix.y=numpix.y, lon.ll=lon.ll, lat.ll=lat.ll,
                lon.res=lon.res, lat.res=lat.res)


     # 'traj' is a vector
     if (existsr(paste("stiltresult", part, sep=""), path=path)) {
        result <- getr(paste("stiltresult", part, sep=""), path=path)
        if (dim(result)[1] != nrows) {
           if (firstflux) print("Trajecmod(): existing stiltresult has wrong dimension; creating new one.")
        } else {
           if (firstflux) print("Trajecmod(): found existing stiltresult, update rows in that.")
           firstflux <- FALSE
        }
     }
     if (firstflux) { # at beginning create result object
        ncols <- length(traj) # all from Trajec(), + 3 from StartInfo (agl, lat, lon)
        result <- matrix(NA, nrow=nrows, ncol=ncols)
        firstflux <- F
     }
     result[rownum, ] <- traj
     dimnames(result) <- list(NULL, c(names(traj)))
     dimnames(result) <- list(NULL, dimnames(result)[[2]])
     # write the object into default database; object names are, e.g., "Crystal.1"
     assignr(paste("stiltresult", part, sep=""), result, path=path)
  }
  rownum <- rownum+1

  ##### calculate footprint, assign in object ########
  if (footprintTF) {
     print(paste("Trajecmod(): ", identname, " running footprint at ", unix("date"), sep=""))
     print(paste("Trajecmod(): memory in use:", memory.size()[1]))
     foot <- Trajecfoot(identname, pathname=path, foottimes=foottimes, zlim=c(zbot, ztop),
                        fluxweighting=NULL, coarse=1, vegpath=vegpath,
                        numpix.x=numpix.x, numpix.y=numpix.y,
                        lon.ll=lon.ll, lat.ll=lat.ll, lon.res=lon.res, lat.res=lat.res)
     assignr(paste("foot", identname, sep=""), foot, path)
     print(paste("Trajecmod(): foot", identname, " assigned", sep=""))
  } # if (exists(foottimes))

  ##### plot footprint ########
  if (footplotTF) { # plot footprints
    foot <- getr(paste("foot", identname, sep=""), path)
    footplot(foot,identname,lon.ll,lat.ll,lon.res,lat.res)
    for (foottimespos in 1:(length(foottimes)-1)) {
    ############# NOT READY YET ###############
    }
  }

  # Specify the function parameters
  if (biomassburnTF) {
     biomassburnoutmat[j, ] <- biomassburn(timesname=StartInfo, burnpath=burnpath, endpath=path, pathname=path, nhrs=nhrs, timesrow=j)
     print(paste("Biomassburning influence calculated to be ", biomassburnoutmat[j,2], " ppbv. Inserted into fireinfluence matrix row ", j, sep=""))
  }

  if(l.remove.Trajfile)unix(paste("rm -f ",paste(path,".RData",identname,sep=""),sep=""))

}                                                           # for (j in 1:nrows)


# Wrap up all of the CO biomassburning calculations
if (biomassburnTF)
  write.table(biomassburnoutmat, file=paste(path, "fireinfluencex", nhrs, "hr_", part, ".txt", sep=""), row.names=F)

##### save mixing ratios at receptor points in file, e.g. stiltresult1.csv for part=1 ########
if (fluxTF) {
   dimnames(result) <- list(NULL, dimnames(result)[[2]])
   # write the object into default database; object names are, e.g., "Crystal.1"
   assignr(paste("stiltresult", part, sep=""), result, path=path)
   print(paste("stiltresult", part, " assigned in ", path, sep=""))
   write.table(result, file=paste(path, "stiltresult", part, ".csv", sep=""), na="", row.names=F)
}

# If evi and lswi maps from vprm calculations is saved to the global environment; it should be removed here

rm(list=objects(pattern="GlobalEvi"), envir=globalenv())
rm(list=objects(pattern="GlobalLswi"), envir=globalenv())
gc(verbose=F)

return(run.info)

}

