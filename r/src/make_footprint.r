Trajecfoot<-function(ident,part=NULL,timelabel=NULL,pathname="",foottimes=c(0,360),zlim=c(0,0),fluxweighting=NULL,coarse=1,dmassTF=T,
                     vegpath="/deas/group/stilt/Vegetation/",numpix.x=376,numpix.y=324,lon.ll=-145,lat.ll=11,lon.res=1/4,lat.res=1/6,landcov="IGBP",wrfinput=NULL){
  #Creates footprint for individual particle runs
  #footrpints are in units of ppm/(micro-moles/m2/s)
  #'ident' is character value specifying the trajectory ensemble to look at 
  #'part' is the object containing particle run; if NULL, then determine from ident and pathname
  #'timelabel' has format of "2002x08x03x10x"; will be used as timestamp for receptor; otherwise, take from 'ident'
  #'pathname' is path where object with particle locations is saved
  #'foottimes' is vector of times between which footprint or influence will be integrated
  #'zlim' (if not default 0,0): vertical interval over which particle distribution is integrated (partial column integrated particle density, or "volume influence")
  #                             default (0,0): footprint calculation ("surface influence")
  #'fluxweighting' if not NULL, but a number (1-31), weighting by that vegetation class;
  #                if not NULL, but a "CO" or "CO2" weighting by these emissions
  #'coarse' degrade resolution (for aggregation error): 0: only 20 km resolution; 
  #        1-16: dynamic resolution, but limited highest resolution
  #	coarse:	  (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
  #	coarsex:c(1,1,2,2,2,4,4,4,8, 8, 8,16,16,16,32,32) #factors by which grids have been made coarser
  #	coarsey:c(1,2,1,2,4,2,4,8,4, 8,16, 8,16,32,16,32)#e.g., '4' means grid is 4 times coarser
  #'dmassTF' if TRUE, weighting by accumulated mass due to violation of mass conservation in met fields
  #'vegpath' is path under wich vegetation and flux grids are stored
  #'numpix.x'			#number of pixels in x directions in grid
  #'numpix.y' 			#number of pixels in y directions in grid
  #'lon.ll'			#lower left corner of grid (longitude of southwest corner of southwest corner gridcell)
  #'lat.ll'			#lower left corner of grid (latitude of southwest corner of southwest corner gridcell)
  #'lon.res'			#resolution in degrees longitude
  #'lat.res'			#resolution in degrees latitude
  #'wrfinput'			#for projection of footprints on WRF grid: absoluzt path to WRF input file (incl. filename)
  #
  #  $Id: Trajecfoot.r,v 1.17 2009/02/11 13:03:13 gerbig Exp $
  #---------------------------------------------------------------------------------------------------
  
  #JCL:(04/19/2004) enable capability to directly pass on the object
  if(is.null(part)){   #when object wasn't passed on
    #Check if object exists
    if(length(unix.shell(paste("if (-e ",pathname,".RData",ident,") echo 'found it!'; endif",sep=""),shell="/bin/csh"))>0){
      print(paste("Trajecfoot(): starting with ident=",ident,sep="")) #found object
      part<-getr(ident,pathname) #get it
    }else{#if not there, break out of function, return NA's
      print(paste("Trajecfoot(): object=",pathname,".RData",ident," NOT FOUND",sep=""))
      return()
    } #if exists or not
  } #if(!is.null(part)){
  
  #check if required particle information is available
  rqdnames<-c("time","lat","lon","agl","zi","index","foot")
  if(dmassTF)rqdnames<-c(rqdnames,"dmass")
  if(any(!(rqdnames%in%dimnames(part)[[2]]))){print(paste("not all columns available for this run, returning NA",sep=""));return(NA)}
  if(!is.null(fluxweighting)){
    if(landcov=="IGBP") {
      veghead<-"veg."
    } else if(landcov=="GLCC") {
      veghead<-"glcc."
    } else if(landcov=="DVN") {
      veghead<-"devanveg."
    } else if(landcov=="SYNMAP") {
      veghead<-"synmap."
    } else if(landcov=="SYNMAP.VPRM8") {
      veghead<-"synvprm8."
    } else {
      stop("Improperly specified Landcover format; Exiting now!")
      veghead<-"veg."
    }
    if(fluxweighting=="CO2"|fluxweighting=="CO"){
      if(is.null(timelabel)){ #if timestamp not passed on, then take from 'ident'
        if(nchar(ident)!=34){stop("need timelabel or proper name (ident) including time stamp for weighting by emissions in Trajecfoot")} 
        pos<-id2pos(ident)
        time<-month.day.year(floor(pos[1]))
      }  
      if(!is.null(timelabel)){pos<-id2pos(timelabel);time<-month.day.year(floor(pos[1]))}  
      yr4<-time$year #4 digit year
      mon<-time$month
      day<-time$day
      hr<-round((pos[1]-floor(pos[1]))*24)
    } #if flux weighting by CO2 or CO fossil fuel emissions
  } #if weighting by flux
  
  #get grid indices
  #For horizontal grids (lower left corner of south-west gridcell: lat.ll, lon.ll; resolution: lon.res, lat.res, x is lon, y is lat)
  if(is.null(wrfinput)){
    #adapt for domain that crosses dateline
    if(lon.ll+numpix.x*lon.res>180)part[part[,"lon"]<0,"lon"]<-part[part[,"lon"]<0,"lon"]+360
    gitx<-floor(1/lon.res*(part[,"lon"]-lon.ll)+1)
    gity<-floor(1/lat.res*(part[,"lat"]-lat.ll)+1)
  } else { #use wrfinput to get projection info
    require(RNetCDF)
    nc      <- open.nc(wrfinput)
    latx    <- var.get.nc(nc,"XLAT_V") #diff columns 1-ny = N-S, diff rows 1-nx = W-E
    longx   <- var.get.nc(nc,"XLONG_U") #using staggered grid, to get boundaries between cells
    moadstandlat_1 <- att.get.nc(nc,"NC_GLOBAL","TRUELAT1")
    moadstandlat_2 <- att.get.nc(nc,"NC_GLOBAL","TRUELAT2")
    moadstandlon   <- att.get.nc(nc,"NC_GLOBAL","STAND_LON")    
    moadknownlat   <- att.get.nc(nc,"NC_GLOBAL","MOAD_CEN_LAT")
    moadknownlon   <- att.get.nc(nc,"NC_GLOBAL","CEN_LON")
    dx     <- att.get.nc(nc,"NC_GLOBAL","DX")
    dy     <- att.get.nc(nc,"NC_GLOBAL","DY")
    close.nc(nc)
    
    projout<-paste("+proj=lcc +lat_1=",moadstandlat_1," +lat_2=",moadstandlat_2," +lat_0=",moadknownlat," +lon_0=",moadstandlon,sep="")
    xy.corner<-project.nice(x.in=longx[1,1], y.in=latx[1,1],proj=projout,inv=FALSE,matTF=FALSE) #inv needs to be false for ll to xy
    xy<-  project.nice(x.in=part[,"lon"], y.in=part[,"lat"],proj=projout,inv=FALSE,matTF=FALSE,matinTF=FALSE) #using staggered grids
    gitx<-floor((xy[[1]]-xy.corner[[1]])/dx+1)
    gity<-floor((xy[[2]]-xy.corner[[2]])/dy+1)
    #set/reset grid definition
    numpix.y<-dim(longx)[2] #staggered, so use lon (or drop last column of lat)
    numpix.x<-dim(latx)[1]
  }
  part<-cbind(part,gitx,gity)
  dimnames(part)<-list(NULL,dimnames(part)[[2]])
  
  #switch to time back when times negative, switch to hours
  if(any(part[,"time"]<0))part[,"time"]<--part[,"time"]/60 else part[,"time"]<-part[,"time"]/60 
  dimnames(part)[[2]][dimnames(part)[[2]]=="time"]<-"btime"
  part<-part[order(part[, "btime"], part[, "index"]),] #order by btime, then by index
  
  if (dmassTF) {
    #remove particles with too strong dmass violation
    ind<-unique(part[part[,"dmass"]>1E3|part[,"dmass"]<1/1E3,"index"])
    if (length(ind) >= length(unique(part[, "index"]))/2){
      message("Trajecvprm(): ", length(ind), ' of ', length(unique(part[, "index"])), ' particles have mass defect; returning NA')
      return(NA)
    }
    part<-part[!part[,"index"]%in%ind,]
    
    # get average dmass to "correct correction" (allow multiplication w/ dmass without changing total mass)
    # i.e. correction for average mass loss of particles, since they get attracted to areas of mass destruction
    mean.dmass <- tapply(part[, "dmass"], part[, "btime"], mean) # this gives for each btime a mean dmass
    # DMM
    # To account for situations where mean.dmass is zero (mass violation total), need to avoid division by zero to
    # avoid downstream problems.
    mean.dmass[which(mean.dmass == 0)] <- 0.00001
    
    # need to "merge" this with part; can't use array since not all paticles are left at very large btime
    nparleft <- rle(part[, "btime"])$length # number of times the same btime is repeated
    mean.dmass <- rep(mean.dmass, nparleft) # long vector of mean dmass
    # need to link this info to each particles dmass: normalize individual dmass by mean.dmass
    part[, "dmass"] <- part[, "dmass"]/mean.dmass             # Dan Matross gets problems with that
  }  
  
  nparstilt<-length(unique(part[,"index"]))
  
  #need to order (might not be...)
  part<-part[order(part[,"btime"],part[,"index"]),]
  
  #remove particles when they cross the longitude -145 West for the first time (that's where the climatology is valid)
  inbgarea<-floor(part[,"gitx"])<1
  part[inbgarea,c("gitx","gity")]<-NA
  dimnames(part)<-list(NULL,dimnames(part)[[2]])
  
  #remove points when they enter the background area for the first time
  sumx<-tapply(part[,"gitx"],part[,"index"],cumsum) #cumsum gives "NA" after first occurence of "NA"
  ordert<-order(part[,"index"],part[,"btime"]) #order first by index, then by time
  ordern<-order(part[ordert,"btime"],part[ordert,"index"]) #to recreate original order
  sumx<-unlist(sumx)[ordern]
  part<-part[!is.na(sumx),]
  dimnames(part)<-list(NULL,dimnames(part)[[2]])
  
  #only keep points, when information is changed:
  if(zlim[2]==0){ #SURFACE INFLUENCE ONLY; need to apply selection (only first, last, or influence)
    inngm<-floor(part[,"gitx"])<=numpix.x&floor(part[,"gitx"])>=1&floor(part[,"gity"])<=numpix.y&floor(part[,"gity"])>=1
    selinf<-part[,"foot"]>0&inngm
    part<-part[selinf>0,]
  }
  if(zlim[2]>0){#want "volume" influence rather than surface influence
    #get residence time (i.e. time used for each time step), convert from hours to minutes
    for(id in unique(part[,"index"]))part[part[,"index"]==id,"foot"]<-c(part[part[,"index"]==id,"btime"][1],diff(part[part[,"index"]==id,"btime"]))*60
    part<-part[part[,"agl"]>zlim[1]&part[,"agl"]<=zlim[2],]
  }
  
  if(is.null(dim(part)))return()
  
  #move x and y position of final position to initialization area (gitx=1, gity= 1 to numpix.y), at least make sure they are not outside domain
  part[part[,"gitx"]>numpix.x,"gitx"]<-numpix.x
  part[part[,"gity"]>numpix.y,"gity"]<-numpix.y
  part[part[,"gitx"]<1,"gitx"]<-1
  part[part[,"gity"]<1,"gity"]<-1
  
  part[,"btime"]<-round(part[,"btime"],2)
  
  #create object for output: 3d array (lat-lon-time)
  foot.arr<-array(0,dim=c(numpix.y,numpix.x,length(foottimes)-1))
  
  for(foottimespos in 1:(length(foottimes)-1)){ #loop over time intervals
    
    # adding ,drop=FALSE to prevent R from converting 1-D matrix to vector-DW,DM
    subpart<-part[(part[,"btime"]>(foottimes[foottimespos]))&(part[,"btime"]<=(foottimes[foottimespos+1])),,drop=FALSE]
    
    if(length(subpart)<=21)next
    
    subpart<-subpart[order(subpart[,"btime"],subpart[,"index"]),,drop=FALSE]
    
    if(!is.null(fluxweighting)){
      fluxweight<-fluxweighting; 
      if(fluxweighting=="CO2"|fluxweighting=="CO"){ #for CO and CO2: time of day and day of week weighting
        #time factor for CO emissions
        ltime<-weekdayhr(yr4,mon,day,hr,-subpart[,"btime"]*60,diffGMT=round(subpart[,"lon"]*24/360)) #last column is weekday, sunday=0, monday=1 etc.
        #hourly factors from NAPAP
        hfac<-c(0.272,0.231,0.214,0.226,0.322,0.705,1.22,1.39,1.35,1.40,1.49,1.54,1.58,1.59,1.65,1.74,1.69,1.43,1.05,0.789,0.688,0.592,0.483,0.370)
        #weekday-factors: sun,mon,tue,wed,thu,fri,sat,sun; napap-grid contains fluxes for summer-weekday, therefore weekday factor of 1
        dfac<-c(0.865,1,1,1,1,1,0.875)
        emfac<-hfac[ltime[,"hr"]+1]*dfac[ltime[,"weekd"]+1]
        if(fluxweighting=="CO2"){
          #use hourly factors from CO, but with reduced amplitude (0.4 reduction)
          hfac<-1+0.4*(hfac-1)
          #weekday-factors: use factors from CO, but w/ reduced amplitude, and with mean==1 (CO2 emissions are for average day, not for weekday)
          dfac<-dfac/mean(dfac)
          dfac<-1+0.4*(dfac-1)
          emfac<-hfac[ltime[,"hr"]+1]*dfac[ltime[,"weekd"]+1]
        }
        if(landcov=="GLCC"|landcov=="IGBP") nBaseVeg<-17
        if(landcov=="SYNMAP") nBaseVeg<-48
        if(landcov=="SYNMAP.VPRM8") nBaseVeg<-8
        if(landcov=="DVN") nBaseVeg<-12
        if(fluxweighting=="CO2")fluxweight<-nBaseVeg+1
        if(fluxweighting=="CO")fluxweight<-nBaseVeg+2
      } #
    } #fluxweighting
    
    
    #get different resolutions for surface grids depending on range in x and y and on particle number for each timestep
    #get selector for first and last row w/ a given btime
    selfirst<-c(T,diff(subpart[,"btime"])>0)
    selast<-c(diff(subpart[,"btime"])>0,T)
    
    max.x<-subpart[order(subpart[,"btime"],subpart[,"gitx"]),,drop=FALSE][selast>0,"gitx"]
    min.x<-subpart[order(subpart[,"btime"],subpart[,"gitx"]),,drop=FALSE][selfirst>0,"gitx"]
    max.y<-subpart[order(subpart[,"btime"],subpart[,"gity"]),,drop=FALSE][selast>0,"gity"]
    min.y<-subpart[order(subpart[,"btime"],subpart[,"gity"]),,drop=FALSE][selfirst>0,"gity"]
    btime<-subpart[order(subpart[,"btime"],subpart[,"gity"]),,drop=FALSE][selfirst>0,"btime"]
    names(max.x)<-NULL;names(min.x)<-NULL;names(max.y)<-NULL;names(min.y)<-NULL
    
    #now get information back in format for all timesteps and index-numbers
    minmax.yx<-cbind(btime,max.x,min.x,max.y,min.y)
    minmax.yx<-merge(subpart[,c("btime","index"),drop=FALSE],minmax.yx,by="btime")
    max.x<-minmax.yx[,"max.x"]
    min.x<-minmax.yx[,"min.x"]
    max.y<-minmax.yx[,"max.y"]
    min.y<-minmax.yx[,"min.y"]
    names(max.x)<-NULL;names(min.x)<-NULL;names(max.y)<-NULL;names(min.y)<-NULL
    
    #Call 'getgrid' to get correct emission grid--necessary b/c emission grid is too large, so divided into several diff objects
    #use getgridp.ssc function: don't allow the resolution to get finer at earlier backtime; use cummax(ran.x)
    
    gridresult<-getgridp(min.x,max.x,min.y,max.y,numpix.x,numpix.y,coarse=coarse)
    emissname<-paste(gridresult[,"xpart"],gridresult[,"ypart"],gridresult[,"gridname"],sep="")
    #Extract appropriate emissions within each emission grid--do one grid at a time b/c reduces # of times grid has to be accessed
    coarsex<-c(1,1,2,2,2,4,4,4,8,8,8,16,16,16,32,32)  #factors by which grids have been made coarser
    coarsey<-c(1,2,1,2,4,2,4,8,4,8,16,8,16,32,16,32)  #e.g., '4' means grid is 4 times coarser
    #loop over different surface grids
    emissgrid.all<-matrix(0,nrow=numpix.y,ncol=numpix.x) #initialize fine grid with "0" everywhere
    for (name in unique(emissname)){
      emissgrid<-matrix(0,nrow=numpix.y,ncol=numpix.x) #initialize fine grid with "0" everywhere
      if(!is.null(fluxweighting))emissgridwt<-getr(paste(veghead,fluxweight,".",name,sep=""),path=vegpath) #get CO emission grid
      sel<-emissname==name
      gridname<-unique(gridresult[sel,"gridname"])   #gridname can be 1~16, representing diff. resolutions of grid
      x<-subpart[sel,"gitx"]
      y<-subpart[sel,"gity"]
      #Convert mins & maxes from coordinate values to rows & columns
      shrink.x<-coarsex[gridname];shrink.y<-coarsey[gridname]
      #grids have NOT been divided
      x<-ceiling(x/shrink.x)
      y<-ceiling(y/shrink.y)
      #get index pairs to extract surface pixels
      yxv<-y*1000+x #separate in order of magnitude
      
      #get index pairs to extract surface pixels, for extracting emission grid (e.g. CO)
      yx<-cbind(y,x)
      
      ##########BUDGET##########
      if(zlim[2]==0){ #surface influence
        #Take emission values at particle positions; multiply by "foot", i.e. sensitivity of mixing ratio changes to fluxes, 
        #in ppm/(micro-mol/m^2/s) 
        influence<-subpart[sel,"foot"]
        #also multiplied by dmass (accumulated weight of particles due to mass violation, normalized by average dmass to conserve total mass over time)
        if(dmassTF)influence<-influence*subpart[sel,"dmass"]
        
        #USE CO EMISSION, TIME OF DAY AND DAY OF WEEK FACTOR AS WEIGHT
        if(!is.null(fluxweighting)){
          if(fluxweighting=="CO2"|fluxweighting=="CO"){
            influence<-influence*emfac[sel]*emissgridwt[yx]
          } else {
            influence<-influence*emissgridwt[yx]
          }
        }
        
      }else{ #different "budget" for volume influence (i.e. zlim[2]>0)
        influence<-subpart[sel,"foot"]*60 #res. time in seconds
        if(dmassTF)influence<-influence*subpart[sel,"dmass"]
      }
      
      for(yxg in unique(yxv)){ #for each gridcell...
        #now need to vary each x and y between coarse grid cells to map it on fine grid cell
        selunix<-yxv==yxg
        xfine<-yxg-floor(yxg/1000)*1000
        yfine<-floor(yxg/1000)
        xfine<-(xfine-1)*shrink.x+1
        yfine<-(yfine-1)*shrink.y+1
        xfine<-xfine:(xfine+shrink.x-1)
        yfine<-yfine:(yfine+shrink.y-1)
        #cut out areas which are not in fine grid
        xfine<-xfine[xfine<=numpix.x]
        yfine<-yfine[yfine<=numpix.y]
        #get number of pixels (every combination if fine pixel indices)
        npixfx<-length(xfine)
        npixfy<-length(yfine)
        
        emissgrid[yfine,xfine]<-sum(influence[selunix])/(npixfx*npixfy) #'dilute' influence over larger area of fine grid
        
        
      } # for each gridcell...
      
      emissgrid.all<-emissgrid.all+emissgrid/nparstilt #uses number of particles used to determine footprint
    }#for different emissname
    
    #JCL(5/23/2004)------- not normalize by 'foottimes',b/c want TIME-INTEGRATED footprint----------------#
    #foot.arr[,,foottimespos]<-emissgrid.all/(foottimes[foottimespos+1]-foottimes[foottimespos])
    #JCL(5/23/2004)------- not normalize by 'foottimes',b/c want TIME-INTEGRATED footprint----------------#
    foot.arr[,,foottimespos]<-emissgrid.all
  } #loop over time intervals
  
  #For horizontal grids (lower left corner of south-west gridcell: 11N,145W; resolution: 1/4 lon, 1/6 lat, 376 (x) times 324 (y))
  if(is.null(wrfinput))dimnames(foot.arr)<-list(seq(lat.ll,(lat.ll+(numpix.y-1)*lat.res),lat.res),seq(lon.ll,(lon.ll+(numpix.x-1)*lon.res),lon.res),foottimes[1:(length(foottimes)-1)])
  if(!is.null(wrfinput))dimnames(foot.arr)<-list(round(latx[round(numpix.x/2),-dim(latx)[2]],2),round(longx[-dim(longx)[1],round(numpix.y/2)],2),foottimes[1:(length(foottimes)-1)])
  return(foot.arr)
}
