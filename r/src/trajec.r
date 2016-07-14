#' Run HYSPLIT dispersion model based on user parameters in SCRIPTNAME.r
#' 
#' @param yr year for testing
#' 
#' 
#' 
#' 
#' 
#' 
Trajec<-function(yr=02,mon=8,day=1,hr=6,lat=42.536,lon=-72.172,agl=30,nhrs=-48,
                 delt=0.0,numpar=1000,ndump=0,random=T,outdt=0.0,veght=0.5,metlib="/deas/group/stilt/Metdata/",
                 metd="edas",doublefiles=F,metfile=NULL,nturb=F,outfrac=0.9,conv=F,ziscale=NULL,mgmin=1000, kmix0 = 250, kmixd = 3,
                 siguverr=NULL,TLuverr=NULL,zcoruverr=NULL,horcoruverr=NULL,hymodelc.exe=NULL,
                 varsout=c("time","index","lat","lon","agl","grdht","foot","temp0","swrad","zi","dens","dmass"),
                 rundir=NULL,nummodel=NULL,outname=NULL,outpath="",overwrite=T,emisshrs=1/100,sourcepath="./",debugTF=F,max.counter=NULL,
                 sigzierr=NULL,TLzierr=NULL,horcorzierr=NULL,zsg.name=NULL,create.X0=FALSE){
  #Function to run HYSPLIT particle dispersion model and to check distribution of particles
  #Written by JCL, modified by DVM on 10/4/13 (added mgmin variable)
  #Modified by DVM 1/7/2016, added KMIX0 and KMIXD variables
  
  
  #INPUT:
  #'yr','mon','day','hr': starting time
  #'lat','lon',&'agl' can be a VECTOR of the same length--to have multiple starting locations; agl in meters above ground
  #'nhrs' is number of hours model would be run--NEGATIVE values mean model is run BACKWARDS
  #'delt' is timestep in minutes (integer...); if 0 then dynamic timestep (depending on 'trat')
  #'numpar' is number of particles emitted over the # of hrs specified in 'emisshrs'
  #'ndump' to dump out all particle/puff points to a file PARDUMP that can be read at start of new simulation to continue prev calc.
  #   valid NDUMP settings: 0 - no I/O, 1- read and write, 2 - read only, 3 - write only. Default value = 0
  #'random' is flag that means that random number generator would generate diff random number sequence each time model is run
  #'outdt' is the interval [min] that elapses before particle results are written out to PARTICLE.DAT
  #     if outdt=0.0, then data at EVERY timestep is written out; outdt should be a POSITIVE number
  #'veght' is height in meters above ground (model ground) below which time is counted as particle seeing the ground
  #     if <1 then interpreted as fraction of zi (mixed layer height as derived from met data)
  #'metd' is character vector with names (descriptors) of met files to be used; possible entries: "edas","fnl","RAMS" (not yet)
  #'doublefiles' should concatenated met files be used? Allows starting times between files
  #       concatenation with "cat file1 file2 > file12"
  #'metfile' specifies the meteorological input file;if not specified, then let 'getmetfile' automatically determine filename based on time and 'metd'
  #'nturb' is NotTURBulence flag that turns turbulence on (FALSE) or off (TRUE)
  #'outfrac' is fraction of particles which are allowed to leave the model area before hysplit stops
  #'conv' turns on convection (RAMS winds: grell convection scheme, EDAS and FNL: simple excessive redistribution within vertical range with CAPE>0)
  #'ziscale' is a vector with which to scale the modelled mixed-layer height
  #    each element specifies scaling factor for each model simulation hour (ziscale can be of length that is smaller than abs(nhrs))
  #'siguverr' & 'TLuverr' refer to the stddev of magnitude in horizontal wind errors [m/s] and their correlation timescale [min]
  #   'zcoruverr' refers to the vertical correlation lengthscale of horizontal winds [m]
  #   'horcoruverr' refers to the horizontal correlation lengthscale of horizontal winds [km]
  #'varsout' specifies output variables from STILT
  #      can be any subset of c("time","sigmaw","TL","lon","lat","agl","grdht","index","cldidx","temp","temp0","sampt","foot","shtf","lhtf","tcld","dmass","dens","rhf","sphu","solw","lcld","zloc","swrad","wbar","zi","totrain","convrain","zconv","pres")
  #'nummodel' specifies copy of directory where fortran executable is executed; needs to be different for different runs running parallel on same filesystem
  #'rundir' specifies main directory where different copy directories are found (see nummodel)
  #'outname' specifies name of the object for output; if not specified, uses default name
  #     based on time and position using pos2id() (e.g. "2002x08x16x06x42.54Nx072.17Wx00030")
  #'outpath' specifies the directory in which the object will be saved
  #'overwrite' if TRUE (default), overwrite existing object with same 'outname' or same default name
  #'emisshrs' specifies the hrs over which the particles will be emitted
  #'sourcepath' specifies the directory in which r scripts are located (needed for ECMWF level info)
  #'sigzierr' & 'TLzierr' refer to the stddev of magnitude in mixed layer height errors [%] and their correlation timescale [min]
  #   'horcorzierr' refers to the horizontal correlation lengthscale of mixed layer height errors [km]
  #'zsg.name' is the name of the file containing prescribed heights (e.g., "ZSG_LEVS.IN.AWRF" is read in by default for AWRF fields)
  #'mgmin' determines the size of the sub domain, set >1000 when working with WRF-ARW met fields as particles may terminate early without this variable
  
  
  #OUTPUT:
  #assigns the output of particle dispersion model in MATRIX format to an object called 'outname' (or default name indicating time & position);
  #object is saved in database at location depending on outpath
  #e.g. for outname="tmp" and outpath="/home/gerbig/modeloutput/" the database will be saved as "/home/gerbig/modeloutput/.RDatatmp"
  #and contains the object "tmp" that can be retrieved with getr("tmp",path="/home/gerbig/modeloutput/")
  #columns are specified using 'varsout' argument
  #returns list of:
  #  defaultname; all input data; metd with times when switched;
  #-status
  #  1: new object assigned, no problem;
  #  2: new object assigned, ended early
  #  3: object already exists, not overwriten
  #  4: no object assigned; failed
  #
  #
  #Note that 'metfile' can be a VECTOR of file names--necessary if the model run time encompasses info from TWO files
  #Returns the output of particle dispersion model in a matrix format
  #Calls 'getmetfile' to get name of necessary meteorological file
  #3/3/2004 by JCL and CHG
  #
  #  $Id: Trajec.r,v 1.42 2009/01/19 12:45:37 gerbig Exp $
  #---------------------------------------------------------------------------------------------------
  
  qcycle<-0;emissrate<-1;sampintv<-2;TLfrac=0.1; #define some former arguments
  #
  if(is.null(rundir))rundir<-"~/STILT/Exe/" #directory in which hymodelc is run (default)
  if(is.null(nummodel)){
    if(substring(unix("hostname"),1,3)=="nod"){ #on grid.deas, node decides which copy number
      nummodel<-as.numeric(substring(unix("hostname"),6,7))
    }else{
      nummodel<-0
    }}
  
  
  if((length(lat)!=length(lon))|(length(lon)!=length(agl))|(length(lat)!=length(agl))){
    stop("lat, lon, & agl have to be the same length!")
  } #if((length(lat)!=length(lon))|(length(lon)!=length(agl))|(length(lat)!=length(agl))){
  
  npos<-length(lat)
  status<-4 #assume worst case: nothin writen, nothing found
  yr4<-(yr+2000)*(yr<50)+(yr+1900)*(yr>=50) #4 digit year
  if(is.null(outname)){
    if(npos==1)outname<-pos2id(julian(mon,day,yr4)+hr/24,lat,lon,agl)
    else stop("'outname' needs to be specified when using multiple starting locations")
  }
  
  returninfo<-c(yr,mon,day,hr,lat,lon,agl,nhrs,delt,numpar,ndump,random,outdt,veght,metlib,metd,doublefiles,metfile,nturb,outfrac,conv,
                ziscale,siguverr,TLuverr,zcoruverr,horcoruverr,sigzierr,TLzierr,horcorzierr,varsout,nummodel,outname,outpath,overwrite,status)
  names.returninfo<-c("yr","mon","day","hr",paste("lat",1:npos,sep=""),paste("lon",1:npos,sep=""),paste("agl",1:npos,sep=""),"nhrs",
                      "delt","numpar","ndump","random","outdt","veght","metlib",paste("metd",1:length(metd),sep=""),"doublemetfiles")
  if(!is.null(metfile))names.returninfo<-c(names.returninfo,paste("metfile",1:length(metfile),sep=""))
  names.returninfo<-c(names.returninfo,"nturb","outfrac","conv")
  if(!is.null(ziscale))names.returninfo<-c(names.returninfo,paste("ziscale",1:length(ziscale),sep=""))
  if(!is.null(siguverr))names.returninfo<-c(names.returninfo,"siguverr")
  if(!is.null(TLuverr))names.returninfo<-c(names.returninfo,"TLuverr")
  if(!is.null(zcoruverr))names.returninfo<-c(names.returninfo,"zcoruverr")
  if(!is.null(horcoruverr))names.returninfo<-c(names.returninfo,"horcoruverr")
  if(!is.null(sigzierr))names.returninfo<-c(names.returninfo,"sigzierr")
  if(!is.null(TLzierr))names.returninfo<-c(names.returninfo,"TLzierr")
  if(!is.null(horcorzierr))names.returninfo<-c(names.returninfo,"horcorzierr")
  names.returninfo<-c(names.returninfo,paste("varsout",1:length(varsout),sep=""),"nummodel","outname","outpath","overwrite","status")
  
  names(returninfo)<-names.returninfo
  
  #JCL:(05/03/2004 this line causes model to crash, if it can't find the object
  #if(!overwrite&existsr(outname,outpath)){
  if(!overwrite){
    if(existsr(outname,outpath)){
      returninfo["status"]<-3
      cat("Trajec(): found object", outname, " in ", outpath, "; use this.\n")
      return(returninfo)
    } #if(existsr(outname,outpath)){
  } #check if there, if so, don't overwrite, just return with status 3
  
  rundir<-paste(rundir,"Copy",nummodel,"/",sep="")
  #cat("Trajec(): directory where STILT is run: ", rundir, "\n")
  #
  #
  #
  input1<-paste(rundir,"CONTROL",sep="")  #general input file for 'chghymodelc'
  input2<-paste(rundir,"SETUP.CFG",sep="")  #namelist file for 'chghymodelc'
  input3<-paste(rundir,"ZICONTROL",sep="") #file for prescribing mixed-layer heights
  input4<-paste(rundir,"WINDERR",sep="")  #wind error covariance
  input5<-paste(rundir,"ZSG_LEVS.IN",sep="")  #file for prescribing heights in met fields to achieve better match with internal levels
  input6<-paste(rundir,"ZIERR",sep="")  #mixed layer height error covariance
  #
  ##First delete any previous files--so not have same run results between runs if run doesn't succeed
  unix(paste("rm -f ",input1,sep=""))
  unix(paste("rm -f ",input2,sep=""))
  unix(paste("rm -f ",input3,sep=""))
  unix(paste("rm -f ",input4,sep=""))
  unix(paste("rm -f ",input5,sep=""))
  unix(paste("rm -f ",input6,sep=""))
  #unix(paste("rm -f ",rundir,"hymodelc.out",sep=""))
  #
  #Write the prescribed heights in met fields for ECMWF fields (hybrid coordinate)
  ecflag<-FALSE;if(!is.null(metfile)){if(length(grep("ec",tolower(metfile)))>0)ecflag<-TRUE}
  if(length(grep("ec",tolower(metd)))>0)ecflag<-TRUE
  if(ecflag){
    metfile1<-getmetfile(yr=yr,mon=mon,day=day,hr=hr,nhrs=nhrs,metd="ECmetF",doublefiles=doublefiles)[1] #get name(s) of met files required to drive model
    if (is.null(zsg.name)) {
      zname<-paste(substring(metfile1,1,nchar(metfile1)-nchar("arl")),"IN",sep="")
    } else {
      zname <- zsg.name
    }
    #print(paste(metlib,zname,sep=""))
    file.copy(from=paste(metlib,zname,sep=""), to=paste(rundir,"ZSG_LEVS.IN",sep=""), overwrite = TRUE) #use correct sigma levels, specific for ECMWF metdata file
  }  #if(ecflag){
  
  #Write the prescribed scaling factors for mixed-layer heights to 'ZICONTROL'
  if(!is.null(ziscale)){
    cat(paste(length(ziscale),"\n",sep=""),file=input3)
    for(j in 1:length(ziscale))cat(paste(ziscale[j],"\n",sep=""),file=input3,append=T)  #scaling factor for mixed-layer height
  }  #if(!is.null(ziscale)){
  
  #Write the stddev of magnitude in horizontal wind errors [m/s] and their correlation timescale [min] & length scales to 'WINDERR'
  if(!is.null(siguverr)&!is.null(TLuverr)&!is.null(zcoruverr)&!is.null(horcoruverr)){
    cat(paste(siguverr,"\n",sep=""),file=input4)
    cat(paste(TLuverr,"\n",sep=""),file=input4,append=T)
    cat(paste(zcoruverr,"\n",sep=""),file=input4,append=T)    #vertical correlation lengthscale [m]
    cat(paste(horcoruverr,"\n",sep=""),file=input4,append=T)  #horizontal correlation lengthscale [km]
  }  #if(!is.null(siguverr)&!is.null(TLuverr)&!is.null(zcoruverr)&!is.null(horcoruverr)){
  
  #Write the stddev of magnitude in mixed layer height errors [%] and their correlation timescale [min] & length scale to 'ZIERR'
  if(!is.null(sigzierr)&!is.null(TLzierr)&!is.null(horcorzierr)){
    cat(paste(sigzierr,"\n",sep=""),file=input6)
    cat(paste(TLzierr,"\n",sep=""),file=input6,append=T)
    cat(paste(horcorzierr,"\n",sep=""),file=input6,append=T)  #horizontal correlation lengthscale [km]
  }  #if(!is.null(sigzierr)&!is.null(TLzierr)&!is.null(horcorzierr)){
  
  #create batch file to run hymodelc in 'rundir'
  #create batch file to run hymodelc in 'rundir'
  batchname<-paste(rundir,"runhymodelc.bat",sep="")  #name for executable
  cat(paste("cd ",rundir,"\n",sep=""),file=batchname)
  if (is.null(hymodelc.exe)) hymodelc.exe <- "hymodelc"
  if (debugTF) {
    cat(paste(hymodelc.exe," >! hymodelc.out","\n",sep=""),file=batchname,append=T)
  } else {
    cat(paste(hymodelc.exe," >! hymodelc.out","\n",sep=""),file=batchname,append=T)
  }
  
  #translate variable names from R to fortran
  r.names<-c("time","sigmaw","TL",  "lon", "lat", "agl", "grdht","index","cldidx","temp","temp0","sampt","foot","shtf","lhtf","tcld","dmass","dens","rhf", "sphu","solw","lcld","zloc","swrad","wbar","zi",  "totrain","convrain","zconv","pres")
  f.names<-c("time","sigw",  "tlgr","long","lati","zagl","zsfc", "indx", "icdx",  "temz","temp", "samt", "foot","shtf","whtf","tcld","dmas", "dens","rhfr","sphu","solw","lcld","zloc","dswf", "wout","mlht","rain",   "crai","zcfx","pres")
  varsout.f<-f.names[match(varsout,r.names)]
  if(sum(is.na(match(varsout,r.names)))>0)stop(paste("wrong names:",varsout[is.na(match(varsout,r.names))]))
  n.col<-length(varsout)
  
  outdat<-NULL
  #timesofar<-0 #time calculated so far, as absolute value in hours
  #
  if(is.null(metfile)){
    for (i in 1:length(metd)){
      metf<-getmetfile(yr=yr,mon=mon,day=day,hr=hr,nhrs=nhrs,metd=metd[i],doublefiles=doublefiles) #get name(s) of met files required to drive model
      metfile<-c(metfile,metf) 
    }
  }
  cat("Trajec(): metfile that will be used: ", metfile, "\n")
  
  #check if metfiles available
  for(mm in metfile){  #loop over the number of meteorological files
    #  if(! my.file.exists(paste(metlib,mm,sep=""))){
    if(! is.element(mm,list.files(metlib))){
      stop(paste("Trajec(): Metfile ",metlib,mm," not found",sep=""))
    }
  }
  
  if(length(grep("awrf",tolower(metd)))>0) {
    #Write the prescribed heights in met fields for AWRF
    if (is.null(zsg.name)) {
      zname<-"ZSG_LEVS.IN.AWRF"
    } else {
      zname <- zsg.name
    }
    file.copy(from=paste(sourcepath,zname,sep=""), to=paste(rundir,"ZSG_LEVS.IN",sep=""), overwrite = TRUE)
  }
  
  trat<-0.75 #ensure Courant condition
  #
  #Generate 'CONTROL' file
  cat(paste(yr,' ',mon,' ',day,' ',hr,'\n',sep=""),file=input1)  #starting time
  if((length(lat)==length(lon))&(length(lon)==length(agl))){
    cat(paste(length(agl),'\n',sep=""),file=input1,append=T)  #print out number of starting locations
    for(i in 1:length(agl)){
      cat(paste(lat[i],' ',lon[i],' ',agl[i],'\n',sep=""),file=input1,append=T)  #print out each starting location
    }
  }else{
    return("need to have same length for lat, lon, & agl!")
  }
  cat(paste(round(nhrs),'\n',0,'\n','25000.0','\n',sep=""),file=input1,append=T)
  #numb of hrs model will be run, vertical motion calc method (0 = vertical velocity from data), top of model domain (m AGL)
  
  cat(paste(length(metfile),'\n',sep=""),file=input1,append=T)  #number of met data files
  for (i in 1:length(metfile)){ 
    cat(paste(metlib,'\n',metfile[i],'\n',sep=""),file=input1,append=T)  #met directory, met filename
  } 
  
  cat(paste(1,'\n','test','\n',emissrate,'\n',emisshrs,'\n',sep=""),file=input1,append=T)  #numb of pollutants, pollutant name, mass units emitted per hour, hours of emission
  cat(paste('00 00 00 00 00','\n',sep=""),file=input1,append=T)  #starting time of emissions--all 0's mean simulation starting time
  #   all default values for concentration grid definition
  cat(paste(1,'\n','0.0 0.0','\n','0.5 0.5','\n','30.0 30.0','\n','./','\n','cdump','\n','1','\n',sep=""),file=input1,append=T)
  cat(paste(100,'\n','00 00 00 00 00','\n','00 00 00 00 00','\n',sep=""),file=input1,append=T)
  cat(paste('00 ',sampintv,' 00','\n',sep=""),file=input1,append=T)   #'sampintv' is time interval(hrs) between which concentration grid  output would be written to file
  #   all default values for deposition definitions
  cat(paste(1,'\n','0.0 0.0 0.0','\n','0.0 0.0 0.0 0.0 0.0','\n','0.0 0.0 0.0','\n','0.0','\n','0.0','\n',sep=""),file=input1,append=T)
  #
  #
  #Generate 'Setup.cfg' file
  cat(paste(" &SETUP","\n",sep=""),file=input2)
  cat(paste(" TRATIO=",trat,",\n",sep=""),file=input2,append=T)  		#fraction of gridcell travelled by particles during a single timestep
  cat(paste(" INITD=0",",\n",sep=""),file=input2,append=T)  		#3-D particle simulation
  cat(paste(" KHMAX=9999",",\n",sep=""),file=input2,append=T)  		#max age a particle is allowed to attain
  cat(paste(" NUMPAR=",numpar,",\n",sep=""),file=input2,append=T)  	#total number of particles released at one time
  cat(paste(" QCYCLE=",qcycle,",\n",sep=""),file=input2,append=T)  	#number of hours between emission cycles
  cat(paste(" KRND=6",",\n",sep=""),file=input2,append=T)  		#at this interval in hrs, enhanced puff merging occurs
  cat(paste(" FRMR=0.0",",\n",sep=""),file=input2,append=T)  		#fraction of mass permitted to be removed at KRND intervals.
  cat(paste(" DELT=",delt,",\n",sep=""),file=input2,append=T)  		#nonzero value sets integration timestep to fixed step
  cat(paste(" ISOT=0",",\n",sep=""),file=input2,append=T)  		#flag used to set isotropic turbulence option
  cat(paste(" MGMIN=",mgmin,",\n",sep=""),file=input2,append=T)           #determines the size of the sub domain, set >1000 when working with high-res WRF-ARW met fields
  cat(paste(" KMIX0=",kmix0,",\n",sep=""),file=input2,append=T)           #minimum mixing depth (abs(kmix0) is used as the minimum mixing depth), 
  #negative values are used to force mixing heights coincident with model levels
  cat(paste(" KMIXD=",kmixd,",\n",sep=""),file=input2,append=T)           #PBL height computation: compute from bulk Ri profile (but see ziscale)
  cat(paste(" OUTFRAC=",outfrac,",\n",sep=""),file=input2,append=T)       #fraction of particles which are allowed to leave the model area before hysplit stops
  cat(paste(" NDUMP=",ndump,",\n",sep=""),file=input2,append=T)  		#can be set to dump out all the particle/puff points at the end of a simulation to a file called
  #PARDUMP. This file can be read from root directory at start of new simulation to continue previous calculation.
  #Valid NDUMP settings: 0 - no I/O, 1- read and write, 2 - read only, 3 - write only. Default value = 0
  i<-0;if(random)i<-1
  cat(paste(" RANDOM=",i,",\n",sep=""),file=input2,append=T)  		#flag used to set isotropic turbulence option
  cat(paste(" OUTDT=",outdt,",\n",sep=""),file=input2,append=T) 		#interval [min] that will determine how often particle data are written out to PARTICLE.DAT
  cat(paste(" VEGHT=",veght,",\n",sep=""),file=input2,append=T) 		#height below which time is counted as particle seeing the ground
  i<-0;if(nturb)i<-1
  cat(paste(" NTURB=",i,",\n",sep=""),file=input2,append=T)  		#flag used to switch off turbulence
  i<-0;if(conv)i<-1
  cat(paste(" ICONVECT=",i,",\n",sep=""),file=input2,append=T)  		#flag used to switch on convection
  if(is.null(ziscale)){  							#'ziscale' is vector of scaling factors used to prescribe mixed-layer height during model run
    cat(paste(" ZICONTROLTF=0,\n",sep=""),file=input2,append=T)
  }else{
    cat(paste(" ZICONTROLTF=1,\n",sep=""),file=input2,append=T)
  }  #if(is.null(ziscale)){
  if((is.null(siguverr)|is.null(TLuverr))&(is.null(sigzierr)|is.null(TLzierr))){
    cat(paste(" WINDERRTF=0,\n",sep=""),file=input2,append=T)
  }else if(!(is.null(siguverr)|is.null(TLuverr))&(is.null(sigzierr)|is.null(TLzierr))){
    cat(paste(" WINDERRTF=1,\n",sep=""),file=input2,append=T)
  }else if((is.null(siguverr)|is.null(TLuverr))&!(is.null(sigzierr)|is.null(TLzierr))){
    cat(paste(" WINDERRTF=2,\n",sep=""),file=input2,append=T)
  }else if(!(is.null(siguverr)|is.null(TLuverr))&!(is.null(sigzierr)|is.null(TLzierr))){
    cat(paste(" WINDERRTF=3,\n",sep=""),file=input2,append=T)
  }  #if((is.null(siguverr)|is.null(TLuverr))&(is.null(sigzierr)|is.null(TLzierr))){
  
  varsouttxt<-NULL; for(i in 1:length(varsout.f)){varsouttxt<-paste(varsouttxt,"'",varsout.f[i],"',",sep="")}
  cat(paste(" IVMAX=",length(varsout.f),",\n",sep=""),file=input2,append=T) #number of output variables
  cat(paste(" VARSIWANT=",varsouttxt,"\n",sep=""),file=input2,append=T) #4-letter code for output variables
  
  cat(paste(" /","\n",sep=""),file=input2,append=T)
  
  #remove old PARTICLE.DAT
  unix(paste("rm -f ",rundir,"PARTICLE.DAT",sep=""))
  unix(paste("rm -f ",rundir,"core",sep=""))
  #debug:
  if (debugTF) {
    for (tmp.file in paste(rundir,c("CONTROL","SETUP.CFG"),sep="")) {
      cat(tmp.file,":\n")
      unix(paste("cat",tmp.file,sep=" "),intern=F)
    }
  }
  #Call the executable 'chghymodelc'
  unix.shell(paste("sync",sep=""),shell="/bin/csh")   #sync before to make sure files are writen before read
  unix.shell(paste("source ",rundir,"runhymodelc.bat",sep=""),shell="/bin/csh")   #output in 'PARTICLE.DAT' in same directory
  unix.shell(paste("sync",sep=""),shell="/bin/csh")   #sync after to make sure files are writen before read
  
  #doesn't crash when PARTICLE.DAT isn't there
  if(! file.exists(paste(rundir, "PARTICLE.DAT",sep='')) ){
    dat<-matrix(nrow=2,ncol=n.col);dimnames(dat)<-list(NULL,varsout);print("Trajec(): PARTICLE.DAT not found");
  }else #doesn't read the data when core was dumped; necessary since PARDUMP file will not be updated
    if(length(unix.shell(paste("if (-e ",rundir,"core) echo 'core dumped'; endif",sep=""),shell="/bin/csh"))>0){
      dat<-matrix(nrow=2,ncol=n.col);dimnames(dat)<-list(NULL,varsout);print("Trajec(): core was dumped, don't use");
    }else #doesn't crash when PARTICLE.DAT is empty or has only one line (checking if larger than 500 b)
      if(as.numeric(unix(paste("cat ", rundir,"PARTICLE.DAT | wc -l",sep=""))) < 2){
        dat<-matrix(nrow=2,ncol=n.col);dimnames(dat)<-list(NULL,varsout);print("Trajec(): PARTICLE.DAT too short");
      }else{
        #DMM Modification to cover downstream to solve ******* in EDAS40 output of PARTICLE.DAT; This is a cluge!
        if(length(grep("edas40",tolower(metd)))>0) {
          datb<-scan(paste(rundir,"PARTICLE.DAT",sep=""),what=character(),skip=1,quiet=T) #now read it
          datb[which(substring(datb,1,2)=="**")]<-NA
          datbb<-as.numeric(datb)
          dat<-matrix(datbb,byrow=T,ncol=n.col)
        } else {
          dat<-matrix(scan(paste(rundir,"PARTICLE.DAT",sep=""),skip=1,quiet=T),byrow=T,ncol=n.col) #now read it
        }  #if(metdat=="edas40") {
        outdat<-rbind(outdat,dat)
        dimnames(outdat)<-list(NULL,varsout)
        
      } #of if can't read Particle.dat or else...
  
  metoutname<-NULL;for(i in 1:length(metfile))metoutname<-paste(metoutname,metfile[i],sep="x")
  #  returninfo<-c(returninfo,timesofar)
  names(returninfo)[length(returninfo)]<-metoutname
  returninfo["status"]<-1 #perfect, all times done
  if (debugTF) {
    for (save.file in c('MESSAGE','JCLmessage','hymodelc.out','PARTICLE.DAT'))
      file.copy(from=paste(rundir,save.file,sep=""), to=paste(rundir,save.file,".",loop.counter,".",outname,sep=""),
                overwrite = TRUE)
  }
  
  assignr(outname,outdat,outpath,printTF=T)
  unix(paste("rm -f ",input5,sep=""))
  if (create.X0) {
    min.time <- min(abs(outdat[,'time']))
    first.rows <- (1:dim(outdat)[1])[abs(outdat[,'time']) == min.time]
    outdat0 <- outdat[first.rows,,drop=FALSE]
    outname0 <- paste(outname,'X0',sep='')
    assignr(outname0,outdat0,outpath,printTF=T)
  }
  return(returninfo)
}

