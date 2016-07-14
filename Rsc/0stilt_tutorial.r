#This is a script that illustrates how to run STILT and how STILT works
#3/17/2004 by JCL
#
#  $Id: 0stilt_tutorial.r,v 1.6 2007-06-27 11:54:03 skoerner Exp $
#---------------------------------------------------------------------------------------------------

###### source necessary scripts and functions ######
#sourcepath<-"./";source("sourceall.r")   #'sourcepath' specifies the directory where all the scripts and functions are found


################################
#set parameters needed for STILT
################################

###### set directories ######
path<-"/home/johnlin/Rdat/"             #path where output gets saved, 
                                        #also input data (Receptor locations and times) 
                                        #and boundary mixing ratio objects are read from 'path'
metpath<-"/deas/group/stilt/Metdata/"        #where met data are stored in ARL format
vegpath<-"/deas/group/stilt/Vegetation/"     #path to get to surface fluxes and vegetation grids at various resolutions
rundir<-"/home/johnlin/STILT/Exe/"      #specifies main directory where different directories are found to run hymodelc


#####################################
####### PARTICLE LOCATIONS ##########
###### control parameters for what STILT should do: Particle locations, mixing ratios or footprints ######
varstrajec<-c("time","index","lat","lon","agl","grdht","foot","sampt","temp0","swrad","zi","dens","dmass","sigmaw","TL")
                        #specifies output variables from STILT can be any subset of 
                        #c("time","sigmaw","TL","lon","lat","agl","grdht","index","cldidx","temp0",
                        #  "foot","sampt","shtf","lhtf","tcld","dmass","dens","rhf","solw","lcld","zloc",
                        #  "swrad","wbar","zi","totrain","convrain")
overwrite=T             #T: rerun hymodelc, even if particle location object found; F: re-use previously calculated particle location object
#----- basic parameters -----#
nhrs<- -48              #for how many hours to run (time-backward runs negative).
nparstilt<-100          #how many particles per receptor? 100 of for receptor oriented modeling with dynamic grid resolution is ok
veght<-0.5              #surface layer (for fluxes); if less than 1 then as fraction of mlhgt; 0.5 is a good value.
convect<-F              #T for convection (RAMS winds: grell convection scheme, EDAS and FNL: simple redistribution within vertical range with CAPE>0)
stepsize<-0             #Enforces Courant criterion. When >0: maximum horizontal distances travelled by a particle in a single timestep, in km. 
                        #For dynamic resolution, choose value 0. First 12 hrs: 20 km; then 60 km
winderrTF<-F            #transport error due to wind error included as stochastic process?
nummodel<-0             #which copy of the model to run?  Will be in directory named paste(rundir,"Copy",nummodel,sep="")
delt<-30                #fixed timestep [min]; set =0 for dynamic timestep

#----- starting times & locations -----#
yr<-2;mon<-8;day<-12;hr<-18
lat<-42.53572;lon<--72.172;agl<-30

###### Run STILT ######
info<-Trajec(yr=yr,mon=mon,day=day,hr=hr,lat=lat,lon=lon,agl=agl,nhrs=nhrs,maxdist=stepsize,
               doublefiles=T,metd=c("edas","fnl"),delt=delt,veght=veght,nummodel=nummodel,
               metlib=metpath,conv=convect,overwrite=overwrite,outpath=path,varsout=varstrajec,rundir=rundir)

#get back the object that has been stored in .RData* file (filename is info["outname"]
dat<-getr(info["outname"],path=path)
###### generate map of particles #####
library(maps);flttrack(lon=dat[,"lon"],lat=dat[,"lat"],col=2,cex=0.5,cities=T);map("state",add=T)


###### Generate 3D plot of particles #####
library(scatterplot3d)
X11()
scatterplot3d(x=dat[,"lon"],y=dat[,"lat"],z=dat[,"agl"],pch=16,xlab="Longitude",ylab="Latitude",zlab="Altitude [m AGL]",
                            cex.symbols=0.5,angle=60,scale.y=1.2)

###### Generate time series of particle vertical distribution  #######
X11();xpar<-par()$mar;xpar[4]<-xpar[4]+3;par(mar=xpar);plot(dat[,c("time","agl")],pch=16,cex=0.4)
#  overlay average mixed-layer height at each timestep
xzi<-tapply(dat[,"zi"],dat[,"time"],mean,na.rm=T);tt<-as.numeric(names(xzi))
lines(tt,xzi,type="o",pch=16,col=2)
#  overlay average sampt at each timestep
#xsampt<-tapply(dat[,"sampt"],dat[,"time"],mean,na.rm=T);tt<-as.numeric(names(xsampt))
#par(new=T);plot(tt,xsampt,type="o",pch=16,col=3,axes=F,xlab="",ylab="");axis(4,col=3)
#  overlay average footprint at each timestep
xfoot<-tapply(dat[,"foot"],dat[,"time"],mean,na.rm=T);tt<-as.numeric(names(xfoot))
par(new=T);plot(tt,xfoot,type="o",pch=16,col=3,axes=F,xlab="",ylab="");axis(4,col=3)
#  overlay mass violation ("dmass") at each timestep
xdmass<-tapply(dat[,"dmass"],dat[,"time"],mean,na.rm=T);tt<-as.numeric(names(xdmass))
par(new=T);plot(tt,xdmass,type="o",pch=16,col=4,axes=F,xlab="",ylab="");axis(4,col=4,pos=par()$usr[2]+0.05*(par()$usr[2]-par()$usr[1]))
legend(x=par()$usr[1]+0.0*(par()$usr[2]-par()$usr[1]),y=par()$usr[3]+0.95*(par()$usr[4]-par()$usr[3]),
              c("agl","zi","foot","dmass"),pch=rep(16,4),col=1:4,ncol=2)


###### Generate vertical profiles of turbulence parameters #########
X11(width=10,height=8);par(mfrow=c(1,2))
sel<-abs(dat[,"time"])==60*36
plot(dat[sel,c("sigmaw","agl")],pch=16) #standard deviation of w' [m/s]
abline(h=dat[sel,"zi"],col=8,lwd=3)     #add horizontal line indicating the mixed-layer height
plot(dat[sel,c("TL","agl")],pch=16)     #Lagrangian timescale [s]
abline(h=dat[sel,"zi"],col=8,lwd=3)     #add horizontal line indicating the mixed-layer height


####### FOOTPRINTS AND INFLUENCE ##########
#Output will be a 3 D array of Influence (or Surface Influence/Footprints)
#     (i.e., sensitivity of atmospheric concentrations to surface fluxes) 
#these objects are given names that reflect starting location and time
#e.g. .RDatafoot2002x08x01x00x42.54Nx072.17Wx00030
foottimes<-c(0,24,48)           #vector of times (backtimes) in hours between which footprint is computed
zbot<-0                         #lower vertical bound for influence projection, in meters agl
ztop<-0                         #upper vertical bound for influence projection, in meters agl
                                #   if ztop set to zero, *surface* influence will be calculated
numpix.x<-376                   #number of pixels in x directions in grid
numpix.y<-324                   #number of pixels in y directions in grid
lon.ll<--145                    #lower left corner of grid
lat.ll<-11                      #lower left corner of grid
lon.res<-1/4                    #resolution in degrees longitude
lat.res<-1/6                    #resolution in degrees latitude

ident<-info["outname"]   #ident is identifier flag for an "2002x08x12x18x42.54Nx072.17Wx00030"
foot<-Trajecfoot(ident=ident,pathname=path,foottimes=foottimes,zlim=c(zbot,ztop),fluxweighting=NULL,coarse=1,vegpath=vegpath,
                 numpix.x=numpix.x,numpix.y=numpix.y,lon.ll=lon.ll,lat.ll=lat.ll,lon.res=lon.res,lat.res=lat.res)
xfoot<-apply(foot,c(1,2),sum)   #(a) generate *time-integrated* footprint
#xfoot<-foot[,,2]               #(b) only select one time-period
#  extract only the box enclosing non-zero footprint values
ilims<-range(row(xfoot)[xfoot>0]);jlims<-range(col(xfoot)[xfoot>0]);xfoot<-xfoot[ilims[1]:ilims[2],jlims[1]:jlims[2]]
xlon<-as.numeric(dimnames(xfoot)[[2]]);xlat<-as.numeric(dimnames(xfoot)[[1]])
xfoot.log<-log10(xfoot);xfoot.log[is.inf(xfoot.log)]<-NA
library(maps)
#colors.jcl<-c("white","purple","lightblue","darkblue","yellow","orange","red")  #color scheme (see output from 'colors()' for complete list of colors)
#  1) Use 'image' to plot
#library(fields)
#X11();image.plot(x=xlon,y=xlat,z=t(xfoot.log),col=col.grads(10),horizontal=T,legend.shrink=0.6,legend.width=0.03)
#map("state",add=T);map("world",c("Canada","Mexico"),add=T)
#  2) Use 'filled.contour' to plot
X11();filled.contour(x=xlon,y=xlat,z=t(xfoot.log),plot.axes={map("state",add=T);map("world",c("Canada","Mexico"),add=T)},
                           asp=1/cos(mean(dat[,"lat"],na.rm=T)*pi/180),color.palette=col.grads,nlevels=10)


