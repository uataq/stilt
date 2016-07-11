#This is a script that illustrates how to run STILT and how STILT works
#3/17/2004 by JCL
#
#  $Id: 0stilt_tutorial.r,v 1.6 2007-06-27 11:54:03 skoerner Exp $
#---------------------------------------------------------------------------------------------------

###### source necessary scripts and functions ######
sourcepath<-"./stiltR/";source("./stiltR/sourceall.r")   #'sourcepath' specifies the directory where all the scripts and functions are found

source("setStiltparam.r")
################################
#set parameters needed for STILT
################################

###### set directories ######
path<-"/Net/Groups/BSY/people/tkoch/STILT_icos_europe_100m_new/Output/RData/Hohenpeissenberg_501/"             #path where output gets saved, 


#####################################
####### PARTICLE LOCATIONS ##########

#----- starting times & locations -----#
yr<-7;mon<-6;day<-26;hr<-15
lat<-47.80;lon<-11.01;agl<-100

#2007x08x10x15x47.80Nx011.01Ex00100
###### Run STILT ######
info<-Trajec(yr=yr,mon=mon,day=day,hr=hr,lat=lat,lon=lon,agl=agl,nhrs=nhrs,
               doublefiles=T,metd=c("ECmetF"),delt=delt,veght=veght,nummodel=nummodel,
               metlib=metpath,conv=convect,overwrite=overwrite,outpath=path,varsout=varstrajec,rundir=rundir)

#get back the object that has been stored in .RData* file (filename is info["outname"]
dat<-getr(info["outname"],path=path)
#dat<-getr("2007x08x10x15x47.80Nx011.01Ex00100",path=path)
###### generate map of particles #####
library(maps);flttrack(lon=dat[,"lon"],lat=dat[,"lat"],col=2,cex=0.5,cities=F,
                    xylims=c(lon.ll,lon.ll+numpix.y*lon.res,lat.ll,lat.ll+numpix.y*lat.res));
map("world",add=T)

par(mfrow=c(2,2))
plot(dat[,c("lon","lat")], type="n")
for(tt in unique(dat[,"time"]))points(dat[dat[,"time"]==tt,c("lon","lat")],pch=19,cex=0.5)

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
sel<-abs(dat[,"time"])==60*2
plot(dat[sel,c("sigmaw","agl")],pch=16) #standard deviation of w' [m/s]
abline(h=dat[sel,"zi"],col=8,lwd=3)     #add horizontal line indicating the mixed-layer height
plot(dat[sel,c("TL","agl")],pch=16)     #Lagrangian timescale [s]
abline(h=dat[sel,"zi"],col=8,lwd=3)     #add horizontal line indicating the mixed-layer height


####### FOOTPRINTS AND INFLUENCE ##########
#Output will be a 3 D array of Influence (or Surface Influence/Footprints)
#     (i.e., sensitivity of atmospheric concentrations to surface fluxes) 
#these objects are given names that reflect starting location and time
#e.g. .RDatafoot2002x08x01x00x42.54Nx072.17Wx00030

#below lines are "done" in setStiltparam.r
#foottimes<-c(0,24,48)           #vector of times (backtimes) in hours between which footprint is computed
#zbot<-0                         #lower vertical bound for influence projection, in meters agl
#ztop<-0                         #upper vertical bound for influence projection, in meters agl
#                                #   if ztop set to zero, *surface* influence will be calculated
#numpix.x<-376                   #number of pixels in x directions in grid
#numpix.y<-324                   #number of pixels in y directions in grid
#lon.ll<--145                    #lower left corner of grid
#lat.ll<-11                      #lower left corner of grid
#lon.res<-1/4                    #resolution in degrees longitude
#lat.res<-1/6                    #resolution in degrees latitude

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

imagell(xfoot,lg=F,zlims=range(xfoot),zoom=NULL,center=c(lat,lon),rays=F,
                  lon.ll,lat.ll,lon.res,lat.res,nlevel=64,col=image.wry(64),citiesTF=F,riversTF=F,
                  part=NULL,dev="X11",type="",res.png=72,map.title=NULL,map.title.sub=NULL)

