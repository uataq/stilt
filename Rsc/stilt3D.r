stilt3D<-function(dat,col="blue",size=2){
#This functions generates an *interactive* 3D plot of STILT particles using the "rgl" package 
#   in R and also adds a map to the ground surface.  
#The "rgl" package is based upon the OpenGL 3D visualization device system and is really cool!
#   It allows you to rotate and zoom in/out of the resulting plot.
#Required R packages:  "rgl" and "maps" (can be downloaded from "http://www.r-project.org/")
#9/13/2008 by JCL (jcl@uwaterloo.ca)

library(rgl);library(maps)
library(mapdata)  #contains hi-res world map

xlims<-range(dat[,"lon"],na.rm=T)
ylims<-range(dat[,"lat"],na.rm=T)

open3d()
par3d(cex=0.6)  #don't want axis labels to be too large...
plot3d(x=dat[,"lon"],y=dat[,"lat"],z=dat[,"agl"],col=col,
       xlab="Lon",ylab="Lat",zlab="AGL [m]",box=F,axes=F,size=size)
axis3d(edge="x-");axis3d(edge="y-");axis3d(edge="z+")

#Add map to particle plot
#  extract map data:  world
mdat<-list(x=NULL,y=NULL)
canflag<-max(dat[,"lat"],na.rm=T)>42&(mean(dat[,"lon"],na.rm=T)>-150&mean(dat[,"lon"],na.rm=T)<(-48))
if(canflag){
  canada<-map("worldHires","Canada",plot=F)
  sel<-canada$x>=min(dat[,"lon"],na.rm=T)
  sel<-sel&canada$x<=max(dat[,"lon"],na.rm=T)
  sel<-sel&canada$y>=min(dat[,"lat"],na.rm=T)
  sel<-sel&canada$y<=max(dat[,"lat"],na.rm=T)
  mdat$x<-c(mdat$x,canada$x[sel],NA)
  mdat$y<-c(mdat$y,canada$y[sel],NA)
} #if(canflag){
usaflag<-ylims[1]<49&(mean(dat[,"lon"],na.rm=T)>-150&mean(dat[,"lon"],na.rm=T)<(-48))
if(usaflag){
  sdat<-map("state",plot=F)
  sel<-sdat$x>=xlims[1]
  sel<-sel&sdat$x<=xlims[2]
  sel<-sel&sdat$y>=ylims[1]
  sel<-sel&sdat$y<=ylims[2]
  mdat$x<-c(mdat$x,sdat$x[sel],NA)
  mdat$y<-c(mdat$y,sdat$y[sel],NA)
} #if(usaflag){
if(!canflag&!usaflag){
  mdat3<-map("worldHires",plot=F,xlim=xlims,ylim=ylims)
  mdat$x<-c(mdat$x,mdat3$x,NA)
  mdat$y<-c(mdat$y,mdat3$y,NA)
} #if(!canflag&!usaflag){
lines3d(x=mdat$x,y=mdat$y,z=0,col="black")

} #stilt3D<-function(dat){





