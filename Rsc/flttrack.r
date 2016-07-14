flttrack<-function(lon=NULL,lat=NULL,col=6,cex=0.8,cities=F,minpop=100000,pts=T,newmotifTF=T,database=NULL,xylims=NULL){
#Plots the flight tracks on top of U.S. state map
#'pts' specifies whether points are plotted or not
#If don't want to specify 'lon' & 'lat' separately, then could just give entire particle object as first argument--e.g., 'flttrack(bdat)'
#'minpop' specifies the minimum population in cities to plot
#'xylims': lower left and upper right corners (x1, x2, y1, y2)
#8/14/2000 by JCL
#
#  $Id: flttrack.r,v 1.13 2009-03-09 08:02:47 gerbig Exp $
#---------------------------------------------------------------------------------------------------

#------for debugging--------#
#lon<-c(-60,-90,-120);lat<-c(48,45,65)
#lon<-c(lon,-117.78662,-97.48695);lat<-c(lat,34.1248,24.71251)
#pts<-T;newmotifTF<-T;cities<-F
#cex<-0.8;col<-6
#------for debugging--------#

library(maps)

if(!is.null(lon)&is.null(lat)){
   if(sum(dimnames(lon)[[2]]=="lat")==1)lat<-lon[,"lat"]
   if(sum(dimnames(lon)[[2]]=="LAT")==1)lat<-lon[,"LAT"]
   if(sum(dimnames(lon)[[2]]=="lon")==1)lon<-lon[,"lon"]
   if(sum(dimnames(lon)[[2]]=="LON")==1)lon<-lon[,"LON"]
}

if(newmotifTF)motif()

xlims<-range(lon,na.rm=T);ylims<-range(lat,na.rm=T)
xlims[1]<-floor(xlims[1]-0.2*diff(xlims));xlims[2]<-ceiling(xlims[2]+0.2*diff(xlims))
ylims[1]<-floor(ylims[1]-0.2*diff(ylims));ylims[2]<-ceiling(ylims[2]+0.2*diff(ylims))
canflag<-F;mexicoflag<-F;franceflag<-F
map.database <- database
if (is.null(database)) {
  map.database <- 'state'
  if(ylims[2]>42&ylims[2]<47&xlims[2]>(-4)&xlims[2]<3)franceflag<-T
  if(ylims[2]>49)canflag<-T
  if(ylims[2]>43&xlims[2]>(-67))canflag<-T
  if(ylims[1]<28)mexicoflag<-T
  if(franceflag){canflag<-F;mexicoflag<-F}
} else {
  canflag <- database == 'canada'
  mexicoflag <- database == 'mexico'
  franceflag <- database == 'france'
}
if(!is.null(xylims)){
  canflag<-F;mexicoflag<-F;franceflag<-F
  xlims<-xylims[1:2];ylims<-xylims[3:4];
  if (is.null(database)) database <- 'world'
  map.database <- database
}

#implement aspect ratio to make pretty map
yctr<-mean(lat,na.rm=T);xctr<-mean(lon,na.rm=T)
dx<-abs(diff(xlims));dy<-abs(diff(ylims))
asp<-cos(pi*yctr/180)
if(is.null(xylims)){
  if(dy<asp*dx){
     dy<-asp*dx;ylims<-c(yctr-dy/2,yctr+dy/2)
  }else{
     dx<-dy/asp;xlims<-c(xctr-dx/2,xctr+dx/2)
  }
}
if(!canflag&!mexicoflag&!franceflag){
   map(map.database,xlim=xlims,ylim=ylims)
}
if(franceflag){
   xlims<-c(-4,2.5);ylims<-c(42.5,46.5)
   #these lines make it possible to plot France on same map
   plot(0,0,xlim=xlims,ylim=ylims,type="n",axes=F,xlab="",ylab="")
   #map("world",region=c("canada"),add=T,err=-1)
   map("france",add=T)
   cdat<-map("world",plot=F)
#   lines(x=cdat$x[(160:245)],y=cdat$y[(160:245)],lwd=1,col="black") #to cover up ugly borders
   lines(x=cdat$x[],y=cdat$y[],lwd=1,col="black") #
}  #if(franceflag){
if(canflag&!mexicoflag){
   #these lines make it possible to plot Canada on same map
   plot(0,0,xlim=xlims,ylim=ylims,type="n",axes=F,xlab="",ylab="")
   #map("world",region=c("canada"),add=T,err=-1)
   map("state",add=T)
   cdat<-map("world",region="canada",plot=F)
   #lines(x=cdat$x[120:280],y=cdat$y[120:280],lwd=3,col="black") #to cover up ugly borders
   lines(x=cdat$x[-1*(120:280)],y=cdat$y[-1*(120:280)],lwd=1,col="black") #to cover up ugly borders
   lines(x=cdat$x[(160:245)],y=cdat$y[(160:245)],lwd=1,col="black") #to cover up ugly borders
   lines(c(-79.5,-79.5,-74.5),c(50.8,47.0,45.0))  #add crude border betw Quebec & Ontario
}  #if(canflag&!mexicoflag){
if(!canflag&mexicoflag){
   #these lines make it possible to plot Mexico on same map
   plot(0,0,xlim=xlims,ylim=ylims,type="n",axes=F,xlab="",ylab="")
   map("state",add=T)
   #map("world",region=c("mexico"),add=T,err=-1)
   mdat<-map("world",region="mexico",plot=F)
   #lines(x=mdat$x[-1*(1:25)],y=mdat$y[-1*(1:25)],lwd=1,col="red") #to cover up ugly borders
   #lines(x=mdat$x[-1*(226:240)],y=mdat$y[-1*(226:240)],lwd=1,col="yellow") #to cover up ugly borders
   lines(x=mdat$x[-1*c(1:25,226:280)],y=mdat$y[-1*c(1:25,226:280)],lwd=1,col="black") #to cover up ugly borders
}  #if(!canflag&mexicoflag){
if(canflag&mexicoflag){
   #these lines make it possible to plot Canada & Mexico on same map
   plot(0,0,xlim=xlims,ylim=ylims,type="n",axes=F,xlab="",ylab="")
   map("state",add=T)
   #map("world",region=c("canada","mexico"),add=T,err=-1)
   lines(c(-79.5,-79.5,-74.5),c(50.8,47.0,45.0))  #add crude border betw Quebec & Ontario
   cdat<-map("world",region="canada",plot=F)
   mdat<-map("world",region="mexico",plot=F)
   #lines(x=cdat$x[120:280],y=cdat$y[120:280],lwd=3,col="black") #to cover up ugly borders
   lines(x=cdat$x[-1*(120:280)],y=cdat$y[-1*(120:280)],lwd=1,col="black") #to cover up ugly borders
   lines(x=cdat$x[(160:245)],y=cdat$y[(160:245)],lwd=1,col="black") #to cover up ugly borders
   lines(c(-79.5,-79.5,-74.5),c(50.8,47.0,45.0))  #add crude border betw Quebec & Ontario

   #lines(x=mdat$x[-1*(1:25)],y=mdat$y[-1*(1:25)],lwd=1,col="black") #to cover up ugly borders
   #lines(x=mdat$x[-1*(226:240)],y=mdat$y[-1*(226:240)],lwd=1,col="black") #to cover up ugly borders
   lines(x=mdat$x[-1*c(1:25,226:280)],y=mdat$y[-1*c(1:25,226:280)],lwd=1,col="black") #to cover up ugly borders
}  #if(canflag&mexicoflag){

if(pts)points(lon,lat,pch=16,col=col,cex=cex)
#if(cities)points(city.x,city.y,pch=16,col=3,cex=0.7,err=-1)
if(cities){
  #remove cities close to one another--otherwise letters cover one another and is ugly
  data(world.cities) #need to first load data
  tmp<-world.cities

  sel<-tmp$name%in%"Longueuil";sel<-sel|(tmp$name%in%"Laval")
  dum<-data.frame(tmp$name[!sel],tmp$lat[!sel],tmp$long[!sel],tmp$pop[!sel],
                tmp$country.etc[!sel],tmp$capital[!sel])
  names(dum)<-c("name","lat","long","pop","country.etc","capital")
  map.cities(x=dum,minpop=minpop)   #plots cities over minimum population specified by 'minpop'
 
} #if(cities){

}  #flttrack<-function(lon=NULL,lat=NULL,col=6,cex=0.8,cities=F,pts=T,newmotifTF=T){



