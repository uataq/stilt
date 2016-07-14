project.nice<-function(x.in,y.in,proj="+proj=lcc +lat_1=30 +lat_2=60 +lat_0=50 +lon_0=10",inv=FALSE,matTF=TRUE,matinTF=FALSE){
#project lat lon to x,y or vice versa
#proj: projection definition according to proj4 library
#inv: TRUE for x,y to lat,lon
#matTF: want to create matrix of all combinations x.in and y.in, and return matrix?
#matinTF: use matrix as input for x.in and y.in, and return matrix

#tmp<-map("france",plot=FALSE)
#plot(tmp$x,tmp$y)
#y.in<-tmp$y;x.in<-tmp$x
#xy<-cbind(x.in,y.in)
#x.in<-ll[,2];y.in<-ll[,1]

if(matTF){
x.mat<-rep(1,length(y.in))%o%x.in
y.in<-y.in%o%rep(1,length(x.in))
x.in<-x.mat
}

# identity if projection is geo
if (!is.na(pmatch("+proj=geo",proj))) {
  return(list(x.in,y.in)) 
}

#print(x.in)
#library(Rmap)
require(proj4)
xy<-cbind(as.vector(x.in),as.vector(y.in))
#proj<-paste(proj," +ellps=clrk66",sep="") #select projection (Clark 1866)
## proj<-paste(proj," +a=6371200 +es=0.0",sep="") #select geoid sphere (used by WRF!!!)
proj<-paste(proj," +a=6370000 +es=0.0",sep="") #select geoid sphere (used by WRF!!!)

lls<-project(xy=xy,proj=proj,inv=inv) #takes about 1 second...
x<-lls[,1]
y<-lls[,2]
if(matTF|matinTF){
  x<-matrix(x,nrow=dim(x.in)[1],byrow=FALSE)
  y<-matrix(y,nrow=dim(y.in)[1],byrow=FALSE)
}
return(list(x,y))
}
