trajwind<-function(yr,mon,day,hr,lat,lon,agl,
                   metlib=paste(unix("echo $HOME"),"/Metdat/",sep=""),
                   metfile=NULL,
                   rundir="~/STILT/Exe/",
                   nummodel=0,
                   metd="edas",
                   varsout=c("time","index","lon","lat","agl","grdht","zi","temp","temp0","rhf",
                     "swrad","totrain","convrain","shtf","lhtf","tcld","dens","pres"),
                   doublefiles=T,
                   outname="tdat",
                   outpath="./",# specify outpath to remove confusion in directory where info is returned to
                   ...){
#Uses differences in positions of mean trajectories to calculate the mean wind at a specified time & position
#'lat','lon',&'agl' can be VECTORS--then return a MATRIX
#Calls 'Trajec'
#2/21/2001 by JCL; Updated on 5/12/2004 by JCL
#
#  $Id: trajwind.r,v 1.7 2007/10/12 13:41:37 tnehrkor Exp $
#---------------------------------------------------------------------------------------------------

if(length(unique(c(length(lat),length(lon),length(agl))))>1)stop("'lat', 'lon', and 'agl' have to all be the same length!")

#run mean trajectory, then take difference in position from original position to derive winds
info<-Trajec(nturb=T,yr=yr,mon=mon,day=day,hr=hr,lat=lat,lon=lon,agl=agl,nhrs=1,numpar=1,outdt=1,delt=0,
             metlib=metlib,metfile=metfile,metd=metd,doublefiles=doublefiles,overwrite=T,
             rundir=rundir,nummodel=nummodel,varsout=varsout,outname=outname,outpath=outpath,...) #pass on the outpath to Trajec.r
tdat<-getr(info["outname"],path=outpath)# search .RDatatdat in outpath
if(length(tdat)==0){return(NA)}
if(length(tdat)==1){if(is.na(tdat))return(NA)}

sel<-abs(tdat[,"time"])==min(unique(abs(tdat[,"time"])))

#check whether get correct number of trajectories--has to be same number as length of vector of lat and lon
if(length(unique(tdat[sel,"index"]))!=length(lat)){
   print("Wrong # of simulated trajectories compared to specified starting locations; NA returned!")
   return(NA)
}  #if(length(unique(tdat[sel,"index"]))!=length(lat)){

nmins<-abs(tdat[sel,"time"])
delx<-1000*distance(x1=tdat[sel,"lon"],x2=lon,y1=lat,y2=lat)    #distance in [m]
ubar<-sign(tdat[sel,"lon"]-lon)*delx/(nmins*60)     #U-velocity [m/s]
dely<-1000*distance(x1=lon,x2=lon,y1=tdat[sel,"lat"],y2=lat)    #distance in [m]
vbar<-sign(tdat[sel,"lat"]-lat)*dely/(nmins*60)     #V-velocity [m/s]
delz<-abs(tdat[sel,"agl"]-agl)    #distance in [m]
wbar<-sign(tdat[sel,"agl"]-agl)*delz/(nmins*60)     #W-velocity [m/s]
if(length(lat)==1){
  result<-cbind(ubar,vbar,wbar,t(tdat[sel,varsout[-c(1:5)]]))
}else{
  result<-cbind(ubar,vbar,wbar,(tdat[sel,varsout[-c(1:5)]]))
}
return(result)

}

