footplot<-function(foot,identname,lon.ll,lat.ll,lon.res,lat.res){
#plots footprints
#uses settings in setStiltparam.r
  zlims<-NULL
  if(all(is.na(foot))){print(paste("no surface influence",sep=""));return()} 
  if(is.null(foot)){print(paste("no surface influence",sep=""));return()} 
  if(all(range(as.vector(foot))==c(0,0))){print(paste("no surface influence",sep=""));return()} 
print(paste("creating footprint",sep=""))
  xfoot<-apply(foot,c(1,2),sum) #(a) generate *time-integrated* footprint
  #get date and time formatted; setting sec to 0.1 avoids skipping the time for midnight
  map.title<-paste("Footprint ",identname,sep="")
  #create directory if not available
  if(!file.exists(paste(path,"Footprints/",sep="")))dir.create(path=paste(path,"Footprints/",sep=""), showWarnings = TRUE, recursive = FALSE)
  try(imagell(xfoot,lg=T,zlims=zlims,zoom=NULL,center=id2pos(identname)[c("lon","lat")],rays=F,lon.ll,lat.ll,lon.res,lat.res,nlevel=64,citiesTF=F,part=NULL,
          dev=paste(path,"Footprints/",identname,".png",sep=""),type="foot",col=rev(heat.colors(64)),res.png=144,map.title=map.title))
  return()
}
