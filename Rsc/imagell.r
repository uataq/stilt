imagell<-function(mat,lg=F,zlims=NULL,zoom=NULL,center=c(-72.172,42.536),rays=F,
                  lon.ll=-145,lat.ll=11,lon.res=1/4,lat.res=1/6,nlevel=64,col=topo.colors(nlevel),citiesTF=F,riversTF=F,
                  part=NULL,dev="X11",type="",res.png=72,map.title=NULL){
#function to plot lat/lon image with maps
#'mat' matrix different rows for different latitudes, diff. columns as diff. longitudes (e.g. surface fluxes or footprint)
#'lg' flag for log
#'zlims' vector c(zmin,zmax) for scaling
#'zoom' number corresponding to about half image size in degrees (e.g. zoom=1: 2*2 dg image)
#       if NULL: no zoom
#       if 0:    auto-zoom (only rectangular area with non-zero values, plus 10% boundary each side)
#'center' receptor location lat lon 
#'rays' flag; if TRUE, 32 rays and 30 circles corresponding to a polar projection are drawn
#'lon.ll' longitude of southwest corner of southwest corner gridcell
#       if NULL uses dimnames of mat as lat/lon
#'lat.ll' latitude of southwest corner of southwest corner gridcell
#'lon.res'			#resolution in degrees longitude
#'lat.res'			#resolution in degrees latitude
#'nlevel' Number of color levels used in legend strip
#'col' Color table to use for image ( see help file on image for details). Default is a pleasing range of 64 divisions on a topographic scale. 
#'citiesT^F' if T, plot locations of cities
#'part' a particle location object can be passed on, it will be plotted on top of footprints
#'dev' device to which graphic is send (currently X11 or file name with .png ending, creates png file)
#'type' determines label for color legend
#       "foot" is footprint, uses ppm/micro-moles/m2/s
#       "infl" is influence, uses 1/m^2 (influence density integrated over altitude)
#       any other: used as label for color legend
#'res.png'			#resolution for png file (in dots per inch)
#'map.title'			#character string, allows for title for plot
#03/02/2004 by CHG
#
#  $Id: imagell.r,v 1.9 2007-06-27 11:54:03 skoerner Exp $
#---------------------------------------------------------------------------------------------------

if(!("package:maps"%in%search())){library(maps);print("loaded library maps")}
if(riversTF&!("package:mapdata"%in%search())){library(mapdata);print("loaded library mapdata")}
if(!("package:fields"%in%search())){library(fields);print("loaded library fields")}
if(".png"==substring(dev,nchar(dev)-3,nchar(dev))){bitmap(dev,res=res.png);print(paste("output in ",dev,sep=""))}
  if(lg)mat<-log10(mat)
  if(is.null(lon.ll)){ #use grid specifications from dimnames
    xs<-as.numeric(dimnames(mat)[[2]]); xs<-xs+0.5*(xs[2]-xs[1]) #centers of cells
    ys<-as.numeric(dimnames(mat)[[1]]); ys<-ys+0.5*(ys[2]-ys[1]) #centers of cells
  } else { #use grid specifications from arguments
    xs<-seq(lon.ll,by=lon.res,length=dim(mat)[2])+0.5*lon.res #centers of cells
    ys<-seq(lat.ll,by=lat.res,length=dim(mat)[1])+0.5*lat.res #centers of cells
  }
  if(!is.null(zoom)){
    if(zoom!=0){ #zoom in over HF
      xsel<-xs>(center[1]-zoom)&xs<(center[1]+zoom)
      ysel<-ys>center[2]-0.7*zoom&ys<center[2]+0.7*zoom #aspect ratio 0.7 deg lat per deg lon
      xs<-xs[xsel]
      ys<-ys[ysel]
      mat<-mat[ysel,xsel]
    } else { #when zoom 0, auto zoom
      if(!lg){
        xsel<-apply(mat,2,sum)>0
        ysel<-apply(mat,1,sum)>0
      } else { #back to normal scale
        xsel<-apply(exp(mat),2,sum)>0
        ysel<-apply(exp(mat),1,sum)>0
      }
      xran<-max(which(xsel))-min(which(xsel))
      yran<-max(which(ysel))-min(which(ysel))
      frac<-0.2
      xsel1<-max(min(which(xsel))-floor(xran*frac),1)
      xsel2<-min(max(which(xsel))+floor(xran*frac),dim(mat)[2])
      ysel1<-max(min(which(ysel))-floor(yran*frac),1)
      ysel2<-min(max(which(ysel))+floor(yran*frac),dim(mat)[1])
      #reset if aspect not ok
      increase.x<-abs(ys[ysel2]-ys[ysel1])/(abs(xs[xsel2]-xs[xsel1])*0.7) #factor by which need to increase x scale
      if(increase.x<1){#need more latitude range
        frac<-frac+(1/increase.x-1)/2 #each side
        ysel1<-max(min(which(ysel))-floor(yran*frac),1)
        ysel2<-min(max(which(ysel))+floor(yran*frac),dim(mat)[1])
      } else { #need more longitude range
        frac<-frac+(increase.x-1)/2 #each side
        xsel1<-max(min(which(xsel))-floor(xran*frac),1)
        xsel2<-min(max(which(xsel))+floor(xran*frac),dim(mat)[2])
      }
      xsel<-xsel1:xsel2
      ysel<-ysel1:ysel2
      xs<-xs[xsel]
      ys<-ys[ysel]
      mat<-mat[ysel,xsel]
    }
  } #when zoom null, don't zoom.

  if(is.null(zlims)){
    zlims<-range(mat,na.rm=T)
    if(is.inf(zlims[1])){#get the first non-infinite value
      mato<-mat[order(mat)]
      zlims[1]<-mato[sum(is.inf(mato))+1]
    }
  }
print(zlims)
  xlims<-range(xs)
  ylims<-range(ys)

  old.par <- par(no.readonly = TRUE)
  par(mar=c(par()$mar[1:3],4))
  cex.set<-1.2
  par(cex=cex.set)
  par(mgp=c(2.0,0.3,0))
  par(tcl=-0.1)
  image.plot.fix(x=xs,y=ys,z=t(mat),zlim=zlims,nlevel=nlevel,col=col,xlab="longitude/[deg]",ylab="latitude/[deg]",legend.width=0.03,
		offset=0.05,legend.only=F,lg=lg)
  if(center[1]< -50 & center[2]>0){#somewhere in Amerika
    if(!is.null(zoom)){if(zoom>1|zoom==0)map("world",region=c("canada","mexico"),xlim=xlims,ylim=ylims,add=T)}
    else if(is.null(zoom)){map("world",region=c("canada","mexico"),xlim=xlims,ylim=ylims,add=T)}
    map("state",xlim=xlims,ylim=ylims,add=T)
  }
  if(center[1]> -50 & center[2]>0){
print(paste("zoom",zoom))
    if(!is.null(zoom)){if(zoom>1|zoom==0)map("world",xlim=xlims,ylim=ylims,add=T)}
    else if(is.null(zoom)){map("world",xlim=xlims,ylim=ylims,add=T)}
    if(xlims[1]>5&xlims[2]<10&ylims[1]>42&ylims[2]<52)map("france",xlim=xlims,ylim=ylims,add=T)
  }
  if(riversTF)map('rivers', add=TRUE, col="blue")
  if(citiesTF)map.cities(minpop=3E5)
  points(center[1],center[2],col=3,pch=3,cex=3)

  if(!is.null(part)){ #plot particle location on top
    part<-part[part[,"lat"]>ylims[1]&part[,"lat"]<ylims[2]&part[,"lon"]>xlims[1]&part[,"lon"]<xlims[2],]
    if("foot"%in%dimnames(part)[[2]]){
      points(part[part[,"foot"]==0,c("lon","lat")],cex=0.2,pch=16,col="gray")
      points(part[part[,"foot"]>0,c("lon","lat")],cex=0.2,pch=16,col="black")
    } else {
      points(part[,c("lon","lat")],cex=0.2,pch=16,col="black")
    }
  }
  

  if(type=="foot")zlab<-expression(paste("surface influence [ppm/(",mu,"mole/(",m^2*s,"))]"))
  if(type=="infl")zlab<-expression(paste("influence [1/",m^2,"]"))
  if(!(type%in%c("infl","foot")))zlab<-type
  mtext(zlab,side=4,line=3.4,cex=cex.set)
  box()
  par(old.par)

#also plot rays for sectors, and circles
  if(rays){
    #set up polar, dr=20 km @ 0, increasing ~ r, with 40km at 400km
    rs<-NULL;r<-0;dr<-20;r.dr<-5 #distance/delta(r)
    d.phi<-32 #sectors
    while(r<8000){r<-r+dr;dr<-max(r/r.dr,20);area<-pi*(r^2-(r-dr)^2)/400/r.dr #area in 20 km pixels, in single sector
      rs<-rbind(rs,c(r,dr,area))
    }
    dimnames(rs)<-list(NULL,c("r","dr","area"))
    
    llray1<-rp2ll(0,0)
    for(i in 0:31){
      p<-i*360/32
      llray2<-rp2ll(1000,p)
      lines(c(llray1[1],llray2[1]),c(llray1[2],llray2[2]))
    }
    r.mat<-rep(1,d.phi)%o%(0.5*(rs[,"r"]+c(0,rs[-length(rs[,"r"]),"r"]))) #center positions in r coordinates
    phi.mat<-(((1:d.phi)-0.5)*360/d.phi)%o%rep(1,length(rs[,"r"])) #center positions in phi coordinates
    gitll<-rp2ll(r.mat,phi.mat)
    gitll[,1:30][gitll[,1:30]<xlims[1]]<-NA;gitll[,1:30][gitll[,1:30]>xlims[2]]<-NA;
    gitll[,31:60][gitll[,31:60]<ylims[1]]<-NA;gitll[,31:60][gitll[,31:60]>ylims[2]]<-NA;
    points(gitll[,1:30],gitll[,31:60],pch=16,cex=0.5) #also cell centers
  } #if rays
if(!is.null(map.title))title(map.title)
if(".png"==substring(dev,nchar(dev)-3,nchar(dev))){dev.off()}
  return(zlims)
}


