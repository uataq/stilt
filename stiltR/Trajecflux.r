Trajecflux<-function(ident,pathname="",tracers=c("CO","CO2"),coarse=1,dmassTF=T,nhrs=NULL,
                     vegpath="/deas/group/stilt/Vegetation/",linveg=FALSE,detailsTF=FALSE){
#calculate tracer concentrations in trajectories
#maps trajectories onto surface fluxes

#'ident' is character value specifying the trajectory ensemble to look at 
#'pathname' is path where object with particle locations is saved
#'tracers' vector of names for which mixing ratios are wanted; any subset of c("co","co2","ch4","h2")
#'coarse' degrade resolution (for aggregation error): 0: only 20 km resolution; 
#        1-16: dynamic resolution, but limited highest resolution
#	coarse:	 (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
#	coarsex:c(1,1,2,2,2,4,4,4,8, 8, 8,16,16,16,32,32) #factors by which grids have been made coarser
#	coarsey:c(1,2,1,2,4,2,4,8,4, 8,16, 8,16,32,16,32)#e.g., '4' means grid is 4 times coarser
#'dmassTF' if TRUE, weighting by accumulated mass due to violation of mass conservation in met fields
#'vegpath' is path under wich vegetation and flux grids are stored
#'linveg' if TRUE, CO2 fluxes are represented as linear functions of temperature and radiation (otherwise GEE is non linear)
#'detailsTF' if TRUE, for each particle the flux contribution is saved in a big object (same name as ident, with "result" attached)
#30/04/01 by CHG
#
#  $Id: Trajecflux.r,v 1.13 2007-06-27 11:54:03 skoerner Exp $
#---------------------------------------------------------------------------------------------------

tracers<-tolower(tracers) #use lower case

#get names for output
tracersini<-c(paste(tracers,"ini",sep=""))
if("co"%in%tracers)tracersini<-c(tracersini,"coinio") #also want initial CO without OH...

names.output<-c("ident","latstart","lonstart","aglstart","btime","late","lone","agle",tracersini,"zi","grdht","nendinarea",paste("sd",tracersini,sep=""))
if("co2"%in%tracers)names.output<-c(names.output,"inflwater","co2ffm")
if("co"%in%tracers&!("co2"%in%tracers))names.output<-c(names.output,"landinfl")
if("co"%in%tracers)names.output<-c(names.output,"coffm")
if("co2"%in%tracers)names.output<-c(names.output,"inflfrst","geefrst","respfrst","inflshrb","geeshrb","respshrb",
	"inflcrop","geecrop","respcrop","inflwetl","geewetl","respwetl")
if("ch4"%in%tracers)names.output<-c(names.output,"ch4ffm")
if("h2"%in%tracers)names.output<-c(names.output,"h2m")

#Check if object exists
if(length(unix.shell(paste("if (-e ",pathname,".RData",ident,") echo 'found it!'; endif",sep=""),shell="/bin/csh"))>0){
  print(paste("Trajecflux(): starting with ident=",ident,sep="")) #found object
  part<-getr(ident,pathname) #get it
}else #check for .gz compression
if(length(unix.shell(paste("if (-e ",pathname,".RData",ident,".gz) echo 'found it!'; endif",sep=""),shell="/bin/csh"))>0){
  print(paste("Trajecflux(): starting with ident=",ident,sep="")) #found object
  part<-getr(ident,pathname,gz=TRUE) #get it
}else{#if not there, break out of function, return NA's
print(paste("object ",pathname,".RData",ident," NOT FOUND",sep=""))
  lastresult<-rep(NA,length(names.output)-1);lastresult<-c(ident,lastresult);names(lastresult)<-names.output
  return(lastresult)
} #if exists or not

#get time and position information from name (ident)
pos<-id2pos(ident)
time<-month.day.year(floor(pos[1]))
yr4<-time$year #4 digit year
yr<-yr4%%100 #2 digit year (or 1 digit...)
mon<-time$month
day<-time$day
hr<-round((pos[1]-floor(pos[1]))*24)

#check if required particle information is available
rqdnames<-c("time","lat","lon","agl","zi","index","temp0","foot")
if("co2"%in%tracers)rqdnames<-c(rqdnames,"swrad")
if(dmassTF)rqdnames<-c(rqdnames,"dmass")
for(nm in rqdnames){if(!(nm%in%dimnames(part)[[2]]))stop(paste("need column '",nm,"' for this run",sep=""))}

part[,"time"]<--part[,"time"]/60 #use time back from now on, and transform to hours
dimnames(part)[[2]][dimnames(part)[[2]]=="time"]<-"btime"
#print("1");print(dim(part))
#get grid indices
#For horizontal grids (lower left corner of south-west gridcell: 11N,145W; resolution: 1/4 lon, 1/6 lat, 376 (x) times 324 (y))
gitx<-floor(4*(part[,"lon"]+145)+1)
gity<-floor(6*(part[,"lat"]-11)+1)
numpix.x<-376;numpix.y<-324 #number of pixels in x & y directions in non-reduced grid
filehead<-"veg."
part<-cbind(part,gitx,gity)
dimnames(part)<-list(NULL,dimnames(part)[[2]])

if(dmassTF){
  #get average dmass to "correct correction" (allow multiplication w/ dmass without changing total mass)
  #i.e. correction for average mass loss of particles, since they get attracted to areas of mass destruction
  mean.dmass<-tapply(part[,"dmass"],part[,"btime"],mean) #this gives for each btime a mean dmass
  #need to "merge" this with part; can't use array since not all paticles are left at very large btime
  nparleft<-rle(part[,"btime"])$length #number of times the same btime is repeated
  mean.dmass<-rep(mean.dmass,nparleft) #long vector of mean dmass
  #need to link this info to each particles dmass: normalize individual dmass by mean.dmass
  #part[,"dmass"]<-part[,"dmass"]/mean.dmass #DMM this normalizes dmass to time-step ensemble average
}  
#To allow short runs, cut off after nhrs
if(!is.null(nhrs))part<-part[abs(part[,"btime"])<=abs(nhrs),]

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
#print("2");print(dim(part));print(max(part[,"btime"]))
#only keep points, when information is changed:
#1. position and times at boundary of NGM grid
#2. position and times when surface influences particles
#3. position and times at start of trajectory

#1. boundary
ordert<-order(part[,"index"],part[,"btime"]) #order first by index, then by time
ordern<-order(part[ordert,"btime"],part[ordert,"index"]) #to recreate original order
delbte<-c(diff(part[ordert,"btime"]),-1000)[ordern] #timestep will be negative at last obs. for each particle
selend<-delbte<0

#2. surface influence over NGM area
inngm<-floor(part[,"gitx"])<=numpix.x&floor(part[,"gitx"])>=1&floor(part[,"gity"])<=numpix.y&floor(part[,"gity"])>=1
selinf<-part[,"foot"]>0&inngm

#3. start of particle trajectory
delbte<-c(-1000,diff(part[ordert,"btime"]))[ordern] #timestep will be negative at first obs. for each particle
selfirst<-delbte<0

####################CO + OH losses########################################
#first get OH at particle position, calculate COrel. loss and CH4 source
#OH from SAS, parameterized as oh=oh0+oh1*p+oh2*p**2, with parameters ohi for each month and for 30 and 60 lat
pmb<-1013*exp(-part[,"agl"]/8000) #8 km scale height
ohm<-oh[oh[,"month"]==mon,]
oh60<-ohm[ohm[,"lat"]==60,"oh0"]+ohm[ohm[,"lat"]==60,"oh1"]*pmb+ohm[ohm[,"lat"]==60,"oh2"]*pmb**2 #OH at 60 north, same altitude
oh30<-ohm[ohm[,"lat"]==30,"oh0"]+ohm[ohm[,"lat"]==30,"oh1"]*pmb+ohm[ohm[,"lat"]==30,"oh2"]*pmb**2 #OH at 30 north, same altitude
oh60[oh60<0]<-0;oh30[oh30<0]<-0
#local OH, interpolated
ohl<-(oh60*(part[,"lat"]-30)+oh30*(60-part[,"lat"]))/(60-30)
ohl[ohl<0]<-0
tair.k<-part[,"temp0"]-part[,"agl"]*6.5/1000
#get CO loss
kohCO<-1.3E-13*(1+(0.6*pmb/1000))*300/tair.k #in 1/cm^3/sec
#integral (k*OH*dt)
#get dbtime in same format as kohCO and ohl
delbt<-c(0,diff(part[ordert,"btime"]))[ordern][!selfirst]*60*60 #timestep in seconds
kohdt<-kohCO[!selfirst]*ohl[!selfirst]*delbt
#integral over dt: sum over timesteps for each particle individually
Ikohdt<-tapply(kohdt,part[!selfirst,"index"],sum)
CO.frac<-exp(-Ikohdt) #preliminary result: factors for each particles initial CO
CO.fact<-part[,"btime"]*0 #initialize
CO.fact[selend]<-CO.frac
####################CO from CH4 + OH######################################
kohCH4<-2.3E-12*exp(-1765/tair.k);
kohCH4dt<-kohCH4[!selfirst]*ohl[!selfirst]*delbt
IkohCH4dt<-tapply(kohCH4dt,part[!selfirst,"index"],sum)
COfrCH4<-1.780*IkohCH4dt #CO from CH4 in ppm
#not all will make it, account for some losses (about half the losses for COini)
COfrCH4<-COfrCH4*exp(-0.5*Ikohdt)
CO.frCH4<-part[,"btime"]*0 #initialize
CO.frCH4[selend]<-COfrCH4
part<-cbind(part,CO.fact,CO.frCH4) #add both to particle location object

#keep only particles where they matter: create flag for first, last, or when surface influence
#apply selection here (only first, last, or when surface influence)
part<-part[(selinf+selend+selfirst)>0,]
#print("3");print(dim(part));print(sum(inngm));print(sum(part[,"foot"]>0))
#move x and y position of final position to initialization area (gitx=1, gity= 1 to numpix.y), at least make sure they are not outside NGM
part[part[,"gitx"]>numpix.x,"gitx"]<-numpix.x
part[part[,"gity"]>numpix.y,"gity"]<-numpix.y
part[part[,"gitx"]<1,"gitx"]<-1
part[part[,"gity"]<1,"gity"]<-1

#get different resolutions for surface grids depending on range in x and y and on particle number for each timestep
#get selector for first and last row w/ a given btime
selfirst<-c(T,diff(part[,"btime"])>0)
selast<-c(diff(part[,"btime"])>0,T)

max.x<-part[order(part[,"btime"],part[,"gitx"]),][selast>0,"gitx"]
min.x<-part[order(part[,"btime"],part[,"gitx"]),][selfirst>0,"gitx"]
max.y<-part[order(part[,"btime"],part[,"gity"]),][selast>0,"gity"]
min.y<-part[order(part[,"btime"],part[,"gity"]),][selfirst>0,"gity"]
btime<-part[order(part[,"btime"],part[,"gity"]),][selfirst>0,"btime"]
names(max.x)<-NULL;names(min.x)<-NULL;names(max.y)<-NULL;names(min.y)<-NULL

#now get information back in format for all timesteps and index-numbers
minmax.yx<-cbind(btime,max.x,min.x,max.y,min.y)
minmax.yx<-merge(part[,c("btime","index")],minmax.yx,by="btime")
max.x<-minmax.yx[,"max.x"]
min.x<-minmax.yx[,"min.x"]
max.y<-minmax.yx[,"max.y"]
min.y<-minmax.yx[,"min.y"]
names(max.x)<-NULL;names(min.x)<-NULL;names(max.y)<-NULL;names(min.y)<-NULL

#Call 'getgrid' to get correct emission grid--necessary b/c emission grid is too large, so divided into several diff objects
#use getgridp.ssc function: don't allow the resolution to get finer at earlier backtime; use cummax(ran.x)

gridresult<-getgridp(min.x,max.x,min.y,max.y,numpix.x,numpix.y,coarse)
emissname<-paste(gridresult[,"xpart"],gridresult[,"ypart"],gridresult[,"gridname"],sep="")
#Extract appropriate emissions within each emission grid--do one grid at a time b/c reduces # of times grid has to be accessed
coarsex<-c(1,1,2,2,2,4,4,4,8,8,8,16,16,16,32,32)  #factors by which grids have been made coarser
coarsey<-c(1,2,1,2,4,2,4,8,4,8,16,8,16,32,16,32)  #e.g., '4' means grid is 4 times coarser

#loop over different surface grids
vegs<-NULL #initialize vector containing all flux grid numbers
if("co2"%in%tracers)vegs<-c(vegs,1:18)
if("co"%in%tracers&!("co2"%in%tracers))vegs<-c(vegs,17) #water influence (1 - land influence...)
if("co"%in%tracers)vegs<-c(vegs,19)
if("ch4"%in%tracers)vegs<-c(vegs,20:31)
nveg<-length(vegs)
result<-cbind(part,matrix(ncol=nveg,nrow=length(part[,"btime"])))
dimnames(result)<-list(NULL,c(dimnames(part)[[2]],paste("v",vegs,sep="")))
for (vegnum in vegs){
for (name in unique(emissname)){
        emissgrid<-getr(paste(filehead,vegnum,".",name,sep=""),path=vegpath)  #select appropriate emission grid--note extension '.mat' for new format
        if(vegnum==17&!("co2"%in%tracers))emissgrid<-1-emissgrid #use land mask for CO emission w/o prior
        sel<-emissname==name
        #time<-part[sel,"btime"]
        xpart<-unique(gridresult[sel,"xpart"])   #xpart can be 0~3
        ypart<-unique(gridresult[sel,"ypart"])   #ypart can be 0~3
        gridname<-unique(gridresult[sel,"gridname"])   #gridname can be 1~16, representing diff. resolutions of grid
        x<-part[sel,"gitx"]
        y<-part[sel,"gity"]
        #Convert mins & maxes from coordinate values to rows & columns
        shrink.x<-coarsex[gridname];shrink.y<-coarsey[gridname]
        if (gridname>3|xpart==0|ypart==0){
                #grids have NOT been divided
                x<-ceiling(x/shrink.x)
                y<-ceiling(y/shrink.y)
        }else{
                #grids have been divided up
                x<-ceiling((x-ceiling((numpix.x/4)*(xpart-1)))/shrink.x) #checked, was wrong before 11/13/03
#                if(xpart>1) x<-floor((x-ceiling((numpix.x/4)*(xpart-1))+1)/shrink.x)
#                if(xpart==1) x<-ceiling(x/shrink.x)
                y<-ceiling((y-ceiling((numpix.y/4)*(ypart-1)))/shrink.y) #checked, was wrong before 11/13/03
#                if(ypart>1) y<-floor((y-ceiling((numpix.y/4)*(ypart-1))+1)/shrink.y)
#                if(ypart==1) y<-ceiling(y/shrink.y)
        }
	#get index pairs to extract surface pixels
	yx<-cbind(y,x)
#BUDGET
        #Take emission values at particle positions; multiply by "foot", i.e. sensitivity of mixing ratio changes to fluxes, 
        #in ppm/(micro-mol/m^2/s) 
        EMCO<-emissgrid[yx]*part[sel,"foot"]
	#also multiplied by dmass (accumulated weight of particles due to mass violation, normalized by average dmass to conserve total mass over time)
        if(dmassTF)EMCO<-EMCO*part[sel,"dmass"] 
        result[sel,paste("v",vegnum,sep="")]<-EMCO

}#for different emissname
}#for different surface grid
dimnames(result)<-list(NULL,dimnames(result)[[2]])
#add factors for emission (CO, CO2(ff), R, GEE)
###################################

#time factor for CO emissions
	ltime<-weekdayhr(yr4,mon,day,hr,-result[,"btime"]*60,diffGMT=round(result[,"lon"]*24/360)) #last column is weekday, sunday=0, monday=1 etc.
        #hourly factors from NAPAP
        hfac<-c(0.272,0.231,0.214,0.226,0.322,0.705,1.22,1.39,1.35,1.40,1.49,1.54,1.58,1.59,1.65,1.74,1.69,1.43,1.05,0.789,0.688,0.592,0.483,0.370)
        #weekday-factors: sun,mon,tue,wed,thu,fri,sat,sun; napap-grid contains fluxes for summer-weekday, therefore weekday factor of 1
	dfac<-c(0.865,1,1,1,1,1,0.875)

	if("co"%in%tracers){
	  emfacco<-hfac[ltime[,"hr"]+1]*dfac[ltime[,"weekd"]+1]
	  result[,"v19"]<-result[,"v19"]*emfacco
#	  if(!("co2"%in%tracers))result[,"v17"]<-result[,"v17"]*emfacco #land mask weighted w/ CO emissions
	}
	
#time factor for CO2 emissions
	if("co2"%in%tracers){
          #use hourly factors from CO, but with reduced amplitude (0.4 reduction)
	  hfac<-1+0.4*(hfac-1)
          #weekday-factors: use factors from CO, but w/ reduced amplitude, and with mean==1 (CO2 emissions are for average day, not for weekday)
	  dfac<-dfac/mean(dfac)
	  dfac<-1+0.4*(dfac-1)
	  emfacco2<-hfac[ltime[,"hr"]+1]*dfac[ltime[,"weekd"]+1]
	  result[,"v18"]<-result[,"v18"]*emfacco2
#factor for GEE and Respiration

	  #get parameters for coarser vegetation classes
	  #Water: 17 (water), 15 (Snow and Ice)
	  #vgroup1: Forrests: 1, 5, 4, 2, 3 (not existent)
	  #vgroup2: shrublands etc.:7,10,8,16,6,9
	  #vgroup3: Crops etc.: 14,12
	  #vgroup4: wetland: 11
	  #Rest:13(urban)
	  #regroup into fewer classes
	  idfrst<-c(1,2,3,4,5) #group 1
	  idshrb<-c(6,7,8,9,10,16) #group 2
	  idcrop<-c(12,14) #group 3
	  idwetl<-c(11) #group 4
	  #first guesses for parameters for coarser classes, wetland=shrubland in parameters
          if(!linveg){
  	    dresp.dT<-c(dlambda.veg[,"drdt"],dlambda.veg[2,"drdt"])
	    a3<-c(dlambda.veg[,"a3"],dlambda.veg[2,"a3"])
	    a4<-c(dlambda.veg[,"a4"],dlambda.veg[2,"a4"])
	    #get parameterized GEE and R
          }else{ #fluxes linear in temp and radiation
            dresp.dT<-c(dlambda.simp.veg[,"drdt"],dlambda.simp.veg[2,"drdt"])
            a3<-c(dlambda.simp.veg[,"a3"],dlambda.simp.veg[2,"a3"])
	  }#of if not linear fluxes as fct of temp and radiation
	  #calculate t and par dependent fluxes
	  v20<-apply(result[,paste("v",as.character(idfrst),sep='')],1,sum) #add all influences from this vegetation class
	  if(!linveg)v21<-a3[1]*result[,"swrad"]/(result[,"swrad"]+a4[1])*v20 #get GEE
	  if(linveg)v21<-a3[1]*result[,"swrad"]*v20 #get GEE
	  v22<-dresp.dT[1]*(result[,"temp0"]-273)*v20 #get respiration

	  v23<-apply(result[,paste("v",as.character(idshrb),sep='')],1,sum) #add all influences from this vegetation class
	  if(!linveg)v24<-a3[2]*result[,"swrad"]/(result[,"swrad"]+a4[2])*v23 #get GEE
	  if(linveg)v24<-a3[2]*result[,"swrad"]*v23 #get GEE
	  v25<-dresp.dT[2]*(result[,"temp0"]-273)*v23 #get respiration

	  v26<-apply(result[,paste("v",as.character(idcrop),sep='')],1,sum) #add all influences from this vegetation class
	  if(!linveg)v27<-a3[3]*result[,"swrad"]/(result[,"swrad"]+a4[3])*v26 #get GEE
	  if(linveg)v27<-a3[3]*result[,"swrad"]*v26 #get GEE
	  v28<-dresp.dT[3]*(result[,"temp0"]-273)*v26 #get respiration

	  v29<-result[,paste("v",as.character(idwetl),sep='')] #idwetl has only one element, so no apply necessary
	  if(!linveg)v30<-a3[4]*result[,"swrad"]/(result[,"swrad"]+a4[4])*v29 #get GEE
	  if(linveg)v30<-a3[4]*result[,"swrad"]*v29 #get GEE
	  v31<-dresp.dT[4]*(result[,"temp0"]-273)*v29 #get respiration

	  #keep only coarse vegetation classes, prescribed emissions (fossil fuel + ch4 prior), and water influence
	  v18<-result[,"v18"]
	  v17<-result[,"v17"]
	} #of if("co2"%in%tracers)

	if("co"%in%tracers)v19<-result[,"v19"]
	if("co"%in%tracers&!("co2"%in%tracers))v17<-result[,"v17"] #as "land mask" (non-water)
	if("ch4"%in%tracers){ #also put ch4 fluxes....
	  #CH4: need month selector
	  monthlong<-weekdayhr(yr4,mon,day,hr,-result[,"btime"]*60)[,"mon"]
	  v32<-result[,"v31"]*0 #initialize
	  for(monch4 in 1:12){ #loop over 12 months
	    v32[monthlong==monch4]<-result[monthlong==monch4,paste("v",19+monch4,sep="")]
	  }
	}

	if("h2"%in%tracers){ #H2: use soil moisture parameterization in presence of vegetation
	  #For volumetric soil moisture 0- 0.15:
	  #-Flux=80*soil moisture
	  #For soil moisture 0.15-0.50
	  #-Flux=-34.3*soil moisture + 17.14
	  #The units are nmol m-2 s-1.  You get max flux of 12 at 0.15 linearly
	  #decreasing to 0 at soil moisture values of 0 and 0.5. 
	  idveg<-c(idfrst,idshrb,idcrop,idwetl)
	  v33<-apply(result[,paste("v",as.character(idveg),sep='')],1,sum) #add all influences from this vegetation class
	  if(!("solw"%in%dimnames(result)[[2]]))stop("need soil moisture for H2")
	  vsm<-result[,"solw"]/1000
	  v33<-v33*( (vsm<0.15)*(80*vsm) + (vsm>=0.15&vsm<0.5)*(17.14-34.3*vsm))/1000 #from nmol to mumol
	}

	result<-result[,substring(dimnames(result)[[2]],1,1)!="v"]
	if("co2"%in%tracers)result<-cbind(result,v17,v18)
	if("co"%in%tracers&!("co2"%in%tracers))result<-cbind(result,v17) #as "land mask" (non-water)
	if("co"%in%tracers)result<-cbind(result,v19)
	if("co2"%in%tracers)result<-cbind(result,v20,v21,v22,v23,v24,v25,v26,v27,v28,v29,v30,v31)
	if("ch4"%in%tracers)result<-cbind(result,v32)
	if("h2"%in%tracers)result<-cbind(result,v33)

	#rename all flux inputs, start with v1
	dimnames(result)[[2]][substring(dimnames(result)[[2]],1,1)=="v"]<-paste("v",1:length(dimnames(result)[[2]][substring(dimnames(result)[[2]],1,1)=="v"]),sep="")
	#this means, v1 is water influence, v2 and v3 are CO2 and CO ff emission inputs, then dco2ini, then infl., gee and R due to different vegetations, 
	# then (if req.) CH4, then H2

	#get selector for last observation for each particle
	ordert<-order(result[,"index"],result[,"btime"]) #order first by index, then by time
	ordern<-order(result[ordert,"btime"],result[ordert,"index"]) #to recreate original order
	delbte<-c(diff(result[ordert,"btime"]),-1000)[ordern] #timestep will be negative at last obs. for each particle
	selend<-delbte<0

	#get selector for first observation for each particle
	delbte<-c(-1000,diff(part[ordert,"btime"]))[ordern] #timestep will be negative at first obs. for each particle
	selfirst<-delbte<0

	#get initial conditions for CO and CO2, read with read<-co2<-bg<-aug00.ssc and read<-co<-bg<-aug00.ssc
	#co2.ini and co..ini objects are 3d arrays, agl*lat*sasdate; sasdate is 0 at 1/1/1960 (i.e. elapsed days since 1/1/1960, same as julian(m,d,y))
	aglg<-round(result[selend,"agl"]/500) #get agl in 500 m intervals
	latg<-round((result[selend,"lat"]-10)/2.5) #get lat in 2.5 deg intervals, starting at 10 deg
	#get 1st boundary field to set things up
	if(!existsr(paste(tracers[1],".ini",sep=""),pathname))stop(paste("need ",tracers[1]," boundary condition",sep=""))
	ini<-getr(paste(tracers[1],".ini",sep=""),pathname)
	startday<-as.numeric(dimnames(ini)[[3]][1]);delday<-as.numeric(dimnames(ini)[[3]][2])-startday
	sasdate<-julian(mon,day,yr4)-result[selend,"btime"]/24-startday #days elapsed since beginning of .ini file
	pointer<-cbind(aglg+1,latg+1,round(sasdate/delday)+1) #array indices must start with 1
	#use constant ini field when no initial data available
	if(any(pointer[,3]>dim(ini)[3]))print(paste("Trajecflux(): extrapolating ",tracers[1],".ini, need later times",sep=""))
	if(any(pointer[,3]<1))print(paste("Trajecflux(): extrapolating ",tracers[1],".ini, need earlier times",sep=""))
	pointer[pointer[,3]>dim(ini)[3],3]<-dim(ini)[3]
	pointer[pointer[,3]<1,3]<-1
	#limit to lat range of initial field
	#limit pointer to valid values for latg and aglg
	pointer[pointer[,1]>dim(ini)[1],1]<-dim(ini)[1]
	pointer[pointer[,2]>dim(ini)[2],2]<-dim(ini)[2]
	pointer[pointer[,1]<1,1]<-1
	pointer[pointer[,2]<1,2]<-1

	for(tr in tracers[1:length(tracers)]){
	  if(tr>1) if(!existsr(paste(tr,".ini",sep=""),pathname))stop(paste("need ",tr," boundary condition",sep=""))
	  if(tr>1) ini<-getr(paste(tr,".ini",sep=""),pathname)
	  cini<-result[,1]*0 #initialize vector with the right length
	  cini[selend]<-ini[pointer]
	  result<-cbind(result,cini)
	}
		
	dimnames(result)[[2]][(dim(result)[2]-length(tracers)+1):dim(result)[2]]<-paste(tracers,"ini",sep="") #give advected boundary mixing ratio names

	if("co"%in%tracers){#CO: need also advected boundary value with no chemistry
	  coinio<-result[,"coini"] #no chemistry yet
	  result[selend,"coini"]<-result[selend,"CO.fact"]*coinio[selend]+result[selend,"CO.frCH4"]
	  result<-cbind(result,coinio)
	}
	dimnames(result)<-list(NULL,dimnames(result)[[2]])

	#add all influences
	v<-result[,substring(dimnames(result)[[2]],1,1)=="v"]
	
	sv<-v*0 #initialize sum of veg. influence
	dimnames(sv)<-list(NULL,paste("s",dimnames(sv)[[2]],sep=""))
	for(veg in 1:length(dimnames(sv)[[2]])){
		cumsamp<-tapply(v[,paste("v",veg,sep="")],result[,"index"],cumsum)
		cumsamp<-unlist(cumsamp)[ordern]
		sv[,paste("sv",veg,sep="")]<-cumsamp
	}#for veg in different surface flux variables
	#################################
	result<-cbind(result,sv)
	dimnames(result)<-list(NULL,dimnames(result)[[2]])

	#ONLY FOR WRITING TIMESERIES OF POLLUTANT INPUT
	timeseriesTF<-F
	if(timeseriesTF){
	  timesera<-tapply(result[,"v3"],result[,"btime"],sum)/100
	  timesera<-timesera[order(-as.numeric(names(timesera)))]
	  timesera<-cumsum(timesera)
	  timesera<-cbind(as.numeric(names(timesera)),timesera)
	  dimnames(timesera)<-list(NULL,c("btime","dCO"))
	  assignr(paste("timeser",ident,sep=""),timesera)
	}

	#drop all vegetation fluxes, keep only the accumulated fluxes
	result<-result[,substring(dimnames(result)[[2]],1,1)!="v"]
	dimnames(result)<-list(NULL,dimnames(result)[[2]])

	#keep only information which is required: last entry position and accumulated fluxes for coarse vegetation classes, and CO and CO2 fossil fuel emissions
	lastresult<-result[selend>0,]
	dimnames(lastresult)<-list(NULL,dimnames(lastresult)[[2]])

	#keep also certain parameters at beginning of particle trajectory
	firstresult<-result[selfirst>0,]

	#give each particle the same weight throughout the run,
	#discard a complete particle run when particle stops before max. btime and before it reaches boundary: done with object "part" higher up
	#get stdeviation for co2ini, coini
	sdini<-apply(lastresult[,tracersini],2,stdev)

	#get number of particles ending insight area before btime=btimemax; allow for ~2.5 degree band at border and 1 day time
	endinarea<-(lastresult[,"gitx"]>10&lastresult[,"gitx"]<(numpix.x-10))&(lastresult[,"gity"]>10&lastresult[,"gity"]<(numpix.y-10))&lastresult[,"btime"]<nhrs-0.05*nhrs
	endinarea<-sum(endinarea)
	#get averages
	lastresult<-apply(lastresult,2,mean)
	firstresult<-apply(firstresult,2,mean)
	#keep only ident, btime, lat, lon, co2ini, coini, mixed layer height and ground height, at start, #ending in area, sdev. co and co2 ini, influence from water, and veg. influences

	output<-c(pos, lastresult[c("btime","lat","lon","agl",tracersini,"zi","grdht")],endinarea,sdini,lastresult[substring(names(lastresult),1,2)=="sv"])

	names(output)<-names.output

#to get trajectory time series data for storage, average all particles, and keep 2-hourly data
#sum over all particles and divide by number of particles initially used (nparstilt-1), sum over timesteps within each 2-hour interval
#all done in get.mean.traj ("get_mean_traj.ssc")
#resulthr<-get.mean.traj(result,hr=2) #only for looking at data...
#assign("resulthr",resulthr,where=1,immediate=T)
if(detailsTF){
  selnam<-substring(dimnames(result)[[2]],1,2)=="sv"
  dimnames(result)[[2]][selnam]<-names.output[(length(names.output)-sum(selnam)+1):length(names.output)]
  assignr(paste(ident,"result",sep=""),result,pathname)
}
return(output)
}
