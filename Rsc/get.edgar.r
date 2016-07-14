get.edgar<-function(ll.x=-12,ll.y=35,numpix.x=376,numpix.y=324,xres=1/4*1/2,yres=1/6*1/2,
                    epath="/Net/Groups/BSY/BSY_3/cgerbig/RData/ROAM/Fluxes/Emission/Global_Edgar32FT2000/",
                    tempvarTF=FALSE,year=2000,mon=NULL,day=NULL,hour=NULL,arrayTF=FALSE,plotTF=FALSE,spec=c("CO2","CO")){
#converts 1x1 deg edgar data to mivro-moles/m2/s for the STILT grid
#call with e.g.:
#edgar.eu<-get.edgar(ll.x=-12,ll.y=35,numpix.x=376,numpix.y=324,xres=1/4*1/2,yres=1/6*1/2,
#                    epath="/group/wofsy/stilt/Global_Edgar32FT2000/",
#                    tempvarTF=TRUE,year=2000,mon=1,day=2,hour=14,arrayTF=FALSE,plotTF=TRUE,spec="CO2")

#
#'numpix.x'			#number of pixels in x directions in output grid (STILT grid)
#'numpix.y' 			#number of pixels in y directions in output grid
#'ll.x'			#lower left corner of output grid (longitude of southwest corner of southwest corner gridcell)
#'ll.y'			#lower left corner of output grid (latitude of southwest corner of southwest corner gridcell)
#'xres'			#resolution in degrees longitude of output grid
#'yres'			#resolution in degrees latitude of output grid
#'epath'		#path where emissions files are located
#'tempvarTF'		#include temporal variation?
                        #uses: factors for Europe, from EDGAR, 
                        #takes into account southern hemisphere offset (6 months) for power, industry, and residential
                        #takes into account no monthly variation near equator (+/- 20 lat) for power, industry, and residential
                        #does not take into account slight changes of monthly variations in "processes" with continent or group of countries
                        #does not take into account changes of weekly variations with religion
                        #details see excel file "temporal-variation-TROTREP_POET_doc_v2_tcm32-14425.xls"
#'year','mon','day','hour'  #time for which emission is calculated; needs to be scalar, not vector
#'arrayTF'              #not used yet
#'plotTF'               #plot different categories
#'spec'                 #"CO2" or "CO", can be vector
#
#chg April 5, 2006
#
#-------------------------------------------------------
#EMISSIONS ON GRID FOR CARBON DIOXIDE IN 2000
#
#Calculated by EDGAR on 22-jul-2005 by laeedg
#calculation name : G: 32FT2000 CO2 BB-AVG EF-32
#(this file created on 25-jul-2005 11:05:42)
#-------------------------------------------------------
#EDGAR Inventory 1x1 EDG 32FT2000 CO2 BB-AVG EF-32
#CARBON DIOXIDE 2000 annual   # cells: 12379 (<>0)
#Values: min: 9.9200e+03 max: 2.9500e+11 sum: 3.6957e+13
#Units : kg CO2 (FMM)/yr       For cells <> 0: avg: 2.9854e+09
#-------------------------------------------------------
#  $Id: get.edgar.r,v 1.3 2007-06-27 11:54:03 skoerner Exp $
#---------------------------------------------------------------------------------------------------



######################################################################################################################
#################################### first without temporal variation ################################################
######################################################################################################################
if(!tempvarTF){
if(any(spec=="CO2")){
edgar<-read.table(paste(epath,"edgar_32ft2000_co2_ant.1x1",sep=""),skip=12,sep=",")
lats<--90:89
lons<--180:179
edgar.mat<-matrix(ncol=length(lons),nrow=length(lats))
edgar.mat[cbind(edgar[,2]-min(lats)+1,edgar[,1]-min(lons)+1)]<-edgar[,3]
edgar.mat[is.na(edgar.mat)]<-0
if(plotTF)image(t(log(edgar.mat)))
}
##########################################################
#make sure have same day, need to use factors for EDGAR
##########################################################


#also CO
#-------------------------------------------------------
#EMISSIONS ON GRID FOR CARBON MONOXIDE IN 2000
#
#Calculated by EDGAR on 27-jul-2005 by laeedg
#calculation name : G: 32FT2000 CO BB-AVG EF-32
#(this file created on 27-jul-2005 16:24:37)
#-------------------------------------------------------
#EDGAR Inventory 1x1 EDGV32FT2000 CO BB-AVG EF-32
#CARBON MONOXIDE 2000 annual   # cells: 12456 (<>0)
#Values: min: 2.5000e+02 max: 4.6200e+09 sum: 1.0410e+12
#Units : kg CO (FMM)/yr       For cells <> 0: avg: 8.3575e+07
#-------------------------------------------------------
if(any(spec=="CO")){
edgar<-read.table(paste(epath,"edgar_32ft2000_co_ant.1x1",sep=""),skip=12,sep=",")
edgar.mat2<-matrix(ncol=length(lons),nrow=length(lats))
edgar.mat2[cbind(edgar[,2]-min(lats)+1,edgar[,1]-min(lons)+1)]<-edgar[,3]
edgar.mat2[is.na(edgar.mat2)]<-0
if(plotTF)image(t(edgar.mat2))
}
#conversion from kg(CO or CO2)/y/gridcell to micromoles/m2-sec
lat1<-lats;lat2<-lats+1;gridareae<-3.1416*6.377e6*6.377e6*(sin(lat2/180*3.1416)-sin(lat1/180*3.1416))/180;
if(any(spec=="CO2"))edgar.mat<-edgar.mat*1000/(3600*24*365)/gridareae/(12+16+16)*1E6
if(any(spec=="CO"))edgar.mat2<-edgar.mat2*1000/(3600*24*365)/gridareae/(12+16)*1E6


#transform to higher resolution

#assume resolution lat-lon finer, so use edgar grid for loop
latbins<-seq(ll.y,ll.y+(numpix.y-1)*yres,by=yres)
lonbins<-seq(ll.x,ll.x+(numpix.x-1)*xres,by=xres)

edgar.x<-floor(lonbins-min(lons)+1) #assume lower left as ref. in edgar
edgar.y<-floor(latbins-min(lats)+1)

x.mat<-rep(1,length(edgar.y))%o%edgar.x
y.mat<-edgar.y%o%rep(1,length(edgar.x))
yx<-cbind(as.vector(y.mat),as.vector(x.mat))

if(any(spec=="CO2"))edgar.eu<-matrix(edgar.mat[yx],nrow=length(latbins),ncol=length(lonbins))
if(any(spec=="CO"))edgar.eu2<-matrix(edgar.mat2[yx],nrow=length(latbins),ncol=length(lonbins))
if(length(spec)==2){
  out<-array(dim=c(dim(edgar.eu),2))
  dimnames(out)[[3]]<-c("CO2","CO")
  out[,,"CO2"]<-edgar.eu
  out[,,"CO"]<-edgar.eu2
  return(out)
}
if(any(spec=="CO2"))out<-edgar.eu
if(any(spec=="CO"))out<-edgar.eu2
return(out)
}#of if(!tempvarTF){

######################################################################################################################
#################################### with temporal variation #########################################################
######################################################################################################################
if(tempvarTF){ #with temporal variation
lats<--90:89
lons<--180:179

#categories from files, CO2:
#substring(dir(paste(epath,"edgar_32ft2000_co2_total_tcm32-21997/",sep="")),1,3)
maps.file<-c("F80","b10","b20","b40","f10","f20","f30","f40","f51","f54","f57","f58","f60","i40")
cats<-c(       4,    5,    1,    3,    5,    1,    4,    3,    7,    8,    9,    4,    4,    4) 
#rem                                                                            no variation, therefore use refineries (4)

#substring(dir(paste(epath,"edgar_32ft2000_co_avg_32_total_tcm32-22060/",sep="")),1,3)
maps.file2<-c("b10","b20","b30","b40","b51","f10","f20","f30","f40","f51","f54","f57","f58","f80","i10","i20","i50","w40")
cats2<-rep(NA,length(maps.file2))
cats2[(substring(maps.file2,1,1)=="b"|substring(maps.file2,1,1)=="f")&substring(maps.file2,2,3)=="10"]<-5
cats2[(substring(maps.file2,1,1)=="b"|substring(maps.file2,1,1)=="f")&substring(maps.file2,2,3)=="20"]<-1
cats2[(substring(maps.file2,1,1)=="b"|substring(maps.file2,1,1)=="f")&substring(maps.file2,2,3)=="30"]<-1
cats2[(substring(maps.file2,1,1)=="b"|substring(maps.file2,1,1)=="f")&substring(maps.file2,2,3)=="40"]<-3
cats2[(substring(maps.file2,1,1)=="b"|substring(maps.file2,1,1)=="f")&substring(maps.file2,2,3)=="51"]<-7
cats2[(substring(maps.file2,1,1)=="b"|substring(maps.file2,1,1)=="f")&substring(maps.file2,2,3)=="54"]<-8
cats2[(substring(maps.file2,1,1)=="b"|substring(maps.file2,1,1)=="f")&substring(maps.file2,2,3)=="57"]<-9
cats2[(substring(maps.file2,1,1)=="b"|substring(maps.file2,1,1)=="f")&substring(maps.file2,2,3)=="58"]<-4
cats2[(substring(maps.file2,1,1)=="b"|substring(maps.file2,1,1)=="f")&substring(maps.file2,2,3)=="80"]<-4
cats2[substring(maps.file2,1,3)=="i10"]<-5
cats2[substring(maps.file2,1,3)=="i20"]<-5
cats2[substring(maps.file2,1,3)=="i50"]<-5
cats2[substring(maps.file2,1,3)=="w40"]<-5

#categories for activities
#   1 Power plants
#   2 Area source combustion "industry"
#   3 Small combustion sources "residential"
#   4 Refineries
#   5 Industrial processes: metallurgy, also food 
#   6 Solvent use
#   7 Transport Road transport (incl. evaporation)
#   8 Transport Transport non-road
#   9 Transport aircraft

monthly_var = matrix(c(                         #  category
      1.2,1.15,1.05,1.0,0.9,0.85,0.8,0.875,0.95,1.0,1.075,1.15,    #  1 Power plants
      1.1,1.075,1.05,1.0,0.95,0.9,0.93,0.95,0.97,1.0,1.025,1.05,   #  2 Area source combustion "industry"
      1.7,1.5,1.3,1.0,0.7,0.4,0.2,0.4,0.7,1.05,1.4,1.65,           #  3 Small combustion sources "residential"
      1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,             #  4 Refineries
      0.98,1.04,1.05,1.04,1.06,1.05,0.97,0.88,1.01,1.04,1.01,0.87, #  5 Industrial processes (EUROPE)
      0.95,0.96,1.02,1.0,1.01,1.03,1.03,1.01,1.04,1.03,1.01,0.91,  #  6 Solvent use
      0.88,0.92,0.98,1.03,1.05,1.06,1.01,1.02,1.06,1.05,1.01,0.93, #  7 Transport
      1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,             #  8 Transport
      1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,             #  9 Transport
#      0.45,1.3,2.35,1.7,0.85,0.85,0.85,1.0,1.1,0.65,0.45,0.45,     # 10 agricult. (NH3-manure)
#      0.45,1.3,2.35,1.7,0.85,0.85,0.85,1.0,1.1,0.65,0.45,0.45,     # 11 agricult. (NH3-manure)
#      0.45,1.3,2.35,1.7,0.85,0.85,0.85,1.0,1.1,0.65,0.45,0.45      # 12 agricult. (NH3-manure)
      ),ncol=9,nrow=12)

weekly_var = matrix(c(    #   category
      1.06,1.06,1.06,1.06,1.06,0.85,0.85,   #   1 Power plants
      1.08,1.08,1.08,1.08,1.08,0.80,0.80,   #   2 Area source combustion
      1.08,1.08,1.08,1.08,1.08,0.80,0.80,   #   3 Small combustion sources
      1.00,1.00,1.00,1.00,1.00,1.00,1.00,   #   4 Refineries
      1.00,1.00,1.00,1.00,1.00,1.00,1.00,   #   5 Industrial processes: metallurgy
#      1.40,1.40,1.40,1.40,1.40,0.00,0.00,   #   5 Industrial processes: food
      1.20,1.20,1.20,1.20,1.20,0.50,0.50,   #   6 Solvent use
      1.02,1.06,1.08,1.10,1.14,0.81,0.79,   #   7 Transport
      1.02,1.06,1.08,1.10,1.14,0.81,0.79,   #   8 Transport
      1.02,1.06,1.08,1.10,1.14,0.81,0.79,   #   9 Transport
#      1.00,1.00,1.00,1.00,1.00,1.00,1.00,   #  10 agricult. (NH3-manure)
#      1.00,1.00,1.00,1.00,1.00,1.00,1.00,   #  11 agricult. (NH3-manure)
#      1.00,1.00,1.00,1.00,1.00,1.00,1.00    #  12 agricult. (NH3-manure)
      ),ncol=9,nrow=7)

hourly_var = matrix(c(                         #  category
      0.79,0.72,0.72,0.71,0.74,0.80,0.92,1.08,1.19,1.22,1.21,1.21, #   1 Power plants
      1.17,1.15,1.14,1.13,1.10,1.07,1.04,1.02,1.02,1.01,0.96,0.88, #
      0.75,0.75,0.78,0.82,0.88,0.95,1.02,1.09,1.16,1.22,1.28,1.30, #   2 Area source combustion
      1.22,1.24,1.25,1.16,1.08,1.01,0.95,0.90,0.85,0.81,0.78,0.75, #
      0.38,0.36,0.36,0.36,0.37,0.50,1.19,1.53,1.57,1.56,1.35,1.16, #   3 Small combustion sources
      1.07,1.06,1.00,0.98,0.99,1.12,1.41,1.52,1.39,1.35,1.00,0.42, #
      1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00, #   4 Refineries
      1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00, #
      1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00, #   5 Industrial processes
      1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00, #
      0.50,0.35,0.20,0.10,0.10,0.20,0.75,1.25,1.40,1.50,1.50,1.50, #   6 Solvent use
      1.50,1.50,1.50,1.50,1.50,1.40,1.25,1.10,1.00,0.90,0.80,0.70, #
      0.19,0.09,0.06,0.05,0.09,0.22,0.86,1.84,1.86,1.41,1.24,1.20, #   7 Transport
      1.32,1.44,1.45,1.59,2.03,2.08,1.51,1.06,0.74,0.62,0.61,0.44, #
      0.19,0.09,0.06,0.05,0.09,0.22,0.86,1.84,1.86,1.41,1.24,1.20, #   8 Transport
      1.32,1.44,1.45,1.59,2.03,2.08,1.51,1.06,0.74,0.62,0.61,0.44, #
      0.19,0.09,0.06,0.05,0.09,0.22,0.86,1.84,1.86,1.41,1.24,1.20, #   9 Transport
      1.32,1.44,1.45,1.59,2.03,2.08,1.51,1.06,0.74,0.62,0.61,0.44, #
#      0.60,0.60,0.60,0.60,0.60,0.65,0.75,0.90,1.10,1.25,1.45,1.60, #  10 agricult. (NH3-manure)
#      1.80,1.75,1.70,1.55,1.35,1.10,0.90,0.75,0.65,0.60,0.60,0.60, #
#      0.60,0.60,0.60,0.60,0.60,0.65,0.75,0.90,1.10,1.25,1.45,1.60, #  11 agricult. (NH3-manure)
#      1.80,1.75,1.70,1.55,1.35,1.10,0.90,0.75,0.65,0.60,0.60,0.60, #
#      0.60,0.60,0.60,0.60,0.60,0.65,0.75,0.90,1.10,1.25,1.45,1.60, #  12 agricult. (NH3-manure)
#      1.80,1.75,1.70,1.55,1.35,1.10,0.90,0.75,0.65,0.60,0.60,0.60
      ),ncol=9,nrow=24)

times<-weekdayhr(rep(year,length(lons)),rep(mon,length(lons)),rep(day,length(lons)),rep(hour,length(lons)),runtt=0,diffGMT=lons/360*24)
year<-times[,"yr"]
mon<-times[,"mon"]
day<-times[,"day"]
hr<-times[,"hr"]
wday<-times[,"weekd"] #monday is first weekday; from local time to GMT
wday[wday==0]<-wday[wday==0]+7 #sunday is 7
if(plotTF)par(mfrow=c(3,5))
if(any(spec=="CO2")){
  if(!existsr("edgar.mat",epath)){
    edgar.mat<-array(dim=c(length(lats),length(lons),length(maps.file))) #for each category one layer, CO2
    edgar.mats<-edgar.mat[,,1]
    for(fnm in maps.file){
      id<-which(maps.file==fnm)
      edgar<-read.table(paste(epath,"edgar_32ft2000_co2_total_tcm32-21997/",fnm,"00co2.1x1",sep=""),skip=12,sep=",")
      edgar.mat[,,id][cbind(edgar[,2]-min(lats)+1,edgar[,1]-min(lons)+1)]<-edgar[,3]
      edgar.mat[is.na(edgar.mat)]<-0
      if(plotTF){image(t(log(edgar.mat[,,id])));title(fnm)}
    }
    assignr("edgar.mat",edgar.mat,epath)
  } else {
    edgar.mat<-getr("edgar.mat",epath)
  }
  if(plotTF){for(fnm in maps.file){id<-which(maps.file==fnm);image(t(log(edgar.mat[,,id])));title(fnm)}}
  edgar.mats<-edgar.mat[,,1]*0;
  for(i in 1:length(cats)){
    time.fact<-monthly_var[mon,cats[i]]*weekly_var[wday,cats[i]]*hourly_var[hour,cats[i]]
    time.fact.sh<-time.fact
    time.fact.eq<-time.fact
    if(any(i==c(1,2,3))){
      time.fact.sh<-monthly_var[(mon+5)%%12+1,cats[i]]*weekly_var[wday,cats[i]]*hourly_var[hour,cats[i]] #southern hemisphere
      time.fact.eq<-1*weekly_var[wday,cats[i]]*hourly_var[hour,cats[i]] #equator
    }
    tmp<-t(t(edgar.mat[,,i])*time.fact)
    tmp[lats<0,]<-t(t(edgar.mat[,,i])*time.fact.sh)[lats<0,]
    tmp[abs(lats)<=20,]<-t(t(edgar.mat[,,i])*time.fact.eq)[abs(lats)<=20,]
    edgar.mats<-edgar.mats+tmp
  }
}
##########################################################
#make sure have same day, need to use factors for EDGAR
##########################################################

#image(t(log(edgar.mats)))
#also CO
#-------------------------------------------------------
#EMISSIONS ON GRID FOR CARBON MONOXIDE IN 2000
#
#Calculated by EDGAR on 27-jul-2005 by laeedg
#calculation name : G: 32FT2000 CO BB-AVG EF-32
#(this file created on 27-jul-2005 16:24:37)
#-------------------------------------------------------
#EDGAR Inventory 1x1 EDGV32FT2000 CO BB-AVG EF-32
#CARBON MONOXIDE 2000 annual   # cells: 12456 (<>0)
#Values: min: 2.5000e+02 max: 4.6200e+09 sum: 1.0410e+12
#Units : kg CO (FMM)/yr       For cells <> 0: avg: 8.3575e+07
#-------------------------------------------------------
if(plotTF)par(mfrow=c(4,5))
if(any(spec=="CO")){
  #fix inconsistent name for one category
  unix(paste("mv ",epath,"/edgar_32ft2000_co_avg_32_total_tcm32-22060/b2000co_.1x1 ",epath,"/edgar_32ft2000_co_avg_32_total_tcm32-22060/b2000co.1x1",sep=""))
  if(!existsr("edgar.mat2",epath)){
    edgar.mat2<-array(dim=c(length(lats),length(lons),length(maps.file2))) #for each category one layer, CO
    for(fnm in maps.file2){
      id<-which(maps.file2==fnm)
      edgar<-read.table(paste(epath,"edgar_32ft2000_co_avg_32_total_tcm32-22060/",fnm,"00co.1x1",sep=""),skip=12,sep=",")
      edgar.mat2[,,id][cbind(edgar[,2]-min(lats)+1,edgar[,1]-min(lons)+1)]<-edgar[,3]
      edgar.mat2[is.na(edgar.mat2)]<-0
    }
    assignr("edgar.mat2",edgar.mat2,epath)
  } else{
    edgar.mat2<-getr("edgar.mat2",epath)
  }
  edgar.mats2<-edgar.mat2[,,1]*0;
  if(plotTF){for(fnm in maps.file2){id<-which(maps.file2==fnm);image(t(log(edgar.mat2[,,id])));title(fnm)}}
  for(i in 1:length(cats2)){
    time.fact<-monthly_var[mon,cats2[i]]*weekly_var[wday,cats2[i]]*hourly_var[hour,cats2[i]]
    time.fact.sh<-time.fact
    time.fact.eq<-time.fact
    if(any(i==c(1,2,3))){
      time.fact.sh<-monthly_var[(mon+5)%%12+1,cats2[i]]*weekly_var[wday,cats2[i]]*hourly_var[hour,cats2[i]] #southern hemisphere, shift 6 months
      time.fact.eq<-1*weekly_var[wday,cats2[i]]*hourly_var[hour,cats2[i]] #equator, set to constant seasonally
    }
    tmp<-t(t(edgar.mat2[,,i])*time.fact)
    tmp[lats<0,]<-t(t(edgar.mat2[,,i])*time.fact.sh)[lats<0,]
    tmp[abs(lats)<=20,]<-t(t(edgar.mat2[,,i])*time.fact.eq)[abs(lats)<=20,]
    edgar.mats2<-edgar.mats2+tmp
  }
#    if(plotTF){image(t(log(edgar.mat2[,,id])))}
}
#conversion from kg(CO or CO2)/y/gridcell to micromoles/m2-sec
lat1<-lats;lat2<-lats+1;gridareae<-3.1416*6.377e6*6.377e6*(sin(lat2/180*3.1416)-sin(lat1/180*3.1416))/180;
if(arrayTF){
  if(any(spec=="CO2"))edgar.mat<-edgar.mat*1000/(3600*24*365)/gridareae/(12+16+16)*1E6
  if(any(spec=="CO"))edgar.mat2<-edgar.mat2*1000/(3600*24*365)/gridareae/(12+16)*1E6
}
if(any(spec=="CO2"))edgar.mats<-edgar.mats*1000/(3600*24*365)/gridareae/(12+16+16)*1E6
if(any(spec=="CO"))edgar.mats2<-edgar.mats2*1000/(3600*24*365)/gridareae/(12+16)*1E6

if(arrayTF){ # return array with different emission maps for different categories
  if(length(spec)==2){
    out<-list(edgar.mat,edgar.mat2)
    names(out)<-c("CO2","CO")
    return(out)
  }
  if(any(spec=="CO2"))out<-edgar.mat
  if(any(spec=="CO"))out<-edgar.mat2
}

# only use sum of emission categories
if(any(spec=="CO2"))edgar.mat<-edgar.mats
if(any(spec=="CO"))edgar.mat2<-edgar.mats2
#transform to higher resolution

#assume resolution lat-lon finer, so use edgar grid for loop
latbins<-seq(ll.y,ll.y+(numpix.y-1)*yres,by=yres)
lonbins<-seq(ll.x,ll.x+(numpix.x-1)*xres,by=xres)

edgar.x<-floor(lonbins-min(lons)+1) #assume lower left as ref. in edgar
edgar.y<-floor(latbins-min(lats)+1)

x.mat<-rep(1,length(edgar.y))%o%edgar.x
y.mat<-edgar.y%o%rep(1,length(edgar.x))
yx<-cbind(as.vector(y.mat),as.vector(x.mat))

if(any(spec=="CO2"))edgar.eu<-matrix(edgar.mat[yx],nrow=length(latbins),ncol=length(lonbins))
if(any(spec=="CO"))edgar.eu2<-matrix(edgar.mat2[yx],nrow=length(latbins),ncol=length(lonbins))
if(length(spec)==2){
  out<-array(dim=c(dim(edgar.eu),2))
  dimnames(out)[[3]]<-c("CO2","CO")
  out[,,"CO2"]<-edgar.eu
  out[,,"CO"]<-edgar.eu2
  return(out)
}
if(any(spec=="CO2"))out<-edgar.eu
if(any(spec=="CO"))out<-edgar.eu2
return(out)
}#of if(tempvarTF){

}



