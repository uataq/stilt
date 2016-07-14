getmetfile<-function(yr=0,mon,day,hr,nhrs,metd="edas",doublefiles=F){
#A function that is smart enough to return the meteorological filename to be used in 
#	driving particle dispersion model when given the yr, mon, day, hr, & nhrs
#Arguments are all SINGLE values
#Note that function returns TWO FILES if magnitude of nhrs is large and run time would overlap two files
#Note that function only works right now for running model BACKWARDS (i.e., NEGATIVE nhrs)
#'yr' 2 digit (!) year (starting time)
#'mon' month (starting time)
#'day' day (starting time)
#'hr' hour (starting time)
#'nhrs' is the number of hours that particle model would be run 
#	if nhrs<0, then running model BACKWARDS
#'metd' meteorological data identifier, such as "edas", or "fnl" or "brams"
#'doublefiles' should concatenated met files be used? Allows starting times between files
#
#1/20/2000 by JCL, modifications by CHG
#
#  $Id: getmetfile.r,v 1.19 2008-08-14 12:18:02 gerbig Exp $
#---------------------------------------------------------------------------------------------------



#Determine whether to use NGM (lower reso) or EDAS (higher reso) met files
#internaly: change yr from numeric to character

if(metd=="brams"){
#get all days involved
#assume daily files
  ndays<-floor((abs(nhrs)+24-hr)/24)
  if(yr<50){yr<-yr+2000} else if(yr<100) {yr<-yr+1900}
  newtime<-weekdayhr(yr,mon,day,hr,sign(nhrs)*60*24*(0:ndays))
  day<-newtime[,"day"]; mon<-newtime[,"mon"]; yr<-newtime[,"yr"]; hr<-newtime[,"hr"];
  day<-paste(substring(100+abs(day),2,3),sep="") #2 digits 
  yr<-paste(substring(10000+abs(yr),2,5),sep="") #2 digits 
  mon<-paste(substring(100+abs(mon),2,3),sep="") #2 digits 
  metfile<-paste("brams_",mon,"_",day,"_",yr,"_1.bin",sep="")
  return(metfile)
}#if rams

if(length(grep('wrf|d[0-9][0-9]',metd)) == 1) {
#get all days involved
#assume daily files
  dnn <- if (metd=="wrf") "d01" else metd
  hroffset <- 6 #time at which daily files start: 06z - 06z the next day
  ndays<-floor((abs(nhrs)+24-(hr-hroffset))/24)
  seq.ndays <- 0:ndays
  if (hr < hroffset && nhrs > 0) seq.ndays <- c(-1,seq.ndays)
  if(yr<50){yr<-yr+2000} else if(yr<100) {yr<-yr+1900}
  newtime<-weekdayhr(yr,mon,day,hr,sign(nhrs)*60*24*seq.ndays)
  day<-newtime[,"day"]; mon<-newtime[,"mon"]; yr<-newtime[,"yr"]; hr<-newtime[,"hr"];
  jdat<-julian(mon,day,yr)
  day<-paste(substring(100+abs(day),2,3),sep="") #2 digits 
  yr<-paste(substring(10000+abs(yr),2,5),sep="") #4 digits 
  mon<-paste(substring(100+abs(mon),2,3),sep="") #2 digits 
  metfile<-paste(dnn,".",yr,mon,day,".arl",sep="")
  return(metfile)
}#if rams


if(metd=="alad"){ #also allow for aladin model
#get all days involved
#assume daily files
  ndays<-floor((abs(nhrs)+24-hr)/24)
  if(yr<50){yr<-yr+2000} else if(yr<100) {yr<-yr+1900}
  newtime<-weekdayhr(yr,mon,day,hr,sign(nhrs)*60*24*(0:ndays))
  day<-newtime[,"day"]; mon<-newtime[,"mon"]; yr<-newtime[,"yr"]; hr<-newtime[,"hr"];
  jdat<-julian(mon,day,yr)
  sel<-jdat>=julian(5,15,2005)&jdat<=julian(6,26,2005)
  day<-day[sel];mon<-mon[sel];yr<-yr[sel];hr<-hr[sel]
  day<-paste(substring(100+abs(day),2,3),sep="") #2 digits 
  yr<-paste(substring(10000+abs(yr),2,5),sep="") #2 digits 
  mon<-paste(substring(100+abs(mon),2,3),sep="") #2 digits 
  metfile<-paste("alad",yr,mon,day,sep="")
  return(metfile)
}#if alad

if(metd=="ecmw"){ #also allow for ECMWF fields (daily patched short term forecasts)
#get all days involved
#assume daily files
  ndays<-floor((abs(nhrs)+24-hr)/24)
  if(yr<50){yr<-yr+2000} else if(yr<100) {yr<-yr+1900}
  newtime<-weekdayhr(yr,mon,day,hr,sign(nhrs)*60*24*(0:ndays))
  day<-newtime[,"day"]; mon<-newtime[,"mon"]; yr<-newtime[,"yr"]; hr<-newtime[,"hr"];
  jdat<-julian(mon,day,yr)
  sel<-jdat>=julian(5,04,2005)&jdat<=julian(6,26,2005)
  day<-day[sel];mon<-mon[sel];yr<-yr[sel];hr<-hr[sel]
  day<-paste(substring(100+abs(day),2,3),sep="") #2 digits 
  yr<-paste(substring(10000+abs(yr),4,5),sep="") #2 digits 
  mon<-paste(substring(100+abs(mon),2,3),sep="") #2 digits 
  metfile<-paste("ecmw.",yr,mon,day,sep="")
  return(metfile)
}#if rams

if(metd=="ECmetF"){ #also allow for ECMWF fields (monthly patched short term forecasts)
#get all days involved
  back<-F
  if(nhrs<0) back<-T
  ndays<-floor((abs(nhrs)+24-hr)/24)
  if(yr<50){yr<-yr+2000} else if(yr<100) {yr<-yr+1900}
  newtime<-weekdayhr(yr,mon,day,hr,sign(nhrs)*60*24*(0:ndays))
  day<-newtime[,"day"]; mon<-newtime[,"mon"]; yr<-newtime[,"yr"]; hr<-newtime[,"hr"];
  jdat<-julian(mon,day,yr)
#fix for FEB 06, change in Resolution starting on 2nd!
#  sel<-yr==2006&mon==2&day==1
#  day[sel]<-31;mon[sel]<-1
  day<-paste(substring(100+abs(day),2,3),sep="") #2 digits 
  yr<-paste(substring(10000+abs(yr),4,5),sep="") #2 digits 
  mon<-paste(substring(100+abs(mon),2,3),sep="") #2 digits 
  yrmns<-unique(paste(yr,mon,sep=""))
  metfile<-paste("ECmetF.",yrmns,"01","00",".arl",sep="")
#  metfile[metfile=="ECmetF.06020100.arl"]<-"ECmetF.06020200.arl" #fix for FEB 06, change in Resolution starting on 2nd!
  return(metfile)
}#if rams

#-----------------------------------------------------------------------------------------------------------------
if(metd=="GDAS"){ #allow NCEP GDAS-arl.files  
#get all days involved
  back<-F
  if(nhrs<0) back<-T
  ndays<-floor((abs(nhrs)+24-hr)/24)
  if(yr<50){yr<-yr+2000} else if(yr<100) {yr<-yr+1900}
  newtime<-weekdayhr(yr,mon,day,hr,sign(nhrs)*60*24*(0:ndays))
  day<-newtime[,"day"]; mon<-newtime[,"mon"]; yr<-newtime[,"yr"]; hr<-newtime[,"hr"];
  jdat<-julian(mon,day,yr)
#fix for FEB 06, change in Resolution starting on 2nd!
#  sel<-yr==2006&mon==2&day==1
#  day[sel]<-31;mon[sel]<-1
  day<-paste(substring(100+abs(day),2,3),sep="") #2 digits 
  yr<-paste(substring(10000+abs(yr),4,5),sep="") #2 digits 
  mon<-paste(substring(100+abs(mon),2,3),sep="") #2 digits 
  yrmns<-unique(paste(yr,mon,sep=""))
  metfile<-paste("ECmetF.",yrmns,"01","00",".arl",sep="")
#  metfile[metfile=="ECmetF.06020100.arl"]<-"ECmetF.06020200.arl" #fix for FEB 06, change in Resolution starting on 2nd!
  return(metfile)
}#if gdas



if(doublefiles){#make sure that gat proper files when starting close to file end or start
#get filename for metdata; add >3h (EDAS) or >6h (FNL) to be able to interpolate between timesteps;
        if(!is.na(match("edas",metd)))
                timejump<-6
        else if(!is.na(match("fnl",metd)))
                timejump<-9
        else if(!is.na(match("edas40",metd)))
                timejump<-3
        else if(!is.na(match("gdas",metd)))
                timejump<-3
        else if(!is.na(match("narr",metd)))
                timejump<-3
        else
                stop("can't use double files for this type")
        newtime<-weekdayhr(yr,mon,day,hr,-sign(nhrs)*60*timejump)[1,]
        day<-newtime["day"]
        mon<-newtime["mon"]
        yr<-newtime["yr"]
        hr<-newtime["hr"]
        nhrs=nhrs+sign(nhrs)*60*(timejump+3)
}

if(yr==0){yr="00"}else if(yr<10){yr<-paste("0",as.character(yr),sep="")}else{ yr<-as.character(yr)}

if(as.numeric(yr)>96){dattype<-"edas.subgrd."}
else if(as.numeric(yr)<50){dattype<-"edas.subgrd."}
	else{dattype<-"ngm."}

metfiles<-paste(dattype,c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec"),yr,sep="")

#number of days in each month
numdays<-c(31,28,31,30,31,30,31,31,30,31,30,31)
if(leap.year(1900+as.numeric(yr)))numdays[2]<-29
	
if(day<=15){suffix<-".001"}else{suffix<-".002"}	#determine suffix to met file name
metfile<-paste(metfiles[mon],suffix,sep="")   #name of meteorological file

#determine if only one, or two meteorological files would be needed in simulation--needed if nhrs is large
numhrs<-24*(day-1)+hr   #number of hrs since beginning of file 
if(suffix==".002")numhrs<-24*(day-15-1)+hr    #if in 2nd half of month
doub<-F;back<-F
if(nhrs<0) back<-T
if(nhrs<0&(numhrs+nhrs)<0) doub<-T	#running model BACKward & need two meteorological files
if(nhrs>0&(numhrs+nhrs)>(15*24)) doub<-T #running model FORward & need two meteorological files

#y2k resistant
yrm1<-as.character(as.numeric(yr)-1); yrp1<-as.character(as.numeric(yr)+1)
if(yr=="00"){yrm1<-"99";yrp1<-"01"}else if(yr=="99"){yrp1<-"00"}else if(yr=="09"){yrp1<-"10"}
else if(substring(yr,1,1)=="0"){yrm1<-paste("0",as.character(as.numeric(yr)-1),sep="");yrp1<-paste("0",as.character(as.numeric(yr)+1),sep="")}
else{  #dummy statement--to satisfy 'else if' construct
}


if(back&doub&(mon==1)&(suffix==".001")){
		metfile<-c(metfile,paste(dattype,"dec",yrm1,".002",sep="")) 	#special case:e.g.,  'ngm_jan96.001' & 'ngm_dec95.002'			
}else if(back&doub&(suffix==".001")){
		metfile<-c(metfile,paste(metfiles[mon-1],".002",sep=""))		#e.g., 'ngm_apr96.001' & 'ngm_mar96.002'
}else if(back&doub&(suffix==".002")){
		metfile<-c(metfile,paste(metfiles[mon],".001",sep="")) 		#e.g., 'ngm_jul96.002' & 'ngm_jul96.001'
}else{  #dummy statement--to satisfy 'else if' construct
}


if(!back&doub&(mon==12)&(suffix==".002")){
		metfile<-c(metfile,paste(dattype,"jan",yrp1,".001",sep="")) 	#special case:e.g.,  'ngm_dec95.002' & 'ngm_jan96.001'			
}else if(!back&doub&(suffix==".002")){
		metfile<-c(metfile,paste(metfiles[mon+1],".001",sep=""))		#e.g., 'ngm_apr96.002' & 'ngm_may96.001'
}else if(!back&doub&(suffix==".001")){
		metfile<-c(metfile,paste(metfiles[mon],".002",sep="")) 		#e.g., 'ngm_jul96.001' & 'ngm_jul96.002'
}else{  #dummy statement--to satisfy 'else if' construct
}

if(doublefiles){ #files pasted together for 2 half month periods; get in single name
  if(length(metfile)>1){
    if(back){metfile<-paste(substring(metfile[2],1,15),substring(metfile[1],13,17),sep="")
    }else{metfile<-paste(substring(metfile[1],1,15),substring(metfile[2],13,17),sep="")}
  }else{metfile<-paste(substring(metfile,1,15),substring(metfile,13,17),sep="")}
}

if(metd=="fnl"|metd=="fnl.nh")metfile<-gsub(dattype,"fnl.nh.",metfile)
if(metd=="fnl.sh")metfile<-gsub(dattype,"fnl.sh.",metfile)
if(metd=="edas40")metfile<-gsub(dattype,"edas40.",metfile)
if(metd=="gdas")metfile<-gsub(dattype,"gdas.",metfile)
if(metd=="narr")metfile<-gsub(dattype,"narr.",metfile)

return(metfile)

}

