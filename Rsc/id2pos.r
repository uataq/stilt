id2pos<-function(id,sep="x"){
#revers of id2pos
#function to create read identifying label for single receptor (location&time)
#returns time as fractional julian day since 1/1/1960
#and alt as altitude above ground in meters
#example:
# id2pos("2002x08x03x10x+45.00x+090.00x00030")
#[1] 15555.42    45.00    90.00    30.00
#3/4/2004 by CHG
#
#  $Id: id2pos.r,v 1.4 2007-06-27 11:54:03 skoerner Exp $
#---------------------------------------------------------------------------------------------------

#get format using seperator location
cpos<-regexpr(sep,id) #position of first sep
yr4<-as.numeric(substring(id,1,cpos-1));id<-substring(id,cpos+1);cpos<-regexpr(sep,id)
mon<-as.numeric(substring(id,1,cpos-1));id<-substring(id,cpos+1);cpos<-regexpr(sep,id)
day<-as.numeric(substring(id,1,cpos-1));id<-substring(id,cpos+1);cpos<-regexpr(sep,id)
hr<-as.numeric(substring(id,1,cpos-1));id<-substring(id,cpos+1);cpos<-regexpr(sep,id)
lat<-as.numeric(substring(id,1,cpos-2));
if(substring(id,cpos-1,cpos-1)=="S")lat<-(-lat)
id<-substring(id,cpos+1);cpos<-regexpr(sep,id)
lon<-as.numeric(substring(id,1,cpos-2));
if(substring(id,cpos-1,cpos-1)=="W")lon<-(-lon)
id<-substring(id,cpos+1) 
if(regexpr("[a-z]",id)-1>0)id<-substring(id,1,regexpr("[a-z]",id)-1) #remove any trailing non-numeric letters
alt<-as.numeric(id)

time<-julian(mon,day,yr4)+hr/24
pos<-c(time,lat,lon,alt)
names(pos)<-c("time","lat","lon","alt")
return(pos)
}
