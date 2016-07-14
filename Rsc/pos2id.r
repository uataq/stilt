pos2id<-function(jultime,lat,lon,alt,sep="x"){
#function to create identifying label for single receptor (location&time)
#expects jultime as fractional julian day since 1/1/1960
#expects alt as altitude above ground in meters
#example:
# pos2id(15555.4166667,45,-90,0.03)
#[1] "2002x08x03x10x+45.00x+090.00x00030"
#
#  $Id: pos2id.r,v 1.5 2007-06-27 11:54:03 skoerner Exp $
#---------------------------------------------------------------------------------------------------

#time
if(is.matrix(jultime)){ #allows Timesname object of starting times and locations to be used as single argument
  lat<-jultime[,"lat"]
  lon<-jultime[,"lon"]
  alt<-jultime[,"agl"]
  jultime<-jultime[,"fjul"]
}
hr<-round((jultime-floor(jultime))*24)
jultime[hr==24]<-ceiling(jultime[hr==24]);hr[hr==24]<-0 #go to next day
jultime<-month.day.year(floor(jultime))
yr4<-jultime$year #4 digit year
mon<-jultime$month
day<-jultime$day

#need leading zeros
hr<-substring(100+hr,2) #2 digit hr
yr4<-substring(10000+yr4,2) #4 digit year
mon<-substring(100+mon,2) #2 digit mon
day<-substring(100+day,2) #2 digit day

#location
lat<-round(lat,2) #1 km roundoff
lon<-round(lon,2) #1 km roundoff
alt<-round(alt) #1 m roundoff
#no below ground...
alt[alt<0]<-0
#need leading zeros and sign
lats<-rep("S",length(lat))
lons<-rep("W",length(lat))
lats[lat>=0]<-"N"
lons[lon>=0]<-"E"
lat<-paste(substring(100.001+abs(lat),2,6),lats,sep="") #2 dig. before decimal + sign
lon<-paste(substring(1000.001+abs(lon),2,7),lons,sep="") #3 dig. before decimal + sign
alt<-substring(100000+alt,2)

paste(yr4,mon,day,hr,lat,lon,alt,sep=sep) #put everything together
}
