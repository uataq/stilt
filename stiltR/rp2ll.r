rp2ll<-function(r,p,lon0=-72.172,lat0=42.536){
#function to convert r, phi (distance and angle) to lon,lat
#'r'	distance from center
#'p'	angle (north=0, east=90)
#'lon0' center longitude (origin)
#'lat0'	center latitude (origin)
#
#returns lat and lon
#
#  $Id: rp2ll.r,v 1.2 2007-06-27 11:54:03 skoerner Exp $
#---------------------------------------------------------------------------------------------------

  Rearth <- 6371 #mean radius of earth [km]
  lat<-lat0+r/(Rearth*pi/180)*cos(p*pi/180)
  lon<-lon0+r/(cos(((lat0+lat)/2 * pi)/180)*Rearth*pi/180)*sin(p*pi/180) #r in km, need distance in deg lon: 116,116
  pnew<-atan((lon-lon0)/(lat-lat0))*180/pi
#  print(cbind(p,pnew))
  #need to iterate...
  #distance is fine, but angle not
  if(length(r)==1)return(c(lon,lat))
  if(length(r)>1)return(cbind(lon,lat))
}
