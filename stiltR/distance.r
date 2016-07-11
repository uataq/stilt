distance<-function(x1, x2, y1, y2)
{
#Calculates the distance [km] between 2 locations, given their lat/lon coordinates
#'x1'&'x2' are locations in deg-LON
#'y1'&'y2' are locations in deg-LAT
#NOTE:  the linear distance given by differences in LON (which is a function of LAT) is
#       approximated by the MEAN of the two LATs
	#NOTE: not identical with Fltplan/distance.r
#4/10/2000 by JCL
#10/30/01 cy CHG: changed from taking mean(lat1,lat2) to (lat1+lat2)/2 to allow vector processing
#
#  $Id: distance.r,v 1.2 2007-06-27 11:54:03 skoerner Exp $
#---------------------------------------------------------------------------------------------------

   Rearth <- 6371
   #mean radius of earth [km]
   rr <- cos(((y1 + y2)/2 * pi)/180)
   #factor by which distance in degs lon has to be multiplied by to be comparable to degs lat
   dx <- rr * (x1 - x2)
   dy <- y1 - y2
   result <- (Rearth * sqrt(dx^2 + dy^2) * pi)/180
   return(result)
}
