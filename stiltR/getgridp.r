getgridp<-function(min.x, max.x, min.y, max.y, numpix.x, numpix.y,coarse.factor=1)
{
#Implements decisions in trajectory model about which emission grid to use
#  by taking in vectors of min & max grid points
#Each element is a time point
#Returns vector of filenames (char) that would give correct file for emission grid
#8/20/99 by JCL
#CHG 5/3/01 changes to use function for particles: ran.x<-cummax(ran.x)
#CHG 4/01/02 changed to allow different resolution: coarse = 1,2,4,8,16,32 for 20, 40, 80, 160, 320, 640 km
#CHG 11/06/03 changed to allow high resolution: coarse = 0 for 20 km
#coarse.factor: to allow different resolution: coarse = 1,2,4,8,16,32 for 20, 40, 80, 160, 320, 640 km; coarse = 0 for 20 km only
#
#  $Id: getgridp.r,v 1.4 2008-10-30 09:38:39 gerbig Exp $
#---------------------------------------------------------------------------------------------------

   leng <- length(min.x)
   #number of elements, or timepoints
   #this will come in handy later
   ran.x <- max.x - min.x
   ran.y <- max.y - min.y
	ran.x<-cummax(ran.x)#don't allow resolution to get finer at earlier btimes
	ran.y<-cummax(ran.y)
	#adapt coarse.factor
   cf<-c(1, 2, 4, 8,16,32)
	mr<-c(0, 6,12,24,48,96)
	mr<-cbind(cf,mr)
	mr<-mr[mr[,"cf"]==coarse.factor,"mr"]
	ran.x[ran.x<mr]<-mr*rep(1,leng)[ran.x<mr]
	ran.y[ran.y<mr]<-mr*rep(1,leng)[ran.y<mr]
   #calculate ranges
   minranx <- c(0, 0, 6, 6, 6, 12, 12, 12, 24, 24, 24, 48, 48, 48, 96, 96)
   minrany <- c(0, 6, 0, 6, 12, 6, 12, 24, 12, 24, 48, 24, 48, 96, 48, 96)
   gridname <- rep(0, leng)
   #Loop through each of 16 elements of 'minranx' & 'maxrany' one at a time
   #'gridname' represents different resolutions of emission grid--16 is coarsest & 1 is finest
   for(i in 1:length(minranx)) {
           minranxTF <- ran.x >= minranx[i]
           minranyTF <- ran.y >= minrany[i]
           gridname[minranxTF & minranyTF] <- i
   }
	if(coarse.factor==0)gridname<-rep(1, leng) #use high resolution for coarse.factor 0
   #xpart & ypart can range betw. 0~3.These tell you which piece of whole grid is needed.
   #If value is 0, then means that ENTIRE grid is used.
   #Note that only when resolution is fine (gridname<=3) is emission grid broken down.
   #Note also that broken down emission grids can OVERLAP.
   xpart <- rep(0, leng)
   ypart <- rep(0, leng)
   first.x <- (min.x >= 1) & (max.x <= 0.5 * numpix.x)
   second.x <- (min.x > 0.25 * numpix.x) & (max.x <= 0.75 * numpix.x)
   third.x <- (min.x > 0.5 * numpix.x) & (max.x <= 1 * numpix.x)
   first.y <- (min.y >= 1) & (max.y <= 0.5 * numpix.y)
   second.y <- (min.y > 0.25 * numpix.y) & (max.y <= 0.75 * numpix.y)
   third.y <- (min.y > 0.5 * numpix.y) & (max.y <= 1 * numpix.y)
   largedx <- !(first.x | second.x | third.x)
   largedy <- !(first.y | second.y | third.y)
   sel <- gridname == 1
   xpart[first.x & sel] <- 1
   xpart[second.x & sel] <- 2
   xpart[third.x & sel] <- 3
   ypart[first.y & sel] <- 1
   ypart[second.y & sel] <- 2
   ypart[third.y & sel] <- 3
   sel <- gridname == 2
   xpart[first.x & sel] <- 1
   xpart[second.x & sel] <- 2
   xpart[third.x & sel] <- 3
   ypart[(!largedx) & sel] <- 1
   #if deltax is small enough, then ypart must = 1
   sel <- gridname == 3
   ypart[first.y & sel] <- 1
   ypart[second.y & sel] <- 2
   ypart[third.y & sel] <- 3
   xpart[(!largedy) & sel] <- 1
   #if deltay is small enough, then xpart must = 1
   #deltax or deltay is too large, so have to use coarse grid
   xpart[largedx | largedy] <- 0
   ypart[largedx | largedy] <- 0
   result <- cbind(xpart, ypart, gridname)
   dimnames(result) <- list(NULL, c("xpart", "ypart", "gridname"))
   return(result)
}

