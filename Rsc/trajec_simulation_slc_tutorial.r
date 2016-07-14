##########################################################################################################################################
##########################################################################################################################################
#####                                                                                                         ############################
#####   This script will be used to generate trajectories for some time peroid and store the RDATA output     ############################
#####   files in the appropriate directory. These trajectories will use WRF output and the ZSG levels         ############################
#####   from this output.                                                                                     ############################
#####                                                                                                         ############################
#####   Created by DVM on 02/25/2014                                                                          ############################
#####                                                                                                         ############################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################



library(chron)
library(ncdf4)
library(maps)


#First we reference other stiltR functions that we may use and the path for the model output and input. Also opens all libraries that are used within this
#script

sourcepath<-"/uufs/chpc.utah.edu/common/home/u0791983/STILT_modeling/stiltR/";source("sourceall.r")

path     <-"/uufs/chpc.utah.edu/common/home/u0791983/STILT_modeling/Output/"             #model output
metpath  <-"/uufs/chpc.utah.edu/common/home/u0791983/met/"                                     #where met data are stored in ARL format
rundir   <-"/uufs/chpc.utah.edu/common/home/u0791983/STILT_modeling/STILT_Exe/"                                         #specifies main directory where different directories are found to run hymodelc
ncdfpath <-path 			      #path for netcdf files

zsg.name <- NULL


#Please enter a starting date and time in the following format: "YY/DD/YYYY", and "HH:MM:SS" along with the end time and time
#Also enter the time increment (for now this NEEDS to be in 6-hourly increments)

start.date = "04/28/2010"
start.hour = "00:00:00"

end.date = "04/28/2010"
end.hour = "06:00:00"


#Use the Chron module which creates an array of numbers that we will loop through to determine the starting date and
#times for the trajec.r portion of the STILT model

start.time <- chron(start.date, start.hour)
end.time   <- chron(end.date, end.hour)

run.times <- seq(start.time, end.time, by = times("01:00:00"))      #Comboing it up with using the Chron function to get us the dates that we want #to create a sequence of!
run.times <- as.character(run.times)	                            #Sets string to a character so we can manipulate it in the next section


# Substrings our time stamps by year, month, day, and hour

month  = substr(run.times, 2, 3)
days   = substr(run.times, 5, 6)
hour   = substr(run.times, 11, 12)
year   = substr(run.times, 8, 9)
year   = as.numeric(year)
year   = year + 2000
year   = as.character(year)

time.matrix = cbind(year,month,days,hour)                           #Time matrix which has all of our time variables (year, month, day, hour)



cat("\n"); cat("\n"); cat("\n")
cat("Starting WRF-STILT simulation loop, this could take some time.")
cat("\n"); cat("\n"); cat("\n")

for(time in 1:nrow(time.matrix)){


	# Prints time that we are working on

	cat("Working on time: "); cat(time.matrix[time,"year"]);
	cat(time.matrix[time,"month"]); cat(time.matrix[time,"days"]);
	cat(time.matrix[time,"hour"]);cat("\n"); cat("\n")



	##############################################################################################################
	####### PARTICLE LOCATIONS ###################################################################################
	####### control parameters for what STILT should do: Particle locations, mixing ratios or footprints #########
	##############################################################################################################


	#Picks the variables that we want to output in the .RData file. Ceck the trajec.r function
	#for addtional variables that can be outputted


        varstrajec<-c("time","index","lat","lon","agl","grdht","foot","sampt","dmass","zi","pres","sphu","temp","temp0","shtf","lhtf","sigmaw","dens","tcld","totrain","convrain","wbar")



	#----- basic parameters -----#

	nhrs<- -12                       #for how many hours to run (time-backward runs negative).
	veght<-0.5              	 #surface layer (for fluxes); if less than 1 then as fraction of mlhgt; 0.5 is a good value.
	convect<-F              	 #T for convection (RAMS winds: grell convection scheme, EDAS and FNL: simple redistribution within vertical range with CAPE>0)
	stepsize<-0            		 #Enforces Courant criterion. When >0: maximum horizontal distances travelled by a particle in a single timestep, in km.
        convect<-F     		         #For dynamic resolution, choose value 0. First 12 hrs: 20 km; then 60 km
	winderrTF<-F           		 #transport error due to wind error included as stochastic process?
        overwrite=T                      #T: rerun hymodelc, even if particle location object found; F: re-use previously calculated particle location object
	nummodel<-0               	 #which copy of the model to run?  Will be in directory named paste(rundir,"Copy",nummodel,sep="")
	delt<-2                  	 #fixed timestep [min]; set =0 for dynamic timestep
        mgmin<-2000                      #Set to 2000 to keep particles from dropping out of high res WRF runs. ONLY USE WITH WRF otherwrise comment this out, and delete from traj function below!
	numpar<-1000                     #number of particles for STILT simulation



	#Sets specific met file name

	metfile<-c("wrfout_d02.arl","wrfout_d01.arl","gdas1.apr10.arl")


	#----- starting times & locations -----#

	yr <-as.numeric(time.matrix[time,"year"])-2000        #Year needs to be in form of years "before" or "after" 2000
	mon<-as.numeric(time.matrix[time,"month"])
	day<-as.numeric(time.matrix[time,"days"])
	hr <-as.numeric(time.matrix[time,"hour"])

	lat<-40.75;lon<--111.87;agl<-5


	#This command runs the STILT trajec.r function. Also Commented out
	#maxdist=stepsize as there was an issue with this particulation function

        cat("Running hymodelC executable..."); cat("\n"); cat("\n")

	info<-Trajec(yr=yr,mon=mon,day=day,hr=hr,lat=lat,lon=lon,agl=agl,nhrs=nhrs,doublefiles=T,metd=c("fnl","awrf"),delt=delt,numpar=numpar,mgmin=mgmin,veght=veght,metfile=metfile,nummodel=nummodel,metlib=metpath,conv=convect,overwrite=overwrite,outpath=path,varsout=varstrajec,rundir=rundir,zsg.name=zsg.name)


	dat<-getr(info["outname"],path=path)


	timestamp <- paste(time.matrix[time,"year"],formatC(time.matrix[time,"month"],width=2,flag="0"),formatC(time.matrix[time,"days"],width=2,flag="0"),formatC(time.matrix[time,"hour"],width=2,flag="0"),sep="")


	##############################################################################################
	##############################################################################################
	###########   FOOTPRINTS AND INFLUENCE:                                                  #####
	###########   Output will be a 3 D array of Influence (or Surface Influence/Footprints)  #####
	###########   (i.e., sensitivity of atmospheric concentrations to surface fluxes)        #####
	###########   these objects are given names that reflect starting location and time      #####
	###########   e.g. .RDatafoot2002x08x01x00x42.54Nx072.17Wx00030                          #####
	##############################################################################################
	##############################################################################################

	cat("Generating footprints..."); cat("\n"); cat("\n")

	#Basic parameters for the trajecfoot.r
	#Foot print for EDGAR emissions

	foottimes<-seq(0,24,1)              #Want to generate a footprint for each hour
	zbot<-0				    #Set equal to 0 if we want the surface influence volume
	ztop<-0				    #Set equal to 0 if we want the surface influence volume
	numpix.x<-600                       #number of grid cells in the x direction (lon)
	numpix.y<-600                       #number of grid cells in the y direction (lat)
	lon.ll<--130.0                      #lower left corner of our footprint grid
	lat.ll<- 20.0                       #lower right corner of our footprint grid
	lon.res<-1/10                      #horizontal resolution of x grid
	lat.res<-1/10                      #horizontal resolution of y grid


	#ident is identifier flag for an "2002x08x12x18x42.54Nx072.17Wx00030"

	ident<-info["outname"]
	foot  <-Trajecfoot(ident=ident,pathname=path,foottimes=foottimes,zlim=c(zbot,ztop),fluxweighting=NULL,coarse=1,vegpath=vegpath,numpix.x=numpix.x,numpix.y=numpix.y,lon.ll=lon.ll,lat.ll=lat.ll,lon.res=lon.res,lat.res=lat.res)


	cat("Saving footprint grid to netcdf output file since foot.createTF = TRUE"); cat("\n"); cat("\n");cat("\n"); cat("\n"); cat("\n")


	#As an additional step lets save the time-integrated footprint as a netcdf file as well!

	xfoot.wfei <-apply(foot,c(1,2),sum)


	foot <- aperm(foot, c(2,1,3))
	netcdf.name = paste(timestamp,"_footprint.nc",sep="")

	foot.lat <- as.numeric(rownames(foot))
	foot.lon <- as.numeric(colnames(foot))


	x <- ncdim_def( "Lat",  "degreesN (row)", foot.lon)                 #Set equal to our lat lon vectors we created earlier
	y <- ncdim_def( "Lon",  "degreesE (column)", foot.lat)


	foot.var <- ncvar_def(name="footprint", units="PPM/(umoles/m^2 s)", list(y,x),longname="footprint")

	ncnew <- nc_create(filename=netcdf.name,vars=foot.var)

	ncvar_put(nc=ncnew,varid=foot.var,vals=xfoot.wfei)            #puts our variable into our netcdf file


	nc_close(ncnew)                                               #Closes our netcdf4 file


	#Move the output file name to our model output directory

	system(paste("mv",netcdf.name,ncdfpath))


}



#Finally lets move our model output

cat("Simulation completed, moving  model output file..."); cat("\n"); cat("\n")


#End of script!
