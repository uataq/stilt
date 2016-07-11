###### README file for stilt model R code ##################
#  $Id: 00README.TXT,v 1.19 2008-08-14 12:43:29 gerbig Exp $
############################################################

This directory contains all R scripts needed to run STILT, the receptor oriented modelling package.
Instructions to install STILT is given below under 'INSTALLING STILT'.
An example on how to run STILT is given below under 'RUNNING STILT'.
A listing of all functions used and what they do is given under 'DIRECTORY LISTING'.

######### INSTALLING STILT ########

1) If not yet done, install R, also libraries "foreign" and "fields".

2) Copy the STILT directory incuding all subdirectories (on grid.deas.harvard.edu under
   /group/wofsy/stilt/STILTdownload/) to your home directory or to any directory where you have
   write permission.

3) Compile an executable "hymodelc" from the sources and Makefile in the "stilt_hysplit" directory.
   Put the "hymodelc" executable file in a directory that is in the PATH environment variable (e.g.
   ~/bin).

4) Compile the R extensions in stiltR/shlib.
   For that, you need a compiled netCDF library with Fortran 90 bindings.
   Also make sure that your computing environment is set up properly to run shared libraries;
   probably they refer themselfs to runtime libraries of the Fortran compiler and need to find them,
   which is often accomplished by setting some environment variables (see the compiler manual).
   E.g. at Jena:
   Please source file /usr/local/etc/sources/Intel_fort.(c)sh (bash or csh) to set Fortran
   environment and libraries.
   Also set netCDF library path, e.g. export LD_LIBRARY_PATH=/Net/Groups/BSY/tpobw/STILT/IER/netCDF/netcdf-3.6.2/fortran/.libs:/Net/Groups/BSY/tpobw/STILT/IER/netCDF/netcdf-3.6.2/libsrc/.libs:$LD_LIBRARY_PATH

5) Obtain meteorological data in ARL format
   Data sources:
        - http://www.arl.noaa.gov/ss/transport/archives.html
        - some are on grid.deas.harvard.edu under /group/wofsy/stilt/Metdata/
        - some are on galactica.bgc-jena.mpg.de under
          /Net/Groups/BSY/STILT_Metdata/
        - or create from wether forecasting model output (RAMS, BRAMS, WRF, ECMWF) using e.g.
          RAMS2arl.exe (ask chg or jcl or Stefan Koerner at Jena)

6) Run the script "setup.sh". (NOTE: this does not exist as of August 08)
   This will set up working directories "Copy0", "Copy1" etc. to run the hymodelc executable within
   a main run-time directory "exe", allowing for multiple runs on different nodes. At least the
   "Copy0" directory is required. At the same place (in "exe") a directory "bdyfiles" containing
   the files ASCDATA.CFG, LANDUSE.ASC, ROUGLEN.ASC is created.
   The script sets the following variables in "setStiltparam.r":
      rundir
      shlibpath
      metpath

Modelling working directory "must have"s:
  - a "stiltR" subdirectory or link containing the R part of STILT
  - a file "setStiltparam.r" with the STILT run parameters
For job sumission it is convenient to have copies of the respective helper scripts "stilt.bsub.sh"
or "stilt.qsub.csh" here.


REST NOT UPDATED YET!
############################################
#### old 4) and 5) and 6) not rqd???

4) Start R, assign the pathname to the r script directory, and source all required scripts with sourceall.r
	Example:  sourcepath<-"~/STILT/Rsc/";source("~/STILT/Rsc/sourceall.r")

5) Only in case of mapping particles to flux grids:
	get boundary condition objects with lateral boundary condition for each desired tracer
        (time dependent latitudinal cross sections) covering desired time period
        A selection of ASCII files are on grid.deas.harvard.edu in /group/stilt/Boundary/
	Example: use function read.bg() (from script read.bg.r)
        	read.bg(spec=c("CO2","CO"),datename="1_1_99_12_31_02",pathin="/group/stilt/Boundary/",pathout="~/Rdat/")
        read.bg reads
        -- ASCII files with 4 columns: agl lat co2m sasdate,
           i.e. altitude (in m), lat (in DEG N), co2 (in ppm), day (since 1/1/1960)
	Details see in "read.bground.r" or below under 'DIRECTORY LISTING'

6) Only in case of mapping particles to flux grids AND not running on grid.deas:
        get surface flux grids (dumpfile vegflux.dmp),
        regrid into varying resolution (for dynamic aggregation)
        -- run script gen_veg_mat.r in R, before need to assign vegpath (e.,g. vegpath<-"/group/stilt/Vegetation/")
        -- quit R without saving: q(save="no")


################################


######### RUNNING STILT ########

STILT is here regarded as the whole package, that consists of the main functions
Trajecmod() ==> loop over receptor locations (e.g. points along flight track, or different times at a ground station)
Trajec() ==> calls fortran executable hymodelc and generates object with particle locations (and other variables from the metfields)
Trajecflux() ==> maps particle location objects to surface flux grids, calculates tracer mixing ratios (GSB model)
Trajecvprm() ==> maps particle location objects to surface flux grids, calculates tracer mixing ratios (VPRM model)
Trajecfoot() ==> calculates footprints or influence functions from particle location objects

A tutorial of how to use the most important functions and how to plot influence functions is in 0stilt_tutorial.r


## ## ## ## ## ## ## ## ## ##
##Preparing input and parameters
## ## ## ## ## ## ## ## ## ##
INPUT
Input to STILT are receptor locations and times. An R-object needs to be created that contains that
information as a matrix with 4 columns: fjul (fractional julian day since 1/1/1960), lat (deg N), lon (deg E) and
altitude (meters above ground). The name of this object and the location where it is stored is then used by STILT. An
example of how to create such an object is given in "create_times.r".

PARAMETERS
The script "setStiltparam.r" contains assignments for many parameters that influence STILT.
To change the parameters, setStiltparam.r has to be edited.

Basic parameters that are set generally once at the first time are: path (path where output gets saved,
and input data (Receptor locations and times) and boundary mixing ratio objects are read from);
metpath (where met data are stored in ARL format) and vegpath (path to get surface flux objects)

Typical parameters that change often are: nhrstilt (for how many hours backward), convect (T to turn on
convection), Timesname (name of object containing receptor information), overwrite (even if particle
location files are there, Trajectories will be recalculated), fluxTF (calculating mixing ratios from mapping to
fluxes), footprintTF (calculating footprints or influences).
For details see documentation within setStiltparam.r

## ## ## ## ## ## ## ## ## ##

## ## ## ## ## ## ## ## ## ##
##Running STILT
## ## ## ## ## ## ## ## ## ##

Change to directory STILT/Rsc (in which all r scripts are found).

From command line call STILT with:

node-09% stilt.bat

stilt.bat contains 1 line: 'R BATCH stilt.r stilt.log', i.e. it calls R in batch mode with the
script stilt.r, and it generates a log file called stilt.log. stilt.r sources all required functions
in the directory (refered to as sourcepath), and calls the function 'Trajecmod()'. 'Trajecmod()' calls
'Trajec()' for each starting time to calculate paricle dirstributions, and (if desired) calls
Trajecflux() to calculate mixing ratios and Trajecfoot() to derive footprints.

## ## ## ## ## ## ## ## ## ##
##Output from STILT
## ## ## ## ## ## ## ## ## ##

Particle location objects:
objects containing time since release (backward with negative time), with names that reflect
time and position of the receptor (e.g. "2002x08x16x06x42.54Nx072.17Wx00030").
These objects are saved in files like "/STILT/Output/.Rdata2002x08x16x06x42.54Nx072.17Wx00030" using assignr() (see below),
and they can be retrieved using getr (see below).

Logfile:
--info for individual receptor runs:
saves logfile to object with date and time in name
example: ./Runs.done/.RDatarun.info.Mar..9.18:36:21.2004
this can be retrieved to R with info<-getr("run.info.Mar..9.18:36:21.2004",path="./Runs.done/")
--stilt.log contains all printed output from stilt.r

Maps/Matrices of Footprints/Influence functions:
objects containing Footprints/Influence, with names that reflect
time and position of the receptor (e.g. "foot2002x08x16x06x42.54Nx072.17Wx00030").
These objects are saved in files like "/STILT/Output/.Rdatafoot2002x08x16x06x42.54Nx072.17Wx00030" using assignr() (see below),
and they can be retrieved using getr (see below).


################################


######### DIRECTORY LISTING ########

###################################
assignr.r
#function assignr(xname, value, path="",printTF=FALSE,gz=FALSE)
#see also getr() and existsr()
#to avoid large databases that take too long to load memory (as R tries...)
#similar to assign(), but object is removed after call to assignr
#assigns object to name xname, saves it with save() under path/.RDataxname
#and REMOVES it from local database (search()[1])
#example: assignr("test",temp,path="/mydir/") creates file "/mydir/.RDatatest"
#that can be attached and contains a single object with name "test"
#printTF:  if TRUE, prints info about file creation
#gz:       if TRUE, compresses file after creation with gzip
###################################

###################################
 combine.met.r
#function combine.met(year=99,startmonth=1,endmonth=12, metd="fnl",metlib="/group/stilt/Metdata/",outlib=metlib)
#function to combine arl met data (to avoid hymodelc crashing at beginning of month...)
#combines met files, using cat
#'year':        2 digit year
#'startmonth':  first month of selected period
#'endmonth':    last month of selected period
#'metd':        type of metfile ("edas" or "fnl")
#'metlib':      where metdata are read from
#'outlib':      where metdata are written to
###################################

###################################
create_times.r
#script
#creates object "Times.hf" w/ starting times for HF
#output is a matrix with 4 columns:
#-fjul (fractional julian day since 1/1/1960)
#-lat (deg N)
#-lon (deg E)
#-altitude (meters above ground)
###################################

###################################
day.of.week.r
#function day.of.week(month, day, year)
#day.of.week(month, day, year) returns day of week as number (0: Sun, 6: Sat)
###################################

###################################
existsr.r
#function existsr(xname, path="")
#see also assignr() and getr()
#similar to exists(), but
#checks for object on stored location (file path/.Rdataxname)
###################################

####################################
flttrack.r
#function(lon,lat,col,cex,cities,pts,newmotifTF)
#Plots the flight tracks on top of map
####################################

###################################
gen_veg_mat.r
#script
#generates matrices with different resolutions for vegetation maps, from 1/6lat*1/4lon to coarser resolution
#saves objects with assignr() in location defined by 'vegpath' (set in script setStiltparam.r)
###################################

###################################
getgridp.r
#function getgridp(min.x, max.x, min.y, max.y, numpix.x, numpix.y,coarse.factor=1)
#Implements decisions in trajectory model about which emission grid to use
#  by taking in vectors of min & max grid points
#Each element is a time point
#Returns vector of filenames (char) that would give correct file for emission grid
###################################

###################################
getmetfile.r
#function getmetfile(yr=0,mon,day,hr,nhrs,metd="edas",doublefiles=F)
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
#'metd' meteorological data identifier, such as "edas", or "fnl" or "brams" or "rams" (RAMS NOT YET!)
#'doublefiles' should concatenated met files be used? Allows starting times between files
###################################

###################################
getr.r
#function getr(xname, path="",gz=FALSE)
#see also assignr() and existsr()
#to avoid large databases that take too long to load memory (as R tries...)
#similar to get()
#gets object from stored location (file path/.Rdataxname)
#name of object is xname if file was writen with assignr()
#gz: flag, if TRUE, file is first gunziped (if uncompressed version does not exist), then uncompressed version is deleted.
###################################

###################################
id2pos.r
#function id2pos(id,sep="x")
#revers of id2pos
#function to create read identifying label for single receptor (location&time)
#returns time as fractional julian day since 1/1/1960
#and alt as altitude above ground in kilometers
#example:
# id2pos("2002x08x03x10x+45.00x+090.00x00030")
#[1] 15555.42    45.00    90.00    30.00
###################################

###################################
imagell.r
#function imagell(mat,lg=F,zlims=NULL,zoom=0,center=c(-72.172,42.536),rays=F){
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
#'lon.res'                      #resolution in degrees longitude
#'lat.res'                      #resolution in degrees latitude
#'nlevel' Number of color levels used in legend strip
#'col' Color table to use for image ( see help file on image for details). Default is a pleasing range of 64 divisions on a topographic scale.
#'citiesTF' if T, plot locations of cities
#'part' a particle location object can be passed on, it will be plotted on top of footprints
#'dev' device to which graphic is send (currently X11 or file name with .png ending, creates png file)
#'type' determines label for color legend
#       "foot" is footprint, uses ppm/micro-moles/m2/s
#       "infl" is influence, uses 1/m^2 (influence density integrated over altitude)
#       any other: used as label for color legend
#
#returns vector c(zmin,zmax) for scaling
###################################

###################################
image.plot.fix.r
#function image.plot.fix()
#fix to image.plot()
###################################

###################################
image.plot.plt.fix.r
#function image.plot.plt.fix()
#fix to image.plot.plt()
###################################

###################################
is.inf.r
#function is.inf(x)
#similar to Splus function; returns T if x is infinite
###################################

###################################
julian.r
#function julian(m, d, y, origin.)
#from Splus
#returns day since 1/1/1960
#'origin.' has format c(month = 1, day = 1, year = 1960)
###################################

###################################
leap.year.r
#function leap.year(y)
#from Splus
###################################

###################################
lsr.r
function lsr(path)
#lists all r objects in 'path' stored with assignr()
###################################

###################################
memory.size.r
#function memory.size()
#instead of memory.size in splus
#for linux only
#returns string indicating memory usage
###################################

###################################
month.day.year.r
#function month.day.year(jul, origin.)
#reverse of julian
#returns month, day, year
#'jul' julian since 'origin.'
#'origin.' has format c(month = 1, day = 1, year = 1960)
###################################

###################################
motif.r
#function(...)
#For backward compatibility with older S scripts
###################################

###################################
pos2id.r
#function pos2id(jultime,lat,lon,alt,sep="x")
#function to create identifying label for single receptor (location&time)
#expects jultime as fractional julian day since 1/1/1960
#expects alt as altitude above ground in meters
#example:
# pos2id(15555.4166667,45,-90,0.03)
#[1] "2002x08x03x10x+45.00x+090.00x00030"
###################################

###################################
read.asc.r
#function read.asc(file)
#Facilitates importing files generated by SAS, which are blank-separated values
#  and may contain '.' as symbol for NAs
###################################

###################################
read.bground.r
#function read.bg(spec,datename,pathin,pathout)
#reads initial field for CO2, CO, CH4, or H2, from ascii data that were created using sas
#
#spec is any subset of c("co","co2","ch4","h2")
#datename is a string indicating time range, e.g. "1_1_99_12_31_02" for 1/1/1999 to 12/31/2002
#pathin is path where to read ASCII data from
#pathout is path where to assign boundary condition objects
#
#OUTPUT
#assigns *.ini objects (3-D array altitude * Latitude * Day since 1/1/1960)
#time is 0 at 1/1/1960 (i.e. elapsed days since 1/1/1960; similar to julian() with default origin.)
# e.g. read.bg(spec=c("CO2","CO"),datename="1_1_99_12_31_02",pathin="/group/stilt/Boundary/",pathout="~johnlin/Rdat/")

###################################

###################################
rmr.r
#function rmr(xname,path)
#similar to exists(), but
#removes object at stored location (file path/.Rdataxname)
###################################

###################################
rp2ll.r
#function rp2ll(r,p,lon0=-72.172,lat0=42.536)
#function to convert from polar coordinates (r, p as distance and angle) to lon,lat
#'r'	distance from center
#'p'	angle (north=0, east=90)
#'lon0' center longitude (origin)
#'lat0'	center latitude (origin)
#
#returns lat and lon
###################################

###################################
setStiltparam.r
#script
#set parameters needed for STILT
##THIS NEEDS TO BE EDITED SO THAT STILT DOES THE RIGHT THINGS
###################################

###################################
sourceall.r
#script
#sources all required R functions found in directory 'sourcepath' (global variable)
###################################

###################################
stdev.r
#function stdev(x,na.rm=F)
#returns standard deviation (sqrt(var))
###################################

###################################
stilt.r
#script
#sources all required R functions for STILT, then calls Trajecmod()
#saves logfile (info for individual receptor runs) to object with date and time in name
#example: ./Runs.done/.RDatarun.info.Mar..9.18:36:21.2004
#this can be retrieved to R with info<-getr("run.info.Mar..9.18:36:21.2004",path="./Runs.done/")
###################################

###################################
trajwind.r
#function (yr,mon,day,hr,lat,lon,agl,metlib=paste(unix("echo $HOME"),"/Metdat/",sep=""),
#                      metfile=NULL,rundir="~/STILT/Exe/",RAMSTF=F,nummodel=0){
#Uses differences in positions of mean trajectories to calculate the mean wind at a specified time & position
###################################

###################################
Trajecflux.r
#function Trajecflux(ident,pathname="",tracers=c("CO","CO2"),coarse=1,dmassTF=T,nhrs=NULL,vegpath="/group/stilt/Vegetation/")
#calculate tracer concentrations for trajectories
#maps trajectories onto surface fluxes
###
#'ident' is character value specifying the trajectory ensemble to look at
#'pathname' is path where object with particle locations is saved
#'tracers' vector of names for which mixing ratios are wanted; any subset of c("co","co2","ch4","h2")
#'coarse' degrade resolution (for aggregation error): 0: only 20 km resolution;
#        1-16: dynamic resolution, but limited highest resolution
#	coarse:	  (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
#	coarsex:c(1,1,2,2,2,4,4,4,8, 8, 8,16,16,16,32,32) #factors by which grids have been made coarser
#	coarsey:c(1,2,1,2,4,2,4,8,4, 8,16, 8,16,32,16,32)#e.g., '4' means grid is 4 times coarser
#'dmassTF' if TRUE, weighting by accumulated mass due to violation of mass conservation in met fields
#'vegpath' is path under wich vegetation and flux grids are stored
###
#returns vector containing mixing ratios etc. at endpoint (receptor)
#for timeseries (e.g. mixing ratio as fct of time back) needs some modifications
###################################

###################################
Trajecfoot.r
#function Trajecfoot(ident,pathname="",foottimes,zlim=c(0,0),fluxweighting=NULL,coarse=1,dmassTF=T,vegpath="/home/gerbig/Rdat/Vegetation/",
		     numpix.x=376,numpix.y=324,lon.ll=-145,lat.ll=11,lon.res=1/4,lat.res=1/6)
#Creates footprint or influence for individual particle runs
#'ident' is character value specifying the trajectory ensemble to look at
#'pathname' is path where object with particle locations is saved
#'foottimes' is vector of times between which footprint or influence will be integrated
#'zlim' (if not default 0,0): vertical interval for which influence is calculated
#'fluxweighting' if not NULL, but a number (1-31), weighting by that vegetation class or emission (e.g. 19 for CO fossil fuel emissions)
#'coarse' degrade resolution (for aggregation error): 0: only 20 km resolution;
#        1-16: dynamic resolution, but limited highest resolution
#	coarse:	  (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
#	coarsex:c(1,1,2,2,2,4,4,4,8, 8, 8,16,16,16,32,32) #factors by which grids have been made coarser
#	coarsey:c(1,2,1,2,4,2,4,8,4, 8,16, 8,16,32,16,32)#e.g., '4' means grid is 4 times coarser
#'dmassTF' if TRUE, weighting by accumulated mass due to violation of mass conservation in met fields
#'vegpath' is path under wich vegetation and flux grids are stored
#'numpix.x'			#number of pixels in x directions in grid
#'numpix.y' 			#number of pixels in y directions in grid
#'lon.ll'			#lower left corner of grid
#'lat.ll'			#lower left corner of grid
#'lon.res'			#resolution in degrees longitude
#'lat.res'			#resolution in degrees latitude
###
#returns: 3-D array (lat, lon, time) containing influence grids or footprints (lat,lon) for different times back
###################################

###################################
Trajecmod.r
#function Trajecmod()
#no arguments; these are set in setStiltparam.r
#Function that loops over all receptors (or only a part of it)
#Calls 'Trajec()' for each starting time to calculate paricle dirstributions
##
#depending on flags set in set in setStiltparam.r:
#calls Trajecflux(), and assigns mixing ratio results for all receptors (e.g. 'stiltresult1' for first part)
#calls Trajecfoot() to derive footprints
##
#returns: object containing info generated by each run of Trajec()
###################################

###################################
Trajec.r
#function Trajec(yr=02,mon=8,day=1,hr=6,lat=42.536,lon=-72.172,agl=30,nhrs=-48,maxdist=20,
           delt=0.0,numpar=100,ndump=0,random=T,outdt=0.0,veght=0.5,metlib="/group/stilt/Metdata/",
           metd="edas",doublefiles=F,metfile=NULL,nturb=F,outfrac=0.9,conv=F,ziscale=NULL,
           siguverr=NULL,TLuverr=NULL,zcoruverr=NULL,horcoruverr=NULL,
           varsout=c("time","index","lat","lon","agl","grdht","foot","temp0","swrad","zi","dens","dmass"),
	   rundir=NULL,nummodel=0,outname=NULL,outpath="",overwrite=T){
#
#to run HYSPLIT particle dispersion model and to check distribution of particles
##
#INPUT:
#'yr','mon','day','hr': starting time
#'lat','lon',&'agl' can be a VECTOR of the same length--to have multiple starting locations; agl in meters above ground
#'nhrs' is number of hours model would be run--NEGATIVE values mean model is run BACKWARDS
#'maxdist' maximum distance (in km) that particles travel between timesteps.
#          assumes a default gridcell size of 80 km (edas: 80 km; fnl: 180 km; rams: 45 km depending on 'metd')
#'delt' is timestep in minutes (integer...); if 0 then dynamic timestep (depending on 'trat')
#'numpar' is number of particles emitted over the # of hrs specified in 'emisshrs'
#'ndump' to dump out all particle/puff points to a file PARDUMP that can be read at start of new simulation to continue prev calc.
#   valid NDUMP settings: 0 - no I/O, 1- read and write, 2 - read only, 3 - write only. Default value = 0
#'random' is flag that means that random number generator would generate diff random number sequence each time model is run
#'outdt' is the interval [min] that elapses before particle results are written out to PARTICLE.DAT
#     if outdt=0.0, then data at EVERY timestep is written out; outdt should be a POSITIVE number
#'veght' is height in meters above ground (model ground) below which time is counted as particle seeing the ground
#     if <1 then interpreted as fraction of zi (mixed layer height as derived from met data)
#'metd' is character vector with names (descriptors) of met files to be used; possible entries: "edas","fnl","brams" or "rams" (not yet)
#'doublefiles' should concatenated met files be used? Allows starting times between files
#       concatenation with "cat file1 file2 > file12"
#'metfile' specifies the meteorological input file;if not specified, then let 'getmetfile' automatically determine filename based on time and 'metd'
#'nturb' is NotTURBulence flag that turns turbulence on (FALSE) or off (TRUE)
#'outfrac' is fraction of particles which are allowed to leave the model area before hysplit stops
#'conv' turns on convection (RAMS winds: grell convection scheme, EDAS and FNL: simple excessive redistribution within vertical range with CAPE>0)
#'ziscale' is a vector with which to scale the modelled mixed-layer height
#    each element specifies scaling factor for each model simulation hour (ziscale can be of length that is smaller than abs(nhrs))
#'siguverr' & 'TLuverr' refer to the stddev of magnitude in horizontal wind errors [m/s] and their correlation timescale [min]
#   'zcoruverr' refers to the vertical correlation lengthscale of horizontal winds [m]
#   'horcoruverr' refers to the horizontal correlation lengthscale of horizontal winds [km]
#'varsout' specifies output variables from STILT
#      can be any subset of c("time","sigmaw","TL","lon","lat","agl","grdht","index","cldidx","temp0","sampt","foot","shtf","lhtf","tcld","dmass","dens","rhf","sphu","solw","lcld","zloc","swrad","wbar","zi","totrain","convrain")
#'nummodel' specifies copy of directory where fortran executable is executed; needs to be different for different runs running parallel on same filesystem
#'rundir' specifies main directory where differend copy directories are found (see nummodel)
#'outname' specifies name of the object for output; if not specified, uses default name
#     based on time and position using pos2id() (e.g. "2002x08x16x06x42.54Nx072.17Wx00030")
#'outpath' specifies the directory in which the object will be saved
#'overwrite' if TRUE (default), overwrite existing object with same 'outname' or same default name
##
#OUTPUT:
#assigns the output of particle dispersion model in MATRIX format to an object called 'outname' (or default name indicating time & position);
#object is saved in database at location depending on outpath
#e.g. for outname="tmp" and outpath="/home/gerbig/modeloutput/" the database will be saved as "/home/gerbig/modeloutput/.RDatatmp"
#and contains the object "tmp" that can be retrieved with getr("tmp",path="/home/gerbig/modeloutput/")
#columns are specified using 'varsout' argument
#returns list of:
#  defaultname; all input data; metd with times when switched;
#-status
#  1: new object assigned, no problem;
#  2: new object assigned, ended early
#  3: object already exists, not overwriten
#  4: no object assigned; failed
###################################

###################################
unix.r
#function unix(command, intern = TRUE, ignore.stderr = FALSE)
#redefines unix as system() with intern=T as default
###################################

###################################
unix.shell.r
#function unix.shell(command, shell = "/bin/sh", ...){
#like in Splus
###################################

###################################
weekdayhr.r
#function weekdayhr(yr,mon,day,hr,runtt,diffGMT=NA)
#Determines day of week and hour of day (needed to determine emission factors used to scale emissions)
#'yr','mon','day','hr' are individual numbers used to specify the starting time
#		let 'hr' run from 0~23
#'runtt' is the runtime from starting time in MINUTES; 'runtt' can be both pos or neg, depending on forward or backward run
#'diffGMT' can be a vector with the same length as 'runtt'--the DIFFERENCE from GMT at each timestep, which can be dependent on
#	the longitude of current timestep--this is important if need this function to return LOCAL TIME
#Returns MATRIX with following columns:   1) yr 2) mon  3) day  4) hour of day  5)day of week
#	Day of week is value of 0~6, denoting Sun~Sat
###################################


