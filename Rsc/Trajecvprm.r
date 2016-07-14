#***************************************************************************************************
# calculate tracer concentrations in trajectories by mapping trajectories onto surface fluxes and
# adding CO2 boundary conditions
#***************************************************************************************************

Trajecvprm <- function(ident, pathname="", tracers=c("CO", "CO2"), coarse=1, dmassTF=T,
                       nhrs=NULL,
                       vegpath="/home/dmatross/ModelLab/VPRM/DevanVegRevise/",
                       evilswipath="/deas/group/cobra/DMMvprmTest/EVILSWI2004FULL/",
                       vprmconstantspath="/Net/Groups/BSY/BSY_3/cgerbig/RData/CarboEurope/",vprmconstantsname="vprmConstants.optCE",
                       nldaspath="/deas/group/cobra/DMMvprmTest/Radiation/NLDAS/",
                       nldasrad=FALSE, nldastemp=FALSE, pre2004=FALSE, keepevimaps=FALSE,
                       detailsTF=FALSE, linveg=FALSE,
                       numpix.x=376, numpix.y=324, lon.ll=-145, lat.ll=11,
                       lon.res=1/4, lat.res=1/6, bios="VPRM", landcov="SYNMAP.VPRM8") {

# --------------------------------------------------------------------------------------------------
# Interface
# =========
#
# ident         character value specifying the trajectory ensemble to look at
# pathname      path where object with particle locations is saved
# tracers       vector of names for which mixing ratios are wanted; any subset of
#               c("co", "co2", "ch4", "h2", "n2o")
# coarse        degrade resolution (for aggregation error): 0: only 20 km resolution;
#               1-16: dynamic resolution, but limited highest resolution
#               coarse:         (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
#               coarsex:c(1,1,2,2,2,4,4,4,8, 8, 8,16,16,16,32,32) # factors by which grids have been made coarser
#               coarsey:c(1,2,1,2,4,2,4,8,4, 8,16, 8,16,32,16,32)# e.g., 4 means grid is 4 times coarser
# dmassTF       if TRUE, weighting by accumulated mass due to violation of mass conservation in met fields
# vegpath       path under which vegetation and flux grids are stored
# evilswipath   path under which evi and lswi grids are stored
# vprmconstantspath     path under which vprm constants are stored (e.g. "/Net/Groups/BSY/BSY_3/cgerbig/RData/CarboEurope/")
# vprmconstantsname     name under which vprm constants are stored (e.g. "vprmconstants.12" or "vprmconstants.optCE" for EU) 
# nldaspath     path where gridded NLDAS temperature and radiation are stored
# nldasrad, nldastemp   use nldas inventory radiation and temperature instead of assimilated meteorology
# pre2004       is the year prior to 2004 or 2004+
# keepevimaps   this assigns evi and lswi maps to the global environment -- useful if this function will be called repeatedly
#               BE CAREFUL when using this, as it will put a number of LARGE objects in your database. If dynamic memory
#               allocation problems are anticipated, this feature should NOT be used!!!
# detailsTF     if TRUE, for each particle the flux contribution is saved in a big object (same name as ident, with "result" attached)
# linveg        if TRUE, CO2 fluxes are represented as linear functions of temperature and radiation (otherwise GEE is non linear)
# numpix.x      number of pixels in x directions in grid
# numpix.y      number of pixels in y directions in grid
# lon.ll        lower left corner of grid (longitude of southwest corner of southwest corner gridcell)
# lat.ll        lower left corner of grid (latitude of southwest corner of southwest corner gridcell)
# lon.res       resolution in degrees longitude
# lat.res       resolution in degrees latitude
# bios          which biosphere to use?
#               GSB -- Gerbig stupid biosphere
#               VPRM -- Devans biosphere
# landcover     which landcover to use?
#               IGBP -- Original IGBP, from Christoph
#               GLCC -- Original VPRM, based on updated IGBP (GLC2.0)
#               DVN -- The updated Devan system
#               SYNMAP Martin Jung synmap product
#               SYNMAP.VPRM8 8 VPRM classes based on synmap product
#               SYNMAP.VPRM9 9 VPRM classes based on synmap product

#global variables (see setStiltparam.r): 
# inikind    trajectory initial values. Possible: "climat", "CT" (CarbonTracker), "TM3". Vector, specifies info for each desired tracer
# inifile    absolute file name. Vector, specifies info for each desired tracer
#
#
# ---------------------------------------------------------------------------------------------------
#
# Algorithm
# =========
#
#    Step 1: Initialize & Set up output
#    Step 2: Retrieve trajectory information
#    Step 3: Set-up Grid information
#    Step 4: Accounting--mass correction and organizing particle info into only what is relevant to
#            our ground
#    Step 5: OH effects on CO
#    Step 6: Further reduce particle information to surface influence
#    Step 7: Define Emission grid wrt dynamic gridding
#    Step 8: Setup vegetation grid loop
#    Step 10: Adjust Fossil Fuel Emissions
#    Step 11: Reclassify vegetation, calculate a priori vegetation flux
#    Step 12: add CO2 boundary values
#
#
# ---------------------------------------------------------------------------------------------------
#
# History
# =======
#
# This is an attempt a rebuild Trajecflux to include the VPRM (Vegetation and Photosynthesis Model)
# into the STILT code. It represents a major departure from the GSB (Gerbig's Stupid Biosphere) originally written into
# the Trajecflux code.
#
# 1) We will now use the Xiao VPM model (Xiao et al. 2004) to get an a priori flux estimate.
#
# NEE = GEE + R
# GEE ~= GPP
# GPP = eta * FAPAR.pav * PAR
#      eta=eta.0*T.sc * W.sc * P.sc * PAR.sc, where eta.0 = 0.044 initially (according to Xiao)
#            T.sc = [(T-Tmin)(T-Tmax)]/[((T-Tmin)(T-Tmax)) - (T - Topt)^2]
#            W.sc = 1+LSWI/1+LSWImax
#           P.sc = 1+LSWI/2 {Bud Burst period} -OR- P.sc = 1 {rest of year}
#            PAR.sc = 1/(1 + PAR/PAR.0)
#      FAPAR.pav = a*EVI where a~=1 initially 
#
#     this translates to:
# GPP = lambda.par * T.sc * P.sc * W.sc * PAR.sc * EVI * PAR
#
# R = alpha*T + beta
#   
# Given an increase in the number of vegetation classes {up to 11} there will be an individual lambda parameter,
# where lambda is a scaling factor. EVI will vary over time, and we will assume a functional form with fitted 
# parameters, based on Julian time. The maximum values and minimum values will be derived from the table
# which comes from that functional form.
#
#
# 2) New capability has been added to make it possible to use Devan's vegetation classes (moniker= "DVN").
# With this, there are 11 vegetation classes, 4 evergreen and 7 others. The new vegetation scheme is as follows:
#
#                 VPRM Class    STILT-VPRM class   
# Evergreen A                1A                1
# Evergreen B                1B                2
# Evergreen C                1C                3
# Evergreen D                1D                4
# Deciduous                   2                5
# Mixed forest                3                6
# Shrubland                   4                7
# Savanna                     5                8
# Cropland-Soy               6A                9
# Cropland-Corn              6B                9
# Grassland                   7               10
# Peat                        8               11
# Others                      9               12
#
# 3) Potential global variable capability for repeated runs, eliminating memory calls
#
# 4) Now include NLDAS temperature as well as NLDAS radiation. The NLDAS radiation fields and temperature
# fields remain our best option for inputs. The STILT-based temp0 and swrad can also be used, depending
# on met input, but NLDAS represents a reanalysis version.
#
# 5) Added 'pre2004' variable. For 2003 and earlier, our radiation conversion factors are different. This is mostly
# used in converting shortwave radiation to par.
#
#---------------------------------------------------------------------------------------------------
# $Id: Trajecvprm.r,v 1.26 2009-11-25 08:51:03 gerbig Exp $
#---------------------------------------------------------------------------------------------------


####################################################################################################
# Step 1: Initialize & Set up output
####################################################################################################
#for each tracer have fields for fluxes and for LBC specified in setStiltparam.r
tracers <- tolower(tracers)                                 # use lower case

#get names for lateral boundary condition tracers for each tracer
tracersini <- c(paste(tracers, "ini", sep=""))
if ("co"%in%tracers) tracersini <- c(tracersini, "coinio")  # also want initial CO without OH...

#get names for flux tracer for each tracer ("ffm" appended from legacy code, "fossil fuel model")
tracersffm <- c(paste(tracers, "ffm", sep=""))
#if ("co2"%in%tracers&bios!="") #add tracer release over water when have landuse available for biosphere
#   tracersffm <- c(tracersffm, "inflwater")

# get names for output
names.output <- c("ident", "latstart", "lonstart", "aglstart", "btime", "late", "lone", "agle",
                  "zi", "grdht", "nendinarea", tracersini, paste("sd", tracersini, sep=""), tracersffm)

#add names for tracers associated with biospheric fluxes: GEE (photosynthesis), respiration, and influence for each vegetation type
out.type <- NULL
output.veg <- NULL
if ("co2"%in%tracers)  {
   if (bios != "GSB" & bios != "VPRM" & bios != "")
      stop("Trajecvprm: parameter bios incorrectly specified. Redefine as GSB or VPRM or quote-unquote.")
   if (bios == "GSB" | bios == "VPRM"){
      out.type <- c("infl", "gee", "resp")
      cat(bios, "\n")
      if (bios == "GSB") {
         output.veg <- c("frst", "shrb", "crop", "wetl")
      }
      if (bios != "GSB"&(landcov == "GLCC"|landcov == "SYNMAP"|landcov == "SYNMAP.VPRM8")){
        output.veg <- c("evergreen", "decid", "mixfrst", "shrb", "savan", "crop", "grass", "peat") #"peat" is replaced by "others" in Jena VPRM preproc.
        if(substring(evilswipath,1,21)=="/Net/Groups/BSY/data/")output.veg <- c("evergreen", "decid", "mixfrst", "shrb", "savan", "crop", "grass", "others") #"peat" is replaced by "others" in Jena VPRM preproc.
      }
      if (bios != "GSB"&(landcov == "DVN"))
         output.veg <- c("evergreenA", "evergreenB", "evergreenC", "evergreenD", "decid", "mixfrst", "shrb", "savan", "crop", "grass", "peat")
      names.output <- c(names.output, as.vector(outer(out.type, output.veg, paste, sep="")))
   }
}

####################################################################################################
# Step 2: Retrieve trajectory information
####################################################################################################

# Check if object exists
if (existsr(ident,pathname)) {
   cat("Trajecvprm(): starting with ident=", ident, "\n")   #found object
   part <- getr(ident,pathname)                             #get it
} else {
   if (paste(pathname,".RData",ident,".gz",sep="") %in%
      dir(pathname,pattern=paste(".RData",ident,sep=""),all.files=TRUE)) {
      cat("Trajecvprm(): starting with .gz ident=",ident, "\n") #found .gz object
      part <- getr(ident,pathname,gz=TRUE)                  #get it
   } else {
    cat("object ", pathname,".RData", ident, " NOT FOUND\n", sep="")
    lastresult <- rep(NA,length(names.output)-1)
    lastresult <- c(ident,lastresult)
    names(lastresult) <- names.output
    return(lastresult)
   }
} #if exists or not


# get time and position information from name (ident)
pos <- id2pos(ident)
time <- month.day.year(floor(pos[1]))
yr4 <- time$year                                            # 4 digit year
yr <- yr4%%100                                              # 2 digit year (or 1 digit...)
mon <- time$month
day <- time$day
hr <- round((pos[1]-floor(pos[1]))*24)

# check if required particle information is available
rqdnames <- c("time", "lat", "lon", "agl", "zi", "index", "temp0", "foot")
if ("co2"%in%tracers) {
   rqdnames <- c(rqdnames, "swrad")
   # need pressure for hybrid coordinates
   if (inikind["co2"] == "CT" || inikind["co2"] == "TM3") rqdnames <- c(rqdnames, "pres")
}
if (dmassTF)
        rqdnames <- c(rqdnames, "dmass")
for (nm in rqdnames) {
   if (!(nm%in%dimnames(part)[[2]])) {
      cat("need column '", nm, "' for this run\n", sep="")
      lastresult <- rep(NA, length(names.output)-1)
      lastresult <- c(ident, lastresult)
      names(lastresult) <- names.output
      return(lastresult)
   }
}


####################################################################################################
# Step 3: Set-up Grid information
####################################################################################################

part[, "time"] <- (-part[, "time"]/60) # use time back from now on, and transform to hours
dimnames(part)[[2]][dimnames(part)[[2]] == "time"] <- "btime"
part<-part[order(part[, "btime"], part[, "index"]),] #order by btime, then by index
# get grid indices
# For horizontal grids (lower left corner of south-west gridcell: 11N,145W; resolution: 1/4 lon, 1/6 lat, 376 (x) times 324 (y))
gitx <- floor(1/lon.res*(part[, "lon"]-lon.ll)+1)
gity <- floor(1/lat.res*(part[, "lat"]-lat.ll)+1)
part <- cbind(part, gitx, gity)
dimnames(part) <- list(NULL, dimnames(part)[[2]])



####################################################################################################
# Step 4: Accounting--mass correction and organizing particle info into only what is relevant to our
# ground area, as defined in step 3.
####################################################################################################

if (dmassTF) {
  #remove particles with too strong dmass violation
  ind<-unique(part[part[,"dmass"]>1E3|part[,"dmass"]<1/1E3,"index"])
  if (length(ind) >= length(unique(part[, "index"]))/2){
    message("Trajecvprm(): ", length(ind), ' of ', length(unique(part[, "index"])), ' particles have mass defect; returning NA')
    lastresult <- rep(NA, length(names.output)-1)
    lastresult <- c(ident, lastresult)
    names(lastresult) <- names.output
    return(lastresult)
  }
  part<-part[!part[,"index"]%in%ind,]

  # get average dmass to "correct correction" (allow multiplication w/ dmass without changing total mass)
  # i.e. correction for average mass loss of particles, since they get attracted to areas of mass destruction
  mean.dmass <- tapply(part[, "dmass"], part[, "btime"], mean) # this gives for each btime a mean dmass
  # DMM
  # To account for situations where mean.dmass is zero (mass violation total), need to avoid division by zero to
  # avoid downstream problems.
  mean.dmass[which(mean.dmass == 0)] <- 0.00001

  # need to "merge" this with part; can't use array since not all paticles are left at very large btime
  nparleft <- rle(part[, "btime"])$length # number of times the same btime is repeated
  mean.dmass <- rep(mean.dmass, nparleft) # long vector of mean dmass
  # need to link this info to each particles dmass: normalize individual dmass by mean.dmass
  part[, "dmass"] <- part[, "dmass"]/mean.dmass             # Dan Matross gets problems with that
}  
# To allow short runs, cut off after nhrs
if (!is.null(nhrs))
  part <- part[abs(part[, "btime"])<=abs(nhrs), ]

# remove particles when they cross the longitude -145 West for the first time (that's where the climatology is valid)
inbgarea <- floor(part[, "gitx"])<1
part[inbgarea, c("gitx", "gity")] <- NA
dimnames(part) <- list(NULL, dimnames(part)[[2]])

# remove points when they enter the background area for the first time
sumx <- tapply(part[, "gitx"], part[, "index"], cumsum) # cumsum gives "NA" after first occurence of "NA"
ordert <- order(part[, "index"], part[, "btime"]) # order first by index, then by time
ordern <- order(part[ordert, "btime"], part[ordert, "index"]) # to recreate original order
sumx <- unlist(sumx)[ordern]
part <- part[!is.na(sumx), ]
dimnames(part) <- list(NULL, dimnames(part)[[2]])

# only keep points, when information is changed:
# 1. position and times at boundary of NGM grid
# 2. position and times when surface influences particles
# 3. position and times at start of trajectory

# 1. boundary
ordert <- order(part[, "index"], part[, "btime"]) # order first by index, then by time
ordern <- order(part[ordert, "btime"], part[ordert, "index"]) # to recreate original order
delbte <- c(diff(part[ordert, "btime"]), -1000)[ordern] # timestep will be negative at last obs. for each particle
selend <- delbte<0

# 2. surface influence over domain area
in.domain <- floor(part[, "gitx"])<=numpix.x&floor(part[, "gitx"])>=1&floor(part[, "gity"])<=numpix.y&floor(part[, "gity"])>=1
selinf <- part[, "foot"]>0&in.domain

# 3. start of particle trajectory
delbte <- c(-1000, diff(part[ordert, "btime"]))[ordern] # timestep will be negative at first obs. for each particle
selfirst <- delbte<0



####################################################################################################
# Step 5: OH effects on CO
####################################################################################################

if ("co"%in%tracers) {
   #################### CO + OH losses########################################
   # first get OH at particle position, calculate COrel. loss and CH4 source
   # OH from SAS, parameterized as oh=oh0+oh1*p+oh2*p**2, with parameters ohi for each month and for 30 and 60 lat
   pmb <- 1013*exp(-part[, "agl"]/8000) # 8 km scale height
   ohm <- oh[oh[, "month"] == mon, ]
   oh60 <- ohm[ohm[, "lat"] == 60, "oh0"]+ohm[ohm[, "lat"] == 60, "oh1"]*pmb+ohm[ohm[, "lat"] == 60, "oh2"]*pmb**2 # OH at 60 north, same altitude
   oh30 <- ohm[ohm[, "lat"] == 30, "oh0"]+ohm[ohm[, "lat"] == 30, "oh1"]*pmb+ohm[ohm[, "lat"] == 30, "oh2"]*pmb**2 # OH at 30 north, same altitude
   oh60[oh60<0] <- 0; oh30[oh30<0] <- 0
   # local OH, interpolated
   ohl <- (oh60*(part[, "lat"]-30)+oh30*(60-part[, "lat"]))/(60-30)
   ohl[ohl<0] <- 0
   tair.k <- part[, "temp0"]-part[, "agl"]*6.5/1000 # average lapse rate
   # get CO loss
   kohCO <- 1.3E-13*(1+(0.6*pmb/1000))*300/tair.k # in 1/cm^3/sec
   tau <- 1/(kohCO*ohl)/3600/24; # in days
   # integral (k*OH*dt)
   # get dbtime in same format as kohCO and ohl
   delbt <- c(0, diff(part[ordert, "btime"]))[ordern][!selfirst]*60*60 # timestep in seconds
   kohdt <- kohCO[!selfirst]*ohl[!selfirst]*delbt
   # integral over dt: sum over timesteps for each particle individually
   Ikohdt <- tapply(kohdt, part[!selfirst, "index"], sum)
   CO.frac <- exp(-Ikohdt) # preliminary result: factors for each particles initial CO
   CO.fact <- part[, "btime"]*0 # initialize
   CO.fact[selend] <- CO.frac

   #################### CO from CH4 + OH ######################################
   kohCH4 <- 2.3E-12*exp(-1765/tair.k);
   kohCH4dt <- kohCH4[!selfirst]*ohl[!selfirst]*delbt
   IkohCH4dt <- tapply(kohCH4dt, part[!selfirst, "index"], sum)
   COfrCH4 <- 1.780*IkohCH4dt # CO from CH4 in ppm
   # not all will make it, account for some losses (about half the losses for COini)
   COfrCH4 <- COfrCH4*exp(-0.5*Ikohdt)
   CO.frCH4 <- part[, "btime"]*0 # initialize
   CO.frCH4[selend] <- COfrCH4
   part <- cbind(part, CO.fact, CO.frCH4) # add chem loss factor for CO and production from CH4 to particle location object
}



####################################################################################################
# Step 6: Further reduce particle information to surface influence
####################################################################################################

# keep only particles where they matter: create flag for first, last, or when surface influence
# apply selection here (only first, last, or when surface influence)
part <- part[(selinf+selend+selfirst)>0, ]
# print("3"); print(dim(part)); print(sum(in.domain)); print(sum(part[, "foot"]>0))
# move x and y position of final position to initialization area (gitx=1, gity= 1 to numpix.y), at least make sure they are not outside NGM
part[part[, "gitx"]>numpix.x, "gitx"] <- numpix.x
part[part[, "gity"]>numpix.y, "gity"] <- numpix.y
part[part[, "gitx"]<1, "gitx"] <- 1
part[part[, "gity"]<1, "gity"] <- 1

# get different resolutions for surface grids depending on range in x and y and on particle number for each timestep
# get selector for first and last row w/ a given btime
selfirst <- c(T, diff(part[, "btime"])>0)
selast <- c(diff(part[, "btime"])>0, T)

max.x <- part[order(part[, "btime"], part[, "gitx"]), ][selast>0, "gitx"]
min.x <- part[order(part[, "btime"], part[, "gitx"]), ][selfirst>0, "gitx"]
max.y <- part[order(part[, "btime"], part[, "gity"]), ][selast>0, "gity"]
min.y <- part[order(part[, "btime"], part[, "gity"]), ][selfirst>0, "gity"]
btime <- part[order(part[, "btime"], part[, "gity"]), ][selfirst>0, "btime"]
names(max.x) <- NULL; names(min.x) <- NULL; names(max.y) <- NULL; names(min.y) <- NULL

# now get information back in format for all timesteps and index-numbers
minmax.yx <- cbind(btime, max.x, min.x, max.y, min.y)
minmax.yx <- merge(part[, c("btime", "index")], minmax.yx, by="btime")
max.x <- minmax.yx[, "max.x"]
min.x <- minmax.yx[, "min.x"]
max.y <- minmax.yx[, "max.y"]
min.y <- minmax.yx[, "min.y"]
names(max.x) <- NULL; names(min.x) <- NULL; names(max.y) <- NULL; names(min.y) <- NULL



####################################################################################################
# Step 7: Define Emission grid wrt dynamic gridding
####################################################################################################

# Call 'getgrid' to get correct emission grid--necessary b/c emission grid is too large, so divided into several diff objects
# use getgridp.ssc function: don't allow the resolution to get finer at earlier backtime; use cummax(ran.x)

gridresult <- getgridp(min.x, max.x, min.y, max.y, numpix.x, numpix.y, coarse)
emissname <- paste(gridresult[, "xpart"], gridresult[, "ypart"], gridresult[, "gridname"], sep="")
# Extract appropriate emissions within each emission grid--do one grid at a time b/c reduces # of times grid has to be accessed
coarsex <- c(1,1,2,2,2,4,4,4,8,8,8,16,16,16,32,32)  # factors by which grids have been made coarser
coarsey <- c(1,2,1,2,4,2,4,8,4,8,16,8,16,32,16,32)  # e.g., '4' means grid is 4 times coarser



####################################################################################################
# Step 8: Setup tracer and vegetation type grid loop
####################################################################################################

# loop over different surface grids
# initialize vector of pointers (tracer ids) for all tracers (not including biospheric tracers yet)
tr.ids <- 1:length(tracers)
result <- cbind(part, matrix(ncol=length(tr.ids), nrow=length(part[, "btime"])))
dimnames(result) <- list(NULL, c(dimnames(part)[[2]], paste("v", tr.ids, sep="")))

# mapping to flux fields
for (speci in tracers) { # all fossil fuel emissions, do the ones with netCDF format
   fnum<-which(tracers==speci)
   if (ncdfTF[tolower(speci)]) { # read surface fluxes from netcdf files
        # BUDGET
              # Take emission values at particle positions; multiply by "foot", i.e. sensitivity of mixing ratio changes to fluxes,
              # in ppm/(micro-mol/m^2/s)
      # need time in YY MM DD hh mm, not local, but GMT!
      gmtime <- weekdayhr(yr4, mon, day, hr, -result[, "btime"]*60, diffGMT=rep(0, dim(result)[1])) # last column is weekday, sunday=0, monday=1 etc.
      fdate <- as.vector(rbind(gmtime[, "yr"], gmtime[, "mon"], gmtime[, "day"], gmtime[, "hr"], rep(0, dim(result)[1])))
      if (sum(is.na(fdate))>0) stop("NA in fdate!")
      # print(fdate)
      # Convert mins & maxes from coordinate values to rows & columns
      shrink.x <- coarsex[gridresult[, "gridname"]]; shrink.y <- coarsey[gridresult[, "gridname"]]
      x <- ceiling(part[, "gitx"]/shrink.x)
      y <- ceiling(part[, "gity"]/shrink.y)
      emiss <- NULL
      # print(pathname)
      # assignr("testncdf", list(gmtime=gmtime, fdate=fdate, i=x, j=y, ires=shrink.x, jres=shrink.y, emiss=emiss), pathname)
      # stop("testing")
      # print("date 1")
      # print(date())
#      if(speci!="cofire")emiss <- get.fossEU.netcdf(fdate, i=x, j=y, ires=shrink.x, jres=shrink.y, spec=speci)
#      if(speci=="cofire")emiss <- get.fireBARCA.netcdf(fdate, i=x, j=y, ires=shrink.x, jres=shrink.y, spec=speci)
      emiss <- get.Emis.netcdf(fdate, i=x, j=y, ires=shrink.x, jres=shrink.y, spec=speci,numpix_x=numpix.x,numpix_y=numpix.y)
      EMCO <- emiss*part[, "foot"]
      # also multiplied by dmass (accumulated weight of particles due to mass violation, normalized by average dmass to conserve total mass over time)
      if (dmassTF) EMCO <- EMCO*part[, "dmass"]
      result[, paste("v", fnum, sep="")] <- EMCO
   } # of if (ncdf) {
} # for different surface grid

# mapping to VPRM fields (EVI and LSWI indices and vegetation fraction)
if ((bios=="VPRM"|bios=="GSB")&"co2"%in%tracers) { # VPRM with fortran call to directly access netcdf files
        # BUDGET
              # Take emission values at particle positions; multiply by "foot", i.e. sensitivity of mixing ratio changes to fluxes,
              # in ppm/(micro-mol/m^2/s)
      # need time in YY MM DD hh mm, not local, but GMT!
      gmtime <- weekdayhr(yr4, mon, day, hr, -result[, "btime"]*60, diffGMT=rep(0, dim(result)[1])) # last column is weekday, sunday=0, monday=1 etc.
      fdate <- as.vector(rbind(gmtime[, "yr"], gmtime[, "mon"], gmtime[, "day"], gmtime[, "hr"], rep(0, dim(result)[1])))
      if (sum(is.na(fdate))>0) stop("NA in fdate!")
      # Convert mins & maxes from coordinate values to rows & columns
      shrink.x <- coarsex[gridresult[, "gridname"]]; shrink.y <- coarsey[gridresult[, "gridname"]]
      x <- ceiling(part[, "gitx"]/shrink.x)
      y <- ceiling(part[, "gity"]/shrink.y)
      if(bios=="VPRM"){
        vprmset <- NULL
        # get EVI parameters, vegetation specific
        vprmset  <- get.modis.netcdf(fdate, i=x, j=y, ires=shrink.x, jres=shrink.y, dpath=evilswipath)
        eviset <- vprmset$EVI
        lswiset <- vprmset$LSWI
        eviMaxVec <- vprmset$EVI_amax
        eviMinVec <- vprmset$EVI_amin
        lswiMaxVec <- vprmset$LSWI_amax
        lswiMinVec <- vprmset$LSWI_amin
        veg.fra<-vprmset$VEG_FRA
      } else {
        veg.fra <- get.vegfrac.netcdf(fdate, i=x, j=y, ires=shrink.x, jres=shrink.y, dpath=evilswipath)

      }
      infl.veg <- veg.fra*part[, "foot"] #influence from specific vegetation type
      # also multiplied by dmass (accumulated weight of particles due to mass violation, normalized by average dmass to conserve total mass over time)
      if (dmassTF) infl.veg <- infl.veg*part[, "dmass"]
} # of if ((bios=="VPRM"|bios=="GSB")&"co2"%in%tracers) {

dimnames(result) <- list(NULL, dimnames(result)[[2]])
###################################


####################################################################################################
# Step 11: Reclassify vegetation, calculate a priori vegetation flux
####################################################################################################

if (landcov == "IGBP"&bios == "GSB") {
   nReclss <- 5 #number of classes for flux calculation
   reclss <- c(1,1,1,1,1,2,2,2,2,2,4,3,0,3,5,2,5) #conversion from landcov class to flux class

   # GSB is based on grouping IGBP classes
   # -------------------------------------
   # vgroup1: Forrests: IGBP 1, 5, 4, 2, 3 (not existent)
   # vgroup2: shrublands etc.:7,10,8,16,6,9
   # vgroup3: Crops etc.: 14,12
   # vgroup4: wetland: 11
   # vgroup5: Water IGBP 17 (water), 15 (Snow and Ice)
   # Rest:13 (urban) -- are assumed to have zero vegetative flux influence
}

if (landcov == "GLCC"&bios != "GSB") {
   nReclss <- 8
   reclss <- c(1,1,2,2,3,4,4,5,5,7,7,6,8,7,8,7,8)
   # CHECK THIS
   # NEEDS WORK TO DEFINE VEGETATION GROUPS
   ## Value        Devan Category
   # 1        Evergreen Forest
   # 2        Deciduous Forest
   # 3        Mixed Forest
   # 4        Shrublands
   # 5        Savannas
   # 6        Croplands
   # 7        Grasslands
   # 8        Water, Snow/ice, other
   # Rest(urban, barren) assumed to have zero vegetative flux influence
}

if (landcov == "SYNMAP.VPRM8"&bios == "GSB") {
   nReclss <- 5
   reclss <- c(1,1,1,2,2,3,2,5)
}
if (landcov == "SYNMAP.VPRM8"&bios == "VPRM") {
   nReclss <- 9
   reclss <- 1:9
}
# Devan vegetation classes, updated
#                 VPRM Class    STILT-VPRM class
# Evergreen A                1A                1
# Evergreen B                1B                2
# Evergreen C                1C                3
# Evergreen D                1D                4
# Deciduous                   2                5
# Mixed forest                3                6
# Shrubland                   4                7
# Savanna                     5                8
# Cropland-Soy               6A                9
# Cropland-Maize             6B                9
# Grassland                   7                10
# Peatlands                   8                11
# Others                      9                12
#

if (landcov == "DVN") {
   nReclss <- 12
   reclss <- 1:12
}

vegresult<-NULL
if ("co2"%in%tracers&bios!="") {
   # Conversion factor of SW radiation (W/m2) to par for use in VPRM. Based on accumulated fit of tower data
   # In NLDAS, EDAS gets ingested to supplement GOES at low sun angles. EDAS is an overestimation,
   if (pre2004) {
      swtoparnldas <- 1.60
      swtoparedas <- 1.5 # placeholder
   } else {
      swtoparnldas <- 1.89
      swtoparedas <- 1.61
   }
        
   # Get light and temperature
   swrad <- result[, "swrad"]
   swradpar <- swtoparedas*swrad

   tempAir <- result[, "temp0"]-273.15
   tempAir[which(tempAir<(-100)|tempAir>(100)|is.nan(tempAir))] <- NA
                
   if (bios == "GSB") {
      if (!linveg) {
         dresp.dT <- c(dlambda.veg[, "drdt"], dlambda.veg[2, "drdt"])
         a3 <- c(dlambda.veg[, "a3"], dlambda.veg[2, "a3"])
         a4 <- c(dlambda.veg[, "a4"], dlambda.veg[2, "a4"])
         # get parameterized GEE and R
      } else { # fluxes linear in temp and radiation
         dresp.dT <- c(dlambda.simp.veg[, "drdt"], dlambda.simp.veg[2, "drdt"])
         a3 <- c(dlambda.simp.veg[, "a3"], dlambda.simp.veg[2, "a3"])
      }# of if not linear fluxes as fct of temp and radiation

      # calculate t and par dependent fluxes
      for (k in 1:(nReclss-1)) {#exclude last flux class, since no biosphere flux assumed (water etc)
         # Assume last reclass category is water
         # Determine total influence for all vegetation from this uberclass
         selCols <- which(reclss == k)
         if (length(selCols)>1)
                 v.tot <- apply(infl.veg[, selCols],1, sum)
         else if (length(selCols) == 1)
                 v.tot <- infl.veg[, selCols]
         else # not all classes are present in all areas
#                 stop("Sel Cols missing--Veg types off")
                 v.tot <- infl.veg[,1]*NA
         if (!linveg)
                 v.gee <- a3[k]*swrad/(swrad+a4[k])*v.tot
         if (linveg)
                 v.gee <- a3[k]*swrad*v.tot

         v.resp <- dresp.dT[k]*tempAir*v.tot

         gsbMat <- cbind(v.tot, v.gee, v.resp)
         dimnames(gsbMat)[[2]]<-paste(dimnames(gsbMat)[[2]],k,sep="")
         vegresult <- cbind(vegresult, gsbMat)
      }
   }

   if (bios == "VPRM") {
      # read in values for all classes
      if (landcov != "DVN") {
              vprmConstants <- getr(vprmconstantsname, path=vprmconstantspath)
      } else {

              vprmConstants <- getr("publishedvprmconstants.12", path=vprmconstantspath)
      }


      for (k in 1:(nReclss-1)) {

         # Determine total influence for all vegetation in this uberclass
         selCols <- which(reclss == k)
         if (length(selCols)>1)
                 v.tot <- apply(infl.veg[, selCols],1, sum)
         else if (length(selCols) == 1)
                 v.tot <- infl.veg[, selCols]
         else # not all classes are present in all areas
                 stop("Sel Cols missing--Veg types off")
         # get parameters
         if(names(vprmConstants)[1]=="lambdaGPP"){ #devan style (precompute PAR from SW, then use parameters for PAR)
            lambdaGPP <- vprmConstants$lambdaGPP[k]
            alphaResp <- vprmConstants$alphaResp[k]
            betaResp <- vprmConstants$betaResp[k]
            parZero <- vprmConstants$parZero[k]
            radScalar <- swradpar/(1 + (swradpar/parZero)) # Include radiation factor
         } else { #Jena style, differentiate between PAR and SW, and use SW
            lambdaGPP <- vprmConstants$lambdaGPP.sw[k]
            alphaResp <- vprmConstants$alphaResp[k]
            betaResp <- vprmConstants$intResp[k]
            parZero <- vprmConstants$swradZero[k]
            radScalar <- swrad/(1 + (swrad/parZero)) # Include radiation scalar
         }
         # get EVI parameters, vegetation specific
         evi <- eviset[, k]
         lswi <- lswiset[, k]

         # get temperature
         # temperature parameters, assume calculations in degrees C
         tempOpt <- vprmConstants$tempOpt[k]
         tempMax <- vprmConstants$tempMax[k]
         tempMin <- vprmConstants$tempMin[k]
         if("tempRLow"%in%names(vprmConstants)) #min temp with resp 
            tempRLow <- vprmConstants$tempRLow[k]
         else
            tempRLow <- 0*vprmConstants$tempMin[k] #set to zero

         tempScalar <- ((tempAir-tempMin)*(tempAir-tempMax))/
                     (((tempAir-tempMin)*(tempAir-tempMax))-(tempAir-tempOpt)^2)
         tempScalar[which(tempAir>tempMax | tempAir < tempMin)] <- 0

         eviMax <- eviMaxVec[, k]
         eviMin <- eviMinVec[, k]
         lswiMax <- lswiMaxVec[, k]
         lswiMin <- lswiMinVec[, k]

         # Need to eliminate ALL noise in evi and lswi--when greater than max values, set to max value, same for min (on evi)
         # with extremely low values.

#remove later
#print(paste("k:",k))
#print(paste("Trajecvprm: checking if evi threshold required",sum(evi>eviMax),sum(evi<eviMin)))
#if(sum(evi>eviMax)>0)print(paste("locations evi>max:",which(evi>eviMax, arr.ind = TRUE)))
#print(paste("Trajecvprm: checking if lswi threshold required",sum(lswi>lswiMax),sum(lswi<lswiMin)))
#end remove later

         evi[which(evi>eviMax)] <- eviMax[which(evi>eviMax)]
         evi[which(evi<eviMin)] <- eviMin[which(evi<eviMin)]
         lswi[which(lswi>lswiMax)] <- lswiMax[which(lswi>lswiMax)]
         lswi[which(lswi<lswiMin)] <- lswiMin[which(lswi<lswiMin)]

         # modification for so-called "xeric systems", comprising shrublands and grasslands
         # these have different dependencies on ground water.

         if ((landcov == "DVN"&is.element(k, c(7,10)))|(landcov != "DVN"&is.element(k, c(4,7)))) {
                 wScalar <- (lswi-lswiMin)/(lswiMax-lswiMin)
         } else {
                 wScalar <- (1+lswi)/(1+lswiMax)
         }


         # Case 1: Evergreens--pScalar always = 1
         # Case 2: Grasslands and Savannas -- pScalar never equals 1
         # Case 3: Others -- pScalar = 1 for growing season only as determined from evi threshold
         #        evithreshold = eviMin  + (0.55*range); range = eviMax - eviMin

         pScalar <- (1+lswi)/2
         if ((landcov == "DVN"&is.element(k,1:4))|(landcov != "DVN"&k == 1)) { # if evergreen
                 phenologyselect <- 1:length(evi)
         }
         if ((landcov == "DVN"&is.element(k, c(5,6,7,9,12)))|(landcov != "DVN"&is.element(k, c(2,3,4,6,8)))) { # if decid, mixed, shrub, crop, or other
                 threshmark <- 0.55
                 evithresh <- eviMin+(threshmark*(eviMax-eviMin))
                 phenologyselect <- which(evi>evithresh)
         }
         # by default, grasslands and savannas and peatlands never have pScale=1

         pScalar[phenologyselect] <- 1

         # Error Check evi, lswi, tempScalar
         evi[which(is.na(evi))] <- 0   #checked for Jena implementation: na check not really required
         lswi[which(is.na(lswi))] <- 0 #checked for Jena implementation: na check not really required
         tempScalar[which(is.na(tempScalar))] <- 0  #checked for Jena implementation: na check not really required
         tempScalar[which(tempScalar<0)] <- 0
         # Determine GPP
         # NOTE UNITS--VPRM outputs GPP and Respiration in umol/m2/s (conveniently, what is needed here); when multiplied by
         #                influence (ppm/(umol/m2/s)) get ppm
         xiaoGPP <- lambdaGPP*tempScalar*wScalar*pScalar*evi*radScalar*v.tot*-1

         # want vegetative uptake to be negative with respect to the atmosphere, so multiply by negative one
         # for symmetry
         xiaoGPP[which(is.na(xiaoGPP))] <- 0

         v.gee <- xiaoGPP

         tempForR <- tempAir
         tempForR[which(tempForR<tempRLow)] <- tempRLow

         devanResp <- alphaResp*tempForR*v.tot+betaResp*v.tot
         devanResp[which(is.na(devanResp))] <- 0
         v.resp <- devanResp


         # Output
         msbMat <- cbind(v.tot, v.gee, v.resp)
         dimnames(msbMat)[[2]]<-paste(dimnames(msbMat)[[2]],k,sep="")
         vegresult <- cbind(vegresult, msbMat)

      } # for k
   } # of if (bios == "VPRM")
}# if co2

#special case H2 fluxes from soil: use soil moisture parameterization in presence of vegetation
if ("h2"%in%tracers&!ncdfTF["h2"]) {
   # For volumetric soil moisture 0- 0.15:
   # -Flux=80*soil moisture
   # For soil moisture 0.15-0.50
   # -Flux=-34.3*soil moisture + 17.14
   # The units are nmol m-2 s-1.  You get max flux of 12 at 0.15 linearly
   # decreasing to 0 at soil moisture values of 0 and 0.5.

   idveg <- which(is.element(reclss,1:(nReclss-1))) # Assume vegetation in all categories except the 'final' landclass

   v.h2 <- apply(infl.veg[, idveg],1, sum) # add all influences from this vegetation class
   if (!("solw"%in%dimnames(result)[[2]]))
          stop("need soil moisture for H2")
   vsm <- result[, "solw"]/1000
   v.h2 <- v.h2*( (vsm<0.15)*(80*vsm) + (vsm>=0.15&vsm<0.5)*(17.14-34.3*vsm))/1000 # from nmol to mumol

   fnum<-which(tracers=="h2")
   cnameAdd <- paste("v", fnum, sep="")
   v.h2 <- matrix(v.h2, ncol=1, dimnames=list(NULL, cnameAdd))

   result <- cbind(result, v.h2)
}

####################################################################################################
# Step 12: add CO2 boundary values
####################################################################################################

        # get selector for last observation for each particle
ordert <- order(result[, "index"], result[, "btime"]) # order first by index, then by time
ordern <- order(result[ordert, "btime"], result[ordert, "index"]) # to recreate original order
delbte <- c(diff(result[ordert, "btime"]), -1000)[ordern] # timestep will be negative at last obs. for each particle
selend <- delbte<0

        # get selector for first observation for each particle
delbte <- c(-1000, diff(part[ordert, "btime"]))[ordern] # timestep will be negative at first obs. for each particle
selfirst <- delbte<0

        # get initial conditions for CO and CO2, read with read <- co2 <- bg <- aug00.ssc and read <- co <- bg <- aug00.ssc
        # co2.ini and co.ini objects are 3d arrays, agl*lat*sasdate; sasdate is 0 at 1/1/1960 (i.e. elapsed days since 1/1/1960, same as julian(m, d, y))
aglg <- round(result[selend, "agl"]/500) # get agl in 500 m intervals
latg <- round((result[selend, "lat"]-10)/2.5) # get lat in 2.5 deg intervals, starting at 10 deg
cini <- NULL
for(tr in tracers) cini <-  cbind(cini,result[,1]*0) #default: no boundary condition specified
colnames(cini) <- paste(tracers, "ini", sep="") # give advected boundary mixing ratio names

tracers.clim<-tracers[inikind[tracers]=="climat"]
if(length(tracers.clim)>0){ #are there climatological boundary files to be used
        # get 1st boundary field to set things up
   if (!existsr(paste(tracers.clim[1], ".ini", sep=""), pathname)){
           print("need to use read.bg() to get boundary condition")
           stop(paste("need ", tracers.clim[1], " boundary condition", sep=""))
   }
   ini <- getr(paste(tracers.clim[1], ".ini", sep=""), pathname)

   startday <- as.numeric(dimnames(ini)[[3]][1]); delday <- as.numeric(dimnames(ini)[[3]][2])-startday

   sasdate <- julian(mon, day, yr4)-result[selend, "btime"]/24-startday # days elapsed since beginning of .ini file
   pointer <- cbind(aglg+1, latg+1, round(sasdate/delday)+1) # array indices must start with 1
           # use constant ini field when no initial data available
   if (any(pointer[,3]>dim(ini)[3]))
           cat("Trajecvprm(): extrapolating ", tracers.clim[1], ".ini, need later times\n", sep="")
   if (any(pointer[,3]<1))
           cat("Trajecvprm(): extrapolating ", tracers.clim[1], ".ini, need earlier times\n", sep="")
   pointer[pointer[,3]>dim(ini)[3],3] <- dim(ini)[3]
   pointer[pointer[,3]<1,3] <- 1

           # limit to lat range of initial field
           # limit pointer to valid values for latg and aglg
   pointer[pointer[,1]>dim(ini)[1],1] <- dim(ini)[1]
   pointer[pointer[,2]>dim(ini)[2],2] <- dim(ini)[2]
   pointer[pointer[,1]<1,1] <- 1
   pointer[pointer[,2]<1,2] <- 1

   for (tr in tracers.clim[1:length(tracers.clim)]) {
      if (tr!=tracers.clim[1]) {
         if (!existsr(paste(tr, ".ini", sep=""), pathname)) {
            stop(paste("warning: need ", tr, " boundary condition", sep=""))
         }
         ini <- getr(paste(tr, ".ini", sep=""), pathname)
      }
      cini[selend,paste(tr, "ini", sep="")] <- ini[pointer]
   }
} #if(length(tracers.clim)>0){ #are there climatological boundary files to be used

#combine flux tracers, LBC tracers and biospheric tracers
result <- cbind(result, cini, vegresult)

if ("co2"%in%tracers & inikind["co2"] == "CT") {
   cat("Trajecvprm: using CarbonTracker initial values.\n")
   result <- get.CarbonTracker.netcdf(yr4=yr4, mon=mon, day=day, hr=hr, co2inifile=inifile["co2"],
                              result=result, result.sel=selend)
} else if ("co2"%in%tracers & inikind["co2"] == "TM3") {
   cat("Trajecvprm: using TM3 initial values.\n")
   ftype<-substring(inifile["co2"],nchar(inifile["co2"])-1,nchar(inifile["co2"]))
   if(ftype==".b")result <- get.TM3.bin(yr4=yr4, mon=mon, day=day, hr=hr, co2inifile=inifile["co2"],
                    result=result, result.sel=selend)
   if(ftype=="nc")result <- get.TM3.netcdf(yr4=yr4, mon=mon, day=day, hr=hr, co2inifile=inifile["co2"],
                    result=result, result.sel=selend)
} else if ("co2"%in%tracers & inikind["co2"] == "LMDZ") {
   cat("Trajecvprm: using LMDZ initial values.\n")
   result <- get.LMDZ.netcdf(yr4=yr4, mon=mon, day=day, hr=hr, co2inifile=inifile["co2"],
                    result=result, result.sel=selend)
}


if ("co"%in%tracers & inikind["co"] != "GEMS") {                                      # CO: need also advected boundary value with no chemistry
   coinio <- result[, "coini"]                              # no chemistry yet
   result[selend, "coini"] <- result[selend, "CO.fact"]*coinio[selend]+result[selend, "CO.frCH4"] #add chemistry (linearized)
   result <- cbind(result, coinio)
} else if ("co"%in%tracers & inikind["co"] == "GEMS") {
   cat("Trajecvprm: using GEMS CO initial values.\n")
   result <- get.GEMS_CO.netcdf(yr4=yr4, mon=mon, day=day, hr=hr, co2inifile=inifile["co"],
                    result=result, result.sel=selend,spec="coini")
   coinio <- result[, "coini"]                              # no chemistry yet
   result[selend, "coini"] <- result[selend, "CO.fact"]*coinio[selend]+result[selend, "CO.frCH4"] #add chemistry (linearized)
   result <- cbind(result, coinio)
}
if ("cofire"%in%tracers & inikind["cofire"] == "GEMS") {                                      # CO: need also advected boundary value with no chemistry
   cat("Trajecvprm: using GEMS CO initial values for cofire.\n")
   result <- get.GEMS_CO.netcdf(yr4=yr4, mon=mon, day=day, hr=hr, co2inifile=inifile["cofire"],
                    result=result, result.sel=selend,spec="cofireini")
}


if ("ch4"%in%tracers & inikind["ch4"] == "TM3") {
   cat("Trajecvprm: using TM3 initial values for CH4.\n")
   result <- get.TM3CH4.bin(yr4=yr4, mon=mon, day=day, hr=hr, ch4inifile=inifile["ch4"],
                    result=result, result.sel=selend)
}

dimnames(result) <- list(NULL, dimnames(result)[[2]])

####################################################################################################
# Step 13: integrate influence / tracer contributions over time and trajectory ensemble
####################################################################################################

# add all influences
v <- result[, substring(dimnames(result)[[2]],1,1) == "v"]
sv <- v*0 # initialize sum of veg. influence

if(is.matrix(v)){ #multiple tracers
   dimnames(sv) <- list(NULL, paste("s", dimnames(sv)[[2]], sep=""))
   for (veg in colnames(v)) {
      cumsamp <- tapply(v[, veg], result[, "index"], cumsum)
      cumsamp <- unlist(cumsamp)[ordern]
      sv[, paste("s", veg, sep="")] <- cumsamp
   }# for veg in different surface flux variables
   result <- cbind(result, sv)
} else { #only single tracer
   cumsamp <- tapply(v, result[, "index"], cumsum)
   cumsamp <- unlist(cumsamp)[ordern]
   sv1 <- cumsamp
   result <- cbind(result, sv1)
}
#################################
dimnames(result) <- list(NULL, dimnames(result)[[2]])

# drop all vegetation fluxes, keep only the accumulated fluxes
resulto<-result
result <- result[, substring(dimnames(result)[[2]],1,1) != "v"]
dimnames(result) <- list(NULL, dimnames(result)[[2]])


# keep only information which is required: last entry position and accumulated fluxes for coarse vegetation classes, and CO and CO2 fossil fuel emissions
lastresult <- result[selend>0, ]
dimnames(lastresult) <- list(NULL, dimnames(lastresult)[[2]])

# keep also certain parameters at beginning of particle trajectory
firstresult <- result[selfirst>0, ]

# give each particle the same weight throughout the run,
# discard a complete particle run when particle stops before max. btime and before it reaches boundary: done with object "part" higher up
# get stdeviation for co2ini, coini
if (length(tracersini) > 1) {
   sdini <- apply(lastresult[, tracersini], 2, stdev, na.rm=T)
} else {
   sdini <- stdev(lastresult[, tracersini])
}

# get number of particles ending insight area before btime=btimemax; allow for ~2.5 degree band
# at border and 1 day time
endinarea <- (lastresult[, "gitx"]>10&lastresult[, "gitx"]<(numpix.x-10))&(lastresult[, "gity"]>10&lastresult[, "gity"]<(numpix.y-10))&lastresult[, "btime"]<nhrs-0.05*nhrs
endinarea <- sum(endinarea)
# get ensemble averages
lastresult <- apply(lastresult,2, mean)
firstresult <- apply(firstresult,2, mean)

# keep only ident, btime, lat, lon, co2ini, coini, mixed layer height and ground height, at start,
# at ending in area, sdev. co and co2 ini, influence from water, and veg. influences
output <- c(pos, lastresult[c("btime", "lat", "lon", "agl")], firstresult[c("zi", "grdht")], endinarea,
            lastresult[tracersini],sdini, lastresult[substring(names(lastresult),1,2) == "sv"])
names(output) <- names.output

if (detailsTF) {
  resulthr<-get.mean.traj(resulto,hr=2) #only for looking at data...
  selnam <- substring(dimnames(resulthr)[[2]],1,1) == "v"
  dimnames(resulthr)[[2]][selnam] <- names.output[(length(names.output)-sum(selnam)+1):length(names.output)]
  assignr(paste(ident,"resulthr", sep=""), resulthr,pathname, printTF=TRUE)

  selnam <- substring(dimnames(result)[[2]],1,2) == "sv"
  dimnames(result)[[2]][selnam] <- names.output[(length(names.output)-sum(selnam)+1):length(names.output)]
  if(ncdfTF["co2"]&bios=="VPRM"){
    vprmstuff<-cbind(eviset,eviMaxVec,eviMinVec)
    dimnames(vprmstuff)[[2]]<-c(paste(output.veg,".evi",sep=""),paste(output.veg,".evimax",sep=""),paste(output.veg,".evimin",sep=""))
    result<-cbind(result,vprmstuff)
  }
  assignr(paste(ident, "result", sep=""), result, pathname, printTF=TRUE)
}

gc(verbose=F)
return(output)

}
