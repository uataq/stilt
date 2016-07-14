# Ben Fasoli

# User inputs ------------------------------------------------------------------
# paths - gives a list of paths to
#  rsc - directory containing R functions
#  met - directory containing arl formatted met files
#  run - directory of desired stilt hymodelc copy
# t_start - simulation start time
# t_end   - simulation end time
# r_lat   - receptor latitude
# r_lon   - receptor longitude
# r_agl   - receptor height above ground level
# n_hour  - number of hours to run each simulation (negative indicates backward)

paths <- list()
paths$stilt <- Sys.getenv('HOME')
paths$r     <- file.path(paths$stilt, 'STILT_modeling/stiltR/')
paths$met   <- file.path(paths$stilt, 'met/')
paths$run   <- file.path(paths$stilt, 'STILT_modeling/STILT_Exe/')
paths$out   <- file.path(paths$stilt, 'STILT_modeling/Output/')

# t_start <- '2010-07-18 16:00:00'
# t_end   <- '2010-07-18 17:00:00'
t_start <- '2010-04-28 00:00:00'
t_end   <- '2010-04-28 06:00:00'

r_lat   <- 40.75
r_lon   <- -111.87
r_agl   <- 5

n_hour     <- -12
z_veg      <- 0.5
convect    <- F
step_size  <- 0
winderr    <- F
overwrite  <- T
delta_t    <- 2
mg_min     <- 2000
n_particle <- 1000


# Source dependencies ----------------------------------------------------------
# setwd(paths$r)
rsc <- dir(file.path(paths$r, 'src'), pattern = '.*\\.r$', full.names = T)
lapply(rsc, source)
load_libs('dplyr', 'uataq')

# Model run timing -------------------------------------------------------------
time <- data_frame(posix = seq(from = t_start %>% as.POSIXct(tz='UTC'),
                               to   = t_end   %>% as.POSIXct(tz='UTC'),
                               by   = 'hour'),
                   yyyy  = format(posix, '%Y') %>% as.numeric,
                   yy    = format(posix, '%y') %>% as.numeric,
                   mm    = format(posix, '%m') %>% as.numeric,
                   dd    = format(posix, '%d') %>% as.numeric,
                   HH    = format(posix, '%H') %>% as.numeric)


# Met path auto symlink --------------------------------------------------------
met_sym <- paste0(paths$stilt, format(time$posix[1], tz = 'UTC', format = '/.met-%Y-%m-%d/'))
met_sym_init <- function() {
  system(paste('ln -s', paths$met, met_sym))
}
met_sym_close <- function() {
  system(paste('unlink', met_sym))
}

# Find met files ---------------------------------------------------------------
met_sym_init()
# met_files <- dir(met_sym, pattern = '\\.')#, pattern = '.*\\.arl', full.names = T)
met_files<-c("wrfout_d02.arl","wrfout_d01.arl","gdas1.apr10.arl")
if (length(met_files) < 1) {
  met_sym_close()
  stop(paste('No met files found in', met_sym))
}


# Time step loop ---------------------------------------------------------------
#uataq::clc()
uataq::br(2)
pb <- txtProgressBar(title = 'Progress:', style=3)

vars_trajec<-c("time","index","lat","lon","agl","grdht","foot","sampt","dmass",
               "zi","pres","sphu","temp","temp0","shtf","lhtf","sigmaw","dens","tcld",
               "totrain","convrain","wbar")

for (i in 1:nrow(time)) {
  print(paste('Time:', format(time$posix[i], '%Y-%m-%d %H:%M', 'UTC')))
  time_step <- time[i, ]
  
  info <- Trajec(yr = time_step$yy, mon = time_step$mm, day = time_step$dd, hr = time_step$HH,
                 lat = r_lat, lon = r_lon, agl = r_agl, nhrs = n_hour, doublefiles = T,
                 metd = c('fnl', 'awrf'), delt = delta_t, numpar = n_particle,
                 mgmin = mg_min, veght = z_veg, metfile = met_files, nummodel = 0,
                 metlib = met_sym, conv = convect, overwrite = overwrite,
                 outpath = paths$out, varsout = vars_trajec, rundir = paths$run,
                 zsg.name = NULL)
  
  dat <- getr(info['outname'], path = paths$out)
  timestamp <- format(time$posix[i], tz='UTC', format = '%Y%m%d%H')
  
  # Generate footprints --------------------------------------------------------
  foottimes<-0:24              #Want to generate a footprint for each hour
  zbot<-0				                     #Set equal to 0 if we want the surface influence volume
  ztop<-0				                       #Set equal to 0 if we want the surface influence volume
  numpix.x<-600                       #number of grid cells in the x direction (lon)
  numpix.y<-600                       #number of grid cells in the y direction (lat)
  lon.ll<--130.0                      #lower left corner of our footprint grid
  lat.ll<- 20.0                       #lower right corner of our footprint grid
  lon.res<-1/10                      #horizontal resolution of x grid
  lat.res<-1/10                      #horizontal resolution of y grid
  
  ident<-info["outname"]
  foot <- Trajecfoot(ident = ident, pathname = paths$out, foottimes = foottimes,
                     zlim = c(zbot, ztop), fluxweighting = NULL, coarse = 1, vegpath = NULL,
                     numpix.x = numpix.x, numpix.y = numpix.y, lon.ll = lon.ll, lat.ll = lat.ll,
                     lon.res = lon.res, lat.res = lat.res)
  
  xfoot.wfei <-apply(foot, c(1,2), sum)
  foot <- aperm(foot, c(2,1,3))
  netcdf.name = paste(paths$out, timestamp,"_footprint.nc",sep="")
  
  foot.lat <- as.numeric(rownames(foot))
  foot.lon <- as.numeric(colnames(foot))
  x <- ncdim_def( "Lon",  "degreesE", foot.lon)
  y <- ncdim_def( "Lat",  "degreesN", foot.lat)
  
  foot.var <- ncvar_def(name="footprint", units="PPM/(umoles/m^2 s)", list(y,x),longname="footprint")
  ncnew <- nc_create(filename=netcdf.name,vars=foot.var)
  ncvar_put(nc=ncnew,varid=foot.var,vals=xfoot.wfei)            #puts our variable into our netcdf file
  nc_close(ncnew)                                               #Closes our netcdf4 file
  
  # Update progress bar---------------------------------------------------------
  j <- (i - 1) / nrow(time)
  setTxtProgressBar(pb, j)
  #uataq::clc()
  uataq::br(2)
}

met_sym_close()
