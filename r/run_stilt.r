# STILT R Executable
# For documentation, see https://github.com/benfasoli/stilt
# Ben Fasoli

# User inputs ------------------------------------------------------------------
project <- '{{project}}'
stilt_wd <- file.path('{{wd}}', project)
lib.loc <- .libPaths()[1]

rm_dat  <- T
n_nodes <- 1
n_cores <- 1

slurm   <- F
slurm_options <- list(
  time = '300:00:00',
  account = 'lin-kp',
  partition = 'lin-kp'
)

# Simulation timing, yyyy-mm-dd HH:MM:SS
t_start <- '2015-06-17 00:00:00'
t_end   <- '2015-06-17 00:00:00'
run_times <- seq(from = as.POSIXct(t_start, tz='UTC'),
                 to   = as.POSIXct(t_end, tz='UTC'),
                 by   = 'hour')

# Receptor locations
lati <- 40.766189
long <- -111.847672
zagl <- 10

# Expand the run times, latitudes, and longitudes to form the unique receptors
# that are used for each simulation
receptors <- expand.grid(run_time = run_times, lati = lati, long = long,
                         zagl = zagl, KEEP.OUT.ATTRS = F, stringsAsFactors = F)

# Meteorological data input
met_directory   <- '/uufs/chpc.utah.edu/common/home/lin-group6/hrrr/data/west'
met_file_format <- '%Y%m%d.%Hz.hrrra'
n_met_min       <- 5

# Model control
run_trajec <- T
n_hours    <- -24
convect    <- F
delt       <- 0
numpar     <- 200
outdt      <- 0
timeout    <- 3600
varsiwant  <- c('time', 'indx', 'long', 'lati', 'zagl', 'sigw', 'tlgr', 'zsfc',
                'icdx', 'temp', 'samt', 'foot', 'shtf', 'tcld', 'dmas',
                'dens', 'rhfr', 'sphu', 'solw', 'lcld', 'zloc', 'dswf', 'wout',
                'mlht', 'rain', 'crai')

# Transport and dispersion settings
iconvect    <- 0
isot        <- 0
mgmin       <- 2000
ndump       <- 0
nturb       <- 0
outfrac     <- 0.9
random      <- 1
tlfrac      <- 0.1
tratio      <- 0.9
veght       <- 0.5
w_option    <- 0
winderrtf   <- F
z_top       <- 25000
zicontroltf <- 0

# Footprint grid settings
xmn <- -114.5
xmx <- -109
ymn <- 37
ymx <- 42
xres <- 0.1
yres <- xres
time_integrate <- F


# Source dependencies ----------------------------------------------------------
setwd(stilt_wd)
source(file.path(stilt_wd,'r/dependencies.r'))


# Met path symlink -------------------------------------------------------------
# Auto symlink the meteorological data path to the working directory to
# eliminate issues with long (>80 char) paths in fortran. Note that this assumes
# that all meteorological data is found in the same directory.
if ((nchar(paste0(met_directory, met_file_format)) + 2) > 80) {
  met_loc <- file.path(path.expand('~'), paste0('m', project))
  if (file.exists(met_loc))
    system(paste('unlink', met_loc))
  system(paste('ln -s', met_directory, met_loc))
} else met_loc <- met_directory


# Run trajectory simulations ---------------------------------------------------
# Gather varsiwant into a single character string and fork the process to apply
# simulation_step() to each receptor across n_cores and n_nodes
validate_varsiwant(varsiwant)
if (!is.null(varsiwant[1]))
  varsiwant <- paste(varsiwant, collapse = '/')

output <- stilt_apply(X = 1:nrow(receptors), FUN = simulation_step,
                      slurm = slurm, slurm_options = slurm_options,
                      n_cores = n_cores, n_nodes = n_nodes, rm_dat = rm_dat,
                      delt = delt, iconvect = iconvect, isot = isot,
                      lib.loc = lib.loc, met_file_format = met_file_format,
                      met_loc = met_loc, mgmin = mgmin, n_hours = n_hours,
                      n_met_min = n_met_min, ndump = ndump, nturb = nturb,
                      numpar = numpar, outdt = outdt, outfrac = outfrac,
                      run_trajec = run_trajec, r_run_time = receptors$run_time,
                      r_lati = receptors$lati, r_long = receptors$long,
                      r_zagl = receptors$zagl, random = random,
                      stilt_wd = stilt_wd, time_integrate = time_integrate,
                      timeout = timeout, tlfrac = tlfrac, tratio = tratio,
                      varsiwant = varsiwant, veght = veght, w_option = w_option,
                      winderrtf = winderrtf, zicontroltf = zicontroltf,
                      z_top = z_top, xmn = xmn, xmx = xmx, xres = xres,
                      ymn = ymn, ymx = ymx, yres = yres)

q('no')
