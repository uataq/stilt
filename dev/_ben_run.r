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
#

paths <- list(rsc = '/uufs/chpc.utah.edu/common/home/u0791983/stilt/Rsc/',
              met = '/uufs/chpc.utah.edu/common/home/u0791983/met/',
              run = '/uufs/chpc.utah.edu/common/home/u0791983/stilt/STILT_Exe/Copy0')

t_start <- '2015-07-18 00:00:00'
t_end   <- '2015-07-18 01:00:00'
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
library(dplyr)
library(uataq)

# Model run timing -------------------------------------------------------------
time <- data_frame(posix = seq(from = t_start %>% as.POSIXct(tz='UTC'),
                               to   = t_end   %>% as.POSIXct(tz='UTC'),
                               by   = 'hour'),
                   yyyy  = format(posix, '%Y'),
                   mm    = format(posix, '%m'),
                   dd    = format(posix, '%d'),
                   HH    = format(posix, '%H'))


# Met path auto symlink --------------------------------------------------------
met_sym <- paste0('~/.met_', gsub(' ', '', t_start))
# system(paste('ln -s', paths$met, met_sym))

# Find met files ---------------------------------------------------------------
met_files <- dir(met_sym, pattern = '.*\\.arl', full.names = T)
# if (length(met_files) < 1)
  # stop(paste('No met files found in', met_sym))


# Time step loop ---------------------------------------------------------------
uataq::clc()
uataq::br(2)
pb <- txtProgressBar(title = 'Progress:', style=3)
for (i in 1:nrow(time)) {
  print(paste('Time:', format(time$posix[i], '%Y-%m-%d %H:%M', 'UTC')))
  # k <- format(time$posix[i], '%Y-%m-%d %H:%M', 'UTC')
  j <- (i - 1) / nrow(time)
  setTxtProgressBar(pb, j)
  Sys.sleep(0.5)
  uataq::clc()
  uataq::br(2)
}




# system(paste('unlink', met_sym))
