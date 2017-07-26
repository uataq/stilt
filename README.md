# Quick-start

![STILT Footprint Example](https://air.utah.edu/~benfasoli/img/stilt-heb.png)

Done this before and just looking to start a new project?

1. Initialize myproject with `Rscript -e "uataq::stilt_init('myproject')"`  
1. Edit `r/run_stilt.r` to set receptor timing/locations and model parameters  
1. Run the model with `Rscript r/run_stilt.r`  

# Installation

STILT has been compiled to run on UNIX platforms (Mac, Linux). Required software includes

- R (v >= 3.25), [https://www.r-project.org/](https://www.r-project.org/)  
    - `dplyr` package, for speed and data manipulation
    - `parallel` package, for parallel computation on a single node
    - `raster` package, for raster-based spatial analysis
    - `rslurm` package, for parallel computation across multiple nodes
    - `uataq` package, for data manipulation
- One of the following Fortran compilers
    - Portland Group Compiler (pgf90)
    - Intel Fortran Compiler (ifort)
    - GNU Fortan Compiler (gfortran)
    - G95 Fortran Compiler (g95)
- Git, [https://git-scm.com/](https://git-scm.com/)

## Install methods

Two options exist to initialize a new STILT project. 

### R (preferred)

This method uses R to initialize a new project. `stilt_init()` in the Utah Atmospheric Trace gas & Air Quality (UATAQ) R package, which includes tools for trace gas data manipulation and analysis, is a wrapper around several system commands that do much of the heavy lifting. The `uataq` R package is available on [Github](https://github.com/benfasoli/uataq/) and can be installed in R using `devtools`.

```r
if (!require('devtools')) install.packages('devtools')
devtools::install_github('benfasoli/uataq')
```

A new STILT project can then be initialized in plain R code with

```r
uataq::stilt_init('myproject')
```

This function  
1. clones the stilt Github repository into a local `stilt` directory which is then renamed `myproject`  
2. builds the `permute.so` dynamic link library used to apply gaussian kernels for footprint output
3. compiles the hymodelc executable and moves to the `exe` directory
4. populates the project name and paths in `myproject/r/run_stilt.r`  

### Manual

While the R method is preferable since it streamlines the process of initializing new projects, the same can be accompolished manually. To reproduce the results above,

Clone the repository and set the name of the project

```
git clone https://github.com/benfasoli/stilt
mv stilt myproject
```

Compile and move the hymodelc executable to the `exe` directory using the provided `setup` script (or manually)

```
cd myproject
chmod +x setup
./setup
```

Finally, edit `r/run_stilt.r` setting `project`, `paths`, receptor timing/locations, and model parameters.


# Repository structure

The repository contains a number of directories that serve as the framework for the model. Initially, only the `r/` and `fortran/` directories are fully populated while `fortran/` and `setup` are used for the project initialization and compilation of `hymodelc`.

```
exe/
  hymodelc
  ...
fortran/
  ...
out/
  yyyymmddHH_lati_long_zagl
    yyyymmddHH_lati_long_zagl_traj.rds
    yyyymmddHH_lati_long_zagl_foot.nc
    hymodelc
    SETUP.CFG
    CONTROL
    ...
  ...
r/
  src/
    ...
  dependencies.r
  run_stilt.r
setup
```

## exe/

Location for the compiled `hymodelc` executable as well as shared model configuration files, such as `ASCDATA.CFG`, `CONC.CFG`, `LANDUSE.ASC`, `ROUGLEN.ASC`, etc. These files are simlinked to each simulation (`out/`) directory.

## fortran/

`hymodelc` fortran source code and make file used for compilation. Compilation is most easily done with the included `setup` script (which runs as part of the standard `uataq::stilt_init()`{r} installation) or manually with

```
cd fortran
make
mv hymodelc ../exe/hymodelc
cd ..
chmod +x exe/hymodelc
```

`permute.f90` contains a spatial permutation subroutine that applies the particle influences using variable bandwidth gaussian kernels. This is designed to be compiled as an R-compatible dynamic link library, which allows it to be used interactively within the R environment. Compilation is most easily done with the included `setup` script (which runs as part of the standard `uataq::stilt_init()` installation) or manually with

```
cd fortran
R CMD SHLIB permute.f90
mv permute.so ../r/src/permute.so
```

## out/

Initially empty, this folder propagates subdirectories with the naming convention `yyyymmddHH_lati_long_zagl` when simulations begin. See Model output.

## r/

`run_stilt.r` is used to adjust model parameters, execute the simulations, and produce upstream influence footprints. See Workflow and run_stilt.r.

`dependencies.r` is used to load the necessary functions and packages for `run_stilt.r` and is loaded on each forked process during parallelization.

The `src/` subdirectory contains the bulk of the R source code but should not need to be modified. Each file contains a single R function with metadata documenting function arguments and usage instructions for making programatic adjustments to STILT's workflow.


# Workflow and run_stilt.r

After turning the desired knobs in the sections below, the `r/run_stilt.r` script symlinks the meteorological data path to the user's home directory, with the default format `paste0('m', project)`. This is done to avoid issues with paths longer than 80 characters, which will result in unsuccessful fortran simulations.

User parameters are then passed to `stilt_apply`, which chooses the appropriate parallel or serial function to dispatch the simulations. If using SLURM for job submission, `stilt_apply` will use the `rslurm` package to submit jobs across `n_nodes` and `n_cores` per node. If running in parallel on a single node without SLURM, `stilt_apply` will use the `parallel` package to run simulations on the current node across `n_cores`. Otherwise, `stilt_apply` will run the simulations serially using `lapply()`.

## Computational arguments

Arg       | Description
----------|-------------------------------------------------------------------------
`project` | Project name
`stilt_wd`| Working directory for the given project
`lib.loc` | Path to R package installations, passed to `library()`
`rm_dat`  | Logical indicating whether to delete `PARTICLE.DAT` after each simulation. Default to TRUE to reduce disk space since all of the trajectory information is also stored in `STILT_OUTPUT.rds` alongside the calculated upstream influence footprint
`n_nodes` | If using SLURM for job submission, number of nodes to utilize
`n_cores` | Number of cores per node to parallelize simulations by receptor locations and times
`slurm`   | Logical indicating the use of rSLURM to submit job(s)
`slurm_options` | Named list of options passed to `sbatch` using `rslurm::slurm_apply`, which typically includes `time`, `account`, and `partition` values

## Simulation timing

Arg               | Description
------------------|-------------------------------------------------------------------------
`t_start/t_end`   | Simulation timing, formatted as `'yyyy-mm-dd HH:MM:SS'` UTC
`run_times`       | Hourly simulations spanning `t_start` through `t_end`

## Receptor locations

Arg               | Description
------------------|-------------------------------------------------------------------------
`lati`            | Receptor latitude(s), in degrees. Can be a single value or a vector with length equal to `length(run_times)` of receptor latitudes
`long`            | Receptor longitude(s), in degrees. Can be a single value or a vector with length equal to `length(run_times)` of receptor longitudes
`zagl`            | Receptor height(s), in meters above ground level. Can be a single value or a vector with length equal to `length(run_times)` of receptor heights

## Meteorological data input

Arg               | Description
------------------|-------------------------------------------------------------------------
`met_directory`   | Full directory path in which ARL compatible meteorological data files can be found
`met_file_format` | `strftime` compatible file naming convention to identify meteorological data files necessary for the timing of the simulation, such as `%Y%m%d.%Hz.hrrra`. Also accepts wildcards that match regular expression patterns passed to `dir()`, such as `wrfout_d0*_jul.arl`
`n_met_min`       | Minimum number of meteorological data files with which to proceed with simulation. Useful for capturing missing periods. For a -24 hour simulation using the 6 hour HRRR met data files, `n_met_min` should be set to 5, since `find_met_files()` ensures that data before and after the simulation is included.

## Model control

Arg             | Description
----------------|-------------------------------------------------------------------------
`delt` | integration timestep [min]; if set to 0.0, then timestep is dynamically determined
`n_hours` | Number of hours to run each simulation; negative indicates backward in time
`numpar` | number of particles to be run; defaults to 200
`outdt` | interval [min] to output data to PARTICLE.DAT; defaults to 0.0, which outputs at every timestep
`run_trajec` | Logical indicating whether to produce new trajectories with `hymodelc`. If FALSE, will load the previous `STILT_OUTPUT.rds` for regridding purposes
`timeout` | number of seconds to allow hymodelc to complete before sending SIGTERM and moving to the next simulation; defaults to 3600 (1 hour)
`varsiwant` | character vector of 4-letter hymodelc variables. Options include `NULL` for all variables, or a vector containing a minimum of `c('time', 'indx', 'lati', 'long', 'zagl', 'foot')` and optionally containing 'sigw', 'tlgr', 'zsfc', 'icdx', 'temp', 'samt', 'shtf', 'tcld', 'dmas', 'dens', 'rhfr', 'sphu', 'solw', 'lcld', 'zloc', 'dswf', 'wout', 'mlht', 'rain', 'crai'

## Transport and dispersion

Arg             | Description
----------------|-------------------------------------------------------------------------
`iconvect` | flag for convection. If set to 1, then runs excessive convection as described in Gerbig et al., 2003. For specialized RAMS output, the particles will be vertically redistributed according to the output convective mass fluxes; defaults to 0
`isot` | flag used to set the isotropic turbulence option; defaults to 0 to compute horizontal turbulence from wind field deformation. Setting to 1 results in the horizontal turbulence to be the same in both the u and v directions
`mgmin` | **Needs documentation**
`ndump` | flag to dump all particle/puff points at the end of a simulation to a file called PARDUMP. This can be read at the start of a new simulation to continue the previous calculation. Valid settings include 0 (no i/o), 1 (read/write), 2 (read only), 3 (write only); defaults to 0
`nturb` | no turbulence flag; defaults to 0, which includes turbulence rather than simulating mean trajectories
`outfrac` | the fraction of the particles that are allowed to leave the model domain (given by met data); defaults to 0.9. If exceeded, the model stops
`random` | flag that tells the random number generator whether to have a different random sequence for each model run (0 - false, 1 - true); defaults to 1
`tlfrac` | the fraction of the lagrangian timescale TL to set as timestep in dispersion subroutine. The smaller this fraction is, the more finely the turbulence is resolved; defaults to 0.1
`tratio` | maximum fraction of gridcell to be travelled by a particle in a single integration timestep. This determines the timestep if DELT is set to be dynamic
`veght` | height below which a particle's time spent is tallied; defaults to 0.5, which specifies half of the PBL. Setting <=1.0 specifies a fraction of the boundary layer height, and setting >1.0 specifies a height above ground in meters
`w_option` | vertical motion calculation method. 0: use vertical velocity from data, 1: isob, 2: isen, 3: dens, 4: sigma; defaults to 0
`winderrtf` | flag that specifies whether to have particle motions be affected by horizontal wind errors; defaults to 0. If set to 1, then STILT looks for a file called "WINDERR" that has four lines: 1. Standard deviation of errors [m/s], 2. Correlation timescale of errors [min], 3. Vertical correlation lengthscale of errors [m], 4. Horizontal correlation lengthscale of errors [km]. All statistical properties are applied equally in the u and v wind components
`z_top` | top of model domain, in meters above ground level; defaults to 25000.0
`zicontroltf` | flag that specifies whether to scale the PBL heights in STILT uniformly in the entire model domain; defaults to 0. If set to 1, then STILT looks for a file called "ZICONTROL" that specifies the scaling for the PBL height. The first line indicates the number of hours that the PBL height will be changed, and each subsequent line indicates the scaling factor for that hour

## Footprint gridding

Arg             | Description
----------------|-------------------------------------------------------------------------
`xmn` | grid start longitude, in degrees from -180 to 180
`xmx` | grid end longitude, in degrees from -180 to 180
`ymn` | grid start latitude, in degrees from -180 to 180
`ymx` | grid end latitude, in degrees from -180 to 180
`xres` | resolution for longitude grid, in degrees
`yres` | resolution for latitude grid, in degrees
`smooth_factor` | factor by which to linearly scale footprint smoothing; 0 to disable, defaults to 1
`time_integrate` | logical indicating whether to integrate footprint over time or retain discrete time steps


# Model output

The `out/` directory contains model initialization files and outputs in subdirectories organized by receptor timing and locations. The directories are named with the convention `yyyymmddHH_long_lati_zagl` where

Abbreviation   | Value
---------------|----------------------------------
`yyyy`         | Year (start)
`mm`           | Month (start)
`dd`           | Day (start)
`HH`           | Hour (start)
`lati`         | Receptor latitude (deg)
`long`         | Receptor longitude (deg)
`zagl`         | Receptor height above ground (m)

These directories also contain symbolic links to the files in `exe/` including the `hymodelc` executable. When the model runs trajectory simulations, diagnostic information (`hymodelc.out`) and model output are also stored here.

## Trajectory .rds

Particle trajectories and run information are packaged and saved in a compressed .rds file with the naming convention with the naming convention `YYYYMMDDHH_LONG_LATI_ZAGL_traj.rds`. This allows regridding of the footprints without recalculating particle trajectories. Load into R with `out <- readRDS('YYYYMMDDHH_LONG_LATI_ZAGL_traj.rds')`{.r}.

This object is structured as

```{r, eval = F}
out <- readRDS('YYYYMMDDHH_LONG_LATI_ZAGL_foot.rds')
str(out)
# List of 4
#  $ runtime : POSIXct[1:1], format: "2015-06-16"
#  $ file    : chr "/uufs/chpc.utah.edu/common/home/u0791983/stilt-sims/test/out/2015061600_-111.847672_40.766189_10/2015061600_-111.847672_40.7661"| __truncated__
#  $ receptor:Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	1 obs. of  4 variables:
#   ..$ run_time: POSIXct[1:1], format: "2015-06-16"
#   ..$ lati    : num 40.8
#   ..$ long    : num -112
#   ..$ zagl    : num 10
#  $ particle:Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	128441 obs. of  26 variables:
#   ..$ time: num [1:128441] -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 ...
#   ..$ indx: num [1:128441] 1 2 3 4 5 6 7 8 9 10 ...
#   ..$ long: num [1:128441] -112 -112 -112 -112 -112 ...
#   ..$ lati: num [1:128441] 40.8 40.8 40.8 40.8 40.8 ...
#   ..$ zagl: num [1:128441] 61.3 100.6 87.3 98.3 88.6 ...
#   ..$ sigw: num [1:128441] 1.15 1.15 1.15 1.15 1.15 ...
#   ..$ tlgr: num [1:128441] 8.05 8.05 8.05 8.05 8.05 ...
#   ..$ zsfc: num [1:128441] 1532 1532 1532 1532 1532 ...
#   ..$ icdx: num [1:128441] 2 2 2 2 2 2 2 2 2 2 ...
#   ..$ temp: num [1:128441] 301 301 301 301 301 ...
#   ..$ samt: num [1:128441] 2 2 2 2 2 ...
#   ..$ foot: num [1:128441] 0.0219 0.0219 0.0219 0.0219 0.0219 ...
#   ..$ shtf: num [1:128441] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ tcld: num [1:128441] 87.9 87.9 87.9 87.9 87.9 ...
#   ..$ dmas: num [1:128441] 1.104 1.104 1.104 0.991 1.119 ...
#   ..$ dens: num [1:128441] 0.968 0.965 0.966 0.965 0.966 ...
#   ..$ rhfr: num [1:128441] 0.277 0.282 0.28 0.282 0.28 ...
#   ..$ sphu: num [1:128441] 0.00222 0.00222 0.00222 0.00222 0.00222 ...
#   ..$ solw: num [1:128441] -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 ...
#   ..$ lcld: num [1:128441] -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 ...
#   ..$ zloc: num [1:128441] -999 -999 -999 -999 -999 -999 -999 -999 -999 -999 ...
#   ..$ dswf: num [1:128441] 25.2 25.2 25.2 25.2 25.2 25.2 25.2 25.2 25.2 25.2 ...
#   ..$ wout: num [1:128441] 0.722 0.722 0.722 0.722 0.722 ...
#   ..$ mlht: num [1:128441] 329 329 329 329 329 ...
#   ..$ rain: num [1:128441] 4.02e-07 4.02e-07 4.02e-07 4.02e-07 4.02e-07 ...
#   ..$ crai: num [1:128441] -0.9 -0.9 -0.9 -0.9 -0.9 -0.9 -0.9 -0.9 -0.9 -0.9 ...
```

Particle trajectory data is stored in a data frame with columns corresponding with varsiwant and can be accessed with `out$particle`{.r}.

## Footprint .nc

Footprints are packaged and saved in a compressed .nc file with the naming convention `YYYYMMDDHH_LONG_LATI_ZAGL_foot.nc`. This object contains information about the model domain, the grid resolution, and footprints per timestep (if `!time_integrate`). For those familiar with raster operations, the default output is compatible with the `raster` package. If `time_integrate`, the data can be loaded with `raster()`. If `!time_integrate`, the data can be loaded with `brick()`, which contains the POSIX time (UTC seconds since 1970-01-01) in the Z dimension (accessed with `getZ()`).

For more information about raster manipulation, the [Raster R package](https://geoscripting-wur.github.io/IntroToRaster/) is a good place to start.
