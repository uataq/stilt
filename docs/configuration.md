## STILT configuration

The following parameters are found in `r/run_stilt.r` and are used to configure STILT. These settings are used to adjust model parameters, execute parallelized simulations, and calculate produce upstream influence footprints.

### System configuration

| Arg         | Description                                                                                     |
| ----------- | ----------------------------------------------------------------------------------------------- |
| `project`   | Project name. Defaults to the name of the directory specified in `uataq::stilt_init()`          |
| `stilt_wd`  | Root directory of the STILT project. Defaults to the directory created by `uataq::stilt_init()` |
| `output_wd` | Directory containing simulation output files. Defaults to `<stilt_wd>/out/`                     |
| `lib.loc`   | Path to installed R packages, passed to `library()`                                             |

### Parallel simulation settings

| Arg             | Description                                                                                                                                                                               |
| --------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `n_nodes`       | If using SLURM for job submission, number of nodes to utilize                                                                                                                             |
| `n_cores`       | Number of cores per node to parallelize simulations by receptor locations and times                                                                                                       |
| `slurm`         | Logical indicating the use of rSLURM to submit job(s). When using SLURM, a `<stilt_wd>/_rslurm` directory is created to contain the SLURM submission scripts and node-specific log files. |
| `slurm_options` | Named list of options passed to `sbatch` using `rslurm::slurm_apply()`. This typically includes `time`, `account`, and `partition` values                                                 |

### Receptor placement

| Arg             | Description                                                                                  |
| --------------- | -------------------------------------------------------------------------------------------- |
| `t_start/t_end` | Receptor time(s) to initialize simulations, formatted as `'yyyy-mm-dd HH:MM:SS'` UTC         |
| `run_times`     | Hourly time increments spanning `t_start` through `t_end` of length _n_                      |
| `lati`          | Receptor latitude(s), in degrees as a single value or a vector of length _n_                 |
| `long`          | Receptor longitude(s), in degrees as a single value or a vector of length _n_                |
| `zagl`          | Receptor height(s), in meters above ground level as a single value or a vector of length _n_ |
| `receptors`     | Expands the combinations of time and location to generate the unqiue receptors               |

> Simulation timing and receptor locations are defined in this way for convenience and then expanded to contain the unique receptors in a _x_, _y_, _z_, _t_ table. To specify the receptors manually, a data frame named `receptors` can be given with column names `run_time` (POSIXct), `long` (numeric), `lati` (numeric), and `zagl` (numeric). If supplying `receptors` directly, the `t_start/t_end`, `run_times`, `lati`, `long`, and `zagl` parameters can be omitted.

```r
str(receptors)
# 'data.frame':	100 obs. of  4 variables:
#   $ run_time: POSIXct, format: "2015-07-02 11:00:00" "2015-07-02 11:00:00" ...
#   $ long    : num  -112 -112 -112 -112 -112 ...
#   $ lati    : num  40.8 40.8 40.8 40.8 40.8 ...
#   $ zagl    : num  5 5 5 5 5 5 5 5 5 5 ...
```

### Footprint calculation methods

| Arg              | Description                                                                                                                                                                                                                                                                                                                                     |
| ---------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `hnf_plume`      | logical indicating whether to apply a vertical gaussian plume model to rescale the effective dilution depth for particles in the hyper near-field. This acts to scale up the influence of hyper-local fluxes on the receptor. If enabled, requires `varsiwant` to include a minimum of `dens, tlgr, sigw, foot, mlht, samt`. Default is enabled |
| `projection`     | [proj4](https://proj4.org/usage/quickstart.html) string defining the map projection of the footprint netCDF output. Defaults to `+proj=longlat`                                                                                                                                                                                                 |
| `smooth_factor`  | factor by which to linearly scale footprint smoothing. Defaults to 1.                                                                                                                                                                                                                                                                           |
| `time_integrate` | logical indicating whether to integrate footprint over time or retain discrete hourly time steps in footprint output                                                                                                                                                                                                                            |
| `xmn`            | grid start longitude, in degrees from -180 to 180                                                                                                                                                                                                                                                                                               |
| `xmx`            | grid end longitude, in degrees from -180 to 180                                                                                                                                                                                                                                                                                                 |
| `ymn`            | grid start latitude, in degrees from -180 to 180                                                                                                                                                                                                                                                                                                |
| `ymx`            | grid end latitude, in degrees from -180 to 180                                                                                                                                                                                                                                                                                                  |
| `xres`           | resolution for longitude grid, in projection units (degrees for lat/lon, meters for other)                                                                                                                                                                                                                                                      |
| `yres`           | resolution for latitude grid, in projection units (degrees for lat/lon, meters for other)                                                                                                                                                                                                                                                       |

### Meteorological data input

| Arg                  | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| -------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `met_path`           | Absolute path to ARL compatible meteorological data files                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| `met_file_format`    | String detailing file naming convention for meteorological data files using a mixture of datetime and regex syntax. The formatting string accepts `grep` compatible regular expressions (`.\*.arl`), `strftime` compatible datetime strings (`%Y%m%d%H`) or any combination of the two. Datetime syntax is expanded to all unique combinations required for the receptor and simulation duration and the intersection between the requested files and files available in `met_path` is determined with `grep`, allowing partial matching and compatible regular expressions to be used to identify the relevant data. Matching does not require the full format to be specified - e.g. `\*.arl`, `%Y`, `%Y%m%d`, `%Y%m%d_d0.*.arl` would all match with a file named `20180130_d01.arl`. |
| `met_subgrid_buffer` | Percent to extend footprint area for meteorological subdomain when using `met_subgrid_enable`. Defaults to 0.1 (10%)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| `met_subgrid_enable` | Enables extraction of spatial subdomains from files in `met_path` using HYSPLIT's `xtrct_grid` binary prior to executing simulations. If enabled, will create files in `<output_wd>/met/`. This can substantially accelerate simulation speed at the cost of increased disk usage. Defaults to disabled                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| `met_subgrid_levels` | If set, extracts the defined number of vertical levels from the meteorological data files to further accelerate simulations. Defaults to `NA`, which includes all vertical levels available                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| `n_met_min`          | Require a minimum number of meteorological data files to be matched by `met_file_format` for the simulation to proceed. Useful for handling periods where meteorological data may be missing. For a -24 hour simulation using the 6 hour HRRR met data files, `n_met_min` should be set to 5. Defaults to 1.                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |

> NOAA publishes High Resolution Rapid Refresh (HRRR) mesoscale model data in the ARL packed format required for STILT at [ftp://arlftp.arlhq.noaa.gov/pub/archives/hrrr/](ftp://arlftp.arlhq.noaa.gov/pub/archives/hrrr/). This is often the easiest place to start but is only available after June 15, 2015. The coupling of the popular Weather Research and Forecasting (WRF) model with STILT is well documented by [Nehrkorn, 2010](https://link.springer.com/article/10.1007%2Fs00703-010-0068-x). You can access various compatible gridded meteorological data products at [https://www.ready.noaa.gov/archives.php](https://www.ready.noaa.gov/archives.php).

### Model control

| Arg             | Description                                                                                                                                                                                                                                        |
| --------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `n_hours`       | Number of hours to run each simulation; negative indicates backward in time                                                                                                                                                                        |
| `numpar`        | number of particles to be run; defaults to 200                                                                                                                                                                                                     |
| `rm_dat`        | Logical indicating whether to delete `PARTICLE.DAT` after each simulation. Default to TRUE to reduce disk space since all of the trajectory information is also stored in `STILT_OUTPUT.rds` alongside the calculated upstream influence footprint |
| `run_foot`      | Logical indicating whether to produce footprints. If FALSE, `run_trajec` must be TRUE. This can be useful when calculating trajectories separate from footprints                                                                                   |
| `run_trajec`    | Logical indicating whether to produce new trajectories with `hycs_std`. If FALSE, will try to load the previous trajectory outputs. This is often useful for regridding purposes                                                                   |
| `simulation_id` | Unique identifier for each simulation; defaults to NA which determines a unique identifier for each simulation by hashing the time and receptor location                                                                                           |
| `timeout`       | number of seconds to allow `hycs_std` to complete before sending SIGTERM and moving to the next simulation; defaults to 3600 (1 hour)                                                                                                              |
| `varsiwant`     | character vector of 4-letter `hycs_std` variables. Defaults to the minimum required variables including `'time', 'indx', 'long', 'lati', 'zagl', 'foot', 'mlht', 'dens', 'samt', 'sigw', 'tlgr'`. Can optionally include options listed below.     |

#### Optional `varsiwant` arguments

- `crai` convective rainfall rate [m/min]
- `dens` air density [kg/m3]
- `dmas` particle weight changes due to mass violation in wind fields [initial value = 1.0]
- `dswf` downward shortwave radiation [W/m2]
- `foot` footprint, or sensitivity of mixing ratio to surface fluxes [ppm/(μmole/m2/s)]
- `icdx` cloud index when using RAMS (Grell scheme) [1=updraft,2=environment,3=downdraft]
- `indx` unique particle identifier
- `lati` latitude position of particle [degrees]
- `lcld` low cloud cover [%]
- `lhtf` latent heat flux [W/m2]
- `long` longitude position of particle [degrees]
- `mlht` mixed-layer height [m]
- `rain` total rainfall rate [m/min]
- `rhfr` relative humidity fraction [0~1.0]
- `samt` amount of time particle spends below VEGHT (see section on SETUP.CFG) [min]
- `shtf` sensible heat flux [W/m2]
- `sigw` standard deviation of vertical velocity; measure of strength of vertical turbulence [m/s]
- `solw` soil moisture
- `sphu` specific humidity [g/g]
- `tcld` total cloud cover [%]
- `temp` air temperature at lowest model layer [K]
- `time` time since start of simulation; negative if going backward in time [min] indx particle index
- `tlgr` Lagrangian decorrelation timescale [s]
- `wout` vertical mean wind [m/s]
- `zagl` vertical position of particle [m above ground level]
- `zloc` limit of convection heights [m]
- `zsfc` terrain height [m above sea level]

### Transport and dispersion

| Arg           | Description                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| ------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `w_option`    | vertical motion calculation method. 0: use vertical velocity from data, 1: isob, 2: isen, 3: dens, 4: sigma; defaults to 0                                                                                                                                                                                                                                                                                                                   |
| `zicontroltf` | flag that specifies whether to scale the PBL heights in STILT uniformly in the entire model domain; defaults to 0. If set to 1, then STILT looks for a file called "ZICONTROL" that specifies the scaling for the PBL height. The first line indicates the number of hours that the PBL height will be changed, and each subsequent line indicates the scaling factor for that hour                                                          |
| `ziscale`     | manually scale the mixed-layer height, with each element specifying a scaling factor for each simulation hour (ziscale can be of length that is smaller than abs(nhrs). A vector can be passed as a list (e.g. `ziscale <- list(rep(0.8, 24))` scales the mixed layer height to 80% for the first 24 hours of all simulations) or a list of vectors specific to each simulation (e.g. `ziscale <- rep(list(rep(0.8, 24)), nrow(receptors))`) |
| `z_top`       | top of model domain, in meters above ground level; defaults to 25000.0                                                                                                                                                                                                                                                                                                                                                                       |

> Additional arguments can be referenced in the [HYSPLIT user's guide](https://www.arl.noaa.gov/documents/reports/hysplit_user_guide.pdf)

### Transport error calculations

| Arg           | Description                                                          |
| ------------- | -------------------------------------------------------------------- |
| `siguverr`    | standard deviation of horizontal wind errors [m/s]                   |
| `tluverr`     | standard deviation of horizontal wind error timescale [min]          |
| `zcoruverr`   | vertical correlation lengthscale [m]                                 |
| `horcoruverr` | horizontal correlation lengthscale [km]                              |
| `sigzierr`    | standard deviation of mixed layer height errors [%]                  |
| `tlzierr`     | standard deviation of mixed layer height timescale [min]             |
| `horcorzierr` | horizontal correlation lengthscale of mixed layer height errors [km] |

### Inject user defined functions

This is an advanced option that can be used to inject user code into the simulation flow. This can be used to add custom logic at key points in the code, such as correcting for [mass violation in the input wind fields](https://github.com/uataq/stilt/issues/41#issuecomment-656174491).

| Arg                | Description                                                                                                                                                   |
| ------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `before_trajec`    | function that returns the `output` object and is executed prior to calculating the ensemble particle trajectories                                             |
| `before_footprint` | function that returns the `output` object and is executed prior to calculating the gridded fotprint but after calculating the ensmble's particle trajectories |

---

## Next steps

- [Execution](execution.md) for details on how to run your STILT simulations
- [Tutorial: Stationary simulations](https://github.com/uataq/stilt-tutorials/tree/main/01-wbb)
