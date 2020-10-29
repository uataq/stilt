## STILT project structure

The STILT project framework contains an initial scaffolding with a number of directories for model configuration, source code, binaries, and outputs. Initially, only the `r/` and `fortran/` directories are fully populated.

```
exe/
  hycs_std
  xtrct_grid
  ...
out/
  by-id/
    yyyymmddHH_lati_long_zagl/
      yyyymmddHH_lati_long_zagl_traj.rds
      yyyymmddHH_lati_long_zagl_foot.nc
      hycs_std
      SETUP.CFG
      CONTROL
      ...
    ...
  footprints/
    yyyymmddHH_lati_long_zagl_foot.nc
    ...
  particles/
    yyyymmddHH_lati_long_zagl_traj.rds
    ...
r/
  src/
    ...
  dependencies.r
  run_stilt.r
```

### exe/

Files here are shared across all model runs. Each file stored within exe/ is symbolically linked to each simulation directory in out/by-id/.

This is where you will find the compiled `hycs_std` and `xtrct_grid` executables as well as global model configuration files such as `ASCDATA.CFG`, `CONC.CFG`, `LANDUSE.ASC`, and `ROUGLEN.ASC`.

### out/

Initially nonexistant and configurable with the `output_wd` parameter, this folder contains subdirectories containing simulation information and outputs. These are organized into the following three subdirectories for convenience.

#### out/by-id/

Contains simulation files by simulation id, with the naming convention `yyyymmddHHMM_lati_long_zagl`.

| Abbreviation | Value                            |
| ------------ | -------------------------------- |
| yyyy         | Year (start)                     |
| mm           | Month (start)                    |
| dd           | Day (start)                      |
| HH           | Hour (start)                     |
| MM           | Minute (start)                   |
| lati         | Receptor latitude (deg)          |
| long         | Receptor longitude (deg)         |
| zagl         | Receptor height above ground (m) |

This becomes the working directory for each unique simulation, containing symbolic links to all of the shared files in `exe/` as well as simulation specific `CONTROL`, `SETUP.CFG`, and output files.

STILT outputs two files for analysis. The trajectories of the particle ensemble are saved to a `<simulation_id>_traj.rds` file. Gridded footprints are saved to a `<simulation_id>_foot.nc` file. For guidance on working with these output files, see [output files](output-files.md).

#### out/footprints/

Contains symbolic links to footprint files found in `out/by-id/<simulation_id>` directories for convenience.

#### out/particles/

Contains symbolic links to trajectory files found in `out/by-id/<simulation_id>` directories for convenience.

#### out/met/

Optional directory used when `met_subgrid_enable == TRUE` which contains extracted subdomains from the input meteorological data files.

### r/

Contains configuration data and source code.

`run_stilt.r` is the primary script users will interact with. It contains settings used to adjust model parameters, execute parallelized simulations, and calculate produce upstream influence footprints. These parameters are documented in [configuration](configuration.md).

`dependencies.r` is used to install and load the necessary functions on each forked parallel process.

#### r/src/

Contains the bulk of the source code for the control layer and footprint gridding algorithm. **The source code found here will not need to be modified by the majority of users**. Each file contains a single R function with metadata documenting function arguments and usage instructions for making programatic adjustments to STILT’s workflow.

---

## Next steps

- [Best practices](best-practices.md) for your STILT workflow
- [Tutorial: Stationary simulations](https://github.com/uataq/stilt-tutorials/tree/main/01-wbb)
