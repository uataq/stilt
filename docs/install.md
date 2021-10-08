## Install STILT

STILT has been compiled to run on UNIX platforms (Mac, Linux). Required software includes

- [R (version >= 3.5.0)](https://www.r-project.org/)
  - `dplyr` package for data manipulation
  - `parallel` package for single node parallelism
  - `rslurm` package for multi node parallelism
  - `raster` package for gridded spatial analysis tools
  - `uataq` package for STILT project initiation workflow
- [Git](https://git-scm.com/) for fetching STILT source code from GitHub
- [NetCDF (version >= 4.0)](https://www.unidata.ucar.edu/software/netcdf/) for storage libraries
- [GDAL](https://gdal.org) for geospatial transformations

## Fair use policy

STILT is freely available and we encourage others to use it. Kindly keep us informed of how you are using the model and of any publication plans. Please acknowledge the source as a citation. STILT is continuously updated and improved by the development consortium, and in some cases (as when new elements are used for the first time) we may suggest that one or more of the developers be offered participation as authors. If your work directly competes with our analysis that uses unpublished features, we may ask that we have the opportunity to submit a manuscript first. The software is updated from time to time and it is your responsibility to ensure that your publication is consistent with the most recent version.

## Installing with R

The preferred method uses R to initialize a new project and requires the `uataq` R package, which is open source and can be installed using `devtools`.

```r
if (!require('devtools')) install.packages('devtools')
devtools::install_github('uataq/uataq')
```

A STILT project can then be initialized in plain R code. A project name other than `stilt` should be chosen for projects to avoid naming conflicts with the repository.

```r
uataq::stilt_init('myproject')
# Cloning into 'myproject'...
# remote: Enumerating objects: 60, done.
# remote: Counting objects: 100% (60/60), done.
# remote: Compressing objects: 100% (56/56), done.
# remote: Total 60 (delta 3), reused 25 (delta 2), pack-reused 0
# Unpacking objects: 100% (60/60), 2.26 MiB | 4.09 MiB/s, done.
# Compiling footprint kernel aggregation subroutine...
#
# STILT installation successful.
#
# We strongly suggest you subscribe to the critical update notifications at
# https://uataq.github.io/stilt/
# to be notified if important STILT model updates updates.
```

This method sources statically compiled `hycs_std`, `xtrct_grid`, `xtrct_time`, and `arw2arl` binaries for 64-bit Intel and AMD systems. These binaries provide the same functionality as those compiled from source but their use with forecast data to model future scenarios has been disabled.

## Installing from source

> Requires one of the `gfortran` or `pgf90` Fortran compilers and NOAA ARL user registration

Compiling from source requires user [registration with NOAA ARL](https://www.ready.noaa.gov/HYSPLIT_register.php) to access the HYSPLIT source code. HYSPLIT version 5.0 or newer is required.

To compile the `hycs_std`, `xtrct_grid`, and `xtrct_time` binaries using `gfortran` (recommended) from the root directory of the HYSPLIT source code -

```bash
cp Makefile.inc.gfortran Makefile.inc
make library/libhysplit.a library/liblbfgsb.a
(cd exec && make hycs_std xtrct_grid xtrct_time)
```

> Versions 10.0 and newer of `gfortran` require modifying `Makefile.inc` to specify `FFLAGS=-fallow-argument-mismatch <other_fflags>`

If compilation is successful, you'll find binaries at `exec/` within the HYSPLIT source code. We now need to create a new STILT project, include these binaries, and compile the `permute.f90` DLL used by the footprint kernel algorithm.

```bash
git clone --depth=1 https://github.com/uataq/stilt myproject
# Cloning into 'myproject'...
# remote: Enumerating objects: 60, done.
# remote: Counting objects: 100% (60/60), done.
# remote: Compressing objects: 100% (56/56), done.
# remote: Total 60 (delta 3), reused 25 (delta 2), pack-reused 0
# Unpacking objects: 100% (60/60), 2.26 MiB | 3.37 MiB/s, done.

ls myproject
# Dockerfile  README.md  bin/  exe/  r/  setup  test/

cp hysplit/exec/hycs_std hysplit/exec/xtrct_grid hysplit/exec/xtrct_time myproject/exe/

R CMD SHLIB myproject/r/src/permute.f90
```

Edit the configuration in `myproject/r/run_stilt.r` and be sure to specify the project name and working directory.

---

## Next steps

- [Project structure](project-structure.md) for an overview of the files in your STILT project
- [Tutorial: Stationary simulations](https://github.com/uataq/stilt-tutorials/tree/main/01-wbb)
