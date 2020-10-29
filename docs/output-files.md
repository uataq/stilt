## STILT outputs

The model outputs can be found in the directory configured with `output_wd` (defaults to `<stilt_wd>/out/`, see [project structure](http://localhost:3000/#/project-structure)). STILT outputs two files for analysis -

- a `<simulation_id>_traj.rds` file containing the trajectories of the particle ensemble
- a `<simulation_id>_foot.nc` file containing gridded footprint values and metadata

Simulation identifiers follow a `yyyymmddHHMM_lati_long_zagl` convention, see [project structure](project-structure.md?id=outby-id).

## Particle trajectories

Particle trajectories and simulation configuration information are packaged and saved in a compressed `.rds` file ([serialized single R object](https://stat.ethz.ch/R-manual/R-devel/library/base/html/readRDS.html)) with the naming convention with the naming convention `<simulation_id>_traj.rds`. Preserving the particle trajectories enables regridding the footprints at a later time without the computational cost of recalculating particle trajectories.

This object can be loaded with `readRDS(<path>)` and is structured as

```r
traj <- readRDS('<simulation_id>_traj.rds')
str(traj)
# List of 4
# $ file    : chr "<stilt_wd>/out/by-id/<simulation_id>/<simulation_id>_traj.rds
# $ receptor:List of 4
# ..$ run_time: POSIXct[1:1], format: "1983-09-18 21:00:00"
# ..$ lati    : num 39.6
# ..$ long    : num -80.4
# ..$ zagl    : num 10
# $ particle:Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	11712136 obs. of  12 variables:
# ..$ time                : num [1:11712136] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ...
# ..$ indx                : num [1:11712136] 1 2 3 4 5 6 7 8 9 10 ...
# ..$ long                : num [1:11712136] -80.4 -80.4 -80.4 -80.4 -80.4 ...
# ..$ lati                : num [1:11712136] 39.6 39.6 39.6 39.6 39.6 ...
# ..$ zagl                : num [1:11712136] 5.63 41.37 30.95 3.32 28.88 ...
# ..$ foot                : num [1:11712136] 0.0626 0.043 0.043 0.0626 0.0627 ...
# ..$ mlht                : num [1:11712136] 1459 1459 1459 1459 1459 ...
# ..$ dens                : num [1:11712136] 1.1 1.1 1.1 1.1 1.1 ...
# ..$ samt                : num [1:11712136] 1 1 1 1 1 1 1 1 1 1 ...
# ..$ sigw                : num [1:11712136] 1.08 1.09 1.09 1.08 1.08 ...
# ..$ tlgr                : num [1:11712136] 1.67 4.98 4.98 1.67 1.67 ...
# ..$ foot_no_hnf_dilution: num [1:11712136] 0.00224 0.00224 0.00224 0.00224 0.00224 ...
```

The `traj$receptor` object is a named list with the time and location of the release point for the simulation. The `traj$particle` object is a data frame containing each particle's position and characteristics over time. 

## Gridded footprints

Footprints are packaged and saved in a compressed NetCDF file using [Climate and Forecast (CF)](http://cfconventions.org) compliant metadata with the naming convention `<simulation_id>_foot.nc`. This object contains information about the model domain, the grid resolution, and footprint values. This object is typically a three dimensional array with dimensions ordered (*x*, *y*, *t*). However, the object will only have dimensions (*x*, *y*) for time integrated footprints.

```bash
ncdump -h <simulation_id>_foot.nc

netcdf <simulation_id>_foot {
dimensions:
	lon = 1320 ;
	lat = 710 ;
	time = 6 ;
variables:
	double lon(lon) ;
		lon:units = "degrees_east" ;
		lon:standard_name = "longitude" ;
		lon:long_name = "longitude at cell center" ;
	double lat(lat) ;
		lat:units = "degrees_north" ;
		lat:standard_name = "latitude" ;
		lat:long_name = "latitude at cell center" ;
	double time(time) ;
		time:units = "seconds since 1970-01-01 00:00:00Z" ;
		time:standard_name = "time" ;
		time:long_name = "utc time" ;
		time:calendar = "standard" ;
	float foot(time, lat, lon) ;
		foot:units = "ppm (umol-1 m2 s)" ;
		foot:_FillValue = -1.f ;
		foot:standard_name = "footprint" ;
		foot:long_name = "stilt surface influence footprint" ;

// global attributes:
		:crs = "+proj=longlat" ;
		:crs_format = "PROJ.4" ;
		:documentation = "github.com/uataq/stilt" ;
		:title = "STILT Footprint" ;
		:time_created = "1992-04-15 00:00:00" ;
}

```

For those familiar with raster datasets, the footprint adheres to the [CF-1.4 metadata conventions](http://cfconventions.org/) which is inherently compatible compatible with the R `raster` package. For more information about raster manipulation, the [Raster R package](https://geoscripting-wur.github.io/IntroToRaster/) is a great place to start.


## Interfacing with Raster R package

The footprint NetCDF files can be read direactly into raster (2d) and brick (3d) objects in R using the `raster` package. The POSIX time (UTC seconds since 1970-01-01) is stored in the Z dimension and can be easily accessed with `getZ(footprint)`.

Using the [Raster R package](https://geoscripting-wur.github.io/IntroToRaster/) the data can be loaded with `brick()` and is structured as -

```r
library(raster)
footprint <- brick('<simulation_id>_foot.nc')
footprint
# class       : RasterBrick
# dimensions  : 710, 1320, 937200, 6  (nrow, ncol, ncell, nlayers)
# resolution  : 0.01, 0.01  (x, y)
# extent      : -84.3, -71.1, 39.6, 46.7  (xmin, xmax, ymin, ymax)
# coord. ref. : +proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0
# data source : <stilt_wd>/out/by-id/<simulation_id>/<simulation_id>_foot.nc
# names       : X432745200, X432748800, X432752400, X432756000, X432759600, X432763200
# z-value     : 432745200, 432748800, 432752400, 432756000, 432759600, 432763200
# varname     : foot

time <- as.POSIXct(getZ(footprint), tz = 'UTC', origin = '1970-01-01')
str(time)
# POSIXct[1:5], format: "2015-07-02 05:00:00" "2015-07-02 06:00:00" ...
```

#### Manual loading with ncdf4 R package

Alternatively the footprint data can be loaded using standard NetCDF methods.

```r
library(ncdf4)
nc <- nc_open('<simulation_id>_foot.nc')
nc
# File <simulation_id>_foot.nc (NC_FORMAT_NETCDF4):
# 
#      1 variables (excluding dimension variables):
#         float foot[lon,lat,time]   (Contiguous storage)  
#             units: ppm (umol-1 m2 s)
#             _FillValue: -1
#             standard_name: footprint
#             long_name: stilt surface influence footprint
# 
#      3 dimensions:
#         lon  Size:1320
#             units: degrees_east
#             standard_name: longitude
#             long_name: longitude at cell center
#         lat  Size:710
#             units: degrees_north
#             standard_name: latitude
#             long_name: latitude at cell center
#         time  Size:6
#             units: seconds since 1970-01-01 00:00:00Z
#             standard_name: time
#             long_name: utc time
#             calendar: standard
# 
#     5 global attributes:
#         crs: +proj=longlat
#         crs_format: PROJ.4
#         documentation: github.com/uataq/stilt
#         title: STILT Footprint
#         time_created: 1992-04-15 00:00:00
```

---

## Next steps

- [Tutorial: Stationary simulations](https://github.com/uataq/stilt-tutorials/tree/master/01-wbb)
