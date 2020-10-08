# STILT Dependency Loader
# For documentation, see https://github.com/uataq/stilt
# Ben Fasoli

if (!'stilt_wd' %in% ls())
  stilt_wd <- getwd()

# Source r/src R scripts
rsc <- dir(file.path(stilt_wd, 'r', 'src'), pattern = '.*\\.r$', full.names = T)
invisible(lapply(rsc, source))

# Load external libraries
if (!'lib.loc' %in% ls()) lib.loc <- NULL
libs <- load_libs('dplyr',
                  'ncdf4',
                  'parallel',
                  'raster',
                  'rslurm',
                  lib.loc = lib.loc)

# Load permute fortran dll for footprint matrix permutation
permute_exe <- file.path(stilt_wd, 'r/src/permute.so')
if (!file.exists(permute_exe))
  stop('calc_footprint(): failed to find permute.so in r/src/')
dyn.load(permute_exe)

# Validate arguments and load dependencies if necessary
if ((!class(projection) == 'function') && ('projection' %in% ls()))
  validate_projection(projection)
if ('varsiwant' %in% ls())
  validate_varsiwant(varsiwant)
if (all(c('xmn', 'xmx', 'ymn', 'ymx') %in% ls()))
  validate_footprint_extent(xmn, xmx, ymn, ymx)

# Disable grouping message from dplyr >=1.0.0
options(dplyr.summarise.inform = F)
