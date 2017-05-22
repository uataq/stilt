# STILT Dependency Loader
# For documentation, see https://github.com/benfasoli/stilt
# Ben Fasoli

if (!'stilt_wd' %in% ls())
  stilt_wd <- getwd()

rsc <- dir(file.path(stilt_wd, 'r', 'src'), pattern = '.*\\.r$', full.names = T)
invisible(lapply(rsc, source))


if (!'lib.loc' %in% ls())
  lib.loc <- NULL
load_libs('dplyr', 'ncdf4', 'parallel', 'raster', 'rslurm', 'uataq',
          lib.loc = lib.loc)


permute_exe <- file.path(stilt_wd, 'r/src/permute.so')
if (!file.exists(permute_exe))
  stop('calc_footprint(): failed to find permute.so in r/src/')
dyn.load(permute_exe)
