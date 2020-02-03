#!/bin/bash
# Ben Fasoli
# Integration testing for STILT R wrapper
set -e

# Fetch the tutorial data
git clone https://github.com/uataq/stilt-tutorials

# Replace {{project}} and {{wd}}
sed -i 's|{{project}}|stilt-test|g' r/run_stilt.r
sed -i "s|file.path('{{wd}}', project)|getwd()|g" r/run_stilt.r

# Set receptor and footprint information
sed -i 's|2015-06-18 22:00:00|2015-12-10 00:00:00|g' r/run_stilt.r
sed -i 's|xmn <- NA|xmn <- -113|g' r/run_stilt.r
sed -i 's|xmx <- NA|xmx <- -111|g' r/run_stilt.r
sed -i 's|ymn <- NA|ymn <- 39.5|g' r/run_stilt.r
sed -i 's|ymx <- NA|ymx <- 41.5|g' r/run_stilt.r

# Set met_directory
sed -i "s|'/uufs/chpc.utah.edu/common/home/lin-group6/hrrr/data/utah'|file.path(stilt_wd, 'stilt-tutorials/01-wbb/met')|g" r/run_stilt.r

# Minimize run duration
sed -i 's|n_hours    <- -24|n_hours <- -6|g' r/run_stilt.r

echo "Running r/run_stilt.r"
Rscript r/run_stilt.r

# Check output
model_output=$(ls out/by-id/2015121000_-112_40.5_5/2015121000_-112_40.5_5* | wc -l)
if [ $model_output -lt 2 ]; then
  echo "Model output not found."

  echo "run_stilt.r configuration:"
  cat r/run_stilt.r
  
  echo "log.txt:"
  cat out/by-id/2015121000_-112_40.5_5/log.txt
  exit 1
fi

echo "out/by-id/<id> contents:"
ls -lh out/by-id/2015121000_-112_40.5_5

echo "Removing model outputs"
rm out/by-id/2015121000_-112_40.5_5/*
rm out/footprints/*
rm out/particles/*
