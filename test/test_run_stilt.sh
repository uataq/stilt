#!/bin/bash
# Ben Fasoli
# Integration testing for STILT R wrapper
set -e

# Fetch the tutorial data
[[ -d stilt-tutorials ]] || git clone https://github.com/uataq/stilt-tutorials

# Replace {{project}} and {{wd}}
sed -i'.bak' -e 's|{{project}}|stilt-test|g' r/run_stilt.r
sed -i'.bak' -e "s|file.path('{{wd}}', project)|getwd()|g" r/run_stilt.r

# Set receptor and footprint information
sed -i'.bak' -e 's|2015-06-18 22:00:00|2015-12-10 00:00:00|g' r/run_stilt.r
sed -i'.bak' -e 's|xmn <- NA|xmn <- -113|g' r/run_stilt.r
sed -i'.bak' -e 's|xmx <- NA|xmx <- -111|g' r/run_stilt.r
sed -i'.bak' -e 's|ymn <- NA|ymn <- 39.5|g' r/run_stilt.r
sed -i'.bak' -e 's|ymx <- NA|ymx <- 41.5|g' r/run_stilt.r

# Set met_path
sed -i'.bak' -e "s|'<path_to_arl_meteorological_data>'|file.path(stilt_wd, 'stilt-tutorials/01-wbb/met')|g" r/run_stilt.r

# Minimize run duration
sed -i'.bak' -e 's|n_hours    <- -24|n_hours <- -6|g' r/run_stilt.r

echo "Running r/run_stilt.r"
Rscript r/run_stilt.r

# Check output
model_output=$(ls out/by-id/201512100000_-112_40.5_5/201512100000_-112_40.5_5* | wc -l)
if [ $model_output -lt 2 ]; then
  echo "Model output not found."

  echo "run_stilt.r configuration:"
  cat r/run_stilt.r

  echo "stilt.log:"
  cat out/by-id/201512100000_-112_40.5_5/stilt.log
  exit 1
fi

echo "out/by-id/<id> contents:"
ls -lh out/by-id/201512100000_-112_40.5_5

echo "Removing model outputs"
rm out/by-id/201512100000_-112_40.5_5/*
rm out/footprints/*
rm out/particles/*
rm r/run_stilt.r.bak

echo "run_stilt.r test successful"
