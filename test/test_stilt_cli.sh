#!/bin/bash
# Ben Fasoli
# Integration testing for STILT R wrapper
set -e

chmod +x r/stilt_cli.r

echo "Running r/stilt_cli.r"
r/stilt_cli.r \
  r_run_time=2015-12-10T00:00:00Z \
  r_lati=40.5 \
  r_long=-112.0 \
  r_zagl=5 \
  met_loc=$(pwd)/stilt-tutorials/01-wbb/met \
  met_file_format=%Y%m%d.%Hz.hrrra \
  n_hours=-12 \
  xmn=-113 \
  xmx=-111 \
  xres=0.01 \
  ymn=39.5 \
  ymx=41.5 \
  yres=0.01

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

echo "Model test successful"
