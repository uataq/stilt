#!/bin/bash
# Ben Fasoli
# Integration testing for STILT R wrapper
set -e

# Fetch the tutorial data
[[ -d stilt-tutorials ]] || git clone https://github.com/uataq/stilt-tutorials

chmod +x r/stilt_cli.r

echo "Running r/stilt_cli.r"
r/stilt_cli.r \
  r_run_time=2015-12-10T00:00:00Z \
  r_lati=40.5,40.51 \
  r_long=-112.0,-112.01 \
  r_zagl=0,1000 \
  met_path=$(pwd)/stilt-tutorials/01-wbb/met \
  met_file_format=%Y%m%d.%Hz.hrrra \
  n_hours=-12 \
  xmn=-113 \
  xmx=-111 \
  xres=0.01 \
  ymn=39.5 \
  ymx=41.5 \
  yres=0.01 \
  simulation_id=test-slant-column

# Check output
model_output=$(ls out/by-id/test-slant-column/test-slant-column* | wc -l)
if [ $model_output -lt 2 ]; then
  echo "Model output not found."

  echo "stilt.log:"
  cat out/by-id/test-slant-column/stilt.log
  exit 1
fi

echo "out/by-id/<id> contents:"
ls -lh out/by-id/test-slant-column

echo "Removing model outputs"
rm out/by-id/test-slant-column/*
rm out/footprints/*
rm out/particles/*

echo "stilt_cli.r test successful"
