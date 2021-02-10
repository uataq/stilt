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
  r_lati=40.5 \
  r_long=-112.0 \
  r_zagl=5 \
  met_path=$(pwd)/stilt-tutorials/01-wbb/met \
  met_file_format=%Y%m%d.%Hz.hrrra \
  n_hours=-12 \
  xmn=-113 \
  xmx=-111 \
  xres=0.01 \
  ymn=39.5 \
  ymx=41.5 \
  yres=0.01

# Check output
model_output=$(ls out/by-id/201512100000_-112_40.5_5/201512100000_-112_40.5_5* | wc -l)
if [ $model_output -lt 2 ]; then
  echo "Model output not found."

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

echo "stilt_cli.r test successful"
