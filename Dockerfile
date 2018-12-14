# STILT Dockerfile
# Ben Fasoli
#
# Builds docker image for dependency free single-shot cli STILT runs by
# supplying run_stilt.r parameters as command args
#
#   docker build -t stilt .
#
# The following mounts are required on call to docker run:
#
#   METDIR=/path/to/met/directory
#   OUTDIR=/path/to/out/directory
#   --mount type=bind,source=$METDIR,destination=/app/met,readonly \
#   --mount type=bind,source=$OUTDIR,destination=/app/out/by-id \
#
#
# Example
#
# Fetch example met data for testing
#   git clone https://github.com/uataq/stilt-tutorials /tmp/stilt-tutorials
# Create host input/output paths
#   METDIR=/tmp/stilt-tutorials/01-wbb/met
#   OUTDIR=/tmp/stilt-out
#   mkdir -p $METDIR $OUTDIR
#   docker run \
#     --rm \
#     --mount type=bind,source=$METDIR,destination=/app/met,readonly \
#     --mount type=bind,source=$OUTDIR,destination=/app/out/by-id \
#     stilt \
#     r_run_time=2015-12-10T00:00:00Z \
#     r_lati=40.5 \
#     r_long=-112.0 \
#     r_zagl=5 \
#     met_loc=/app/met \
#     met_file_format=%Y%m%d.%H \
#     xmn=-112.3 \
#     xmx=-111.52 \
#     xres=0.01 \
#     ymn=40.39 \
#     ymx=40.95 \
#     yres=0.01

FROM debian:stretch-slim

WORKDIR /app

COPY . /app

ENV TZ UTC
ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
                build-essential \
                git \
                libhdf5-serial-dev \
                libnetcdf-dev \
                libssl-dev \
                locales \
                netcdf-bin \
                procps \
                r-base \
                r-base-dev \
                unzip \
                wget \
        && locale-gen en_US.UTF-8 \
        && update-locale \
        && bash setup 3 \
        && Rscript r/dependencies.r \
        && apt-get remove --purge -y \
                build-essential \
                git \
                locales \
                wget \
        && apt-get autoremove -y \
        && rm -rf /var/lib/apt/lists/*

VOLUME ["/app/met", "/app/out"]

ENTRYPOINT ["/app/r/stilt_cli.r", \
                "met_loc=/app/met", \
                "output_wd=/app/out"]
