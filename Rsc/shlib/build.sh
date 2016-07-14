#! /bin/bash
#
#  Builds the shared libraries of the stiltR part
#
# at Jena : 
# -  login at pc027
# -  source /usr/local/etc/sources/Intel_fort.sh

NETCDFLIB=/Net/Groups/BSY/tpobw/STILT/IER/netCDF/netcdf-3.6.2
source /usr/local/etc/sources/Intel_fort.sh

SOURCES=(getfossEUnetCDF.f90 getTM3bin.f90 getTM3CH4bin.f90 getvegfracNetCDF.f90 getMODISnetCDF.f90 getEmisnetCDF.f90 )
#SOURCES=(getMODISnetCDF.f90)
#SOURCES=(getMODISnetCDF_dyn.f90)
#SOURCES=(getfossBARCAnetCDF.f90)
#SOURCES=(getfossO2BARCAnetCDF.f90)
#SOURCES=(getfireBARCAnetCDF.f90)
#SOURCES=(getMODISnetCDF_dynV1.f90)
for M in ${SOURCES[@]}; do
   ifort -fPIC -shared -assume byterecl -warn all,nodec,interfaces,unused -ftrapuv -check all -nomixed-str-len-arg\
   -traceback -O2 -fp-stack-check -I$NETCDFLIB/f90 \
   ${M} -o ${M/.f90/.so} \
   $NETCDFLIB/fortran/.libs/libnetcdff.so $NETCDFLIB/libsrc/.libs/libnetcdf.so
   echo $M done $'=====\n'
done
