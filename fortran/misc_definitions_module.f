!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODULE MISC_DEFINITIONS_MODULE
!
! This module defines various non-meteorological constants that are used
!   by other modules for readability.
!
!-------------------------------------------------------------------------------
!  $Id: misc_definitions_module.f,v 1.1 2009/10/26 15:36:54 jel Exp $
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module misc_definitions_module

   real, parameter :: NAN=1.E20

   real, parameter :: NOT_MASKED   = -2.,  &
                      MASKED_BOTH  = -1.,  &
                      MASKED_WATER =  0.,  &
                      MASKED_LAND  =  1.

   integer, parameter :: OUTSIDE_DOMAIN=1E8, NOT_PROCESSED=1E9, INVALID=1E9

   integer, parameter :: SIXTEEN_POINT=1, FOUR_POINT=2, N_NEIGHBOR=3, &
                         AVERAGE4=4, AVERAGE16=5, W_AVERAGE4=6, W_AVERAGE16=7, &
                         SEARCH=8

   integer, parameter :: BOTTOM_TOP=1, TOP_BOTTOM=2

   integer, parameter :: CONTINUOUS=0, CATEGORICAL=1, SP_CONTINUOUS=2

   integer, parameter :: M=1, U=2, V=3, HH=4, VV=5

   integer, parameter :: ONETWOONE=1, SMTHDESMTH=2, SMTHDESMTH_SPECIAL=3

   integer, parameter :: BINARY=1, NETCDF=2, GRIB1=3, HDF=4

   ! Projection codes for proj_info structure:
   INTEGER, PUBLIC, PARAMETER  :: PROJ_LATLON = 0
   INTEGER, PUBLIC, PARAMETER  :: PROJ_LC = 1
   INTEGER, PUBLIC, PARAMETER  :: PROJ_PS = 2
   INTEGER, PUBLIC, PARAMETER  :: PROJ_PS_WGS84 = 102
   INTEGER, PUBLIC, PARAMETER  :: PROJ_MERC = 3
   INTEGER, PUBLIC, PARAMETER  :: PROJ_GAUSS = 4
   INTEGER, PUBLIC, PARAMETER  :: PROJ_ROTLL = 203

end module misc_definitions_module
