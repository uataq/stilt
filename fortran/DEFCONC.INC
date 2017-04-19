!-------------------------------------------------------------------------------
! DEFCONC.INC 
! Defines concentration grid  structure, the pollutant identification, and
! release information.  Required as an INCLUDE for most CON routines.
! LAST REVISED: 29 Jan 1997 (RRD)
!               20 Sep 2000 (RRD) - fortran90 upgrade
!               27 Aug 2003 (RRD) - volume units conversion flag 
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! time identification 
!-------------------------------------------------------------------------------

  TYPE SSET
     SEQUENCE

     INTEGER    ::  yr      ! year
     INTEGER    ::  mo      ! month
     INTEGER    ::  da      ! day
     INTEGER    ::  hr      ! hour
     INTEGER    ::  mn      ! minute
     INTEGER    ::  ic      ! forecast hour
     INTEGER    ::  macc    ! accumulated minutes

  END TYPE

!-------------------------------------------------------------------------------
! define the geographic GRID parameters
!-------------------------------------------------------------------------------

  TYPE CSET
     SEQUENCE

     CHARACTER(80) :: dir           ! directory location
     CHARACTER(80) :: file          ! file name
     INTEGER       :: unit          ! associated unit number 
     INTEGER       :: numb_lat      ! number of latitude points
     INTEGER       :: numb_lon      ! number of longitude points
     REAL          :: delt_lat      ! latitude spacing in degrees
     REAL          :: delt_lon      ! longitude spacing in degrees
     REAL          :: size          ! grid size (km)
     REAL          :: x1y1_lat      ! lower left latitude
     REAL          :: x1y1_lon      ! lower left longitude
     INTEGER       :: snap          ! type of averaging (0:std 1:snapshot)
     INTEGER       :: levels        ! number of levels
     INTEGER       :: height(99)    ! heights of each level (m)
     TYPE(sset)    :: start         ! sample start time
     TYPE(sset)    :: stop          ! sample stop time
     TYPE(sset)    :: delta         ! sample integration interval
     TYPE(sset)    :: now           ! current time of sample

  END TYPE

!-------------------------------------------------------------------------------
! define the pollutant parameters
!-------------------------------------------------------------------------------

  TYPE PSET
     SEQUENCE

     CHARACTER(4)  :: ident   ! identifcation label
     REAL          :: qrate   ! emission rate in units per hour
     REAL          :: qhrs    ! hours of emission
     TYPE(sset)    :: start   ! emission starting time
     LOGICAL       :: dowet   ! wet removal flag
     LOGICAL       :: dodry   ! dry removal flag
     LOGICAL       :: dorad   ! radioactive decay flag
     LOGICAL       :: dores   ! resistance method flag
     LOGICAL       :: dogas   ! gas definition flag
     LOGICAL       :: dogrv   ! gravitational settling flag
     LOGICAL       :: dosus   ! resuspension flag
     LOGICAL       :: dovol   ! convert mass to volume units (emsgrd)
     REAL          :: wetin   ! wet removal ratio for within-cloud (by volume)
     REAL          :: wetlo   ! wet removal coefficient below-cloud (1/sec)
     REAL          :: wetgas  ! henry's consant for wet removal of gases
     REAL          :: henry   ! effective henry's constant for dry deposition
     REAL          :: dryvl   ! dry removal deposition velocity (m/s)
     REAL          :: pdens   ! particle density (g/m3)
     REAL          :: pdiam   ! particle diameter (um)
     REAL          :: shape   ! particle shape factor (0-1)
     REAL          :: rhalf   ! radioactive half life (days)
     REAL          :: gpmol   ! gram molecular weight
     REAL          :: acvty   ! pollutant activity ratio for dry deposition
     REAL          :: difty   ! ratio pollutant diffusivity to water vapor
     REAL          :: srate   ! deposition resuspension constant

  END TYPE
