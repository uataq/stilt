!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  CYL2XY           CYLindrical equidistant to XY 
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:01-11-07
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   CONVERTS LAT LON POSITION TO X Y GRID COORDINATES BASED UPON
!   UPON THE CYCLINDRICAL EQUIDISTANT PROJECTION
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 12 Jan 2007 (RRD) - initial version
!
! USAGE:  CALL CYL2XY(PLAT,PLON,X,Y)  
! 
!   INPUT ARGUMENT LIST:     see below
!   OUTPUT ARGUMENT LIST:    see below
!   INPUT FILES:             none
!   OUTPUT FILES:            none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE CYL2XY(PLAT,PLON,X,Y)

  IMPLICIT NONE

  REAL,    INTENT(IN)  :: plat,plon   ! location to compute x,y
  REAL,    INTENT(OUT) :: x,y         ! location of plat,plon

  REAL                 :: xypdeg,coslat,rlat,rlon,xr,yr,tlon 

  COMMON /cylproj/ xypdeg,coslat,rlat,rlon,xr,yr

  TLON=PLON
! internal system always 0-360
  IF(TLON.LT.0.0)   TLON=360.0+TLON
  IF(TLON.GT.360.0) TLON=TLON-360.0

! compute difference from reference longitude
  TLON=TLON-RLON

! difference cannot be greater than 180
  IF(TLON.LT.-180.0) TLON=360.0+TLON
  IF(TLON.GT. 180.0) TLON=TLON-360.0

  X=XYPDEG*TLON*COSLAT+XR
  Y=XYPDEG*(PLAT-RLAT)+YR

END SUBROUTINE cyl2xy

!------------------------------------------------------------------

SUBROUTINE CYL2LL(X,Y,PLAT,PLON)

  IMPLICIT NONE

  REAL,    INTENT(IN)  :: x,y         ! location of plat,plon
  REAL,    INTENT(OUT) :: plat,plon   ! location to compute x,y

  REAL                 :: xypdeg,coslat,rlat,rlon,xr,yr,tlon 

  COMMON /cylproj/ xypdeg,coslat,rlat,rlon,xr,yr

  PLAT=RLAT+(Y-YR)/XYPDEG
  PLON=RLON+(X-XR)/COSLAT/XYPDEG

! return with 0-360 system
! IF(PLON.LT.0.0)   PLON=360.0+PLON
! IF(PLON.GT.360.0) PLON=PLON-360.0

! return with -180 to +180 system
  IF(PLON.LT.-180.0) PLON=360.0+PLON
  IF(PLON.GT. 180.0) PLON=PLON-360.0

END SUBROUTINE cyl2ll

!------------------------------------------------------------------

SUBROUTINE CYLSET(DISTXY,CLAT,PLAT,PLON,XP,YP)

  IMPLICIT NONE

  REAL,     INTENT(IN) :: distxy      ! distance / grid-point
  REAL,     INTENT(IN) :: clat        ! projection center latitude  
  REAL,     INTENT(IN) :: plat,plon   ! reference lat,lon position  
  REAL,     INTENT(IN) :: xp,yp       ! reference position in x,y 

  REAL                 :: xypdeg      ! grid-point / degree-lat
  REAL                 :: coslat
  REAL                 :: rlat,rlon,xr,yr

  REAL,      PARAMETER :: REARTH = 6371.2    ! radius of earth in km
  REAL,      PARAMETER :: PI     = 3.14159265358979
  REAL,      PARAMETER :: DEGPRD = 180.0/PI  ! deg per radian

  COMMON /cylproj/ xypdeg,coslat,rlat,rlon,xr,yr

! position in x-y at reference lat-lon
  RLAT=PLAT
  RLON=PLON
  XR=XP
  YR=YP

! internal system always 0-360
  IF(RLON.LT.0.0)   RLON=360.0+RLON
  IF(RLON.GT.360.0) RLON=RLON-360.0

! latitude scale factor
  COSLAT = COS(CLAT/DEGPRD)

! gp/deg = (km/deg) / (km/gp) 
  XYPDEG = REARTH/DEGPRD/DISTXY 

END SUBROUTINE cylset

!------------------------------------------------------------------

REAL FUNCTION CYLZLL(PLAT,PLON)

  IMPLICIT NONE

  REAL,     INTENT(IN) :: plat,plon   ! reference lat,lon position  

  REAL                 :: distxy,xypdeg,coslat,rlat,rlon,xr,yr

  REAL,      PARAMETER :: REARTH = 6371.2    ! radius of earth in km
  REAL,      PARAMETER :: PI     = 3.14159265358979
  REAL,      PARAMETER :: DEGPRD = 180.0/PI  ! deg per radian

  COMMON /cylproj/ xypdeg,coslat,rlat,rlon,xr,yr

  DISTXY = REARTH/DEGPRD/XYPDEG 
  CYLZLL = SQRT(DISTXY*DISTXY*COSLAT)

END FUNCTION cylzll
