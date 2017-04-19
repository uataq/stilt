!-----------------------------------------------------------
! Module for internal plotting grid for map optimization
!-----------------------------------------------------------
! Last Revised: 24 Jul 2002 (RRD) - initial version DEFPLOT
!               03 Dec 2003 (RRD) - converted to module 
!-----------------------------------------------------------

MODULE mapbox
   
  INTEGER, ALLOCATABLE :: LATLON(:,:)    ! counter array (klat,klon)     
                                         ! 0 to 180 and 0 to 360

  REAL                 :: GINC = 1.0     ! grid spacing degrees
  INTEGER              :: KLAT = 181     ! number of latitudes
  INTEGER              :: KLON = 360     ! number of longitudes
  REAL                 :: GCLAT = -90.0  ! grid corner latitude
  REAL                 :: GCLON = 0.0    ! grid corner longitude

!dwen(20090825)  SAVE latlon,ginc,klat,klon,gclat,gclon

END MODULE mapbox
