!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  GBLDLL           GloBaL grid size by LAT LON coord
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:01-09-03
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   DETERMINES THE GRID SPACING IN KM FOR A LAT LON POSITION FOR
!   LAT LON GRID SYSTEMS.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 20 Apr 2001 (RRD) - initial version
!                 29 Aug 2001 (RRD) - simultaneous multiple meteo
!                 09 Sep 2002 (RRD) - fortran coding standards
!
! USAGE:  CALL GBLDLL(KG,KT,CLAT,GSX,GSY)
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

SUBROUTINE GBLDLL(KG,KT,CLAT,GSX,GSY)

  use module_defgrid ! meteorology grid and file
  IMPLICIT NONE


  INTEGER, INTENT(IN)  :: kg          ! active grid number    
  INTEGER, INTENT(IN)  :: kt          ! active time number    
  REAL,    INTENT(IN)  :: clat        ! grid position         
  REAL,    INTENT(OUT) :: gsx,gsy     ! grid size at that location

  REAL, PARAMETER :: REARTH = 6371.2    ! radius of earth in km
  REAL, PARAMETER :: PI     = 3.14159265358979
  REAL, PARAMETER :: DEGPRD = 180.0/PI  ! deg per radian


!-------------------------------------------------------------------------------

  IF(.NOT.GRID(KG,KT)%LATLON) RETURN  

! latitude grid spacing 
  GSY = REARTH*GRID(KG,KT)%REF_LAT/DEGPRD

! longitude grid spacing 
  GSX = COS(CLAT/DEGPRD)*REARTH*GRID(KG,KT)%REF_LON/DEGPRD 

END SUBROUTINE gbldll
