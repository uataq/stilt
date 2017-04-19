!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  GBL2LL           GloBaL position to Lat Lon      
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:01-09-03
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   CONVERTS X,Y GRID POSITION TO LAT LON COORDINATES BASED UPON
!   UPON THE GRID SPACING SPECIFIED.          
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 09 Mar 2001 (RRD) - initial version
!                 29 Aug 2001 (RRD) - simultaneous multiple meteo
!                 09 Sep 2002 (RRD) - fortran coding standards
!
! USAGE:  CALL GBL2LL(KG,KT,X,Y,CLAT,CLON)  
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

SUBROUTINE GBL2LL(KG,KT,X,Y,CLAT,CLON)

  use module_defgrid ! meteorology grid and file
  IMPLICIT NONE


  INTEGER, INTENT(IN)  :: kg          ! active grid number    
  INTEGER, INTENT(IN)  :: kt          ! active time number    
  REAL,    INTENT(IN)  :: x,y         ! grid position         
  REAL,    INTENT(OUT) :: clat,clon   ! latlon location       


  IF(.NOT.GRID(KG,KT)%LATLON) RETURN  

! Grid system is simply defined as the number of grid points
! from the corner point at 1,1 using an even lat-lon increment
! for the x and y directions. Grid distances are computed
! where needed according to the latitude of the grid point

  CLAT=GRID(KG,KT)%SYNC_LAT+(Y-1.0)*GRID(KG,KT)%REF_LAT
  IF(CLAT.GT. 90.0)CLAT= 180.0-CLAT
  IF(CLAT.LT.-90.0)CLAT=-180.0-CLAT

  CLON=GRID(KG,KT)%SYNC_LON+(X-1.0)*GRID(KG,KT)%REF_LON
  CLON=MOD(CLON,360.0)
  IF(CLON.GT.180.0)CLON=CLON-360.0

END SUBROUTINE gbl2ll
