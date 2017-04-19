!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  GBL2XY           GloBaL position to XY cordinates
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:01-09-03
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   CONVERTS LAT LON POSITION TO X Y GRID COORDINATES BASED UPON
!   UPON THE GRID SPACING SPECIFIED.          
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 12 Mar 2001 (RRD) - initial version
!                 29 Aug 2001 (RRD) - simultaneous multiple meteo
!                 18 Dec 2001 (RRD) - removed cyclic BC test
!                 21 Feb 2002 (RRD) - prime meridian test
!                 10 Aug 2006 (RRD) - refined prime meridian test
!
! USAGE:  CALL GBL2XY(KG,KT,CLAT,CLON,X,Y)  
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

SUBROUTINE GBL2XY(KG,KT,CLAT,CLON,X,Y)

  use module_defgrid ! meteorology grid and file
  IMPLICIT NONE


  INTEGER, INTENT(IN)  :: kg          ! active grid number    
  INTEGER, INTENT(IN)  :: kt          ! active time number    
  REAL,    INTENT(IN)  :: clat,clon   ! latlon location       
  REAL,    INTENT(OUT) :: x,y         ! grid position         
  REAL                 :: tlat,tlon 


  IF(.NOT.GRID(KG,KT)%LATLON) RETURN  

! Grid system is simply defined as the number of grid points
! from the corner point at 1,1 using an even lat-lon increment
! for the x and y directions. Grid distances are computed
! where needed according to the latitude of the grid point

  TLAT=CLAT
  IF(TLAT.GT. 90.0)TLAT= 180.0-TLAT
  IF(TLAT.LT.-90.0)TLAT=-180.0-TLAT
  Y=1.0+(TLAT-GRID(KG,KT)%SYNC_LAT)/GRID(KG,KT)%REF_LAT

  TLON=CLON
  IF(.NOT.GRID(KG,KT)%PRIME)THEN
!    use 0-360 system except near prime
     IF(TLON.LT.0.0)  TLON=360.0+TLON
     IF(TLON.GT.360.0)TLON=TLON-360.0
  END IF
  TLON=TLON-GRID(KG,KT)%SYNC_LON
  IF(TLON.LT.0.0)TLON=TLON+360.0
  X=1.0+TLON/GRID(KG,KT)%REF_LON

END SUBROUTINE gbl2xy
