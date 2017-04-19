!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  GBLDXY           GloBaL grid size by x,y coordinate
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:01-09-03
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   DETERMINES THE GRID SPACING IN KM FOR A X,Y COORDINATE POSITION.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 12 Mar 2001 (RRD) - initial version
!                 29 Aug 2001 (RRD) - simultaneous multiple meteo
!                 09 Sep 2002 (RRD) - fortran coding standards
!
! USAGE:  CALL GBLDXY(KG,KT,X,Y,GSX,GSY)
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

SUBROUTINE GBLDXY(KG,KT,X,Y,GSX,GSY)

  use module_defgrid ! meteorology grid and file
  IMPLICIT NONE


  INTEGER, INTENT(IN)  :: kg          ! active grid number    
  INTEGER, INTENT(IN)  :: kt          ! active time number    
  REAL,    INTENT(IN)  :: x,y         ! grid position         
  REAL,    INTENT(OUT) :: gsx,gsy     ! grid size at that location

  REAL            :: CLAT,CLON
  REAL, PARAMETER :: REARTH = 6371.2    ! radius of earth in km
  REAL, PARAMETER :: PI     = 3.14159265358979
  REAL, PARAMETER :: DEGPRD = 180.0/PI  ! deg per radian


!-------------------------------------------------------------------------------
  INTERFACE
    SUBROUTINE GBL2LL(KG,KT,X,Y,CLAT,CLON)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: kg          ! active grid number    
    INTEGER, INTENT(IN)  :: kt          ! active time number    
    REAL,    INTENT(IN)  :: x,y         ! grid position         
    REAL,    INTENT(OUT) :: clat,clon   ! latlon location       
    END SUBROUTINE gbl2ll
  END INTERFACE
!-------------------------------------------------------------------------------

  IF(.NOT.GRID(KG,KT)%LATLON) RETURN  

  CALL GBL2LL(KG,KT,X,Y,CLAT,CLON)

! latitude grid spacing 
  GSY = REARTH*GRID(KG,KT)%REF_LAT/DEGPRD

! longitude grid spacing 
  GSX = COS(CLAT/DEGPRD)*REARTH*GRID(KG,KT)%REF_LON/DEGPRD 

END SUBROUTINE gbldxy
