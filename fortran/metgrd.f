!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  METGRD           computes METeorological GRiD size
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   METEOROLOGICAL GRID COMPUTES THE GRID SIZE FOR THE CURRENT GRID
!   IN ADDITION LOADS OTHER GRID SENSITIVE PARAMETERS SUCH AS LANDUSE
!   AND ROUGHNESS LENGTH.  VARIABLES PASSED THROUGH AS REQUIRED.
!   THE ROUTINE MUST BE CALLED EACH TIME THE SUBGRID CHANGES
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 11 May 1998 (RRD)
!                  14 Apr 1999 (RRD) - refined grid size = 0 test
!                  29 Sep 2000 (RRD) - fortran90 upgrade
!                  16 Mar 2001 (RRD) - optional global lat lon grid
!                  24 Sep 2001 (RRD) - simultaneous multiple meteorology
!                  05 Dec 2001 (RRD) - option to read terrain file
!                  08 Feb 2002 (RRD) - grid spacing near pole
!                  09 Sep 2002 (RRD) - fortran coding standards
!                  02 Apr 2004 (RRD) - generic file unit numbers
!
! USAGE:  CALL METGRD(KG,KT,LX1,LY1,NXS,NYS,GX,GY,Z0,LU,ZT)
!
!   INPUT ARGUMENT LIST:      see below
!   OUTPUT ARGUMENT LIST:     see below
!   INPUT FILES:              none
!   OUTPUT FILES:             unit KF21 diagnostic MESSAGE file
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$
   
SUBROUTINE METGRD(KG,KT,LX1,LY1,NXS,NYS,GX,GY,Z0,LU,ZT)

  USE funits
  use module_defgrid ! meteorology grid and file

  IMPLICIT NONE


!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  INTEGER,    INTENT(IN)     :: kg            ! grid selection index
  INTEGER,    INTENT(IN)     :: kt            ! time selection index
  INTEGER,    INTENT(IN)     :: lx1,ly1       ! subgrid lower left position
  INTEGER,    INTENT(IN)     :: nxs,nys       ! subgrid dimensions

  REAL,       INTENT(OUT)    :: gx (:,:)      ! grid size array (m)
  REAL,       INTENT(OUT)    :: gy (:,:)      ! grid size array (m)
  REAL,       INTENT(OUT)    :: z0 (:,:)      ! aerodynamic roughness length (m)
  INTEGER,    INTENT(OUT)    :: lu (:,:)      ! land-use category (1-11)
  REAL,       INTENT(INOUT)  :: zt (:,:)      ! terrain height (m)

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  LOGICAL                    :: hset
  INTEGER                    :: ii,jj
  REAL                       :: xi,yj,clat,clon 
  REAL                       :: gdx,gdy,cgszxy_wps

!-------------------------------------------------------------------------------
! external variables
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
  INTERFACE
  SUBROUTINE SFCINP(CLAT,CLON,ZNOT,LUSE,HSET,HGTS)
  IMPLICIT NONE
  REAL,     INTENT(IN)    :: clat,clon        ! Lat/Lon of required point
  REAL,     INTENT(OUT)   :: znot             ! Aerodynamic rougness length (m)
  INTEGER,  INTENT(OUT)   :: luse             ! Land-use categories (1-11)
  LOGICAL,  INTENT(IN)    :: hset             ! Read terrain file flag
  REAL,     INTENT(OUT)   :: hgts             ! terrain height value (m)
  END SUBROUTINE sfcinp
!-------------------------------------------------------------------------------
  SUBROUTINE GBL2LL(KG,KT,X,Y,CLAT,CLON)
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: kg          ! active grid number    
  INTEGER, INTENT(IN)  :: kt          ! active time number    
  REAL,    INTENT(IN)  :: x,y         ! grid position         
  REAL,    INTENT(OUT) :: clat,clon   ! latlon location       
  END SUBROUTINE GBL2LL
!-------------------------------------------------------------------------------
  SUBROUTINE GBLDXY(KG,KT,X,Y,GSX,GSY)
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: kg          ! active grid number    
  INTEGER, INTENT(IN)  :: kt          ! active time number    
  REAL,    INTENT(IN)  :: x,y         ! grid position         
  REAL,    INTENT(OUT) :: gsx,gsy     ! grid size at that location
  END SUBROUTINE GBLDXY
  END INTERFACE
!-------------------------------------------------------------------------------

  HSET=.FALSE.
! only read terrain file for pressure-sigma file with no terrain height
  IF(.NOT.DREC(KG,KT)%SHGT.AND.DREC(KG,KT)%Z_FLAG.EQ.1)HSET=.TRUE.

  DO JJ=1,NYS
  DO II=1,NXS

!    convert index from subgrid to full grid
     XI=FLOAT(II+LX1-1)
     YJ=FLOAT(JJ+LY1-1)

     IF(GRID(KG,KT)%LATLON)THEN
!       find position of grid node
        CALL GBL2LL(KG,KT,XI,YJ,CLAT,CLON)

!       meters per grid cell
        CALL GBLDXY(KG,KT,XI,YJ,GDX,GDY)

!       at pole set spacing to 0.25 grid point below pole (2/8/02)
        IF(CLAT.EQ. 90.0) CALL GBLDXY(KG,KT,XI,(YJ-0.25),GDX,GDY)
        IF(CLAT.EQ.-90.0) CALL GBLDXY(KG,KT,XI,(YJ+0.25),GDX,GDY)

!       convert to meters
        GX(II,JJ)=GDX*1000.0
        GY(II,JJ)=GDY*1000.0

     ELSE
!       meters per grid cell
        GX(II,JJ)=CGSZXY_wps(GRID(KG,KT)%GBASE,XI,YJ,GRID(KG,KT))*1000.0

!       accounts for error in grid conversion at pole points
        IF(GX(II,JJ).LT.1.0)GX(II,JJ)=CGSZXY_wps(GRID(KG,KT)%GBASE,XI+0.5,YJ+0.5,GRID(KG,KT))*1000.0
!       conformal projection
        GY(II,JJ)=GX(II,JJ)

!       find position of grid node
        CALL CXY2LL_wps(GRID(KG,KT)%GBASE,XI,YJ,CLAT,CLON,GRID(KG,KT)%proj)
     END IF

!    uses units 60 (land-use) and 62 (roughness length)
!    to load data into each node from lat/lon based input file
     CALL SFCINP(CLAT,CLON,Z0(II,JJ),LU(II,JJ),HSET,ZT(II,JJ))

  END DO
  END DO

  WRITE(KF21,*)' NOTICE metgrd: (kg, xyr,xy1) - ',KG,NXS,NYS,LX1,LY1

END SUBROUTINE metgrd  
