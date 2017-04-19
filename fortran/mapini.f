!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  MAPINI           OPTIMIZE MAP DIMENSIONS
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   TO GENERATE THE MAP REQUIREMENTS FOR PLOTTING
!   OPTIMIZES MAP DIMENSIONS FOR DISPLAY DEPENDING UPON THE PROJECTION
!   AND THE DIMENSIONS OF THE CONCENTRATION PLUME OR TRAJECTORY
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 07 Dec 1998 (RRD)
!                 28 Mar 2000 (RRD) - revised based upon centroid
!                 14 Apr 2000 (RRD) - map limits reflon patch
!                 03 Nov 2000 (RRD) - internal grid spacing
!                 20 Nov 2000 (RRD) - fortran90 upgrade
!                 29 Oct 2001 (RRD) - variable map margin from argument
!                 14 Jan 2002 (RRD) - single point map specification
!                 25 Feb 2002 (RRD) - revised scaling optimization scheme
!                 17 Jul 2002 (RRD) - major revision      
!                 09 Sep 2002 (RRD) - small map grid size minimum
!                 04 Dec 2003 (RRD) - converted include to module
!                 12 May 2004 (RRD) - restrict lambert projections
!                 12 Jan 2007 (RRD) - add cylindrical equidistant
!
! USAGE: CALL MAPINI(PARMAP,SCALE,NXP,NYP,DLAT,DLON,FRAC,IPROJ,
!                    OLAT,OLON,ALONL,ALONR,ALATT,ALATB)
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

SUBROUTINE MAPINI(PARMAP,SCALE,NXP,NYP,DLAT,DLON,FRAC,IPROJ,   &
                  OLAT,OLON,ALONL,ALONR,ALATT,ALATB)

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  REAL,    INTENT(OUT)   :: parmap (9)    ! conformal map projection parameters
  REAL,    INTENT(IN)    :: dlat,dlon     ! resolution of input data grid
  REAL,    INTENT(IN)    :: frac          ! margin of map outside plume
  INTEGER, INTENT(INOUT) :: iproj         ! map projection type  
  REAL,    INTENT(IN)    :: olat,olon     ! computation origin
  REAL,    INTENT(OUT)   :: alonl,alatb   ! limits of final optimized map
  REAL,    INTENT(OUT)   :: alonr,alatt   ! limits of final optimized map
  INTEGER, INTENT(OUT)   :: nxp,nyp       ! grid size of grid for input data
  REAL,    INTENT(IN)    :: scale         ! map scale factor (delta-x / delta-y)

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  LOGICAL, PARAMETER :: DIAG = .FALSE.
! LOGICAL, PARAMETER :: DIAG = .TRUE.  

  INTEGER            :: kproj
  REAL               :: SLAT,SLON,GLAT,GLON
  REAL               :: qlat,qlon,clat,clon,delx,dely
  REAL               :: xc,yc,x1,y1,x2,y2
  REAL               :: grid,reflon,tnglat,cgszll 

!-------------------------------------------------------------------------------

  INTERFACE
  SUBROUTINE MAPOPT(DIAG,IPROJ,KPROJ,PARMAP,QLAT,QLON,X1,Y1,X2,Y2,REFLON,  &
                    TNGLAT,SLAT,SLON,GLAT,GLON,GRID)
  IMPLICIT NONE
  LOGICAL, INTENT(IN)    :: DIAG
  INTEGER, INTENT(IN)    :: IPROJ
  INTEGER, INTENT(OUT)   :: KPROJ
  REAL,    INTENT(OUT)   :: PARMAP(9)
  REAL,    INTENT(INOUT) :: QLAT,QLON
  REAL,    INTENT(OUT)   :: X1,Y1,X2,Y2
  REAL,    INTENT(IN)    :: GRID
  REAL,    INTENT(OUT)   :: REFLON,TNGLAT,SLAT,SLON,GLAT,GLON
  END SUBROUTINE mapopt
  END INTERFACE

!-------------------------------------------------------------------------------

! grid spacing as 0.5 of input grid spacing (rrd: 10/3/00)
  GRID=0.5*MIN(DLAT,DLON)*100.0
  IF(DIAG)WRITE(*,*)'Spacing: ',GRID 

! arbitrary projection to determine location center
  QLAT=OLAT
  QLON=OLON

! find best projection
  CALL MAPOPT(DIAG,IPROJ,KPROJ,PARMAP,QLAT,QLON,X1,Y1,X2,Y2,REFLON,TNGLAT, &
              SLAT,SLON,GLAT,GLON,GRID)

! fine tune projection
! CALL MAPOPT(DIAG,IPROJ,KPROJ,PARMAP,QLAT,QLON,X1,Y1,X2,Y2,REFLON,TNGLAT, &
!             SLAT,SLON,GLAT,GLON,GRID)

! lambert grids not permitted to encompass the poles
  IF(KPROJ.EQ.2)THEN
     CALL CLL2XY(PARMAP,SIGN(90.0,QLAT),0.0,XC,YC)
     IF(DIAG)WRITE(*,*)'Pole xy: ',XC,YC
!    force polar sterographic projection
     IF(XC.GE.X1.AND.XC.LE.X2.AND.YC.GE.Y1.AND.YC.LE.Y2) THEN
        IF(DIAG)WRITE(*,*)'Force polar stereographic!'
        IPROJ=1 
        CALL MAPOPT(DIAG,IPROJ,KPROJ,PARMAP,QLAT,QLON,X1,Y1,X2,Y2,REFLON, &
                    TNGLAT,SLAT,SLON,GLAT,GLON,GRID)
     END IF
  END IF

! new map center
  XC=0.5*(X1+X2)
  YC=0.5*(Y1+Y2)

! scale map according to aspect ratio
  IF(ABS((X2-X1)/(Y2-Y1)).LE.SCALE)THEN
!    expand in x-direction
     DELX=(Y2-Y1)*SCALE
     X1=XC-DELX/2.0
     X2=XC+DELX/2.0
  ELSE
!    expand in y-direction
     DELY=(X2-X1)/SCALE
     Y1=YC-DELY/2.0
     Y2=YC+DELY/2.0
  END IF
  IF(DIAG)WRITE(*,*)'X,Y Asp: ',X1,Y1,X2,Y2

! projection zoom factor
  DELX=ABS(X2-X1)
  DELY=ABS(Y2-Y1)
  Y1=Y1-SIGN(FRAC*DELY,(Y2-Y1))
  X1=X1-SIGN(FRAC*DELX,(X2-X1))
  Y2=Y2+SIGN(FRAC*DELY,(Y2-Y1))
  X2=X2+SIGN(FRAC*DELX,(X2-X1))
  IF(DIAG)WRITE(*,*)'X,Y Zum: ',X1,Y1,X2,Y2

! round map corners to match even grid index for plotting
  Y1=NINT(Y1)
  Y2=NINT(Y2)
  DELX=(Y2-Y1)*SCALE
  X1=NINT(X1)
  X2=X1+NINT(DELX)
  IF(DIAG)WRITE(*,*)'X,Y Adj: ',X1,Y1,X2,Y2

! compute new map corners
  IF(KPROJ.EQ.4)THEN
     CALL CYL2LL(X1,Y1,ALATB,ALONL)
     CALL CYL2LL(X2,Y2,ALATT,ALONR)
  ELSE
     CALL CXY2LL(PARMAP,X1,Y1,ALATB,ALONL)
     CALL CXY2LL(PARMAP,X2,Y2,ALATT,ALONR)
  END IF
  IF(DIAG)WRITE(*,*)'Corners: ',ALATB,ALONL,ALATT,ALONR

! map exceeds limits
  IF(ALATT.GT.90.0.OR.ALATB.LT.-90.0)THEN 
     WRITE(*,*)'Warning: map projection exceeds limits!'
     WRITE(*,*)'Increase zoom or change/force projection'
  END IF

! rescale map by defining 1,1 at lower left corner
  X1=1.0
  Y1=1.0
  IF(KPROJ.EQ.4)THEN
     CALL CYLSET(GRID,TNGLAT,ALATB,ALONL,X1,Y1)
     CALL CYL2XY(ALATT,ALONR,X2,Y2)
  ELSE
     CALL STLMBR(PARMAP,TNGLAT,REFLON)
     CALL STCM1P(PARMAP,X1,Y1,ALATB,ALONL,GLAT,GLON,GRID,0.0)
     CALL CLL2XY(PARMAP,ALATT,ALONR,X2,Y2)
  END IF
  IF(DIAG)WRITE(*,*)'Final  : ',X1,Y1,X2,Y2

! number of points
  NXP=NINT(X2)  
  NYP=NINT(Y2) 

! return updated projection
  IPROJ=KPROJ

END SUBROUTINE mapini

!-------------------------------------------------------------------------------
! map optimization   
!-------------------------------------------------------------------------------

SUBROUTINE MAPOPT(DIAG,IPROJ,KPROJ,PARMAP,QLAT,QLON,X1,Y1,X2,Y2,REFLON,   &
                  TNGLAT,SLAT,SLON,GLAT,GLON,GRID)

  USE mapbox

  IMPLICIT NONE

  LOGICAL, INTENT(IN)    :: DIAG
  INTEGER, INTENT(IN)    :: IPROJ
  INTEGER, INTENT(OUT)   :: KPROJ
  REAL,    INTENT(OUT)   :: PARMAP(9)
  REAL,    INTENT(INOUT) :: QLAT,QLON
  REAL,    INTENT(OUT)   :: X1,Y1,X2,Y2
  REAL,    INTENT(IN)    :: GRID
  REAL,    INTENT(OUT)   :: REFLON,TNGLAT,SLAT,SLON,GLAT,GLON

  INTEGER :: ig,jg
  REAL    :: xc,yc,plat,plon,alatb,alatt,alonl,alonr

  IF(IPROJ.EQ.0)THEN
!    determine projection type (1-polar 2-lambert 3-mercator 4-cycl equid)
     KPROJ=2
     IF(QLAT.GT.55.0.OR. QLAT.LT.-55.0)KPROJ=1
     IF(QLAT.LT.25.0.AND.QLAT.GT.-25.0)KPROJ=3
  ELSE
     KPROJ=IPROJ
  END IF

! projection cut
  REFLON=QLON

! set tangent latitude
  IF(KPROJ.EQ.1)THEN
     IF(QLAT.GE.0.0)THEN
        TNGLAT= 90.0
     ELSE
        TNGLAT=-90.0
     END IF
     SLAT=TNGLAT
     SLON=0.0
     GLAT=QLAT
     GLON=REFLON
  ELSEIF(KPROJ.EQ.2)THEN
     TNGLAT=QLAT
     SLAT=QLAT
     SLON=QLON
     GLAT=QLAT
     GLON=QLON
  ELSEIF(KPROJ.EQ.3.OR.KPROJ.EQ.4)THEN
     TNGLAT=0.0
     SLAT=QLAT
     SLON=QLON
     GLAT=0.0  
     GLON=REFLON
  END IF

  IF(DIAG)WRITE(*,*)'Map Prj: ',KPROJ
  IF(DIAG)WRITE(*,*)'Center : ',QLAT,QLON
  IF(DIAG)WRITE(*,*)'Tangent: ',TNGLAT,'  Reflon: ',REFLON

! initial point
  XC=500.0
  YC=500.0 

  IF(KPROJ.EQ.4)THEN 
     CALL CYLSET(GRID,TNGLAT,SLAT,SLON,XC,YC)
     CALL CYL2XY(QLAT-GINC/2.0,QLON-GINC/2.0,X1,Y1)
     CALL CYL2XY(QLAT+GINC/2.0,QLON+GINC/2.0,X2,Y2)
  ELSE
!    create projection
     CALL STLMBR(PARMAP,TNGLAT,REFLON)
     CALL STCM1P(PARMAP,XC,YC,SLAT,SLON,GLAT,GLON,GRID,0.0)
!    find new map center and map corners
     CALL CLL2XY(PARMAP,QLAT-GINC/2.0,QLON-GINC/2.0,X1,Y1)
     CALL CLL2XY(PARMAP,QLAT+GINC/2.0,QLON+GINC/2.0,X2,Y2)
  END IF
  IF(DIAG)WRITE(*,*)'X,Y Set : ',X1,Y1,X2,Y2
  IF(DIAG)WRITE(*,*)'Grid corner, increment from mapbox or conndx :',GCLAT,GCLON,GINC
  IF(DIAG)WRITE(*,*)'1st guess center : ',QLAT,QLON

  DO JG=1,KLON
  DO IG=1,KLAT
     IF(LATLON(IG,JG).GT.0)THEN
        PLAT=(IG-1)*GINC+GCLAT
        PLON=(JG-1)*GINC+GCLON
        IF(PLON.GT.180.0)PLON=PLON-360.0
        IF(KPROJ.EQ.4)THEN
           CALL CYL2XY(PLAT,PLON,XC,YC)
        ELSE
           CALL CLL2XY(PARMAP,PLAT,PLON,XC,YC)
        END IF
        X1=MIN(X1,XC)
        Y1=MIN(Y1,YC)
        X2=MAX(X2,XC)
        Y2=MAX(Y2,YC)
     END IF
  END DO
  END DO
  IF(DIAG)WRITE(*,*)'X,Y Ini: ',X1,Y1,X2,Y2

  XC=0.5*(X1+X2)
  YC=0.5*(Y1+Y2)
  IF(KPROJ.EQ.4)THEN
     CALL CYL2LL(XC,YC,QLAT,QLON)
  ELSE
     CALL CXY2LL(PARMAP,XC,YC,QLAT,QLON)
  END IF
  IF(DIAG)WRITE(*,*)'Center : ',QLAT,QLON

! compute new map corners
  IF(DIAG)THEN
     IF(KPROJ.EQ.4)THEN
        CALL CYL2LL(X1,Y1,ALATB,ALONL)
        CALL CYL2LL(X2,Y2,ALATT,ALONR)
     ELSE
        CALL CXY2LL(PARMAP,X1,Y1,ALATB,ALONL)
        CALL CXY2LL(PARMAP,X2,Y2,ALATT,ALONR)
     END IF
     WRITE(*,*)'Corners: ',ALATB,ALONL,ALATT,ALONR
  END IF

END SUBROUTINE mapopt
