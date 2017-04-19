!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  ADVRNG           ADVection RaNGe to determine subgrid
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:99-03-03
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   ADVECTION RANGE COMPUTES THE RANGE OF PARTICLE POSITIONS ON THE
!   COMPUTATIONAL GRID (METEO GRID) WHICH IS USED TO DETERMINE THE
!   OPTIMUM LOCATION AND SIZE OF THE METEOROLOGICAL SUBGRID. THE
!   SUBGRID SIZE IS ALSO DETERMINED BY THE FREQUENCY OF THE METEO DATA
!   IN THAT IT IS NECESSARY TO CONSIDER HOW LONG A PARTICLE REMAINS ON
!   THE SUBGRID BEFORE MORE DATA IS REQUIRED TO BE LOADED.  IT IS BETTER
!   TO AVOID LOADING DATA MORE THAN ONCE FOR ANY TIME PERIOD.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 16 Mar 1999 (RRD) - Initial version of the code
!                 01 Dec 1999 (RRD) - added minimum grid size from namelist
!                 03 Sep 2000 (RRD) - fortran90 upgrade
!                 05 Mar 2001 (RRD) - full grid option adjustment
!                 15 Mar 2001 (RRD) - global cyclic boundary conditions
!                 02 Oct 2001 (RRD) - simultaneous multiple meteorology
!                 17 Jan 2002 (RRD) - umax units & index change in kg loop
!                 09 Feb 2002 (RRD) - global when either at 75%
!                 13 Aug 2002 (RRD) - max subgrid limited to min fullgrid
!                 09 Sep 2002 (RRD) - fortran coding standards
!                 29 May 2003 (RRD) - test for various grid combinations
!                 02 Apr 2004 (RRD) - generic file unit numbers
!                 30 Apr 2006 (RRD) - removed global test (conflict with metsub)
!
! USAGE:  CALL ADVRNG(NGRD,MGMIN,UMAX,KPM,XPOS,YPOS,PGRD)
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

SUBROUTINE ADVRNG(NGRD,MGMIN,UMAX,KPM,XPOS,YPOS,PGRD)

  USE funits
  use module_defgrid ! meteorology grid and file

  IMPLICIT NONE


!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  INTEGER,    INTENT(IN)    :: ngrd     ! minimum subgrid size (from namelist)
!dwen(20090825)  INTEGER,    INTENT(IN)    :: mgmin    ! minimum subgrid size (from namelist)
  INTEGER,    INTENT(IN)    :: mgmin    ! minimum subgrid size (from namelist)
  REAL,       INTENT(IN)    :: umax     ! maximum wind speed (km / min)
  INTEGER,    INTENT(IN)    :: kpm      ! number of particles
  REAL,       INTENT(IN)    :: xpos (:) ! particle center positions (grid units)
  REAL,       INTENT(IN)    :: ypos (:) ! particle center positions (grid units)
  INTEGER,    INTENT(IN)    :: pgrd (:) ! particle meteorological grid 

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  INTEGER, ALLOCATABLE  :: lymax(:),lymin(:),lxmax(:),lxmin(:),mgmin_save(:)

  INTEGER               :: kt = 1          ! analysis the same for all times
                                           ! such that the spatial extent of 
                                           ! all grids under kt index identical

  INTEGER               :: k2,kg,kp,kret
  REAL                  :: xpt,ypt,dist,tlat,tlon

!-------------------------------------------------------------------------------
  SAVE lymax,lymin,lxmax,lxmin,mgmin_save
!-------------------------------------------------------------------------------

! check for sufficient number
  IF(KPM.LE.0)RETURN

  IF(.NOT.ALLOCATED(lymax))THEN
     ALLOCATE(lymax(ngrd), STAT=kret)
     ALLOCATE(lymin(ngrd), STAT=kret)
     ALLOCATE(lxmax(ngrd), STAT=kret)
     ALLOCATE(lxmin(ngrd), STAT=kret)
     ALLOCATE(mgmin_save(ngrd), STAT=kret)

!    initialize grid domain range
     GRID%LXR=0
     GRID%LYR=0
!dwen(20090815) *************************
! CHG&JCL (03/10/2004) change minimum subgrid dimension from 10
!     a value depending on grid size
!     at 40 km and larger still keep 10, at smaller grid sizes keep
!     physical dimension constant
!dwen(20090825)      MGMIN=MAX(MGMIN,IDNINT(MGMIN*40.0/GRID(KG)%SIZE))
     do kg=1,ngrd
        MGMIN_save(kg)=MAX(MGMIN,NINT(MGMIN*40.0/GRID(KG,kt)%SIZE))
     end do
     
  END IF

! set maximimum and minimum limits for each grid
  DO KG=1,NGRD
     LXMAX(KG)=1
     LXMIN(KG)=GRID(KG,KT)%NX
     LYMAX(KG)=1
     LYMIN(KG)=GRID(KG,KT)%NY
  END DO

! For all positions deterime max/min limits, but for multiple grids       
! convert limits from one grid to all other grids to determine if particles
! should influence subgrid limits on other grids besides their current 
! computational grid.

  DO KP=1,KPM
     KG=PGRD(KP)

     DO K2=1,NGRD
        IF(KG.NE.K2)THEN
!          for multiple grids convert position from current grid to true
!          then convert true back to new grid coordinates (mod 05/29/2003)
           IF(GRID(KG,1)%LATLON.AND.GRID(K2,1)%LATLON)THEN
              CALL GBL2LL(KG,KT,XPOS(KP),YPOS(KP),TLAT,TLON)
              CALL GBL2XY(K2,KT,TLAT,TLON,XPT,YPT)
           ELSEIF(GRID(KG,1)%LATLON.AND.(.NOT.GRID(K2,1)%LATLON))THEN
              CALL GBL2LL(KG,KT,XPOS(KP),YPOS(KP),TLAT,TLON)
              CALL CLL2XY_wps(GRID(K2,KT)%GBASE,TLAT,TLON,XPT,YPT,GRID(K2,KT)%proj)
           ELSEIF(.NOT.GRID(KG,1)%LATLON.AND.GRID(K2,1)%LATLON )THEN
              CALL CXY2LL_wps(GRID(KG,KT)%GBASE,XPOS(KP),YPOS(KP),TLAT,TLON,GRID(KG,KT)%proj)
              CALL GBL2XY(K2,KT,TLAT,TLON,XPT,YPT)
           ELSE
              CALL CXY2LL_wps(GRID(KG,KT)%GBASE,XPOS(KP),YPOS(KP),TLAT,TLON,GRID(KG,KT)%proj)
              CALL CLL2XY_wps(GRID(K2,KT)%GBASE,TLAT,TLON,XPT,YPT,GRID(K2,KT)%proj)
           END IF

        ELSE
           XPT=XPOS(KP)
           YPT=YPOS(KP)
        END IF

!       particle must be within the main grid to be valid (1/16/2002)
        IF(INT(XPT).GE.1.AND.INT(XPT).LT.GRID(K2,KT)%NX.AND.   &
           INT(YPT).GE.1.AND.INT(YPT).LT.GRID(K2,KT)%NY)THEN     
           LXMAX(K2)=MAX(LXMAX(K2),INT(XPT)+1)
           LXMIN(K2)=MIN(LXMIN(K2),INT(XPT))
           LYMAX(K2)=MAX(LYMAX(K2),INT(YPT)+1)
           LYMIN(K2)=MIN(LYMIN(K2),INT(YPT))
        END IF
     END DO
  END DO

!-------------------------------------------------------------------------------
! determine subgrid corner, range, and center position for each grid using kt=1
! as the test grid and writing the same result for that grid to all kt indicies
!-------------------------------------------------------------------------------

  DO KG=1,NGRD

!    Max transport distance = (grid_units/min)*(meteo time interval)
!    subgrid should be "dist" larger to contain particles on edge

     DIST=MAX(1.0,UMAX*FLOAT(DREC(KG,KT)%DELTA)/GRID(KG,KT)%SIZE)

!    The subgrid range should be at least MGMIN grid units or the size
!    determined from the advection distance plus particle distribution. 
!    Subgrid is not permitted to shrink during a simulation


!****************************************

     GRID(KG,:)%LXR=MAX(MGMIN_save(kg),(NINT(2.0*DIST)+LXMAX(KG)-LXMIN(KG)+3),&
                     GRID(KG,KT)%LXR)
     GRID(KG,:)%LYR=MAX(MGMIN_save(kg),(NINT(2.0*DIST)+LYMAX(KG)-LYMIN(KG)+3),&
                     GRID(KG,KT)%LYR)
 
!    when subgrid in either direction reaches 75% of the full grid dimension
!    then set limits to maximum or previously set subgrid to global
     IF((GRID(KG,KT)%LXR.GT.NINT(0.75*FLOAT(GRID(KG,KT)%NX))).OR. & 
        (GRID(KG,KT)%LYR.GT.NINT(0.75*FLOAT(GRID(KG,KT)%NY))))THEN
        
!        .OR.GRID(KG,KT)%GLOBAL)THEN -- removed from above if test (30 Apr 2006)

        GRID(KG,:)%LXR=GRID(KG,KT)%NX
        GRID(KG,:)%LXC=FLOAT(GRID(KG,KT)%NX+1)/2.0
        GRID(KG,:)%LX1=1

        GRID(KG,:)%LYR=GRID(KG,KT)%NY
        GRID(KG,:)%LYC=(GRID(KG,KT)%NY+1.0)/2.0
        GRID(KG,:)%LY1=1

!       input grid is global then set subgrid to global
!       comment following line due to conflict with metsub global test (30 Apr 2006)
!       IF(GRID(KG,KT)%GBLDAT) GRID(KG,KT)%GLOBAL=.TRUE.

     ELSE
!       compute grid center and corner position
        GRID(KG,:)%LXC=FLOAT(LXMAX(KG)+LXMIN(KG))/2.0
        GRID(KG,:)%LX1=MAX(1,NINT(GRID(KG,KT)%LXC-FLOAT(GRID(KG,KT)%LXR)/2.0))
!       check upper end to avoid subgrid compression
        IF(GRID(KG,KT)%LX1+GRID(KG,KT)%LXR.GT.GRID(KG,KT)%NX) &
           GRID(KG,: )%LX1=GRID(KG,KT)%NX-GRID(KG,KT)%LXR

!       compute grid center and corner position
        GRID(KG,:)%LYC=(LYMAX(KG)+LYMIN(KG))/2.0
        GRID(KG,:)%LY1=MAX(1,NINT(GRID(KG,KT)%LYC-GRID(KG,KT)%LYR/2.0))
!       check upper end to avoid subgrid compression
        IF(GRID(KG,KT)%LY1+GRID(KG,KT)%LYR.GT.GRID(KG,KT)%NY) &
           GRID(KG,: )%LY1=GRID(KG,KT)%NY-GRID(KG,KT)%LYR

     END IF

!    Subgrid data load flag turned off. When particle moves off the subgrid as
!    determined in metsub, the flag is set to true. Subsequent offgrid particles
!    result in an automatic expansion of the subgrid until the flag is reset
!    in this routine resulting in the definition of a new optimal subgrid.

     GRID(KG,:)%DATLOAD=.FALSE.

!    optional diagnostic message
     WRITE(KF21,*)' NOTICE advrng: (kg , mgmin, xyr,xy1) - ',KG, mgmin_save(kg), &
          GRID(KG,KT)%LXR,GRID(KG,KT)%LYR,GRID(KG,KT)%LX1,GRID(KG,KT)%LY1 

  END DO

END SUBROUTINE advrng
