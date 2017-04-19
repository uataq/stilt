!$$$  SUBPROGRAM DOCUMENTATION BLOCK   
!
! SUBPROGRAM:  METSUB           METeorological POSitioning finds record
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:01-09-07
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   METEOROLOGICAL POSITIONING ROUTINE USES THE CURRENT CALCULATION
!   POSITION TO CHECK IF THE POSITION FALLS WITHIN THE METEOROLOGICAL
!   SUBGRID. IF NOT THE NEW SUBGRID LIMITS ARE CALCULATED. NOTE THAT 
!   A SUBGRID CHANGE REQUIRES DATA AT THE OLD AND NEW TIMES BE LOADED 
!   FOR TEMPORAL INTERPOLATION.
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 15 May 1998 (RRD)
!                  05 Oct 2001 (RRD) - initial version from metpos
!                  02 Dec 2001 (RRD) - full grid test
!                  17 Jan 2002 (RRD) - disable kret=1 test  
!                  09 Feb 2002 (RRD) - global subgrid option
!                  28 Mar 2002 (RRD) - global switch for all files
!                  13 Aug 2002 (RRD) - test for grid limits exceeded
!                  09 Sep 2002 (RRD) - fortran coding standards
!                  14 Mar 2003 (RRD) - no subgrid expansion at edge
!                  26 Mar 2003 (RRD) - set range when grid exceeds limit
!                  17 Sep 2003 (RRD) - correction to 26 Mar correction
!                  02 Apr 2004 (RRD) - generic file unit numbers
!                  12 May 2005 (RRD) - subgrid expansion after shift fails
!
! USAGE:  CALL METSUB(XP,YP,KG,KT,OFFG,LX1,LY1,LXR,LYR,KRET)
!
!   INPUT ARGUMENT LIST:     see below
!   OUTPUT ARGUMENT LIST:    see below
!   INPUT FILES:             none
!   OUTPUT FILES:            unit KF21 diagnostic ouput MESSAGE file
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

 SUBROUTINE METSUB(XP,YP,KG,KT,OFFG,LX1,LY1,LXR,LYR,KRET)

  USE funits
  use module_defgrid ! meteorology grid and file

  IMPLICIT NONE


!-------------------------------------------------------------------------------
! argument list definitions
!-------------------------------------------------------------------------------

  REAL,      INTENT(IN)    :: xp,yp          ! x,y grid particle position
  INTEGER,   INTENT(IN)    :: kg             ! grid number being loaded
  INTEGER,   INTENT(IN)    :: kt             ! grid number being loaded
  LOGICAL,   INTENT(OUT)   :: offg           ! load new subgrid flag  
  INTEGER,   INTENT(INOUT) :: lx1,ly1        ! subgrid corner point
  INTEGER,   INTENT(INOUT) :: lxr,lyr        ! subgrid range        
  INTEGER,   INTENT(OUT)   :: kret           ! termination flag

!-------------------------------------------------------------------------------
! internal definitions
!-------------------------------------------------------------------------------

  LOGICAL   :: expand,keep_expanding
  INTEGER   :: ii,jj
  REAL      :: lxc,lyc
  integer :: max_offset

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

! start default assumption that advection point is on the grid
  OFFG=.FALSE.
  KRET=0
  IF(GRID(KG,KT)%GLOBAL) RETURN 

! convert position to index units
  II=INT(XP)
  JJ=INT(YP)

!-------------------------------------------------------------------------------
! Determine if data load is required because position has moved off sub-grid.
! Previous subroutine determined point was somewhere on the full grid. Now the
! position is checked against the current subgrid and if off the subgrid, the
! meteorological data are reloaded into a repositioned or enlarged subgrid.
! Particle must be within the outer ring of grid cells as the transport distance
! over one advection step may otherwise move the particle off the subgrid. 
!-------------------------------------------------------------------------------

  max_offset = 2
!       For WRF input, stay away from one additional row/column of grid edge
!       (because rows/columns for non-staggered field are padded)
  if (GRID(KG,kt)%MODEL_ID(2:4) .eq. 'WRF') max_offset = max_offset+1
!cg(20100804) same for ECMWF fields
  if (GRID(KG,kt)%MODEL_ID(1:2) .eq. 'EC') max_offset = max_offset+1
  
  IF(II .LT. LX1+2 .OR. II .GE. (LX1+LXR-max_offset) .OR. &
       & JJ .LT. LY1+2 .OR. JJ .GE. (LY1+LYR-max_offset)) THEN
!    If the subgrid is already set to full grid dimensions then terminate the
!    particle calculation.

     IF(LX1.EQ.1.AND.LY1.EQ.1.AND.LXR.GE.GRID(KG,KT)%NX       &
                             .AND.LYR.GE.GRID(KG,KT)%NY) THEN
        IF(GRID(KG,KT)%GBLDAT)THEN
           OFFG=.TRUE.
           GRID(KG,:)%GLOBAL=.TRUE.  
           WRITE(KF21,*)' NOTICE metsub: switch to global grid'
        END IF
!       KRET=1 (termination not required, particles terminate from metpos)
        RETURN
     END IF

!    Check to see if particle near the full grid edge.  In that situation
!    no subgrid adjustment is required (14 Mar 2003 - RRD) unless the full 
!    grid is global at which point the subgrid should be expanded (17 Sep 2003) 

	
     IF(((II .LT. 3 .OR. II .GE. GRID(KG,KT)%NX-max_offset) .AND. &
          & (JJ .LT. 3 .OR. JJ .GE. GRID(KG,KT)%NY-max_offset)) &
        .AND. .NOT. GRID(KG,KT)%GBLDAT) RETURN	

!    When the data load flag is false (set with each entry to advrng), then
!    meteorological data have not been reloaded since the last call to advrng,
!    which determines the optimal subgrid location for the current particle 
!    distribution.  At this point it is assumed that the  particle is 
!    approaching the edge of the grid and a new subgrid position
!    is required.  The flag is then set to true, which if other particles
!    approach the edge of the grid, the subgrid is now expanded in size to
!    avoid multiple subgrid reloads before the next call to advrng

     IF(.NOT.GRID(KG,KT)%DATLOAD)THEN

!       The first time that an off grid particle is found after a call to
!       advrng, datload is false, and the subgrid is shifted in position and/or
!       expanded depending upon advrng results stored in structure variable.  
!       Subsequent offgrid particles, before another call to advrng result 
!       in a subgrid expansion.

        EXPAND=.FALSE.
        GRID(KG,KT)%DATLOAD=.TRUE. 

!       Values derived from advrng and where the maximum range is the maximum
!       for any defined grid to avoid reloading data into memory when the grid 
!       switches and the dimensions are not the same. Note that corner point
!       is not optimized when grid dimensions set to larger value, but it
!       will include all computational points.  

        LX1=GRID(KG,KT)%LX1
        LY1=GRID(KG,KT)%LY1
        LXR=GRID(KG,KT)%LXR
        LYR=GRID(KG,KT)%LYR
        WRITE(KF21,*)' NOTICE metsub: (kg ,xyr,xy1) - ',KG,LXR,LYR,LX1,LY1,' shift'

!       after grid shift the point is still offgrid then expand (5/12/2005)
        IF(II.LT.LX1+2.OR.II.GE.(LX1+LXR-max_offset).OR.JJ.LT.LY1+2.OR.JJ.GE.(LY1+LYR-max_offset)) &
           EXPAND=.TRUE. 

     ELSE
!       second pass into the routine go straight to grid expansion
        EXPAND=.TRUE.
     END IF

     IF(EXPAND)THEN
       ! prevend an inifitive loop  (tk 2010.09.30)
        if (lxr .lt. 2) then
            lxr=2
            endif
        if (lyr .lt. 2) then
            lyr=2
           endif
  
        keep_expanding=.TRUE.
        do while(keep_expanding)
!       compute center of existing subgrid
        LXC=FLOAT(LX1)+FLOAT(LXR)/2.0
        LYC=FLOAT(LY1)+FLOAT(LYR)/2.0

!       subgrid loaded from previous entry now requires larger subgrid
        LXR=MIN(NINT(LXR*1.5),GRID(KG,KT)%NX)
        LYR=MIN(NINT(LYR*1.5),GRID(KG,KT)%NY)
 
!       Compute sub-grid corner (in full grid units) based on current position
        IF(LXR.EQ.GRID(KG,KT)%NX)THEN
!          subgrid set to full grid dimensions
           LX1 = 1
        ELSE
!          compute new corner position of shifted (and perhaps) expanded subgrid
           LX1 = MIN(MAX(1,NINT(LXC-FLOAT(LXR)/2.0)),GRID(KG,KT)%NX)
        END IF

        IF(LYR.EQ.GRID(KG,KT)%NY)THEN
           LY1 = 1
        ELSE
           LY1 = MIN(MAX(1,NINT(LYC-FLOAT(LYR)/2.0)),GRID(KG,KT)%NY)
        END IF

        WRITE(KF21,*)' NOTICE metsub: (kg ,xyr,xy1) - ',KG,LXR,LYR,LX1,LY1,' expand'
        IF(LXR.GE.GRID(KG,KT)%NX) THEN
           LX1=1
           LXR=GRID(KG,KT)%NX
        END IF
        IF(LYR.GE.GRID(KG,KT)%NY) THEN
           LY1=1
           LYR=GRID(KG,KT)%NY
        END IF
!     check high end to avoid subgrid compression (RRD - 13 Aug 2002)
        IF(LX1+LXR.GT.GRID(KG,KT)%NX) LX1=MAX(1,GRID(KG,KT)%NX-LXR)
        IF(LY1+LYR.GT.GRID(KG,KT)%NY) LY1=MAX(1,GRID(KG,KT)%NY-LYR)

        keep_expanding=(II.LT.LX1+2.OR.II.GE.(LX1+LXR-max_offset).OR.JJ.LT.LY1+2.OR.JJ.GE.(LY1+LYR-max_offset)) &
        .and. .not.(lxr >= GRID(KG,KT)%NX.and. lyr >= GRID(KG,KT)%NY)
     enddo
     END IF

!    if global data available and subgrid bounds on edge then switch to global
     IF(GRID(KG,KT)%GBLDAT)THEN
        IF(LX1.EQ.1.OR.LY1.EQ.1.OR.       &
           LX1+LXR.GE.GRID(KG,KT)%NX.OR.  &
           LY1+LYR.GE.GRID(KG,KT)%NY)  THEN   

           LX1=1
           LY1=1
           LXR=GRID(KG,KT)%NX
           LYR=GRID(KG,KT)%NY
           GRID(KG,:)%GLOBAL=.TRUE.  
           WRITE(KF21,*)' NOTICE metsub: switch to global grid'
        END IF
     END IF
!    set flag   
     OFFG=.TRUE.
  END IF

! Section below is required in the situation where a particle shifts
! to a new grid but is still on the data subgrid. In this case it is
! possible for the subgrid to exceed the new full data grid bounds.  

! no subgrid requirements for large array space reset corner point
! added range set (26 Mar 2003 - RRD) 
  IF(LXR.GE.GRID(KG,KT)%NX) THEN
     LX1=1
     LXR=GRID(KG,KT)%NX
  END IF
  IF(LYR.GE.GRID(KG,KT)%NY) THEN
     LY1=1
     LYR=GRID(KG,KT)%NY
  END IF

! check high end to avoid subgrid compression (RRD - 13 Aug 2002)
  IF(LX1+LXR.GT.GRID(KG,KT)%NX) LX1=MAX(1,GRID(KG,KT)%NX-LXR)
  IF(LY1+LYR.GT.GRID(KG,KT)%NY) LY1=MAX(1,GRID(KG,KT)%NY-LYR)

  IF(OFFG)THEN
!    new values available to advrng
     GRID(KG,KT)%LX1=LX1
     GRID(KG,KT)%LY1=LY1
     GRID(KG,KT)%LXR=LXR
     GRID(KG,KT)%LYR=LYR
  END IF

END SUBROUTINE metsub
