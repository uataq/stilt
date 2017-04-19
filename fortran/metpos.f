!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  METPOS           METeorological POSitioning finds record
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   METEOROLOGICAL POSITIONING ROUTINE USES THE CURRENT CALCULATION
!   POSITION AND TIME TO CHECK IF THE POSITION FALLS WITHIN THE METEO
!   DATA ALREADY LOADED INTO MEMORY. IF NOT IS IS DETERMINED IF DATA
!   ON A NEW GRID OR TIME ARE REQUIRED. THE ROUTINE RETURNS A POSITIVE 
!   RECORD NUMBER IF DATA ARE TO BE READ FROM THE INPUT FILE.
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 15 May 1998 (RRD)
!                  08 Mar 1999 (RRD) - optimized subgrid selection
!                                    - removed grid number dependence
!                  17 Jun 1999 (RRD) - corrected newx flag at first entry
!                  23 Nov 1999 (RRD) - do not set newx at times between files
!                                    - terminate when no more data fix
!                  06 Mar 2000 (RRD) - mtime between files when delta-t differs
!                  14 Aug 2000 (XUE) - permits meteo data at non-synoptic time
!                  22 Sep 2000 (RRD) - fortran90 upgrade
!                  16 Mar 2001 (RRD) - subgrid point selection optimized
!                  23 Apr 2001 (RRD) - global grid testing
!                  23 May 2001 (RRD) - full grid test switched from OR to AND
!                  04 Oct 2001 (RRD) - simultaneous multiple meteorology   
!                  17 Jan 2002 (RRD) - fixed problem with more than 2 grids
!                  09 Sep 2002 (RRD) - switched test from global to gbldat 
!                  02 Apr 2004 (RRD) - generic file unit numbers
!                  18 Oct 2007 (RRD) - improved diagnostic message
!
! USAGE:  CALL METPOS(BACK,XP,YP,JET,NGRD,NTIM,FTIME,MTIME,POINT,OFFG,
!                     KGC,KGX,KT1,KT2)
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

 SUBROUTINE METPOS(BACK,XP,YP,JET,NGRD,NTIM,FTIME,MTIME,POINT,OFFG,   &
                   KGC,KGX,KT1,KT2)

  USE funits
  use module_defgrid ! meteorology grid and file

  IMPLICIT NONE


!-------------------------------------------------------------------------------
! argument list definitions
!-------------------------------------------------------------------------------

  LOGICAL,   INTENT(IN)    :: back           ! defines backward integration
  REAL,      INTENT(INOUT) :: xp,yp          ! x,y grid particle position
  INTEGER,   INTENT(IN)    :: jet            ! elapsed time (minutes)
  INTEGER,   INTENT(IN)    :: ngrd           ! number of meteo grids   
  INTEGER,   INTENT(IN)    :: ntim           ! number of meteo times   
  INTEGER,   INTENT(IN)    :: ftime(:,:)     ! time of current array data

  INTEGER,   INTENT(OUT)   :: mtime (2)      ! time of requested input data
  INTEGER,   INTENT(OUT)   :: point (2)      ! index pointer to array 
  LOGICAL,   INTENT(OUT)   :: offg           ! off entire grid flag
  INTEGER,   INTENT(IN)    :: kgc            ! current position grid number 
  INTEGER,   INTENT(OUT)   :: kgx            ! new position grid number 
  INTEGER,   INTENT(OUT)   :: kt1            ! grid number time last    
  INTEGER,   INTENT(OUT)   :: kt2            ! grid number time next    

!-------------------------------------------------------------------------------
! internal definitions
!-------------------------------------------------------------------------------

  LOGICAL :: diag = .true.
  LOGICAL :: met1, met2
  REAL    :: tlat,tlon,xpt(mgrd),ypt(mgrd)
  INTEGER :: ii,jj,kg,kt,k1,k2
  integer :: max_offset

!-------------------------------------------------------------------------------
! external definitions
!-------------------------------------------------------------------------------

  SAVE DIAG

!-------------------------------------------------------------------------------
! First determine which meteorological grid is the optimum advection grid.
! Assume that the finest grid is always defined as #1 and particles will be 
! advected on that grid if possible. The search is required each entry because 
! a particle may move from a coarser grid to a finer grid.  Time is defined
! such that the particle moves away from index #1 toward index #2.
!-------------------------------------------------------------------------------

  MET1  = .FALSE. ! time flag if outside of temporal domain
  MET2  = .FALSE.
  OFFG  = .TRUE.  ! start with the assumption that particle within domain
  KGX=0           ! grid selection

  gloop : DO KG=NGRD,1,-1

!    initialize time pointers for each grid
     K1=0  
     K2=0

     tloop : DO KT=1,NTIM

!       check if this time period contains a valid grid
        IF(GRID(KG,KT)%NUMBER.LT.0) CYCLE tloop

!       if test grid not equal to current then remap position
        IF(KG.NE.KGC)THEN  

!          for multiple grids convert position from current grid to true
!          then convert true back to new grid coordinates
           IF(GRID(KGC,1)%LATLON.AND.GRID(KG,1)%LATLON)THEN
              CALL GBL2LL(KGC,1,XP,YP,TLAT,TLON)
              CALL GBL2XY(KG,KT,TLAT,TLON,XPT(KG),YPT(KG))
           ELSEIF(     GRID(KGC,1)%LATLON.AND.(.NOT.GRID(KG,1)%LATLON))THEN
              CALL GBL2LL(KGC,1,XP,YP,TLAT,TLON)
              CALL CLL2XY_wps(GRID(KG,KT)%GBASE,TLAT,TLON,XPT(KG),YPT(KG),GRID(KG,KT)%proj)
           ELSEIF(.NOT.GRID(KGC,1)%LATLON.AND.      GRID(KG,1)%LATLON )THEN
              CALL CXY2LL_wps(GRID(KGC,1)%GBASE,XP,YP,TLAT,TLON,GRID(KGC,1)%proj)
              CALL GBL2XY(KG,KT,TLAT,TLON,XPT(KG),YPT(KG))
           ELSEIF(.NOT.GRID(KGC,1)%LATLON.AND.(.NOT.GRID(KG,1)%LATLON))THEN
              CALL CXY2LL_wps(GRID(KGC,1)%GBASE,XP,YP,TLAT,TLON,GRID(KGC,1)%proj)
              CALL CLL2XY_wps(GRID(KG,KT)%GBASE,TLAT,TLON,XPT(KG),YPT(KG),GRID(KG,KT)%proj)
           END IF

!          convert particle position to index units
           II=INT(XPT(KG))
           JJ=INT(YPT(KG))

        ELSE
!          convert particle position to index units
           II=INT(XP)
           JJ=INT(YP)
        END IF


!       particle must be within the spatial domain of one of the files
!       within a 2 grid cell external band
        IF(.NOT.GRID(KG,KT)%GBLDAT)THEN
           max_offset = 2
!       For WRF input, stay away from one additional row/column of grid edge
!       (because rows/columns for non-staggered field are padded)
           if (GRID(KG,kt)%MODEL_ID(2:4) .eq. 'WRF') max_offset = max_offset+1
!cg(20100804) same for ECMWF fields
           if (GRID(KG,kt)%MODEL_ID(1:2) .eq. 'EC') max_offset = max_offset+1
           IF(II .LT. 2 .OR. II .GE. GRID(KG,KT)%NX-max_offset .OR.             &
              JJ .LT. 2 .OR. JJ .GE. GRID(KG,KT)%NY-max_offset)    CYCLE tloop
        END IF
        OFFG=.FALSE.

!       particle must be within temporal domain of one of the files

        IF(.NOT.BACK)THEN
           IF(JET.GE.FILE(KG,KT)%FIRST%MACC.AND.                &
              JET.LT.FILE(KG,KT)%LAST%MACC)THEN
!             data point falls within the file
              K1=KT
              K2=KT
           ELSEIF(JET.GE.FILE(KG,KT)%LAST%MACC.AND.             &
                  JET.LT.FILE(KG,KT)%LAST%MACC+DREC(KG,KT)%DELTA)THEN
!             data point potentially falls after the end of the file
              IF(K1.EQ.0)K1=KT
              MET2=.TRUE.
           ELSEIF(JET.LT.FILE(KG,KT)%FIRST%MACC.AND.            &
                  JET.GE.FILE(KG,KT)%FIRST%MACC-DREC(KG,KT)%DELTA)THEN
!             data point potentially falls before the start of the file
              IF(K2.EQ.0)K2=KT
              MET1=.TRUE.
           ELSE
              IF(JET.LT.FILE(KG,KT)%FIRST%MACC) MET1=.TRUE.
              IF(JET.GT.FILE(KG,KT)%LAST%MACC ) MET2=.TRUE.
           END IF

        ELSE
           IF(JET.LE.FILE(KG,KT)%LAST%MACC.AND.                 &
              JET.GT.FILE(KG,KT)%FIRST%MACC)THEN
!             data point falls within the file
              K1=KT
              K2=KT
           ELSEIF(JET.LE.FILE(KG,KT)%FIRST%MACC.AND.            &
                  JET.GT.FILE(KG,KT)%FIRST%MACC-DREC(KG,KT)%DELTA)THEN
!             data point potentially falls after the end of the file
              IF(K1.EQ.0)K1=KT
              MET2=.TRUE.
           ELSEIF(JET.GT.FILE(KG,KT)%LAST%MACC.AND.             &
                  JET.LE.FILE(KG,KT)%LAST%MACC+DREC(KG,KT)%DELTA)THEN
!             data point potentially falls before the start of the file
              IF(K2.EQ.0)K2=KT
              MET1=.TRUE.
           ELSE
              IF(JET.LT.FILE(KG,KT)%FIRST%MACC) MET2=.TRUE.
              IF(JET.GT.FILE(KG,KT)%LAST%MACC ) MET1=.TRUE.
           END IF
        END IF

     END DO tloop

!    when both times are defined then the grid is valid for computations
     IF(K1.NE.0.AND.K2.NE.0)THEN
        KGX=KG
        KT1=K1
        KT2=K2
     END IF

  END DO gloop

!-------------------------------------------------------------------------------
! The computational point must have been on one of the meteorological grids
! such that the grid number KGX is not equal to zero.
!-------------------------------------------------------------------------------

  IF(KGX.EQ.0)THEN
!    terminate calculation if position no on any computational grid 
     IF(SUM(FTIME).EQ.0) THEN   
!       internal arrays still empty if data had not been previously loaded
        WRITE(*,*)'*ERROR* metpos: start point not within (x,y,t) any data file'
        IF(MET1)WRITE(*,*)' - start time before start of meteorology data'
        IF(MET2)WRITE(*,*)' - start time after end of meteorology data'
        IF(OFFG)WRITE(*,*)' - start location outside of domain - ',XP,YP
        STOP 900  
     ELSE
        IF(DIAG)THEN
!          particle termination diagnostic message occurs only once
!dwen(20090901)           DIAG=.FALSE.

           IF(OFFG)THEN
              WRITE(KF21,*)'WARNING metpos: off spatial domain of all grids'

           ELSE
              IF(NGRD.EQ.1)THEN
                 IF(K1.EQ.K2.OR.NTIM.EQ.1)THEN
                    WRITE(KF21,*)  &
                   'WARNING metpos: no meteo data at current time - ',JET
                 ELSE
                    WRITE(KF21,*)  &
                   'WARNING metpos: too much time between files - ',JET
                 END IF
              ELSE
!                multiple grids not sure which grid has the time constraint 
                 WRITE(KF21,*)     &
                'WARNING metpos: no data at current time on any grid - ',JET
              END IF
           END IF
        END IF

        OFFG=.TRUE.
        RETURN
     END IF

  ELSEIF(KGX.NE.KGC)THEN
!    position remapped to the new grid's coordinate system
     XP=XPT(KGX)
     YP=YPT(KGX)
  END IF

! advection point is on a meteorological grid
  OFFG=.FALSE.

!-------------------------------------------------------------------------------
! Compute the meteo data times required for interpolation based upon the time
! of the computational point.  These times may already be loaded in memory and
! then a new meteorological data read and processing is not required.
!-------------------------------------------------------------------------------

  IF(BACK)THEN
     MTIME(1)=FILE(KGX,KT1)%LAST%MACC-INT((FILE(KGX,KT1)%LAST%MACC-JET)/    &
              DREC(KGX,KT1)%DELTA)*DREC(KGX,KT1)%DELTA
     MTIME(2)=MTIME(1)-DREC(KGX,KT1)%DELTA

  ELSE
     MTIME(1)=FILE(KGX,KT1)%FIRST%MACC+INT((JET-FILE(KGX,KT1)%FIRST%MACC)/  &
              DREC(KGX,KT1)%DELTA)*DREC(KGX,KT1)%DELTA
     MTIME(2)=MTIME(1)+DREC(KGX,KT1)%DELTA
  END IF

!-------------------------------------------------------------------------------
! Determine data array pointers such that when the request time #2 reaches the
! end of the interpolation time period such that it equals the loaded time #2
! then the array pointer is switched and new data would be loaded into the 
! other array element, alternating between #1 and #2
!-------------------------------------------------------------------------------

  POINT=0  

! note that point(1) indicates the first data time array element
! and point(2) indicates the last time data array

  IF(MTIME(1).EQ.FTIME(KGX,1)) POINT(1)=1
  IF(MTIME(2).EQ.FTIME(KGX,2)) POINT(2)=2
  IF(MTIME(1).EQ.FTIME(KGX,2)) POINT(1)=2
  IF(MTIME(2).EQ.FTIME(KGX,1)) POINT(2)=1

  IF(POINT(1).EQ.0.AND.POINT(2).EQ.0)THEN
!    at initial time no data are yet loaded
     POINT(1)=1
     POINT(2)=2
  ELSEIF(POINT(1).EQ.0)THEN
!    valid data in second element, new data loaded into first element
     POINT(1)=MOD(POINT(2),2)+1
  ELSEIF(POINT(2).EQ.0)THEN
!    valid data in first element, new data loaded into second element
     POINT(2)=MOD(POINT(1),2)+1
  END IF     

!-------------------------------------------------------------------------------
! Print diagnostics when times do not match
!-------------------------------------------------------------------------------

  IF(MTIME(1).NE.FTIME(KGX,POINT(1))) THEN
     IF(KT1.NE.KT2) WRITE(KF21,*)' NOTICE metpos: (kgx,kt1,kt2) - ',KGX,KT1,KT2
     WRITE(KF21,*)' NOTICE metpos: (mtime,ftime) - ',MTIME(1),FTIME(KGX,POINT(1))
  END IF

  IF(MTIME(2).NE.FTIME(KGX,POINT(2))) THEN
     IF(KT1.NE.KT2) WRITE(KF21,*)' NOTICE metpos: (kgx,kt1,kt2) - ',KGX,KT1,KT2
     WRITE(KF21,*)' NOTICE metpos: (mtime,ftime) - ',MTIME(2),FTIME(KGX,POINT(2))
  END IF

END SUBROUTINE metpos
