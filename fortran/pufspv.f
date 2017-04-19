!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PUFSPV           vertical PUFf SPLitting
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PUFF SPLITTING VERTICAL OCCURS WHEN PUFF EXCEEDS VERTICAL LIMITS
!   IS THEN SPLIT INTO SMALLER UNITS.  VERTICAL DISTRIBUTION IS ALWAYS
!   ASSUMED TO BE UNIFORM (TOP-HAT). NUMBER OF SPLITS IS DETERMINED
!   THE VERTICAL PUFF STANDARD DEVIATION IS THEN EVENLY DIVIDED BY THE
!   NUMBER AND THE NEW VERTICAL POSITION IS AT THE CENTRAL POINT.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 14 Feb 1997 (RRD)
!                 10 Oct 2000 (RRD) - fortran90 upgrade
!                 09 Sep 2002 (RRD) - fortran coding standards
!                 12 Aug 2004 (RRD) - variable name change
!                 18 Jan 2008 (RRD) - added splitting shutdown message
!
! USAGE:  CALL PUFSPV(KPM,NLVL,ZSG,MASS,XPOS,YPOS,ZPOS,SIGH,SIGW,HDWP,    
!                     PAGE,PTYP,PGRD,NSORT,NUMPAR,KRET)
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

SUBROUTINE PUFSPV(KPM,NLVL,ZSG,MASS,XPOS,YPOS,ZPOS,SIGH,SIGW,HDWP,PAGE,  &
                  PTYP,PGRD,NSORT,NUMPAR,KRET)

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  INTEGER,  INTENT(INOUT) :: kpm          ! total number of puffs or particles
  INTEGER,  INTENT(IN)    :: nlvl         ! number of levels in subgrid
  REAL,     INTENT(IN)    :: zsg   (:)    ! internal model sigma levels
  REAL,     INTENT(INOUT) :: mass  (:,:)  ! mass of pollutant (arbitrary units)
  REAL,     INTENT(INOUT) :: xpos  (:)    ! center positions (grid units)
  REAL,     INTENT(INOUT) :: ypos  (:)    ! center positions (grid units)
  REAL,     INTENT(INOUT) :: zpos  (:)    ! puff center height (sigma)
  REAL,     INTENT(INOUT) :: sigh  (:)    ! horiz sigma (sigma)
  REAL,     INTENT(INOUT) :: sigw  (:)    ! vert sigma (sigma)
  INTEGER,  INTENT(INOUT) :: hdwp  (:)    ! Horizontal distribution pollutant
  INTEGER,  INTENT(INOUT) :: page  (:)    ! pollutant age since release (min)
  INTEGER,  INTENT(INOUT) :: ptyp  (:)    ! pollutant type index number
  INTEGER,  INTENT(INOUT) :: pgrd  (:)    ! meteorological grid puff position
  INTEGER,  INTENT(INOUT) :: nsort (:)    ! sorted array index values
  INTEGER,  INTENT(IN)    :: numpar       ! number of particles in calculation
  INTEGER,  INTENT(OUT)   :: kret         ! status of split 0=yes 1=no 2=stop   

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  LOGICAL         :: split         ! logical flag to determine split
  REAL, PARAMETER :: sigr  = 1.54  ! vertical extent

  INTEGER         :: mm,kk,kp,kb,kt,kpt,khrs,nsplit,maxdim,hdwpx
  REAL            :: sgt,sgb,frac,znew,zbot,ztop,radius,sigma,depth 

!-------------------------------------------------------------------------------

  MAXDIM = SIZE(mass,1)  ! number of pollutants on single particle
  KPT=KPM

! determine status of splitting (NLVL = max value for NSPLIT)
  IF((KPT+NLVL).LT.NUMPAR)THEN
!    room in the array for at least one split
     KRET=0
  ELSE
!    particle array at limit
     KRET=1
     RETURN
  END IF

! go through particles again to split if needed
  DO KP=1,KPM

!    split at less frequent intervals after 24h
     SPLIT=.FALSE.
     KHRS=ABS(PAGE(KP))/60
     IF(KHRS.LE.24.OR.MOD(KHRS,3).EQ.0)SPLIT=.TRUE.

     HDWPX=HDWP(KP)                            ! complex mode
     IF(HDWPX.GE.100)HDWPX=MOD(HDWP(KP)/10,10) ! simple mode

!    check for time, valid puff on grid, and not vertical particle
     IF(SPLIT.AND.PGRD(KP).GT.0.AND.(HDWPX.EQ.1.OR.HDWPX.EQ.2))THEN

!       scan extent
        RADIUS=SIGR*SIGW(KP)

!       position at bottom and top of puff
        SGT=MAX(ZPOS(KP)-RADIUS, 0.0)
        SGB=MIN(ZPOS(KP)+RADIUS, 1.0)

!       compute vertical index
        KB=1
        KT=NLVL
        DO KK=1,(NLVL-1)
           IF(ZSG(KB+1).GE.SGB)KB=KB+1
           IF(ZSG(KT-1).LE.SGT)KT=KT-1
        END DO

!       check for vertical split (when extent = grid size) and
        NSPLIT=KT-KB
        IF(NSPLIT.GE.3.AND.(KPT+NSPLIT).LE.NUMPAR)THEN

!          top-hat new vertical sigma and layers
           DEPTH=(SGB-SGT)/NSPLIT
           SIGMA=DEPTH/SIGR/2.0
           FRAC=1.0/NSPLIT

!          new layer extent
           ZBOT=SGB
           ZTOP=ZBOT-DEPTH

!          simple bifurcation split across the old sigma
           DO KK=1,NSPLIT
              XPOS(KK+KPT)=XPOS(KP)
              YPOS(KK+KPT)=YPOS(KP)

!             new vertical position
              ZNEW=0.5*(ZBOT+ZTOP)
              ZPOS(KK+KPT)=MAX(MIN(ZNEW,1.0),0.0)

              PAGE(KK+KPT)=ABS(PAGE(KP))
              SIGH(KK+KPT)=SIGH(KP)
              SIGW(KK+KPT)=SIGMA
              HDWP(KK+KPT)=HDWP(KP)
              PTYP(KK+KPT)=PTYP(KP)
              PGRD(KK+KPT)=PGRD(KP)
              NSORT(KK+KPT)=KK+KPT

!             mass is reduced to fraction of previous total
              MASS(1,KK+KPT)=FRAC*MASS(1,KP)
              MM=MAXDIM            
              DO WHILE(MM.GT.1)
                 MASS(MM,KK+KPT)=FRAC*MASS(MM,KP)
                 MM=MM-1
              END DO

!             update layer
              ZBOT=ZTOP
              ZTOP=ZBOT-DEPTH

           END DO
           KPT=KPT+NSPLIT

!          set meteo grid to (0) for old elements
           PGRD(KP)=0

        ELSEIF((KPT+NSPLIT).GT.NUMPAR)THEN
!          special return code for splitting shutdown within loop
           KRET=2
           PAGE(KP)=-ABS(PAGE(KP))

!       vertical split test
        END IF

!    particle loop
     END IF 

  END DO 

! reset particle counter
  KPM=KPT

END SUBROUTINE pufspv
