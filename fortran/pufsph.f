!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PUFSPH           horizontal PUFf SPLitting
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PUFF SPLITTING HORIZONTAL OCCURS WHEN PUFF EXCEEDS PREDEFINED LIMIT
!   THEN SPLIT INTO SMALLER UNITS.  GAUSSIAN PUFFS ARE SPLIT INTO
!   FIVE PUFFS, WHILE TOP-HAT PUFFS ALWAYS SPLIT INTO FOUR PUFFS.
!   MASS IS EQUALLY DIVIDED FOR TOP-HAT PUFFS WHILE FOR THE GAUSSIAN
!   PUFF 60% OF THE MASS IS PUT INTO THE CENTRAL PUFF AND 10% GOES
!   INTO EACH OF THE FOUR CORNER PUFFS.  THE NEW PUFF POSITIONS ARE
!   SET BY 0.5 OF THE HORIZONTAL EXTENT (RADIUS).
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 11 Jun 1997 (RRD)
!                 10 Oct 2000 (RRD) - fortran90 upgrade
!                 09 Sep 2002 (RRD) - fortran coding standards
!                 15 Sep 2003 (RRD) - more focused test on hdwp
!                 02 Apr 2004 (RRD) - generic file unit numbers
!                 12 Aug 2004 (RRD) - variable name change
!                 18 Jan 2008 (RRD) - added splitting shutdown message
!                 28 Jan 2008 (RRD) - dynamica split distance
!
! USAGE:  CALL PUFSPH(KPM,HGD,MASS,XPOS,YPOS,ZPOS,SIGH,SIGW,HDWP,PAGE,PTYP, 
!                     PGRD,NSORT,NUMPAR,CGSIZE,SPLITF,KRET)
!
!   INPUT ARGUMENT LIST:
!   OUTPUT ARGUMENT LIST:
!   INPUT FILES:
!   OUTPUT FILES:
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE PUFSPH(KPM,HGD,MASS,XPOS,YPOS,ZPOS,SIGH,SIGW,HDWP,PAGE,PTYP,   &
                  PGRD,NSORT,NUMPAR,CGSIZE,SPLITF,KRET)

  USE funits

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  INTEGER,  INTENT(INOUT) :: kpm          ! total number of puffs or particles
  REAL,     INTENT(IN)    :: hgd   (:)    ! meteo grid spacing (km)
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
  INTEGER,  INTENT(IN)    :: numpar       ! number of particle permitted 
  REAL,     INTENT(IN)    :: cgsize       ! minimum size (conc/meteo) in km
  REAL,     INTENT(IN)    :: splitf       ! horizontal splitting adjustment
  INTEGER,  INTENT(OUT)   :: kret         ! status of split 0=yes 1=no 2=stop

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

!                    horizontal split distance factors (fraction of extent)
  REAL, PARAMETER :: xf (5) = (/0.5,0.0,-0.5,0.0,0.0/) 
  REAL, PARAMETER :: yf (5) = (/0.0,-0.5,0.0,0.5,0.0/)

  INTEGER         :: kk,mm,kp,kpt,maxdim,nsplit,hdwpx
  REAL            :: hgdm,radius,frac,sigr

!-------------------------------------------------------------------------------

  MAXDIM = SIZE(mass,1)
  KPT=KPM

! determine status of splitting (5 = max value for NSPLIT)
  IF((KPT+5).LT.NUMPAR)THEN
!    room in the array for at least one split
     KRET=0
  ELSE
!    particle array at limit
     KRET=1
     RETURN
  END IF

! go through particles again to split if needed
  DO KP=1,KPM

     HDWPX=HDWP(KP)                            ! complex mode
     IF(HDWPX.GE.100)HDWPX=MOD(HDWP(KP)/10,10) ! simple mode

!    check for on-grid and puff particle
     IF(PGRD(KP).GT.0.AND.(HDWPX.GT.0.AND.HDWPX.LT.5))THEN

!       distribution type determines split procedure
        IF(HDWPX.EQ.1.OR.HDWPX.EQ.3)THEN
!          gaussian
           NSPLIT=5
           SIGR=3.0
           FRAC=0.10
        ELSEIF(HDWPX.EQ.2.OR.HDWPX.EQ.4)THEN
!          top-hat
           NSPLIT=4
           SIGR=1.54
           FRAC=0.25
        ELSE
           WRITE(KF21,*)'*ERROR* pufsph: Invalid horizontal type - ', HDWPX
           WRITE(*,*)   '*ERROR* pufsph: see message file for more information'
           STOP 900
        END IF
!       horizontal extent
        RADIUS=SIGR*SIGH(KP)

!       check for horizontal split (when extent = 1.54 grid size)
        HGDM=MAX(CGSIZE,HGD(PGRD(KP)))*ABS(SPLITF)*1000.0
        IF(RADIUS.GE.(1.54*HGDM).AND.(KPT+NSPLIT).LT.NUMPAR)THEN

           DO KK=1,NSPLIT
!             split into new puffs offset by the sigma
              XPOS(KK+KPT)=XPOS(KP)+XF(KK)*RADIUS/HGDM
              YPOS(KK+KPT)=YPOS(KP)+YF(KK)*RADIUS/HGDM

!             vertical position and time remain the same
              ZPOS(KK+KPT)=ZPOS(KP)

!             puffs that split are permitted to grow again
              PAGE(KK+KPT)=ABS(PAGE(KP))

!             horizontal sigma reduced by half
              SIGH(KK+KPT)=SIGH(KP)/2.0

!             vertical sigma is used as a marker for puffs
!             using mixed mode dispersion
              IF(HDWPX.EQ.3.OR.HDWPX.EQ.4)THEN
!                2d puffs reset vertical turbulence after split
                 SIGW(KK+KPT)=0.0
              ELSE
!                3d puffs or particles maintain value
                 SIGW(KK+KPT)=SIGW(KP)
              END IF

!             distribution, pollutant, and meteo grid
              HDWP(KK+KPT)=HDWP(KP)
              PTYP(KK+KPT)=PTYP(KP)
              PGRD(KK+KPT)=PGRD(KP)

!             save index value until new sort
              NSORT(KK+KPT)=KK+KPM

!             central point (5th) of Gaussian gets 60% of the mass
              IF(KK.EQ.5)FRAC=0.60

!             mass is reduced depending upon distribution
              MASS(1,KK+KPT)=FRAC*MASS(1,KP)
              MM=MAXDIM
              DO WHILE(MM.GT.1)
                 MASS(MM,KK+KPT)=FRAC*MASS(MM,KP)
                 MM=MM-1
              END DO
           END DO
           KPT=KPT+NSPLIT

!          set meteo grid to (0) for old elements
           PGRD(KP)=0

        ELSEIF((KPT+NSPLIT).GE.NUMPAR)THEN
!          special return code for splitting shutdown within loop
           KRET=2
           PAGE(KP)=-ABS(PAGE(KP))
        END IF

!    particle loop
     END IF
  END DO 

! reset particle counter
  KPM=KPT

END SUBROUTINE pufsph
