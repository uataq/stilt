!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PARINP           PARticle INPut from file
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:05-08-03
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   READS A FILE OF POLLUTANTSD PARTICLE POSITIONS WHICH IS USED 
!   TO INITIALIZE THE MODEL.  THE FILE MAY HAVE BEEN GENERATED
!   FROM A PREVIOUS SIMULATION OR CREATED FROM MEASUREMENT DATA.
!   THE FILE MAY CONTAIN ONE OR MORE TIME PERIODS.  WHEN THERE
!   ARE MULTIPLE TIME PERIODS, THESE PARTICLES ARE ADDED TO THE
!   SIMULATION WHEN THE MODEL TIME MATCHES THE VALID TIME OF THE 
!   PARTICLE ARRAY.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 08 May 2003 (RRD) - Initial version
!                 02 Apr 2004 (RRD) - Generic file unit numbers
!                 10 Aug 2004 (RRD) - Argument list change
!                 13 Dec 2004 (RRD) - Extended error testing
!                 12 May 2005 (RRD) - Dropped single grid restriction
!                 25 Jan 2007 (RRD) - additional iostat test
!
! USAGE:  CALL PARINP(JET,KG,KT,KPM,MASS,TLAT,TLON,XPOS,YPOS,
!              ZPOS,SIGH,SIGU,SIGV,SIGW,HDWP,PAGE,PTYP,PGRD,NSORT,
!              ZMDL,MAXPAR,NINIT)
!
!   INPUT ARGUMENT LIST:    see below
!   OUTPUT ARGUMENT LIST:   see below
!   INPUT FILES:            none
!   OUTPUT FILES:           none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE PARINP(JET,KG,KT,KPM,MASS,TLAT,TLON,XPOS,YPOS,ZPOS,         &
                  SIGH,SIGU,SIGV,SIGW,HDWP,PAGE,PTYP,PGRD,NSORT,       &
                  ZMDL,MAXPAR,NINIT)

  USE funits
  use module_defgrid

  IMPLICIT NONE


!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  INTEGER, INTENT(IN)    :: jet        ! current elapsed time
  INTEGER, INTENT(IN)    :: kg         ! current meteo grid in space
  INTEGER, INTENT(IN)    :: kt         ! current meteo grid in time
  INTEGER, INTENT(INOUT) :: kpm        ! total number of puffs or particles
  REAL,    INTENT(INOUT) :: mass (:,:) ! mass of pollutant (arbitrary units)
  REAL,    INTENT(INOUT) :: tlat (:)   ! puff center positions (lat)
  REAL,    INTENT(INOUT) :: tlon (:)   ! puff center positions (lon)
  REAL,    INTENT(INOUT) :: xpos (:)   ! puff center positions (grid units)
  REAL,    INTENT(INOUT) :: ypos (:)   ! puff center positions (grid units)
  REAL,    INTENT(INOUT) :: zpos (:)   ! puff center height (sigma)
  REAL,    INTENT(INOUT) :: sigh (:)   ! horizontal puff sigma 
  REAL,    INTENT(INOUT) :: sigu (:)   ! turbulence u'2     
  REAL,    INTENT(INOUT) :: sigv (:)   ! turbulence v'2
  REAL,    INTENT(INOUT) :: sigw (:)   ! turbulence w'2 or vertical puff sigma
  INTEGER, INTENT(INOUT) :: hdwp (:)   ! Horizontal distribution pollutant
  INTEGER, INTENT(INOUT) :: page (:)   ! pollutant age since release (min)
  INTEGER, INTENT(INOUT) :: ptyp (:)   ! pollutant type index number
  INTEGER, INTENT(INOUT) :: pgrd (:)   ! meteorological grid of puff position
  INTEGER, INTENT(INOUT) :: nsort(:)   ! sort index array by position
  REAL,    INTENT(IN)    :: zmdl       ! vertical index scaling height
  INTEGER, INTENT(IN)    :: maxpar     ! maximum number of particles  
  INTEGER, INTENT(IN)    :: ninit      ! particle initialization type 

!-------------------------------------------------------------------------------

  LOGICAL  :: LIMIT = .FALSE. ! exceeding particle number limit

  REAL     :: ZTMP       ! particle height
  INTEGER  :: KP,KPAR    ! particle counter
  INTEGER  :: KRET       ! IO return code
  INTEGER  :: NUMTYP     ! number of different pollutants
  INTEGER  :: NUMPOL     ! number of different pollutants
  INTEGER  :: IYR        ! file year
  INTEGER  :: IMO        ! file month
  INTEGER  :: IDA        ! file day 
  INTEGER  :: IHR        ! file hour 
  INTEGER  :: IMN        ! file minutes
  INTEGER  :: MTIME      ! file time accumulated minutes
  INTEGER  :: I,J        ! local indicies


  SAVE LIMIT

!-------------------------------------------------------------------------------

! set pollutant dimension
  NUMTYP=SIZE(mass,1)

  IMN=0
  KRET=0
  DO WHILE(KRET.EQ.0)

!    check for valid header record for time period
     READ (KF23,IOSTAT=KRET)KPAR,NUMPOL,IYR,IMO,IDA,IHR
     IF(KRET.NE.0)THEN
!       end of file or error, shutdown further reads
        KPAR=0
        CLOSE(KF23)
        WRITE(KF21,*)'WARNING parinp: particle initialization file closed'
     ELSE
        CALL TM2MIN(IYR,IMO,IDA,IHR,IMN,MTIME)
!       reset and wait until time matches 
        IF(MTIME.GT.JET)THEN
           WRITE(KF21,*)'WARNING parinp: initialization in time wait mode'
           BACKSPACE(KF23)
           RETURN
        END IF
     END IF

!    depending upon particle initialization procedure, set pointer
!    ninit=2 adds new particles to existing
!    ninit=3 replaces existing particles with new particles
     IF(MTIME.EQ.JET.AND.KPAR.GT.0.AND.NINIT.EQ.3) KPM=0

     DO J=1,KPAR
        KP=KPM+J
        IF(KP.GT.MAXPAR)THEN
           KP=MAXPAR
           IF(.NOT.LIMIT)THEN
              WRITE(KF21,*)'WARNING parinp: exceeding particle number limit'
              LIMIT=.TRUE.
           END IF
        END IF

        READ (KF23,IOSTAT=KRET)(MASS(I,KP),I=1,NUMPOL)
        IF(KRET.NE.0)THEN
           KPAR=0
           CLOSE(KF23)
           WRITE(KF21,*)'WARNING parinp: error in particle file record #1'
        END IF
        IF(NUMPOL.NE.NUMTYP)THEN
           WRITE(KF21,*)'*ERROR* parinp: mismatch in number of pollutants'
           WRITE(KF21,*)'Mode: ',NUMTYP,'     File: ',NUMPOL
           STOP 900
        END IF

!       data records 
!       additional variables not yet supported in supplemental programs
!       READ (KF23,IOSTAT=KRET)TLAT(KP),TLON(KP),ZTMP,SIGH(KP),SIGW(KP),SIGV(KP),SIGU(KP)

        READ (KF23,IOSTAT=KRET)TLAT(KP),TLON(KP),ZTMP,SIGH(KP),SIGW(KP),SIGV(KP)
        IF(KRET.NE.0)THEN
           KPAR=0
           CLOSE(KF23)
           WRITE(KF21,*)'WARNING parinp: error in particle file record #2'
        END IF

        READ (KF23,IOSTAT=KRET)PAGE(KP),HDWP(KP),PTYP(KP),PGRD(KP),NSORT(KP)
        IF(KRET.NE.0)THEN
           KPAR=0
           CLOSE(KF23)
           WRITE(KF21,*)'WARNING parinp: error in particle file record #3'
        END IF

!       additional variables not yet supported in supplemental programs
        SIGU(KP)=SIGV(KP)

!       process if time period matches and grid ID is not missing
        IF(MTIME.EQ.JET)THEN            
        IF(PGRD(KP).NE.0)THEN

! disabled, but run requires same meteo files as run that created particles
!          force meteo grid index back to default start grid
! 5/12/05  PGRD(KP)=KG

!          convert positions from lat/lon to internal grid 
           IF(GRID(PGRD(KP),KT)%LATLON)THEN
              CALL GBL2XY(PGRD(KP),KT,TLAT(KP),TLON(KP),               &
                          XPOS(KP),YPOS(KP))
           ELSE
              CALL CLL2XY_wps(GRID(PGRD(KP),KT)%GBASE,TLAT(KP),TLON(KP),   &
                          XPOS(KP),YPOS(KP),GRID(PGRD(KP),KT)%proj)
           END IF

!          convert height to sigma assuming terrain height = 0
           ZPOS(KP)=(ZMDL-ZTMP)/ZMDL

! disabled, but run requires same meteo files as run that created particles
!          flag any puffs that might be offgrid, or in the previous run
!          may have created particles on a different grid
! 5/12/05  IF(XPOS(KP).LE.1.0.OR.XPOS(KP).GE.FLOAT(GRID(KG,KT)%NX).OR.      &
! 5/12/05     YPOS(KP).LE.1.0.OR.YPOS(KP).GE.FLOAT(GRID(KG,KT)%NY))PGRD(KP)=0

        END IF
        END IF
     END DO 

     IF(MTIME.EQ.JET)THEN
!       set particle count to last valid input
        KPM=KPM+KPAR
        WRITE(KF21,*)' NOTICE parinp: particle initialization - ',IYR,IMO,IDA,IHR
        WRITE(KF21,*)' Initialized with particles: ',KPM
        KRET=1       ! forces loop exit
     ELSE
!       go on to next record until times match or exceeds current model time
        IF(KRET.EQ.0)  &
        WRITE(KF21,*)'WARNING parinp: particle record skipped - ',IYR,IMO,IDA,IHR
     END IF

  END DO

END SUBROUTINE parinp
