!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  EMSGRD           EMiSsion GRiDded input starts puff/particle
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   EMISSION GRIDDED INPUT - STARTS A NEW PUFF OR PARTICLE FROM POINT
!   DEFINED IN AN INPUT FILE. OTHER RELEASE CRITERIA ARE THE SAME.
!   INPUT FILE READ IN SUBROUTINE EMSINP.
!   EMISSIONS CAN BE STARTED AGAIN AT QCYCLE INTERVALS.
!   THE ROUTINE IS CALLED EACH TIME STEP.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 26 Apr 1998 (RRD)
!                 22 Dec 1998 (RRD) - default comments for chemistry
!                 22 Jul 1999 (RRD) - variable name change
!                 14 Apr 2000 (RRD) - generalized start stop tests
!                 16 Mar 2001 (RRD) - global lat lon grid option
!                 13 Jul 2001 (RRD) - zero mass prior to setting value
!                 04 Oct 2001 (RRD) - multiple simultaneous meteorology
!                 09 Sep 2002 (RRD) - fortran coding standards
!                 22 Jul 2003 (RRD) - test for pollutant ID match
!                 27 Aug 2003 (RRD) - test for volume units conversion
!                 02 Apr 2004 (RRD) - generic file unit numbers
!                 16 Jul 2004 (RRD) - emission on/off message
!                 22 Jul 2004 (RRD) - exceeding puff limit no exit
!                 10 Aug 2004 (RRD) - enhanced particle-puff conversion
!                 24 Oct 2005 (RRD) - emission less than time step
!                 06 Jan 2006 (RRD) - backward emission fix
!                 17 Oct 2007 (RRD) - revised emission start/stop tests 
!                 05 Nov 2007 (RRD) - short time release correction
!                 07 Aug 2008 (RRD) - refined INITD test
!
! USAGE:  CALL EMSGRD(SPOT,DIRT,NUMTYP,KPM,INITD,DT,JET,NSORT,MASS,XPOS,
!              YPOS,ZPOS,SIGH,SIGU,SIGV,SIGW,HDWP,PAGE,PTYP,PGRD,QCYCLE,
!              NUMPAR,NPVAL,NQLON,NQLAT,QDLON,QDLAT,QAREA,POLID,MAXPAR)
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

SUBROUTINE EMSGRD(SPOT,DIRT,NUMTYP,KPM,INITD,DT,JET,NSORT,MASS,XPOS,          &
                  YPOS,ZPOS,SIGH,SIGU,SIGV,SIGW,HDWP,PAGE,PTYP,PGRD,QCYCLE,   &
                  NUMPAR,NPVAL,NQLON,NQLAT,QDLON,QDLAT,QAREA,POLID,MAXPAR)

  USE funits
  use module_defgrid ! meteorology file

  IMPLICIT NONE

  INCLUDE 'DEFCONC.INC' ! pollutant and concentration grid
  INCLUDE 'DEFSPOT.INC' ! multiple source information

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  TYPE(rset),INTENT(IN)    :: spot(:)     ! source location characteristics
  TYPE(pset),INTENT(INOUT) :: dirt(:)     ! for each pollutant type 
  INTEGER,   INTENT(IN)    :: numtyp      ! number of pollutant types
  INTEGER,   INTENT(INOUT) :: kpm         ! number of puffs or particles
  INTEGER,   INTENT(IN)    :: initd       ! initial distribution type
  REAL,      INTENT(IN)    :: dt          ! time step (min)
  INTEGER,   INTENT(IN)    :: jet         ! current elapsed time (min)
  REAL,      INTENT(IN)    :: qcycle      ! optional emission cycle time in hours
  INTEGER,   INTENT(IN)    :: numpar      ! maximum number of particles permitted
  INTEGER,   INTENT(INOUT) :: nsort (:)   ! index of sorted elements
  REAL,      INTENT(INOUT) :: mass  (:,:) ! mass of pollutant (arbitrary units)
  REAL,      INTENT(INOUT) :: xpos  (:)   ! horizontal position (grid units)
  REAL,      INTENT(INOUT) :: ypos  (:)   ! horizontal position (grid units)
  REAL,      INTENT(INOUT) :: zpos  (:)   ! puff center height (sigma)
  REAL,      INTENT(INOUT) :: sigh (:)    ! horizontal puff sigma 
  REAL,      INTENT(INOUT) :: sigu (:)    ! turbulence u'2    
  REAL,      INTENT(INOUT) :: sigv (:)    ! turbulence v'2
  REAL,      INTENT(INOUT) :: sigw (:)    ! turbulence w'2 or vertical puff sigma
  INTEGER,   INTENT(INOUT) :: hdwp  (:)   ! Horizontal distribution pollutant
  INTEGER,   INTENT(INOUT) :: page  (:)   ! pollutant age since release (min)
  INTEGER,   INTENT(INOUT) :: ptyp  (:)   ! pollutant type index number
  INTEGER,   INTENT(INOUT) :: pgrd  (:)   ! meteorological grid of puff position

  INTEGER,      INTENT(IN) :: npval            ! numb pollutants
  INTEGER,      INTENT(IN) :: nqlon, nqlat     ! number of grid points
  REAL,         INTENT(IN) :: qdlon, qdlat     ! grid spacing
  REAL,         INTENT(IN) :: qarea (:,:,:,:)  ! gridded emissions array 
  CHARACTER(4), INTENT(IN) :: polid (:)        ! emissions pollutant id
  INTEGER,      INTENT(IN) :: maxpar           ! maximum number of particles

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  LOGICAL         :: emit
  LOGICAL         :: pmit   = .FALSE.
  REAL, PARAMETER :: sigr   = 1.5       ! top-hat radius
  REAL, PARAMETER :: PI     = 3.14159265358979
  REAL, PARAMETER :: DEGPRD = 180.0/PI  ! deg per radian

  INTEGER  :: iyr,imo,ida,ihr,imn
  INTEGER  :: nrpts,ih,nphr,jj,ii,npar,kk,kcycl,ktp,kd,kp,kg,kt
  REAL     :: yp,area,clat,clon,xp,delz,qtot,qval,qstep,qsum,qmax
  INTEGER  :: ksb,ksd,initk

  REAL,ALLOCATABLE :: tfact(:)

!-------------------------------------------------------------------------------
! external variables
!-------------------------------------------------------------------------------

  SAVE   PMIT

!-------------------------------------------------------------------------------

  INTERFACE 
  SUBROUTINE GBL2XY(KG,KT,CLAT,CLON,X,Y)
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: kg          ! active grid number    
  INTEGER, INTENT(IN)  :: kt          ! active time number    
  REAL,    INTENT(IN)  :: clat,clon   ! latlon location       
  REAL,    INTENT(OUT) :: x,y         ! grid position         
  END SUBROUTINE GBL2XY
  END INTERFACE

!-------------------------------------------------------------------------------
! check to determine if emissions required for any pollutant
!-------------------------------------------------------------------------------

! compute current hour and convert to array index (depends on #q per day)
  CALL TM2DAY(JET,IYR,IMO,IDA,IHR,IMN)
  IH=IHR+1  ! index 1 = average 0000 to 0100 UTC

  IF(.NOT.ALLOCATED(tfact)) ALLOCATE (tfact(numtyp))

  TFACT=1.0
  EMIT=.FALSE.
  DO KK=1,NUMTYP
     IF(INT(DT).GT.0)THEN
        KSB=JET-DIRT(KK)%START%MACC
        KSD=MAX(INT(DIRT(KK)%QHRS*60.0),INT(DT))
     ELSE
        KSB=DIRT(KK)%START%MACC-JET
        KSD=INT(MAX(ABS(DIRT(KK)%QHRS*60.0),ABS(DT)))
     END IF
     IF(KSB.GE.0.AND.KSB.LT.KSD)THEN
        EMIT=.TRUE.
!       adjust for short duration emissions relative to time step
        IF(KSD-KSB.LT.ABS(INT(DT))) TFACT(KK)=(KSD-KSB)/ABS(DT)
     END IF
  END DO

  IF(EMIT.AND..NOT.PMIT)WRITE(KF21,*)' NOTICE emsgrd: emissions started'
  IF(.NOT.EMIT.AND.PMIT)WRITE(KF21,*)' NOTICE emsgrd: emissions terminated'
  PMIT=EMIT
  IF(.NOT.EMIT)RETURN

!-------------------------------------------------------------------------------
! check to see if any of the pollutants in control file match table
!-------------------------------------------------------------------------------

  EMIT=.FALSE.
  DO KK=1,NUMTYP
!    match pollutant to release with emissions table
     DO KD=1,NPVAL
        IF(DIRT(KK)%IDENT.EQ.POLID(KD))THEN
           EMIT=.TRUE.
!###       WRITE(KF21,*)' NOTICE emsgrd: emitting - ',DIRT(KK)%IDENT
        END IF
     END DO
  END DO

  IF(.NOT.EMIT)THEN
     WRITE(KF21,*)'WARNING emsgrd: No matching pollutants found'
     WRITE(KF21,*)'  File content: ',POLID
     WRITE(KF21,*)'  Control file: ',(DIRT(KK)%IDENT,KK=1,NUMTYP)
     STOP 900
  END IF

! number of pollutants on single particle
  IF(SIZE(mass,1).LT.NUMTYP)THEN
     WRITE(KF21,*)'*ERROR* emsgrd: Number of pollutants in CONTROL file'
     WRITE(KF21,*)' exceeds internal array. Increase MAXDIM in SETUP.CFG' 
     WRITE(KF21,*)' namelist file to at least: ',NUMTYP
     WRITE(*,*)   '*ERROR* emsgrd: see message file for more information'
     STOP 900
  END IF

!-------------------------------------------------------------------------------
! set the number of units to emit
!-------------------------------------------------------------------------------

  INITK=INITD
  IF(INITK.GE.100) INITK=MOD(INITD/10,10)

  IF(INITK.EQ.1.OR.INITK.EQ.2)THEN
!    gaussian or top-hat emissions
     NPAR=1
  ELSE
     QMAX=0.0
     NRPTS=NQLAT*NQLON
!    find maximum emission duration
     DO KK=1,NUMTYP
        QMAX=MAX(QMAX,ABS(DT/60.0),ABS(DIRT(KK)%QHRS))
     END DO
!    particle emission rate: number particles per hour per source
     NPHR=INT(FLOAT(NUMPAR/NRPTS)/QMAX)
!    particle emissions per time step
     NPAR=MAX(1, INT(FLOAT(NPHR)*ABS(DT/60.0)))
  END IF

!-------------------------------------------------------------------------------
! loop through the number of independent source locations
!-------------------------------------------------------------------------------

  KT=1                   ! grids the same at all times
  KG=SPOT(1)%KG          ! initial grid index set in main program
  EMIT=.TRUE.

  DO JJ=1,NQLAT
  DO II=1,NQLON

!     check for non-zero emissions
      QSUM=0.0
      DO KK=1,NPVAL
         QSUM=QSUM+QAREA(II,JJ,KK,IH)
      END DO
      IF(QSUM.GT.0.0)THEN

!     compute position from lower left corner point (1,1 point)
!     assume emissions at center of grid cell (from SW corner)
      CLAT=FLOAT(JJ-1)*QDLAT+SPOT(1)%OLAT+QDLAT/2.0
      CLON=FLOAT(II-1)*QDLON+SPOT(1)%OLON+QDLON/2.0

!     convert to grid units
      IF(GRID(KG,KT)%LATLON)THEN
         CALL GBL2XY(KG,KT,CLAT,CLON,XP,YP)
      ELSE
         CALL CLL2XY_wps(GRID(KG,KT)%GBASE,CLAT,CLON,XP,YP,GRID(KG,KT)%proj)
      END IF

!     check for location relative to meteo grid
      IF(XP.GT.2.0.AND.XP.LT.FLOAT(GRID(KG,KT)%NX-1).AND.YP.GT.2.0.AND.   &
         YP.LT.FLOAT(GRID(KG,KT)%NY-1))THEN

!     multiple emissions at each cell only for particles
      ploop : DO KP=1,NPAR

!        increment particle/puff counter
         KPM=KPM+1
         IF(KPM.GT.MAXPAR)THEN
            KPM=MAXPAR
            IF(EMIT)WRITE(KF21,*)'Warning: emsgrd - exceeding puff limit'
            EMIT=.FALSE.
            EXIT ploop
         END IF
         NSORT(KPM)=KPM

!        initial position
         XPOS(KPM)=XP
         YPOS(KPM)=YP

!        set equivalent horizontal sigma regardless of distribution 
!        assume area source (111000 m / deg - squared)
         AREA=1.2E+10*QDLAT*QDLON*COS(CLAT/DEGPRD)
!        compute sigma for uniform radius
         SIGH(KPM)=SQRT(AREA/PI)/SIGR

!        initial depth defined from release height
!        input height assumed to define layer from ground
         DELZ=1.0-SPOT(1)%ZP

!        particles are distributed uniformly in height from the
!        ground to the specified release height spot%zp
         IF(INITK.EQ.1.OR.INITK.EQ.2)THEN
!           puff variance set to layer depth
            SIGW(KPM)=ABS(DELZ/SIGR/2.0)
            ZPOS(KPM)=1.0-DELZ/2.0
         ELSE
            SIGW(KPM)=0.0
!           particles get distributed in the layer
            ZPOS(KPM)=1.0-FLOAT(KP)*DELZ/FLOAT(NPAR)
         END IF

!        horizontal variances start at zero
         SIGU(KPM)=0.0
         SIGV(KPM)=0.0

!        initial distribution (see main for definitions)
         HDWP(KPM)=INITD
!        initial age at zero
         PAGE(KPM)=0
!        variable not used in this configuration (all in single element)
         PTYP(KPM)=1      
!        initial grid is the default startup grid from main
         PGRD(KPM)=KG

!-------------------------------------------------------------------------------
!        loop through pollutants
!-------------------------------------------------------------------------------

         DO KK=1,NUMTYP

!           make sure no residual mass left over from previous use (7/13/01-RRD)
            MASS(KK,KPM)=0.0

!           match pollutant to release with emissions table
            KTP=0
            DO KD=1,NPVAL
               IF(DIRT(KK)%IDENT.EQ.POLID(KD))KTP=KD
            END DO

!           check if this pollutant requires a start
!           gridded emissions are assume to be continuous after start
!           qhrs variable is the time that gives the emission rate
            IF(INT(DT).GT.0)THEN
               KSB=JET-DIRT(KK)%START%MACC
               KSD=MAX(INT(DIRT(KK)%QHRS*60.0),INT(DT))
            ELSE
               KSB=DIRT(KK)%START%MACC-JET
               KSD=INT(MAX(ABS(DIRT(KK)%QHRS*60.0),ABS(DT)))
            END IF

            IF(KSB.GE.0.AND.KSB.LT.KSD)THEN           
 
!              number of time steps in emission period
               QSTEP=MAX(1.0,60.0*ABS(DIRT(KK)%QHRS/DT))
               IF(KTP.NE.0)THEN
!                 total emission (qrate from input becomes a multiplier)
                  QTOT=ABS(DIRT(KK)%QHRS)*QAREA(II,JJ,KTP,IH)*DIRT(KK)%QRATE
               ELSE
!                 set rate as specified in the control file
                  QTOT=ABS(DIRT(KK)%QHRS)*DIRT(KK)%QRATE
               END IF

!              emission per time step
               QVAL=QTOT/QSTEP
!              divide amount over the number of units emitted
               MASS(KK,KPM)=TFACT(KK)*QVAL/NPAR

!              optional conversion to ppm (volume factor not included)
!              assume emission in kg, molecular weight in grams, flag is
!              true when the molecular weight is set to a negative value
               IF(DIRT(KK)%DOVOL)  &
                  MASS(KK,KPM)=22.4E+06*MASS(KK,KPM)/DIRT(KK)%GPMOL

!           start time test
            END IF

!        pollutant type loop
         END DO

!     number of particles loop
      END DO ploop

!     within meteo grid test
      END IF

!     non-zero emissions test
      END IF

! number of sources loop
  END DO
  END DO

!-------------------------------------------------------------------------------
! check for emission cycling
!-------------------------------------------------------------------------------

! to avoid the emission of an excessive number of particles or puffs
! (ie the whole grid each time step), gridded emission simulations should
! be made by cycling the emissions.  This is accomplished by emitting a
! short burst at the cycle interval (qcycle).  The cycle interval is set
! in SETUP.CFG.  For instance if the emission duration were set to 0.1 h,
! then the emission rate multiplier (source term in the CONTROL file)
! should be set to 10.0 if the cycle time is 1 h.  This would then emit the
! correct amount of mass over the 1 hour cycle period.

  DO KK=1,NUMTYP
!    test for end of emission cycle
     IF(INT(DT).GT.0)THEN
        KSB=JET+INT(DT)-DIRT(KK)%START%MACC
        KSD=MAX(INT(DIRT(KK)%QHRS*60.0),INT(DT))
     ELSE
        KSB=DIRT(KK)%START%MACC-JET-INT(DT)
        KSD=INT(MAX(ABS(DIRT(KK)%QHRS*60.0),ABS(DT)))
     END IF

     IF(KSB.GE.KSD)THEN        
!       optional restart of emissions at some later/earlier time
        KCYCL=NINT(SIGN(QCYCLE*60.0,DIRT(KK)%QHRS))
        DIRT(KK)%START%MACC=DIRT(KK)%START%MACC+KCYCL
     END IF
  END DO

END SUBROUTINE emsgrd
