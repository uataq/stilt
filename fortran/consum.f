!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  CONSUM           CONcentration SUMmation puffs or particles
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   CONCENTRATION SUMMATION VALUES FOR PUFFS OR PARTICLES IS CALLED
!   EACH TIME STEP. PARTICLE CONCENTRATIONS ARE COMPUTED FROM THE
!   CELL VOLUME WHILE PUFF CONCENTRATIONS ARE COMPUTED FROM THE
!   PUFF VOLUME OVER EACH CELL CENTER POINT. THE MAXIMUM VALUE OR
!   THE SUMMATION IS RETAINED OVER THE AVERAGING INTERVAL. NO
!   INTERPOLATION IS PERFORMED AND THE TIME STEP SHOULD BE SUFFICIENTLY 
!   SMALL WITH RESPECT TO THE MINIMUM CONCENTRATION GRID CELL SIZE. THIS
!   ROUTINE CALLED FOR EACH PARTICLE OR PUFF AFTER ITS NEW ADVECTION 
!   POSITION IS COMPUTED.
!
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 24 Nov 1997 (RRD)
!                 20 Apr 1999 (RRD) - added terrain compression adjustment
!                 20 Oct 1999 (RRD) - correction for dateline
!                 08 Nov 1999 (RRD) - grid scan limits test simplified
!                 14 Apr 2000 (RRD) - generalized start-stop test
!                 05 Jul 2000 (RRD) - msl or agl option
!                 21 Sep 2000 (RRD) - fortran90 upgrade
!                 21 Feb 2001 (RRD) - further longitude grid correction
!                 02 Nov 2001 (RRD) - shortened averaging time options
!                 11 Jan 2002 (RRD) - time setp grid size set in loop 
!                 12 Apr 2002 (RRD) - set minimum particle height
!                 21 Jun 2002 (RRD) - refined dateline span test
!                 24 Jun 2002 (RRD) - converted coordinate from Z* to Z
!                 09 Sep 2002 (RRD) - fortran coding standards
!                 15 Sep 2003 (RRD) - surface particle option (hdwp=5)
!                 02 Apr 2004 (RRD) - generic file unit numbers
!                 12 Aug 2004 (RRD) - variable name change
!                 12 Oct 2005 (RRD) - lagrangian sampling option test
!                 17 Jan 2007 (RRD) - maximum value option SNAP=2
!                 17 Oct 2007 (RRD) - revised start/stop test  
!
! USAGE:  CALL CONSUM(CONC,NUMGRD,PLAT,PLON,DT,JET,ZMDL,ZSFC,KMSL,
!              CGSIZE,MASS,DEPT,ZPOS,SIGH,SIGW,HDWP,PTYP,CSUM)
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

SUBROUTINE CONSUM(conc,numgrd,plat,plon,dt,jet,zmdl,zsfc,kmsl,                &
                  cgsize,mass,dept,zpos,sigh,sigw,hdwp,ptyp,csum)

  USE funits

  IMPLICIT NONE 

  INCLUDE 'DEFCONC.INC'             ! pollutant and concentration grid

!-------------------------------------------------------------------------------
! argument list definitions
!-------------------------------------------------------------------------------

  TYPE(cset), INTENT(IN) :: conc(:)           ! for each concentration grid 
  INTEGER, INTENT(IN)    :: numgrd            ! number of concentration grids
  REAL,    INTENT(IN)    :: plat,plon         ! particle position lat & long
  REAL,    INTENT(IN)    :: dt                ! time step (min)
  INTEGER, INTENT(IN)    :: jet               ! elapsed time (min)
  REAL,    INTENT(IN)    :: zmdl              ! model domain top
  REAL,    INTENT(IN)    :: zsfc              ! height of ground surface (m)
  INTEGER, INTENT(IN)    :: kmsl              ! agl=0 or msl=1 height units
  REAL,    INTENT(INOUT) :: cgsize            ! minimum sampling grid size (km)
  REAL,    INTENT(IN)    :: mass (:)          ! mass of pollutant 
  REAL,    INTENT(IN)    :: dept (:)          ! deposition amount of mass
  REAL,    INTENT(IN)    :: zpos              ! puff center height (sigma)
  REAL,    INTENT(IN)    :: sigh,sigw         ! horiz (m) and vert sigma (sigma)
  INTEGER, INTENT(IN)    :: hdwp              ! horiz distribution in pollutant
  INTEGER, INTENT(IN)    :: ptyp              ! pollutant type index number
  REAL,    INTENT(INOUT) :: csum (:,:,:,:,:)  ! concentration summation matrix

!-------------------------------------------------------------------------------
! internal variable definitions
!-------------------------------------------------------------------------------

  REAL,   PARAMETER      :: sigr   = 1.54         ! distribution
  REAL,   PARAMETER      :: ympd   = 111198.323   ! meters/deg 
  REAL,   PARAMETER      :: twopi  = 6.283185     ! Gaussian constants (2 PI)
  REAL,   PARAMETER      :: PI     = 3.14159265358979
  REAL,   PARAMETER      :: DGPR   = 180.0/PI  ! deg per radian

  INTEGER                :: ksb,ksd
  REAL                   :: xmpd,zpar,radius
  INTEGER                :: kt,kg,ii,jj,i1,i2,j1,j2
  REAL                   :: dfxy0,dfxyz,zbot,ztop
  REAL                   :: vdist,sfact,glon,gfact
  REAL                   :: xp,yp,xi1,yj1,xi2,yj2
  REAL                   :: cval,clat,clon,delx,dely,dist
  INTEGER                :: kl,kk,nlvl,maxdim
  REAL                   :: pbot,ptop

!-------------------------------------------------------------------------------

  IF(HDWP.EQ.6) RETURN        ! 10/12/05 lagrangian sampling option uses CONPAR
  XMPD=YMPD*COS(PLAT/DGPR)    ! meters per degree longitude
  KT=PTYP                     ! pollutant type (over-ride if MAXDIM>1)
  MAXDIM = SIZE (mass,1)      ! number of species per single particle

!-------------------------------------------------------------------------------
! main loop through the number of concentration grids and set grid factors
!-------------------------------------------------------------------------------

  gloop : DO KG=1,NUMGRD

!   sample start stop 
    IF(INT(DT).GT.0)THEN
       KSB=JET-CONC(KG)%START%MACC
    ELSE
       KSB=CONC(KG)%START%MACC-JET
    END IF
    KSD=ABS(CONC(KG)%STOP%MACC-CONC(KG)%START%MACC)

!   test for time within sampling interval
    IF(KSB.LT.0.OR.KSB.GE.KSD) CYCLE gloop

    IF(CONC(KG)%SNAP.EQ.0)THEN
!      factors for integrations per sampling interval
       SFACT=ABS(DT)/CONC(KG)%DELTA%MACC
    ELSEIF(CONC(KG)%SNAP.EQ.1)THEN
!      no factor required for snapshot maps
       SFACT=1.0
    ELSEIF(CONC(KG)%SNAP.EQ.2)THEN
!      maximum value set in calling program
       CYCLE gloop
    ELSEIF(CONC(KG)%SNAP.LT.0)THEN
!      shorter averaging time options
       SFACT=ABS(DT)/ABS(CONC(KG)%SNAP)
    ELSE
!      unsupported option default to averaging
       SFACT=ABS(DT)/CONC(KG)%DELTA%MACC
    END IF

!   compute longitude domain to test for dateline span
    DELY=CONC(KG)%DELT_LON*CONC(KG)%NUMB_LON

!   longitude correction to avoid dateline problems (RRD - 10/20/99)
    IF(CONC(KG)%X1Y1_LON.GE.0.0.AND.PLON.LT.0.0)THEN
       GLON=PLON+360.0
!   large grids may span dateline (RRD - 06/21/02) 
    ELSEIF((DELY.GT.180.0).AND.(PLON-CONC(KG)%X1Y1_LON.LE.0.0))THEN
       GLON=PLON+360.0
    ELSE
       GLON=PLON
    END IF

!-------------------------------------------------------------------------------
!   horizontal plume factors
!-------------------------------------------------------------------------------

    IF(HDWP.EQ.0.OR.HDWP.EQ.5)THEN
       RADIUS=1.0E+25   ! particle model - no distribution
!      normalized dispersion factor = cell area
       DFXY0=1.0/(XMPD*CONC(KG)%DELT_LON*YMPD*CONC(KG)%DELT_LAT)

    ELSEIF(HDWP.EQ.1.OR.HDWP.EQ.3)THEN
       RADIUS=3.0*SIGH  ! horizontal gaussian model - scan to 3 sigma
       DFXY0=1.0/(TWOPI*SIGH*SIGH) ! dispersion factor for area deposition

    ELSEIF(HDWP.EQ.2.OR.HDWP.EQ.4)THEN
       RADIUS=SIGR*SIGH ! top hat distribution - scan 1.54 sigma
       DFXY0=1.0/(PI*RADIUS*RADIUS)

    ELSE
       WRITE(KF21,*)'*ERROR* consum: Invalid horizontal type - ', HDWP
       WRITE(*,*)   '*ERROR* consum: see message file for more information'
       STOP 900
    END IF

!-------------------------------------------------------------------------------
!   vertical plume factors
!-------------------------------------------------------------------------------

    IF(HDWP.EQ.1.OR.HDWP.EQ.2)THEN
!      plume vertical extent
       PBOT=(ZMDL-ZSFC)*(1.0-MAX(ZPOS+SIGR*SIGW,0.0))
       PTOP=(ZMDL-ZSFC)*(1.0-MIN(ZPOS-SIGR*SIGW,1.0))
!      vertical distribution uniform over layer
       VDIST=MAX(PTOP-PBOT,1.0)  
       DFXYZ=DFXY0/VDIST           ! uniform concentration within volume
!      adjust for vertical units MSL rather than AGL (rrd - 7/5/00)
       IF(KMSL.EQ.1)THEN
          PBOT=PBOT+ZSFC
          PTOP=PTOP+ZSFC
       END IF

    ELSEIF(HDWP.EQ.0.OR.HDWP.EQ.3.OR.HDWP.EQ.4)THEN
       ZPAR=(ZMDL-ZSFC)*(1.0-ZPOS) ! convert particle position to meters
       ZPAR=MAX(ZPAR,1.0)          ! correction 12 Apr 2002
       IF(KMSL.EQ.1)ZPAR=ZPAR+ZSFC ! adjust vertical units MSL rather than AGL

    ELSEIF(HDWP.EQ.5)THEN
       ZPAR=0
       IF(KMSL.EQ.1)ZPAR=ZSFC 

    ELSE
       WRITE(KF21,*)'*ERROR* consum: Invalid horizontal type - ', HDWP
       WRITE(*,*)   '*ERROR* consum: see message file for more information'
       STOP 900
    END IF

!-------------------------------------------------------------------------------
!   find the grid position on current concentration grid
!-------------------------------------------------------------------------------

    XP=1.0+(GLON-CONC(KG)%X1Y1_LON)/CONC(KG)%DELT_LON
    YP=1.0+(PLAT-CONC(KG)%X1Y1_LAT)/CONC(KG)%DELT_LAT

!   for particles (no distribution) insure no scan
    IF(HDWP.EQ.0.OR.HDWP.EQ.5)THEN
!      round position so that cell centered over point
       XI1=FLOAT(NINT(XP))
       YJ1=FLOAT(NINT(YP))
       XI2=XI1
       YJ2=YJ1

    ELSE
!      find lower left corner based upon scan radius
       CLAT=PLAT-RADIUS/YMPD
       CLON=GLON-RADIUS/XMPD

!      convert corner to grid index number
       XI1=1.0+(CLON-CONC(KG)%X1Y1_LON)/CONC(KG)%DELT_LON
       YJ1=1.0+(CLAT-CONC(KG)%X1Y1_LAT)/CONC(KG)%DELT_LAT

!      opposite corner from puff position for scan
       XI2=2.0*XP-XI1
       YJ2=2.0*YP-YJ1
    END IF

!   convert real grid to nearest integer value
    I1=INT(XI1)
    I2=NINT(XI2)
    J1=INT(YJ1)
    J2=NINT(YJ2)

!-------------------------------------------------------------------------------
! second nested loop for species
!-------------------------------------------------------------------------------

  sloop: DO KK=1,MAXDIM

!   grid pollutant species (kt) must match particle pollutant species (kk)
    IF(MAXDIM.GT.1)KT=KK  ! defaults (kt=ptyp) if multiple species not defined

!-------------------------------------------------------------------------------
! third nested loop for levels  
!-------------------------------------------------------------------------------

  ZBOT=0.0
  NLVL=CONC(KG)%LEVELS               ! sum through all levels within plume
  lloop : DO KL=1,NLVL

!   determine vertical cell sizes (input defines top)
    ZTOP=FLOAT(CONC(KG)%HEIGHT(KL))
    VDIST=MAX(ZTOP-ZBOT,1.0)         ! particle type

!-------------------------------------------------------------------------------
! loop through grid (sampling) points within plume
!-------------------------------------------------------------------------------

  jloop : DO JJ=J1,J2
     IF(JJ.GT.CONC(KG)%NUMB_LAT.OR.JJ.LT.1) CYCLE jloop
  iloop : DO II=I1,I2
     IF(II.GT.CONC(KG)%NUMB_LON.OR.II.LT.1) CYCLE iloop

!    keep track of active minimum for delt-t calculation
     CGSIZE=MIN(CONC(KG)%SIZE,CGSIZE)

!    distance (meters) of sampling point to plume center
     CLON=(II-1.0)*CONC(KG)%DELT_LON+CONC(KG)%X1Y1_LON
     CLAT=(JJ-1.0)*CONC(KG)%DELT_LAT+CONC(KG)%X1Y1_LAT
     DELY=(PLAT-CLAT)*YMPD
     DELX=(GLON-CLON)*XMPD
     DIST=SQRT(DELY*DELY+DELX*DELX)

!    grid point within plume radius
     IF(DIST.GE.RADIUS) CYCLE iloop

!-------------------------------------------------------------------------------
!    vertical particle distributions
!-------------------------------------------------------------------------------

     IF(HDWP.EQ.0.OR.HDWP.GE.3)THEN

        GFACT=1.0
!       horizontal gaussian - vertical particle
        IF(HDWP.EQ.3)GFACT=EXP(-0.5*DIST*DIST/SIGH/SIGH)

!       deposition output requested
        IF(ZBOT.EQ.0.AND.ZTOP.EQ.0.0)THEN
           CSUM(II,JJ,KL,KT,KG) = CSUM(II,JJ,KL,KT,KG)+DEPT(KK)*DFXY0*GFACT

!       particle falls within vertical extent of cell
        ELSEIF(ZPAR.GT.ZBOT.AND.ZPAR.LE.ZTOP)THEN
!          air concentration option
           CVAL=MASS(KK)*SFACT*DFXY0*GFACT/VDIST
           CSUM(II,JJ,KL,KT,KG) = CSUM(II,JJ,KL,KT,KG) + CVAL 

        ELSE
           CYCLE iloop
        END IF

!-------------------------------------------------------------------------------
!    vertical puff distribution 
!-------------------------------------------------------------------------------

     ELSE

        GFACT=1.0
!       compute horizontal gaussian factor if required
        IF(HDWP.EQ.1)GFACT=EXP(-0.5*DIST*DIST/SIGH/SIGH)

!       deposition output requested (level=0 defined)
        IF(ZTOP.EQ.0.0)THEN
           CSUM(II,JJ,KL,KT,KG)=CSUM(II,JJ,KL,KT,KG)+DEPT(KK)*DFXY0*GFACT

!       output grid node within vertical plume
        ELSEIF(ZTOP.GE.PBOT.AND.ZTOP.LT.PTOP)THEN
           CVAL=MASS(KK)*SFACT*DFXYZ*GFACT
           CSUM(II,JJ,KL,KT,KG)= CSUM(II,JJ,KL,KT,KG) + CVAL 

        ELSE
           CYCLE iloop
        END IF

!   distribution test
    END IF

!-------------------------------------------------------------------------------

  END DO iloop
  END DO jloop

  ZBOT=ZTOP
  END DO lloop

  END DO sloop
  END DO gloop

END SUBROUTINE consum
