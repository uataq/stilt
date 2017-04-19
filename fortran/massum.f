!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  MASSUM           MASSsUMmation puffs or particles
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:05-05-25
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   MASS SUMMATION FOR PUFFS OR PARTICLES IS CALLED EACH TIME STEP
!   THE MASS OF PARTICLES OR PUFFS ARE SUMMED TO A CELL ACCORDING 
!   TO ITS CENTER POSITION.  NO INTERPOLATION IS PERFORMED AND THE
!   TIME STEP SHOULD BE SUFFICIENTLY SMALL WITH RESPECT TO THE
!   MINIMUM CONCENTRATION GRID CELL SIZE. ROUTINE CALLED FOR EACH
!   PARTICLE OR PUFF AFTER ITS NEW ADVECTION POSITION IS COMPUTED.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 25 May 2005 (RRD) - orignal version from consum
!                 10 Jan 2007 (RRD) - support SNAP=2 option
!                 17 Oct 2007 (RRD) - revised start/stop test 
!
! USAGE:  CALL MASSUM(CONC,NUMGRD,PLAT,PLON,DT,JET,ZMDL,ZSFC,KMSL,
!              CGSIZE,MASS,ZPOS,PTYP,CSUM)
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

SUBROUTINE MASSUM(conc,numgrd,plat,plon,dt,jet,zmdl,zsfc,kmsl,                &
                  cgsize,mass,zpos,ptyp,csum)

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
  REAL,    INTENT(IN)    :: zpos              ! puff center height (sigma)
  INTEGER, INTENT(IN)    :: ptyp              ! pollutant type index number
  REAL,    INTENT(INOUT) :: csum (:,:,:,:,:)  ! concentration summation matrix

!-------------------------------------------------------------------------------
! internal variable definitions
!-------------------------------------------------------------------------------

  INTEGER     :: ksb,ksd
  REAL        :: zpar,zbot,ztop,sfact,glon,dely
  INTEGER     :: kl,kk,nlvl,maxdim,kt,kg,ii,jj

!-------------------------------------------------------------------------------

  KT=PTYP                     ! pollutant type (over-ride if MAXDIM>1)
  MAXDIM = SIZE (mass,1)      ! number of species per single particle

!-------------------------------------------------------------------------------
! main loop through the number of concentration grids 
!-------------------------------------------------------------------------------

  gloop : DO KG=1,NUMGRD

!    sample start stop
     IF(INT(DT).GT.0)THEN
        KSB=JET-CONC(KG)%START%MACC
     ELSE
        KSB=CONC(KG)%START%MACC-JET
     END IF
     KSD=ABS(CONC(KG)%STOP%MACC-CONC(KG)%START%MACC)

!    test for time within sampling interval
     IF(KSB.LT.0.OR.KSB.GE.KSD) CYCLE gloop

     IF(CONC(KG)%SNAP.EQ.0)THEN
!       factors for integrations per sampling interval
        SFACT=ABS(DT)/CONC(KG)%DELTA%MACC
     ELSEIF(CONC(KG)%SNAP.EQ.1)THEN
!       no factor required for snapshot maps
        SFACT=1.0
     ELSEIF(CONC(KG)%SNAP.EQ.2)THEN
!       maximum value set in calling program
        CYCLE gloop
     ELSEIF(CONC(KG)%SNAP.LT.0)THEN
!       shorter averaging time options
        SFACT=ABS(DT)/ABS(CONC(KG)%SNAP)
     ELSE
        SFACT=ABS(DT)/CONC(KG)%DELTA%MACC
     END IF

!    compute longitude domain to test for dateline span
     DELY=CONC(KG)%DELT_LON*CONC(KG)%NUMB_LON

!    longitude correction to avoid dateline problems (RRD - 10/20/99)
     IF(CONC(KG)%X1Y1_LON.GE.0.0.AND.PLON.LT.0.0)THEN
       GLON=PLON+360.0
!    large grids may span dateline (RRD - 06/21/02) 
     ELSEIF((DELY.GT.180.0).AND.(PLON-CONC(KG)%X1Y1_LON.LE.0.0))THEN
        GLON=PLON+360.0
     ELSE
        GLON=PLON
     END IF

!    find the particle position on current concentration grid
     II=NINT(1.0+(GLON-CONC(KG)%X1Y1_LON)/CONC(KG)%DELT_LON)
     JJ=NINT(1.0+(PLAT-CONC(KG)%X1Y1_LAT)/CONC(KG)%DELT_LAT)
     IF(JJ.GT.CONC(KG)%NUMB_LAT.OR.JJ.LT.1) CYCLE gloop
     IF(II.GT.CONC(KG)%NUMB_LON.OR.II.LT.1) CYCLE gloop

!    keep track of active minimum grid size for delt-t calculation
     CGSIZE=MIN(CONC(KG)%SIZE,CGSIZE)

!    vertical position      
     ZPAR=(ZMDL-ZSFC)*(1.0-ZPOS) ! convert particle position to meters
     ZPAR=MAX(ZPAR,1.0)          ! correction 12 Apr 2002
     IF(KMSL.EQ.1)ZPAR=ZPAR+ZSFC ! adjust vertical units MSL rather than AGL

!-------------------------------------------------------------------------------
!    second nested loop for species
!-------------------------------------------------------------------------------

     sloop: DO KK=1,MAXDIM

!       grid pollutant species (kt) must match particle pollutant species (kk)
!       which defaults to kt=ptyp when multiple species are not defined
        IF(MAXDIM.GT.1)KT=KK 

!-------------------------------------------------------------------------------
!       third nested loop for levels  
!-------------------------------------------------------------------------------

        ZBOT=0.0
        NLVL=CONC(KG)%LEVELS    
        lloop : DO KL=1,NLVL
 
!          determine vertical cell sizes (input defines top)
           ZTOP=FLOAT(CONC(KG)%HEIGHT(KL))

!          particle falls within vertical extent of cell
           IF(ZPAR.GT.ZBOT.AND.ZPAR.LE.ZTOP)THEN
!             air concentration option
              CSUM(II,JJ,KL,KT,KG) = CSUM(II,JJ,KL,KT,KG)+MASS(KK)*SFACT
              EXIT lloop
           END IF

           ZBOT=ZTOP
        END DO lloop

     END DO sloop
  END DO gloop

END SUBROUTINE massum
