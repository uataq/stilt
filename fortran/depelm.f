!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  DEPELM           DEPosition of a pollutant ELeMent
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   DEPOSITION OF A POLLUTANT ELEMENT COMPUTES GRAVITATIONAL SETTLING,
!   DRY DEPOSITION EITHER EXPLICIT OR VIA THE RESISTANCE METHOD,
!   WET REMOVAL, AND RADIOACTIVE DECAY AS APPLIED TO ONE POLLUTANT
!   PARTICLE OR PUFF EACH TIME STEP.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 17 Nov 1997 (RRD)
!                 27 Oct 1998 (MDC) - Below cloud scavenging correction
!                                   - Added explicit gravitational settling
!                 20 Apr 1999 (RRD) - terrain compression factor
!                 18 Jun 1999 (RRD) - consistent gas definition
!                 03 Sep 2000 (RRD) - fortran90 upgrade
!                 09 Sep 2002 (RRD) - fortran coding standards
!                 27 Aug 2003 (RRD) - resistance method returns dry dep vel
!                 15 Sep 2003 (RRD) - compute probability of dry depositon
!                 12 Aug 2004 (RRD) - variable name change
!                 12 Oct 2005 (RRD) - lagrangian sampling option test
!                 17 Oct 2007 (RRD) - positive time step forced
!
! USAGE:  CALL DEPELM(DIRT,OLAT,IBMO,NLVL,ZX,DT,ZMDL,ZSFC,ZSG,MASS,DEPT,ZPOS,
!              SIGW,ICHEM,KTYP,LAND,ROUG,SFCL,USTR,PSI,SFLX,HDWP,RAIN,
!              DD,TT,QQ,DRYD)
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

SUBROUTINE DEPELM(DIRT,OLAT,IBMO,NLVL,ZX,DT,ZMDL,ZSFC,ZSG,MASS,DEPT,ZPOS,  &
                  SIGW,ICHEM,KTYP,LAND,ROUG,SFCL,USTR,PSI,SFLX,HDWP,RAIN,  &
                  DD,TT,QQ,DRYD)

  IMPLICIT NONE

  INCLUDE 'DEFCONC.INC'         ! concentration and pollutant structure

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  TYPE(pset), INTENT(IN)   :: dirt (:) ! for each pollutant type 
  REAL,     INTENT(IN)     :: olat     ! origin location
  INTEGER,  INTENT(IN)     :: ibmo     ! computational month
  INTEGER,  INTENT(IN)     :: nlvl     ! number of vertical levels
  REAL,     INTENT(IN)     :: zx       ! vertical array position (index)
  REAL,     INTENT(IN)     :: dt       ! time step (min)
  REAL,     INTENT(IN)     :: zmdl     ! model domain top (meters)
  REAL,     INTENT(IN)     :: zsfc     ! height of terrain surface (m)
  REAL,     INTENT(IN)     :: zsg  (:) ! model sigma levels
  REAL,     INTENT(INOUT)  :: mass (:) ! mass (arbitrary units)
  REAL,     INTENT(INOUT)  :: zpos     ! puff center height (sigma)
  REAL,     INTENT(INOUT)  :: sigw     ! vert sigma (sigma)
  INTEGER,  INTENT(IN)     :: ichem    ! special depositon options   
  INTEGER,  INTENT(IN)     :: ktyp     ! pollutant type index number
  INTEGER,  INTENT(IN)     :: land     ! land use category (1-11)
  REAL,     INTENT(IN)     :: roug     ! aerodynamic roughness length (m)
  REAL,     INTENT(IN)     :: sfcl     ! height of the surface layer (m)
  REAL,     INTENT(IN)     :: ustr     ! friction velocity (m/s)
  REAL,     INTENT(IN)     :: psi      ! integrated stability function for heat
  REAL,     INTENT(IN)     :: sflx     ! incident short wave flux (w/m2)
  INTEGER,  INTENT(INOUT)  :: hdwp     ! pollutant distribution type (index)
  REAL,     INTENT(IN)     :: rain     ! precipitation value (m/min)
  REAL,     INTENT(IN)     :: dd   (:) ! air density profile (kg/m3)
  REAL,     INTENT(IN)     :: tt   (:) ! temperature profile
  REAL,     INTENT(IN)     :: qq   (:) ! humidity profile (fraction 0-1)
  REAL,     INTENT(OUT)    :: dept (:) ! deposition total (mass units)
  REAL,     INTENT(OUT)    :: dryd (:) ! deposition velocity as computed

!-------------------------------------------------------------------------------
! internally defined variables
!-------------------------------------------------------------------------------

  REAL,  PARAMETER :: grav  = 9.801        ! GRAVITY (m/s2) 
  REAL,  PARAMETER :: dmvc  = 1.789E-02    ! DYNAMIC VISCOSITY (g m-1 s-1)    
  REAL,  PARAMETER :: frep  = 6.53E-08     ! MEAN FREE PATH (m at stp)
  REAL,  PARAMETER :: dstp  = 1.2E+03      ! STP DENSITY (g/m3)     
  REAL,  PARAMETER :: rgas  = 0.082057     ! gas constant (atm-liter / deg-mole)
  REAL,  PARAMETER :: sigr  = 1.54         ! vertical puff scan factor
 
  INTEGER          :: k,kt,kk,krh,klvl,kbot,ktop,maxdim
  REAL             :: sc,vb,vd,vg,fracb,frbct,rtc,depv,zlvl,beta,frea
  REAL             :: drop,aird,dens,pdia,pbot,ptop,pdepth,cbot,ctop,rvalue

!-------------------------------------------------------------------------------

  INTERFACE
  SUBROUTINE DEPDRY(DIRT,OLAT,IBMO,KPOL,LAND,ROUG,SFCL,USTR,PSI,SFLX,AIRD,   &
                    TEMP,PDIA,VG,VD)
  IMPLICIT NONE
  INCLUDE 'DEFCONC.INC'         ! concentration and pollutant structure
  TYPE(pset), INTENT(IN)  :: dirt(:)    ! for each pollutant type 
  REAL,     INTENT(IN)    :: olat       ! origin latitude
  INTEGER,  INTENT(IN)    :: ibmo       ! computational month
  INTEGER,  INTENT(IN)    :: kpol       ! polluant index number
  INTEGER,  INTENT(IN)    :: land       ! land-use category
  REAL,     INTENT(IN)    :: roug       ! aerodynamic roughness length (m)
  REAL,     INTENT(IN)    :: sfcl       ! height of constant flux layer (m)
  REAL,     INTENT(IN)    :: ustr       ! friction velocity (m/s)
  REAL,     INTENT(IN)    :: psi        ! integrated stability function heat
  REAL,     INTENT(IN)    :: sflx       ! solar irradiation at sfc (watts/m2)
  REAL,     INTENT(IN)    :: aird       ! ambient air density (g/m3)
  REAL,     INTENT(IN)    :: temp       ! canopy air temperature (deg K)
  REAL,     INTENT(IN)    :: pdia       ! particle diameter (meters)
  REAL,     INTENT(IN)    :: vg         ! gravitational settling velocity (m/s)
  REAL,     INTENT(OUT)   :: vd         ! total deposition velocity (m/s)
  END SUBROUTINE depdry
  END INTERFACE

!-------------------------------------------------------------------------------

  IF(HDWP.EQ.5) RETURN   ! deposited particles cannot deposit again
  IF(HDWP.EQ.6) RETURN   ! 10/12/05 lagrangian sampling option uses CONPAR

  MAXDIM = SIZE(mass,1)  ! number of pollutants on single particle
 
! rounded vertical index position for meteorology profiles
  KLVL=NINT(ZX)

  IF(HDWP.EQ.1.OR.HDWP.EQ.2)THEN
!    height value at bottom and top of puff in meters
     PBOT=(ZMDL-ZSFC)*(1.0- MIN(ZPOS+SIGR*SIGW, 1.0)  )
     PTOP=(ZMDL-ZSFC)*(1.0- MAX(ZPOS-SIGR*SIGW, 0.0)  )
     PDEPTH=MAX(PTOP-PBOT,1.0)
  ELSE
!    particle distribution for vertical
     PBOT=(ZMDL-ZSFC)*(1.0-ZPOS)
     PTOP=PBOT
     KK=MIN(INT(ZX)+1,NLVL)
!    assume depth equal to meteorological cell size
     PDEPTH=(ZMDL-ZSFC)*(ZSG(KK-1)-ZSG(KK))
     PDEPTH=MAX(SFCL,PDEPTH)
  END IF

! set default pollutant type (over-ride if MAXDIM>1)
  KT=KTYP

!-------------------------------------------------------------------------------
! check for simultaneous species at position
!-------------------------------------------------------------------------------

  polnum : DO KK=1,MAXDIM

!    use default if multiple species not defined
     IF(MAXDIM.GT.1)KT=KK

     IF(DIRT(KT)%DODRY)THEN
!       explicit definition of the dry deposition velocity
        VD=DIRT(KT)%DRYVL
!       set gravitational settling if defined as particle
        IF(DIRT(KT)%DOGAS)THEN
           VG=0.0
        ELSE
           VG=VD
        END IF

     ELSEIF(DIRT(KT)%DOGRV.OR.DIRT(KT)%DORES)THEN
!       local air density
        AIRD=DD(KLVL)
!       convert (kg/m3) --> (g/m3)
        AIRD=AIRD*1000.0

!       compute gravitational settling for particles
        IF(DIRT(KT)%DOGRV)THEN
!          particle density from g/cc g/m3
           DENS=DIRT(KT)%PDENS*1.0E+06
!          particle diameter (um) to (m)
           PDIA=DIRT(KT)%PDIAM*1.0E-06
!          base settling velocity
           VB=PDIA*PDIA*GRAV*(DENS-AIRD)/(18.0*DMVC)
!          slip correction (mean free path = particle diameter)
           FREA=FREP*(DSTP/AIRD)
           SC=1.0+(2.0*FREA/PDIA)*(1.26+0.4*EXP(-0.55*PDIA/FREA))
!          final value apply shape correction factor
           VG=VB*SC/DIRT(KT)%SHAPE
        ELSE
           VG=0.0
        END IF

        IF(DIRT(KT)%DORES)THEN
!          compute resistance based VD
!**************************
!dwen(20090818):HYSPLIT4.5 removed DIRT from COMMON bloc and added to arguement list
           CALL DEPDRY(DIRT,OLAT,IBMO,KT,LAND,ROUG,SFCL,USTR,PSI,SFLX, &
                       AIRD,TT(1),PDIA,VG,VD)
        ELSE
!          without resistance computation assume Vd = settling
           VD=VG
        END IF

     ELSE
        VD=0.0
        VG=0.0
     END IF
     DRYD(KT)=VD

!    change in vertical position due to settling
     IF(VG.GT.0.0)THEN
        DROP=60.0*VG*ABS(DT)
        PTOP=MAX(0.0,PTOP-DROP)
        PBOT=MAX(0.0,PBOT-DROP)
        ZPOS=1.0-0.5*(PTOP+PBOT)/(ZMDL-ZSFC)
!       adjust values for puffs - don't overwrite for particles
        IF(HDWP.EQ.1.OR.HDWP.EQ.2)THEN
           SIGW=(PTOP-PBOT)/SIGR/2.0/(ZMDL-ZSFC)
           PDEPTH=MAX(PTOP-PBOT,1.0)
        END IF
     END IF

!    zero out removal rate constant (1/min)
     BETA=0.0

!    if puff within first layer compute dry removal time constant
     IF(VD.GT.0.0.AND.PBOT.LT.SFCL)THEN
        BETA=60.0*VD/PDEPTH

        IF((ICHEM.EQ.5.OR.ICHEM.EQ.7).AND.(BETA.GT.0.0))THEN
!          The probability deposition option implies that the particle
!          stays on the surface and drops all its mass rather than
!          just losing a fraction of its mass. The probability that
!          the particle will deposit on the surface is determined to
!          be true if a random number (0-1) is less than beta*dt. This
!          can only be used for single pollutant particles (maxdim=1)

           CALL RANDOM_NUMBER(RVALUE)

!          apply exponential for removal > 1%
           IF(ABS(DT)*BETA.GE.0.01) BETA=(1.0-EXP(-ABS(DT)*BETA))/ABS(DT)

           IF(RVALUE.LT.BETA*ABS(DT))THEN
!             random value within dry deposit range
              IF(ICHEM.EQ.7.AND.LAND.EQ.7)THEN
!                convert to deposited particle over water surfaces
!                and hold on to mass for future transport
                 HDWP=5 
              ELSE
!                over solid surfaces particles drop mass
                 DEPT(KK)=MASS(KK)
                 MASS(KK)=0.0                  
              END IF
              EXIT polnum
           ELSE
!             if no dry deposit then it may still wet deposit
              BETA=0.0
           END IF
        END IF

     END IF

!    test for wet removal processes
     IF(DIRT(KT)%DOWET.AND.RAIN.GT.0.0)THEN

!       determine bottom and top of the precip layer (80% to 60%)
        KBOT=0
        KTOP=NLVL
        DO K=1,NLVL
           KRH=INT(QQ(K)*100.0+0.5)
           IF(KBOT.EQ.0.AND.KRH.GE.80)KBOT=K
           IF(KBOT.NE.0.AND.KTOP.EQ.NLVL.AND.KRH.LE.60)KTOP=K
        END DO

!       default layer (to be consistent with rain even if no cloud)
        IF(KBOT.EQ.0)THEN
           DO K=1,NLVL
              ZLVL=(ZMDL-ZSFC)*(1.0-ZSG(K))
              IF(ZLVL.LT. 500.0)KBOT=K
              IF(ZLVL.LT.3000.0)KTOP=K
           END DO
        END IF

!       check if any part of pollutant within cloud layer
        CTOP=(ZMDL-ZSFC)*(1.0-ZSG(KTOP))
        CBOT=(ZMDL-ZSFC)*(1.0-ZSG(KBOT))

        IF(PBOT.LT.CTOP)THEN
!          fraction of pollutant below cloud top
           FRBCT=1.0
           IF(PTOP.GT.CTOP)FRBCT=1.0-(PTOP-CTOP)/PDEPTH

!          fraction of pollutant above cloud bottom
           FRACB=1.0
           IF(PBOT.LT.CBOT)FRACB=1.0-(CBOT-PBOT)/PDEPTH
           IF(PTOP.LT.CBOT)FRACB=0.0

           IF(DIRT(KT)%DOGAS)THEN
!             equilibrium concentration (mass units) for gases
!             applies as long as material is below cloud top
              IF(PBOT.LT.CTOP)THEN
!                deposition velocity
                 DEPV=DIRT(KT)%WETGAS*RGAS*TT(KLVL)*RAIN
!                rate constant
                 BETA=BETA+DEPV*FRBCT/PDEPTH
              END IF

           ELSE
!             for particles below cloud use scavenging coefficient
!             only fraction of mass below cloud removed
!             correction 10/27/98
              IF(PBOT.LT.CBOT)                                             &
                 BETA=BETA+DIRT(KT)%WETLO*(1.0-FRACB)*60.0

!             for particles within cloud use scavenging ratio
              IF(PTOP.GT.CBOT)THEN
!                deposition velocity
                 DEPV=DIRT(KT)%WETIN*RAIN
!                rate constant
                 BETA=BETA+DEPV*FRBCT*FRACB/PDEPTH
              END IF

!          gas or particle
           END IF

!       layer test
        END IF

!    rain test
     END IF

!    test for radioactive decay
     IF(DIRT(KT)%DORAD)THEN
!       convert half-life (days) to time constant (1/min)
        RTC=LOG(0.5)/(DIRT(KT)%RHALF*1440.0)
!       apply immediately to mass since it doesn't deposit
        MASS(KK)=MASS(KK)*EXP(ABS(DT)*RTC)
     END IF

     IF(BETA.EQ.0.0)THEN
!       no deposition
        DEPT(KK)=0.0
     ELSEIF(ABS(DT)*BETA.LT.0.01)THEN
!       small removal values assume linear approximation
        DEPT(KK)=MASS(KK)*ABS(DT)*BETA
        MASS(KK)=MASS(KK)-DEPT(KK)
     ELSE
!       apply exponential for removal > 1%
        DEPT(KK)=MASS(KK)*(1.0-EXP(-ABS(DT)*BETA))
!       can't remove more than exists
        DEPT(KK)=MIN(MASS(KK),DEPT(KK))
        MASS(KK)=MASS(KK)-DEPT(KK)
     END IF

! species loop
  END DO polnum

END SUBROUTINE depelm
