!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  DEPDRY           DEPosition DRY via the resistance method.
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   DEPOSITION DRY VIA THE RESISTANCE METHOD. REQUIRES INPUT FILES
!   FOR ROUGHNESS LENGTH AND LANDUSE.  PREVIOUS ROUTINES USED TO
!   CALCULATE SOLAR INSOLATION FROM TIME, POSITION, AND HUMIDITY.
!   INPUT DATA INITIALIZED IN SUBROUTINE SFCINP.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 03 Mar 1997 (RRD)
!                 20 Apr 1999 (RRD) - variable name change
!                 03 Sep 2000 (RRD) - fortran90 upgrade
!                 09 Sep 2002 (RRD) - fortran coding standards
!
! USAGE:  CALL DEPDRY(DIRT,OLAT,IBMO,KPOL,LAND,ROUG,SFCL,USTR,PSI,SFLX,AIRD,
!                     TEMP,PDIA,VG,VD)
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

SUBROUTINE DEPDRY(DIRT,OLAT,IBMO,KPOL,LAND,ROUG,SFCL,USTR,PSI,SFLX,AIRD,   &
                  TEMP,PDIA,VG,VD)

  IMPLICIT NONE

  INCLUDE 'DEFCONC.INC'       ! concentration and pollutant structure

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

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

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  REAL, PARAMETER :: VONK  = 0.4         ! VonKarman 
  REAL, PARAMETER :: GRAV  = 9.801       ! gravity (m/s2)
  REAL, PARAMETER :: DMVC  = 1.789E-02   ! dynamic viscosity (g m-1 s-1)
  REAL, PARAMETER :: DSTP  = 1.2E+03     ! stp density (g/m3)
  REAL, PARAMETER :: FREP  = 6.53E-08    ! mean free path at stp (m)
  REAL, PARAMETER :: AMWT  = 28.9644     ! air molec wt (g/mole)
  REAL, PARAMETER :: BOLT  = 1.3804E-20  ! boltzmann (g-m2 /K-s2)
  REAL, PARAMETER :: PR    = 0.923       ! Prandtl

! scaling resistances (sec/meter), 11 categories, summer and winter
! (urban, agricul, dry range, deciduous, coniferous, mixed forest
! water, desert, wetlands, mixed range, rocky)

  REAL, DIMENSION (11,2) :: ri,rlu,rac,rg1,rg2,rc1,rc2 

  INTEGER :: klim
  REAL    :: rl,rd,ru,rm,rs,rz,rg,rx,ra,rb,rc
  REAL    :: cimp,tfact,c3,c2,c4,c1,tadj,stke,dhx,f0,z0
  REAL    :: gama,tsfc,hnry,pmwt,zlog,shmt,sc,diff,frea

!-------------------------------------------------------------------------------

  DATA RI     /1E25, 60.,120., 70.,130.,100.,1E25,1E25, 80.,100.,150.,         &
               1E25,1E25,1E25,1E25,400.,800.,1E25,1E25,1E25,1E25,1E25/
  DATA RLU    /1E25, 2E3, 2E3, 2E3, 2E3, 2E3,1E25,1E25,25E2, 2E3, 4E3,         &
               1E25,1E25,1E25,1E25, 6E3, 9E3,1E25,1E25, 9E3, 9E3, 9E3/
  DATA RAC    /100.,200.,100., 2E3, 2E3, 2E3,  0.,  0.,300.,150.,200.,         &
               100., 10., 10., 1E3, 2E3,15E2,  0.,  0., 50., 10., 50./
  DATA RG1    /400.,150.,350.,500.,500.,100.,  0., 1E3,  0.,220.,400.,         &
               100.,100.,100.,100.,100.,100.,  0., 1E3,100.,100., 50./
  DATA RG2    /300.,150.,200.,200.,200.,300., 2E3,400., 1E3,180.,200.,         &
               600.,35E2,35E2,35E2,35E2,35E2, 2E3,400.,35E2,35E2,35E2/
  DATA RC1    /1E25, 2E3, 2E3, 2E3, 2E3, 2E3,1E25,1E25,25E2, 2E3, 4E3,         &
               1E25,1E25,1E25, 9E3,200.,400.,1E25,1E25, 9E3,1E25, 9E3/
  DATA RC2    /1E25, 1E3, 1E3, 1E3, 1E3, 1E3,1E25,1E25, 1E3, 1E3, 1E3,         &
               1E25, 1E3, 1E3,400.,15E2,600.,1E25,1E25,800., 1E3,800./

  SAVE ri,rlu,rac,rg1,rg2,rc1,rc2 

!-------------------------------------------------------------------------------
! unit conversions
!-------------------------------------------------------------------------------

! canopy surface temperature in deg C
  TSFC=TEMP-273.16
! pollutant molecular weight
  PMWT=DIRT(KPOL)%GPMOL
! henry's constant for gases
  HNRY=DIRT(KPOL)%HENRY
! activity parameter
  F0=DIRT(KPOL)%ACVTY
! diffusivity ratio of pollutant to that of H2O
  DHX=DIRT(KPOL)%DIFTY
! kinematic viscosity
  GAMA=DMVC/AIRD

!-------------------------------------------------------------------------------
! atmospheric resistance
!-------------------------------------------------------------------------------

  Z0=ROUG
! roughness length over land then over water ...
! Smith (1988) correction for Charnock (1955)
  IF(LAND.EQ.7)Z0=(0.011*USTR*USTR/GRAV)+GAMA/(9.1*USTR)

! stability correction
  ZLOG=LOG(SFCL/Z0)
  RA=PR*(ZLOG-PSI)/(VONK*USTR)

!-------------------------------------------------------------------------------
! quasilaminar sublayer resistance
!-------------------------------------------------------------------------------

  IF(PDIA.EQ.0.0)THEN
!    diffusivity of gas pollutant from grahams law
     DIFF=GAMA*SQRT(AMWT/PMWT)
  ELSE
!    particles slip correction (mean free path = particle diameter)
     FREA=FREP*(DSTP/AIRD)
     SC=1.0+(2.0*FREA/PDIA)*(1.26+0.4*EXP(-0.55*PDIA/FREA))
     DIFF=BOLT*TEMP*SC/(9.42*DMVC*PDIA)
  END IF

! schmidt number
  SHMT=GAMA/DIFF

! brownian diffusion resistance component for gasses
  RB=2.0*PR*SHMT**0.67/(VONK*USTR)

! stokes number based correction for particle impaction
! most appropriate for homogeneous uniform surfaces
  IF(PDIA.GT.0.0)THEN
!    stokes number = inertial stopping / length scale
!    where the length scale is defined by laminar layer==>
     STKE=2.0*VG*USTR/(GRAV*SC*GAMA/USTR)
     CIMP=USTR*(STKE/(STKE+0.8))**2
     RB=PR/(1.0/RB+CIMP)
  END IF

! over-water correction
  IF(LAND.EQ.7)RB=10.0

!-------------------------------------------------------------------------------
! canopy resistance for gases only
!-------------------------------------------------------------------------------

  IF(PDIA.LE.0.0)THEN

     KLIM=1
!    determine season according to hemisphere (maxveg=1)
!    default max vegetation - should be customized for application
     IF(OLAT.GT.0.0)THEN
        IF(IBMO.GE.10.OR.IBMO.LE.3)KLIM=2
     ELSE
        IF(IBMO.GE.4.AND.IBMO.LE.9)KLIM=2
     END IF

!    stomatal resistance
     RS=1.0E+25
     IF(TSFC.GT.0.0.AND.TSFC.LT.40.0)THEN
        TFACT=400.0/TSFC/(40.0-TSFC)
        RS=RI(LAND,KLIM)*DHX*(1.+(200./(SFLX+0.1))**2)*TFACT
     END IF

!    mesophyll
     RM=1.0/(HNRY/3000.0+100.0*F0)
!    upper canopy leaf cuticles
     RU=RLU(LAND,KLIM)/(HNRY*1.0E-05+F0)
!    gas-phase convection
     RD=100.0*(1.0+1000.0/(SFLX+10.0))
!    lower canopy surfaces
     RL=1.0/(HNRY*1.0E-05/RC1(LAND,KLIM)+F0/RC2(LAND,KLIM))
!    canopy height factor
     RZ=MAX( 1.0,RAC(LAND,KLIM) )
!    ground resistance
     RX=MAX( 1.0,RG1(LAND,KLIM) )
     RG=1.0/(HNRY*1.0E-05/RX+F0/RG2(LAND,KLIM))

!    temperature adjustment
     TADJ=1000.0*EXP(-TSFC-4.0)
     C1=1.0/(RS+RM)
     C2=1.0/(RU+TADJ)
     C3=1.0/(RD+RL+TADJ)
     C4=1.0/(RZ+RG+TADJ)

!    total bulk canopy resistance
     RC=1.0/(C1+C2+C3+C4)
!    limit unrealistic high and low values
     RC=MAX(10.0,RC)
     RC=MIN(9999.0,RC)

  ELSE
!    no resistance to particle deposition
     RC=0.0
  END IF

! FINAL DRY DEPOSITION
  VD=1.0/(RA+RB+RC+RA*RB*VG)+VG

END SUBROUTINE depdry
