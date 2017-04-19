!-------------------------------------------------------------------------------
! PRECIPITATION CHEMISTRY MODULE 
! REVISED: 12/01/94 
!          09/20/96 - corrected SO2 to SO4 conversion (molecular weights) 
!          08/27/03 - restructured to fit within Hysplit 4.6 
!          04/02/04 - generic file unit numbers

!-------------------------------------------------------------------------------
SUBROUTINE PRCHEM(TPPT,MASS,DEPT,QQ,ZDEP,GFLX,SSA,                             &
                  KA,JSTEP,LBOT,LTOP,NLVL,                                     &
                  TOTSO2,DEPVEL)

  USE funits

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! ARGUMENT DEFINITIONS
 
  REAL, INTENT(IN)    :: TPPT       ! PRECIPITATION RATE (M/HR)
  REAL, INTENT(INOUT) :: MASS(:) 
                                    ! MASS(1)- MASS OF SO2 (KG)
                                    ! MASS(2)- MASS OF SO4 DRY (KG)
                                    ! MASS(3)- MASS OF SO4 WET (KG)
  REAL, INTENT(OUT)   :: DEPT(:)
                                    ! DEPT(1)- DRY DEPOSIT SO2 (KG)
                                    ! DEPT(2)- DRY DEPOSIT SO4 (KG)
                                    ! DEPT(3)- DRY DEPOSIT SO4W IN PRECIP (KG)

  REAL, INTENT(IN)    :: QQ(:)      ! RELATIVE HUMIDITY FRACTION (0.70=70%) 
  REAL, INTENT(IN)    :: ZDEP       ! TOP OF SFC DEPOSITION LAYER (M) 
  REAL, INTENT(IN)    :: GFLX       ! SOLAR IRRADIATION (W/M2)
  REAL, INTENT(IN)    :: SSA        ! SIN OF THE SOLAR ANGLE

  INTEGER, INTENT(IN) :: KA         ! SIGMA LAYER INDEX NUMBER OF PUFF
  REAL,    INTENT(IN) :: JSTEP      ! TIME STEP (HOURS)
  REAL,    INTENT(IN) :: LBOT       ! BOTTOM SIGMA LAYER OF PUFF (M)
  REAL,    INTENT(IN) :: LTOP       ! TOP SIGMA LAYER OF PUFF (M)
  INTEGER, INTENT(IN) :: NLVL       ! NUMBER OF METEO LEVELS

  REAL, INTENT(IN)    :: TOTSO2     ! LOCAL AIR CONCENTRATION OF SO2 (KG/M^3)
  REAL, INTENT(IN)    :: DEPVEL(:)  ! DEPOSITION VELOCITY (M/S) BY SPECIES
                                    ! From main program deposition subroutine

!-------------------------------------------------------------------------------
! INTERNAL DEFINITIONS

  REAL :: SO4PPT  ! DEPOSITION OF SO4 IN PRECIPITATION (KG)
  REAL :: SO2DEP  ! DRY DEPOSITION OF SO2 (KG)
  REAL :: SO4DEP  ! DRY DEPOSITION OF SO4 (KG)

  LOGICAL :: DEFVAL=.TRUE. ! FLAG TO ENTER CONSTANTS ONE TIME
  LOGICAL :: DAYTM         ! FLAG TO INDICATE DAYTIME
  LOGICAL :: WETCNV        ! DO INCLOUD CONVERSION PROCESSES
  LOGICAL :: INCRMV        ! DO IN-CLOUD SCAVENGING BY RAIN
  LOGICAL :: EVAP          ! CONVERT WET TO DRY (.65) IF CLOUD EVAPORATED

  REAL :: KSO2DY     ! SO2 DAYTIME GAS OXIDATION RATE (/HR)
  REAL :: KSO2NT     ! SO2 NIGHTTIME GAS OXIDATION RATE (/HR)
  REAL :: KSO2W      ! SO2 WET OXIDATION RATE (/HR)
  REAL :: KSO2DD     ! DRY DEPOSITION VELOCITY OF SO2 DIVIDED BY DEPTH (/HR)
  REAL :: KSO4DD     ! DRY DEPOSITION VELOCITY OF SO4 DIVIDED BY DEPTH (/HR)
  REAL :: KSO2       ! SO2 DRY OXIDATION RATE (/HR)
  REAL :: SIV        ! AMOUNT OF SO2 DISSOLVED IN CLOUD WATER
  REAL :: LWC        ! LIQUID WATER CONTENT (G/M^3)
  REAL :: D2WFR      ! FRACTION OF SO4 ACTIVATED IN CLOUD
  REAL :: SCAV       ! FRACTION OF MASS REMOVED BY RAIN
  REAL :: H2O2       ! HYDROGEN PEROXIDE CONCENTRATION (PPB)
  REAL :: RHCRIT     ! MINIMUM RH FOR WET OXIDATION TO OCCUR (FRACTION)
  REAL :: CLDEVP=0.7 ! RELATIVE HUMIDITY WHEN CLOUD EVAPORATES

  REAL :: so2,so2i,so2dd,so2cnv,so2ddp,so2ppt
  REAL :: so4c,so4w,so4d,so4di,so4dd,so4ddp
  REAL :: betta,scavic,sumbal,rain2,depth,h2o2i 

  INTEGER :: kldtop,kldbot,kl 

!-------------------------------------------------------------------------------
! Set up constants 

  IF(DEFVAL)THEN
!    enter constants first time only
     DEFVAL=.FALSE.
     OPEN(KF33, FILE='PRCHEM.DAT')
     READ(KF33,*)RHCRIT
     READ(KF33,*)KSO2DY
     READ(KF33,*)KSO2NT
     READ(KF33,*)D2WFR
     READ(KF33,*)H2O2I
     READ(KF33,*)LWC
     CLOSE (KF33)
     WRITE(KF21,*)'*** Precipitation Chemistry Enabled ***'
     WRITE(KF21,*)'   --- Sulfate Chemistry Enabled ---'
  END IF

!-------------------------------------------------------------------------------
! Re-initialize precipitation depositions for new puff
  SO2DEP=0.0
  SO4DEP=0.0
  SO4PPT=0.0
  SO2PPT=0.0

! set puff depth variableI
  DEPTH=MAX(1.0,(LTOP-LBOT))

! check for day or night
  IF(GFLX.GT.0.0)THEN
     DAYTM=.TRUE.
  ELSE
     DAYTM=.FALSE.
  END IF

! set up rate constants for day or night
  IF(.NOT.DAYTM)THEN
     KSO2=KSO2NT
  ELSE
     KSO2=KSO2DY
  END IF

! Solar Sun Angle (SSA) from main program
! latitude variation of SO2 conversion
  KSO2 = KSO2DY * SSA
  IF(KSO2.LT.KSO2NT)KSO2=KSO2NT

! latitude variation of H2O2 concentration
  H2O2 = H2O2I * SSA
  IF(H2O2.LT.0.05)H2O2=0.05

! Compute SO2 K constant (dependent on Total SO2 air conc.)
  KSO2W=41.5*LWC*H2O2*EXP(-0.233E+09*TOTSO2)

! Convert deposition velocity (m/s) to m/hr & divide by layer depth
  KSO2DD=3600.*DEPVEL(1)/DEPTH
  KSO4DD=3600.*DEPVEL(2)/DEPTH

! make sure rain is >= 0.0 and multiply it by the time step in hours
  RAIN2=MAX(0.0,TPPT)*JSTEP

! Determine if particle is in cloud from relative humidity at
! sigma level of pollutant (q replaced by rh) and determine
! if in-cloud or below-cloud scavenging is to be done and
! whether or not a cloud has evaported.
  EVAP=.TRUE.
  WETCNV=.FALSE.
  INCRMV=.FALSE.

! determine top and bottom of cloud
  KLDBOT=0
  KLDTOP=NLVL
  DO KL=1,NLVL
     IF(KLDBOT.EQ.0.AND.QQ(KL).GE.RHCRIT) KLDBOT=KL
     IF(KLDBOT.NE.0.AND.KLDTOP.EQ.NLVL.AND.QQ(KL).LT.CLDEVP)KLDTOP=KL
  END DO     

! Do we have a cloud?
  IF(KLDBOT.GT.0)THEN
!    Is particle in cloud?
     IF(KA.GE.KLDBOT.AND.KA.LE.KLDTOP)THEN
        WETCNV=.TRUE.
        EVAP=.FALSE.
!       Is it raining?
        IF(RAIN2.GT.0.0)INCRMV=.TRUE.
     ELSE
!       Don't evaoporate cloud if RH > CLDEVP and < RHCRIT
        IF(QQ(KA).GE.CLDEVP)EVAP=.FALSE.
     END IF
  ELSE
!    Don't evaoporate cloud if RH > CLDEVP and < RHCRIT
     IF(QQ(KA).GE.CLDEVP)EVAP=.FALSE.
  END IF

! if it's raining don't evaporate
  IF(RAIN2.GT.0.0)EVAP=.FALSE.

! read in initial mass values (KG)
  SO2I= MASS(1)
  SO4DI=MASS(2)
  SO4W =MASS(3)

!*******************************************************************
!                    GAS-PHASE CHEMISTRY
!*******************************************************************

! no-cloud
  IF(.NOT.WETCNV)THEN

!    first sulfate dry deposition constants
     IF(LBOT.LE.ZDEP) THEN
        SO2DD=KSO2DD
        SO4DD=KSO4DD
     ELSE
        SO2DD=0.
        SO4DD=0.
     END IF

!--->fraction of SO2 remaining after conversion to dry SO4
     SO2CNV= EXP(-KSO2 * JSTEP)
!--->fraction of SO2 remaining after dry deposition
     SO2DDP= EXP(-SO2DD * JSTEP)
!--->total SO2 remaining after conversion and dry deposition
     SO2   = SO2I * SO2CNV * SO2DDP
!--->amount of SO2 dry deposited
     SO2DEP= SO2I*(1. - SO2DDP)
!--->total dry SO4 created from SO2
     SO4C  = SO2I*(1. - SO2CNV)
!--->account for molecular weights of SO4 and SO2 (96/64)
     SO4C = SO4C*1.5
!--->fraction of SO4 remaining after dry deposition
     SO4DDP= EXP(-SO4DD * JSTEP)
!--->total dry SO4 remaining after dry deposition
     SO4D  = (SO4DI + SO4C) * SO4DDP
!--->amount of SO4 dry deposited
     SO4DEP= (SO4DI + SO4C) * (1. - SO4DDP)
     SUMBAL = SO2I-(SO2+SO2DEP+SO4C)


  ELSE
!*******************************************************************
!                    IN-CLOUD CHEMISTRY
!*******************************************************************

!    If there is a cloud we do in-cloud chemistry.
     IF(LBOT.LE.ZDEP) THEN
        SO2DD=KSO2DD
        SO4DD=KSO4DD
     ELSE
        SO2DD=0.
        SO4DD=0.
     END IF

!--->fraction of SO2 remaining after conversion to wet SO4
     SO2CNV= EXP(-KSO2W*JSTEP)
!--->fraction of SO2 remaining after dry deposition
     SO2DDP= EXP(-SO2DD*JSTEP)
!--->total SO2 remaining after conversion and dry deposition
     SO2   = SO2I*SO2CNV*SO2DDP
!--->amount of SO2 dry deposited
     SO2DEP= SO2I*(1.-SO2DDP)
!--->total wet SO4 created from SO2
     SO4C  = SO2I*(1.-SO2CNV)
!--->account for molecular weights of SO4 and SO2 (96/64)
     SO4C = SO4C*1.5
!--->fraction of SO4 remaining after dry deposition
     SO4DDP= EXP(-SO4DD*JSTEP)
!--->total SO4 dry after deposition
     SO4D = SO4DI*SO4DDP
!--->amount of SO4 dry deposited
     SO4DEP = SO4DI*(1.-SO4DDP)
!--->total wet SO4 remaining including amount converted from dry SO4
     SO4W  = SO4W + SO4C + D2WFR*SO4D
!--->Amount of dry SO4 remaining
     SO4D   =(1.0 - D2WFR)*SO4D
     SUMBAL= SO2I-(SO2+SO2DEP+SO4C)
  END IF

!*******************************************************************
!                    EVAPORATION OF CLOUD
!*******************************************************************

! when RH is less than CLDEVP we consider any cloud to have evaporated,
! and therefore some or all of the wet component goes back to dry.

! in dry air all wet SO4 evaporates to dry SO4
  IF(EVAP)THEN
     SO4D=SO4D+SO4W
     SO4W=0.0
  END IF

!*******************************************************************
!                    WET DEPOSITION
!*******************************************************************

! wet deposition in kilograms from puff from in-cloud

  IF(RAIN2.GT.0.0) THEN
     IF(INCRMV)THEN
!------>Compute repetative function
        BETTA=RAIN2/DEPTH

!------>additional so2 removed by dissolving in cloud water
!------>In-cloud scavenging ratio
        SCAVIC = 2.5E+05
        SIV= 0.053*LWC*(SO2/(1.0+0.053*LWC))
        SO2PPT = SIV * SCAVIC * BETTA
        SO2PPT = MIN( SIV, SO2PPT)
        SO2 =   SO2 - SO2PPT

!------>now SO4
!------>In-cloud scavenging ratio
        SCAVIC = 2.5E+05
        SO4PPT = SO4W * SCAVIC * BETTA
        SO4PPT = MIN( SO4W, SO4PPT)
        SO4W =   SO4W - SO4PPT

     ELSEIF(.NOT.EVAP)THEN
!------>wet deposition in kilograms if not in cloud
!------>Below-cloud scavenging coefficient in units of /s

!------>SO4
        SCAVIC=2.5E-05
        SCAV =   EXP(-3600. * JSTEP * SCAVIC)
        SCAV =   MIN( SCAV , 1.0)
        SO4PPT = (SO4D + SO4W) * (1.0 - SCAV)
        SO4D   = SO4D * SCAV
        SO4W   = SO4W * SCAV
     END IF
  END IF

!*******************************************************************
!                    WRITE FINAL VALUES
!*******************************************************************

  MASS(1)=SO2
  MASS(2)=SO4D
  MASS(3)=SO4W

  DEPT(1)=SO2DEP+SO2PPT
  DEPT(2)=SO4DEP
  DEPT(3)=SO4PPT

END SUBROUTINE prchem
