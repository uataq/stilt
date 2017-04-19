!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  LAGSUM           LAGrangian sampler SUMmation  
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:05-10-14
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   WRITES ALL THE CONCENTRATIONS COLLECTED BY THE LAGRANGIAN    
!   SAMPLERS AS ONE OUTPUT FILE PER SAMPLER
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 17 Oct 2005 (RRD) - Initial version from parout
!                 17 Oct 2007 (RRD) - revised grid start/stop test 
!                 18 Jan 2008 (RRD) - absolute value page
!
! USAGE:  CALL LAGSUM(LAGS,CONC,NUMGRD,PLAT,PLON,ZPOS,DELT,JET,ZMDL,ZSFC,  &
!                     KMSL,MASS,PAGE,PTYP,CSUM)
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

SUBROUTINE LAGSUM(LAGS,CONC,NUMGRD,PLAT,PLON,ZPOS,DELT,JET,ZMDL,ZSFC,  &
                  KMSL,MASS,PAGE,PTYP,CSUM)

  USE funits

  IMPLICIT NONE

  INCLUDE 'DEFLAGS.INC'                   ! lagrangian sampling configuration
  INCLUDE 'DEFCONC.INC'                   ! pollutant and concentration grid

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  TYPE(lset), INTENT(INOUT) :: lags(:)    ! lagrangian sampler configuration
  TYPE(cset), INTENT(IN)    :: conc(:)    ! for each concentration grid 
  INTEGER,    INTENT(IN)    :: numgrd     ! number of concentration grids
  REAL,       INTENT(IN)    :: plat,plon  ! particle position lat & long
  REAL,       INTENT(IN)    :: zpos       ! center height (sigma)
  REAL,       INTENT(IN)    :: delt       ! current integration time step
  INTEGER,    INTENT(IN)    :: jet        ! current elapsed time
  REAL,       INTENT(IN)    :: zmdl       ! model domain top
  REAL,       INTENT(IN)    :: zsfc       ! height of ground surface (m)
  INTEGER,    INTENT(IN)    :: kmsl       ! agl=0 or msl=1 height units
  REAL,       INTENT(INOUT) :: mass (:)   ! mass of pollutant (arbitrary units)
  INTEGER,    INTENT(INOUT) :: page       ! averaging counter 
  INTEGER,    INTENT(IN)    :: ptyp       ! lagrangian sampler number  
  REAL,       INTENT(IN)    :: csum (:,:,:,:,:)  ! concentration summation 

!-------------------------------------------------------------------------------

  INTEGER  :: KSB,KSD                     ! time in sample and duration
  REAL     :: ZPAR,ZBOT,ZTOP              ! vertical positions in meters
  REAL     :: GLAT,GLON,DELY              ! final corrected particle position
  INTEGER  :: IYR,IMO,IDA,IHR,IMN         ! date/time
  INTEGER  :: KUNIT                       ! file unit number  
  INTEGER  :: KG                          ! snapshot concentration grid
  INTEGER  :: II,JJ,KK                    ! position on the concentration grid
  LOGICAL  :: DIAG = .TRUE.               ! diagnostic message only once

  INTEGER  :: j,km                        ! indicies

  SAVE DIAG

  KM=1  ! current version supports only mass index=1 

!-------------------------------------------------------------------------------
! zero the summation if required (avrg=0 for snapshot output)
!-------------------------------------------------------------------------------

  IF(ABS(PAGE).GE.LAGS(PTYP)%AVRG.OR.LAGS(PTYP)%AVRG.EQ.0) THEN
     MASS(KM)=0.0
     PAGE=0
  END IF 

!-------------------------------------------------------------------------------
! find the snapshot concentration grid
!-------------------------------------------------------------------------------

  KG=0
  DO J=1,NUMGRD
     IF(CONC(J)%SNAP.EQ.1) KG=J
  END DO

  IF(KG.EQ.0)THEN
     WRITE(KF21,*)'*ERROR* lagsum: snapshot concentration grid not found'
     STOP 900
  END IF

! sample start stop 
  IF(INT(DELT).GT.0)THEN
     KSB=JET-CONC(KG)%START%MACC
  ELSE
     KSB=CONC(KG)%START%MACC-JET
  END IF
  KSD=ABS(CONC(KG)%STOP%MACC-CONC(KG)%START%MACC)

! test for time within sampling interval
  IF(KSB.LT.0.OR.KSB.GE.KSD) THEN 
     IF(DIAG)THEN
        WRITE(KF21,*)'WARNING lagsum: snapshot concentration grid turned off'
        WRITE(KF21,*)'... outside of temporal start-stop limits'
        DIAG=.FALSE.
     END IF
     RETURN
  END IF

!-------------------------------------------------------------------------------
! find the location of the lagrangian sampler on the concentration grid
!-------------------------------------------------------------------------------

! compute longitude domain to test for dateline span
  DELY=CONC(KG)%DELT_LON*CONC(KG)%NUMB_LON

! longitude correction to avoid dateline problems 
  IF(CONC(KG)%X1Y1_LON.GE.0.0.AND.PLON.LT.0.0)THEN
     GLON=PLON+360.0
! large grids may span dateline 
  ELSEIF((DELY.GT.180.0).AND.(PLON-CONC(KG)%X1Y1_LON.LE.0.0))THEN
     GLON=PLON+360.0
  ELSE
     GLON=PLON
  END IF

! find the particle position on current concentration grid
  II=NINT(1.0+(GLON-CONC(KG)%X1Y1_LON)/CONC(KG)%DELT_LON)
  JJ=NINT(1.0+(PLAT-CONC(KG)%X1Y1_LAT)/CONC(KG)%DELT_LAT)
  IF(JJ.GT.CONC(KG)%NUMB_LAT.OR.JJ.LT.1.OR.             &
     II.GT.CONC(KG)%NUMB_LON.OR.II.LT.1) THEN   
     IF(DIAG)THEN
        WRITE(KF21,*)'WARNING lagsum: lagrangian sampler moved off the grid'  
        DIAG=.FALSE.
     END IF
     RETURN
  END IF

! vertical position in meters
  ZPAR=MAX(1.0,(ZMDL-ZSFC)*(1.0-ZPOS)) 
  IF(KMSL.EQ.1)ZPAR=ZPAR+ZSFC             ! vertical units MSL rather than AGL

!-------------------------------------------------------------------------------
! find the vertical position of the sampler on the concentration grid
! currently only one mass species is valid for sampling (first index)
!-------------------------------------------------------------------------------

  ZBOT=0.0
  DO KK=1,CONC(KG)%LEVELS
 
!    determine vertical cell sizes (input defines top)
     ZTOP=FLOAT(CONC(KG)%HEIGHT(KK))

!    particle falls within vertical extent of cell
     IF(ZPAR.GT.ZBOT.AND.ZPAR.LE.ZTOP)THEN
        MASS(KM)=MASS(KM)+CSUM(II,JJ,KK,KM,KG)*ABS(DELT)
        PAGE=SIGN( (ABS(PAGE)+INT(ABS(DELT))), PAGE)
        RETURN 
     END IF
     ZBOT=ZTOP
  END DO 

  IF(DIAG)THEN
     WRITE(KF21,*)'WARNING lagsum: lagrangian sampler above conc grid'
     DIAG=.FALSE.
  END IF

END SUBROUTINE lagsum
