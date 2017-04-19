!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  LAGOUT           LAGrangian sampler file OUTput
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:05-10-14
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   WRITES ALL THE CONCENTRATIONS COLLECTED BY THE LAGRANGIAN    
!   SAMPLERS AS ONE OUTPUT FILE PER SAMPLER
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 18 Oct 2005 (RRD) - Initial version from parout
!                 18 Jan 2008 (RRD) - absolute value for page
!
! USAGE:  CALL LAGOUT(LAGS,JET,KPM,DELT,MASS,XPOS,YPOS,ZPOS,HDWP,PAGE,
!                     PTYP,PGRD,ZMDL)
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

SUBROUTINE LAGOUT(LAGS,JET,KPM,DELT,MASS,XPOS,YPOS,ZPOS,HDWP,PAGE,    &
                  PTYP,PGRD,ZMDL)

  USE funits
  use module_defgrid                   ! meteorology grid information

  IMPLICIT NONE

  INCLUDE 'DEFLAGS.INC'                   ! lagrangian sampling configuration

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  TYPE(lset), INTENT(INOUT) :: lags(:)    ! lagrangian sampler confuration
  INTEGER,    INTENT(IN)    :: jet        ! current elapsed time
  INTEGER,    INTENT(IN)    :: kpm        ! particle count
  REAL,       INTENT(IN)    :: delt       ! current integration time step
  REAL,       INTENT(IN)    :: mass (:,:) ! mass of pollutant (arbitrary units)
  REAL,       INTENT(IN)    :: xpos (:)   ! center positions (grid units)
  REAL,       INTENT(IN)    :: ypos (:)   ! center positions (grid units)
  REAL,       INTENT(IN)    :: zpos (:)   ! center height (sigma)
  INTEGER,    INTENT(IN)    :: hdwp (:)   ! always=6 as lagrangian marker    
  INTEGER,    INTENT(IN)    :: page (:)   ! averaging counter 
  INTEGER,    INTENT(IN)    :: ptyp (:)   ! lagrangian sampler number  
  INTEGER,    INTENT(IN)    :: pgrd (:)   ! meteorological grid of puff position
  REAL,       INTENT(IN)    :: zmdl       ! vertical index scaling height

!-------------------------------------------------------------------------------

  REAL     :: AVRG                        ! average concentration
  REAL     :: TLAT,TLON,THGT              ! particle position
  INTEGER  :: IYR,IMO,IDA,IHR,IMN         ! date/time
  INTEGER  :: J,KUNIT                     ! index and file unit number  


!-------------------------------------------------------------------------------
! find each lagrangian sampler
!-------------------------------------------------------------------------------

  CALL TM2DAY(JET,IYR,IMO,IDA,IHR,IMN)

  jloop : DO J=1,KPM    

!    test for lagrangian sampler
     IF(HDWP(J).NE.6) CYCLE jloop

!    has the sampler been shutdown (perhaps exceeded max number permitted)
     IF(LAGS(PTYP(J))%SACM.EQ.0) CYCLE jloop      

!    is it less that the output time for this sampler number
     IF(JET.LT.LAGS(PTYP(J))%SACM) CYCLE jloop      

!    update the time for the next output
     IF(LAGS(PTYP(J))%DISK.EQ.0)THEN
        LAGS(PTYP(J))%SACM=LAGS(PTYP(J))%SACM+INT(DELT)
     ELSE
        LAGS(PTYP(J))%SACM=LAGS(PTYP(J))%SACM+LAGS(PTYP(J))%DISK
     END IF

!    valid output units range from KF51 to KF52 (see funits.f)
     KUNIT=KF51+PTYP(J)-1

!    convert positions to lat/lon before dump
     IF(GRID(PGRD(J),1)%LATLON)THEN
        CALL GBL2LL(PGRD(J),1,XPOS(J),YPOS(J),TLAT,TLON)
     ELSE
        CALL CXY2LL(GRID(PGRD(J),1)%GBASE,XPOS(J),YPOS(J),TLAT,TLON)
     END IF

!    convert sigma to agl assuming terrain height = 0
     THGT=(1.0-ZPOS(J))*ZMDL

     AVRG=0.0
     IF(PAGE(J).NE.0) AVRG=MASS(1,J)/ABS(PAGE(J))

     WRITE(KUNIT,'(5I3,2F10.3,F10.1,E15.4 )')    &
           IYR,IMO,IDA,IHR,IMN,TLAT,TLON,THGT,AVRG

  END DO jloop

END SUBROUTINE lagout
