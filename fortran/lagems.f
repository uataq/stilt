!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  LAGEMS           EMiSsion of LAGrangian samplers
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:05-10-13
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   EMISSION INITIALIZATION FOR THE SPECIAL CASE OF LAGRANGIAN  
!   SAMPLERS (E.G. A BALLOON). MULTIPLE LOCATIONS, EACH WITH
!   ITS OWN RELEASE TIME MAY BE DEFINED IN FILE LAGSET.CFG.
!   EACH SAMPLER IS RELEASED ONLY ONCE.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 17 Oct 2005 (RRD) - initial version
!
! USAGE:  CALL LAGEMS(NUMLAG,LAGS,MAXPAR,JET,KPM,NSORT,MASS,XPOS,YPOS,
!                     ZPOS,SIGH,SIGU,SIGV,SIGW,HDWP,PAGE,PTYP,PGRD)
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

SUBROUTINE LAGEMS(NUMLAG,LAGS,MAXPAR,JET,KPM,NSORT,MASS,XPOS,YPOS,ZPOS, &
                  SIGH,SIGU,SIGV,SIGW,HDWP,PAGE,PTYP,PGRD)

  USE funits

  IMPLICIT NONE

  INCLUDE 'DEFLAGS.INC'                    ! meteorology grid and file

  INTEGER,    INTENT(IN)    :: numlag      ! total number of lagrangian samp
  TYPE(lset), INTENT(INOUT) :: lags(:)     ! lagrangian sampler confuration
  INTEGER,    INTENT(IN)    :: maxpar      ! maximum particle number
  INTEGER,    INTENT(IN)    :: jet         ! current elapsed time (min)
  INTEGER,    INTENT(INOUT) :: kpm         ! number of puffs or particles

  INTEGER,    INTENT(INOUT) :: nsort (:)   ! index of sorted elements
  REAL,       INTENT(INOUT) :: mass  (:,:) ! mass of pollut (arbitrary units)
  REAL,       INTENT(INOUT) :: xpos  (:)   ! horizontal position (grid units)
  REAL,       INTENT(INOUT) :: ypos  (:)   ! horizontal position (grid units)
  REAL,       INTENT(INOUT) :: zpos  (:)   ! puff center height (sigma)
  REAL,       INTENT(INOUT) :: sigh  (:)   ! horizontal puff sigma 
  REAL,       INTENT(INOUT) :: sigu  (:)   ! turbulence u'2    
  REAL,       INTENT(INOUT) :: sigv  (:)   ! turbulence v'2
  REAL,       INTENT(INOUT) :: sigw  (:)   ! turbulence w'2 or vertical  sigma
  INTEGER,    INTENT(INOUT) :: hdwp  (:)   ! Horizontal distribution pollutant
  INTEGER,    INTENT(INOUT) :: page  (:)   ! lagrangian samp averaging counter
  INTEGER,    INTENT(INOUT) :: ptyp  (:)   ! lagrangian samp sequence number
  INTEGER,    INTENT(INOUT) :: pgrd  (:)   ! meteorological grid puff position

  INTEGER                   :: n

!-------------------------------------------------------------------------------
! loop through each lagrangian sampler
!-------------------------------------------------------------------------------

nloop : DO N=1,NUMLAG

    IF(JET.LT.LAGS(N)%RACM.OR.LAGS(N)%RACM.EQ.0) CYCLE nloop      

    LAGS(N)%RACM=0 ! only one emission permitted, then shutdown

    KPM=KPM+1
    IF(KPM.GT.MAXPAR)THEN
       KPM=MAXPAR
       WRITE(KF21,*)'WARNING lagems: exceeding source array limit'
       WRITE(KF21,*)'... increase the value of MAXPAR in SETUP.CFG' 
       EXIT nloop
    END IF
    NSORT(KPM)=KPM

!   initial position 
    XPOS(KPM)=LAGS(N)%XP
    YPOS(KPM)=LAGS(N)%YP
    ZPOS(KPM)=LAGS(N)%ZP

!   variances all start at zero
    SIGH(KPM)=0.0
    SIGU(KPM)=0.0
    SIGV(KPM)=0.0
    SIGW(KPM)=0.0

!   initial distribution (see main for definitions)
    HDWP(KPM)=6      
!   initial age at zero is the averaging counter
    PAGE(KPM)=0
!   pollutant type is the sampler number
    PTYP(KPM)=N
!   initial meteorology grid 
    PGRD(KPM)=LAGS(N)%KG

!   mass summation array
    MASS(:,KPM)=0.0

END DO nloop

END SUBROUTINE lagems
