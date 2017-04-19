!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PARPUF           PARticle to PUFf conversion routine
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:04-07-23
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   IS USED TO CONVERT PUFFS TO PARTICLES WHEN THE PUFF SIZE GROWS TO
!   AN ARBITRARY RADIUS, USUALLY DEFINED AS EQUAL TO THE CONCENTRATION
!   GRID SIZE.  THE HORIZINTAL DISTRIBUTION VARIABLE IS USED TO 
!   IDENTIFY PUFFS THAT ARE TO BE CONVERTED.  VALUES OF HDWP GREATER
!   THAN 100 HAVE A SPECIAL MEANING:
!   --------------------------------------------------------
!   HDWP     Initial Distribution        Conversion Property
!   0        3D particle                 None
!   1        Gaussian Puff               None
!   2        Top-Hat Puff                None
!   3        Gaussian Puff/Particle      None
!   4        Top-Hat Puff/Particle       None
!   103      3D particle                 Gaussian Puff/Particle
!   104      3D particle                 Top-Hat Puff/Particle
!   130      Gaussian Puff/Particle      3D particle
!   140      Top-Hat Puff/Particle       3D particle
!   --------------------------------------------------------
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 23 Jul 2004 (RRD) - initial version
!                 15 Apr 2008 (RRD) - global model consistency
!                 01 Jul 2008 (RRD) - changed conversion test to age
!
! USAGE:  CALL PARPUF(KPM,PAGE,SIGU,SIGV,SIGW,HDWP,PGRD,CONAGE)
!
!   INPUT ARGUMENT LIST:
!   OUTPUT ARGUMENT LIST:
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE PARPUF(KPM,PAGE,SIGU,SIGV,SIGW,HDWP,PGRD,CONAGE)

  USE funits

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  INTEGER,  INTENT(IN)    :: kpm        ! total number of puffs or particles
  INTEGER,  INTENT(IN)    :: page (:)   ! particle age (minutes)        
  REAL,     INTENT(INOUT) :: sigu (:)   ! horizontal u-component turbulence
  REAL,     INTENT(INOUT) :: sigv (:)   ! horizontal v-component turbulence
  REAL,     INTENT(INOUT) :: sigw (:)   ! vertical   w-component turbulence
  INTEGER,  INTENT(INOUT) :: hdwp (:)   ! horizontal distribution of pollutant
  INTEGER,  INTENT(IN)    :: pgrd (:)   ! meteorological grid index 
  INTEGER,  INTENT(IN)    :: conage     ! conversion age (minutes)          

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  INTEGER :: kp,kn

!-------------------------------------------------------------------------------

  KN=0

! go through particles again to split if needed
  ploop : DO KP=1,KPM

!    check for on-grid and puff property 
     IF(HDWP(KP).NE.103.AND.HDWP(KP).NE.104) CYCLE ploop

!    counter number to be converted 
     KN=KN+1

!    determine if large enough for conversion
     IF(PAGE(KP).LT.CONAGE) CYCLE ploop

!    velocity variances all set to zero
     SIGU(KP)=0.0
     SIGV(KP)=0.0
     SIGW(KP)=0.0

!    particle distribution becomes puff 
     HDWP(KP)=MOD(HDWP(KP),100)

!    reduce number
     KN=KN-1

! particle loop
  END DO ploop

  IF(KN.GT.0)WRITE(KF21,*)' NOTICE parpuf: ',KN,' unconverted particles'

END SUBROUTINE parpuf
