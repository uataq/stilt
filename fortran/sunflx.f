!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  SUNFLX           SUN FLuX incident solar radiation at sfc
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   SUN FLUX RETURNS THE INCIDENT SOLAR IRRADIATION AT THE SURFACE 
!   BASED UPON THE AVERAGE RH AND SOLAR ELEVATION ANGLE.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 02 Mar 1998 (RRD) - corrected 80% limit on Rh
!                 05 Sep 2000 (RRD) - fortran90 upgrade
!
! USAGE:  CALL SUNFLX(NLVL,SEA,QQ,SWF,TR)
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

SUBROUTINE SUNFLX(NLVL,SEA,QQ,SWF,TR)

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  INTEGER,  INTENT(IN)   :: nlvl      ! number of vertical levels in profile
  REAL,     INTENT(IN)   :: sea       ! sine of the solar elevation angle 
  REAL,     INTENT(IN)   :: qq(:)     ! RH fraction profile 
  REAL,     INTENT(OUT)  :: swf       ! incident short wave flux (w/m2)
  REAL,     INTENT(OUT)  :: tr        ! transmissivity

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  REAL, PARAMETER :: tcrit    = 0.5        ! critical transmissivity
  REAL, PARAMETER :: solc     = 1104.0     ! cloud top solar const (w/m2)

  REAL    :: fcld,ravg,rcnt,rhmax
  INTEGER :: k,kmax

!-------------------------------------------------------------------------------

! find max RH

  RHMAX=0.0
  DO K=1,NLVL
     IF(QQ(K).GT.RHMAX)THEN
        KMAX=K
        RHMAX=QQ(K)
     END IF
  END DO

! find 3 layer average about max

  RAVG=0.0
  RCNT=0.0
  DO K=MAX(1,KMAX-1),MIN(NLVL,KMAX+1)
     RAVG=RAVG+QQ(K)
     RCNT=RCNT+1.0
  END DO
  RAVG=MAX(0.80,RAVG/RCNT)

! fractional cloud cover when Rh >0.8

  FCLD=(5.0*(RAVG-0.8))**2
  FCLD=MAX(0.0,MIN(1.0,FCLD))

! transmissivity based upon critical RH (TR=0.5 when FCLD=1)

  TR=1.0-FCLD*(1.0-TCRIT)

! adjust solar constant for elevation angle
  SWF=MAX(0.0,SEA*SOLC*TR)

END SUBROUTINE sunflx
