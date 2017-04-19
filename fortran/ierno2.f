!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  IERNO2           Computes NO2 photolysis rate coefficient
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!            uses the solar angle and surface short wave flux to compute
!            the NO2 photolysis rate coefficient.  The rate coefficient
!            is used in the integration of VOC for smog produced
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 15 Jan 1997 (RRD)
!                 24 Jul 2003 (RRD) - fortran90 upgrade
!
! USAGE:  CALL IERNO2(EA,SWF,RK1)
!   INPUT ARGUMENT LIST: 	see below
!   OUTPUT ARGUMENT LIST:	see below
!   INPUT FILES:		none
!   OUTPUT FILES:		none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE IERNO2(EA,SWF,RK1)

  IMPLICIT NONE

  REAL, INTENT(IN)  :: EA      ! solar elevation angle (deg)
  REAL, INTENT(IN)  :: SWF     ! short-wave flux (w/m2)
  REAL, INTENT(OUT) :: RK1     ! rate coefficient for NO2 photolysis

  REAL, PARAMETER   :: RPDG = 0.01745329  ! radians per degree

  REAL  :: zag,zar

! zenith angle in radians from elevation angle

  ZAG=90.0-EA
  ZAR=ZAG*RPDG

! rate coefficient for NO2 photolysis

  IF(ZAG.GE.0.0.AND.ZAG.LT.47.0)THEN
     RK1=(4.23E-04+(1.09E-04/COS(ZAR)))*SWF
  ELSEIF(ZAG.GE.47.0.AND.ZAG.LT.64.0)THEN
     RK1=5.82E-04*SWF
  ELSEIF(ZAG.GE.64.0.AND.ZAG.LE.90.0)THEN
     RK1=(-0.997E-04+1.2E-03*(1.0-COS(ZAR)))*SWF
  ELSE
     RK1=0.0
  END IF

END SUBROUTINE ierno2
