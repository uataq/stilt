!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  EMRISE           EMission plume RISE from Briggs 
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:06-03-08
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   EMISSION PLUME RISE COMPUTES AN AIR PARCEL'S RISE THROUGH
!   BOUANCY ONLY BASED UPON THE HEAT RELEASED AT THE SURFACE.
!   EQUATIONS ARE BASED UPON BRIGGS (1969) WITH UPDATES BY 
!   SP ARYA (1999).
!
! PROGRAM HISTORY LOG:
!   Last Revision: 08 Mar 2006 (RRD) - initial version
!                  19 Mar 2007 (AS)  - added duration test
!
! USAGE:  CALL EMRISE (KDUR,HEAT,STAB,UBAR,USTR,MIXD,RISE)
!
!   INPUT ARGUMENT LIST:      see below
!   OUTPUT ARGUMENT LIST:     see below
!   INPUT FILES:              none
!   OUTPUT FILES:             none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE EMRISE(KDUR,HEAT,STAB,UBAR,USTR,MIXD,RISE)

  USE funits

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: KDUR   ! emission duration (min)
  REAL(4), INTENT(IN)  :: HEAT   ! heat release (watts)
  REAL(4), INTENT(IN)  :: STAB   ! static stability (1/s2)
  REAL(4), INTENT(IN)  :: UBAR   ! wind speed (m/s)
  REAL(4), INTENT(IN)  :: USTR   ! friction velocity (m/s)
  REAL(4), INTENT(IN)  :: MIXD   ! mixed layer depth (m)  
  REAL(4), INTENT(OUT) :: RISE   ! plume rise (m)

  REAL(4)              :: FB     ! bouyancy flux (m4/s3)

! bouancy flux based upon heat release
  FB=7.6E-07*HEAT

  IF(STAB.LE.0.0)THEN
!    neutral unstable conditions
     RISE=1.3*FB/UBAR/USTR/USTR

!    check for limits
     IF(KDUR.LT.1440)THEN
!       short-duration fires cannot penetrate PBL
        RISE=MIN(0.75*MIXD, RISE)
     ELSE
!       three kilometer limit typical daytime rise
        RISE=MIN(3000.0,RISE)
     END IF

  ELSE
     IF(UBAR.GT.0.5)THEN
!       stable windy
        RISE=2.6*(FB/UBAR/STAB)**0.333333
     ELSE
!       stable calm    
        RISE=5.3*(FB**0.25)/(STAB**0.375)
     END IF

!    check for limits
     IF(KDUR.LT.1440)THEN
        RISE=MIN(2.0*MIXD, RISE)
     ELSE   
        RISE=MIN(3000.0,RISE)
     END IF 
  END IF

! WRITE(KF21,*)' NOTICE: emrise (mixd,rise) - ',MIXD,RISE,UBAR,USTR,STAB

END SUBROUTINE emrise
