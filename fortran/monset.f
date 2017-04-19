!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  MONSET           RETURNS CHARACTER MONTH FROM NUMERIC VALUE
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   RETURNS CHARACTER MONTH FROM NUMERIC VALUE
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 18 Feb 1997 (RRD) 
!                 12 Dec 2000 (RRD) - fortran90 upgrade
!                 02 Nov 2001 (RRD) - uppercase option
!
! USAGE:  CALL MONSET(MON,MONTH)
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

SUBROUTINE MONSET(MON,MONTH,KASE)

  IMPLICIT NONE 

  INTEGER,       INTENT(IN) :: mon         ! numeric month 1-12
  INTEGER,       INTENT(IN) :: kase        ! mixed (1) all upper (2)
  CHARACTER(3), INTENT(OUT) :: month       ! character month

  CHARACTER                 :: MONID(12,2)*3

  DATA MONID/'Jan','Feb','Mar','Apr','May','Jun',  &
             'Jul','Aug','Sep','Oct','Nov','Dec',  &
             'JAN','FEB','MAR','APR','MAY','JUN',  &
             'JUL','AUG','SEP','OCT','NOV','DEC'/  

  MONTH=MONID(MON,KASE)

END SUBROUTINE monset
