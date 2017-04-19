!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  TM2JUL           TiMe2JUL converts time to day since 1 Jan
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL      DATE:00-07-06
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   TM2JUL CONVERTS THE CURRENT TIME TO SEQUENTIAL DAY AS COUNTED FROM
!   THE BEGINNING OF EACH YEAR.
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 06 Jul 2000 (RRD)
!                  02 Sep 2000 (RRD) - fortran90 upgrade
!                  16 Mar 2001 (RRD) - argument list change
!
! USAGE:  CALL TM2JUL(IY,IM,ID,JULD)
!
!   INPUT ARGUMENT LIST:     see below
!   OUTPUT ARGUMENT LIST:    see below
!   INPUT FILES:             none
!   OUTPUT FILES:            none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE TM2JUL(IY,IM,ID,JULD)

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  INTEGER,   INTENT(IN)  :: iy,im,id         ! date  
  INTEGER,   INTENT(OUT) :: juld             ! days since Jan 1 st

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

!             default number of accumulated days in each month (non leap year)
  INTEGER  :: nadpm(12) = (/0,31,59,90,120,151,181,212,243,273,304,334/)
  INTEGER  :: iyear

!-------------------------------------------------------------------------------

!     compute four digit year
      IYEAR=MOD(IY+99,100)+1901

!     add accumulated days for this year
      JULD=NADPM(IM)+ID

!     adjust accumulated days for leap year
!     does not account for centuries where mod(iyear,400)=0
      IF(MOD(IYEAR,4).EQ.0.AND.IM.GT.2)JULD=JULD+1

END SUBROUTINE tm2jul
