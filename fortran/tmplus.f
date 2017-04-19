!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  TMPLUS           TiMePLUS is used to add time
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   TIMEPLUS CONVERTS A DATE PLUS HOURS TO THE NEW DATE
!
! PROGRAM HISTORY LOG:
!   Last Revised: 14 Feb 1997 (RRD)
!                 06 Aug 1999 (RRD) - century should be 00 not 100
!                 05 Sep 2000 (RRD) - fortran90 upgrade
!                 09 Sep 2002 (RRD) - fortran coding standards
!
! USAGE:  CALL TMPLUS(IY,IM,ID,IH,IC)
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

SUBROUTINE TMPLUS(IY,IM,ID,IH,IC)

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  INTEGER,  INTENT(INOUT) :: iy,im,id      ! current date
  INTEGER,  INTENT(INOUT) :: ih            ! current hour
  INTEGER,  INTENT(IN)    :: ic            ! increment in hours

!-------------------------------------------------------------------------------

  INTEGER   :: NDM(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
  INTEGER   :: loop, iyear

!-------------------------------------------------------------------------------

! new century year should be 00 not 100
  IY=MOD(IY,100)

! convert to 4 digit year
  IYEAR=MOD(IY+99,100)+1901

! check for leap year (doesn't account for century leap years)
  NDM(2)=28
  IF(MOD(IYEAR,4).EQ.0)NDM(2)=29

!-------------------------------------------------------------------------------
! current time (min) incremented by IC (+/-)
!-------------------------------------------------------------------------------

  IF(IC.EQ.0)RETURN

! loop through by hour
  DO LOOP=1,ABS(IC)
     IH=IH+SIGN(1,IC)

!-------------------------------------------------------------------------------
!    forward clock
!-------------------------------------------------------------------------------

     IF(IC.GT.0)THEN
        IF(IH.GE.24)THEN
           IH=IH-24
           ID=ID+1
           IF(ID.GT.NDM(IM))THEN
              ID=1
              IM=IM+1
              IF(IM.GT.12)THEN
                 IM=1
                 IY=IY+1
                 IY=MOD(IY,100)
              END IF
           END IF
        END IF

!-------------------------------------------------------------------------------
!    backward clock
!-------------------------------------------------------------------------------

     ELSE
        IF(IH.LT.0)THEN
           IH=24+IH
           ID=ID-1
           IF(ID.LT.1)THEN
              IM=IM-1
              IF(IM.LT.1)THEN
                 IM=12
                 IY=IY-1
              END IF
              ID=NDM(IM)
           END IF
        END IF
     END IF

  END DO

END SUBROUTINE tmplus
