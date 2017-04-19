!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  TM2MIN           TiMe2MINute converts time to min since 1970
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   TIME2MINUTE CONVERTS TIME AS EXPRESSED IN DATE FORM:
!   YEAR, MONTH, DAY, HOUR, MINUTE  TO ACCUMULATED MINUTES
!   SINCE THE BEGINNING OF 1900
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 10 Apr 1998 (RRD) 
!                  07 May 1998 (RRD) - y2k mod
!                  03 Sep 2000 (RRD) - fortran90 upgrade
!                  02 Feb 2001 (RRD) - expanded calendar 1900-2100
!                  27 May 2003 (RRD) - enter with 2 or 4 digit year
!
! USAGE:  CALL TM2MIN(IY,IM,ID,IH,MN,MACC)
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

SUBROUTINE TM2MIN(IY,IM,ID,IH,MN,MACC)

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  INTEGER,  INTENT(IN)   :: iy,im,id,ih,mn       ! date and time
  INTEGER,  INTENT(OUT)  :: macc                 ! minutes accumulated 

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  INTEGER  :: iyear,idt

! default number of accumulated days in each month (non leap year)
  INTEGER  :: NADPM (12) = (/0,31,59,90,120,151,181,212,243,273,304,334/)

! accumulated days table from Jan 1, 1900 where the first value is Jan 1, 1901
  INTEGER  :: NADPY (200)
  COMMON /TMVALS/ NADPY

  SAVE NADPM

!-------------------------------------------------------------------------------

! check for initialization
  IF(NADPY(1).NE.365)THEN
     WRITE(*,*)'*ERROR* tm2min: calendar routines not initialized (tminit)'
     STOP 900
  END IF

! compute four digit year
  IF(IY.LE.100)THEN
     IYEAR=MOD(IY+99,100)+1901
  ELSEIF(IY.LT.1941)THEN
     WRITE(*,*)'*ERROR* tm2min: year out of valid range (1941-2040) - ',IY
     STOP 900
  ELSE
     IYEAR=IY
  END IF

! assume that no data prior to 1940 hence those years represents 2000-2040
  IF(IYEAR.LT.1940)IYEAR=IYEAR+100

! number of accumulated days until this year
  IDT=NADPY(IYEAR-1900)

! add accumulated days for this year
  IDT=IDT+NADPM(IM)+(ID-1)

! adjust accumulated days for leap year
! does not account for centuries where mod(iyear,400)=0
  IF(MOD(IYEAR,4).EQ.0.AND.IM.GT.2)IDT=IDT+1

! convert to minutes
  MACC=MN+(IH+IDT*24)*60

END SUBROUTINE tm2min
