!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  TM2DAY           TiMe2DAY converts min since 1970 to time
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   TIME2DAY CONVERTS ACCUMULATED MINUTES SINCE BEGINNING OF
!   1970 TO TIME AS YEAR, MONTH, DAY, HOUR, MINUTE.
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 10 Apr 1998 (RRD)
!                  07 May 1998 (RRD) - y2k mod
!                  02 Sep 2000 (RRD) - fortran90 upgrade
!                  02 Feb 2001 (RRD) - expanded calendar 1900-2100
!
! USAGE:  CALL TM2DAY(MACM,IY,IM,ID,IH,MN)
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

SUBROUTINE TM2DAY(MACM,IY,IM,ID,IH,MN)

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  INTEGER,   INTENT(IN)    :: macm            ! accumulate minutes
  INTEGER,   INTENT(OUT)   :: iy,im,id,ih,mn  ! current date/time

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

!             default number of accumulated days in each month (non leap year)
  INTEGER  :: nadpm (12) = (/0,31,59,90,120,151,181,212,243,273,304,334/)

  INTEGER  :: k,kk,macc,iyear,mhr,mda,mmo,idt

! accumulated days table from Jan 1, 1900 where the first value is Jan 1, 1901
  INTEGER  :: NADPY (200)  
  COMMON /TMVALS/ NADPY

!-------------------------------------------------------------------------------

! check for initialization
  IF(NADPY(1).NE.365)THEN
     WRITE(*,*)'*ERROR* tm2day: calendar routines not initialized (tminit)'
     STOP 900
  END IF

! compute four digit year
  KK=1
  DO WHILE (NADPY(KK)*1440.LT.MACM)
     KK=KK+1
  END DO
  IYEAR=1900+(KK-1)

! compute minutes in this year
  MACC=MACM-NADPY(KK-1)*1440

! convert back to two digit year
  IY=MOD(IYEAR-1800,100)

! current minute
  MHR=MACC/60
  MN=MACC-MHR*60

! current hour
  MDA=MHR/24
  IH=MHR-MDA*24

! current month and day
  DO K=1,12
     IDT=NADPM(K)

!    adjust accumulated days for leap year
     IF(MOD(IYEAR,4).EQ.0.AND.K.GT.2)IDT=IDT+1

     IF(IDT.LE.MDA)THEN
        MMO=IDT
        IM=K
     END IF
  END DO
  ID=MDA-MMO+1

! check for end-of-year
  IF(ID.GT.31)THEN
     ID=ID-31
     IM=1
     IY=MOD(IY+1,100)
  END IF

END SUBROUTINE tm2day
