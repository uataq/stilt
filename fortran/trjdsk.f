!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  TRJDSK           TRaJectory endpoints to DiSK
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   TRAJECTORY OUTPUT TO DISK
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 19 Jul 1996 (RRD)
!                  05 Sep 2000 (RRD) - fortran90 upgrade
!                  18 Mar 2002 (RRD) - expanded diagnostic format
!                  02 Apr 2004 (RRD) - generic file unit numbers
!                  31 May 2005 (RRD) - updated trajectory file format
!
! USAGE:  CALL TRJDSK(KP,KG,JET,MC0,IFHR,TLAT,TLON,TZHT,NDIAG,TDIAG)
!
!   INPUT ARGUMENT LIST:     see below
!   OUTPUT ARGUMENT LIST:    see below
!   INPUT FILES:             none
!   OUTPUT FILES:            unit KF11 to file as defined in input CONTROL file
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE TRJDSK(KP,KG,JET,MC0,IFHR,TLAT,TLON,TZHT,NDIAG,TDIAG)

  USE funits

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  INTEGER,  INTENT(IN) :: kp            ! trajectory id sequence number
  INTEGER,  INTENT(IN) :: kg            ! current meteo grid number
  INTEGER,  INTENT(IN) :: jet           ! elapsed minutes
  INTEGER,  INTENT(IN) :: mc0           ! initial minutes
  INTEGER,  INTENT(IN) :: ifhr          ! current forecast hour
  REAL,     INTENT(IN) :: tlat,tlon     ! trajectory position
  REAL,     INTENT(IN) :: tzht          ! trajectory height AGL
  INTEGER,  INTENT(IN) :: ndiag         ! number of diagnostic meteo variables
  REAL,     INTENT(IN) :: tdiag (ndiag) ! diagnostic meteo array

!-------------------------------------------------------------------------------

  REAL                 :: hours
  INTEGER              :: iyr,imo,ida,ihr,imn,kd

!-------------------------------------------------------------------------------
  INTERFACE
  SUBROUTINE TM2DAY(MACM,IY,IM,ID,IH,MN)
  IMPLICIT NONE
  INTEGER,   INTENT(IN)    :: macm            ! accumulate minutes
  INTEGER,   INTENT(OUT)   :: iy,im,id,ih,mn  ! current date/time
  END SUBROUTINE tm2day
  END INTERFACE
!-------------------------------------------------------------------------------

! convert accumulated minutes to clock time
  HOURS=(JET-MC0)/60.0
  CALL TM2DAY(JET,IYR,IMO,IDA,IHR,IMN)

  WRITE(KF11,'(8I6,F8.1,2F9.3,9(1X,F8.1))')   &
     KP,KG,IYR,IMO,IDA,IHR,IMN,IFHR,HOURS,TLAT,TLON,TZHT,(TDIAG(KD),KD=1,NDIAG)

END SUBROUTINE trjdsk
