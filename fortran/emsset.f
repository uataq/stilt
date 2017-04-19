!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  EMSSET           EMiSsion SET data entry
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   EMISSION TERM SET IS THE DATA ENTRY FOR EMISSION RATES, POLLUTANT
!   INFORMATION, STARTING TIME AND RELEASE DURATION.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 01 Apr 1997 (RRD)
!                 05 Jul 2000 (RRD) - fixed relative emission start
!                 04 Sep 2000 (RRD) - fortran90 upgrade
!                 30 Jul 2001 (RRD) - added option for temporal emission file
!                 24 Oct 2001 (RRD) - simplified backward dispersion options
!                 11 Feb 2002 (RRD) - formatted STARTUP file
!                 12 Jul 2002 (RRD) - input data decoder
!                 15 Dec 2003 (RRD) - space/comma delimted format on emit file
!                 02 Apr 2004 (RRD) - generic file unit numbers
!                 29 Nov 2004 (RRD) - four digit year in efile tm call
!                 07 Mar 2006 (RRD) - moved file I/O to emsprt
!
! USAGE:  CALL EMSSET(DIRT,NUMTYP,IUNIT,IBYR,IBMO,IBDA,IBHR,OLAT,OLON,BACK)
!
!   INPUT ARGUMENT LIST:     see below
!   OUTPUT ARGUMENT LIST:    see below
!   INPUT FILES:             unit 5 or unit KF21 if input from file CONTROL
!   OUTPUT FILES:            unit KF22 if input from unit 5
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE EMSSET(DIRT,NUMTYP,IUNIT,IBYR,IBMO,IBDA,IBHR,OLAT,OLON,BACK)

  USE funits

  IMPLICIT NONE

  INCLUDE 'DEFCONC.INC'         ! pollutant and concentration grid

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  TYPE(pset),    INTENT(OUT)   :: dirt(:)   ! for each pollutant type 
  INTEGER,       INTENT(IN)    :: numtyp    ! number of pollutant types
  INTEGER,       INTENT(IN)    :: iunit     ! unit number for input data
  INTEGER,       INTENT(IN)    :: ibyr,ibmo ! starting date
  INTEGER,       INTENT(IN)    :: ibda,ibhr ! starting time
  REAL,          INTENT(INOUT) :: olat,olon ! starting location
  LOGICAL,       INTENT(IN)    :: back      ! backward integration flag

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  INTEGER       :: kk,macc,kret

!-------------------------------------------------------------------------------

  INTERFACE
  SUBROUTINE TM2DAY(MACM,IY,IM,ID,IH,MN)
  IMPLICIT NONE
  INTEGER,   INTENT(IN)    :: macm            ! accumulate minutes
  INTEGER,   INTENT(OUT)   :: iy,im,id,ih,mn  ! current date/time
  END SUBROUTINE tm2day
!-------------------------------------------------------------------------------
  SUBROUTINE TM2MIN(IY,IM,ID,IH,MN,MACC)
  IMPLICIT NONE
  INTEGER,  INTENT(IN)   :: iy,im,id,ih,mn       ! date and time
  INTEGER,  INTENT(OUT)  :: macc                 ! minutes since 1 Jan 1970
  END SUBROUTINE tm2min
!-------------------------------------------------------------------------------
  SUBROUTINE DECODI(IUNIT,VAR1,VAR2,VAR3,VAR4,VAR5)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: IUNIT   ! unit number
  INTEGER, OPTIONAL, INTENT(INOUT) :: VAR1
  INTEGER, OPTIONAL, INTENT(INOUT) :: VAR2
  INTEGER, OPTIONAL, INTENT(INOUT) :: VAR3
  INTEGER, OPTIONAL, INTENT(INOUT) :: VAR4
  INTEGER, OPTIONAL, INTENT(INOUT) :: VAR5
  END SUBROUTINE decodi
  END INTERFACE

!-------------------------------------------------------------------------------
! generic defaults

  DIRT(1)%IDENT='????'
  DIRT(1)%QRATE=1.0
  DIRT(1)%QHRS=1.0
  DIRT(1)%START%YR=IBYR
  DIRT(1)%START%MO=IBMO
  DIRT(1)%START%DA=IBDA
  DIRT(1)%START%HR=IBHR
  DIRT(1)%START%MN=0

  DO KK=1,NUMTYP

!    multiple pollutants copy over default values
     IF(KK.GT.1)THEN
        DIRT(KK)%IDENT=DIRT(KK-1)%IDENT
        DIRT(KK)%QRATE=DIRT(KK-1)%QRATE
        DIRT(KK)%QHRS=DIRT(KK-1)%QHRS
        DIRT(KK)%START%YR=DIRT(KK-1)%START%YR
        DIRT(KK)%START%MO=DIRT(KK-1)%START%MO
        DIRT(KK)%START%DA=DIRT(KK-1)%START%DA
        DIRT(KK)%START%HR=DIRT(KK-1)%START%HR
        DIRT(KK)%START%MN=DIRT(KK-1)%START%MN
     END IF

!    pollutant ID used for map labels and optional chemistry
     IF(IUNIT.EQ.5)THEN
        WRITE(*,*)'Pollutant 4-Character Identification'
        WRITE(*,'(A4)')DIRT(KK)%IDENT
     END IF
     READ(IUNIT,'(A4)')DIRT(KK)%IDENT
     IF(IUNIT.EQ.5)WRITE(KF22,'(A4)')DIRT(KK)%IDENT

!    emission rate can be replaced if value set with lat/lon
     IF(IUNIT.EQ.5)THEN
        WRITE(*,*)'Emission rate (per hour)'
        WRITE(*,*)DIRT(KK)%QRATE
     END IF
     READ(IUNIT,*)DIRT(KK)%QRATE
     IF(IUNIT.EQ.5)WRITE(KF22,'(E10.3)')DIRT(KK)%QRATE

!    duration of emission applies to one emission cycle
     IF(IUNIT.EQ.5)THEN
        WRITE(*,*)'Hours of emission'
        WRITE(*,*)DIRT(KK)%QHRS
     END IF
     READ(IUNIT,*)DIRT(KK)%QHRS
     IF(IUNIT.EQ.5)WRITE(KF22,'(F10.3)')DIRT(KK)%QHRS
     IF(BACK)DIRT(KK)%QHRS=-ABS(DIRT(KK)%QHRS)

!    emission start can be absolute or relative to file start time
     IF(IUNIT.EQ.5)THEN
        WRITE(*,*)'Release start time: year month day hour minute'
        WRITE(*,*)DIRT(KK)%START%YR, DIRT(KK)%START%MO,                    &
                  DIRT(KK)%START%DA, DIRT(KK)%START%HR, DIRT(KK)%START%MN
     END IF


!    READ(IUNIT,*) DIRT(KK)%START%YR, DIRT(KK)%START%MO,                   &
!                  DIRT(KK)%START%DA, DIRT(KK)%START%HR, DIRT(KK)%START%MN
!    replace fixed format read with variable read 
     CALL DECODI(IUNIT,DIRT(KK)%START%YR, DIRT(KK)%START%MO,               &
                 DIRT(KK)%START%DA, DIRT(KK)%START%HR, DIRT(KK)%START%MN)

!    if month=0 then convert relative time to absolute time
     IF(DIRT(KK)%START%MO.EQ.0)THEN
        DIRT(KK)%START%YR=IBYR
        DIRT(KK)%START%MO=IBMO
        DIRT(KK)%START%DA=DIRT(KK)%START%DA+IBDA
        DIRT(KK)%START%HR=DIRT(KK)%START%HR+IBHR

!       adjust relative date for potential month crossing error
        CALL TM2MIN(DIRT(KK)%START%YR,DIRT(KK)%START%MO,DIRT(KK)%START%DA,     &
                    DIRT(KK)%START%HR, DIRT(KK)%START%MN, MACC)
        CALL TM2DAY(MACC,DIRT(KK)%START%YR,DIRT(KK)%START%MO,DIRT(KK)%START%DA,&
                    DIRT(KK)%START%HR, DIRT(KK)%START%MN)
     END IF

     IF(IUNIT.EQ.5)WRITE(KF22,'(5I3)') DIRT(KK)%START%YR, DIRT(KK)%START%MO,   &
                   DIRT(KK)%START%DA, DIRT(KK)%START%HR, DIRT(KK)%START%MN

     CALL TM2MIN(DIRT(KK)%START%YR, DIRT(KK)%START%MO, DIRT(KK)%START%DA,      &
                 DIRT(KK)%START%HR, DIRT(KK)%START%MN, DIRT(KK)%START%MACC)
  END DO

END SUBROUTINE emsset
