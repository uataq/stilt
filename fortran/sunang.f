!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  SUNANG           SUN ANGle returns the solar angle
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   SUN ANGLE - RETURNS THE SOLAR
!   ELEVATION ANGLE (90-ZENITH ANGLE) AND THE SINE OF THE ANGLE,
!   GIVEN THE LAT/LON AND TIME.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 08 Apr 1997 (RRD)
!                 05 Sep 2000 (RRD) - fortran90 upgrade
!                 27 Feb 2002 (RRD) - dropped year from argument
!                 09 Sep 2002 (RRD) - fortran coding standards
!
! USAGE:  CALL SUNANG(JET,OLAT,OLON,EA,SEA)
!
!   INPUT ARGUMENT LIST:   see below
!   OUTPUT ARGUMENT LIST:  see below
!   INPUT FILES:           none
!   OUTPUT FILES:          none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE SUNANG(JET,OLAT,OLON,EA,SEA)

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER,  INTENT(IN)   :: jet       ! elapsed minutes 
  REAL,     INTENT(IN)   :: olat      ! latitude (+ = north)
  REAL,     INTENT(IN)   :: olon      ! longitude (- = west)
  REAL,     INTENT(OUT)  :: ea        ! solar elevation angle in deg at time
  REAL,     INTENT(OUT)  :: sea       ! sine of EA

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  REAL,  PARAMETER  :: rdpdg  = 0.0174532925199433     ! radians per degree
  REAL    :: d,sd,sig,ha,hr,snoon,df,dtyr
  INTEGER :: jet1,jet0,jyr,jmo,jda,jhr,jmn

  INTERFACE
!-------------------------------------------------------------------------------
  SUBROUTINE TM2MIN(IY,IM,ID,IH,MN,MACC)
  IMPLICIT NONE
  INTEGER,  INTENT(IN)   :: iy,im,id,ih,mn       ! date and time
  INTEGER,  INTENT(OUT)  :: macc                 ! minutes since 1 Jan 1970
  END SUBROUTINE tm2min
!-------------------------------------------------------------------------------
  END INTERFACE

! determine current year (02/27/2002)
  CALL TM2DAY(JET,JYR,JMO,JDA,JHR,JMN)

! adjust for leap year
  DTYR=365.242
  IF(MOD(JYR+1900,4).EQ.0)DTYR=366.242

! compute elapsed time for this year
  CALL TM2MIN(JYR,1,1,0,0,JET0)
  JET1=JET-JET0

! derive fractional part of year
  DF=360.0*FLOAT(JET1)/(DTYR*24.0*60.0)

  SIG=279.9348+DF+1.914827*SIN(RDPDG*DF)-0.079525*COS(RDPDG*DF)            &
      +0.019938*SIN(2.0*RDPDG*DF)-0.00162*COS(2.0*RDPDG*DF)

! solar declination at noon
  SD=SIN(23.4438333*RDPDG)*SIN(RDPDG*SIG)
  D=ASIN(SD)

! time of meridian passage (solar noon time)
  SNOON=12.0+0.12357*SIN(RDPDG*DF)-0.004289*COS(RDPDG*DF)                  &
        +0.153809*SIN(2.0*RDPDG*DF)+0.060783*COS(2.0*RDPDG*DF)

! current hour
  HR=MOD(JET1/60.0,24.0)

! hour angle at current time
  HA=(HR-SNOON)*15.0+OLON

! sine of the solar elevation angle
  SEA=COS(HA*RDPDG)*COS(OLAT*RDPDG)*COS(D)+SIN(OLAT*RDPDG)*SD

! solar elevation angle in degrees
  EA=ASIN(SEA)/RDPDG

END SUBROUTINE sunang
