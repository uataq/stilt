!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  SUNANG           SUN ANGle returns the solar angle
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   SUN ANGLE - RETURNS THE SOLAR
!   ELEVATION ANGLE (90-ZENITH ANGLE) AND THE SINE OF THE ANGLE,
!   GIVEN THE LAT/LON AND TIME.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 08 Apr 1997 - RRD
!
! USAGE:  CALL SUNCLR(JYR,JET,OLAT,OLON,EA,SEA)
!   INPUT ARGUMENT LIST:
!     JET   - int elapsed minutes since January 1st 1970
!     OLAT  - real      latitude in degrees and fraction (+ = north)
!     OLON  - real      longitude in degrees anb fraction (- = west)
!   OUTPUT ARGUMENT LIST:
!     SWF   - clear sky short wave radiative flux
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: sunclr.f,v 1.1 2009/10/26 15:36:55 jel Exp $
!
!$$$
! CHG:(9/24/02) copied from sunang.f, to compute clear sky DSWF
      SUBROUTINE SUNCLR(JET,OLAT,OLON,SWF)

!dwen      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT none

!dwen(20090916)********************
      integer,  intent(in):: jet
      real,     intent(in):: olat
      real,     intent(in):: olon
      real,     intent(out):: swf

      real       :: dtyr,df,snoon,ha,hr,sig,sd,d,  &
                    sea,ea
      integer    :: jyr,jet1,jet0,jmo,jda,jhr,jmn 
!     radians per degree
      real,parameter:: rdpdg = 0.0174532925199433
!     cloud top solar const (w/m2)
      real,parameter:: solc  = 1104.0

      interface
      subroutine tm2min(iy,im,id,ih,mn,macc)
      implicit none
      integer,  intent(in) :: iy,im,id,ih,mn
      integer,  intent(out):: macc
      end subroutine tm2min
      end interface
!************************************

!     radians per degree
!      DATA RDPDG/0.0174532925199433d0/

! CHG(9/24/02)
!     cloud top solar const (w/m2)
!      DATA SOLC/1104.0/

! CHG(9/24/02) get current year from JET
      CALL TM2DAY(JET,JYR,jmo,jda,jhr,jmn)
!     adjust for leap year
      DTYR=365.242
      IF(MOD(JYR+1900,4).EQ.0)DTYR=366.242

!     compute elapsed time for this year
      CALL TM2MIN(JYR,1,1,0,0,JET0)
      JET1=JET-JET0

!     derive fractional part of year
      DF=360.0*JET1/(DTYR*24.0*60.0)

      SIG=279.9348+DF+1.914827*SIN(RDPDG*DF)-0.079525*COS(RDPDG*DF)     &
     &    +0.019938*SIN(2.0*RDPDG*DF)-0.00162*COS(2.0*RDPDG*DF)

!     solar declination at noon
      SD=SIN(23.4438333*RDPDG)*SIN(RDPDG*SIG)
      D=ASIN(SD)

!     time of meridian passage (solar noon time)
      SNOON=12.0+0.12357*SIN(RDPDG*DF)-0.004289*COS(RDPDG*DF)           &
     &      +0.153809*SIN(2.0*RDPDG*DF)+0.060783*COS(2.0*RDPDG*DF)

!     current hour
      HR=AMOD(JET1/60.0,24.0)

!     hour angle at current time
      HA=(HR-SNOON)*15.0+OLON

!     sine of the solar elevation angle
      SEA=COS(HA*RDPDG)*COS(OLAT*RDPDG)*COS(D)+SIN(OLAT*RDPDG)*SD

!     solar elevation angle in degrees
      EA=ASIN(SEA)/RDPDG

! CHG(9/24/02) get clear sky swf
      SWF=aMAX1(0.0,SEA*SOLC)

      RETURN
      END
