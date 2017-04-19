!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  SUNAVE           SUN ANGle returns the solar angle
!   PRGMMR:    Christoph Gerbig ORG: MPI-BGC      DATE:16 Jan 09
!
! ABSTRACT:  THIS CODE WRITTEN AT MPI-BGC ...
!   SUN AVERAGING - RETURNS THE TIME AVERAGED SHORTWAVE SOLAR 
!   RADIATION
!   GIVEN THE LAT/LON, TIME, AND INTEGRATION TIME.
!
!
! USAGE:  CALL SUNAVE(JET,JETDEL,OLAT,OLON,SWF)
!   INPUT ARGUMENT LIST:
!     JET   - int time at end of window, elapsed minutes since January 1st 1970
!     JETDEL   - int time interval in minutes
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
! $Id: sunave.f,v 1.1 2009/10/26 15:36:55 jel Exp $
!
!$$$
      SUBROUTINE SUNAVE(JET,JETDEL,OLAT,OLON,SWF)

!dwen      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT none
!dwen      INTEGER JETDEL
!dwen(20090917)********************
      integer, intent(in) :: jet,jetdel
      real,    intent(in) :: olat,olon
      real,    intent(out):: swf

      real    :: dtyr,df,snoon,ha,hr,sig,sd,d,  &
                 sea,ea,has,seas
      integer :: jyr,jet1,jet0,jmo,jda,jhr,jmn

!     radians per degree
!dwen      DATA RDPDG/0.0174532925199433d0/
      real,parameter:: RDPDG=0.0174532925199433

!     cloud top solar const (w/m2)
!dwen      DATA SOLC/1104.0/
      real,parameter:: SOLC=1104.0

      interface
      subroutine tm2min(iy,im,id,ih,mn,macc)
      implicit none
      integer,intent(in) :: iy,im,id,ih,mn
      integer,intent(out):: macc
      end subroutine tm2min
      end interface
!************************************

!dwen      CALL TM2DAY(JET,JYR,IM,ID,IH,MN)
      CALL TM2DAY(JET,JYR,jMo,jDa,jHr,jMN)

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

!     hour angle at end time
      HA=(HR-SNOON)*15.0+OLON

!     hour angle at start time
      HAS=(HR-JETDEL/60.0-SNOON)*15.0+OLON

!     sine of the solar elevation angle
!      SEA=COS(HA*RDPDG)*COS(OLAT*RDPDG)*COS(D)+SIN(OLAT*RDPDG)*SD

!     time integral of solar elevation angle at end time
      SEA=SIN(HA*RDPDG)*COS(OLAT*RDPDG)*COS(D)/RDPDG+SIN(OLAT*RDPDG)*SD*HA

!     solar elevation angle in degrees
!     time integral of solar elevation angle at end time
      SEAS=SIN(HAS*RDPDG)*COS(OLAT*RDPDG)*COS(D)/RDPDG+SIN(OLAT*RDPDG)*SD*HAS

!     time averaged sine of the solar elevation angle
      SEA=(SEA-SEAS)/(JETDEL*15.0/60.0)

      EA=ASIN(SEA)/RDPDG

      SWF=aMAX1(0.,SEA*SOLC)

      RETURN
      END
