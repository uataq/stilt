SUBROUTINE P10ADJ(P10F,ZWND,UBAR,RAIN,AREA,MASS,PGRD)

!---------------------------------------------------------------
! PM10 emission adjustment based upon local meteorology
! Algorithms based upon Draxler et al.,2001,Atm Environ,35:4315-4330
! Set control file: diameter (3.0 um), density (2.5 g/cc), shape (1.0)
!---------------------------------------------------------------
! Last Revised: 01 Nov 1999
!               10 Apr 2000 (RRD) - emission adjustment factors
!               22 Mar 2002 (RRD) - updated links with hysplit
!               25 Oct 2005 (RRD) - threshold velocity factor
!---------------------------------------------------------------

  IMPLICIT NONE

  REAL,   INTENT(IN)    :: P10F ! threshold emission sensitivity factor
  REAL,   INTENT(IN)    :: ZWND ! release height or ref meteo height
  REAL,   INTENT(IN)    :: UBAR ! wind speed in km/min
  REAL,   INTENT(IN)    :: RAIN ! precip rate
  REAL,   INTENT(IN)    :: AREA ! emission area in sq meters
  REAL,   INTENT(INOUT) :: MASS ! particle/puff mass
  INTEGER,INTENT(INOUT) :: PGRD ! grid identification index

  INTEGER :: LT
  REAL    :: US,UT,UVEL,ZLOG,RCD,FLUX,ROW,GRAV,VONK

  REAL :: Z0NS(7)   ! flux equation roughess
  REAL :: USTR(7)   ! flux equation friction-v 
  REAL :: TKON(7)   ! flux equation emission factor

! assigned roughness length (m) by land use class
  DATA Z0NS/.0002,.0004,.00002,.0005,.003,.00002,.0002/
! assigned threshold friction velocity (m/s) by land use class
  DATA USTR/1.,0.622,0.28,0.685,3.49,3.,0.75/
! soil texture relates soil type to pm10 emission flux
  DATA TKON/2*5.6E-04,  5.6E-04,  4*5.6E-04/

! constants: density(g/m^3), gravity (m/s^2), VonKarman
  DATA       ROW/1225.0/,    GRAV/9.8/,       VONK/0.4/

!---------------------------------------------------------------------
! set all cells to one soil type (active sand sheet)
  LT=3

! convert wind from km/min to m/sec
  UVEL=1000.0*UBAR/60.0
 
! threshold wind speed
  ZLOG=ALOG(ZWND/Z0NS(LT))/VONK
  UT=P10F*USTR(LT)*ZLOG

! drag coefficient
  IF(UVEL.LE.UT)THEN  
     RCD=1.0/ZLOG 
  ELSE
     RCD=1.0/ZLOG+0.003*(1.0-UT/UVEL)
  END IF

!---------------------------------------------------------------------
! determine friction velocity

! local friction velocity based on windspeed
  US=RCD*UVEL
! local friction velocity from advmet subroutine
! US=METO%USTR

!---------------------------------------------------------------------
! pm10 emission flux equations

  IF(UVEL.GT.UT)THEN
!    Gillette equation: vertical mass flux (g m-2 sec-1)
!    FLUX=AMAX1(0.0,TKON(LT)*ROW*US*(US*US-USTR(LT)*USTR(LT))/GRAV)

!    Eq. 14 - Westphal, Toon, & Carlson, 1987, J. Geophys Res, 92:3027-3049
!    FLUX[g/cm2]=1.0E-14*US[cm/s]**4.0
!    FLUX[g/m2] =1.0E-14*(US[m/s]*100.0)**4.0*1.0E+04

     FLUX=0.01*US**4.0
  ELSE
     FLUX=0.0
  END IF

! Precipitation adjustment
  IF(RAIN.GT.0.0)FLUX=0.0

! max limit at 1 mg/m2/s
  FLUX=AMIN1(0.001,FLUX)

! Control file should be configured with .\bdyfiles\dustbdy
! mass is set in emspnt to equal (dt/60)*(qtrm/npar)
! emission duration is not a factor in the computation

  IF(FLUX.GT.0.0)THEN
!    convert to mass = (s/h) (g/m2-s) (m2) (h/part)
     MASS=3600.0*FLUX*AREA*ABS(MASS)

  ELSE
!    zero emissions and index, particle removed in main
     MASS=0.0
!    particle removed in main when grid# is zero
     PGRD=0
  END IF
  
END SUBROUTINE p10adj
