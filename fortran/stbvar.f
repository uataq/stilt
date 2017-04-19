!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  STBVAR           STaBility VARiance methods for mixing
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   STABILITY VARIANCE COMPUTES THE VERTICAL AND HORIZONTAL VELOCITY
!   VARIANCES FROM SURFACE STABILITY PARAMETERS AND BULK RICHARDSON
!   NUMBER GIVEN MONIN-OBUKHOV LENGTH AND FRICTION VELOCITY. 
!   FOR MORE INFORMATION ON EQUATIONS IN THIS ROUTINE SEE ...
!   KANTHA AND CLAYSON, 2000, SMALL SCALE PROCESSES IN GEOPHYSICAL
!   FLUID FLOWS, VOL 67, INTERNATIONAL GEOPHYSICS SERIES, ACADEMIC 
!   PRESS, 883 PP; GARRATT, 1992, THE ATMOSPHERIC BOUNDARY LAYER,
!   CAMBRIDGE UNIVERSITY PRESS, 316 PP.
!
!   The November 2006 revision replaced stable/neutral ratios with those
!   from Lumely & Panofsky: 3.0,4.0,4.5 replaced with 1.7,4.0,5.0)
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 02 Dec 2003 (RRD) - initial version from stbsnd
!                  11 Dec 2003 (RRD) - zero velocity test
!                  21 Nov 2006 (RRD) - lumely & panofsky
!                  04 Jun 2008 (RRD) - moved variables to common
!                  15 Aug 2008 (RRD) - enhanced options
!                  01 Oct 2008 (RRD) - mixing adjustments
!
! USAGE:  CALL STBVAR(TKERD,TKERN,KZMIX,TVMIX,KSFC,NL,UU,VV,TT,ZZ,EE,HH,XX)
!
!   INPUT ARGUMENT LIST:       see below
!   OUTPUT ARGUMENT LIST:      see below
!   INPUT FILES:               none
!   OUTPUT FILES:              none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE STBVAR(TKERD,TKERN,KZMIX,TVMIX,KSFC,NL,UU,VV,TT,ZZ,EE,HH,XX)

  USE stbcon

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  INTEGER,   INTENT(IN)    :: tkerd     ! day vert to horiz turb ratio
  INTEGER,   INTENT(IN)    :: tkern     ! night vert to horiz turb ratio
  INTEGER,   INTENT(IN)    :: kzmix     ! averaged vertical mixing         
  REAL,      INTENT(IN)    :: tvmix     ! tropospheric mixing scale factor
  INTEGER,   INTENT(IN)    :: ksfc      ! index of top of surface layer
  INTEGER,   INTENT(IN)    :: nl        ! number of sigma levels
  REAL,      INTENT(IN)    :: uu (:)    ! horizontal wind component
  REAL,      INTENT(IN)    :: vv (:)    ! horizontal wind component
  REAL,      INTENT(IN)    :: tt (:)    ! virtual potential temperature (pot K)
  REAL,      INTENT(IN)    :: zz (:)    ! height at levels (m)
  REAL,      INTENT(INOUT) :: ee (:)    ! in TKE (J/kg), output  - V'2 (m2/s2)
  REAL,      INTENT(INOUT) :: hh (:)    ! U component turbulence - U'2 (m2/s2) 
  REAL,      INTENT(INOUT) :: xx (:)    ! vertical turbulence    - W'2 (m2/s2)  

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  INTEGER    :: k,kpbl
  REAL       :: wstr2,ustr2,dtdz,delz,delu,zsfc,vmix,etrm,phih 
  REAL       :: uvp2,uvm2,xl,vl,ri,sl,vscale,hscale

!                             stability analysis results each grid point
  REAL       :: fvel         ! scalar friction velocity (m/s)
  REAL       :: ustr         ! u- friction velocity (m/s)
  REAL       :: vstr         ! v- friction velocity (m/s)
  REAL       :: tstr         ! friction temperature (deg K)
  REAL       :: wstr         ! convective velocity scale (m/s)
  REAL       :: slen         ! Obukhov stability length (m)
  REAL       :: zmix         ! mixed layer depth (m)
  REAL       :: psi          ! integrated stability function heat

  COMMON /stbcom/ fvel,ustr,vstr,tstr,wstr,slen,zmix,psi
  COMMON /stblen/ vscale,hscale

! top of surface layer (minimum defined for data)
  ZSFC=MAX(ZZ(KSFC),0.1*ZMIX)

! friction velocity squared
  USTR2=FVEL*FVEL

! convective velocity scale squared
  WSTR2=WSTR*WSTR

  DO K=1,NL

!-------------------------------------------------------------------------------
!    level within the surface layer

     IF(NINT(ZZ(K)).LE.NINT(ZSFC))THEN
        VMIX=-1.0               ! mixing coeffience not required

        IF(SLEN.GE.0.0)THEN
!          stable layer 
           XX(K)=1.7*USTR2      ! w'2
           HH(K)=4.0*USTR2      ! u'2
           EE(K)=5.0*USTR2      ! v'2

        ELSE
!          unstable layer 
           SL=ZZ(K)/SLEN
           XX(K)=1.74*USTR2*(1.0-3.0*SL)**0.67
           HH(K)=0.36*WSTR2
           EE(K)=0.36*WSTR2
        END IF
        KPBL=K

!-------------------------------------------------------------------------------
!    level is within the PBL

     ELSEIF(ZZ(K).LT.ZMIX)THEN
        VMIX=-1.0               ! mixing coeffience not required

        SL=ZZ(K)/ZMIX
        IF(SLEN.GE.0.0)THEN
!          stable layer (to the 3/2 power)
           XL=SQRT((1.0-SL)*(1.0-SL)*(1.0-SL))
           XX(K)=1.7*USTR2*XL   ! w'2
           HH(K)=4.0*USTR2*XL   ! u'2
           EE(K)=5.0*USTR2*XL   ! v'2

        ELSE
!          unstable layer [1.17 = (1+0.5*R**0.67), where R=0.2] 
           XX(K)=1.17*WSTR2*(SL-SL*SL)**0.67
           HH(K)=0.36*WSTR2
           EE(K)=0.36*WSTR2
        END IF
        KPBL=K

!-------------------------------------------------------------------------------
!    compute mixing through the inversion layer (only for convective case)

     ELSEIF(NINT(ZZ(K)).EQ.NINT(ZMIX).AND.WSTR.GT.0.0)THEN

!       Betts and Beljaars inversion layer jump model
        DTDZ=(TT(K)-TT(K-1))/(ZZ(K)-ZZ(K-1))
        DTDZ=MAX(DTDZ,0.1)
        VMIX=-0.4*TSTR*FVEL/DTDZ

!-------------------------------------------------------------------------------
!    level is within the free troposphere

     ELSE

!       bulk Richardson number
        DELZ=ZZ(K)-ZZ(K-1)
        DELU=(UU(K)-UU(K-1))**2+(VV(K)-VV(K-1))**2
        DELU=MAX(DELU,0.1)
        RI=GRAV*DELZ*(TT(K)-TT(K-1))/DELU/TT(K-1)
        RI=MAX(0.0, MIN(20.0, RI))

!       vertical length scale (l)
        VL=1.0/(1.0/VONK/ZZ(K)+1.0/150.0)

!       stability scale (l/lo)
        IF(RI.GE.0.0.AND.RI.LE.0.001)THEN
           XL=1.0893*RI
        ELSE
           XL=A1+RI*(A2+RI*(A3+RI*(A4+A5*RI)))
        END IF

!       local stability function
        ETRM=B*EXP(-D*XL)*(1.0+C-D*XL)
        PHIH=PRN*(1.0+(A*SQRT(1.0+A*B*XL)+ETRM)*XL)
        VMIX=VL*VL*ABS(SQRT(DELU)/DELZ)/PHIH
     END IF
  
!-------------------------------------------------------------------------------
!    proportion horizontal turbulence according to the mean velocity components

     IF(VMIX.GE.0.0)THEN
!       check limits
        VMIX=MAX(VMIN, MIN(VMAX, VMIX))
!       Convert mixing coefficient to turbulent velocity variance 
        XX(K)=VMIX/VSCALE
!       Assume W'2 = U'2 + V'2
        HH(K)=0.5*XX(K)
        EE(K)=0.5*XX(K) 
     END IF

!    horizontal turbulent values computed from one of the sections above
     UVP2=HH(K)+EE(K)
     UVM2=UU(K)*UU(K)+VV(K)*VV(K)

!    turbulent vectors proportional to velocity vectors U'=U2(U'+V')/(U2+V2)
     IF(UVM2.GT.0.0)THEN
        HH(K)=UU(K)*UU(K)*UVP2/UVM2  
        EE(K)=VV(K)*VV(K)*UVP2/UVM2  
     ELSE
        HH(K)=0.0
        EE(K)=0.0
     END IF

! level loop
  END DO

!-------------------------------------------------------------------------------
! PBL mixing profile adjustments

  IF(KZMIX.EQ.1)THEN
!    replace PBL with vertical average
     VMIX=SUM(XX(1:KPBL))/KPBL
     DO K=1,KPBL
        XX(K)=VMIX
     END DO
  ELSEIF(KZMIX.EQ.2)THEN
!    boundary layer adjustment
     DO K=1,KPBL
        XX(K)=XX(K)*TVMIX
     END DO
  ELSEIF(KZMIX.EQ.3)THEN
!    free troposphere adjustment
     KPBL=MIN(KPBL+1,NL)
     DO K=KPBL,NL
        XX(K)=XX(K)*TVMIX
     END DO
  ELSE
     CONTINUE
  END IF

!-------------------------------------------------------------------------------
! urban area enhanced turbulence is defined when more of the turbulence
! is assigned to the vertical component at night (tkern > tkerd) than during
! the daytime. The enhancement is determined from DCNet versus NAM TKE values
! and only applies during stable-neutral conditions.

  IF(TKERN.GT.TKERD.AND.SLEN.GE.0.0) THEN
     EE=1.40*EE
     HH=1.40*HH
     XX=1.40*XX
  END IF

END SUBROUTINE stbvar
