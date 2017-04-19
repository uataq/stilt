!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  STBTKE           STaBility Turbulent Kinetic Energy
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   STABILITY TKE COMPUTES THE VERTICAL AND HORIZONTAL VELOCITY
!   VARIANCES FROM THE TURBULENT KINETC ENERGY FIELD PROVIDED 
!   WITH THE OTHER METEOROLOGICAL MODEL OUTPUT FIELS.
!   FOR MORE INFORMATION ON EQUATIONS IN THIS ROUTINE SEE ...
!   KANTHA AND CLAYSON, 2000, SMALL SCALE PROCESSES IN GEOPHYSICAL
!   FLUID FLOWS, VOL 67, INTERNATIONAL GEOPHYSICS SERIES, ACADEMIC 
!   PRESS, 883 PP; GARRATT, 1992, THE ATMOSPHERIC BOUNDARY LAYER,
!   CAMBRIDGE UNIVERSITY PRESS, 316 PP. 
!
!   The November 2006 revision replaced stable/neutral ratios with those 
!   from Lumely & Panofsky: 0.52,0.70,0.78 replaced with 0.32,0.74,0.85)
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 02 Dec 2003 (RRD) - initial version from stbsnd
!                  11 Dec 2003 (RRD) - zero velocity test
!                  22 Jul 2004 (RRD) - option to apply TKER at all levels
!                  21 Nov 2006 (RRD) - day night tke partition
!                  04 Jun 2008 (RRD) - moved variables to common
!                  15 Aug 2008 (RRD) - enhanced options selection
!                  01 Oct 2008 (RRD) - mixing adjustments
!
! USAGE:  CALL STBTKE(TKERD,TKERN,KZMIX,TVMIX,KSFC,NL,UU,VV,ZZ,EE,HH,XX)
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

SUBROUTINE STBTKE(TKERD,TKERN,KZMIX,TVMIX,KSFC,NL,UU,VV,ZZ,EE,HH,XX)

  USE stbcon

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  REAL,      INTENT(IN)    :: tkerd     ! day ratio of vert to horiz turb
  REAL,      INTENT(IN)    :: tkern     ! night ratio of vert to horiz turb
  INTEGER,   INTENT(IN)    :: kzmix     ! vertical mixing averaging flag   
  REAL,      INTENT(IN)    :: tvmix     ! tropospheric mixing scale factor
  INTEGER,   INTENT(IN)    :: ksfc      ! index of top of surface layer
  INTEGER,   INTENT(IN)    :: nl        ! number of sigma levels
  REAL,      INTENT(IN)    :: uu (:)    ! horizontal wind component
  REAL,      INTENT(IN)    :: vv (:)    ! horizontal wind component
  REAL,      INTENT(IN)    :: zz (:)    ! height at levels (m)
  REAL,      INTENT(INOUT) :: ee (:)    ! in TKE (J/kg), output  - V'2 (m2/s2)
  REAL,      INTENT(INOUT) :: hh (:)    ! U component turbulence - U'2 (m2/s2) 
  REAL,      INTENT(INOUT) :: xx (:)    ! vertical turbulence    - W'2 (m2/s2)  

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  INTEGER    :: k,kpbl
  REAL       :: wstr2,ustr2,zsfc,vmix 
  REAL       :: uvp2,uvm2,xl,sl,tker
  
  REAL, PARAMETER :: TRATIO = 0.5  ! above pbl vertical to horizontal ratio

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

!-------------------------------------------------------------------------------
! urban area enhanced turbulence is defined when more of the turbulence
! is assigned to the vertical component at night (tkern > tkerd) than during
! the daytime. The enhancement is determined from DCNet versus NAM TKE values 
! and only applies at night during stable conditions.

  IF(TKERN.GT.TKERD.AND.SLEN.GE.0.0) EE=1.40*EE

  IF(TKERD.GT.0.0.AND.TKERN.GE.0.0)THEN
!    force the ratio of vertical to horizontal turbulence from setup.cfg
     IF(SLEN.GE.0.0)THEN
        TKER=TKERN
     ELSE
        TKER=TKERD
     END IF

     DO K=1,NL
!       from the basic definition TKE = 0.5 ( U'2 + V'2 + W'2 )
!       and TKER is defined as the ratio w'2 / (u'2 + v'2)
        XL=1.0/TKER         
        XX(K)=2.0*EE(K)/(1.0+XL)
        HH(K)=0.50*XL*XX(K)
        EE(K)=0.50*XL*XX(K)
     END DO

  ELSE
!    compute the ratio according to similarity equations

!    top of surface layer (minimum defined for data)
     ZSFC=MAX(ZZ(KSFC),0.1*ZMIX)
!    friction velocity squared
     USTR2=FVEL*FVEL
!    convective velocity scale squared
     WSTR2=WSTR*WSTR

     DO K=1,NL

!       ----------------------------------------------
!       level within the surface layer
        IF(NINT(ZZ(K)).LE.NINT(ZSFC))THEN

           IF(SLEN.GE.0.0)THEN
!             stable layer 
              XX(K)=0.32*EE(K)     ! w'2
              HH(K)=0.74*EE(K)     ! u'2
              EE(K)=0.85*EE(K)     ! v'2

           ELSE
!             unstable layer 
              SL=ZZ(K)/SLEN
              XL=0.41*(WSTR2/USTR2)/(1.0-3.0*SL)**0.67

              XX(K)=2.0*EE(K)/(1.0+XL)
              HH(K)=0.50*XL*XX(K)
              EE(K)=0.50*XL*XX(K)
           END IF
           KPBL=K

!       -----------------------------------------------
!       level is within the PBL
        ELSEIF(ZZ(K).LT.ZMIX)THEN

           IF(SLEN.GE.0.0)THEN
!             stable layer 
              XX(K)=0.32*EE(K)     ! w'2
              HH(K)=0.74*EE(K)     ! u'2
              EE(K)=0.85*EE(K)     ! v'2

           ELSE
!             unstable layer         [1.17 = (1+0.5*R**0.67), where R=0.2] 
!                                    [0.61 = 0.72/1.17] 
              SL=ZZ(K)/ZMIX
              XL=0.61/(SL**0.67)/(1.0-SL)**0.67

              XX(K)=2.0*EE(K)/(1.0+XL)
              HH(K)=0.50*XL*XX(K)
              EE(K)=0.50*XL*XX(K)
           END IF
           KPBL=K

!       ---------------------------------------------------
!       level is within the free troposphere
        ELSE
           XL=1.0/TRATIO       
           XX(K)=2.00*EE(K)/(1.0+XL)
           HH(K)=0.50*XL*XX(K)
           EE(K)=0.50*XL*XX(K)
        END IF
  
!    level loop
     END DO

  END IF

!-------------------------------------------------------------------------------
! proportion horizontal turbulence according to the mean velocity components

  DO K=1,NL

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

END SUBROUTINE stbtke
