!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  STBSND           STaBility SouNDing computes vertical mixing
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   STABILITY SOUNDING COMPUTES THE VERTICAL MIXING COEFFICIENT
!   PROFILE FROM SURFACE STABILITY PARAMETERS AND BULK RICHARDSON
!   NUMBER GIVEN MONIN-OBUKHOV LENGTH AND FRICTION VELOCITY. USING
!   SIMILARITY FLUX PROFILE RELATIONSHIPS, THE VERTICAL DIFFUSIVITY
!   IS COMPUTED WITHIN THE PBL AND IN THE REMAINDER OF THE ATMOSPHERE
!   USING STANDARD MIXING LENGTH THEORY. DIFFUSIVITY IS CONVERTED
!   TO VELOCITY VARIANCE.
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 13 Mar 1997 (RRD)
!                  18 Aug 1998 (RRD) - isotropic turbulence no BL average
!                  02 Sep 2000 (RRD) - fortran90 upgrade
!                  09 Sep 2002 (RRD) - fortran coding standards
!                  10 Nov 2003 (RRD) - TKE and velocity variances 
!                  02 Dec 2003 (RRD) - simplified code, new routines TKE
!                  11 May 2005 (RRD) - revised argument list, limits test
!                  04 Jun 2008 (RRD) - moved variables to common block
!                  15 Aug 2008 (RRD) - enhanced options
!                  01 Oct 2008 (RRD) - mixing adjustments
!
! USAGE:  CALL STBSND(TKERD,TKERN,KZMIX,TVMIX,KSFC,NL,UU,VV,TT,ZZ,XX)
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

!dwen(20090806) SUBROUTINE STBSND(TKERD,TKERN,KZMIX,TVMIX,KSFC,NL,UU,VV,TT,ZZ,EE,HH,XX)
SUBROUTINE STBSND(TKERD,TKERN,KZMIX,TVMIX,KSFC,NL,UU,VV,TT,ZZ,EE,HH,XX,  &
                  xv,ttl,sigmaw,z0)

  USE stbcon

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  INTEGER,   INTENT(IN)    :: tkerd     ! day vert to horiz turb ratio     
  INTEGER,   INTENT(IN)    :: tkern     ! night vert to horiz turb ratio     
  INTEGER,   INTENT(IN)    :: kzmix     ! mixing profile averaging flag     
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

!dwen(20090806) ***********************
!        add xv as vertical mixing coefficient
 real,      intent(out)    :: xv(:)

 real,      intent(out)    :: ttl(:)
 real,      intent(out)    :: sigmaw(:)
 real,      intent(in)     :: z0
!dwen **********************************
!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  INTEGER    :: k,kpbl 
  REAL       :: dtdz,delz,phih0,delu,zsfc,phim,vmix,phim0,etrm,phih 
  REAL       :: uvp2,uvm2,vsum,wght,prd,xl,vl,ri,sl,wm,wh,vscale,hscale 

!                             stability analysis results each grid point
  REAL       :: fvel         ! scalar friction velocity (m/s)
  REAL       :: ustr         ! u- friction velocity (m/s)
  REAL       :: vstr         ! v- friction velocity (m/s)
  REAL       :: tstr         ! friction temperature (deg K)
  REAL       :: wstr         ! convective velocity scale (m/s)
  REAL       :: slen         ! Obukhov stability length (m)
  REAL       :: zmix         ! mixed layer depth (m)
  REAL       :: psi          ! integrated stability function heat

!dwen(20090806) *********************
! JCL:typical value for Coriolis parameter in mid-lats (s-1)
!dwen(20090823)      DATA FF/0.0001/   
  real,  parameter :: ff = 0.0001
  real             :: zh
!dwen ******************************

  COMMON /stbcom/ fvel,ustr,vstr,tstr,wstr,slen,zmix,psi
  COMMON /stblen/ vscale,hscale

!-------------------------------------------------------------------------------

! summation for average mixing
  VSUM   = 0.0
  WGHT   = 0.0

! total flux for diabatic Prandtl number calculation
  WM=(FVEL*FVEL*FVEL+0.6*WSTR*WSTR*WSTR)**0.3333

! top of surface layer (minimum defined for data)
  ZSFC=MAX(ZZ(KSFC),0.1*ZMIX)

!-------------------------------------------------------------------------------
! Analyze each level, where the methodology depends upon the turbulence 
! method selected in the namelist file and the available input data. The default
! analysis depends upon whether the level is in the boundary layer, the
! interfacial area, or the free troposphere. Other options already have some of
! turbulence data available in the input meteorological file.
!-------------------------------------------------------------------------------

  DO K=1,NL
 
!-------------------------------------------------------------------------------
!    level is within the PBL

    IF(ZZ(K).LE.ZMIX)THEN

!       recompute stablility for level height
        SL=ZZ(K)/SLEN
        SL=MAX(-2.0, MIN(10.0,SL))

!dwen(20090806) ************************
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! JCL:      use more sophisticated parameterization for Lagrangian
!           timescale, and directly calculate sigmaw, instead of
!           deriving it from eddy diffusivity
!           based on Hanna(1982) & Stohl's FLEXPART Ver3.1
! JCL:z/h, the alt divided by mixed-layer height--will be used often
            ZH=ZZ(K)/ZMIX
!***********STABLE***********
            IF(SL.GE.0.0)THEN
! JCL:        parameterization for stable conditions
!dwen(20090903): use fvel instead of USTR according HYSPLIT 4.9
!              SIGMAW(K)=1.3*USTR*(1.0-ZH)
              SIGMAW(K)=1.3*fvel*(1.0-ZH)
! CHG(10/30/02): avoid SIGMAW becoming zero
!dwen(20090825)              SIGMAW(K)=DMAX1(SIGMAW(K),DBLE(0.00001))
              SIGMAW(K)=MAX(SIGMAW(K),0.00001)
! JCL:(10/30/02) wrong exponent used previously from FLEXPART V3.1 Manual
!             TTL(K)=0.1*(ZMIX/SIGMAW(K))*(ZH**0.5)
              TTL(K)=0.1*(ZMIX/SIGMAW(K))*(ZH**0.8)
! JCL:(11/1/02) TL has singularity at ZZ(K)==ZMIX that pumps particles out of ML, so to prevent this assign TL of 100.0 at this leve
              IF(ZH.EQ.1.0)TTL(K)=100.0

!***********UNSTABLE*********
            ELSE

! JCL:      parameterization for unstable conditions
! CHG(10/30/02): use real Hanna, not the modified one found in FLEXPART
!              SIGMAW(K)=WSTR*(((1.2*(1.0-0.9*ZH)*(ZH**0.6666))
!     &           +((1.8-1.4*ZH)*(USTR/WSTR)*(USTR/WSTR)))**0.5)
              IF(ZH.LT.0.03)THEN
                 SIGMAW(K)=WSTR*0.96*(3.0*ZH-SLEN/ZMIX)**0.3333
              ELSEIF(ZH.LT.0.4)THEN
!dwen(20090825)                 SIGMAW(K)=WSTR*                                        &
!dwen(20090825)     &                DMIN1(0.96*(3.0*ZH-SLEN/ZMIX)**0.3333,            &
!dwen(20090825)     &                0.763*ZH**0.175)
                 SIGMAW(K)=WSTR*                                        &
                     MIN(0.96*(3.0*ZH-SLEN/ZMIX)**0.3333,            &
                     0.763*ZH**0.175)
              ELSEIF(ZH.LT.0.96)THEN
                 SIGMAW(K)=WSTR*0.722*(1-ZH)**0.207
              ELSE
                 SIGMAW(K)=WSTR*0.37
              END IF

! CHG(10/30/02): avoid SIGMAW becoming zero
!dwen(20090825)              SIGMAW(K)=DMAX1(SIGMAW(K),DBLE(0.00001))
              SIGMAW(K)=MAX(SIGMAW(K),0.00001)

              IF((ZH.LT.0.1).AND.((ZZ(K)-Z0).LT.ABS(SLEN)))THEN
                TTL(K)=0.1*ZZ(K)/(SIGMAW(K)*(0.55+0.38*(ZZ(K)-Z0)/SLEN))
              ELSEIF((ZH.LT.0.1).AND.((ZZ(K)-Z0).GT.ABS(SLEN)))THEN
                TTL(K)=0.59*(ZZ(K)/SIGMAW(K))
              ELSE
                TTL(K)=0.15*(ZMIX/SIGMAW(K))*(1-EXP(-5.0*ZH))
              END IF
!           IF(SL.GE.0.0)THEN
            END IF

!          WRITE(45,*)"stbs K: ",K," SIGMAW ",SIGMAW(K)," ZH ",ZH,
!     :                "  TL ",TTL(K)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!dwen ******************************************

        IF(SL.GE.0.0)THEN
!          stable surface layer Beljaars-Holtslag
           ETRM=B*EXP(-D*SL)*(1.0-D*SL+C)
           PHIM=1.0+(A+ETRM)*SL
           PHIH=PRN*(1.0+(A*SQRT(1.0+A*B*SL)+ETRM)*SL)

        ELSE
!          unstable surface layer Betchov-Yaglom / Kadar-Perepelkin
           PHIM=((1.0+0.625*SL*SL)/(1.0-7.5*SL))**0.3333
           PHIH=0.64*((3.0-2.5*SL)/(1.0-10.0*SL+50.0*SL*SL))**0.3333
        END IF

!       level within the surface layer
        IF(NINT(ZZ(K)).LE.NINT(ZSFC))THEN
           WH=FVEL/PHIH
           VMIX=VONK*WH*ZZ(K)*(1.0-ZZ(K)/ZMIX)*(1.0-ZZ(K)/ZMIX)

!          save the last value of PHI to be valid at z/Zi=0.1
           PHIM0=PHIM
           PHIH0=PHIH

!       level within boundary layer
        ELSE
!          diabatic corrected Prandtl number
           PRD=(PHIH0/PHIM0)+(7.2*VONK*ZZ(K)*WSTR/ZMIX/WM)
           WH=WM/PRD
           VMIX=VONK*WH*ZZ(K)*(1.0-ZZ(K)/ZMIX)*(1.0-ZZ(K)/ZMIX)
        END IF

!       sum values for average (save last index for pbl)
        KPBL=K
        VSUM=VSUM+VMIX
        WGHT=WGHT+1.0

!-------------------------------------------------------------------------------
!    compute mixing through the inversion layer (only for convective case)

!   ELSEIF(NINT(ZZ(K)).EQ.NINT(ZMIX).AND.WSTR.GT.0.0)THEN

!      Betts and Beljaars inversion layer jump model
!      DTDZ=(TT(K)-TT(K-1))/(ZZ(K)-ZZ(K-1))
!      DTDZ=MAX(DTDZ,0.1)
!      VMIX=-0.4*TSTR*FVEL/DTDZ

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
!       new limit test to match RI=20 (5/11/2005)
        XL=MIN(256.0,XL)

!       local stability function
        ETRM=B*EXP(-D*XL)*(1.0+C-D*XL)
        PHIH=PRN*(1.0+(A*SQRT(1.0+A*B*XL)+ETRM)*XL)
        VMIX=VL*VL*ABS(SQRT(DELU)/DELZ)/PHIH
     END IF
  
!    check limits
     VMIX=MAX(VMIN, MIN(VMAX, VMIX))

!dwen(20090806) ********************
!             assign vertical mixing coefficient to xv before converting to turbulent velocity variance
                     xv(k)=vmix
!dwen *******************************

!    Convert mixing coefficient to turbulent velocity variance
      XX(K)=VMIX/VSCALE


!dwen(20090806) ************************
! JCL:   if particle in free troposphere or in inversion layer,
!           then use old parameterization, with constant TL & SIGMAW
!           determined from eddy diffusitivity
!dwen(20090903)         IF(NINT(ZZ(K)).GT.NINT(ZMIX))THEN
         IF(NINT(ZZ(K)).GT.NINT(ZMIX))THEN
            TTL(K)=100.0
!dwen(20090806)            SIGMAW(K)=SQRT(XX(K)/TTL(K))
            SIGMAW(K)=SQRT(Xv(K)/TTL(K))
         END IF
!dwen **********************************

!-------------------------------------------------------------------------------
!    proportion horizontal turbulence according to the mean velocity components
!    Assume W'2 = U'2 + V'2

     HH(K)=0.5*XX(K)
     EE(K)=0.5*XX(K)

!    horizontal turbulent values 
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
     DO K=1,KPBL
!dwen(20090806) **********
        xv(k)=vsum/wght
!dwen ********************
        XX(K)=VSUM/WGHT/VSCALE
     END DO
  ELSEIF(KZMIX.EQ.2)THEN
!    boundary layer adjustment
     DO K=1,KPBL
!dwen(20090806) **********
        xv(k)=xv(k)*tvmix
!dwen ********************
        XX(K)=XX(K)*TVMIX      
     END DO
  ELSEIF(KZMIX.EQ.3)THEN
!    free troposphere adjustment
     KPBL=MIN(KPBL+1,NL)
     DO K=KPBL,NL
!dwen(20090806) **********
        xv(k)=xv(k)*tvmix
!dwen *******************
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


  IF(TKERN.GT.TKERD.AND.SLEN.GE.0.0) XX=1.40*XX

END SUBROUTINE stbsnd
