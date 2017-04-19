!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  STBANL           STaBility ANaLysis from meteo profiles
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   STABILITY ANALYSIS FROM MODEL SURFACE FLUX PARAMETERS
!   OR FROM PROFILE IF THOSE OUTPUTS ARE NOT AVAILABLE.
!   FOR COMPUTATIONAL PURPOSES LEVEL 2 IS THE TOP OF SURFACE
!   SURFACE LAYER.  ROUTINE RETURNS THE HEIGHT NORMALIZED
!   MONIN-OBUKHOV LENGTH, FRICTION TERMS, AND MIXED LAYER DEPTH.
!
! PROGRAM HISTORY LOG:
!   Last Revision: 01 Jul 1997 (RRD)
!                  07 Jul 2000 (RRD) - abs on scalar exchange coeff
!                  28 Jul 2000 (RRD) - generalized flux computation
!                  02 Sep 2000 (RRD) - fortran90 upgrade
!                  16 Mar 2001 (RRD) - argument list change
!                  09 Sep 2002 (RRD) - fortran coding standards
!                  16 Dec 2002 (RRD) - mixed layer top at lower index
!                  07 Feb 2003 (RRD) - added more diagnostic information
!                  22 Sep 2003 (RRD) - vector component for friction vel
!                  07 Nov 2003 (RRD) - mixed layer from TKE profile
!                  02 Dec 2003 (RRD) - defined module for constants
!                  02 Apr 2004 (RRD) - generic file unit numbers
!                  02 Jul 2004 (RRD) - added heat flux logical test
!                  08 Mar 2006 (RRD) - static stability parameter (ss)
!                  19 Jul 2007 (RRD) - compute mixing depth from Tmin
!                  04 Jun 2008 (RRD) - moved variables to common block
!                  12 Jun 2008 (RRD) - mixed layer test variable
!                  07 Aug 2008 (RRD) - replaced isot with kbls
!
! USAGE:  CALL STBANL(KBLS,I,J,KS,KMIX0,KMIXD,MIXD,TKEN,HFLX,EFLX,UFLX,USTAR,
!                     TSTAR,Z0,FMU,FMV,FHS,NL,UU,VV,TT,ZZ,EE,DEN,PBLH,PSIR,
!                     SSP)
!
!   INPUT ARGUMENT LIST:      see below
!   OUTPUT ARGUMENT LIST:     see below
!   INPUT FILES:              none
!   OUTPUT FILES:             none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

!dwen(20090805) SUBROUTINE STBANL(KBLS,I,J,KS,KMIX0,KMIXD,MIXD,TKEN,HFLX,EFLX,UFLX,USTAR, &
!dwen(20090805)                   TSTAR,Z0,FMU,FMV,FHS,NL,UU,VV,TT,ZZ,EE,DEN,PBLH,PSIR,SSP)


subroutine STBANL(KBLS,I,J,KS,KMIX0,KMIXD,MIXD,TKEN,HFLX,EFLX,UFLX,USTAR,   &
                 TSTAR,Z0,FMU,FMV,FHS,NL,UU,VV,TT,ZZ,EE,DEN,PBLH,PSIR,SSP,  &
                 rh,zloc,iconvect,ziscale)


  USE funits
  USE stbcon

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  INTEGER,   INTENT(IN)    :: kbls         ! stability method  
  INTEGER,   INTENT(IN)    :: i,j          ! horizontal index of processed pt
  INTEGER,   INTENT(IN)    :: ks           ! index top of surface layer
  INTEGER,   INTENT(IN)    :: kmix0        ! minimum mixing depth
  INTEGER,   INTENT(IN)    :: kmixd        ! mixed layer depth options
  LOGICAL,   INTENT(IN)    :: mixd         ! mixed layer depth in data file
  LOGICAL,   INTENT(IN)    :: tken         ! tke field avaialble   
  LOGICAL,   INTENT(IN)    :: eflx         ! momentum flux as scalar exchange
  LOGICAL,   INTENT(IN)    :: uflx         ! fluxes available (momentum)
  LOGICAL,   INTENT(IN)    :: hflx         ! fluxes available (heat)      
  LOGICAL,   INTENT(IN)    :: ustar,tstar  ! friction veloc and temp available
  REAL,      INTENT(IN)    :: z0           ! roughness length (m)

                           ! Variable contains momentum flux:
                           ! exchange coefficient (kg/m2-s) when EFLX=TRUE
                           ! friction velocity (m/s)        when USTAR=TRUE
  REAL,      INTENT(INOUT) :: fmu          ! u-component momentum flux (N/m2) 
  REAL,      INTENT(INOUT) :: fmv          ! v-component momentum flux (N/m2) 

                           ! Variable contains sensible heat flux:
                           ! friction temperature           when TSTAR=TRUE
  REAL,      INTENT(IN)    :: fhs          ! sensible heat flux (W/m2)

  INTEGER,   INTENT(IN)    :: nl           ! number of output sigma levels
  REAL,      INTENT(IN)    :: uu  (:)      ! horizontal wind component (m/s)
  REAL,      INTENT(IN)    :: vv  (:)      ! horizontal wind component (m/s)
  REAL,      INTENT(IN)    :: tt  (:)      ! virtual potential temp (deg K)
  REAL,      INTENT(INOUT) :: zz  (:)      ! height at levels (m)
  REAL,      INTENT(IN)    :: ee  (:)      ! turbulent kinetic energy (J/kg)
  REAL,      INTENT(IN)    :: den (:)      ! air density (kg/m3)
  REAL,      INTENT(INOUT) :: pblh         ! mixed layer depth (m)
  REAL,      INTENT(OUT)   :: psir         ! integrated stability function heat
  REAL,      INTENT(OUT)   :: ssp          ! static stability parameter (1/s2) 
!dwen(20090805) ****************
! CHG:(11/19/01) add humidity as input argument, to allow calculation of LCL, LFC and LOC
  real,      intent(in)    :: rh(:)
! CHG:(12/04/01) add ZLOC (limit of convection) as output argument
  real,      intent(out)   :: zloc
! CHG:(9/17/02) add 'ICONVECT' as convection flag 
  integer,   intent(in)    :: iconvect
  real,   intent(in)    :: ziscale
!dwen ******************************

  LOGICAL    :: zflg
  REAL       :: tfhs,phim,phih,etrm,delt,delu,delz,tbar,ubar
  REAL       :: s,z0h,ri,dz,t,v,zl 
  INTEGER    :: k,kdat 

!                             stability analysis results each grid point
  REAL       :: fvel         ! scalar friction velocity (m/s)
  REAL       :: ustr         ! u- friction velocity (m/s)
  REAL       :: vstr         ! v- friction velocity (m/s)
  REAL       :: tstr         ! friction temperature (deg K)
  REAL       :: wstr         ! convective velocity scale (m/s)
  REAL       :: slen         ! Obukhov stability length (m)
  REAL       :: zmix         ! mixed layer depth (m)
  REAL       :: psi          ! integrated stability function heat

!dwen(20090805) ********************
! CHG:(12/01/01) needed for convective parameterization 
! CHG:(9/13/02) out of bound error for local arrays, use max size 
  real,  allocatable :: tk(:),prs(:),e(:),es(:),r(:),rs(:),  &
                        tp(:),tdry(:),tpv(:),tev(:)
! CHG:(11/19/01)
!     dry air (J/Kg-K)   mb to j/m3
!dwen(20090823)      DATA RDRY/287.04/,      P2JM/100.0/
  real, parameter    :: rdry = 287.04
  real, parameter    :: p2jm = 100.0
  integer            :: kret,kmx,lcl,lfc,loc
  real               :: zmixnew,tstk,tex,tt1,riold,   &
                        tlcl,plcl,ttes1,tmax,t_acc,rtsafe,rsat_t
  real               :: l_kmixd
!dwen ******************************


  COMMON /stbcom/ fvel,ustr,vstr,tstr,wstr,slen,zmix,psi

!*******************************************************
!dwen(20090805):allocate variables
  ALLOCATE (tk(nl),prs(nl),e(nl),es(nl),r(nl),          &
            rs(nl),tp(nl),tdry(nl),tpv(nl),tev(nl),     &
            STAT=kret)
  IF(kret.ne.0)THEN
     WRITE(*,*)'ERROR stbanl: memory allocation - ',KRET
     STOP 900
  END IF
!******************************************************

!-------------------------------------------------------------------------------
! analyze sounding to exclude extrapolated data
! go up the sounding to find the first level of real data
! prf??? routines set height to -height if internal model levels
! are below the first data level.  Prevents using interpolated
! levels for stability analysis
!-------------------------------------------------------------------------------

  KDAT=1
  DO WHILE (ZZ(KDAT).LT.0.0)
!    convert back to positive value and save index
     ZZ(KDAT)=-ZZ(KDAT)
     KDAT=KDAT+1
  END DO

  IF(KDAT.GE.NL)THEN
     WRITE(*,*)   'Error: stbanl - no observed data, see MESSAGE file'
     WRITE(KF21,*)'Error: stbanl - no observed data to match internal levels'
     WRITE(KF21,*)'Number of vertical levels (internal, data): ',NL,KDAT
     WRITE(KF21,*)'At horizontal grid position (i,j): ',I,J 
     WRITE(KF21,*)'Level   Height  Temperature'
     DO K=NL,1,-1
        WRITE(KF21,*)K,ZZ(K),TT(K)
     END DO
     STOP 900
  END IF

! level used as top point for sfc layer stability calculations
  KDAT=MAX(2,KS,KDAT)

!-------------------------------------------------------------------------------
! determine z/L from surface flux fields
!-------------------------------------------------------------------------------

  IF((KBLS.EQ.1).AND.(USTAR.OR.UFLX.OR.EFLX).AND.(TSTAR.OR.HFLX))THEN
     IF(USTAR)THEN
!       wind speed
        UBAR=SQRT(UU(KDAT)*UU(KDAT)+VV(KDAT)*VV(KDAT))
!       friction velocity available as input field
        FVEL=MAX(FMU,0.0001)
!       define vector components to give the same scalar
        USTR=FVEL*UU(KDAT)/UBAR
        VSTR=FVEL*VV(KDAT)/UBAR

     ELSEIF(UFLX)THEN
!       both u and v momentum fluxes are available
        USTR=SIGN(MAX(0.01,SQRT(ABS(FMU)/DEN(KDAT))),FMU)
        VSTR=SIGN(MAX(0.01,SQRT(ABS(FMV)/DEN(KDAT))),FMV)
!       scalar friction velocity
        FVEL=SQRT(USTR*USTR+VSTR*VSTR)

     ELSE
!       only scalar exchange coefficient remains available
        USTR=MAX(0.01,SQRT(ABS(FMU)*ABS(UU(KDAT))/DEN(KDAT)))
        VSTR=MAX(0.01,SQRT(ABS(FMU)*ABS(VV(KDAT))/DEN(KDAT)))
        USTR=SIGN(USTR,UU(KDAT))
        VSTR=SIGN(VSTR,VV(KDAT))
!       scalar friction velocity
        FVEL=SQRT(USTR*USTR+VSTR*VSTR)
     END IF

     IF(TSTAR)THEN
!       friction temperature available
        TSTR=SIGN(MAX(0.0001,ABS(FHS)),FHS)
     ELSE
!       friction temperature from sensible heat flux
        TFHS=MAX(-25.0,FHS)
        TSTR=SIGN(MAX(0.0001,ABS(TFHS/(DEN(KS)*CP*FVEL))),-TFHS)
     END IF

!    normalized Monin-Obukhov length (for top of sfc layer)
     ZL=ZZ(KS)*VONK*GRAV*TSTR/(FVEL*FVEL*TT(KS))

!-------------------------------------------------------------------------------
! determine z/L from wind and temperature sounding
!-------------------------------------------------------------------------------

  ELSE
!    Bulk Richardson number uses only non-extrapolated levels
     TBAR=(TT(KDAT)+TT(1))/2.0
     DELZ=ZZ(KDAT)-ZZ(1)
     DELT=TT(KDAT)-TT(1)
     DELU=(UU(KDAT)-UU(1))**2+(VV(KDAT)-VV(1))**2
     DELU=MAX(0.0001,DELU)
     RI=GRAV*DELT*DELZ/TBAR/DELU

!    correction for excessive height (see Golder Ri~Z^2)
!    when data level above assumed surface layer top
     DZ=ZZ(KDAT)*ZZ(KDAT)/ZZ(KS)/ZZ(KS)
     RI=MAX(-1.0, MIN(RI/DZ,1.0))

!------------------------------------------------------------------------------
! convert bulk Ri to z/L at level KS using Hess formula
!-------------------------------------------------------------------------------

     S=LOG(ZZ(KS)/Z0+1.0)

!    roughness length for heat
     Z0H=0.1*Z0
     T=LOG(ZZ(KS)/Z0H +1.0)
     V=LOG(Z0/Z0H)

     IF(RI.GT.0.0.AND.RI.LT.0.08)THEN
        ZL=(-T+10.0*S*RI+SQRT(T*T-20.0*S*T*RI+20.0*S*S*RI))/(10.0*(1.0-5.0*RI))
     ELSEIF(RI.GE.0.08)THEN
        ZL=(B1*S+B2)*RI*RI+(B3*S-B4*V-B5)*RI
     ELSE
        ZL=RI*(S*S/T-B6)
     END IF
     ZL=SIGN(ZL,RI)

     IF(ZL.GE.0.0)THEN
!       stable surface layer Beljaars-Holtslag
        ETRM=B*EXP(-D*ZL)*(1.0-D*ZL+C)
        PHIM=1.0+(A+ETRM)*ZL
        PHIH=PRN*(1.0+(A*SQRT(1.0+A*B*ZL)+ETRM)*ZL)
     ELSE
!       unstable surface layer Betchov-Yaglom / Kadar-Perepelkin
        PHIM=((1.0+0.625*ZL*ZL)/(1.0-7.5*ZL))**0.3333
        PHIH=0.64*((3.0-2.5*ZL)/(1.0-10.0*ZL+50.0*ZL*ZL))**0.3333
     END IF

!    compute friction terms
     FVEL=VONK*ZZ(KS)*SQRT(DELU)/PHIM/DELZ
     TSTR=VONK*ZZ(KS)*DELT/PHIH/DELZ

!    recompute Z/L from friction terms to be consistent
     ZL=ZZ(KS)*VONK*GRAV*TSTR/(FVEL*FVEL*TT(KS))

!    define vector components to give the same scalar
     USTR=0.707107*FVEL
     VSTR=USTR
  END IF

!-------------------------------------------------------------------------------
! check limits on Z/L (-2 <= z/L <= 15 )
!-------------------------------------------------------------------------------

  ZL=MAX(-2.0, MIN(15.0, ZL))
  IF(ZL.EQ.0.0)ZL=0.001
  SLEN=ZZ(KS)/ZL

!-------------------------------------------------------------------------------
! compute integral (psi) required for deposition calculations
!-------------------------------------------------------------------------------

  IF(ZL.LT.-0.001)THEN
!    unstable
     PSI=P1+ZL*(P2+ZL*(P3+ZL*(P4+P5*ZL)))
  ELSEIF(ZL.GE.-0.001.AND.ZL.LT.0.0)THEN
!    neutral
     PSI=-2.7283*ZL
  ELSE
!    stable
     PSI=-(1.0+A*B*ZL)**1.5-B*(ZL-C/D)*EXP(-D*ZL)-B*C/D+1.0
  END IF

!-------------------------------------------------------------------------------
! determine mixed layer depth
! pblh - input and output, initially comes from meteo data file if available     
! zmix - internally computed value depends upon kmixd, saved in common block
!-------------------------------------------------------------------------------

  ! set pbl height computation options based on availability of data:
  ! for kmixd=0 or 2, default to method 1 if data not available
  ! ziscale can be used to select a different default for kmixd=0
  l_kmixd = kmixd
  if (KMIXD .eq. 0 .and. .not. MIXD) l_kmixd=1
  if (KMIXD .eq. 2 .and. .not. TKEN) l_kmixd=1
  IF (l_kmixd .EQ. 0) THEN
!    mixed layer depth from meteo file
     ZMIX=PBLH
  ELSEif (l_kmixd .EQ. 1) THEN
!    find mixed layer depth as the height T(Zi) > 2 + Tmin
     ZMIX=ZZ(NL)
     ZFLG=.FALSE.
     TBAR=MINVAL(TT)
     DO K=NL,1,-1
        IF(.NOT.ZFLG)ZMIX=ZZ(K)
        IF(.NOT.ZFLG.AND.TT(K).LE.TBAR+2.0)ZFLG=.TRUE.
     END DO
  ELSEIF (l_kmixd .EQ. 2) THEN
!    use tke profile to determine mixed layer depth
     K=1
     ZMIX=ZZ(NL)
     DO WHILE (K.LT.NL)
!       mixed layer when tke drops a factor of two or reaches 0.21
        IF((2.0*EE(K+1).LT.EE(K)).OR.(EE(K+1).LE.0.21))THEN
           ZMIX=ZZ(K)
           K=NL
        END IF
        K=K+1
     END DO
  ELSEIF (l_kmixd .EQ. 3) THEN
     !!dwen(20090805) **********************
     ! CHG(8/11/01): replace by more detailed algorithm
     !      use excess temperature as function of virt. pot. tmp
     !      flux, Ustr and Wstr after Holtslag&Boville
     
     ZFLG=.FALSE.
     ZMIX=ZZ(1)
     RI=0.0
     DO K=2,NL
        ! CHG(02/18/04) back to Ri=0.25 as crit. Rich. number
        IF (.NOT.ZFLG.AND.RI.LT.0.25) THEN
           ! CHG(11/26/02) use Ri=0.5 as crit. Rich. number
           !         IF(.NOT.ZFLG.AND.RI.LT.0.5)THEN
           ZMIX=ZZ(K)
           DELZ=ZZ(K)-ZZ(1)
           !           calculate excess temp for convective cases
           IF(ZL.LT.0.0)THEN
              !             surface layer temperature scale
              TSTK=TSTR*DEN(KS)/DEN(K)
              !dwen(20090903):replace ustr with fvel. 
              !        ustr=fvel in STILT codes, but ustr=0.707107*fvel in HYSPLIT4.9 
              !             WSTR=ABS(GRAV*USTR*TSTK*ZMIX/TT(1))**0.3333
              WSTR=ABS(GRAV*fvel*TSTK*ZMIX/TT(1))**0.3333

              !dwen(20090903):replace ustr with fvel. 
              !        ustr=fvel in STILT codes, but ustr=0.707107*fvel in HYSPLIT4.9 
              !             TEX=DABS(8.5*USTR*TSTK/                                   &
              !             (USTR*USTR*USTR+0.6*WSTR*WSTR*WSTR)**0.3333)
              TEX=ABS(8.5*fvel*TSTK/                                   &
                   (fvel*fvel*fvel+0.6*WSTR*WSTR*WSTR)**0.3333)
           ELSE
              TEX=0.0
           END IF
           !dwen(20090904):according to other appears of RI, tt1 should be average 
           !            need to check if this modificaiton is right
           TT1=TT(1)+TEX
           !dwen(20090906):test this            TT1=(TT(1)+TEX)/2.0
           DELT=TT(K)-TT1
           ! CHG(02/18/04) back to version with surface U,V = UU(1) etc.
           ! CHG(11/1902) Try using surface with U,V==0
           !dwen(20090903):no 100*ustr*ustr for DELU calculation according to the 
           !             other appears of DELU in this subroutine
           !            need to check if this modificaiton is right
           !            DELU=(UU(K)-UU(1))*(UU(K)-UU(1))+                           &
           !                (VV(K)-VV(1))*(VV(K)-VV(1))+                           &
           !                100*USTR*USTR
           !dwen                100*fvel*fvel
           DELU=(UU(K)-UU(1))*(UU(K)-UU(1))+                           &
                (VV(K)-VV(1))*(VV(K)-VV(1))                           
           !dwen(20090825)            DELU=DMAX1(DBLE(0.0001),DELU)
           DELU=MAX(0.0001,DELU)
           RIOLD=RI
           RI=GRAV*DELT*DELZ/TT1/DELU
           !           linear interpolation to altitude where Ri=0.25
           ! CHG(02/18/04) back to Ri=0.25 as crit. Rich. number
           ZMIX=(ZZ(K)-ZZ(K-1))/(RI-RIOLD)*(0.25-RIOLD)+ZZ(K-1)
           ! CHG(11/26/02) CHANGE TO Ri=0.5
           !            ZMIX=(ZZ(K)-ZZ(K-1))/(RI-RIOLD)*(0.5-RIOLD)+ZZ(K-1)
           ! CHG(10/30/02): no linear interpolation to altitude where Ri=0.25,
           ! but use NEAREST level
        ELSE
           ZFLG=.TRUE.
        END IF
     END DO

     !  check wether a 'good' mixing layer height has been found
     IF (ZMIX >= ZZ(NL-1)) THEN
        ZMIX=ZZ(NL-1)
        WRITE(*,*) 'Could not find Rib Zmix, reset to ZMIX=',ZMIX
     END IF

  ELSEIF (l_KMIXD .GE. 10) THEN
!    special option where mixing depth defined in namelist file
     ZMIX=float(KMIXD)
  else
     write (*,*) 'STBANL: unexpected value of l_kmixd=',l_kmixd
     STOP 'STBANL: unexpected value of l_kmixd'
  END IF

  if (ziscale .lt. 0.0) then
     ! Special value to denote use of hpbl from met file, if available
     if (mixd) zmix=pblh
  else
     ! JCL:(9/16/02) scaling factor used to prescribe mixed-layer height from file (ZICONTROL)
     !dwen(20090804)               ZMIX=ZMIX*ZISCALE
     ZMIX=zmix*ZISCALE
  end if
! minimum mixed layer height (day or night); make sure it is always above zz(1)
! (this is needed in stb(snd,var,tke) routines where otherwise KPBL would be undefined
  ZMIX=MAX(ZMIX,abs(FLOAT(KMIX0)),zz(1)+(zz(2)-zz(1))/1000.)

! Force zmix to be at a model level if KMIX0 < 0
  if (KMIX0 .lt. 0) then
     ! Find closest model level above k=1
     zmixnew=zz(2)
     do k=3,nl
        IF(ABS(ZZ(K)-ZMIX) .LT. ABS(ZMIXNEW-ZMIX))THEN
           ZMIXNEW=ZZ(K)
        END IF
     enddo
     zmix=zmixnew
  endif

!-------------------------------------------------------------------------------
! compute convective velocity scale
!-------------------------------------------------------------------------------

  IF(ZL.GE.0.0)THEN
     WSTR=0.0
  ELSE
     WSTR=ABS(GRAV*FVEL*TSTR*ZMIX/TT(KS))**0.3333
  END IF

!dwen(20090805) ********************
! CHG:(9/17/02) add 'ICONVECT' as convection flag
      IF(ICONVECT.EQ.0)THEN
         ZLOC=-999.0
      ELSE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCC    FOR CONVECTIVE REDISTRIBUTION   CCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCC
! CHG:(11/19/01) compute lifting condensation level
!     First: get TLCL, then get PLCL, then first level with P<PLCL

!     get TLCL from surface temperature and humidity, after STULL
!     from TT (potential temperature) to TK (temperature in K)
!     (There is a bug in prfsig.f etc, no virtual temperatures
!     calculated for met data with humidity in relative humidity,
!     only for spec. humidity)
!     use parcel starting at lowest level (1), around 10m, has
!     temperature closest to T02M (EDAS and FNL)

!     initialize
!     lifting condensation level
      LCL=NL
!     level of free convection
      LFC=NL
!     limit of convection
      LOC=1
!     index for level just above zmix
      KMX=NL

!     loop over all levels, from bottom to top
      DO K=1,NL
!CCCC get p, T(K), e, esat, and rsat
!     get pressure (mbar), from TT (K) and DEN (kg/m3)
!     p=DEN*TK*RDRY/P2JM; TK=TT*(p/1000)**0.286
!     i.e. p=DEN*TT*(p/1000)**0.286*RDRY/P2JM
         PRS(K)=(DEN(K)*TT(K)*(1.0/1000.0)**0.286*RDRY/P2JM)            &
     &        **(1/(1.0-0.286))
!     get T (K) from pot. temp. (K)
         TK(K)=TT(K)*(PRS(K)/1000)**0.286

!     get saturation vap. pres. (mbar) from T (K)
         IF(TK(K).GE.273.15)THEN
!     liquid
            ES(K)=6.1078*exp(17.29694*(TK(K)-273.16)/(TK(K)-35.86))
         ELSE
!     ice
            ES(K)=exp(23.33086-6111.72784/TK(K)+0.15215*log(TK(K)))
         END IF

!     from RH (rel.humidity fraction) to vapor pressure (mbar)
         E(K)=RH(K)*ES(K)

!     get spec. hum. r (g/g) from e (mbar) and p (mbar)
         R(K)=0.622*E(K)/(PRS(K)-E(K))

!     get sat. spec. hum. rsat (g/g) from e (mbar) and p (mbar)
         RS(K)=0.622*ES(K)/(PRS(K)-ES(K))

!     get dry adiabate ( T in K of dry adiab. lifted parcel...)
         TDRY(K)=TK(1)* (PRS(K)/PRS(1))**0.286

!CCCC end get p, TK(K), e, esat, and rsat

!CCCCCCC get lifting condensation level
         IF(K.EQ.1)THEN
!     get lifting condensation Temperature (K) from level (1) vap. pres. (mbar) and T (K)
!dwen(20090825)            TLCL=2840.0/(3.5*DLOG(TK(1))-DLOG(E(1)/10.0)-7.108)+55.0
            TLCL=2840.0/(3.5*LOG(TK(1))-LOG(E(1)/10.0)-7.108)+55.0
!     get lifting condensation pressure (mbar) from p (mbar) and TLCL (K)
            PLCL=PRS(1)*(TLCL/TK(1))**3.5
         END IF
!     find level LCL just above condensation
         IF(TDRY(K).GE.TLCL)THEN
            LCL=K+1
         END IF
!CCCCCCC end get lifting condensation level

!CCCCCC get parcel temperature TP (cloud temperature) (in K)
!     starting with parcel from lowest model layer (1);

!     below tlcl: dry adiabat;
         IF(TDRY(K).GT.TLCL)THEN
            TP(K)=TDRY(K)
!     Use specific humidity of parcel to get virt. temp
            TPV(K)=TP(K)*(1+0.61*R(1))
         END IF

!     above TLCL: wet adiabat; use pressure levels from model
         IF(TDRY(K).LE.TLCL)THEN

!     get saturation equivalent potential temperature
!     1. of parcel at true lifting condensation level
            TTES1=TT(1)+R(1)*2500.0*TT(1)/TLCL
!     2. of parcel at level K above LCL, as function of TP
!     (which will be the temperature needed to match the
!     saturation equivalent potential temperature in cloud)
!     rtsafe() is used to find TP above LCL
!     Max. value, needed for iteration (min. = TDRY)
            TMAX=TDRY(K)+30.0
            T_ACC=0.01
            TP(K)=rtsafe(TDRY(K)-1,TMAX,PRS(K),TTES1,T_ACC)
!     Use sat. specific humidity of parcel to get virt. temp
            TPV(K)=TP(K)*(1+0.61*RSAT_T(PRS(K),TP(K)))
         END IF
!CCCCCC end get parcel temperature

!     Get virtual temperature of environment
         TEV(K)=TK(K)*(1+0.61*R(K))

!     Check whether parcel is buoyant:
!     first time: LFC (level of free convection)
!     last time (of buoyant layer): LOC (limit of convection)
!     level of free convection
         IF(TPV(K).GT.TEV(K).AND.LFC.EQ.NL)LFC=K
!     level of free convection
         IF(TPV(K).LT.TEV(K).AND.LFC.LT.NL &
     &      .AND.K.GT.LFC.AND.LOC.EQ.1) LOC=K

!     index for level just above zmix
         IF(ZZ(K).GT.ZMIX.AND.KMX.EQ.NL)KMX=K+1

      END DO
!     end loop over all levels, from bottom to top

!     assign altitude for limit of convection
      ZLOC=ZZ(LOC)

!     check if positively buoyant area is not too far above zi,
!     otherwise change ZLOC to neg. value
      IF(ZZ(LFC).GT.ZMIX+2000.0)LOC=1

!     check if no positively buoyant area was found
      IF(LOC.EQ.1)ZLOC=-999.0

!     check if ZLOC is just above ZMIX (when levels are equal)
      IF(LOC.LE.KMX)ZLOC=-999.0

!     TEST ONLY
!     calculate enhancement factor for conv. redistribution area
!     multiplied with cloud fraction, this gives the fraction to be vertically redistributed
!     Get cloud base mass flux in kg/m2
!     Assume 2 m/s updraft averaged over cloudbase
!     Assume 3 h duration (repetition cycle for met fields)
!      FMCB=2.0*3.0*3600.0*DEN(LCL)
!     get integrated density up to cloud base, in kg/m2
!      DENT=0.0
!      DO K=1,LCL-1,1
!         DENT=DENT+DEN(K)*(ZZ(K+1)-ZZ(K))
!      END DO
!     get enhancement factor as ratio FMCB/DENT
!      FENH=FMCB/DENT
!            WRITE(45,*)'FMCB DENT FENH'
!            WRITE(45,*)FMCB,DENT,FENH
! CHG:(12/05/01) calculated FENH are in the range of 8-60, so redistribute whole gridcell independent of cloud cover

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCC END FOR CONVECTIVE REDISTRIBUTION   CCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CHG: (9/17/02) end of IF(ICONVECT.EQ.1)
      END IF
!dwen *****************************

!-------------------------------------------------------------------------------
! low-level static stability parameter
!-------------------------------------------------------------------------------

  TBAR=(TT(KDAT)+TT(1))/2.0
  DELZ=ZZ(KDAT)-ZZ(1)
  DELT=TT(KDAT)-TT(1)
  SSP=GRAV*DELT/TBAR/DELZ

!-------------------------------------------------------------------------------
! remaining variables to return
!-------------------------------------------------------------------------------

  FMU=USTR
  FMV=VSTR
  PSIR=PSI
  PBLH=ZMIX

  DEALLOCATE (tk,prs,e,es,r,rs,tp,tdry,tpv,tev)

END SUBROUTINE stbanl
