!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PRFPRS           processes meteo PRoFile on PReSsure surface
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PROFILE PRESSURE CONVERTS INPUT DATA ON AN ABSOLUTE PRESSURE
!   COORDINATE SYSTEM TO INTERNAL MODEL TERRAIN FOLLOWING COORDINATE
!   INPUT DATA HEIGHTS ASSUMED TO BE RELATIVE TO MSL, INPUT TEMPERATURE
!   ARE ASSUME TO BE VIRTUAL.  ADDITIONAL DIAGNOSTIC LOCAL VARIABLES
!   COMPUTED ARE DENSITY AND POTENTIAL TEMPERATURE.  VERTICAL DATA
!   LINEARLY INTERPOLATED.  REQUIRED VALUES OUTSIDE OF THE INPUT DATA
!   HEIGHT RANGE ARE EXTRAPOLATED AS APPROPRIATE.
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 31 Mar 1998 (RRD)
!                  20 Apr 1999 (RRD) - option to use terrain height field
!                  02 Oct 2000 (RRD) - fortran90 upgrade
!                  02 Jul 2001 (RRD) - ambient temperature saved
!                  24 Jun 2002 (RRD) - converted coordinate from Z* to Z
!                  09 Sep 2002 (RRD) - fortran coding standards
!                  16 Dec 2002 (RRD) - eliminated virtual temp correction
!                  07 Feb 2003 (RRD) - additional diagnostic information
!                  14 Oct 2003 (RRD) - added turbulent kinetic energy
!                  10 Nov 2003 (RRD) - added velocity variance
!                  02 Apr 2004 (RRD) - generic file unit numbers
!
! USAGE:  CALL PRFPRS(VMIX,TKEN,VELV,QFLG,UFLG,TFLG,PFLG,SFLG,ZSFC,P0,U0,V0,
!         T0,Z0,NZ,PSG,Z,U,V,W,T,R,E,H,X,NL,ZMDL,ZSG,PP,UU,VV,WW,TT,ZZ,RH,
!         DEN,AA,EE,HH,XX)
!
!   INPUT ARGUMENT LIST:     see below
!   OUTPUT ARGUMENT LIST:    see below
!   INPUT FILES:             none
!   OUTPUT FILES:            none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE PRFPRS(VMIX,TKEN,VELV,QFLG,UFLG,TFLG,PFLG,SFLG,ZSFC,         &
                  P0,U0,V0,T0,Z0,NZ,PSG,Z,U,V,W,T,R,E,H,X,NL,ZMDL,ZSG,  &
                  PP,UU,VV,WW,TT,ZZ,RH,DEN,AA,EE,HH,XX)

   USE funits

   IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables 
!-------------------------------------------------------------------------------

  LOGICAL,    INTENT(IN)    :: vmix     ! vertical mixing flag
  LOGICAL,    INTENT(IN)    :: tken     ! turbulent kinetic energy
  LOGICAL,    INTENT(IN)    :: velv     ! velocity variance
  LOGICAL,    INTENT(IN)    :: qflg     ! specific humidity indicator
  LOGICAL,    INTENT(IN)    :: uflg     ! low level wind
  LOGICAL,    INTENT(IN)    :: tflg     ! low level temp 
  LOGICAL,    INTENT(IN)    :: pflg     ! surface pressure
  LOGICAL,    INTENT(IN)    :: sflg     ! surface terrain 
  REAL,       INTENT(INOUT) :: zsfc     ! terrain height (m)
  REAL,       INTENT(INOUT) :: p0       ! surface pressure at data terrain (mb)
  REAL,       INTENT(IN)    :: u0       ! low level horizontal wind component
  REAL,       INTENT(IN)    :: v0       ! low level horizontal wind component
  REAL,       INTENT(IN)    :: t0       ! low level temperaure (deg K)
  REAL,       INTENT(IN)    :: z0       ! roughness length (m)
  INTEGER,    INTENT(IN)    :: nz       ! number of input levels
  REAL,       INTENT(IN)    :: psg(:)   ! data sigma-p profile
  REAL,       INTENT(IN)    :: z  (:)   ! pressure data (non-hydrostatic)
  REAL,       INTENT(IN)    :: u  (:)   ! horizontal wind component
  REAL,       INTENT(IN)    :: v  (:)   ! horizontal wind component
  REAL,       INTENT(IN)    :: w  (:)   ! vertical motion component (dp/dt)
  REAL,       INTENT(IN)    :: t  (:)   ! temperature profile (deg K)
  REAL,       INTENT(IN)    :: r  (:)   ! specific humidity (kg/kg)
  REAL,       INTENT(IN)    :: e  (:)   ! turbulent kinetic energy (m2/s2)
  REAL,       INTENT(IN)    :: h  (:)   ! u-component velocity var (m2/s2)
  REAL,       INTENT(IN)    :: x  (:)   ! w-component velocity var (m2/s2)
  INTEGER,    INTENT(IN)    :: nl       ! number of output sigma levels
  REAL,       INTENT(IN)    :: zmdl     ! internal model top (meters)
  REAL,       INTENT(IN)    :: zsg(:)   ! internal model output sigma levels
  REAL,       INTENT(OUT)   :: pp (:)   ! pressure at sigma level (mb)
  REAL,       INTENT(OUT)   :: uu (:)   ! horizontal wind component
  REAL,       INTENT(OUT)   :: vv (:)   ! horizontal wind component
  REAL,       INTENT(OUT)   :: ww (:)   ! vertical motion term (sigma/time)
  REAL,       INTENT(OUT)   :: tt (:)   ! virtual potential temperature (pot K)
  REAL,       INTENT(OUT)   :: zz (:)   ! internal model sigma height (meters)
  REAL,       INTENT(OUT)   :: rh (:)   ! relative humidity fraction (0-1)
  REAL,       INTENT(OUT)   :: den(:)   ! air density (kg/m3)
  REAL,       INTENT(OUT)   :: aa (:)   ! ambient temperature (deg K)
  REAL,       INTENT(OUT)   :: ee (:)   ! turbulent kinetic energy (m2/s2)
  REAL,       INTENT(OUT)   :: hh (:)   ! u-component velocity var (m2/s2)
  REAL,       INTENT(OUT)   :: xx (:)   ! w-component velocity var (m2/s2)

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  REAL,    PARAMETER :: grav  = 9.80616    ! gravity (m/s2)
  REAL,    PARAMETER :: rdry  = 287.04     ! dry air (J/Kg-K)
  REAL,    PARAMETER :: p2jm  = 100.0      ! mb to j/m3
  INTEGER            :: kflag  = 0         ! diagnostic  

  REAL               :: zlvl
  REAL               :: pbot,zbot,rbot,tbot,ubot,vbot,wbot,ebot,hbot,xbot
  REAL               :: ptop,ztop,rtop,ttop,utop,vtop,wtop,etop,htop,xtop
  REAL               :: frac,abot,omega,atop,esat,delz
  INTEGER            :: k,kk,kz,kl

!-------------------------------------------------------------------------------
! external variables
!-------------------------------------------------------------------------------

  SAVE KFLAG

!-------------------------------------------------------------------------------
! compute height of ground surface if not otherwise available
!-------------------------------------------------------------------------------

  IF(.NOT.SFLG)THEN
     IF(P0.GT.PSG(1))THEN
!       build down to surface using lowest upper level temp
        ZSFC=Z(1) - LOG( P0/PSG(1) ) *RDRY*T(1)/GRAV
     ELSE
!       just interpolate between data levels
        K=1
        DO WHILE (PSG(K).GE.P0)
           K=K+1
        END DO
        FRAC=(PSG(K-1)-P0)/(PSG(K-1)-PSG(K))
        ZSFC=FRAC*(Z(K)-Z(K-1))+Z(K-1)
     END IF
  END IF

!-------------------------------------------------------------------------------
! compute the log of the surface pressure
!-------------------------------------------------------------------------------

  IF(PFLG)THEN
!    surface pressure from input data
     PBOT=LOG(P0)
  ELSE
!    use lowest upper level temp to estimate surface pressure
     PBOT=(Z(1)-ZSFC)*GRAV/RDRY/T(1)+LOG(PSG(1))
     P0=EXP(PBOT)
  END IF

!-------------------------------------------------------------------------------
! use adiabatic profile to estimate surface temperature
!-------------------------------------------------------------------------------

  IF(TFLG)THEN
     TBOT=T0
  ELSE
     TBOT=T(1)*(PSG(1)/P0)**(-0.286)
  END IF
! convert temperature to virtual (rrd - 12/16/2002)
! IF(QFLG)TBOT=TBOT*(1.0+0.61*R(1))

!-------------------------------------------------------------------------------
! convert specific humidity to fraction of RH% to fraction
!-------------------------------------------------------------------------------

  IF(QFLG)THEN
     ESAT=EXP(21.4-(5351.0/T(1)))
     RBOT=R(1)*P0/(0.622*ESAT)
  ELSE
     RBOT=R(1)/100.0
  END IF

!-------------------------------------------------------------------------------
! estimate remaining surface values
!-------------------------------------------------------------------------------

  ZBOT=0.0
  WBOT=0.0

! initial vertical index of internal grid
  KL=1
! first output height level on internal grid system
  ZLVL=(ZMDL-ZSFC)*(1.0-ZSG(KL))

!-------------------------------------------------------------------------------
! go through each input level and interpolate to model level
!-------------------------------------------------------------------------------

  DO KZ=1,NZ

!    set low level wind data if available
     IF(KZ.EQ.1)THEN
        IF(UFLG)THEN
           UBOT=U0
           VBOT=V0
        ELSE
!          use neutral log-law when output below data level
           IF(ZLVL.LT.(Z(1)-ZSFC))THEN
              ATOP=LOG((Z(1)-ZSFC)/Z0)
              ABOT=LOG(ZLVL/Z0)
              UBOT=U(1)*ABOT/ATOP
              VBOT=V(1)*ABOT/ATOP
           ELSE
              UBOT=U(1)
              VBOT=V(1)
           END IF
        END IF

        IF(TKEN.OR.VELV)THEN
           EBOT=E(1)
           IF(VELV)THEN
              HBOT=H(1)
              XBOT=X(1)
           END IF
        END IF
     END IF

!    log of pressure at level
     PTOP=LOG(PSG(KZ))
     ZTOP=Z(KZ)-ZSFC
     UTOP=U(KZ)
     VTOP=V(KZ)
     WTOP=W(KZ)

     IF(TKEN.OR.VELV)THEN
        ETOP=E(KZ)
        IF(VELV)THEN
           HTOP=H(KZ)
           XTOP=X(KZ)
        END IF
     END IF

     TTOP=T(KZ)
!    virtual temperature (rrd - 12/16/2002)
!    IF(QFLG)THEN
!       TTOP=(1.0+0.61*R(KZ))*T(KZ)
!    ELSE
!       TTOP=T(KZ)
!    END IF

!    convert to rh fraction
     IF(QFLG)THEN
        ESAT=EXP(21.4-(5351.0/T(KZ)))
        RTOP=R(KZ)*PSG(KZ)/(0.622*ESAT)
     ELSE
        RTOP=R(KZ)/100.0
     END IF

     DO WHILE (ZLVL.LE.ZTOP)
!       height <0 flags extrapolated levels not to be
!       used in stability calculations
        IF(KZ.EQ.1.AND.(.NOT.TFLG).AND.VMIX)THEN
           ZZ(KL)=-ZLVL
        ELSE
           ZZ(KL)=ZLVL
        END IF

!       basic linear interpolation
        FRAC=(ZLVL-ZBOT)/(ZTOP-ZBOT)
        TT(KL)=FRAC*(TTOP-TBOT)+TBOT
        RH(KL)=FRAC*(RTOP-RBOT)+RBOT
        UU(KL)=FRAC*(UTOP-UBOT)+UBOT
        VV(KL)=FRAC*(VTOP-VBOT)+VBOT

        IF(TKEN.OR.VELV)THEN
           EE(KL)=FRAC*(ETOP-EBOT)+EBOT
           IF(VELV)THEN
              HH(KL)=FRAC*(HTOP-HBOT)+HBOT
              XX(KL)=FRAC*(XTOP-XBOT)+XBOT
           END IF
        END IF

!       linear interpolation of log of pressure
        PP(KL)=EXP(FRAC*(PTOP-PBOT)+PBOT)

!       density and potential temperature from local value
        AA(KL)=TT(KL)
        DEN(KL)=P2JM*PP(KL)/(AA(KL)*RDRY)
        TT(KL)=AA(KL)*(1000.0/PP(KL))**0.286

!       vertical velocity term converted to sigma/time
        OMEGA=FRAC*(WTOP-WBOT)+WBOT
        WW(KL)=P2JM*OMEGA/(DEN(KL)*GRAV*(ZMDL-ZSFC))

        KL=KL+1
        IF(KL.GT.NL)RETURN
        ZLVL=(ZMDL-ZSFC)*(1.0-ZSG(KL))
     END DO

!    update bottom definition when above surface
     IF(KL.GT.1)THEN
        ZBOT=ZTOP
        PBOT=PTOP
        TBOT=TTOP
        UBOT=UTOP
        VBOT=VTOP
        WBOT=WTOP
        RBOT=RTOP
        IF(TKEN.OR.VELV)THEN
           EBOT=ETOP
           IF(VELV)THEN
              HBOT=HTOP
              XBOT=XTOP
           END IF
        END IF
     END IF

! input data level loop
  END DO

!-------------------------------------------------------------------------------
! sounding ends but levels remain to be filled
!-------------------------------------------------------------------------------

  IF(KL.LE.NL)THEN
     IF(KFLAG.EQ.0)THEN
        ZLVL=(ZMDL-ZSFC)*(1.0-ZSG(KL))
        WRITE(KF21,*)'--------------------------------------------------------'
        WRITE(KF21,*)' NOTICE prfprs: extrapolation from level (k,m): ',KL,ZLVL
        WRITE(KF21,*)'Input data levels: ',NZ,'    Internal Sigma levels: ',NL
        WRITE(KF21,*)'--------------------------------------------------------'
        KFLAG=1
     END IF

     DO KK=KL,NL
        ZLVL=(ZMDL-ZSFC)*(1.0-ZSG(KK))
        DELZ=ZLVL-ZZ(KK-1)
        ZZ(KK)=ZLVL
        TT(KK)=TT(KK-1)
        AA(KK)=AA(KK-1)
        RH(KK)=RH(KK-1)

!       use previous pressure and density to find pressure
        PP(KK)=PP(KK-1)-DEN(KK-1)*GRAV*DELZ/P2JM
        DEN(KK)=P2JM*PP(KK)/(AA(KK)*RDRY)
        UU(KK)=UU(KK-1)
        VV(KK)=VV(KK-1)

        IF(TKEN.OR.VELV)THEN
           EE(KK)=EE(KK-1)
           IF(VELV)THEN
              HH(KK)=HH(KK-1)
              XX(KK)=XX(KK-1)
           END IF
        END IF

!       diminish magnitude when no data
        WW(KK)=WW(KK-1)/2.0
     END DO
  END IF

END SUBROUTINE prfprs
