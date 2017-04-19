!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PRFECM           PRoFile ECMwf processes ECMWF data
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PROFILE ECMWF CONVERTS A SOUNDING ON HYBRID ECMWF PRESSURE-SIGMA
!   COORDINATES (TYPE 4) TO MODEL TERRAIN FOLLOWING SIGMA LEVEL.
!   INTERPOLATES ALL METEO VARIABLES TO MODEL VERTICAL GRID
!   HYBRID COORDINATES CONSIST OF A PRESSURE OFFSET AND SIGMA VALUE.
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 31 Mar 1998 (RRD)
!                  16 Apr 1999 (RRD) - added terrain variable
!                  20 Apr 1999 (RRD) - option to use terrain height field
!                  02 Oct 2000 (RRD) - fortran90 upgrade
!                  02 Jul 2001 (RRD) - ambient temperature saved
!                  22 Feb 2002 (RRD) - correction to terrain / sfc press
!                  24 Jun 2002 (RRD) - converted coordinate from Z* to Z
!                  09 Sep 2002 (RRD) - fortran coding standards
!                  16 Dec 2002 (RRD) - eliminated virtual temp correction
!                  14 Oct 2003 (RRD) - added turbulent kinetic energy
!                  10 Nov 2003 (RRD) - added velocity variances
!                  02 Apr 2004 (RRD) - generic file unit numbers
!
! USAGE:  CALL PRFECM(VMIX,TKEN,VELV,QFLG,UFLG,TFLG,PFLG,SFLG,ZSFC,P0,U0,
!         V0,T0,Z0,NZ,PSG,ESG,U,V,W,T,Q,E,H,X,NL,ZMDL,ZSG,PP,UU,VV,WW,TT,
!         ZZ,RH,DEN,AA,EE,HH,XX)
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

SUBROUTINE PRFECM(VMIX,TKEN,VELV,QFLG,UFLG,TFLG,PFLG,SFLG,ZSFC,           & 
                  P0,U0,V0,T0,Z0,NZ,PSG,ESG,U,V,W,T,Q,E,H,X,NL,ZMDL,ZSG,  &
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
  REAL,       INTENT(OUT)   :: esg(:)   ! data sigma-p profile (sfc pres corr)
  REAL,       INTENT(IN)    :: u  (:)   ! horizontal wind component
  REAL,       INTENT(IN)    :: v  (:)   ! horizontal wind component
  REAL,       INTENT(IN)    :: w  (:)   ! vertical motion component (dp/dt)
  REAL,       INTENT(IN)    :: t  (:)   ! temperature profile (deg K)
  REAL,       INTENT(IN)    :: q  (:)   ! specific humidity (kg/kg)
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
  REAL               :: sigma,frac,abot,omega,atop,esat,delz,tbar
  INTEGER            :: k,kpres,kk,kz,kl

!-------------------------------------------------------------------------------
! external variables
!-------------------------------------------------------------------------------

  SAVE KFLAG

!-------------------------------------------------------------------------------
! compute the log of the surface pressure
!-------------------------------------------------------------------------------

  IF(PFLG)THEN
!    surface pressure from input data
     PBOT=LOG(P0)
     IF(.NOT.SFLG)THEN
!       Compute height of ground surface if not otherwise available
        ZSFC=LOG(1013.0/P0)*RDRY*T(1)/GRAV
     END IF

  ELSEIF(SFLG)THEN   
!    use lowest upper level temp and terrain to estimate surface pressure
     PBOT=LOG(1013.0)-ZSFC*GRAV/RDRY/T(1)
     P0=EXP(PBOT)

  ELSE
     WRITE(KF21,*)'*ERROR* prfecm: either surface pressure or terrain required'
     WRITE(*,*)   '*ERROR* prfecm: see message file for more information' 
     STOP 900
  END IF

!-------------------------------------------------------------------------------
! remap sigma coordinates based upon surface pressure
! zsg composed of integer offset (10's of mb) plus fraction sigma
! use esg variable to hold corrected sigma from zsg
!-------------------------------------------------------------------------------

  DO K=1,NZ
     KPRES=INT(PSG(K))
     SIGMA=PSG(K)-FLOAT(KPRES)
     ESG(K)=(FLOAT(KPRES)+P0*SIGMA)/P0
  END DO

!-------------------------------------------------------------------------------
! use adiabatic profile to estimate surface temperature
!-------------------------------------------------------------------------------

  IF(TFLG)THEN
!    available use low level value
     TBOT=T0
  ELSE
!    drop adiabatically to surface
     TBOT=T(1)*ESG(1)**(-0.286)
  END IF
! convert all temperatures to virtual (rrd - 12/16/2002)
! IF(QFLG)TBOT=TBOT*(1.0+0.61*Q(1))

!-------------------------------------------------------------------------------
! convert specific humidity to fraction of RH% to fraction
!-------------------------------------------------------------------------------

  IF(QFLG)THEN
     ESAT=EXP(21.4-(5351.0/T(1)))
     RBOT=Q(1)*P0/(0.622*ESAT)
  ELSE
     RBOT=Q(1)/100.0
  END IF

!-------------------------------------------------------------------------------
! estimate remaining surface values
!-------------------------------------------------------------------------------

! initial vertical index
  KL=1
! initial output sigma level
  ZLVL=(ZMDL-ZSFC)*(1.0-ZSG(KL))

! estimate surface values
  ZBOT=0.0
  WBOT=0.0

!-------------------------------------------------------------------------------
! go through each input level and interpolate to model level
!-------------------------------------------------------------------------------

  DO KZ=1,NZ

!    log of pressure at level
     PTOP=LOG(P0*ESG(KZ))

     TTOP=T(KZ)
!    virtual temperature (rrd - 12/16/2002)
!    IF(QFLG)THEN
!       TTOP=(1.0+0.61*Q(KZ))*T(KZ)
!    ELSE
!       TTOP=T(KZ)
!    END IF

!    use layer average for hypsometric equation
     TBAR=0.5*(TTOP+TBOT)
     DELZ=(PBOT-PTOP)*RDRY*TBAR/GRAV
     ZTOP=ZBOT+DELZ

!    only for the first input data level
     IF(KZ.EQ.1)THEN
!       set low level (ground) wind data if available
        IF(UFLG)THEN
           UBOT=U0
           VBOT=V0
        ELSE
!          use neutral log-law when output below data level
           IF(ZLVL.LT.ZTOP)THEN
              ATOP=LOG(ZTOP/Z0)
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

!    convert to rh fraction
     IF(QFLG)THEN
        ESAT=EXP(21.4-(5351.0/T(KZ)))
        RTOP=Q(KZ)*ESG(KZ)*P0/(0.622*ESAT)
     ELSE
        RTOP=Q(KZ)/100.0
     END IF

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

     DO WHILE (ZLVL.LE.ZTOP)
!       height <0 flags extrapolated levels not to be
!       used in stability calculations
        IF(KZ.EQ.1.AND.(.NOT.TFLG).AND.VMIX)THEN
           ZZ(KL)=-ZLVL
        ELSE
           ZZ(KL)=ZLVL
        END IF

!       basic linear interpolation
        FRAC=(ZLVL-ZBOT)/DELZ
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

     PBOT=PTOP
     ZBOT=ZTOP
     RBOT=RTOP
     TBOT=TTOP
     UBOT=UTOP
     VBOT=VTOP
     WBOT=WTOP

     IF(TKEN.OR.VELV)THEN
        EBOT=ETOP
        IF(VELV)THEN
           HBOT=HTOP
           XBOT=XTOP
        END IF
     END IF

  END DO

!-------------------------------------------------------------------------------
! sounding ends but levels remain to be filled
!-------------------------------------------------------------------------------

  IF(KL.LE.NL)THEN
     IF(KFLAG.EQ.0)THEN
        ZLVL=(ZMDL-ZSFC)*(1.0-ZSG(KL))
        WRITE(KF21,*)'--------------------------------------------------------'
        WRITE(KF21,*)' NOTICE prfecm: extrapolation from level (k,m): ',KL,ZLVL
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
        RH(KK)=RTOP

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

END SUBROUTINE prfecm
