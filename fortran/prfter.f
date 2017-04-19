!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PRFTER           processes PRoFile on TERrain surface
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PROFILE TERRAIN CONVERTS METEO INPUT DATA ON TERRAIN FOLLOWING
!   (TYPE 3) TO MODEL TERRAIN FOLLOWING SIGMA LEVELS
!   SEE DOCBLOCK OF PRFPRS FOR A MORE COMPLETE DESCRIPTION.
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 20 Jan 1999 (RRD) - coamps vertical velocity in m/s
!                  20 Apr 1999 (RRD) - added terrain variable
!                  01 May 2000 (RRD) - top of data scaling height
!                  02 Sep 2000 (RRD) - fortran90 upgrade
!                  02 Jul 2001 (RRD) - saved ambient temperature
!                  24 Jun 2002 (RRD) - hysplit coordinate from Z* to Z
!                  09 Sep 2002 (RRD) - fortran coding standards
!                  16 Dec 2002 (RRD) - eliminated virtual temp correction
!                  14 Oct 2003 (RRD) - added turbulent kinetic energy
!                  10 Nov 2003 (RRD) - added velocity variances
!                  02 Apr 2004 (RRD) - generic file unit numbers
!
! USAGE:  CALL PRFTER(VMIX,TKEN,VELV,QFLG,UFLG,TFLG,PFLG,SFLG,ZSFC,ZMDLT,
!              P0,U0,V0,T0,Z0,NZ,PSG,P,U,V,W,T,Q,E,H,X,NL,ZMDL,ZSG,
!              PP,UU,VV,WW,TT,ZZ,RH,DEN,AA,EE,HH,XX)
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

!dwen(20090805) SUBROUTINE PRFTER(VMIX,TKEN,VELV,QFLG,UFLG,TFLG,PFLG,SFLG,ZSFC,ZMDLT,    &
!dwen(20090805)                   P0,U0,V0,T0,Z0,NZ,PSG,P,U,V,W,T,Q,E,H,X,NL,ZMDL,ZSG,   &
!dwen(20090805)                   PP,UU,VV,WW,TT,ZZ,RH,DEN,AA,EE,HH,XX)

 SUBROUTINE PRFTER(VMIX,TKEN,VELV,QFLG,UFLG,TFLG,PFLG,SFLG,ZSFC,ZMDLT,    &
                   P0,U0,V0,T0,Z0,NZ,PSG,P,U,V,W,T,Q,E,H,X,NL,ZMDL,ZSG,   &
                   PP,UU,VV,WW,TT,ZZ,RH,DEN,AA,EE,HH,XX,ramsflg)

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
!dwen(20090823)
!  REAL,       INTENT(IN)    :: zmdlt    ! scaling height for Z* coordinates
  REAL,       INTENT(INout)    :: zmdlt    ! scaling height for Z* coordinates
  REAL,       INTENT(INOUT) :: p0       ! surface pressure at data terrain (mb)
  REAL,       INTENT(IN)    :: u0       ! low level horizontal wind component
  REAL,       INTENT(IN)    :: v0       ! low level horizontal wind component
  REAL,       INTENT(IN)    :: t0       ! low level temperaure (deg K)
  REAL,       INTENT(IN)    :: z0       ! roughness length (m)
  INTEGER,    INTENT(IN)    :: nz       ! number of input levels
  REAL,       INTENT(IN)    :: psg(:)   ! data sigma-p profile
  REAL,       INTENT(IN)    :: p  (:)   ! pressure data (non-hydrostatic)
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

!dwen(20090805) *********************
! JCL(03/27/03): set flag: whether data from RAMS or not 
 logical,     intent(in)    :: ramsflg
!dwen *********************************

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

!dwen(20090823)
!  REAL,    PARAMETER :: grav  = 9.80616    ! gravity (m/s2)
!  REAL,    PARAMETER :: rdry  = 287.04     ! dry air (J/Kg-K)
  real                :: temp
  REAL                :: grav  
  REAL                :: rdry 
  REAL,    PARAMETER :: p2jm  = 100.0      ! mb to j/m3
  INTEGER            :: kflag  = 0         ! diagnostic  

  REAL               :: zlvl
  REAL               :: pbot,zbot,rbot,tbot,ubot,vbot,wbot,ebot,hbot,xbot
  REAL               :: ptop,ztop,rtop,ttop,utop,vtop,wtop,etop,htop,xtop
  REAL               :: sigz,zdat,frac,abot,omega,atop,esat,delz
  INTEGER            :: kk,kz,kl

!-------------------------------------------------------------------------------
! external variables
!-------------------------------------------------------------------------------

  SAVE KFLAG

!-------------------------------------------------------------------------------
!  model top of input data used for scaling terrain surfaces, the value
!  is obtained from meteorological model -  RAMS: 20000  COAMPS: 34800
!  Mesoscale models define a Z* coordinate system where the height
!  above ground of the models levels varies according to terrain height
!  such that  Z* = Ztop ( Zmsl - Zsfc ) / ( Ztop - Zsfc ).  Here we need
!  to convert back from Z* given with the input data to the actual
!  Zagl at each grid point before interpolation of data to internal
!  hysplit grid.  Above equation is then: Zagl = Z* ( 1 - Zsfc / Ztop )
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! compute height of ground surface if not otherwise available
!-------------------------------------------------------------------------------

!dwen(20090805) ******************
! CHG&JCL(10/07/03): synchronize constants with values used in RAMS (see rconstants.h)
      IF(RAMSFLG)THEN
         RDRY=287.0
         GRAV=9.80
      else
         rdry=287.04
         grav=9.80616
      END IF
!dwen      ZMDLT is assigned in RUNSET
! CHG(09/10/03) use RAMS USP model top, now in ZMDL
      IF(RAMSFLG)ZMDLT=ZMDL
!dwen ***************************

!=>compute height of ground surface if not otherwise available
  IF(.NOT.SFLG)THEN
!    computation assumes surface pressure field
!    estimate height of 1013 hPa using lowest upper level temp
     ZSFC=LOG(1013.0/P0)*RDRY*T(1)/GRAV
  END IF

!-------------------------------------------------------------------------------
! compute the surface pressure if needed
!-------------------------------------------------------------------------------

  IF(.NOT.PFLG)THEN
!    use lowest upper level temp to estimate surface pressure
     P0=EXP(LOG(1013.0)-ZSFC*GRAV/RDRY/T(1))
  END IF

!-------------------------------------------------------------------------------
! use adiabatic profile to estimate surface temperature
!-------------------------------------------------------------------------------

  IF(TFLG)THEN
!    available use low level value
     TBOT=T0
  ELSE
!    first compute lowest data level in agl units
     ZDAT=PSG(1)*(1.0-ZSFC/ZMDLT)
!    convert data level to sigma units
     SIGZ=1.0-ZDAT/(ZMDLT-ZSFC)
!    drop adiabatically to surface (ratio in sigma units)
     TBOT=T(1)*SIGZ**(-0.286)
  END IF
! convert all temperatures to virtual (rrd - 12/16/2002)
! IF(QFLG)TBOT=TBOT*(1.0+0.61*Q(1))
!dwen(20090805) ********************
!               copied the following codes from STILT, but I 
!               commented them out. I am not sure what I did is right.
!               HYSPLIT4.9 has eliminated 
!               virtual temperature correction. 
! CHG&JCL(10/07/03) rams: have T in K and moisture as sphu (want: T in K as non-virtual?)
!dwen      IF(QFLG.AND..NOT.RAMSFLG)TBOT=TBOT*(1.0+0.61*Q(1))
!dwen      IF(QFLG.AND.RAMSFLG)TBOT=TBOT*(1.0+0.61*Q(1))
!dwen **************************************************
!
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
! interpolate the input data to the internal model grid
!-------------------------------------------------------------------------------

! initial vertical index for output
  KL=1
! initial internal sigma level for output
  ZLVL=(ZMDL-ZSFC)*(1.0-ZSG(KL))

! set surface values
  ZBOT=0.0
  PBOT=P0
  WBOT=0.0

  DO KZ=1,NZ

!    set low level wind data if available
     IF(KZ.EQ.1)THEN
        IF(UFLG)THEN
           UBOT=U0
           VBOT=V0
        ELSE
!          use neutral log-law when output below data level
!          convert normalized units to agl for level
           ZTOP=PSG(1)*(1.0-ZSFC/ZMDLT)
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

     TTOP=T(KZ)
!    virtual temperature (rrd - 12/16/2002)
!    IF(QFLG)THEN
!       TTOP=(1.0+0.61*Q(KZ))*T(KZ)
!    ELSE
!       TTOP=T(KZ)
!    END IF
!dwen(20090805) ***************************
!       STILT has following codes, I copied them and commented 
!       them out. HYSPLIT4.9 has eliminated virtual temperature
!       correction
! CHG&JCL(10/07/03) rams: have T in K and moisture as sphu (want: T in K as non-virtual?)
!dwen         IF(QFLG.AND..NOT.RAMSFLG)THEN
!dwen            TTOP=(1.0+0.61*Q(KZ))*T(KZ)
!dwen         ELSEIF(QFLG.AND.RAMSFLG)THEN
!dwen            TTOP=(1.0+0.61*Q(KZ))*T(KZ)
!dwen!            TTOP=T(KZ)
!dwen         ELSE
!dwen            TTOP=T(KZ)
!dwen         END IF
!dwen ************************************

     PTOP=P(KZ)
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

!    convert to rh fraction
     IF(QFLG)THEN
        ESAT=EXP(21.4-(5351.0/T(KZ)))
        RTOP=Q(KZ)*PTOP/(0.622*ESAT)
     ELSE
        RTOP=Q(KZ)/100.0
     END IF

!    convert the normalized height to actual agl of level
     ZTOP=PSG(KZ)*(1.0-ZSFC/ZMDLT)

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
        PP(KL)=FRAC*(PTOP-PBOT)+PBOT

        IF(TKEN.OR.VELV)THEN
           EE(KL)=FRAC*(ETOP-EBOT)+EBOT
           IF(VELV)THEN
              HH(KL)=FRAC*(HTOP-HBOT)+HBOT
              XX(KL)=FRAC*(XTOP-XBOT)+XBOT
           END IF
        END IF

!       density and potential temperature from local value
        AA(KL)=TT(KL)
        DEN(KL)=P2JM*PP(KL)/(AA(KL)*RDRY)
!dwen        TT(KL)=AA(KL)*(1000.0/PP(KL))**0.286
!dwen(20090805) ***********************
!              copied from STILT codes
! CHG(10/07/03) use different constants for RAMS (consistent with rconstants.h)
            IF(.NOT.RAMSFLG)TT(KL)=AA(KL)*(1000.0/PP(KL))**0.286
            IF(RAMSFLG)TT(KL)=AA(KL)*(1000.0/PP(KL))**(RDRY/1004.0)
!dwen **************************************



!       vertical velocity term converted to sigma/time
        OMEGA=FRAC*(WTOP-WBOT)+WBOT

        IF(ZMDLT.LE.20000.0)THEN
!          assume data from RAMS where omega in mb/sec
!          and then ds/dt = omega / (row g Ztop)
           WW(KL)=P2JM*OMEGA/(DEN(KL)*GRAV*(ZMDL-ZSFC))
        ELSE
!          assume data from coamps where omega in m/sec
!          and then ds/dt = - dz/dt / Ztop
           WW(KL)=-OMEGA/(ZMDL-ZSFC)
        END IF
!dwen(20090805) ************************
!               copied from STILT codes
! JCL(03/27/03) simple linear interpolation; convert to dsigma/dt
! CHG but don't apply density normalization yet, do this when using winds locally
! WWND corresponds to "w*-flux", not dZ/dt!
! w* = dz*/dt, z*=Zt * (z - Zsurf)/(Zt-Zsurf)= zagl*Zt/(Zt-Zsurf)
! i.e. z* = Zt * (1.0 - sigma_z)
!            IF(RAMSFLG)WW(KL)=-1.0*(FRAC*(WTOP-WBOT)+WBOT)/(ZMDL-ZSFC)
            IF(RAMSFLG)WW(KL)=-1.0*(FRAC*(WTOP-WBOT)+WBOT)/ZMDL
!dwen ********************************
!
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
        WRITE(KF21,*)' NOTICE prfter: extrapolation from level (k,m): ',KL,ZLVL
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

!dwen(20090805) ***********************
!          copied from STILT codes
!           convert potential back to local temperature
! CHG(10/07/03) use different constants for RAMS (consistent with rconstants.h)
            IF(.NOT.RAMSFLG)TEMP=TT(KK)*(PP(KK)/1000.0)**0.286
            IF(RAMSFLG)TEMP=TT(KK)*(PP(KK)/1000.0)**(RDRY/1004.0)
!dwen ******************************************
 

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

END SUBROUTINE prfter
